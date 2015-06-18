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

/**@file   reopt.c
 * @brief  data structures and methods for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/mem.h"
#include "scip/event.h"
#include "scip/scip.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/var.h"
#include "scip/misc.h"
#include "scip/reopt.h"
#include "scip/tree.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/cons.h"
#include "scip/cons_logicor.h"
#include "scip/clock.h"
#include "scip/heur_reoptsols.h"
#include "blockmemshell/memory.h"

#define DEFAULT_MEM_VARAFTERDUAL    10
#define DEFAULT_MEM_VAR             10
#define DEFAULT_MEM_NODES         1000
#define DEFAULT_MEM_RUN            200
#define DEFAULT_MEM_DUALCONS        10

/* event handler properties */
#define EVENTHDLR_NAME         "Reopt"
#define EVENTHDLR_DESC         "node event handler for reoptimization"

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecReopt)
{/*lint --e{715}*/
   SCIP_NODE*          eventnode;
   SCIP_Real           oldbound;
   SCIP_Real           newbound;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(SCIPvarGetType(SCIPeventGetVar(event)) == SCIP_VARTYPE_BINARY);

   eventnode = SCIPgetCurrentNode(scip);
   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);

   assert( eventnode != NULL );

   /* skip if the node is not the focus nodes */
   if( SCIPnodeGetType(eventnode) != SCIP_NODETYPE_FOCUSNODE
    || SCIPnodeGetDepth(eventnode) != SCIPgetEffectiveRootDepth(scip) )
      return SCIP_OKAY;

   SCIPdebugMessage("catch event for node %lld: <%s>: %g -> %g\n", SCIPnodeGetNumber(eventnode),
         SCIPvarGetName(SCIPeventGetVar(event)), SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));

   assert(SCIPisFeasLT(scip, newbound, oldbound) || SCIPisFeasGT(scip, newbound, oldbound));

   SCIP_CALL( SCIPaddReoptDualBndchg(scip, eventnode, SCIPeventGetVar(event), newbound, oldbound) );

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolReopt)
{
   SCIP_VAR** vars;
   int varnr;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(SCIPisReoptEnabled(scip));

   vars = SCIPgetVars(scip);
   for(varnr = 0; varnr < SCIPgetNVars(scip); ++varnr)
   {
      if( SCIPvarGetType(vars[varnr]) == SCIP_VARTYPE_BINARY )
      {
         SCIP_CALL(SCIPcatchVarEvent(scip, vars[varnr], SCIP_EVENTTYPE_GBDCHANGED, eventhdlr, NULL, NULL));
      }
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolReopt)
{
   SCIP_VAR** vars;
   int varnr;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(SCIPisReoptEnabled(scip));

   vars = SCIPgetVars(scip);

   for(varnr = 0; varnr < SCIPgetNVars(scip); ++varnr)
   {
      if( SCIPvarGetType(vars[varnr]) == SCIP_VARTYPE_BINARY )
      {
         SCIP_CALL(SCIPdropVarEvent(scip, vars[varnr], SCIP_EVENTTYPE_GBDCHANGED , eventhdlr, NULL, -1));
      }
   }
   return SCIP_OKAY;
}

/* ---------------- Callback methods of reoptimization methods ---------------- */

/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that sols[pos] array can store at least num entries */
static
SCIP_RETCODE ensureSolsSize(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   num,                /**< minimum number of entries to store */
   int                   runidx              /**< run index for which the memory should checked */
   )
{
   assert(runidx >= 0);
   assert(runidx <= reopt->runsize);

   if( num > reopt->soltree->solssize[runidx] )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->soltree->sols[runidx], reopt->soltree->solssize[runidx],
             newsize) ); /*lint !e866 */
      reopt->soltree->solssize[runidx] = newsize;
   }
   assert(num <= reopt->soltree->solssize[runidx]);

   return SCIP_OKAY;
}

/** ensures, that sols array can store at least num entries */
static
SCIP_RETCODE ensureRunSize(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   num,                /**< minimum number of entries to store */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   if( num >= reopt->runsize )
   {
      int newsize;
      int s;

      newsize = 2*num;
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->soltree->sols, reopt->runsize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->soltree->nsols, reopt->runsize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->soltree->solssize, reopt->runsize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->prevbestsols, reopt->runsize, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&reopt->objs, newsize) );

      for(s = reopt->runsize; s < newsize; s++)
      {
         reopt->prevbestsols[s] = NULL;
         reopt->objs[s] = NULL;
         reopt->soltree->solssize[s] = 0;
         reopt->soltree->nsols[s] = 0;
         reopt->soltree->sols[s] = NULL;
      }

      reopt->runsize = newsize;
   }
   assert(num < reopt->runsize);

   return SCIP_OKAY;
}

/** check the memory of the reoptimization tree and if necessary reallocate */
static
SCIP_RETCODE reopttreeCheckMemory(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopttree != NULL);
   assert(blkmem != NULL);

   if( SCIPqueueIsEmpty(reopttree->openids) )
   {
      unsigned int id;

      assert(reopttree->nreoptnodes == (int)(reopttree->reoptnodessize));

      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopttree->reoptnodes, reopttree->reoptnodessize,
            2*reopttree->reoptnodessize) ); /*lint !e647*/

      for( id = reopttree->reoptnodessize; id < 2*reopttree->reoptnodessize; id++ )
      {
         SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void*) (size_t) id) ); /*lint !e571*/
         reopttree->reoptnodes[id] = NULL;
      }

      reopttree->reoptnodessize *= 2;
   }

   return SCIP_OKAY;
}

/** check allocated memory of a node within the reoptimization tree and if necessary reallocate */
static
SCIP_RETCODE reoptnodeCheckMemory(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   var_mem,            /**< memory for variables */
   int                   child_mem,          /**< memory for child nodes */
   int                   conss_mem           /**< memory for constraints */
   )
{
   assert(reoptnode != NULL);
   assert(blkmem != NULL);
   assert(var_mem >= 0);
   assert(child_mem >= 0);
   assert(conss_mem >= 0);

   /* check allocated memory for variable and bound information */
   if( var_mem > 0 )
   {
      if( reoptnode->varssize == 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->vars, var_mem) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->varbounds, var_mem) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->varboundtypes, var_mem) );
         reoptnode->varssize = var_mem;
      }
      else if( reoptnode->varssize < var_mem )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->vars, reoptnode->varssize, var_mem) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->varbounds, reoptnode->varssize, var_mem) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->varboundtypes, reoptnode->varssize, var_mem) );
         reoptnode->varssize = var_mem;
      }
   }

   /* check allocated memory for child node information */
   if( child_mem > 0 )
   {
      if( reoptnode->allocchildmem == 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->childids, child_mem) );
         reoptnode->nchilds = 0;
         reoptnode->allocchildmem = child_mem;
      }
      else if( reoptnode->allocchildmem < child_mem )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->childids, reoptnode->allocchildmem, child_mem) );
         reoptnode->allocchildmem = child_mem;
      }
   }

   /* check allocated memory for add constraints */
   if( conss_mem > 0 )
   {
      if( reoptnode->consssize == 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->conss, conss_mem) );
         reoptnode->nconss = 0;
         reoptnode->consssize = conss_mem;
      }
      else if( reoptnode->consssize < conss_mem )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->conss, reoptnode->consssize, conss_mem) );
         reoptnode->consssize = conss_mem;
      }
   }

   return SCIP_OKAY;
}

/*
 * local methods
 */

/** returns the number of stored solutions in the subtree induced by @p solnode */
static
int soltreeNInducedSols(
   SCIP_SOLNODE*         solnode             /**< node within the solution tree */
   )
{
   assert(solnode != NULL);

   if( solnode->father == NULL && solnode->rchild == NULL && solnode->lchild == NULL )
      return 0;
   else if( solnode->rchild == NULL && solnode->lchild == NULL )
      return 1;
   else
   {
      if( solnode->rchild == NULL )
         return soltreeNInducedSols(solnode->lchild);
      else if( solnode->lchild == NULL )
         return soltreeNInducedSols(solnode->rchild);
      else
         return soltreeNInducedSols(solnode->rchild) + soltreeNInducedSols(solnode->lchild);
   }
}

/** returns the similarity of the objective functions of two given iterations */
static
SCIP_Real reoptSimilarity(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   obj1_id,            /**< id of one objective function */
   int                   obj2_id,            /**< id of the other objective function */
   SCIP_VAR**            transvars,          /**< transformed problem variables */
   int                   ntransvars          /**< number of transformed problem variables */
   )
{
   SCIP_Real similarity;
   SCIP_Bool onediffertozero;
   int v;

   assert(reopt != NULL);
   assert(transvars != NULL);
   assert(ntransvars >= 0);

   onediffertozero = FALSE;

   /* calc similarity */
   similarity = 0.0;
   for(v = 0; v < ntransvars; v++)
   {
      SCIP_Real c1;
      SCIP_Real c2;
      SCIP_Real lb;
      SCIP_Real ub;

      assert(SCIPvarIsActive(transvars[v]));
      assert(!SCIPvarIsOriginal(transvars[v]));

      lb = SCIPvarGetLbLocal(transvars[v]);
      ub = SCIPvarGetUbLocal(transvars[v]);

      if( SCIPsetIsFeasLT(set, lb, ub) )
      {
         int idx;

         idx = SCIPvarGetProbindex(transvars[v]);
         assert(0 <= idx && idx < ntransvars);

         c1 = reopt->objs[obj1_id][idx];
         c2 = reopt->objs[obj2_id][idx];

         if( c1 != 0 || c2 != 0 )
            onediffertozero = TRUE;

         /* vector product */
         similarity += c1*c2;
      }
   }

   if( !onediffertozero )
      return -2.0;
   else
      return similarity;
}

/** delete the given reoptimization node */
static
SCIP_RETCODE reoptnodeDelete(
   SCIP_REOPTNODE**      reoptnode,          /**< node of the reoptimization tree */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert((*reoptnode) != NULL );
   assert(blkmem != NULL );

   /* delete data for constraints */
   if((*reoptnode)->consssize > 0 )
   {
      int c;

      assert((*reoptnode)->conss != NULL);

      for(c = 0; c < (*reoptnode)->nconss; c++)
      {
         BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->conss[c]->vals, (*reoptnode)->conss[c]->varssize);
         BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->conss[c]->vars, (*reoptnode)->conss[c]->varssize);
         BMSfreeBlockMemory(blkmem, &(*reoptnode)->conss[c]); /*lint !e866*/
      }
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->conss, (*reoptnode)->consssize);
      (*reoptnode)->nconss = 0;
      (*reoptnode)->consssize = 0;
      (*reoptnode)->conss = NULL;
   }

   /* free list of children */
   if( (*reoptnode)->childids != NULL )
   {
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->childids, (*reoptnode)->allocchildmem);
      (*reoptnode)->nchilds = 0;
      (*reoptnode)->allocchildmem = 0;
      (*reoptnode)->childids = NULL;
   }

   /* delete dual constraint */
   if( (*reoptnode)->dualconscur != NULL )
   {
      assert((*reoptnode)->dualconscur->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconscur->vals, (*reoptnode)->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconscur->vars, (*reoptnode)->dualconscur->varssize);
      BMSfreeBlockMemory(blkmem, &(*reoptnode)->dualconscur);
      (*reoptnode)->dualconscur = NULL;
   }

   if( (*reoptnode)->dualconsnex != NULL )
   {
      assert((*reoptnode)->dualconsnex->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconsnex->vals, (*reoptnode)->dualconsnex->varssize);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconsnex->vars, (*reoptnode)->dualconsnex->varssize);
      BMSfreeBlockMemory(blkmem, &(*reoptnode)->dualconsnex);
      (*reoptnode)->dualconsnex = NULL;
   }

   /* free boundtypes */
   if ((*reoptnode)->varboundtypes != NULL )
   {
      assert((*reoptnode)->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->varboundtypes, (*reoptnode)->varssize);
      (*reoptnode)->varboundtypes = NULL;
   }

   /* free bounds */
   if ((*reoptnode)->varbounds != NULL )
   {
      assert((*reoptnode)->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->varbounds, (*reoptnode)->varssize);
      (*reoptnode)->varbounds = NULL;
   }

   /* free variables */
   if ((*reoptnode)->vars != NULL )
   {
      assert((*reoptnode)->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->vars, (*reoptnode)->varssize);
      (*reoptnode)->vars = NULL;
   }

   (*reoptnode)->varssize = 0;

   /* free afterdual-boundtypes */
   if ((*reoptnode)->afterdualvarboundtypes != NULL )
   {
      assert((*reoptnode)->afterdualvarssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->afterdualvarboundtypes, (*reoptnode)->afterdualvarssize);
      (*reoptnode)->afterdualvarboundtypes = NULL;
   }

   /* free afterdual-bounds */
   if ((*reoptnode)->afterdualvarbounds != NULL )
   {
      assert((*reoptnode)->afterdualvarssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->afterdualvarbounds, (*reoptnode)->afterdualvarssize);
      (*reoptnode)->afterdualvarbounds = NULL;
   }

   /* free afterdual-variables */
   if ((*reoptnode)->afterdualvars != NULL )
   {
      assert((*reoptnode)->afterdualvarssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->afterdualvars, (*reoptnode)->afterdualvarssize);
      (*reoptnode)->afterdualvars = NULL;
   }

   (*reoptnode)->afterdualvarssize = 0;

   BMSfreeBlockMemory(blkmem, reoptnode);
   (*reoptnode) = NULL;

   return SCIP_OKAY;
}

/** reset the given reoptimization node */
static
SCIP_RETCODE reoptnodeReset(
   SCIP_REOPTNODE*       reoptnode,          /**< reoptimization node */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reoptnode != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* remove and delete all constraints */
   if( reoptnode->nconss > 0 )
   {
      int c;

      assert(reoptnode->conss != NULL);
      assert(reoptnode->consssize > 0);

      for(c = 0; c < reoptnode->nconss; c++)
      {
         BMSfreeBlockMemoryArray(blkmem, &reoptnode->conss[c]->vals, reoptnode->conss[c]->varssize);
         BMSfreeBlockMemoryArray(blkmem, &reoptnode->conss[c]->vars, reoptnode->conss[c]->varssize);
         BMSfreeBlockMemory(blkmem, &reoptnode->conss[c]); /*lint !e866 */
      }
      reoptnode->nconss = 0;
   }

   /* remove all children */
   if (reoptnode->childids != NULL )
   {
      reoptnode->nchilds = 0;
   }

   /* delete dual constraint */
   if( reoptnode->dualconscur != NULL )
   {
      assert(reoptnode->dualconscur->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->vals, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem ,&reoptnode->dualconscur->vars, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualconscur);
      reoptnode->dualconscur = NULL;
   }

   if( reoptnode->dualconsnex != NULL )
   {
      assert(reoptnode->dualconsnex->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconsnex->vals, reoptnode->dualconsnex->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconsnex->vars, reoptnode->dualconsnex->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualconsnex);
      reoptnode->dualconsnex = NULL;
   }

   reoptnode->parentID = 0;
   reoptnode->nvars = 0;
   reoptnode->nafterdualvars = 0;
   reoptnode->dualfixing = FALSE;
   reoptnode->reopttype = (unsigned int)SCIP_REOPTTYPE_NONE;
   reoptnode->lowerbound = -SCIPsetInfinity(set);

   return SCIP_OKAY;
}

/** delete the node stored at position @p nodeID of the reoptimization tree */
static
SCIP_RETCODE reopttreeDeleteNode(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          id,                 /**< id of a node */
   SCIP_Bool             softreset           /**< delete at the end of the solving process */
   )
{
   assert(reopttree != NULL );
   assert(id < reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id] != NULL );

   if( softreset )
   {
      SCIP_CALL( reoptnodeReset(reopttree->reoptnodes[id], set, blkmem) );
   }
   else
   {
      SCIP_CALL( reoptnodeDelete(&reopttree->reoptnodes[id], blkmem) );
   }

   assert(softreset || reopttree->reoptnodes[id] == NULL);
   assert(reopttree->reoptnodes[id] == NULL || reopttree->reoptnodes[id]->conss == NULL || reopttree->reoptnodes[id]->nconss == 0);
   assert(reopttree->reoptnodes[id] == NULL || reopttree->reoptnodes[id]->childids == NULL || reopttree->reoptnodes[id]->nchilds == 0);

   --reopttree->nreoptnodes;

   return SCIP_OKAY;
}

/** constructor of the solution tree */
static
SCIP_RETCODE createSolTree(
   SCIP_SOLTREE*         soltree,            /**< solution tree */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int s;

   assert(soltree != NULL);

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &soltree->sols, DEFAULT_MEM_RUN) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &soltree->nsols, DEFAULT_MEM_RUN) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &soltree->solssize, DEFAULT_MEM_RUN) );
   for(s = 0; s < DEFAULT_MEM_RUN; s++)
   {
      soltree->nsols[s] = 0;
      soltree->solssize[s] = 0;
      soltree->sols[s] = NULL;
   }

   /* allocate the root node */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &soltree->root) );
   soltree->root->sol = NULL;
   soltree->root->updated = FALSE;
   soltree->root->father = NULL;
   soltree->root->rchild = NULL;
   soltree->root->lchild = NULL;

   return SCIP_OKAY;
}

/** free the given solution node */
static
SCIP_RETCODE soltreefreeNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< the primal */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SOLNODE**        solnode             /**< node within the solution tree */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(primal != NULL || set->stage == SCIP_STAGE_INIT);
   assert(solnode != NULL);
   assert(blkmem != NULL);

   /* free recursive right subtree */
   if( (*solnode)->rchild != NULL )
   {
      SCIP_CALL( soltreefreeNode(reopt, set, primal, blkmem, &(*solnode)->rchild) );
   }
   assert((*solnode)->rchild == NULL);

   /* free recursive left subtree */
   if( (*solnode)->lchild != NULL )
   {
      SCIP_CALL( soltreefreeNode(reopt, set, primal, blkmem, &(*solnode)->lchild) );
   }
   assert((*solnode)->lchild == NULL);

   if( (*solnode)->sol != NULL )
   {
      assert(set->stage == SCIP_STAGE_PROBLEM);

      SCIP_CALL( SCIPsolFree(&(*solnode)->sol, blkmem, primal) );
   }

   /* free this nodes */
   BMSfreeBlockMemoryNull(blkmem, solnode);

   return SCIP_OKAY;
}

/** free the solution tree */
static
SCIP_RETCODE freeSolTree(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          origprimal,         /**< the origprimal */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(reopt->soltree != NULL);
   assert(reopt->soltree->root != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* free all nodes recursive */
   SCIP_CALL( soltreefreeNode(reopt, set, origprimal, blkmem, &reopt->soltree->root) );

   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->sols, reopt->runsize);
   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->nsols, reopt->runsize);
   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->solssize, reopt->runsize);

   BMSfreeMemory(&reopt->soltree);

   return SCIP_OKAY;
}

/** creates and adds a solution node to the solution tree */
static
SCIP_RETCODE solnodeAddChild(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SOLNODE*         father,             /**< father of the node to add */
   SCIP_Bool             rchild,             /**< 0-branch? */
   SCIP_Bool             lchild              /**< 1-branch? */
   )
{
   SCIP_SOLNODE* solnode;

   assert(father != NULL);
   assert(rchild == !lchild);
   assert((rchild && father->rchild == NULL) || (lchild && father->lchild == NULL));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &solnode) );
   solnode->sol = NULL;
   solnode->updated = FALSE;
   solnode->father = father;
   solnode->rchild = NULL;
   solnode->lchild = NULL;

   if( rchild )
      father->rchild = solnode;
   else
      father->lchild = solnode;

   return SCIP_OKAY;
}

/** add a solution to the solution tree */
static
SCIP_RETCODE soltreeAddSol(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,         /**< orig primal */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            vars,               /**< array of original variables */
   SCIP_SOL*             sol,                /**< solution to add */
   SCIP_SOLNODE**        solnode,            /**< current solution node */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             bestsol,            /**< is the solution an optimal (best found) solution */
   SCIP_Bool*            added               /**< pointer to store the result */
   )
{
   SCIP_SOLNODE* cursolnode;
   int varid;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(origprimal != NULL);
   assert(blkmem != NULL);
   assert(vars != NULL);
   assert(sol != NULL);
   assert(solnode != NULL);

   cursolnode = reopt->soltree->root;
   (*added) = FALSE;

   if( set->reopt_savesols > 0 )
   {
      for( varid = 0; varid < nvars; varid++ )
      {
         if( SCIPvarGetType(vars[varid]) == SCIP_VARTYPE_BINARY
          || SCIPvarGetType(vars[varid]) == SCIP_VARTYPE_INTEGER
          || SCIPvarGetType(vars[varid]) == SCIP_VARTYPE_IMPLINT )
         {
            SCIP_Real objval;

            objval = SCIPsolGetVal(sol, set, stat, vars[varid]);
            if( SCIPsetIsFeasEQ(set, objval, 0.0) )
            {
               if( cursolnode->rchild == NULL )
               {
                  SCIP_CALL( solnodeAddChild(blkmem, cursolnode, TRUE, FALSE) );
                  assert(cursolnode->rchild != NULL);
                  (*added) = TRUE;
               }
               cursolnode = cursolnode->rchild;
            }
            else
            {
               assert(SCIPsetIsFeasEQ(set, objval, 1.0));
               if( cursolnode->lchild == NULL )
               {
                  SCIP_CALL( solnodeAddChild(blkmem, cursolnode, FALSE, TRUE) );
                  assert(cursolnode->lchild != NULL);
                  (*added) = TRUE;
               }
               cursolnode = cursolnode->lchild;
            }
         }
      }

      /* the solution was added or is an optimal solution */
      if( *added || bestsol )
      {
         SCIP_SOL* copysol;

         assert(cursolnode->lchild == NULL && cursolnode->rchild == NULL);

         if( *added )
         {
            SCIP_CALL( SCIPsolCopy(&copysol, blkmem, set, stat, origprimal, sol) );
            cursolnode->sol = copysol;
         }
         else
            /* this is a pseudo add; we do not want to save this solution
             * more than once, but we will link this solution to the solution
             * storage of this round */
            (*added) = TRUE;

         if( bestsol )
         {
            assert(reopt->prevbestsols != NULL);
            assert(cursolnode->sol != NULL);

            reopt->prevbestsols[reopt->run-1] = cursolnode->sol;
         }

         (*solnode) = cursolnode;
      }
   }

   return SCIP_OKAY;
}

/** reset all marks 'updated' to FALSE */
static
void soltreeResetMarks(
   SCIP_SOLNODE*         node                /**< node within the solution tree */
   )
{
   assert(node != NULL);

   if( node->rchild != NULL || node->lchild != NULL )
   {
      /* the node is no leaf */
      assert(node->sol == NULL);
      assert(!node->updated);

      if( node->rchild != NULL )
         soltreeResetMarks(node->rchild);
      if( node->lchild != NULL )
         soltreeResetMarks(node->lchild);
   }
   else
   {
      /* the node is a leaf */
      assert(node->father != NULL);
      assert(node->sol != NULL);
      node->updated = FALSE;
   }
}

/** allocate memory for a node within the reoptimization tree */
static
SCIP_RETCODE createReoptnode(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          id                  /**< id of the node to create */
   )
{
   assert(reopttree != NULL );
   assert(id < reopttree->reoptnodessize);

   SCIPdebugMessage("create a reoptnode at ID %u\n", id);

   if(reopttree->reoptnodes[id] == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopttree->reoptnodes[id]) ); /*lint !e866*/

      reopttree->reoptnodes[id]->conss = NULL;
      reopttree->reoptnodes[id]->nconss = 0;
      reopttree->reoptnodes[id]->consssize = 0;
      reopttree->reoptnodes[id]->childids = NULL;
      reopttree->reoptnodes[id]->allocchildmem = 0;
      reopttree->reoptnodes[id]->nchilds = 0;
      reopttree->reoptnodes[id]->nvars = 0;
      reopttree->reoptnodes[id]->nafterdualvars = 0;
      reopttree->reoptnodes[id]->parentID = 0;
      reopttree->reoptnodes[id]->dualfixing = FALSE;
      reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_NONE;
      reopttree->reoptnodes[id]->varssize = 0;
      reopttree->reoptnodes[id]->afterdualvarssize = 0;
      reopttree->reoptnodes[id]->vars = NULL;
      reopttree->reoptnodes[id]->varbounds = NULL;
      reopttree->reoptnodes[id]->varboundtypes = NULL;
      reopttree->reoptnodes[id]->afterdualvars = NULL;
      reopttree->reoptnodes[id]->afterdualvarbounds = NULL;
      reopttree->reoptnodes[id]->afterdualvarboundtypes = NULL;
      reopttree->reoptnodes[id]->dualconscur = NULL;
      reopttree->reoptnodes[id]->dualconsnex = NULL;
      reopttree->reoptnodes[id]->lowerbound = -SCIPsetInfinity(set);
   }
   else
   {
      assert(reopttree->reoptnodes[id]->nvars == 0);
      assert(reopttree->reoptnodes[id]->nafterdualvars == 0);
      reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_NONE;
      reopttree->reoptnodes[id]->lowerbound = -SCIPsetInfinity(set);
   }

   /* increase the counter */
   ++reopttree->nreoptnodes;

   return SCIP_OKAY;
}

/** constructor of the reoptimization tree */
static
SCIP_RETCODE createReopttree(
   SCIP_REOPTTREE*       reopttree,          /**< pointer to the reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   unsigned int id;

   assert(reopttree != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* allocate memory */
   reopttree->reoptnodessize = DEFAULT_MEM_NODES;
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopttree->reoptnodes, reopttree->reoptnodessize) );

   /* initialize the queue of open IDs */
   SCIP_CALL( SCIPqueueCreate(&reopttree->openids, (int)reopttree->reoptnodessize, 2.0) );

   /* fill the queue, but reserve the 0 for the root */
   for( id = 1; id < reopttree->reoptnodessize; id++ )
   {
      reopttree->reoptnodes[id] = NULL;
      SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void*) (size_t) id) ); /*lint !e571*/
   }
   assert(SCIPqueueNElems(reopttree->openids) == (int)(reopttree->reoptnodessize)-1);

   reopttree->nreoptnodes = 0;
   reopttree->ninfsubtrees = 0;
   reopttree->ntotalfeasnodes = 0;
   reopttree->nfeasnodes = 0;
   reopttree->ninfnodes = 0;
   reopttree->ntotalinfnodes= 0;
   reopttree->nprunednodes = 0;
   reopttree->ntotalprunednodes= 0;
   reopttree->ncutoffreoptnodes = 0;
   reopttree->ntotalcutoffreoptnodes = 0;

   /* initialize the root node */
   reopttree->reoptnodes[0] = NULL;
   SCIP_CALL( createReoptnode(reopttree, set, blkmem, 0) );

   return SCIP_OKAY;
}

/** clears the reopttree, e.g., to restart and solve the next problem from scratch */
static
SCIP_RETCODE clearReoptnodes(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             softreset           /**< delete nodes before exit the solving process */
   )
{
   unsigned int id;

   assert(reopttree != NULL );

   /* clear queue with open IDs */
   SCIPqueueClear(reopttree->openids);
   assert(SCIPqueueNElems(reopttree->openids) == 0);

   /* delete all data about nodes */
   for( id = 0; id < reopttree->reoptnodessize; id++ )
   {
      if( reopttree->reoptnodes[id] != NULL )
      {
         SCIP_CALL( reopttreeDeleteNode(reopttree, set, blkmem, id, softreset) );
         assert(reopttree->reoptnodes[id] == NULL || reopttree->reoptnodes[id]->nvars == 0);
      }

      if( id > 0 )
      {
         SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void* ) (size_t ) id) ); /*lint !e571*/
      }
   }
   assert(SCIPqueueNElems(reopttree->openids) == (int)(reopttree->reoptnodessize)-1);

   reopttree->nreoptnodes = 0;

   return SCIP_OKAY;
}

/** free the reoptimization tree */
static
SCIP_RETCODE freeReoptTree(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree data */
   SCIP_SET*             set,                /**< global SCIP settings  */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopttree != NULL);
   assert(blkmem != NULL);

   /* free nodes */
   SCIP_CALL( clearReoptnodes(reopttree, set, blkmem, FALSE) );

   /* free the data */
   BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes, reopttree->reoptnodessize);
   SCIPqueueFree(&reopttree->openids);

   /* free the tree itself */
   BMSfreeMemory(&reopttree);

   return SCIP_OKAY;
}

/** check memory for the constraint to handle bound changes based on dual information */
static
SCIP_RETCODE checkMemDualCons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   size                /**< size which need to be allocated */
   )
{
   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(size > 0);

   if( reopt->dualcons == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopt->dualcons) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->vars, size) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->vals, size) );
      reopt->dualcons->varssize = size;
      reopt->dualcons->nvars = 0;
   }
   else
   {
      if( reopt->dualcons->varssize > 0 )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->dualcons->vars, reopt->dualcons->varssize, size) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->dualcons->vals, reopt->dualcons->varssize, size) );
      }
      else
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->vars, size) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->vals, size) );
         reopt->dualcons->nvars = 0;
      }

      reopt->dualcons->varssize = size;
   }

   return SCIP_OKAY;
}

/** check the memory to store global constraints */
static
SCIP_RETCODE checkMemGlbCons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   mem                 /**< memory which has to be allocated */
   )
{
   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(mem >= 0);

   if( mem > 0 )
   {
      if( reopt->glbconss == NULL )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss, mem) );
         reopt->nglbconss = 0;
         reopt->allocmemglbconss = mem;
      }
      else if( reopt->allocmemglbconss < mem )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->glbconss, reopt->allocmemglbconss, mem) );

         reopt->allocmemglbconss = mem;
      }
   }

   return SCIP_OKAY;
}

/** update the bound changes made by constraint propagations during current iteration; stop saving the bound changes if
 *  we reach a branching decision based on a dual information.
 */
static
SCIP_RETCODE updateConstraintPropagation(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int          id,                 /**< id of the node */
   SCIP_Bool*            transintoorig       /**< transform variables into originals */
   )
{
   int nvars;
   int nconsprops;
   int naddedbndchgs;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(0 < id && id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL );

   /* get the number of all stored constraint propagations */
   SCIPnodeGetNDomchg(node, NULL, &nconsprops, NULL);
   nvars = reopt->reopttree->reoptnodes[id]->nvars;

   if( nconsprops > 0 )
   {
      /* check the memory */
      SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[id], blkmem, nvars + nconsprops, 0, 0) );

      SCIPnodeGetConsProps(node,
            &reopt->reopttree->reoptnodes[id]->vars[nvars],
            &reopt->reopttree->reoptnodes[id]->varbounds[nvars],
            &reopt->reopttree->reoptnodes[id]->varboundtypes[nvars],
            &naddedbndchgs,
            reopt->reopttree->reoptnodes[id]->varssize-nvars);

      assert(nvars + naddedbndchgs <= reopt->reopttree->reoptnodes[id]->varssize);

      reopt->reopttree->reoptnodes[id]->nvars += naddedbndchgs;

      *transintoorig = TRUE;
   }

   return SCIP_OKAY;
}

/** save bound changes made after the first bound change based on dual information, e.g., mode by strong branching.
 *
 *  this method is can be used during reoptimization. if we want to reconstruct a node containing dual bound changes we
 *  have to split the node into the original one and at least one node representing the pruned part. all bound changes,
 *  i.e., (constraint) propagation, made after the first bound change based on dual information are still valid for
 *  the original node after changing the objective function. thus, we can store them for the following iterations.
 *
 *  it should be noted, that these bound change will be found by (constraint) propagation methods anyway after changing
 *  the objective function. do not saving these information and find them again might be useful for conflict analysis.
 */
static
SCIP_RETCODE saveAfterDualBranchings(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int          id,                 /**< id of the node */
   SCIP_Bool*            transintoorig       /**< transform variables into originals */
   )
{
   int nbranchvars;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(0 < id && id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL );

   nbranchvars = 0;

   /* allocate memory */
   if (reopt->reopttree->reoptnodes[id]->afterdualvarssize == 0)
   {
      assert(reopt->reopttree->reoptnodes[id]->afterdualvars == NULL );
      assert(reopt->reopttree->reoptnodes[id]->afterdualvarbounds == NULL );
      assert(reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes == NULL );

      /* allocate block memory for node information */
      reopt->reopttree->reoptnodes[id]->afterdualvarssize = DEFAULT_MEM_VARAFTERDUAL;
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvars), reopt->reopttree->reoptnodes[id]->afterdualvarssize) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvarbounds), reopt->reopttree->reoptnodes[id]->afterdualvarssize) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes), reopt->reopttree->reoptnodes[id]->afterdualvarssize) );
   }

   assert(reopt->reopttree->reoptnodes[id]->afterdualvarssize > 0);
   assert(reopt->reopttree->reoptnodes[id]->nafterdualvars >= 0);

   SCIPnodeGetBdChgsAfterDual(node,
         reopt->reopttree->reoptnodes[id]->afterdualvars,
         reopt->reopttree->reoptnodes[id]->afterdualvarbounds,
         reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes,
         reopt->reopttree->reoptnodes[id]->nafterdualvars,
         &nbranchvars,
         reopt->reopttree->reoptnodes[id]->afterdualvarssize);

   if( nbranchvars > reopt->reopttree->reoptnodes[id]->afterdualvarssize )
   {
      int newsize;
      newsize = nbranchvars + 1;
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvars), reopt->reopttree->reoptnodes[id]->afterdualvarssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvarbounds), reopt->reopttree->reoptnodes[id]->afterdualvarssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes), reopt->reopttree->reoptnodes[id]->afterdualvarssize, newsize) );
      reopt->reopttree->reoptnodes[id]->afterdualvarssize = newsize;

      SCIPnodeGetBdChgsAfterDual(node,
            reopt->reopttree->reoptnodes[id]->afterdualvars,
            reopt->reopttree->reoptnodes[id]->afterdualvarbounds,
            reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes,
            reopt->reopttree->reoptnodes[id]->nafterdualvars,
            &nbranchvars,
            reopt->reopttree->reoptnodes[id]->afterdualvarssize);
   }

   /* the stored variables of this node need to be transformed into the original space */
   if( nbranchvars > 0 )
      *transintoorig = TRUE;

   SCIPdebugMessage(" -> save %d bound changes after dual reductions\n", nbranchvars);

   assert(nbranchvars <= reopt->reopttree->reoptnodes[id]->afterdualvarssize); /* this should be the case */

   reopt->reopttree->reoptnodes[id]->nafterdualvars = nbranchvars;

   return SCIP_OKAY;
}

/** transform variable and bounds back to the original space */
static
SCIP_RETCODE transformIntoOrig(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   unsigned int          id                  /**< id of the node */
   )
{
   int varnr;

   assert(reopt != NULL );
   assert(0 < id && id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL );

   /* transform branching variables and bound changes applied before the first dual reduction */
   for(varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nvars; varnr++)
   {
      SCIP_Real constant;
      SCIP_Real scalar;

      scalar = 1;
      constant = 0;

      if(!SCIPvarIsOriginal(reopt->reopttree->reoptnodes[id]->vars[varnr]))
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->reopttree->reoptnodes[id]->vars[varnr], &scalar, &constant)) ;
         reopt->reopttree->reoptnodes[id]->varbounds[varnr] = (reopt->reopttree->reoptnodes[id]->varbounds[varnr] - constant) / scalar;
      }
      assert(SCIPvarIsOriginal(reopt->reopttree->reoptnodes[id]->vars[varnr]));
   }

   /* transform bound changes affected by dual reduction */
   for(varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nafterdualvars; varnr++)
   {
      SCIP_Real constant;
      SCIP_Real scalar;

      scalar = 1;
      constant = 0;

      if(!SCIPvarIsOriginal(reopt->reopttree->reoptnodes[id]->afterdualvars[varnr]))
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->reopttree->reoptnodes[id]->afterdualvars[varnr], &scalar, &constant)) ;
         reopt->reopttree->reoptnodes[id]->afterdualvarbounds[varnr] = (reopt->reopttree->reoptnodes[id]->afterdualvarbounds[varnr] - constant) / scalar;
      }
      assert(SCIPvarIsOriginal(reopt->reopttree->reoptnodes[id]->afterdualvars[varnr]));
   }

   return SCIP_OKAY;
}

/** search the next node along the root path that was saved by reoptimization */
static
SCIP_RETCODE getLastSavedNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_NODE**           parent,             /**< parent node within the search tree */
   unsigned int*         parentid,           /**< id of the parent node */
   int*                  nbndchgs            /**< number of bound changes */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(reopt->reopttree->reoptnodes != NULL);

   (*nbndchgs) = 0;
   (*parent) = node;

   /* look for a saved parent along the root-path */
   while( SCIPnodeGetDepth(*parent) != 0 )
   {
      int nbranchings;
      int nconsprop;

      nbranchings = 0;
      nconsprop = 0;

      if( set->reopt_saveconsprop )
         SCIPnodeGetNDomchg((*parent), &nbranchings, &nconsprop, NULL);
      else
         SCIPnodeGetNDomchg((*parent), &nbranchings, NULL, NULL);

      (*nbndchgs) = (*nbndchgs) + nbranchings + nconsprop;
      (*parent) = SCIPnodeGetParent(*parent);

      if( SCIPnodeGetDepth(*parent) == 0)
      {
         (*parentid) = 0;
         break;
      }
      else if( SCIPnodeGetReopttype((*parent)) >= SCIP_REOPTTYPE_TRANSIT )
      {
         assert(SCIPnodeGetReoptID((*parent)) < reopt->reopttree->reoptnodessize);
         (*parentid) = SCIPnodeGetReoptID((*parent));
         assert((*parentid) && (*parentid) < reopt->reopttree->reoptnodessize);
         break;
      }
   }

   return SCIP_OKAY;
}

/** adds the id @p childid to the array of child nodes of @p parentid */
static
SCIP_RETCODE reoptAddChild(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          parentid,           /**< id of the parent node */
   unsigned int          childid             /**< id of the child node */
   )
{
   int nchilds;

   assert(reopttree != NULL);
   assert(blkmem != NULL);
   assert(parentid < (unsigned int)reopttree->reoptnodessize);
   assert(childid < (unsigned int)reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[parentid] != NULL);

   nchilds = reopttree->reoptnodes[parentid]->nchilds;

   /* ensure that the array is large enough */
   SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[parentid], blkmem, 0, nchilds+1, 0) );
   assert(reopttree->reoptnodes[parentid]->allocchildmem > nchilds);

   /* add the child */
   reopttree->reoptnodes[parentid]->childids[nchilds] = childid;
   ++reopttree->reoptnodes[parentid]->nchilds;

   SCIPdebugMessage("add ID %u as a child of ID %u.\n", childid, parentid);

   return SCIP_OKAY;
}

/** move all children to the next node (along the root path) stored in the reoptimization tree */
static
SCIP_RETCODE moveChildrenUp(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          nodeid,             /**< id of the node */
   unsigned int          parentid            /**< id of the parent node */
   )
{
   unsigned int childid;
   int varnr;
   int nvars;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(0 < nodeid && nodeid < reopt->reopttree->reoptnodessize);
   assert(parentid < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[nodeid]->childids != NULL);

   /* ensure that enough memory at the parentID is available */
   SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[parentid], blkmem, 0, reopt->reopttree->reoptnodes[parentid]->nchilds + reopt->reopttree->reoptnodes[nodeid]->nchilds, 0) );

   while( reopt->reopttree->reoptnodes[nodeid]->nchilds > 0 )
   {
      int nchilds;

      nchilds = reopt->reopttree->reoptnodes[nodeid]->nchilds;
      childid = reopt->reopttree->reoptnodes[nodeid]->childids[nchilds-1];
      assert(0 < childid && childid < reopt->reopttree->reoptnodessize);

      /* check the memory */
      SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[childid], blkmem, reopt->reopttree->reoptnodes[childid]->nvars + reopt->reopttree->reoptnodes[nodeid]->nvars, 0, 0) );
      assert(reopt->reopttree->reoptnodes[childid]->varssize >= reopt->reopttree->reoptnodes[childid]->nvars + reopt->reopttree->reoptnodes[nodeid]->nvars);

      /* save branching information */
      for(varnr = 0; varnr < reopt->reopttree->reoptnodes[nodeid]->nvars; varnr++)
      {
         nvars = reopt->reopttree->reoptnodes[childid]->nvars;
         reopt->reopttree->reoptnodes[childid]->vars[nvars] = reopt->reopttree->reoptnodes[nodeid]->vars[varnr];
         reopt->reopttree->reoptnodes[childid]->varbounds[nvars] = reopt->reopttree->reoptnodes[nodeid]->varbounds[varnr];
         reopt->reopttree->reoptnodes[childid]->varboundtypes[nvars] = reopt->reopttree->reoptnodes[nodeid]->varboundtypes[varnr];
         ++reopt->reopttree->reoptnodes[childid]->nvars;
      }

      /* update the ID of the parent node */
      reopt->reopttree->reoptnodes[childid]->parentID = parentid;

      /* insert the node as a child */
      SCIP_CALL( reoptAddChild(reopt->reopttree, blkmem, parentid, childid) );

      /* reduce the number of child nodes by 1 */
      --reopt->reopttree->reoptnodes[nodeid]->nchilds;
   }

   return SCIP_OKAY;
}

/** delete all nodes in the subtree induced by nodeID */
static
SCIP_RETCODE deleteChildrenBelow(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          id,                 /**< id of the node */
   SCIP_Bool             delnodeitself,      /**< should the node deleted after deleting the induced subtree? */
   SCIP_Bool             exitsolve           /**< will the solving process end after deletion */
   )
{
   assert(reopttree != NULL );
   assert(blkmem != NULL);
   assert(id < reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id] != NULL);

   /* delete all children below */
   if( reopttree->reoptnodes[id]->childids != NULL && reopttree->reoptnodes[id]->nchilds > 0 )
   {
      SCIPdebugMessage("-> delete subtree induced by ID %u (hard remove = %u)\n", id, exitsolve);

      while( reopttree->reoptnodes[id]->nchilds > 0 )
      {
         int nchilds;
         unsigned int childid;

         nchilds = reopttree->reoptnodes[id]->nchilds;
         childid = reopttree->reoptnodes[id]->childids[nchilds-1];
         assert(0 < childid && childid < reopttree->reoptnodessize);

         SCIP_CALL( deleteChildrenBelow(reopttree, set, blkmem, childid, TRUE, exitsolve) );

         --reopttree->reoptnodes[id]->nchilds;
      }
   }

   /* delete node data*/
   if( delnodeitself )
   {
      SCIP_CALL( reopttreeDeleteNode(reopttree, set, blkmem, id, exitsolve) );
      SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void*) (size_t) id) );
   }

   return SCIP_OKAY;
}

/** replaces a reoptimization nodes by its stored child nodes */
static
SCIP_RETCODE shrinkNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int          id,                 /**< id of the node */
   SCIP_Bool*            shrank,             /**< pointer to store if the node was shrank */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_NODE* parent;
   int ndomchgs;
   unsigned int parentid;

   assert(reopt != NULL);
   assert(node != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   if( reopt->reopttree->reoptnodes[id]->childids != NULL && reopt->reopttree->reoptnodes[id]->nchilds > 0 )
   {
      ndomchgs = 0;
      parentid = 0;
      parent = NULL;

      SCIP_CALL( getLastSavedNode(reopt, set, node, &parent, &parentid, &ndomchgs) );

      assert(parentid != id);
      assert(reopt->reopttree->reoptnodes[parentid] != NULL );
      assert(reopt->reopttree->reoptnodes[parentid]->childids != NULL && reopt->reopttree->reoptnodes[parentid]->nchilds);

      /* check if we want move all children to the next saved node above
       * we want to shrink the path if either
       * - the maximal number of bound changes fix and the number of bound changes is
       *   less than the given threshold set->reopt_maxdiffofnodes
       * or
       * - the number is calculated dynamically and the number of bound changes
       *   is less than log2(SCIPgetNBinVars - (#vars of parent))
       * */
      if( ndomchgs <= set->reopt_maxdiffofnodes )
      {
         int c;

         SCIPdebugMessage(" -> shrink node %lld at ID %u, replaced by %d child nodes.\n", SCIPnodeGetNumber(node), id, reopt->reopttree->reoptnodes[id]->nchilds);

         /* copy the references of child nodes to the parent*/
         SCIP_CALL( moveChildrenUp(reopt, blkmem, id, parentid) );

         /* delete the current node */
         c = 0;
         while( reopt->reopttree->reoptnodes[parentid]->childids[c] != id && c < reopt->reopttree->reoptnodes[parentid]->nchilds )
            ++c;

         assert(reopt->reopttree->reoptnodes[parentid]->childids[c] == id);

         /* replace the childid at position c by the last one */
         reopt->reopttree->reoptnodes[parentid]->childids[c] = reopt->reopttree->reoptnodes[parentid]->childids[reopt->reopttree->reoptnodes[parentid]->nchilds-1];
         --reopt->reopttree->reoptnodes[parentid]->nchilds;

         SCIP_CALL( reopttreeDeleteNode(reopt->reopttree, set, blkmem, id, TRUE) );
         SCIP_CALL( SCIPqueueInsert(reopt->reopttree->openids, (void*) (size_t) id) );

         *shrank = TRUE;

         /* set the reopttype to none */
         SCIPnodeSetReopttype(node, SCIP_REOPTTYPE_NONE);
      }
   }

   return SCIP_OKAY;
}

/** change all reopttypes in the subtree induced by @p nodeID */
static
SCIP_RETCODE changeReopttypeOfSubtree(
   SCIP_REOPTTREE*       reopttree,          /**< reopttree */
   unsigned int          id,                 /**< id of the node */
   SCIP_REOPTTYPE        reopttype           /**< reopttype */
   )
{
   assert(reopttree != NULL);
   assert(id < reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id] != NULL);

   if( reopttree->reoptnodes[id]->childids != NULL && reopttree->reoptnodes[id]->nchilds > 0 )
   {
      unsigned int childid;
      int nchildids;
      int seenids;

      nchildids = reopttree->reoptnodes[id]->nchilds;
      seenids = 0;

      while( seenids < nchildids )
      {
         /* get childID */
         childid = reopttree->reoptnodes[id]->childids[seenids];
         assert(childid < reopttree->reoptnodessize);
         assert(reopttree->reoptnodes[childid] != NULL);

         /* change the reopttype of the node iff the node is neither infeasible nor induces an
          * infeasible subtree and if the node contains no bound changes based on dual decisions */
         if( reopttree->reoptnodes[childid]->reopttype != SCIP_REOPTTYPE_STRBRANCHED
          && reopttree->reoptnodes[childid]->reopttype != SCIP_REOPTTYPE_INFSUBTREE ) /*lint !e641*/
            reopttree->reoptnodes[childid]->reopttype = reopttype; /*lint !e641*/

         /* change reopttype of subtree */
         SCIP_CALL( changeReopttypeOfSubtree(reopttree, childid, reopttype) );

         ++seenids;
      }
   }

   return SCIP_OKAY;
}

/** delete the constraint handling dual information for the current iteration and replace it with the dual constraint
 *  for the next iteration
 */
static
SCIP_RETCODE reoptnodeUpdateDualConss(
   SCIP_REOPTNODE*       reoptnode,          /**< reoptimization node */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reoptnode != NULL);
   assert(blkmem != NULL);

   if( reoptnode->dualconscur != NULL )
   {
      SCIPdebugMessage("reset dual (1) information\n");

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->vals, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->vars, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualconscur);
      reoptnode->dualconscur = NULL;
   }

   if( reoptnode->dualconsnex != NULL )
   {
      reoptnode->dualconscur = reoptnode->dualconsnex;
      reoptnode->dualconsnex = NULL;
   }

   reoptnode->dualfixing = (reoptnode->dualconscur != NULL ? 1 : 0);

   return SCIP_OKAY;
}

/** calculates a (local) similarity of a given node and returns if the subproblem should be solved from scratch */
static
SCIP_RETCODE reoptCheckLocalRestart(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_VAR**            transvars,          /**< transformed variables */
   int                   ntransvars,         /**< number of transformed variables */
   SCIP_Bool*            localrestart        /**< pointer to store if we want to restart solving the (sub)problem */
   )
{
   unsigned int id;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(transvars != NULL);

   /* node == NULL is equivalent to node == root, this case should be handled by SCIPreoptCheckReopt */
   assert(node != NULL);

   *localrestart = FALSE;

   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* set the id to -1 if the node is not part of the reoptimization tree */
   if( SCIPnodeGetDepth(node) > 0 && id == 0 )
      return SCIP_OKAY;

   if( set->reopt_objsimdelay > -1 )
   {
      SCIP_Real sim;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real oldcoef;
      SCIP_Real newcoef;
      int v;
      int idx;

      if( id == 0 )
         reopt->nlocrestarts = 0;

      sim = 0.0;

      /* since the stored objective functions are already normalize the dot-product is equivalent to the similarity */
      for(v = 0; v < ntransvars; v++)
      {
         lb = SCIPvarGetLbLocal(transvars[v]);
         ub = SCIPvarGetUbLocal(transvars[v]);

         /* skip already fixed variables */
         if( SCIPsetIsFeasLT(set, lb, ub) )
         {
            idx = SCIPvarGetProbindex(transvars[v]);
            assert(0 <= idx && idx < ntransvars);

            oldcoef = SCIPreoptGetOldObjCoef(reopt, reopt->run-1, idx);
            newcoef = SCIPreoptGetOldObjCoef(reopt, reopt->run, idx);

            sim += (oldcoef * newcoef);
         }
      }

      /* delete the stored subtree and information about bound changes
       * based on dual information */
      if( SCIPsetIsLT(set, sim, set->reopt_objsimdelay) )
      {
         /* set the flag */
         *localrestart = TRUE;

         ++reopt->nlocrestarts;
         ++reopt->ntotallocrestarts;

         /* delete the stored subtree */
         SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, id, FALSE, FALSE) );

         /* delete the stored constraints; we do this twice in a row because we want to delete both constraints */
         SCIP_CALL( reoptnodeUpdateDualConss(reopt->reopttree->reoptnodes[id], blkmem) );
         SCIP_CALL( reoptnodeUpdateDualConss(reopt->reopttree->reoptnodes[id], blkmem) );
      }

      SCIPdebugMessage(" -> local similarity: %.4f%s\n", sim, *localrestart ? " (solve subproblem from scratch)" : "");
   }

   return SCIP_OKAY;
}

/** save ancestor branching information up to the next stored node */
static
SCIP_RETCODE saveAncestorBranchings(
   SCIP_REOPTTREE*       reopttree,          /**< reoptimization tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the branch and bound tree */
   SCIP_NODE*            parent,             /**< parent node */
   unsigned int          id,                 /**< id of the node */
   unsigned int          parentid            /**< id of the parent node */
   )
{
   int nbranchvars;

   assert(reopttree != NULL );
   assert(node != NULL );
   assert(parent != NULL );
   assert(1 <= id && id < reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id] != NULL );
   assert(parentid < reopttree->reoptnodessize);
   assert(parentid == 0 || reopttree->reoptnodes[parentid] != NULL ); /* if the root is the next saved node, the nodedata can be NULL */

   SCIPdebugMessage(" -> save ancestor branchings\n");

   /* allocate memory */
   if (reopttree->reoptnodes[id]->varssize == 0)
   {
      assert(reopttree->reoptnodes[id]->vars == NULL );
      assert(reopttree->reoptnodes[id]->varbounds == NULL );
      assert(reopttree->reoptnodes[id]->varboundtypes == NULL );

      /* allocate memory for node information */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], blkmem, DEFAULT_MEM_VAR, 0, 0) );
   }

   assert(reopttree->reoptnodes[id]->varssize > 0);
   assert(reopttree->reoptnodes[id]->nvars == 0);

   SCIPnodeGetAncestorBranchingsPart(node, parent,
         reopttree->reoptnodes[id]->vars,
         reopttree->reoptnodes[id]->varbounds,
         reopttree->reoptnodes[id]->varboundtypes,
         &nbranchvars,
         reopttree->reoptnodes[id]->varssize);

   if( nbranchvars >  reopttree->reoptnodes[id]->varssize )
   {
      /* reallocate memory */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], blkmem, nbranchvars, 0, 0) );

      SCIPnodeGetAncestorBranchingsPart(node, parent,
            reopttree->reoptnodes[id]->vars,
            reopttree->reoptnodes[id]->varbounds,
            reopttree->reoptnodes[id]->varboundtypes,
            &nbranchvars,
            reopttree->reoptnodes[id]->varssize);
   }

   assert(nbranchvars <= reopttree->reoptnodes[id]->varssize); /* this should be the case */

   reopttree->reoptnodes[id]->nvars = nbranchvars;

   assert(nbranchvars <= reopttree->reoptnodes[id]->varssize);
   assert(reopttree->reoptnodes[id]->vars != NULL );

   return SCIP_OKAY;
}

/** save additional all constraints that were additionally added to @p node */
static
SCIP_RETCODE saveLocalConssData(
   SCIP_REOPTTREE*       reopttree,          /**< reopttree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the branch and bound tree */
   unsigned int          id                  /**< id of the node*/
   )
{
   SCIP_CONS** addedcons;
   SCIP_Real constant;
   SCIP_Real scalar;
   int var;
   int consnr;
   int naddedconss;
   int addedconsssize;
   int nconss;

   assert(node != NULL );
   assert(reopttree != NULL);
   assert(id < reopttree->reoptnodessize);

   /* save the added pseudo-constraint */
   if(SCIPnodeGetNAddedConss(node) > 0)
   {
      addedconsssize = SCIPnodeGetNAddedConss(node);

      SCIPdebugMessage(" -> save %d locally added constraints\n", addedconsssize);

      /* get memory */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &addedcons, addedconsssize) );
      SCIPnodeGetAddedConss(node, addedcons, &naddedconss, addedconsssize);

      nconss = reopttree->reoptnodes[id]->nconss;

      /* check memory for added constraints */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], blkmem, 0, 0, nconss+naddedconss) );

      for(consnr = 0; consnr < naddedconss; consnr++)
      {
         SCIP_Bool success;

         SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopttree->reoptnodes[id]->conss[nconss]) ); /*lint !e866*/

         success = FALSE;
         SCIP_CALL( SCIPconsGetNVars(addedcons[consnr], set, &reopttree->reoptnodes[id]->conss[nconss]->nvars, &success) );
         assert(success);
         reopttree->reoptnodes[id]->conss[nconss]->varssize = reopttree->reoptnodes[id]->conss[nconss]->nvars;

         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->conss[nconss]->vars,
               reopttree->reoptnodes[id]->conss[nconss]->nvars) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->conss[nconss]->vals,
               reopttree->reoptnodes[id]->conss[nconss]->nvars) );

         success = FALSE;
         SCIP_CALL( SCIPconsGetVars(addedcons[consnr], set, reopttree->reoptnodes[id]->conss[nconss]->vars, reopttree->reoptnodes[id]->conss[nconss]->nvars, &success) );
         assert(success);

         if( strcmp("sepasol", SCIPconsGetName(addedcons[consnr])) == 0 )
            reopttree->reoptnodes[id]->conss[nconss]->constype = REOPT_CONSTYPE_SEPASOLUTION;
         else if( strcmp("infsubtree", SCIPconsGetName(addedcons[consnr])) == 0 )
            reopttree->reoptnodes[id]->conss[nconss]->constype = REOPT_CONSTYPE_INFSUBTREE;
         else if( strcmp("splitcons", SCIPconsGetName(addedcons[consnr])) == 0 )
            reopttree->reoptnodes[id]->conss[nconss]->constype = REOPT_CONSTYPE_STRBRANCHED;

         assert(reopttree->reoptnodes[id]->conss[nconss]->constype == REOPT_CONSTYPE_SEPASOLUTION
             || reopttree->reoptnodes[id]->conss[nconss]->constype == REOPT_CONSTYPE_INFSUBTREE
             || reopttree->reoptnodes[id]->conss[nconss]->constype == REOPT_CONSTYPE_STRBRANCHED);

         for(var = 0; var < reopttree->reoptnodes[id]->conss[nconss]->nvars; var++)
         {
            constant = 0;
            scalar = 1;

            if(!SCIPvarIsOriginal(reopttree->reoptnodes[id]->conss[nconss]->vars[var]))
            {
               if(SCIPvarIsNegated(reopttree->reoptnodes[id]->conss[nconss]->vars[var]))
               {
                  SCIP_CALL(SCIPvarGetOrigvarSum(&reopttree->reoptnodes[id]->conss[nconss]->vars[var], &scalar, &constant));
                  reopttree->reoptnodes[id]->conss[nconss]->vals[var] = 1;
               }
               else
               {
                  SCIP_CALL(SCIPvarGetOrigvarSum(&reopttree->reoptnodes[id]->conss[nconss]->vars[var], &scalar, &constant));
                  reopttree->reoptnodes[id]->conss[nconss]->vals[var] = 0;
               }
               assert(reopttree->reoptnodes[id]->conss[nconss]->vars[var] != NULL );
            }
            assert(SCIPvarIsOriginal(reopttree->reoptnodes[id]->conss[nconss]->vars[var]));
         }

         /* increase the counter for added constraints */
         ++reopttree->reoptnodes[id]->nconss;
         ++nconss;
      }

      assert(reopttree->reoptnodes[id]->nconss == naddedconss);
      SCIPsetFreeBufferArray(set, &addedcons);
   }

   return SCIP_OKAY;
}

/** collect all bound changes based on dual information
 *
 *  if the bound changes are global, all information are already stored because they were caught by the event handler.
 *  otherwise, we have to use SCIPnodeGetDualBoundchgs.
 *
 *  afterwards, we check if the constraint will be added in the next iteration or after splitting the node.
 */
static
SCIP_RETCODE collectDualInformation(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int          id,                 /**< id of the node */
   SCIP_REOPTTYPE        reopttype           /**< reopttype */
   )
{
   SCIP_Real constant;
   SCIP_Real scalar;
   SCIP_Bool cons_is_next;
   int nbndchgs;
   int v;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id]->dualfixing);
   assert(node != NULL);
   assert(blkmem != NULL);

   cons_is_next = TRUE;

   /* first case, all bound changes were global */
   if( reopt->currentnode == SCIPnodeGetNumber(node) && reopt->dualcons != NULL && reopt->dualcons->nvars > 0 )
   {
      nbndchgs = reopt->dualcons->nvars;
   }
   else
   {
      assert(reopt->currentnode == SCIPnodeGetNumber(node));

      /* get the number of bound changes based on dual information */
      nbndchgs = SCIPnodeGetNDualBndchgs(node);

      /* ensure that enough memory is allocated */
      SCIP_CALL( checkMemDualCons(reopt, blkmem, nbndchgs) );

      /* collect the bound changes */
      SCIPnodeGetDualBoundchgs(node,
            reopt->dualcons->vars,
            reopt->dualcons->vals,
            &nbndchgs,
            reopt->dualcons->varssize);

      assert(nbndchgs <= reopt->dualcons->varssize);

      reopt->dualcons->nvars = nbndchgs;

      /* transform the variables into the original space */
      for(v = 0; v < nbndchgs; v++)
      {
         constant = 0.0;
         scalar = 1.0;

         SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->dualcons->vars[v], &scalar, &constant) );
         reopt->dualcons->vals[v] = (reopt->dualcons->vals[v] - constant) / scalar;

         assert(SCIPvarIsOriginal(reopt->dualcons->vars[v]));
      }
   }

   assert(nbndchgs > 0);

   /* due to the strong branching initialization it can be possible that two
    * constraints handling dual information are stored at the same time.
    * during reoptimizing a node we add the constraint stored at dualconscur only,
    * i.e, if dualconscur is not NULL, we need to store the constraint the
    * constraint for the next iteration at dualconsnex because the constraint
    * stored at dualconscur is needed to split the constraint in the current
    * iteration.
    */
   if( reopt->reopttree->reoptnodes[id]->dualconscur != NULL )
   {
      assert(reopt->reopttree->reoptnodes[id]->dualconsnex == NULL);
      cons_is_next = FALSE;
   }
   assert((cons_is_next && reopt->reopttree->reoptnodes[id]->dualconscur == NULL)
       || (!cons_is_next && reopt->reopttree->reoptnodes[id]->dualconsnex == NULL));

   /* the constraint will be added next */
   if( cons_is_next )
   {
      assert(reopt->reopttree->reoptnodes[id]->dualconscur == NULL);
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopt->reopttree->reoptnodes[id]->dualconscur) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconscur->vars,
            reopt->dualcons->vars, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconscur->vals,
            reopt->dualcons->vals, nbndchgs) );

      reopt->reopttree->reoptnodes[id]->dualconscur->nvars = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconscur->varssize = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconscur->constype = reopttype == SCIP_REOPTTYPE_STRBRANCHED ? REOPT_CONSTYPE_STRBRANCHED : REOPT_CONSTYPE_INFSUBTREE;

      SCIPdebugMessage(" -> save dual information of type 1: node %lld, nvars %d, constype %d\n",
            SCIPnodeGetNumber(node), reopt->reopttree->reoptnodes[id]->dualconscur->nvars,
            reopt->reopttree->reoptnodes[id]->dualconscur->constype);
   }
   /* the constraint will be added after next */
   else
   {
      assert(reopt->reopttree->reoptnodes[id]->dualconsnex == NULL);
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopt->reopttree->reoptnodes[id]->dualconsnex) );
      reopt->reopttree->reoptnodes[id]->dualconsnex->nvars = -1;

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconsnex->vars, reopt->dualcons->vars, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconsnex->vals, reopt->dualcons->vals, nbndchgs) );
      reopt->reopttree->reoptnodes[id]->dualconsnex->nvars = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconsnex->varssize = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconsnex->constype = reopttype == SCIP_REOPTTYPE_STRBRANCHED ? REOPT_CONSTYPE_STRBRANCHED : REOPT_CONSTYPE_INFSUBTREE;

      SCIPdebugMessage(" -> save dual information of type 2: node %lld, nvars %d, constype %d\n",
            SCIPnodeGetNumber(node), reopt->reopttree->reoptnodes[id]->dualconsnex->nvars,
            reopt->reopttree->reoptnodes[id]->dualconsnex->constype);
   }


   return SCIP_OKAY;
}

/** adds a node of the branch and bound tree to the reoptimization tree */
static
SCIP_RETCODE addNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< current node */
   SCIP_REOPTTYPE        reopttype,          /**< reason for storing the node*/
   SCIP_Bool             saveafterdual,      /**< save branching decisions after the first dual */
   SCIP_Bool             isrootnode,         /**< node is the root node */
   SCIP_Real             lowerbound          /**< lower bound of the node */
   )
{
   SCIP_NODE* parent;
   SCIP_Bool shrank;
   unsigned int id;
   unsigned int parentid;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);

   parentid = 0;
   parent = NULL;
   shrank = FALSE;

   if( set->reopt_maxsavednodes == 0 )
      return SCIP_OKAY;

   assert(reopttype == SCIP_REOPTTYPE_TRANSIT
       || reopttype == SCIP_REOPTTYPE_INFSUBTREE
       || reopttype == SCIP_REOPTTYPE_STRBRANCHED
       || reopttype == SCIP_REOPTTYPE_LOGICORNODE
       || reopttype == SCIP_REOPTTYPE_LEAF
       || reopttype == SCIP_REOPTTYPE_PRUNED
       || reopttype == SCIP_REOPTTYPE_FEASIBLE);

   /* start clock */
   SCIPclockStart(reopt->savingtime, set);

   /* the node was created by reoptimization, i.e., we need to update the
    * stored data */
   if( SCIPnodeGetReoptID(node) >= 1 )
   {
      SCIP_Bool transintoorig;

      assert(reopttype != SCIP_REOPTTYPE_LEAF);
      assert(!isrootnode);

      id = SCIPnodeGetReoptID(node);
      assert(id < reopt->reopttree->reoptnodessize);
      assert(reopt->reopttree->reoptnodes[id] != NULL);

      SCIPdebugMessage("update node %lld at ID %u:\n", SCIPnodeGetNumber(node), id);

      transintoorig = FALSE;

      /* store in*/
      if( saveafterdual )
      {
         SCIP_CALL( saveAfterDualBranchings(reopt, blkmem, node, id, &transintoorig) );
      }

      /* update constraint propagations */
      if( set->reopt_saveconsprop )
      {
         SCIP_CALL( updateConstraintPropagation(reopt, blkmem, node, id, &transintoorig) );
      }

      /* ensure that all variables are original */
      if( transintoorig )
      {
         SCIP_CALL( transformIntoOrig(reopt, id) );
      }

      /* update the lowerbound */
      if( !SCIPsetIsEQ(set, REALABS(lowerbound), SCIPsetInfinity(set)) )
         reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

#ifdef SCIP_DEBUG
         int varnr;

         SCIPdebugMessage(" -> nvars: %d, ncons: %d, parentID: %d, reopttype: %d\n",
               reopt->reopttree->reoptnodes[id]->nvars,
               reopt->reopttree->reoptnodes[id]->nconss,
               reopt->reopttree->reoptnodes[id]->parentID, reopttype);
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage(" -> saved variables:\n");
         for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nvars; varnr++)
         {
            SCIPdebugMessage("  <%s> %s %g\n", SCIPvarGetName(reopt->reopttree->reoptnodes[id]->vars[varnr]),
                  reopt->reopttree->reoptnodes[id]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                  "=>" : "<=", reopt->reopttree->reoptnodes[id]->varbounds[varnr]);
         }
         for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nafterdualvars; varnr++)
         {
            SCIPdebugMessage("  <%s> %s %g (after dual red.)\n", SCIPvarGetName(reopt->reopttree->reoptnodes[id]->afterdualvars[varnr]),
                  reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                  "=>" : "<=", reopt->reopttree->reoptnodes[id]->afterdualvarbounds[varnr]);
         }
#endif
#endif

      /* update LPI state if node is pseudobranched or feasible */
      switch( reopttype ) {
         case SCIP_REOPTTYPE_TRANSIT:
            assert(reopt->reopttree->reoptnodes[id]->nconss == 0);

            if( set->reopt_shrinkinner )
            {
               SCIP_CALL( shrinkNode(reopt, set, node, id, &shrank, blkmem) );
            }

            goto TRANSIT;

            break; /*lint !e527*/

         case SCIP_REOPTTYPE_LOGICORNODE:
         case SCIP_REOPTTYPE_LEAF:
            goto TRANSIT;
            break; /*lint !e527*/

         case SCIP_REOPTTYPE_INFSUBTREE:
            /* delete the whole subtree induced be the current node */
            SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, id, FALSE, FALSE) );
            goto PSEUDO;
            break; /*lint !e527*/

         case SCIP_REOPTTYPE_STRBRANCHED:
            goto PSEUDO;
            break; /*lint !e527*/

         case SCIP_REOPTTYPE_FEASIBLE:
            /* delete the subtree */
            if( set->reopt_reducetofrontier )
            {
               SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, id, FALSE, FALSE) );
               SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );
            }
            /* dive through all children and change the reopttype to PRUNED */
            else
            {
               SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, id, SCIP_REOPTTYPE_PRUNED) );
            }
            goto FEASIBLE;
            break; /*lint !e527*/

         case SCIP_REOPTTYPE_PRUNED:
            /* delete the subtree */
            if( set->reopt_reducetofrontier )
            {
               SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, id, FALSE, FALSE) );
               SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );
            }
            /* dive through all children and change the reopttype to LEAF */
            else
            {
               SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, id, SCIP_REOPTTYPE_PRUNED) );
            }

            ++reopt->reopttree->ncutoffreoptnodes;
            ++reopt->reopttree->ntotalcutoffreoptnodes;

            goto PRUNED;
            break; /*lint !e527*/

         default:
            break;
      } /*lint !e788*/

      /* stop clock */
      SCIPclockStart(reopt->savingtime, set);

      return SCIP_OKAY;
   }

   /* get new IDs */
   SCIP_CALL( reopttreeCheckMemory(reopt->reopttree, blkmem) );

   /* the current node is the root node */
   if( isrootnode )
   {
      id = 0;

      switch( reopttype ) {
         case SCIP_REOPTTYPE_TRANSIT:
            /* ensure that no dual constraints are stored */
            SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );

            /* update the lowerbound */
            if( !SCIPsetIsEQ(set, REALABS(lowerbound), SCIPsetInfinity(set)) )
               reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

            goto TRANSIT;
            break; /*lint !e527*/

         case SCIP_REOPTTYPE_INFSUBTREE:
         case SCIP_REOPTTYPE_STRBRANCHED:
            reopt->reopttree->reoptnodes[0]->reopttype = (unsigned int)reopttype;
            reopt->reopttree->reoptnodes[0]->dualfixing = TRUE;
            reopt->reopttree->reoptnodes[0]->nvars = 0;

            if( reopttype == SCIP_REOPTTYPE_INFSUBTREE )
            {
               /* delete the whole subtree induced be the current node */
               SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, 0, FALSE, FALSE) );
            }

            SCIPdebugMessage("update node %d at ID %d:\n", 1, 0);
            SCIPdebugMessage(" -> nvars: 0, ncons: 0, parentID: -, reopttype: %u\n", reopttype);

            /* update the lowerbound */
            if( !SCIPsetIsEQ(set, REALABS(lowerbound), SCIPsetInfinity(set)) )
               reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

            goto PSEUDO;
            break; /*lint !e527*/

         case SCIP_REOPTTYPE_FEASIBLE:
            ++reopt->reopttree->ntotalfeasnodes;
            ++reopt->reopttree->nfeasnodes;
            reopt->reopttree->reoptnodes[0]->reopttype = (unsigned int)SCIP_REOPTTYPE_FEASIBLE;
            reopt->reopttree->reoptnodes[0]->dualfixing = FALSE;

            if( reopt->reopttree->reoptnodes[0]->childids != NULL && reopt->reopttree->reoptnodes[0]->nchilds > 0 )
            {
              /* delete the subtree */
               if( set->reopt_reducetofrontier )
               {
                  SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, 0, FALSE, FALSE) );
                  SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );
               }
               /* dive through all children and change the reopttype to LEAF */
               else
               {
                  SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, 0, SCIP_REOPTTYPE_PRUNED) );
               }
            }
            else
               SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );

            /* update the lowerbound */
            if( !SCIPsetIsEQ(set, REALABS(lowerbound), SCIPsetInfinity(set)) )
               reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

            SCIPdebugMessage("update node %d at ID %d:\n", 1, 0);
            SCIPdebugMessage(" -> nvars: 0, ncons: 0, parentID: -, reopttype: %u\n", reopttype);

            break;

         case SCIP_REOPTTYPE_PRUNED:
            ++reopt->reopttree->nprunednodes;
            ++reopt->reopttree->ntotalprunednodes;
            reopt->reopttree->reoptnodes[0]->reopttype = (unsigned int)SCIP_REOPTTYPE_PRUNED;
            reopt->reopttree->reoptnodes[0]->dualfixing = FALSE;

            if( reopt->reopttree->reoptnodes[0]->childids != NULL && reopt->reopttree->reoptnodes[0]->nchilds > 0 )
            {
               /* delete the subtree */
               if( set->reopt_reducetofrontier )
               {
                  SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, 0, FALSE, FALSE) );
                  SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );
               }
               /* dive through all children and change the reopttype to LEAF */
               else
               {
                  SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, 0, SCIP_REOPTTYPE_PRUNED) );
               }
            }
            else
               SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );

            /* update the lowerbound if it was not set */
            if( !SCIPsetIsEQ(set, REALABS(lowerbound), SCIPsetInfinity(set)) )
               reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

            SCIPdebugMessage("update node %d at ID %d:\n", 1, 0);
            SCIPdebugMessage(" -> nvars: 0, ncons: 0, parentID: -, reopttype: %u\n", reopttype);

            break;

         default:
            assert(reopttype == SCIP_REOPTTYPE_TRANSIT
                || reopttype == SCIP_REOPTTYPE_INFSUBTREE
                || reopttype == SCIP_REOPTTYPE_STRBRANCHED
                || reopttype == SCIP_REOPTTYPE_PRUNED
                || reopttype == SCIP_REOPTTYPE_FEASIBLE);
            break;
      }/*lint !e788*/

      /* reset the information of dual bound changes */
      reopt->currentnode = -1;
      if( reopt->dualcons != NULL )
         reopt->dualcons->nvars = 0;

      /* stop clock */
      SCIPclockStop(reopt->savingtime, set);

      return SCIP_OKAY;
   }
   else
   {
      int nbndchgdiff;
      SCIP_Bool transintoorig;

      SCIPdebugMessage("try to add node #%lld to the reopttree\n", SCIPnodeGetNumber(node));
      SCIPdebugMessage(" -> reopttype = %u\n", reopttype);

      /*
       *  check if we really want to save this node:
       *  1. save the node if reopttype is at least LOGICORNODE
       *  2. save the node if the number of bound changes of this node
       *     and the last saved node is at least a given number n
       */

      /* get the ID of the last saved node or 0 for the root */
      SCIP_CALL( getLastSavedNode(reopt, set, node, &parent, &parentid, &nbndchgdiff) );

      if( reopttype < SCIP_REOPTTYPE_INFSUBTREE && nbndchgdiff <= set->reopt_maxdiffofnodes)
      {
         SCIPdebugMessage(" -> skip saving\n");

         /* stop clock */
         SCIPclockStop(reopt->savingtime, set);

         return SCIP_OKAY;
      }

      /* check if there are free slots to store the node */
      SCIP_CALL( reopttreeCheckMemory(reopt->reopttree, blkmem) );

      id = (unsigned int) (size_t) SCIPqueueRemove(reopt->reopttree->openids);

      SCIPdebugMessage(" -> save at ID %u\n", id);

      assert(reopt->reopttree->reoptnodes[id] == NULL
         || (reopt->reopttree->reoptnodes[id]->nvars == 0 && reopt->reopttree->reoptnodes[id]->nconss == 0));
      assert(id >= 1 && id < reopt->reopttree->reoptnodessize);
      assert(!isrootnode);

      /* get memory for nodedata */
      assert(reopt->reopttree->reoptnodes[id] == NULL || reopt->reopttree->reoptnodes[id]->nvars == 0);
      SCIP_CALL( createReoptnode(reopt->reopttree, set, blkmem, id) );
      reopt->reopttree->reoptnodes[id]->parentID = parentid;

      assert(parent != NULL );
      assert((SCIPnodeGetDepth(parent) == 0 && parentid == 0) || (SCIPnodeGetDepth(parent) >= 1 && parentid > 0));
      assert(id >= 1);

      /* create the array of "child nodes" if they not exist */
      if( reopt->reopttree->reoptnodes[parentid]->childids == NULL
       || reopt->reopttree->reoptnodes[parentid]->allocchildmem == 0 )
      {
         SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[parentid], blkmem, 0, 10, 0) );
      }

      /* add the "child node" */
      SCIP_CALL( reoptAddChild(reopt->reopttree, blkmem, parentid, id) );

      /* save branching path */
      SCIP_CALL( saveAncestorBranchings(reopt->reopttree, blkmem, node, parent, id, parentid) );

      /* save bound changes after some dual reduction */
      if( saveafterdual )
      {
         SCIP_CALL( saveAfterDualBranchings(reopt, blkmem, node, id, &transintoorig) );
      }
      else
      {
         SCIPdebugMessage(" -> skip saving bound changes after dual reductions.\n");
      }

      /* transform all bounds of branched variables and ensure that they are original. */
      SCIP_CALL( transformIntoOrig(reopt, id) );

      /* save pseudo-constraints (if one exists) */
      if (SCIPnodeGetNAddedConss(node) >= 1)
      {
         assert(reopt->reopttree->reoptnodes[id]->nconss == 0);

         SCIP_CALL( saveLocalConssData(reopt->reopttree, set, blkmem, node, id) );
      }

      /* update the lowerbound if it was not set */
      if( !SCIPsetIsEQ(set, REALABS(lowerbound), SCIPsetInfinity(set)) )
         reopt->reopttree->reoptnodes[id]->lowerbound = lowerbound;

      /* set ID */
      SCIPnodeSetReoptID(node, id);

      /* set the REOPTTYPE */
      SCIPnodeSetReopttype(node, reopttype);

#ifdef SCIP_DEBUG
      int varnr;
      SCIPdebugMessage("save node #%lld successful\n", SCIPnodeGetNumber(node));
      SCIPdebugMessage(" -> ID %d, nvars %d, ncons %d, reopttype %d\n",
            id, reopt->reopttree->reoptnodes[id]->nvars + reopt->reopttree->reoptnodes[id]->nafterdualvars,
            reopt->reopttree->reoptnodes[id]->nconss,
            reopttype);
#ifdef SCIP_MORE_DEBUG
      for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nvars; varnr++)
      {
         SCIPdebugMessage("  <%s> %s %g\n", SCIPvarGetName(reopt->reopttree->reoptnodes[id]->vars[varnr]),
               reopt->reopttree->reoptnodes[id]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                     "=>" : "<=", reopt->reopttree->reoptnodes[id]->varbounds[varnr]);
      }
      for (varnr = 0; varnr < reopt->reopttree->reoptnodes[id]->nafterdualvars; varnr++)
      {
         SCIPdebugMessage("  <%s> %s %g (after dual red.)\n",
               SCIPvarGetName(reopt->reopttree->reoptnodes[id]->afterdualvars[varnr]),
               reopt->reopttree->reoptnodes[id]->afterdualvarboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                     "=>" : "<=", reopt->reopttree->reoptnodes[id]->afterdualvarbounds[varnr]);
      }
#endif
#endif
   }

   switch( reopttype ) {
      case SCIP_REOPTTYPE_TRANSIT:
      case SCIP_REOPTTYPE_LOGICORNODE:
      case SCIP_REOPTTYPE_LEAF:
         TRANSIT:

         if( !shrank )
         {
            reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)reopttype;
         }
         else
         {
            SCIPnodeSetReoptID(node, 0);
            SCIPnodeSetReopttype(node, SCIP_REOPTTYPE_NONE);
         }
         break;

      case SCIP_REOPTTYPE_INFSUBTREE:
      case SCIP_REOPTTYPE_STRBRANCHED:
         PSEUDO:

         assert(reopt->currentnode == SCIPnodeGetNumber(node));

         reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)reopttype;
         reopt->reopttree->reoptnodes[id]->dualfixing = TRUE;

         /* get all the dual information and decide if the constraint need
          * to be added next or after next */
         SCIP_CALL( collectDualInformation(reopt, blkmem, node, id, reopttype) );

         break;

      case SCIP_REOPTTYPE_FEASIBLE:
         FEASIBLE:
         reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_FEASIBLE;
         reopt->reopttree->reoptnodes[id]->dualfixing = FALSE;
         ++reopt->reopttree->nfeasnodes;
         ++reopt->reopttree->ntotalfeasnodes;

         break;

      case SCIP_REOPTTYPE_PRUNED:
         PRUNED:

         reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_PRUNED;
         reopt->reopttree->reoptnodes[id]->dualfixing = FALSE;
         ++reopt->reopttree->nprunednodes;
         ++reopt->reopttree->ntotalprunednodes;

         break;

      default:
         assert(reopttype == SCIP_REOPTTYPE_TRANSIT
             || reopttype == SCIP_REOPTTYPE_LOGICORNODE
             || reopttype == SCIP_REOPTTYPE_LEAF
             || reopttype == SCIP_REOPTTYPE_INFSUBTREE
             || reopttype == SCIP_REOPTTYPE_STRBRANCHED
             || reopttype == SCIP_REOPTTYPE_FEASIBLE
             || reopttype == SCIP_REOPTTYPE_PRUNED);
         break;
   } /*lint !e788*/

   /* stop clock */
   SCIPclockStop(reopt->savingtime, set);

   /* reset the information of dual bound changes */
   reopt->currentnode = -1;
   if( reopt->dualcons != NULL )
      reopt->dualcons->nvars = 0;

   return SCIP_OKAY;
}

/** delete the stored information about dual bound changes of the last focused node */
static
void deleteLastDualBndchgs(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   if( reopt->dualcons != NULL && reopt->dualcons->nvars > 0 )
   {
      SCIPdebugMessage("delete %d dual variable information about node %lld\n", reopt->dualcons->nvars,
            reopt->currentnode);
      reopt->dualcons->nvars = 0;
      reopt->currentnode = -1;
   }
}

/** delete the stored constraints that dual information at the given reoptimization node */
static
SCIP_RETCODE reoptnodeResetDualConss(
   SCIP_REOPTNODE*       reoptnode,          /**< reoptimization node */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reoptnode != NULL);
   assert(blkmem != NULL);

   if( reoptnode->dualconscur != NULL )
   {
      SCIPdebugMessage("reset dual (1) information\n");

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->vals, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->vars, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualconscur);
      reoptnode->dualconscur = NULL;
   }

   if( reoptnode->dualconsnex != NULL )
   {
      SCIPdebugMessage("reset dual (2) information\n");

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconsnex->vals, reoptnode->dualconsnex->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconsnex->vars, reoptnode->dualconsnex->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualconsnex);
      reoptnode->dualconsnex = NULL;
   }

   reoptnode->dualfixing = FALSE;

   return SCIP_OKAY;
}

/** generate a global constraint to separate an infeasible subtree */
static
SCIP_RETCODE saveGlobalCons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   REOPT_CONSTYPE        consttype           /**< reopttype of the constraint */
   )
{
   assert(reopt != NULL);
   assert(node != NULL);

   if( consttype == REOPT_CONSTYPE_INFSUBTREE )
   {
      SCIP_BOUNDTYPE* boundtypes;
      int nbranchvars;
      int nvars;
      int nglbconss;
      int v;

      nglbconss = reopt->nglbconss;
      nvars = SCIPnodeGetDepth(node)+1;

      /* check if enough memory to store the global constraint is available */
      SCIP_CALL( checkMemGlbCons(reopt, blkmem, nglbconss+1) );

      /* allocate memory to store the infeasible path
       * we use the permanent allocated array consbounds to store the boundtypes */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopt->glbconss[nglbconss]) ); /*lint !e866*/
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vars, nvars) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vals, nvars) );

      /* allocate buffer */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &boundtypes, nvars) );

      reopt->glbconss[nglbconss]->varssize = nvars;
      reopt->glbconss[nglbconss]->constype = REOPT_CONSTYPE_INFSUBTREE;

      SCIPnodeGetAncestorBranchings(node,
            reopt->glbconss[nglbconss]->vars,
            reopt->glbconss[nglbconss]->vals,
            boundtypes,
            &nbranchvars,
            nvars);

      if( nvars < nbranchvars )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vars, nvars, nbranchvars) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vals, nvars, nbranchvars) );
         nvars = nbranchvars;
         reopt->glbconss[nglbconss]->varssize = nvars;

         SCIPnodeGetAncestorBranchings(node, reopt->glbconss[nglbconss]->vars, reopt->glbconss[nglbconss]->vals,
               boundtypes, &nbranchvars, nvars);
      }

      /* transform into original variables */
      for(v = 0; v < nbranchvars; v++)
      {
         SCIP_Real constant;
         SCIP_Real scalar;

         constant = 0;
         scalar = 1;

         SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->glbconss[nglbconss]->vars[v], &scalar, &constant) );
         reopt->glbconss[nglbconss]->vals[v] = (reopt->glbconss[nglbconss]->vals[v] - constant)/scalar;

         assert(SCIPsetIsFeasEQ(set, reopt->glbconss[nglbconss]->vals[v], 0.0) || SCIPsetIsFeasEQ(set,
               reopt->glbconss[nglbconss]->vals[v], 1.0));
      }

      /* free buffer */
      SCIPsetFreeBufferArray(set, &boundtypes);

      /* increase the number of global constraints */
      ++reopt->nglbconss;
   }

   return SCIP_OKAY;
}


/** move all id of child nodes from reoptimization node stored at @p id1 to the node stored at @p id2 */
static
SCIP_RETCODE reoptMoveIDs(
   SCIP_REOPTTREE*       reopttree,          /**< reopttree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          id1,                /**< source id */
   unsigned int          id2                 /**< target id */
   )
{
   int c;
   int nchilds_id1;
   int nchilds_id2;

   assert(reopttree != NULL);
   assert(blkmem != NULL);
   assert(id1 < reopttree->reoptnodessize);
   assert(id2 < reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id1] != NULL);
   assert(reopttree->reoptnodes[id2] != NULL);

   nchilds_id1 = reopttree->reoptnodes[id1]->nchilds;
   nchilds_id2 = reopttree->reoptnodes[id2]->nchilds;

   /* ensure that the array storing the child id's is large enough */
   SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id2], blkmem, 0, nchilds_id1+nchilds_id2, 0) );
   assert(reopttree->reoptnodes[id2]->allocchildmem >= nchilds_id1+nchilds_id2);

   SCIPdebugMessage("move %d IDs: %u -> %u\n", nchilds_id1, id1, id2);

   /* move the ids */
   for(c = 0; c < nchilds_id1; c++)
   {

#ifdef SCIP_DEBUG
      /* check that no id is added twice */
      int k;
      for(k = 0; k < nchilds_id2; k++)
         assert(reopttree->reoptnodes[id2]->childids[k] != reopttree->reoptnodes[id1]->childids[c]);
#endif

      reopttree->reoptnodes[id2]->childids[nchilds_id2+c] = reopttree->reoptnodes[id1]->childids[c];
   }

   /* update the number of childs */
   reopttree->reoptnodes[id1]->nchilds = 0;
   reopttree->reoptnodes[id2]->nchilds += nchilds_id1;

   return SCIP_OKAY;
}

/** change all bound changes along the root path */
static
SCIP_RETCODE changeAncestorBranchings(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< search tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the branch and bound tree */
   unsigned int          id,                 /**< id of stored node */
   SCIP_Bool             afterdualintobranching /**< apply and convert bound changes made after the first based on dual information into branchings */
   )
{
   SCIP_REOPTTREE* reopttree;
   SCIP_REOPTNODE* reoptnode;
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(cliquetable != NULL);
   assert(node != NULL);
   assert(blkmem != NULL);

   reopttree = reopt->reopttree;
   assert(reopttree != NULL);
   assert(id < reopttree->reoptnodessize);

   reoptnode = reopttree->reoptnodes[id];
   assert(reoptnode != NULL);

   /* copy memory to ensure that only original variables are saved */
   if( reoptnode->nvars == 0 && reoptnode->nafterdualvars == 0)
      return SCIP_OKAY;

   /* change the bounds along the branching path */
   for(v = 0; v < reoptnode->nvars; v++)
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_BOUNDTYPE boundtype;
      SCIP_Real oldlb;
      SCIP_Real oldub;
      SCIP_Real newbound;

      var = reoptnode->vars[v];
      val = reoptnode->varbounds[v];
      boundtype = reoptnode->varboundtypes[v];

      assert(SCIPvarIsOriginal(var));
      SCIP_CALL( SCIPvarGetProbvarBound(&var, &val, &boundtype) );
      assert(SCIPvarIsTransformed(var));
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

      oldlb = SCIPvarGetLbLocal(var);
      oldub = SCIPvarGetUbLocal(var);
      newbound = val;

      assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

      if(boundtype == SCIP_BOUNDTYPE_LOWER
      && SCIPsetIsGT(set, newbound, oldlb)
      && SCIPsetIsFeasLE(set, newbound, oldub))
      {
         SCIPvarAdjustLb(var, set, &newbound);

         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      }
      else if(boundtype == SCIP_BOUNDTYPE_UPPER
           && SCIPsetIsLT(set, newbound, oldub)
           && SCIPsetIsFeasGE(set, newbound, oldlb))
      {
         SCIPvarAdjustUb(var, set, &newbound);

         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      }
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("  (path) <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=",
            newbound);
#endif
   }

   if( afterdualintobranching && reoptnode->nafterdualvars > 0 )
   {
      /* check the memory to convert this bound changes into 'normal' */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], blkmem, reoptnode->nvars + reoptnode->nafterdualvars,
            0, 0) );

      /* change the bounds */
      for(v = 0; v < reoptnode->nafterdualvars; v++)
      {
         SCIP_VAR* var;
         SCIP_Real val;
         SCIP_BOUNDTYPE boundtype;
         SCIP_Bool bndchgd;
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Real newbound;

         var = reoptnode->afterdualvars[v];
         val = reoptnode->afterdualvarbounds[v];
         boundtype = reoptnode->afterdualvarboundtypes[v];

         assert(SCIPvarIsOriginal(var));
         SCIP_CALL( SCIPvarGetProbvarBound(&var, &val, &boundtype) );
         assert(SCIPvarIsTransformed(var));
         assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

         bndchgd = FALSE;

         oldlb = SCIPvarGetLbLocal(var);
         oldub = SCIPvarGetUbLocal(var);
         newbound = val;

         if(boundtype == SCIP_BOUNDTYPE_LOWER
         && SCIPsetIsGT(set, newbound, oldlb)
         && SCIPsetIsFeasLE(set, newbound, oldub))
         {
            SCIPvarAdjustLb(var, set, &newbound);
            SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );

            bndchgd = TRUE;
         }
         else if(boundtype == SCIP_BOUNDTYPE_UPPER
              && SCIPsetIsLT(set, newbound, oldub)
              && SCIPsetIsFeasGE(set, newbound, oldlb))
         {
            SCIPvarAdjustUb(var, set, &newbound);
            SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
                  tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );

            bndchgd = TRUE;
         }

         assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("   (prop) <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=",
               newbound);
#endif
         if( bndchgd )
         {
            int nvars;

            nvars = reoptnode->nvars;
            reoptnode->vars[nvars] = reoptnode->afterdualvars[v];
            reoptnode->varbounds[nvars] = reoptnode->afterdualvarbounds[v];
            reoptnode->varboundtypes[nvars] = reoptnode->afterdualvarboundtypes[v];
            ++reoptnode->nvars;
         }
      }

      /* free the afterdualvars, -bounds, and -boundtypes */
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->afterdualvarboundtypes, reoptnode->afterdualvarssize);
      reoptnode->afterdualvarboundtypes = NULL;

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->afterdualvarbounds, reoptnode->afterdualvarssize);
      reoptnode->afterdualvarbounds = NULL;

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->afterdualvars, reoptnode->afterdualvarssize);
      reoptnode->afterdualvars = NULL;

      reoptnode->nafterdualvars = 0;
      reoptnode->afterdualvarssize = 0;
   }

   return SCIP_OKAY;
}

/** add a constraint to ensure that at least one variable bound gets different */
static
SCIP_RETCODE addSplitcons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< search tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_NODE*            node,               /**< node corresponding to the pruned part */
   unsigned int          id                  /**< id of stored node */
   )
{
   SCIP_CONS* cons;
   const char* name;
   int v;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);
   assert(reopt->reopttree->reoptnodes[id]->dualfixing);
   assert(reopt->reopttree->reoptnodes[id]->dualconscur != NULL);
   assert(scip != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(node != NULL);

   if( reopt->reopttree->reoptnodes[id]->dualconscur->constype == REOPT_CONSTYPE_STRBRANCHED )
   {
      SCIPdebugMessage(" create a split-node #%lld\n", SCIPnodeGetNumber(node));
   }
   else if( reopt->reopttree->reoptnodes[id]->dualconscur->constype == REOPT_CONSTYPE_INFSUBTREE )
   {
      SCIPdebugMessage(" separate an infeasible subtree\n");
   }

   /* if the constraint consists of exactly one variable it can be interpreted
    * as a normal branching step, i.e., we can fix the variable to the negated bound */
   if( reopt->reopttree->reoptnodes[id]->dualconscur->nvars == 1 )
   {
      SCIP_VAR* var;
      SCIP_BOUNDTYPE boundtype;
      SCIP_Real oldlb;
      SCIP_Real oldub;
      SCIP_Real newbound;

      var = reopt->reopttree->reoptnodes[id]->dualconscur->vars[0];
      newbound = reopt->reopttree->reoptnodes[id]->dualconscur->vals[0];
      boundtype = SCIPsetIsFeasEQ(set, newbound, 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;

      assert(SCIPvarIsOriginal(var));
      SCIP_CALL( SCIPvarGetProbvarBound(&var, &newbound, &boundtype) );
      assert(SCIPvarIsTransformed(var));

      oldlb = SCIPvarGetLbLocal(var);
      oldub = SCIPvarGetUbLocal(var);

      /* negate the bound */
      newbound = 1 - newbound;
      boundtype = (SCIP_BOUNDTYPE) (1 - (int)boundtype);

      if(boundtype == SCIP_BOUNDTYPE_LOWER
      && SCIPsetIsGT(set, newbound, oldlb)
      && SCIPsetIsFeasLE(set, newbound, oldub))
      {
         SCIPvarAdjustLb(var, set, &newbound);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      }
      else if(boundtype == SCIP_BOUNDTYPE_UPPER
           && SCIPsetIsLT(set, newbound, oldub)
           && SCIPsetIsFeasGE(set, newbound, oldlb))
      {
         SCIPvarAdjustUb(var, set, &newbound);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      }

      assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

      SCIPdebugMessage("  -> constraint consists of only one variable: <%s> %s %g\n", SCIPvarGetName(var),
            boundtype == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);
   }
   else
   {
      SCIP_VAR** consvars;
      SCIP_Real consval;
      SCIP_BOUNDTYPE consboundtype;

      /* allocate buffer */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, reopt->reopttree->reoptnodes[id]->dualconscur->nvars) );

      for(v = 0; v < reopt->reopttree->reoptnodes[id]->dualconscur->nvars; v++)
      {
         consvars[v] = reopt->reopttree->reoptnodes[id]->dualconscur->vars[v];
         consval = reopt->reopttree->reoptnodes[id]->dualconscur->vals[v];
         consboundtype = SCIPsetIsFeasEQ(set, consval, 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;

        assert(SCIPvarIsOriginal(consvars[v]));
        SCIP_CALL( SCIPvarGetProbvarBound(&consvars[v], &consval, &consboundtype) );
        assert(SCIPvarIsTransformed(consvars[v]));
        assert(SCIPvarGetStatus(consvars[v]) != SCIP_VARSTATUS_MULTAGGR);

        if ( SCIPsetIsFeasEQ(set, consval, 1.0) )
        {
           SCIP_CALL( SCIPvarNegate(consvars[v], blkmem, set, stat, &consvars[v]) );
           assert(SCIPvarIsNegated(consvars[v]));
        }
      }

      if( reopt->reopttree->reoptnodes[id]->dualconscur->constype == REOPT_CONSTYPE_INFSUBTREE )
         name = "infsubtree";
      else
      {
         assert(reopt->reopttree->reoptnodes[id]->dualconscur->constype == REOPT_CONSTYPE_STRBRANCHED);
         name = "splitcons";
      }

      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, reopt->reopttree->reoptnodes[id]->dualconscur->nvars, consvars,
            FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

      SCIPdebugMessage(" -> add constraint in node #%lld:\n", SCIPnodeGetNumber(node));
      SCIPdebugPrintCons(scip, cons, NULL);

      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* free buffer */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** fix all bounds ad stored in dualconscur at the given node @p node_fix */
static
SCIP_RETCODE fixBounds(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< search tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node corresponding to the fixed part */
   unsigned int          id,                 /**< id of stored node */
   SCIP_Bool             updatedualconss     /**< update constraint representing dual bound changes */
   )
{
   SCIP_REOPTTREE* reopttree;
   SCIP_REOPTNODE* reoptnode;
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(cliquetable != NULL);
   assert(node != NULL);
   assert(blkmem != NULL);

   reopttree = reopt->reopttree;
   assert(0 < id && id < reopttree->reoptnodessize);
   assert(reopttree != NULL);

   reoptnode = reopttree->reoptnodes[id];
   assert(reoptnode != NULL);
   assert(reoptnode->dualfixing);
   assert(reoptnode->dualconscur != NULL);

   /* ensure that the arrays to store the bound changes are large enough */
   SCIP_CALL( reoptnodeCheckMemory(reoptnode, blkmem, reoptnode->nvars + reoptnode->dualconscur->nvars, 0, 0) );

   for(v = 0; v < reoptnode->dualconscur->nvars; v++)
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_BOUNDTYPE boundtype;
      SCIP_Bool bndchgd;

      var = reoptnode->dualconscur->vars[v];
      val = reoptnode->dualconscur->vals[v];
      boundtype = SCIPsetIsFeasEQ(set, val, 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;

      SCIP_CALL(SCIPvarGetProbvarBound(&var, &val, &boundtype));
      assert(SCIPvarIsTransformedOrigvar(var));

      bndchgd = FALSE;

      if(boundtype == SCIP_BOUNDTYPE_LOWER
      && SCIPsetIsGT(set, val, SCIPvarGetLbLocal(var))
      && SCIPsetIsFeasLE(set, val, SCIPvarGetUbLocal(var)))
      {
         SCIPvarAdjustLb(var, set, &val);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, val, SCIP_BOUNDTYPE_LOWER, FALSE) );

         bndchgd = TRUE;
      }
      else if(boundtype == SCIP_BOUNDTYPE_UPPER
           && SCIPsetIsLT(set, val, SCIPvarGetUbLocal(var))
           && SCIPsetIsFeasGE(set, val, SCIPvarGetLbLocal(var)))
      {
         SCIPvarAdjustUb(var, set, &val);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, val, SCIP_BOUNDTYPE_UPPER, FALSE) );

         bndchgd = TRUE;
      }
      else if(boundtype != SCIP_BOUNDTYPE_LOWER && boundtype != SCIP_BOUNDTYPE_UPPER)
      {
         printf("** Unknown boundtype: %d **\n", boundtype);
         assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);
      }
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("  (dual) <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", val);
#endif
      /* add variable and bound to branching path information, because we don't want to delete this data */
      if( bndchgd )
      {
         int pos;
         SCIP_Real constant;
         SCIP_Real scalar;

         pos = reoptnode->nvars;

         reoptnode->vars[pos] = var;
         scalar = 1.0;
         constant = 0.0;
         SCIP_CALL( SCIPvarGetOrigvarSum(&reoptnode->vars[pos], &scalar, &constant) );
         assert(SCIPvarIsOriginal(reoptnode->vars[pos]));

         reoptnode->varbounds[pos] = reoptnode->dualconscur->vals[v];
         reoptnode->varboundtypes[pos] = (SCIPsetIsFeasEQ(set, reoptnode->varbounds[pos], 0.0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
         ++reoptnode->nvars;
      }
   }

   if( updatedualconss )
   {
      /* delete dualconscur and move dualconsnex -> dualconscur */
      SCIP_CALL( reoptnodeUpdateDualConss(reoptnode, blkmem) );
   }

   return SCIP_OKAY;
}

/** fix all bounds corresponding to dual bound changes in a previous iteration in the fashion of interdiction branching;
 *  keep the first negbndchg-1 bound changes as stored in dualconscur and negate the negbndchg-th bound.
 */
static
SCIP_RETCODE fixInterdiction(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< search tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< child node */
   unsigned int          id,                 /**< id of the node */
   SCIP_VAR**            vars,               /**< variables in permuted order */
   SCIP_Real*            vals,               /**< bounds in permuted order */
   int                   nvars,              /**< number of variables */
   int                   negbndchg           /**< index of the variable that should negated */
   )
{
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_BOUNDTYPE boundtype;
   int nbndchgs;
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(cliquetable != NULL);
   assert(node != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(nvars >= 0);
   assert(blkmem != NULL);
   assert(0 < id && id < reopt->reopttree->reoptnodessize);

#ifndef NDEBUG
   {
      SCIP_REOPTTREE* reopttree;
      SCIP_REOPTNODE* reoptnode;

      reopttree = reopt->reopttree;
      assert(reopttree != NULL);

      reoptnode = reopttree->reoptnodes[id];
      assert(reoptnode != NULL);
      assert(reoptnode->dualfixing);
   }
#endif

   nbndchgs = MIN(negbndchg, nvars);

   /* change the first negbndchg-1 bounds as stored in dualconscur and negate the negbndchg-th bound */
   for(v = 0; v < nbndchgs; v++)
   {
      var = vars[v];
      val = vals[v];
      boundtype = SCIPsetIsFeasEQ(set, val, 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;

      SCIP_CALL(SCIPvarGetProbvarBound(&var, &val, &boundtype));
      assert(SCIPvarIsTransformedOrigvar(var));

      /* negate the negbndchg-th bound */
      if( v == nbndchgs-1 )
      {
         val = 1-val;
         boundtype = (SCIP_BOUNDTYPE)(SCIP_BOUNDTYPE_UPPER - boundtype); /*lint !e656*/
      }

      if(boundtype == SCIP_BOUNDTYPE_LOWER
      && SCIPsetIsGT(set, val, SCIPvarGetLbLocal(var))
      && SCIPsetIsFeasLE(set, val, SCIPvarGetUbLocal(var)))
      {
         SCIPvarAdjustLb(var, set, &val);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, val, SCIP_BOUNDTYPE_LOWER, FALSE) );
      }
      else if(boundtype == SCIP_BOUNDTYPE_UPPER
           && SCIPsetIsLT(set, val, SCIPvarGetUbLocal(var))
           && SCIPsetIsFeasGE(set, val, SCIPvarGetLbLocal(var)))
      {
         SCIPvarAdjustUb(var, set, &val);
         SCIP_CALL( SCIPnodeAddBoundchg(node, blkmem, set, stat, transprob, origprob,
               tree, reopt, lp, branchcand, eventqueue, cliquetable, var, val, SCIP_BOUNDTYPE_UPPER, FALSE) );
      }
      else if(boundtype != SCIP_BOUNDTYPE_LOWER && boundtype != SCIP_BOUNDTYPE_UPPER)
      {
         printf("** Unknown boundtype: %d **\n", boundtype);
         assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);
      }
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("  (dual) <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", val);
#endif
   }

   return SCIP_OKAY;
}

/** add all constraints stored at @p id to the given nodes @p node_fix and @p node_cons */
static
SCIP_RETCODE addLocalConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the branch and bound tree*/
   unsigned int          id                  /**< id of stored node */
   )
{
   int c;
   const char* name;

   assert(scip != NULL);
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(0 < id && id < reopt->reopttree->reoptnodessize);

   if( reopt->reopttree->reoptnodes[id]->nconss == 0 )
      return SCIP_OKAY;

   SCIPdebugMessage(" -> add %d constraint(s) to node #%lld:\n", reopt->reopttree->reoptnodes[id]->nconss, SCIPnodeGetNumber(node));

   for( c = 0; c < reopt->reopttree->reoptnodes[id]->nconss; c++ )
   {
      SCIP_CONS* cons;
      LOGICORDATA* consdata;
      SCIP_VAR** consvars;
      SCIP_Real consval;
      SCIP_BOUNDTYPE consboundtype;
      int v;

      consdata = reopt->reopttree->reoptnodes[id]->conss[c];
      assert(consdata != NULL);
      assert(consdata->nvars > 0);
      assert(consdata->varssize >= consdata->nvars);

      /* allocate buffer */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, consdata->nvars) );

      /* iterate over all variable and transform them */
      for(v = 0; v < consdata->nvars; v++)
      {
         consvars[v] = consdata->vars[v];
         consval= consdata->vals[v];
         consboundtype = SCIPsetIsFeasEQ(set, consval, 0.0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER;

         assert(SCIPvarIsOriginal(consvars[v]));
         SCIP_CALL( SCIPvarGetProbvarBound(&consvars[v], &consval, &consboundtype) );
         assert(SCIPvarIsTransformed(consvars[v]));
         assert(SCIPvarGetStatus(consvars[v]) != SCIP_VARSTATUS_MULTAGGR);

         if( SCIPsetIsFeasEQ(set, consval, 1.0) )
         {
            SCIP_CALL( SCIPvarNegate(consvars[v], blkmem, set, stat, &consvars[v]) );
            assert(SCIPvarIsNegated(consvars[v]));
         }
      }

      assert(consdata->constype == REOPT_CONSTYPE_INFSUBTREE || consdata->constype == REOPT_CONSTYPE_STRBRANCHED);

      if( consdata->constype == REOPT_CONSTYPE_INFSUBTREE )
         name = "infsubtree";
      else
         name = "splitcons";


      /* create the constraints and add them to the corresponding nodes */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, consdata->nvars, consvars,
            FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

      SCIPdebugPrintCons(scip, cons, NULL);

      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* free buffer */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** reset the internal statistics at the beginning of a new iteration */
static
void resetStats(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   reopt->lastbranched = -1;
   reopt->currentnode = -1;
   reopt->lastseennode = -1;
   reopt->reopttree->nfeasnodes = 0;
   reopt->reopttree->ninfnodes = 0;
   reopt->reopttree->nprunednodes = 0;
   reopt->reopttree->ncutoffreoptnodes = 0;

   return;
}

/** check the stored bound changes of all child nodes for redundancy and infeasibility.
 *
 *  due to strongbranching initialization at node stored at @p id it can happen, that some bound changes stored in the
 *  child nodes of the reoptimization node stored at @p id become redundant or make the subproblem infeasible. in this
 *  method we remove all redundant bound changes and delete infeasible child nodes.
 */
static
SCIP_RETCODE dryBranch(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool*            runagain,           /**< pointer to store of this method should run again */
   unsigned int          id                  /**< id of stored node */
   )
{
   SCIP_REOPTNODE* reoptnode;
   unsigned int* cutoffchilds;
   int ncutoffchilds;
   unsigned int* redchilds;
   int nredchilds;
   int c;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes != NULL);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   reoptnode = reopt->reopttree->reoptnodes[id];

   *runagain = FALSE;
   ncutoffchilds = 0;
   nredchilds = 0;

   SCIPdebugMessage("start dry branching of node at ID %u\n", id);

   /* allocate buffer arrays */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cutoffchilds, reoptnode->nchilds) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &redchilds, reoptnode->nchilds) );

   /* iterate over all child nodes and check each bound changes
    * for redundancy and conflict */
   for( c = 0; c < reoptnode->nchilds; c++ )
   {
      SCIP_REOPTNODE* child;
      SCIP_Bool cutoff;
      SCIP_Bool redundant;
      int* redundantvars;
      int nredundantvars;
      int v;
      unsigned int childid;

      cutoff = FALSE;
      redundant = FALSE;
      nredundantvars = 0;

      childid = reoptnode->childids[c];
      assert(childid < reopt->reopttree->reoptnodessize);
      child = reopt->reopttree->reoptnodes[childid];
      assert(child != NULL);
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("-> check child at ID %d (%d vars, %d conss):\n", childid, child->nvars, child->nconss);
#endif
      if( child->nvars > 0 )
      {
         /* allocate buffer memory to store the redundant variables */
         SCIP_CALL( SCIPsetAllocBufferArray(set, &redundantvars, child->nvars) );

         for(v = 0; v < child->nvars && !cutoff; v++)
         {
            SCIP_VAR* transvar;
            SCIP_Real transval;
            SCIP_BOUNDTYPE transbndtype;
            SCIP_Real ub;
            SCIP_Real lb;

            transvar = child->vars[v];
            transval = child->varbounds[v];
            transbndtype = child->varboundtypes[v];

            /* transform into the transformed space */
            SCIP_CALL( SCIPvarGetProbvarBound(&transvar, &transval, &transbndtype) );

            lb = SCIPvarGetLbLocal(transvar);
            ub = SCIPvarGetUbLocal(transvar);

            /* check for infeasibility */
            if( SCIPsetIsFeasEQ(set, lb, ub) && !SCIPsetIsFeasEQ(set, lb, transval) )
            {
               SCIPdebugMessage(" -> <%s> is fixed to %g, can not change bound to %g -> cutoff\n",
                  SCIPvarGetName(transvar), lb, transval);

               cutoff = TRUE;
               break;
            }

            /* check for redundancy */
            if( SCIPsetIsFeasEQ(set, lb, ub) && SCIPsetIsFeasEQ(set, lb, transval) )
            {
               SCIPdebugMessage(" -> <%s> is already fixed to %g -> redundant bound change\n",
                  SCIPvarGetName(transvar), lb);

               redundantvars[nredundantvars] = v;
               ++nredundantvars;
            }
         }

         if( !cutoff && nredundantvars > 0 )
         {
            for(v = 0; v < nredundantvars; v++)
            {
               /* replace the redundant variable by the last stored variable */
               child->vars[redundantvars[v]] = child->vars[child->nvars-1];
               child->varbounds[redundantvars[v]] = child->varbounds[child->nvars-1];
               child->varboundtypes[redundantvars[v]] = child->varboundtypes[child->nvars-1];
               --child->nvars;
            }
         }

         /* free buffer memory */
         SCIPsetFreeBufferArray(set, &redundantvars);
      }
      else if( child->nconss == 0 )
      {
         redundant = TRUE;
         SCIPdebugMessage(" -> redundant node found.\n");
      }

      if( cutoff )
      {
         cutoffchilds[ncutoffchilds] = childid;
         ++ncutoffchilds;
      }
      else if( redundant )
      {
         redchilds[nredchilds] = childid;
         ++nredchilds;
      }
   }

   SCIPdebugMessage("-> found %d redundant and %d infeasible nodes\n", nredchilds, ncutoffchilds);

   c = 0;

   /* delete all nodes that can be cut off */
   while( ncutoffchilds > 0 )
   {
      /* delete the node and the induced subtree */
      SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, cutoffchilds[ncutoffchilds-1], TRUE, TRUE) );

      /* find the position in the childid array */
      c = 0;
      while( reoptnode->childids[c] != cutoffchilds[ncutoffchilds-1] && c < reoptnode->nchilds )
         ++c;
      assert(reoptnode->childids[c] == cutoffchilds[ncutoffchilds-1]);

      /* replace the ID at position c by the last ID */
      reoptnode->childids[c] = reoptnode->childids[reoptnode->nchilds-1];
      --reoptnode->nchilds;

      /* decrease the number of nodes to cutoff */
      --ncutoffchilds;
   }

   c = 0;

   /* replace all redundant nodes their child nodes or cutoff the node if it is a leaf */
   while( nredchilds > 0 )
   {
      /* find the position in the childid array */
      c = 0;
      while( reoptnode->childids[c] != redchilds[nredchilds-1] && c < reoptnode->nchilds )
         ++c;
      assert(reoptnode->childids[c] == redchilds[nredchilds-1]);

      /* the node is a leaf and we can cutoff them  */
      if( reopt->reopttree->reoptnodes[redchilds[nredchilds-1]]->nchilds == 0 )
      {
         /* delete the node and the induced subtree */
         SCIP_CALL( deleteChildrenBelow(reopt->reopttree, set, blkmem, redchilds[nredchilds-1], TRUE, TRUE) );

         /* replace the ID at position c by the last ID */
         reoptnode->childids[c] = reoptnode->childids[reoptnode->nchilds-1];
         --reoptnode->nchilds;

         /* decrease the number of redundant nodes */
         --nredchilds;
      }
      else
      {
         int cc;
         int ncc;

         /* replace the ID at position c by the last ID */
         reoptnode->childids[c] = reoptnode->childids[reoptnode->nchilds-1];
         --reoptnode->nchilds;

         ncc = reopt->reopttree->reoptnodes[redchilds[nredchilds-1]]->nchilds;

         /* check the memory */
         SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[id], blkmem, 0, reoptnode->nchilds+ncc, 0) );

         /* add all IDs of child nodes to the current node */
         for(cc = 0; cc < ncc; cc++)
         {
            reoptnode->childids[reoptnode->nchilds] = reopt->reopttree->reoptnodes[redchilds[nredchilds-1]]->childids[cc];
            ++reoptnode->nchilds;
         }

         /* delete the redundant node */
         SCIP_CALL( reopttreeDeleteNode(reopt->reopttree, set, blkmem, redchilds[nredchilds-1], TRUE) );

         /* decrease the number of redundant nodes */
         --nredchilds;

         /* update the flag to rerun this method */
         *runagain = TRUE;
      }
   }

   /* free buffer arrays */
   SCIPsetFreeBufferArray(set, &cutoffchilds);
   SCIPsetFreeBufferArray(set, &redchilds);

   return SCIP_OKAY;
}

/** return the number of all nodes in the subtree induced by the reoptimization node stored at @p id */
static
int reopttreeGetNNodes(
   SCIP_REOPTTREE*       reopttree,          /**< reopttree */
   unsigned int          id                  /**< id of stored node */
   )
{
   int nnodes;
   int i;

   assert(reopttree != NULL);
   assert(id < reopttree->reoptnodessize);

   nnodes = 0;

   for(i = 0; i < reopttree->reoptnodes[id]->nchilds; i++)
   {
      nnodes += reopttreeGetNNodes(reopttree, reopttree->reoptnodes[id]->childids[i]);
   }

   return nnodes + 1;
}

/** returns the number of leaf nodes of the induced subtree */
static
int reoptGetNLeaves(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   unsigned int          id                  /**< id of stored node */
   )
{
   int i;
   int nleaves;

   assert(reopt != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   nleaves = 0;

   /* iterate over all child nods and check whether they are leaves or not */
   for(i = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++)
   {
      unsigned int childid;

      childid = reopt->reopttree->reoptnodes[id]->childids[i];
      assert(childid < reopt->reopttree->reoptnodessize);

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
         ++nleaves;
      else
         nleaves += reoptGetNLeaves(reopt, childid);
   }

   return nleaves;
}

/** returns all leaves of the subtree induced by the node stored at @p id*/
static
SCIP_RETCODE reoptGetLeaves(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure*/
   unsigned int          id,                 /**< id of stored node */
   unsigned int*         leaves,             /**< array of leave nodes */
   int                   leavessize,         /**< size of leaves array */
   int*                  nleaves             /**< pointer to store the number of leave nodes */
   )
{
   int i;
   int l;

   assert(reopt != NULL);
   assert(leavessize > 0 && leaves != NULL);
   assert((*nleaves) >= 0);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   for(i = 0, l = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++)
   {
      unsigned int childid;

      assert(*nleaves <= leavessize);

      childid = reopt->reopttree->reoptnodes[id]->childids[i];
      assert(childid < reopt->reopttree->reoptnodessize);

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
      {
         leaves[l] = reopt->reopttree->reoptnodes[id]->childids[i];
         ++l;
         ++(*nleaves);
      }
      else
      {
         int nleaves2;

         nleaves2 = 0;
         SCIP_CALL( reoptGetLeaves(reopt, childid, &leaves[l], leavessize - l, &nleaves2) );
         l += nleaves2;
         (*nleaves) += nleaves2;
      }
   }

   return SCIP_OKAY;
}

/** after restarting the reoptimization and an after compressing the search tree we have to delete all stored information */
static
SCIP_RETCODE reoptResetTree(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             softreset           /**< mark the nodes to overwriteable (TRUE) or delete them completely (FALSE) */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* clear the tree */
   SCIP_CALL( clearReoptnodes(reopt->reopttree, set, blkmem, softreset) );
   assert(reopt->reopttree->nreoptnodes == 0);

   /* reset the dual constraint */
   if( reopt->dualcons != NULL )
      reopt->dualcons->nvars = 0;

   reopt->currentnode = -1;

   return SCIP_OKAY;
}

/** restart the reoptimization by removing all stored information about nodes and increase the number of restarts */
static
SCIP_RETCODE reoptRestart(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* clear the tree */
   SCIP_CALL( reoptResetTree(reopt, set, blkmem, FALSE) );
   assert(reopt->reopttree->nreoptnodes == 0);

   /* allocate memory for the root node */
   SCIP_CALL( createReoptnode(reopt->reopttree, set, blkmem, 0) );

   reopt->nglbrestarts += 1;

   if( reopt->firstrestart == -1 )
      reopt->firstrestart = reopt->run;

   reopt->lastrestart = reopt->run;

   return SCIP_OKAY;
}

/** save the new objective function */
static
SCIP_RETCODE reoptSaveNewObj(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            transvars,          /**< transformed problem variables */
   int                   ntransvars          /**< number of transformed problem variables */
   )
{
   SCIP_Real norm;
   int v;
   int idx;

   assert(reopt != NULL);

   /* check memory */
   SCIP_CALL( ensureRunSize(reopt, reopt->run, blkmem) );

   norm = 0;

   /* get memory */
   SCIP_ALLOC( BMSallocMemoryArray(&reopt->objs[reopt->run-1], ntransvars) ); /*lint !e866*/

   /* save coefficients */
   for( v = 0; v < ntransvars; v++ )
   {
      SCIP_Real glblb;
      SCIP_Real glbub;

      assert(SCIPvarIsActive(transvars[v]));
      assert(!SCIPvarIsOriginal(transvars[v]));

      idx = SCIPvarGetProbindex(transvars[v]);
      assert(idx < ntransvars);
      assert(0 <= idx);

      reopt->objs[reopt->run-1][idx] = SCIPvarGetObj(transvars[v]);

      /* we skip global fixed variables */
      glblb = SCIPvarGetLbGlobal(transvars[v]);
      glbub = SCIPvarGetUbGlobal(transvars[v]);

      if( SCIPsetIsFeasLT(set, glblb, glbub) )
         norm += SQR(reopt->objs[reopt->run-1][idx]);

      /* mark this objective as the first non empty */
      if( reopt->firstobj == -1 && reopt->objs[reopt->run-1][idx] != 0 )
         reopt->firstobj = reopt->run-1;
   }
   assert(norm >= 0);
   norm = SQRT(norm);

   /* normalize the coefficients */
   for(idx = 0; idx < ntransvars && norm > 0; idx++)
      reopt->objs[reopt->run-1][idx] /= norm;

   /* calculate similarity to last objective */
   if( reopt->run-1 > 1 )
   {
      /* calculate similarity to first objective */
      if( reopt->run-1 > 1 && reopt->firstobj < reopt->run-1 && reopt->firstobj >= 0 )
         reopt->simtofirstobj = reoptSimilarity(reopt, set, reopt->run-1, reopt->firstobj, transvars, ntransvars);

      /* calculate similarity to last objective */
      reopt->simtolastobj = reoptSimilarity(reopt, set, reopt->run-1, reopt->run-2, transvars, ntransvars);

      SCIPdebugMessage("new objective has similarity of %g/%g compared to first/previous.\n", reopt->simtofirstobj,
            reopt->simtolastobj);
      printf("new objective has similarity of %g/%g compared to first/previous.\n", reopt->simtofirstobj,
            reopt->simtolastobj);
   }

   SCIPdebugMessage("saved obj for run %d.\n", reopt->run);

   return SCIP_OKAY;
}

/** permute the variable and bound array randomly */
static
void permuteRandom(
   SCIP_VAR**            vars,               /**< variable array to permute */
   SCIP_Real*            vals,               /**< bound array to permute in the same order */
   int                   nvars,              /**< number of variables */
   unsigned int*         randseed            /**< seed value for the random generator */
   )
{
   SCIP_VAR* tmpvar;
   SCIP_Real tmpval;
   int end;
   int i;

   end = nvars;

   /* loop backwards through all variables (and bounds) and always swap the current last element to a random position */
   while( end > 1 )
   {
      --end;

      /* get a random position into which the last variable should be shuffled */
      i = SCIPgetRandomInt(0, end, randseed);

      /* swap the last variable and the random variable */
      tmpvar = vars[i];
      vars[i] = vars[end];
      vars[end] = tmpvar;

      /* swap the last variable and the random variable */
      tmpval = vals[i];
      vals[i] = vals[end];
      vals[end] = tmpval;
   }
}

/*
 * public methods
 */

/* ---------------- methods of general reoptimization ---------------- */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPreoptGetNRestartsGlobal
#undef SCIPreoptGetNRestartsLocal
#undef SCIPreoptGetNTotalRestartsLocal
#undef SCIPreoptGetFirstRestarts
#undef SCIPreoptGetLastRestarts
#undef SCIPreoptGetNFeasNodes
#undef SCIPreoptGetNTotalFeasNodes
#undef SCIPreoptGetNPrunedNodes
#undef SCIPreoptGetNTotalPrunedNodes
#undef SCIPreoptGetNCutoffReoptnodes
#undef SCIPreoptGetNTotalCutoffReoptnodes
#undef SCIPreoptGetNInfNodes
#undef SCIPreoptGetNTotalInfNodes
#undef SCIPreoptGetNInfSubtrees


/** returns the number of global restarts */
int SCIPreoptGetNRestartsGlobal(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->nglbrestarts;
}

/** returns the number of local restarts in the current run */
int SCIPreoptGetNRestartsLocal(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->nlocrestarts;
}

/** returns the number of local restarts over all runs */
int SCIPreoptGetNTotalRestartsLocal(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->ntotallocrestarts;
}

/** returns the number of iteration with the first global restarts */
int SCIPreoptGetFirstRestarts(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->firstrestart;
}

/** returns the number of iteration with the last global restarts */
int SCIPreoptGetLastRestarts(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->lastrestart;
}

/** returns the number of stored nodes providing an improving feasible LP solution in the current run */
int SCIPreoptGetNFeasNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->nfeasnodes;
}

/** returns the number of stored nodes providing an improving feasible LP solution over all runs */
int SCIPreoptGetNTotalFeasNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ntotalfeasnodes;
}

/** returns the number of stored nodes that exceeded the cutoff bound in the current run */
int SCIPreoptGetNPrunedNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->nprunednodes;
}

/** returns the number of stored nodes that exceeded the cutoff bound over all runs */
int SCIPreoptGetNTotalPrunedNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ntotalprunednodes;
}

/** rerturns the number of reoptimized nodes that were cutoff in the same iteration in the current run */
int SCIPreoptGetNCutoffReoptnodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ncutoffreoptnodes;
}

/** rerturns the number of reoptimized nodes that were cutoff in the same iteration over all runs */
int SCIPreoptGetNTotalCutoffReoptnodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ntotalcutoffreoptnodes;
}

/** returns the number of stored nodes with an infeasible LP in the current run */
int SCIPreoptGetNInfNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ninfnodes;
}

/** returns the number of stored nodes with an infeasible LP over all runs */
int SCIPreoptGetNTotalInfNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ntotalinfnodes;
}

/** returns the number of found infeasible subtrees */
int SCIPreoptGetNInfSubtrees(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ninfsubtrees;
}


/** constructor for the reoptimization data */
SCIP_RETCODE SCIPreoptCreate(
   SCIP_REOPT**          reopt,              /**< pointer to reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   int i;

   assert(reopt != NULL);

   SCIP_ALLOC( BMSallocMemory(reopt) );
   (*reopt)->runsize = DEFAULT_MEM_RUN;
   (*reopt)->run = 0;
   (*reopt)->simtolastobj = -2.0;
   (*reopt)->simtofirstobj = -2.0;
   (*reopt)->firstobj = -1;
   (*reopt)->currentnode = -1;
   (*reopt)->lastbranched = -1;
   (*reopt)->dualcons = NULL;
   (*reopt)->glbconss = NULL;
   (*reopt)->nglbconss = 0;
   (*reopt)->allocmemglbconss = 0;
   (*reopt)->ncheckedsols = 0;
   (*reopt)->nimprovingsols = 0;
   (*reopt)->noptsolsbyreoptsol = 0;
   (*reopt)->nglbrestarts = 0;
   (*reopt)->nlocrestarts = 0;
   (*reopt)->ntotallocrestarts = 0;
   (*reopt)->firstrestart = 0;
   (*reopt)->lastrestart = 0;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*reopt)->prevbestsols, (*reopt)->runsize) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*reopt)->objs, (*reopt)->runsize) );

   for(i = 0; i < (*reopt)->runsize; i++)
   {
      (*reopt)->objs[i] = NULL;
      (*reopt)->prevbestsols[i] = NULL;
   }

   /* clocks */
   SCIP_CALL( SCIPclockCreate(&(*reopt)->savingtime, SCIP_CLOCKTYPE_DEFAULT) );

   /* create and initialize SCIP_SOLTREE */
   SCIP_ALLOC( BMSallocMemory(&(*reopt)->soltree) );
   SCIP_CALL( createSolTree((*reopt)->soltree, blkmem) );

   /* create and initialize SCIP_REOPTTREE */
   SCIP_ALLOC( BMSallocMemory(&(*reopt)->reopttree) );
   SCIP_CALL( createReopttree((*reopt)->reopttree, set, blkmem) );

   /* create event handler for node events */
   eventhdlr = NULL;

   /* include event handler into SCIP */
   SCIP_CALL( SCIPeventhdlrCreate(&eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, NULL, NULL, NULL, NULL, eventInitsolReopt,
         eventExitsolReopt, NULL, eventExecReopt, NULL) );
   SCIP_CALL( SCIPsetIncludeEventhdlr(set, eventhdlr) );
   assert(eventhdlr != NULL);

   return SCIP_OKAY;
}

/** frees reoptimization data */
SCIP_RETCODE SCIPreoptFree(
   SCIP_REOPT**          reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          origprimal,         /**< original primal */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int p;

   assert(reopt != NULL);
   assert(*reopt != NULL);
   assert(set != NULL);
   assert(origprimal != NULL || set->stage == SCIP_STAGE_INIT);
   assert(blkmem != NULL);

   /* free reopttree */
   SCIP_CALL( freeReoptTree((*reopt)->reopttree, set, blkmem) );

   /* free solutions */
   if( set->stage >= SCIP_STAGE_PROBLEM )
   {
      for( p = (*reopt)->run-1; p >= 0; p-- )
      {
         if( (*reopt)->soltree->sols[p] != NULL )
         {
            BMSfreeBlockMemoryArray(blkmem, &(*reopt)->soltree->sols[p], (*reopt)->soltree->solssize[p]); /*lint !e866*/
            (*reopt)->soltree->sols[p] = NULL;
         }

         if( (*reopt)->objs[p] != NULL )
         {
            BMSfreeMemoryArray(&(*reopt)->objs[p]);
         }
      }
   }

   /* free solution tree */
   SCIP_CALL( freeSolTree((*reopt), set, origprimal, blkmem) );

   if( (*reopt)->dualcons != NULL )
   {
      if( (*reopt)->dualcons->varssize > 0 )
      {
         BMSfreeBlockMemoryArray(blkmem, &(*reopt)->dualcons->vals, (*reopt)->dualcons->varssize);
         BMSfreeBlockMemoryArray(blkmem, &(*reopt)->dualcons->vars, (*reopt)->dualcons->varssize);
         BMSfreeBlockMemory(blkmem, &(*reopt)->dualcons);
         (*reopt)->dualcons = NULL;
      }
   }

   if( (*reopt)->glbconss != NULL && (*reopt)->allocmemglbconss > 0 )
   {
      --(*reopt)->nglbconss;

      /* free all constraint */
      while( (*reopt)->nglbconss > 0 )
      {
         int c;
         c = (*reopt)->nglbconss;

         if( (*reopt)->glbconss[c] != NULL )
         {
            if( (*reopt)->glbconss[c]->varssize > 0 )
            {
               BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss[c]->vals, (*reopt)->glbconss[c]->varssize);
               BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss[c]->vars, (*reopt)->glbconss[c]->varssize);
               (*reopt)->glbconss[c]->varssize = 0;
            }
            BMSfreeBlockMemory(blkmem, &(*reopt)->glbconss[c]); /*lint !e866*/
         }

         --(*reopt)->nglbconss;
      }
      assert((*reopt)->nglbconss == 0);

      BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss, (*reopt)->allocmemglbconss);
      (*reopt)->allocmemglbconss = 0;
   }

   /* clocks */
   SCIPclockFree(&(*reopt)->savingtime);

   BMSfreeBlockMemoryArray(blkmem, &(*reopt)->prevbestsols, (*reopt)->runsize);
   BMSfreeMemoryArray(&(*reopt)->objs);
   BMSfreeMemory(reopt);

   return SCIP_OKAY;
}

/** returns the number of constraints added by the reoptimization plug-in */
int SCIPreoptGetNAddedConss(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   )
{
   unsigned int id;

   assert(reopt != NULL);
   assert(node != NULL);

   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* set the id to -1 if the node is not part of the reoptimization tree */
   if( SCIPnodeGetDepth(node) > 0 && id == 0 )
      return SCIPnodeGetNAddedConss(node);

   if( id >= 1 && reopt->reopttree->reoptnodes[id]->nconss > 0 )
      return MAX(SCIPnodeGetNAddedConss(node), reopt->reopttree->reoptnodes[id]->nconss); /*lint !e666*/
   else
      return SCIPnodeGetNAddedConss(node);
}

/** add a solution to the solution tree */
SCIP_RETCODE SCIPreoptAddSol(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,         /**< original primal */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SOL*             sol,                /**< solution to add */
   SCIP_Bool             bestsol,            /**< is the current solution an optimal solution? */
   SCIP_Bool*            added,              /**< pointer to store the information if the soltion was added */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   int                   run                 /**< number of the current run (1,2,...) */
   )
{
   SCIP_SOLNODE* solnode;
   SCIP_HEUR* heur;
   int insertpos;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(sol != NULL);
   assert(run > 0);

   assert(reopt->soltree->sols[run-1] != NULL);

   /* if the solution was found by reoptsols the solutions is already stored */
   heur = SCIPsolGetHeur(sol);
   if( heur != NULL && strcmp(SCIPheurGetName(heur), "reoptsols") == 0 && bestsol )
      ++reopt->noptsolsbyreoptsol;
   else if( bestsol )
      reopt->noptsolsbyreoptsol = 0;

   /* check memory */
   SCIP_CALL( ensureSolsSize(reopt, set, blkmem, reopt->soltree->nsols[run-1], run-1) );

   solnode = NULL;

   /* add solution to solution tree */
   SCIP_CALL( soltreeAddSol(reopt, set, stat, origprimal, blkmem, vars, sol, &solnode, nvars, bestsol, added) );

   if( (*added) )
   {
      assert(solnode != NULL);

      /* add solution */
      insertpos = reopt->soltree->nsols[run-1];
      reopt->soltree->sols[run-1][insertpos] = solnode;
      ++reopt->soltree->nsols[run-1];
      assert(reopt->soltree->nsols[run-1] <= set->reopt_savesols);
   }

   return SCIP_OKAY;
}

/** we want to store the optimal solution of each run in a separate array */
SCIP_RETCODE SCIPreoptAddOptSol(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SOL*             sol,                /**< solution to add */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal          /**< original primal */
   )
{
   SCIP_SOL* solcopy;

   assert(reopt != NULL);
   assert(reopt->run-1 >= 0);
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(origprimal != NULL);

   SCIP_CALL( SCIPsolCopy(&solcopy, blkmem, set, stat, origprimal, sol) );
   reopt->prevbestsols[reopt->run-1] = solcopy;

   return SCIP_OKAY;
}

/** add a new iteration after changing the objective function */
SCIP_RETCODE SCIPreoptAddRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data sturcture */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            transvars,          /**< transformed variables */
   int                   ntransvars,         /**< number of transformed variables */
   int                   size                /**< number of expected solutions */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem !=  NULL);
   assert(transvars != NULL);

   /* increase number of runs */
   ++reopt->run;

   /* check memory */
   SCIP_CALL( ensureRunSize(reopt, reopt->run, blkmem) );

   /* allocate memory */
   reopt->soltree->solssize[reopt->run-1] = size;
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->soltree->sols[reopt->run-1], size) ); /*lint !e866*/
   /* save the objective function */
   SCIP_CALL( reoptSaveNewObj(reopt, set, blkmem, transvars, ntransvars) );

   resetStats(reopt);

   return SCIP_OKAY;
}

/** get the number of checked during the reoptimization process */
int SCIPreoptGetNCheckedSols(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->ncheckedsols;
}

/** update the number of checked during the reoptimization process */
void SCIPreoptSetNCheckedSols(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   ncheckedsols        /**< number of updated solutions */
   )
{
   assert(reopt != NULL);

   reopt->ncheckedsols += ncheckedsols;
}

/** get the number of checked during the reoptimization process */
int SCIPreoptGetNImprovingSols(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->nimprovingsols;
}

/** update the number of checked during the reoptimization process */
void SCIPreoptSetNImprovingSols(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   nimprovingsols      /**< number of improving solutions */
   )
{
   assert(reopt != NULL);

   reopt->nimprovingsols += nimprovingsols;
}

/** returns number of solution of a given run */
int SCIPreoptGetNSolsRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   run                 /**< number of the run (1,2,..) */
   )
{
   assert(reopt != NULL);
   assert(0 < run && run <= reopt->runsize);

   if( reopt->soltree->sols[run-1] == NULL )
      return 0;
   else
      return reopt->soltree->nsols[run-1];
}

/** returns number of all solutions of all runs */
int SCIPreoptGetNSols(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   int nsols;
   int r;

   assert(reopt != NULL);

   nsols = 0;

   for( r = 0; r < reopt->run; r++)
      nsols += reopt->soltree->nsols[r];

   return nsols;
}

/** return the stored solutions of a given run */
SCIP_RETCODE SCIPreoptGetSolsRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   run,                /**< number of the run (1,2,...) */
   SCIP_SOL**            sols,               /**< array of solutions to fill */
   int                   solssize,           /**< length of the array */
   int*                  nsols               /**< pointer to store the number of added solutions */
   )
{
   int s;

   assert(reopt != NULL);
   assert(run > 0 && run <= reopt->run);
   assert(sols != NULL);

   assert(solssize > 0);
   *nsols = 0;

   for( s = 0; s < reopt->soltree->nsols[run-1]; s++ )
   {
      if( !reopt->soltree->sols[run-1][s]->updated )
         ++(*nsols);
   }

   if( solssize < (*nsols) )
      return SCIP_OKAY;

   (*nsols) = 0;
   for( s = 0; s < reopt->soltree->nsols[run-1]; s++ )
   {
      if( !reopt->soltree->sols[run-1][s]->updated )
      {
         sols[*nsols] = reopt->soltree->sols[run-1][s]->sol;
         reopt->soltree->sols[run-1][s]->updated = TRUE;
         ++(*nsols);
      }
   }

   return SCIP_OKAY;
}

/** returns the number of saved solutions overall runs */
int SCIPreoptGetNSavedSols(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   int nsavedsols;

   assert(reopt != NULL);
   assert(reopt->soltree->root != NULL);

   nsavedsols = 0;

   if( reopt->soltree->root->lchild != NULL
    || reopt->soltree->root->rchild != NULL)
      nsavedsols = soltreeNInducedSols(reopt->soltree->root);

   return nsavedsols;
}

/** check if the reoptimization process should be (locally) restarted.
 *
 *  first, we check whether the current node is the root node, e.g., node == NULL. in this case, we do not need to calculate
 *  the similarity again. we trigger a restart if
 *    1. the objective function has changed too much
 *    2. the number of stored nodes is exceeded
 *    3. the last n optimal solutions were found by heur_reoptsols (in this case, the stored tree was only needed to
 *    prove the optimality and this can be probably faster by solving from scratch)
 *
 *  if the current node is different to the root node we calculate the local similarity, i.e., exclude all variable
 *  that are already fixed by bounding.
 */
SCIP_RETCODE SCIPreoptCheckRestart(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< current node of the branch and bound tree (or NULL) */
   SCIP_VAR**            transvars,          /**< transformed problem variables */
   int                   ntransvars,         /**< number of transformed problem variables */
   SCIP_Bool*            restart             /**< pointer to store if the reoptimization process should be restarted */
   )
{
   SCIP_Real sim;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(transvars != NULL);
   assert(ntransvars >= 0);

   sim = 1.0;
   *restart = FALSE;

   /* check if the whole reoptimization process should start from scratch */
   if( node == NULL )
   {
      if( reopt->run > 0 && set->reopt_objsimdelay > -1.0 )
      {
         sim = reopt->simtolastobj;
      }

      /* check similarity */
      if( SCIPsetIsFeasLT(set, sim, set->reopt_objsimdelay) )
      {
         SCIPdebugMessage("-> restart reoptimization (objective functions are not similar enough)\n");
         *restart = TRUE;
      }
      /* check size of the reoptimization tree */
      else if( reopt->reopttree->nreoptnodes > set->reopt_maxsavednodes )
      {
         SCIPdebugMessage("-> restart reoptimization (node limit reached)\n");
         *restart = TRUE;
      }
      /* check if the tree was only needed to prove optimality */
      else if( reopt->noptsolsbyreoptsol >= set->reopt_forceheurrestart )
      {
         SCIPdebugMessage("-> restart reoptimization (found last %d optimal solutions by <reoptsols>)\n",
               reopt->noptsolsbyreoptsol);
         reopt->noptsolsbyreoptsol = 0;
         *restart = TRUE;
      }

      if( *restart )
      {
         /* trigger a restart */
         SCIP_CALL( reoptRestart(reopt, set, blkmem) );
      }
   }
   /* check for a local restart, ie, start the solving process of an inner node from scatch */
   else
   {
      SCIP_CALL( reoptCheckLocalRestart(reopt, set, blkmem, node, transvars, ntransvars, restart) );
   }
   return SCIP_OKAY;
}

/** returns the similarity to the previous objective function, if no exist return -2.0 */
SCIP_Real SCIPreoptGetSimToPrevious(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);
   return reopt->simtolastobj;
}

/**returns the similarity to the first objective different to the zero-function function, if no exist return -2.0 */
SCIP_Real SCIPreoptGetSimToFirst(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);
   return reopt->simtofirstobj;
}

/** return the similarity between two of objective functions of two given runs */
SCIP_Real SCIPreoptGetSimilarity(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   run1,               /**< number of the first run */
   int                   run2,               /**< number of the second run */
   SCIP_VAR**            transvars,          /**< original problem variables */
   int                   ntransvars          /**< number of original problem variables */
   )
{
   assert(reopt != NULL);
   assert(run1 > 0 && run1 <= reopt->run);
   assert(run2 > 0 && run2 <= reopt->run);
   assert(transvars != NULL);
   assert(ntransvars >= 0);

   return reoptSimilarity(reopt, set, run1-1, run2-1, transvars, ntransvars);
}

/** returns the best solution of the last run */
SCIP_SOL* SCIPreoptGetLastBestSol(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);
   assert(reopt->prevbestsols != NULL);

   if( reopt->run-2 < 0 )
      return NULL;
   else
   {
      assert(reopt->prevbestsols[reopt->run-2] != NULL);
      return reopt->prevbestsols[reopt->run-2];
   }
}

/** returns the node of the reoptimization tree corresponding to the unique @p id */
SCIP_REOPTNODE* SCIPreoptGetReoptnode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   unsigned int          id                  /**< unique id */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   return reopt->reopttree->reoptnodes[id];
}

/** returns the coefficient of variable with index @p idx in run @p run */
SCIP_Real SCIPreoptGetOldObjCoef(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   run,                /**< number of the run (1,2,...) */
   int                   idx                 /**< index of variable */
   )
{
   assert(reopt != NULL);
   assert(0 < run && run <= reopt->runsize);

   return reopt->objs[run-1][idx];
}

/** return the best solution of a given run.
 *
 *  @note the return solution is part of the original space.
 */
SCIP_SOL* SCIPreoptGetBestSolRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   run                 /**< number of the run (1,2,...) */
   )
{
   assert(reopt != NULL);
   assert(0 < run && run <= reopt->run);

   return reopt->prevbestsols[run-1];
}

/** reset marks of stored solutions to not updated */
void SCIPreoptResetSolMarks(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);
   assert(reopt->soltree != NULL);
   assert(reopt->soltree->root != NULL);

   if( reopt->soltree->root->rchild != NULL )
      soltreeResetMarks(reopt->soltree->root->rchild);
   if( reopt->soltree->root->lchild )
      soltreeResetMarks(reopt->soltree->root->lchild);
}

/** returns the number of stored nodes in the subtree induced by @p node */
int SCIPreoptGetNNodes(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   )
{
   unsigned int id;

   assert(reopt != NULL);

   if( node == NULL || SCIPnodeGetDepth(node) == 0 )
      return reopt->reopttree->nreoptnodes;

   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* set the id to -1 if the node is not part of the reoptimization tree */
   if( SCIPnodeGetDepth(node) > 0 && id == 0 )
      return 0;

   assert(0 < id && id < reopt->reopttree->reoptnodessize);

   return reopttreeGetNNodes(reopt->reopttree, id);
}

/* ---------------- methods of general reoptimization nodes ---------------- */

/** In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPreoptnodeGetNVars
#undef SCIPreoptnodeGetNConss
#undef SCIPreoptnodeGetNDualBoundChgs
#undef SCIPreoptnodeGetNChildren
#undef SCIPreoptnodeGetLowerbound
#undef SCIPreoptnodeGetType

/** returns the number of bound changes stored in the reopttree at ID id */
int SCIPreoptnodeGetNVars(
   SCIP_REOPTNODE*       reoptnode           /**< node of the roepttree */
   )
{
   assert(reoptnode != NULL);

   return reoptnode->nvars + reoptnode->nafterdualvars;
}

/** returns the number of bound changes at the node stored at ID id */
int SCIPreoptnodeGetNConss(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   )
{
   assert(reoptnode != NULL);

   return reoptnode->nconss;
}

/** returns the number of stored bound changes based on dual information in the reopttree at ID id */
int SCIPreoptnodeGetNDualBoundChgs(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   )
{
   assert(reoptnode != NULL);

   if( reoptnode->dualconscur == NULL )
      return 0;
   else
      return reoptnode->dualconscur->nvars;
}

/** returns the number of child nodes of @p reoptnode */
int SCIPreoptnodeGetNChildren(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimizzation tree */
   )
{
   assert(reoptnode != NULL);

   return reoptnode->nchilds;
}

/** return the lower bound stored at @p ID id */
SCIP_Real SCIPreoptnodeGetLowerbound(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   )
{
   assert(reoptnode != NULL);

   return reoptnode->lowerbound;
}

/** returns the type of the @p reoptnode */
SCIP_REOPTTYPE SCIPreoptnodeGetType(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   )
{
   assert(reoptnode != NULL);

   return (SCIP_REOPTTYPE)reoptnode->reopttype;
}

/** create the constraint which splits the node stored at ID id on the basis of
 *  the stored dual information.
 */
void SCIPreoptnodeGetSplitCons(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR**            vars,               /**< array to store the variables of the constraint */
   SCIP_Real*            vals,               /**< array to store the coefficients of the variables */
   REOPT_CONSTYPE*       constype,           /**< type of the constraint */
   int                   conssize,           /**< size of the arrays */
   int*                  nvars               /**< pointer to store the size of the constraints */
   )
{
   int v;

   assert(reoptnode != NULL);
   assert(vars != NULL);
   assert(vals != NULL);

   (*nvars) = reoptnode->dualconscur == NULL ? 0 : reoptnode->dualconscur->nvars;

   if( *nvars == 0 || *nvars > conssize )
      return;

   /* copy the variable information */
   for(v = 0; v < *nvars; v++)
   {
      vars[v] = reoptnode->dualconscur->vars[v];
      vals[v] = reoptnode->dualconscur->vals[v];
   }

   *constype = reoptnode->dualconscur->constype;

   return;
}

/** returns all added constraints at ID id */
void SCIPreoptnodeGetConss(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR***           vars,               /**< 2-dim array of variables */
   SCIP_Real**           vals,               /**< 2-dim array of values */
   int                   mem,                /**< allocated memory for constraints */
   int*                  nconss,             /**< pointer to store the number of constraints */
   int*                  nvars               /**< pointer to store the number of variables */
   )
{
   int c;

   assert(reoptnode != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(nvars != NULL);


   (*nconss) = reoptnode->nconss;

   if( mem < *nconss )
      return;

   for(c = 0; c < *nconss; c++)
   {
      assert(vars[c] != NULL);
      assert(vals[c] != NULL);

      vars[c] = reoptnode->conss[c]->vars;
      vals[c] = reoptnode->conss[c]->vals;
      nvars[c] = reoptnode->conss[c]->nvars;
   }

   return;
}

/** set the parent id */
void SCIPreoptnodeSetParentID(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   unsigned int          parentid            /**< id of the parent node */
   )
{
   assert(reoptnode != NULL);
   assert(parentid <= 536870912); /* id can be at least 2^29 */

   reoptnode->parentID = parentid;
}

/** returns the number of leaf nodes of the subtree induced by @p node (of the whole tree if node == NULL) */
int SCIPreoptGetNLeaves(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node                /**< node of the search tree (or NULL) */
   )
{
   int nleaves;
   unsigned int id;
   int i;

   assert(reopt != NULL);

   nleaves = 0;
   id = (node == NULL) ? 0 : SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* return if the node is not part of the reoptimization tree */
   if( node != NULL && SCIPnodeGetDepth(node) > 0 && id == 0 )
      return nleaves;

   for(i = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++)
   {
      unsigned int childid;

      childid = reopt->reopttree->reoptnodes[id]->childids[i]; /*lint !e713*/
      assert(childid < reopt->reopttree->reoptnodessize);

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
         ++nleaves;
      else
         nleaves += reoptGetNLeaves(reopt, childid);
   }

   return nleaves;
}

/** save information that given node is infeasible */
SCIP_RETCODE SCIPreoptAddInfNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node                /**< node of the search tree */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);

   if( set->reopt_sepaglbinfsubtrees )
   {
      SCIP_CALL( saveGlobalCons(reopt, set, blkmem, node, REOPT_CONSTYPE_INFSUBTREE) );
   }

   ++reopt->reopttree->ninfnodes;
   ++reopt->reopttree->ntotalinfnodes;

   return SCIP_OKAY;
}

/** check the reason for cut off a node and if necessary store the node */
SCIP_RETCODE SCIPreoptCheckCutoff(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_EVENTTYPE        eventtype,          /**< eventtype */
   SCIP_LPSOLSTAT        lpsolstat,          /**< solution status of the LP */
   SCIP_Bool             isrootnode,         /**< the node is the root */
   SCIP_Bool             isfocusnode,        /**< the node is the current focus node */
   SCIP_Real             lowerbound,         /**< lower bound of the node */
   int                   effectiverootdepth  /**< effective root depth */
   )
{
   SCIP_Bool strongbranched;

   assert(reopt != NULL);
   assert(node != NULL);
   assert(eventtype == SCIP_EVENTTYPE_NODEBRANCHED
       || eventtype == SCIP_EVENTTYPE_NODEFEASIBLE
       || eventtype == SCIP_EVENTTYPE_NODEINFEASIBLE);

   if( reopt->lastseennode == SCIPnodeGetNumber(node) )
      return SCIP_OKAY;

   reopt->lastseennode = SCIPnodeGetNumber(node);

   SCIPdebugMessage("catch event %x for node %lld\n", eventtype, SCIPnodeGetNumber(node));

   /* case 1: the current node is the root node
    * we can skip if the root is (in)feasible or branched w/o bound
    * changes based on dual information.
    *
    * case 2: we need to store the current node if it contains
    * bound changes based on dual information or is a leave node
    */

   if( isrootnode )
   {
      if( SCIPreoptGetNDualBndchgs(reopt, node) > 0 )
      {
         goto CHECK;
      }
      else if( eventtype == SCIP_EVENTTYPE_NODEBRANCHED )
      {
         /* store or update the information */
         SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_TRANSIT, TRUE, isrootnode, lowerbound) );
      }
      else if( eventtype == SCIP_EVENTTYPE_NODEFEASIBLE )
      {
         /* delete saved dual information which would lead to split the node in a further iteration */
         SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );

         /* store or update the information */
         SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_FEASIBLE, FALSE, isrootnode, lowerbound) );
      }
      else if( eventtype == SCIP_EVENTTYPE_NODEINFEASIBLE )
      {
         /* delete saved dual information which would lead to split the node in a further iteration */
         SCIP_CALL( SCIPreoptResetDualBndchgs(reopt, node, blkmem) );

         /* store or update the information */
         SCIP_CALL( addNode(reopt, set, blkmem, node, reopt->currentnode == 1 ? SCIP_REOPTTYPE_INFSUBTREE : SCIP_REOPTTYPE_PRUNED, FALSE,
               isrootnode, lowerbound) );
      }

      assert(reopt->currentnode == -1);
      assert(reopt->dualcons == NULL || reopt->dualcons->nvars == 0);

      return SCIP_OKAY;
   }

  CHECK:

   if( effectiverootdepth == SCIPnodeGetDepth(node) )
   {
      strongbranched = SCIPreoptGetNDualBndchgs(reopt, node) > 0 ? TRUE : FALSE;
   }
   else
   {
      strongbranched = SCIPnodeGetNDualBndchgs(node) > 0 ? TRUE : FALSE;
   }

   SCIPdebugMessage("check the reason of cutoff for node %lld:\n", SCIPnodeGetNumber(node));
   SCIPdebugMessage(" -> focusnode       : %s\n", isfocusnode ? "yes" : "no");
   SCIPdebugMessage(" -> depth           : %d (eff. %d)\n", SCIPnodeGetDepth(node), effectiverootdepth);
   SCIPdebugMessage(" -> strong branched : %s\n", strongbranched ? "yes" : "no");
   SCIPdebugMessage(" -> LP lpsolstat    : %d\n", lpsolstat);

   switch( eventtype ) {
      case SCIP_EVENTTYPE_NODEFEASIBLE:
         /* current node has to be the eventnode */
         assert(isfocusnode);

         SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_FEASIBLE);

         /* delete strong branching information of some exists */
         deleteLastDualBndchgs(reopt);

         SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_FEASIBLE, FALSE, isrootnode, lowerbound) );
         break;

      case SCIP_EVENTTYPE_NODEINFEASIBLE:
         /* we have to check if the current node is the event node.
          * if the current node is not the event node, we have to save this node, else we have to
          * look at LP lpsolstat and decide.
          */
         if( isfocusnode )
         {
            /* an after-branch heuristic says NODEINFEASIBLE, maybe the cutoff bound is reached.
             * because the node is already branched we have all children and can delete this node. */
            if( SCIPnodeGetNumber(node) == reopt->lastbranched )
            {
               deleteLastDualBndchgs(reopt);
               break;
            }

            /*
             * if the node is strong branched we possible detect an infeasible subtree, if not,
             * the whole node is either infeasible or exceeds the cutoff bound.
             */
            if( strongbranched )
            {
               /*
                * 1. the LP is not solved or infeasible: the subnode is infeasible and can be discarded
                *    because either the LP proves infeasibility or a constraint handler.
                *    We have to store an infeasible subtree constraint
                * 2. the LP exceeds the objective limit, we have to store the node and can delete the
                *    strong branching information
                */
               if( lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE )
               {
                  /* add a dummy variable, because the bound changes were not global in the
                   * sense of effective root depth
                   */
                  if( SCIPnodeGetDepth(node) > effectiverootdepth )
                  {
                     SCIP_CALL( SCIPreoptAddDualBndchg(reopt, set, blkmem, node, NULL, 0.0, 1.0) );
                  }

                  SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_INFSUBTREE);
                  SCIPdebugMessage(" -> new constype    : %d\n", REOPT_CONSTYPE_INFSUBTREE);

                  /* save the node as a strong branched node */
                  SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_INFSUBTREE, FALSE, isrootnode, lowerbound) );
               }
               else
               {
                  assert(SCIP_LPSOLSTAT_OBJLIMIT || SCIP_LPSOLSTAT_OPTIMAL || SCIP_LPSOLSTAT_NOTSOLVED);

                  SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_PRUNED);

                  /* delete strong branching information of some exists */
                  deleteLastDualBndchgs(reopt);

                  SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_PRUNED, FALSE, isrootnode, lowerbound) );
               }
            }
            else
            {
               /*
                * 1. the LP is not solved or infeasible: the whole node is infeasible and can be discarded
                *    because either the LP proves infeasibility or a constraint handler.
                * 2. the LP exceeds the objective limit, we have to store the node and can delete the
                *    strong branching information
                */
               if( lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE )
               {
                  /* save the information of an infeasible node */
                  SCIPdebugMessage(" -> new reopttype   : infeasible\n");
                  SCIP_CALL( SCIPreoptAddInfNode(reopt, set, blkmem, node) );
               }
               else
               {
                  SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_PRUNED);

                  /* store the node */
                  SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_PRUNED, TRUE, isrootnode, lowerbound) );
               }
            }
         }
         else
         {
            SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_PRUNED);

            /* if the node was created by branch_nodereopt, nothing happens */
            SCIP_CALL(addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_PRUNED, TRUE, isrootnode, lowerbound) );

         }
         break;

      case SCIP_EVENTTYPE_NODEBRANCHED:
         /* current node has to be the eventnode */
         assert(isfocusnode);

         reopt->lastbranched = SCIPnodeGetNumber(node);

         /* we have to check the depth of the current node. if the depth is equal to the effective
          * root depth, then all information about bound changes based on dual information already exists,
          * else we have to look at the domchg-data-structure.*/
         if (SCIPnodeGetDepth(node) == effectiverootdepth)
         {
            /* Save the node if there are added constraints, because this means the node is a copy create by the
             * reoptimization plug-in and contains at least one logic-or-constraint */
            if( strongbranched )
            {
               SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_STRBRANCHED);
               SCIPdebugMessage(" -> new constype    : %d\n", REOPT_CONSTYPE_STRBRANCHED);
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_STRBRANCHED, TRUE, isrootnode, lowerbound) );
            }
            else if( SCIPreoptGetNAddedConss(reopt, node) > 0 )
            {
               SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_LOGICORNODE);
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_LOGICORNODE, TRUE, isrootnode, lowerbound) );
            }
            else
            {
               SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_TRANSIT);
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_TRANSIT, TRUE, isrootnode, lowerbound) );
            }
         }
         else
         {
            /* we only branch on binary variables and var == NULL indicates memory allocation w/o saving information.
             *
             * we have to do this in the following order:
             * 1) all bound-changes are local, thats way we have to mark the node to include bound changes based
             *    on dual information.
             * 2) save or update the node.
             */
            if( strongbranched )
            {
               SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_STRBRANCHED);
               SCIPdebugMessage(" -> new constype    : %d\n", REOPT_CONSTYPE_STRBRANCHED);
               SCIP_CALL( SCIPreoptAddDualBndchg(reopt, set, blkmem, node, NULL, 0.0, 1.0) );
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_STRBRANCHED, TRUE, isrootnode, lowerbound) );
            }
            else if( SCIPreoptGetNAddedConss(reopt, node) > 0 )
            {
               SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_LOGICORNODE);
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_LOGICORNODE, TRUE, isrootnode, lowerbound) );
            }
            else
            {
               SCIPdebugMessage(" -> new reopttype   : %d\n", SCIP_REOPTTYPE_TRANSIT);
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_TRANSIT, TRUE, isrootnode, lowerbound) );
            }
         }
         break;

      default:
         break;
   }

   assert(reopt->currentnode == -1);
   assert(reopt->dualcons == NULL || reopt->dualcons->nvars == 0);

   return SCIP_OKAY;
}

/** store bound change based on dual information */
SCIP_RETCODE SCIPreoptAddDualBndchg(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_VAR*             var,                /**< variables */
   SCIP_Real             newval,             /**< new bound */
   SCIP_Real             oldval              /**< old bound */
   )
{
   SCIP_Real constant;
   SCIP_Real scalar;

   assert(reopt != NULL);
   assert(node != NULL);

   constant = 0.0;
   scalar = 1.0;

   /* If var == NULL, we save all information by calling SCIPreoptNodeFinished().
    * In that case, all bound changes were not global and we can find them within the
    * domchg data structure.
    * Otherwise, we allocate memory and store the information.
    */
   if( var != NULL )
   {
      int allocmem;

      assert(SCIPsetIsFeasEQ(set, newval, 0.0) || SCIPsetIsFeasEQ(set, newval, 1.0));

      allocmem = (reopt->dualcons == NULL || reopt->dualcons->varssize == 0) ? DEFAULT_MEM_DUALCONS : reopt->dualcons->varssize+2;

      /* allocate memory of necessary */
      SCIP_CALL( checkMemDualCons(reopt, blkmem, allocmem) );

      assert(reopt->dualcons->varssize > 0);
      assert(reopt->dualcons->nvars >= 0);
      assert(reopt->currentnode == -1 || reopt->dualcons->nvars > 0);
      assert((reopt->dualcons->nvars > 0 && reopt->currentnode == SCIPnodeGetNumber(node))
           || reopt->dualcons->nvars == 0);

      reopt->currentnode = SCIPnodeGetNumber(node);

      /* transform into the original space and then save the bound change */
      SCIP_CALL(SCIPvarGetOrigvarSum(&var, &scalar, &constant));
      newval = (newval - constant) / scalar;
      oldval = (oldval - constant) / scalar;

      assert(SCIPvarIsOriginal(var));

      reopt->dualcons->vars[reopt->dualcons->nvars] = var;
      reopt->dualcons->vals[reopt->dualcons->nvars] = newval;
      ++reopt->dualcons->nvars;

      SCIPdebugMessage(">> store bound change of <%s>: %g -> %g\n", SCIPvarGetName(var), oldval, newval);
   }
   else
   {
      assert(reopt->currentnode == -1);
      assert(reopt->dualcons == NULL || reopt->dualcons->nvars == 0);

      reopt->currentnode = SCIPnodeGetNumber(node);
   }

   return SCIP_OKAY;
}

/** returns the number of bound changes based on dual information */
int SCIPreoptGetNDualBndchgs(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   )
{
   int ndualbndchgs;

   assert(reopt != NULL);
   assert(node != NULL);

   ndualbndchgs = 0;

   if( SCIPnodeGetNumber(node) == reopt->currentnode )
   {
      assert(reopt->dualcons != NULL);
      ndualbndchgs = reopt->dualcons->nvars;
   }

   return ndualbndchgs;
}

/** returns the child nodes of @p node that need to be reoptimized next or NULL if @p node is a leaf */
SCIP_RETCODE SCIPreoptGetChildIDs(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int*         childs,             /**< array to store the child ids */
   int                   childssize,         /**< size of the childs array */
   int*                  nchilds             /**< pointer to store the number of child nodes */
   )
{
   SCIP_Bool runagain;
   unsigned int id;

   assert(reopt != NULL);
   assert(childssize > 0 && childs != NULL);

   (*nchilds) = 0;

   if( node == NULL )
      id = 0;
   else
      id = SCIPnodeGetReoptID(node);

   assert(id >= 1 || SCIPnodeGetDepth(node) == 0);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   /* check if there are redundant bound changes or infeasible nodes */
   runagain = TRUE;
   while( runagain && reopt->reopttree->reoptnodes[id]->nchilds > 0 )
   {
      SCIP_CALL( dryBranch(reopt, set, blkmem, &runagain, id) );
   }

   /* return the list of child nodes if some exists; otherwise return NULL */
   if( reopt->reopttree->reoptnodes[id]->childids != NULL && reopt->reopttree->reoptnodes[id]->nchilds > 0 )
   {
      int c;

      (*nchilds) = reopt->reopttree->reoptnodes[id]->nchilds;

      if( childssize < *nchilds )
         return SCIP_OKAY;

      for( c = 0; c < *nchilds; c++ )
      {
         childs[c] = reopt->reopttree->reoptnodes[id]->childids[c];
      }
   }

   return SCIP_OKAY;
}

/** returns all leaves of the subtree induced by @p node */
SCIP_RETCODE SCIPreoptGetLeaves(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int*         leaves,             /**< array to the the ids */
   int                   leavessize,         /**< size of leaves array */
   int*                  nleaves             /**< pointer to store the number of leav node */
   )
{
   unsigned int id;
   int i;

   assert(reopt != NULL);
   assert(leavessize > 0 && leaves != NULL);
   assert((*nleaves) >= 0);

   if( node == NULL )
      id = 0;
   else
      id = SCIPnodeGetReoptID(node);

   /* return if the node is not part of the reoptimization tree */
   if( id == 0 && node != NULL )
   {
      (*nleaves) = 0;
      return SCIP_OKAY;
   }

   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   for( i = 0; i < leavessize; i++ )
      leaves[i] = 0;

   for( i = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++ )
   {
      unsigned int childid;

      assert(*nleaves + 1 <= leavessize);

      childid = reopt->reopttree->reoptnodes[id]->childids[i];
      assert(childid < reopt->reopttree->reoptnodessize);

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
      {
         leaves[(*nleaves)] = reopt->reopttree->reoptnodes[id]->childids[i];
         ++(*nleaves);
      }
      else
      {
         int nleaves2;

         nleaves2 = 0;
         SCIP_CALL( reoptGetLeaves(reopt, childid, &leaves[*nleaves], leavessize - (*nleaves), &nleaves2) );
         (*nleaves) += nleaves2;
      }
   }

   return SCIP_OKAY;
}

/** add all unprocessed nodes to the reoptimization tree */
SCIP_RETCODE SCIPreoptSaveOpenNodes(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE**           leaves,             /**< array of open leave nodes */
   int                   nleaves,            /**< number of open leave nodes */
   SCIP_NODE**           childs,             /**< array of open children nodes */
   int                   nchilds,            /**< number of open leave nodes */
   SCIP_NODE**           siblings,           /**< array of open sibling nodes */
   int                   nsiblings           /**< number of open leave nodes */
   )
{
   int n;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(nleaves >= 0);
   assert(nleaves == 0 || leaves != NULL);
   assert(nchilds >= 0);
   assert(nchilds == 0 || childs != NULL);
   assert(nsiblings >= 0);
   assert(nsiblings == 0 || siblings != NULL);

   SCIPdebugMessage("save unprocessed nodes (%d leaves, %d children, %d siblings)\n", nleaves, nchilds, nsiblings);

   /* save open leaves */
   for( n = 0; n < nleaves; n++ )
   {
      SCIP_CALL( addNode(reopt, set, blkmem, leaves[n], SCIP_REOPTTYPE_PRUNED, FALSE, FALSE,
            SCIPnodeGetLowerbound(leaves[n])) );
   }

   /* save open children */
   for( n = 0; n < nchilds; n++ )
   {
      SCIP_CALL( addNode(reopt, set, blkmem, childs[n], SCIP_REOPTTYPE_PRUNED, FALSE, FALSE,
            SCIPnodeGetLowerbound(childs[n])) );
   }

   /* save open siblings */
   for( n = 0; n < nsiblings; n++ )
   {
      SCIP_CALL( addNode(reopt, set, blkmem, siblings[n], SCIP_REOPTTYPE_PRUNED, FALSE, FALSE,
            SCIPnodeGetLowerbound(siblings[n])) );
   }

   return SCIP_OKAY;
}

/** reset the complete tree and set the given search frontier */
SCIP_RETCODE SCIPreoptApplyCompression(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPTNODE**      representatives,    /**< array of representatives */
   int                   nrepresentatives,   /**< number of representatives */
   SCIP_Bool*            success             /**< pointer to store if the method was successful */
   )
{
   SCIP_REOPTTREE* reopttree;
   unsigned int id;
   int r;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(representatives != NULL);
   assert(nrepresentatives > 0);

   reopttree = reopt->reopttree;

   /* reset the current search tree */
   SCIP_CALL( reoptResetTree(reopt, set, blkmem, FALSE) );
   assert(reopttree->nreoptnodes == 0);

   /* create a new root node */
   id = 0;
   SCIP_CALL( createReoptnode(reopttree, set, blkmem, id) );

   /* set the reopttype */
   reopttree->reoptnodes[0]->reopttype = (unsigned int)SCIP_REOPTTYPE_TRANSIT;

   /* add all representatives */
   for( r = 0; r < nrepresentatives; r++ )
   {
      /* get an empty slot*/
      id = (unsigned int) (size_t) SCIPqueueRemove(reopttree->openids);
      assert(1 <= id && id < reopttree->reoptnodessize);
      assert(reopttree->reoptnodes[id] == NULL);

      SCIP_CALL( createReoptnode(reopttree, set, blkmem, id) );
      assert(reopttree->reoptnodes[id] != NULL);

      /* set the new node
       * 1. copy all variables, bounds, and boundtypes
       * 2. copy all constraints
       * 3. set the parent relation */
      if( representatives[r]->nvars > 0 )
      {
         assert(representatives[r]->nvars <= representatives[r]->varssize);
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->vars, representatives[r]->vars,
               representatives[r]->nvars) );
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->varbounds,
               representatives[r]->varbounds, representatives[r]->nvars) );
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->varboundtypes,
               representatives[r]->varboundtypes, representatives[r]->nvars) );
         reopttree->reoptnodes[id]->varssize = representatives[r]->varssize;
         reopttree->reoptnodes[id]->nvars = representatives[r]->nvars;
      }

      if( representatives[r]->nconss > 0 )
      {
         int c;

         assert(representatives[r]->nconss <= representatives[r]->consssize);

         for( c = 0; c < representatives[r]->nconss; c++ )
         {
            SCIP_CALL( SCIPreoptnodeAddCons(reopttree->reoptnodes[id], blkmem, representatives[r]->conss[c]->vars,
                  representatives[r]->conss[c]->vals, representatives[r]->conss[c]->nvars,
                  representatives[r]->conss[c]->constype) );
         }
      }

      reopttree->reoptnodes[id]->parentID = representatives[r]->parentID; /*lint !e732*/

      assert(reopttree->reoptnodes[id]->parentID == 0);
      assert(reopttree->reoptnodes[id]->nvars >= 0);
      assert(reopttree->reoptnodes[id]->nvars <= reopttree->reoptnodes[id]->varssize);
      assert(reopttree->reoptnodes[id]->nconss >= 0);

      /* set the reopttype */
      if( reopttree->reoptnodes[id]->nconss == 0 )
         reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_LEAF;
      else
         reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_LOGICORNODE;

      /* add the representative as a child of the root */
      SCIP_CALL( reoptAddChild(reopttree, blkmem, 0, id) );
   }

   SCIPdebugMessage("-> new tree consists of %d nodes, the root has %d child nodes.\n",
         reopttree->nreoptnodes, reopttree->reoptnodes[0]->nchilds);

   (*success) = TRUE;

   return SCIP_OKAY;
}

/** splits the root into several nodes and moves the child nodes of the root to one of the created nodes */
SCIP_RETCODE SCIPreoptSplitRoot(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          randseed,           /**< seed value for random generator */
   int*                  ncreatedchilds,     /**< pointer to store the number of created nodes */
   int*                  naddedconss         /**< pointer to store the number added constraints */
   )
{
   SCIP_REOPTTREE* reopttree;
   LOGICORDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   unsigned int id;
   int nbndchgs;
   int nchilds;
   int nvars;
   int v;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(reopt->reopttree->reoptnodes[0] != NULL);
   assert(reopt->reopttree->reoptnodes[0]->dualfixing);
   assert(reopt->reopttree->reoptnodes[0]->reopttype == (unsigned int)SCIP_REOPTTYPE_STRBRANCHED);
   assert(set != NULL);
   assert(blkmem != NULL);

   reopttree = reopt->reopttree;

   nchilds = reopttree->reoptnodes[0]->nchilds;

   assert(reopttree->reoptnodes[0]->dualconscur != NULL);
   nbndchgs = reopttree->reoptnodes[0]->dualconscur->nvars;

   vars = NULL;
   vals = NULL;
   nvars = 0;

   (*ncreatedchilds) = 0;
   (*naddedconss) = 0;

   if( !set->reopt_usesplitcons )
   {
      vars = reopttree->reoptnodes[0]->dualconscur->vars;
      vals = reopttree->reoptnodes[0]->dualconscur->vals;
      nvars = reopttree->reoptnodes[0]->dualconscur->nvars;

      /* calculate the order of the variables */
      switch (set->reopt_varorderinterdiction) {
         case 'd':
            break;

         case 'r':
            permuteRandom(vars, vals, nvars, &randseed);
            break;

         default:
            return SCIP_INVALIDDATA;
      }
   }

   /* create a node with all variables fixed, i.e., reconstruct the root of the last iteration */

   /* ensure that two free slots are available  */
   SCIP_CALL( reopttreeCheckMemory(reopttree, blkmem) );
   id = (unsigned int) (size_t) SCIPqueueRemove(reopttree->openids);

   assert(0 < id && id < reopt->reopttree->reoptnodessize);
   assert(reopttree->reoptnodes[id] == NULL || reopttree->reoptnodes[id]->nvars == 0);

   /*   1. create the node
    *   2. add all bound changes
    *   3. move all child nodes to id
    *   4. add id as a child of the root node
    */
   SCIP_CALL( createReoptnode(reopttree, set, blkmem, id) );
   reopttree->reoptnodes[id]->parentID = 0;
   reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_TRANSIT;

   /* check memory */
   SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], blkmem, nbndchgs, nchilds, 0) );
   assert(reopttree->reoptnodes[id]->varssize >= nbndchgs);
   assert(reopttree->reoptnodes[id]->nvars == 0);
   assert(reopttree->reoptnodes[id]->vars != NULL);
   assert(reopttree->reoptnodes[id]->varbounds != NULL);
   assert(reopttree->reoptnodes[id]->varboundtypes != NULL);

   /* copy bounds */
   for(v = 0; v < nbndchgs; v++)
   {
      reopttree->reoptnodes[id]->vars[v] = reopttree->reoptnodes[0]->dualconscur->vars[v];
      reopttree->reoptnodes[id]->varbounds[v] = reopttree->reoptnodes[0]->dualconscur->vals[v];
      reopttree->reoptnodes[id]->varboundtypes[v] = SCIPsetIsFeasEQ(set, reopttree->reoptnodes[0]->dualconscur->vals[v], 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
      ++reopttree->reoptnodes[id]->nvars;
   }
   assert(reopttree->reoptnodes[id]->nvars == reopttree->reoptnodes[0]->dualconscur->nvars);

   /* move the children */
   SCIP_CALL( reoptMoveIDs(reopttree, blkmem, 0, id) );
   assert(reopttree->reoptnodes[0]->nchilds == 0);

   /* add the new reoptimization node as a child of the root node */
   SCIP_CALL( reoptAddChild(reopttree, blkmem, 0, id) );

   ++(*ncreatedchilds);

   if( set->reopt_usesplitcons )
   {
      assert(*ncreatedchilds == 1);

      /* ensure that there is a free slots */
      SCIP_CALL( reopttreeCheckMemory(reopttree, blkmem) );
      id = (unsigned int) (size_t) SCIPqueueRemove(reopttree->openids);
      assert(0 < id && id < reopt->reopttree->reoptnodessize);

      /* 1. create the node
       * 2. add the constraint to ensure that at least one
       *    variable gets different
       * 3. add id as a child of the root node
       */
      SCIP_CALL( createReoptnode(reopttree, set, blkmem, id) );
      reopttree->reoptnodes[id]->parentID = 0;
      reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_LOGICORNODE;

      /* create the constraint */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &consdata) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &consdata->vars, reopttree->reoptnodes[0]->dualconscur->vars, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &consdata->vals, reopttree->reoptnodes[0]->dualconscur->vals, nbndchgs) );

      consdata->varssize = nbndchgs;
      consdata->nvars = nbndchgs;
      consdata->constype = REOPT_CONSTYPE_STRBRANCHED;

      ++(*naddedconss);

      /* check memory for added constraints */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], blkmem, 0, 0, 1) );

      /* add the constraint */
      reopttree->reoptnodes[id]->conss[reopttree->reoptnodes[id]->nconss] = consdata;
      ++reopttree->reoptnodes[id]->nconss;

      /* add id as a child of the root node */
      SCIP_CALL( reoptAddChild(reopttree, blkmem, 0, id) );

      ++(*ncreatedchilds);
   }
   else
   {
      int c;

      assert(*ncreatedchilds == 1);
      assert(vars != NULL);
      assert(vals != NULL);
      assert(nvars > 0);

      /* create nvars nodes in the fashion of interdiction branching */
      for( c = 0; c < nvars; c++ )
      {
         /* ensure that two free slots are available  */
         SCIP_CALL( reopttreeCheckMemory(reopttree, blkmem) );
         id = (unsigned int) (size_t) SCIPqueueRemove(reopttree->openids);

         assert(0 < id && id < reopt->reopttree->reoptnodessize);
         assert(reopttree->reoptnodes[id] == NULL || reopttree->reoptnodes[id]->nvars == 0);

         /*   1. create the node
          *   2. fix the first v bound changes to vals[v] and v+1 to 1-vals[v]
          *   4. add the ID id as a child of the root node
          */
         SCIP_CALL( createReoptnode(reopttree, set, blkmem, id) );
         reopttree->reoptnodes[id]->parentID = 0;
         reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_TRANSIT;

         /* check memory */
         SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], blkmem, c+1, 0, 0) );
         assert(reopttree->reoptnodes[id]->varssize >= c+1);
         assert(reopttree->reoptnodes[id]->nvars == 0);
         assert(reopttree->reoptnodes[id]->vars != NULL);
         assert(reopttree->reoptnodes[id]->varbounds != NULL);
         assert(reopttree->reoptnodes[id]->varboundtypes != NULL);

         /* copy first v bound changes */
         for( v = 0; v < c; v++ )
         {
            reopttree->reoptnodes[id]->vars[v] = vars[v];
            reopttree->reoptnodes[id]->varbounds[v] = vals[v];
            reopttree->reoptnodes[id]->varboundtypes[v] = SCIPsetIsFeasEQ(set, vals[v], 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
            ++reopttree->reoptnodes[id]->nvars;
         }

         /* set bound change v+1 (= c) to 1-vals[c] */
         assert(v == c);
         reopttree->reoptnodes[id]->vars[c] = vars[c];
         reopttree->reoptnodes[id]->varbounds[c] = 1-vals[c];
         reopttree->reoptnodes[id]->varboundtypes[c] = SCIPsetIsFeasEQ(set, 1-vals[c], 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
         ++reopttree->reoptnodes[id]->nvars;

         /* add dummy1 as a child of the root node */
         SCIP_CALL( reoptAddChild(reopttree, blkmem, 0, id) );

         ++(*ncreatedchilds);
      }

      assert(*ncreatedchilds == nvars+1);
   }

   /* free the current dualconscur and assign dualconsnex */
   assert(reopttree->reoptnodes[0]->dualconscur->vars != NULL);
   assert(reopttree->reoptnodes[0]->dualconscur->vals != NULL);

   /* free the current dualconscur and assign dualconsnex */
   SCIP_CALL( reoptnodeUpdateDualConss(reopttree->reoptnodes[0], blkmem) );

   /* change the reopttype of the root node */
   SCIPnodeSetReopttype(SCIPtreeGetRootNode(tree), SCIP_REOPTTYPE_TRANSIT);

   return SCIP_OKAY;
}

/** reset the stored information abound bound changes based on dual information */
SCIP_RETCODE SCIPreoptResetDualBndchgs(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node,               /**< node of the search tree */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   unsigned int id;

   assert(reopt != NULL);
   assert(node != NULL);

   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* return if the node ist not part of the reoptimization tree */
   if( SCIPnodeGetDepth(node) > 0 && id == 0 )
      return SCIP_OKAY;

   /* reset the dual constraint */
   SCIP_CALL( reoptnodeResetDualConss(reopt->reopttree->reoptnodes[id], blkmem) );

   return SCIP_OKAY;
}

/** return the branching path stored of the given node in the reoptimization tree */
void SCIPreoptnodeGetPath(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR**            vars,               /**< array for variables */
   SCIP_Real*            vals,               /**< array for values */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array for bound types */
   int                   varssize,           /**< size of arrays vars, vals, and boundtypes */
   int*                  nbndchgs,           /**< pointer to store the number of bound changes */
   int*                  nbndchgsafterdual   /**< pointer to store the number of bound changes applied after
                                              *  the first dual reduction at the given node */
   )
{
   int v;
   int nvars2;
   int nafterdualvars2;

   assert(reopt != NULL);
   assert(reoptnode != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(boundtypes != NULL);

   (*nbndchgs) = reoptnode->nvars;
   (*nbndchgsafterdual) = reoptnode->nafterdualvars;

   if( varssize == 0 || varssize < *nbndchgs + *nbndchgsafterdual )
      return;

   for(v = 0; v < *nbndchgs; v++)
   {
      vars[v] = reoptnode->vars[v];
      vals[v] = reoptnode->varbounds[v];
      boundtypes[v] = reoptnode->varboundtypes[v];
   }

   for(; v < *nbndchgs + *nbndchgsafterdual; v++)
   {
      vars[v] = reoptnode->afterdualvars[v-(*nbndchgs)];
      vals[v] = reoptnode->afterdualvarbounds[v-(*nbndchgs)];
      boundtypes[v] = reoptnode->afterdualvarboundtypes[v-(*nbndchgs)];
   }

   if( reoptnode->parentID != 0 )
   {
      SCIP_REOPTNODE* parent;

      parent = reopt->reopttree->reoptnodes[reoptnode->parentID];
      SCIPreoptnodeGetPath(reopt, parent, &vars[v], &vals[v], &boundtypes[v], varssize, &nvars2, &nafterdualvars2);

      (*nbndchgs) += nvars2;
      (*nbndchgsafterdual) += nafterdualvars2;
   }

   return;
}

/** delete a node stored in the reoptimization tree */
SCIP_RETCODE SCIPreoptDeleteNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   unsigned int          id,                 /**< id of a stored node */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);
   assert(reopt->reopttree->reoptnodes[id] != NULL);
   assert(blkmem != NULL);

   SCIP_CALL( reopttreeDeleteNode(reopt->reopttree, set, blkmem, id, TRUE) );

   return SCIP_OKAY;
}

/** reactivate the given @p reoptnode and split them into several nodes if necessary */
SCIP_RETCODE SCIPreoptApply(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branching tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          randseed,           /**< seed value for random generator */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree to reactivate */
   unsigned int          id,                 /**< id of the node to reactivate */
   SCIP_Real             estimate,           /**< estimate of the child nodes that should be created */
   SCIP_NODE**           childnodes,         /**< array to store the created child nodes */
   int*                  ncreatedchilds,     /**< pointer to store number of created child nodes */
   int*                  naddedconss,        /**< pointer to store number of generated constraints */
   int                   childnodessize,     /**< available size of childnodes array */
   SCIP_Bool*            success             /**< pointer store the result */
   )
{
   assert(reopt != NULL);
   assert(scip != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(cliquetable != NULL);
   assert(blkmem != NULL);
   assert(reoptnode != NULL);
   assert(childnodes != NULL);
   assert(reopt->reopttree != NULL);
   assert(id < reopt->reopttree->reoptnodessize);

   SCIPdebugMessage("reactivating node at id %u:\n", id);

   *success = FALSE;

   /* check if we need to split the node */
   if( reoptnode->reopttype == (unsigned int)SCIP_REOPTTYPE_STRBRANCHED
    || reoptnode->reopttype == (unsigned int)SCIP_REOPTTYPE_INFSUBTREE )
   {
      int c;

      assert(reoptnode->dualfixing);

      /* we want use a constraint to split the node into two disjoint node */
      if( set->reopt_usesplitcons )
      {
         if( reoptnode->reopttype == (unsigned int)SCIP_REOPTTYPE_INFSUBTREE )
         {
            assert(reoptnode->dualconscur != NULL);
            assert(reoptnode->dualconscur->constype == REOPT_CONSTYPE_INFSUBTREE);
            (*ncreatedchilds) = 1;
         }
         else
         {
            assert(reoptnode->dualconscur != NULL);
            assert(reoptnode->dualconscur->constype == REOPT_CONSTYPE_STRBRANCHED);
            (*ncreatedchilds) = 2;
         }

         /* in both cases we add exactly one constraint */
         (*naddedconss) = 1;

         if( childnodessize < *ncreatedchilds )
            return SCIP_OKAY;

         /* generate the nodes */
         for( c = 0; c < *ncreatedchilds; c++ )
         {
            /* create the child node */
            SCIP_CALL( SCIPnodeCreateChild(&childnodes[c], blkmem, set, stat, tree, 1.0, estimate) );

            /* change all bounds; convert the bound changes after the first based on dual reductions into branching
             * for second node only. if we generate only one node, i.e., the pruned part, we do not need this
             * changes anyway.
             */
            SCIP_CALL( changeAncestorBranchings(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue,
                  cliquetable, blkmem, childnodes[c], id, c == 1) );

            /* add all local constraints */
            SCIP_CALL( addLocalConss(scip, reopt, set, stat, blkmem, childnodes[c], id) );

            if( c == 0 )
            {
               /* in both cases the node generated first represents the pruned is currently not part of the reoptimization tree */
               SCIPnodeSetReopttype(childnodes[c], SCIP_REOPTTYPE_NONE);

               /* add the constraint to the node */
               assert(reopt->reopttree->reoptnodes[id]->dualconscur != NULL);
               SCIP_CALL( addSplitcons(reopt, scip, set, stat, blkmem, transprob, origprob, tree, lp, branchcand,
                     eventqueue, cliquetable, childnodes[c], id) );

               /* fixBounds() does the same, but in this case we go not into it */
               if( reoptnode->dualconscur->constype == REOPT_CONSTYPE_INFSUBTREE )
               {
                  assert(reoptnode->dualconscur->nvars > 0);
                  assert(reoptnode->dualconscur->varssize > 0);

                  /* delete dualconscur and move dualconsnex -> dualconscur */
                  SCIP_CALL( reoptnodeUpdateDualConss(reoptnode, blkmem) );
               }
            }
            else
            {
               /* if we reach this lines of code, the current node represents the original node including all bound
                * changes based in dual information.
                */
               assert(reoptnode->dualconscur->constype == REOPT_CONSTYPE_STRBRANCHED);
               if( reoptnode->nconss == 0 )
                  SCIPnodeSetReopttype(childnodes[c], SCIP_REOPTTYPE_TRANSIT);
               else
                  SCIPnodeSetReopttype(childnodes[c], SCIP_REOPTTYPE_LOGICORNODE);

               /* fix all bound changes based on dual information and convert them into branchings */
               assert(reopt->reopttree->reoptnodes[id]->dualconscur != NULL);
               SCIP_CALL( fixBounds(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue, cliquetable,
                     blkmem, childnodes[c], id, TRUE) );

               /* set the unique id the id of the original node */
               SCIPnodeSetReoptID(childnodes[c], id);
            }

            /* set the estimate */
            if( !SCIPsetIsInfinity(set, REALABS(reoptnode->lowerbound)) )
            {
               if( SCIPsetIsRelGE(set, reoptnode->lowerbound, SCIPnodeGetLowerbound(childnodes[c])) )
                  SCIPnodeSetEstimate(childnodes[c], set, reoptnode->lowerbound);
            }
         }

         /* reset the stored dual constraints */
         SCIP_CALL( reoptnodeUpdateDualConss(reopt->reopttree->reoptnodes[id], blkmem) );

         /* set the reoptimization type */
         if( reopt->reopttree->reoptnodes[id]->dualfixing )
            reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_STRBRANCHED;
         else
            reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_TRANSIT;

         *success = TRUE;
      }
      else
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int nvars;

         vars = reoptnode->dualconscur->vars;
         vals = reoptnode->dualconscur->vals;
         nvars = reoptnode->dualconscur->nvars;

         /* calculate the order of the variables */
         switch (set->reopt_varorderinterdiction)
         {
         case 'd':
            break;

         case 'r':
            permuteRandom(vars, vals, nvars, &randseed);
            break;

         default:
            return SCIP_INVALIDDATA;
         }

         *ncreatedchilds = nvars+1;
         *naddedconss = 0;

         if( childnodessize < *ncreatedchilds )
            return SCIP_OKAY;

         assert(reopt->reopttree->reoptnodes[id] != NULL);
         reoptnode = reopt->reopttree->reoptnodes[id];

         /* enough that the node need to split */
         assert(reoptnode->dualfixing);

         /* iterate over all nodes and change the necessary bounds (nodes[0] corresponds to the original one)
          * we need to do this in the reverse order because we want to transform the bound changes based on dual information
          * into branching decisions at nodes[0].
          */
         for( c = nvars; c >= 0; c-- )
         {
            /* create the child node */
            SCIP_CALL( SCIPnodeCreateChild(&childnodes[c], blkmem, set, stat, tree, 1.0, estimate) );

   #ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage(" change bounds at node %lld\n", SCIPnodeGetNumber(childnodes[c]));
   #endif

            /* change all bounds */
            SCIP_CALL( changeAncestorBranchings(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue,
                  cliquetable, blkmem, childnodes[c], id, FALSE) );

            /* reconstruct the original node and the pruned part, respectively */
            if( c == 0 )
            {
               /* fix bound changes based on dual information and convert all these bound changes to normal bound changes */
               SCIP_CALL( fixBounds(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue, cliquetable,
                     blkmem, childnodes[c], id, TRUE) );

               /* set the reopttype of the node */
               SCIPnodeSetReopttype(childnodes[c], SCIP_REOPTTYPE_TRANSIT);

               /* set the unique id */
               SCIPnodeSetReoptID(childnodes[c], id);
            }
            else
            {
               /* fix the first c bound changes and negate the (c+1)th */
               SCIP_CALL( fixInterdiction(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue, cliquetable,
                     blkmem, childnodes[c], id, vars, vals, nvars, c) );
            }

            /* add all local constraints */
            SCIP_CALL( addLocalConss(scip, reopt, set, stat, blkmem, childnodes[c], id) );

            /* set estimates */
            if( !SCIPsetIsInfinity(set, REALABS(reopt->reopttree->reoptnodes[id]->lowerbound)) )
            {
               if( SCIPsetIsRelGE(set, reoptnode->lowerbound, SCIPnodeGetLowerbound(childnodes[c])))
                  SCIPnodeSetEstimate(childnodes[c], set, reoptnode->lowerbound);
            }
         }

         /* reset the stored dual constraints */
         SCIP_CALL( reoptnodeUpdateDualConss(reopt->reopttree->reoptnodes[id], blkmem) );

         /* set the reoptimization type to transit */
         if( reopt->reopttree->reoptnodes[id]->dualfixing )
            reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_STRBRANCHED;
         else
            reopt->reopttree->reoptnodes[id]->reopttype = (unsigned int)SCIP_REOPTTYPE_TRANSIT;

         *success = TRUE;
      }
   }
   else
   {
      /* we need the create exactly one node to reconstruct the node itself and no additional constraint */
      (*ncreatedchilds) = 1;
      (*naddedconss) = 0;

      if( childnodessize < *ncreatedchilds )
         return SCIP_OKAY;

      /* create the child node */
      SCIP_CALL( SCIPnodeCreateChild(&childnodes[0], blkmem, set, stat, tree, 1.0, estimate) );

      /* change all bounds */
      assert(reoptnode->nafterdualvars == 0);
      SCIP_CALL( changeAncestorBranchings(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue,
            cliquetable, blkmem, childnodes[0], id, FALSE) );

      /* add all local constraints */
      SCIP_CALL( addLocalConss(scip, reopt, set, stat, blkmem, childnodes[0], id) );

      /* set the estimate */
      if( !SCIPsetIsInfinity(set, REALABS(reopt->reopttree->reoptnodes[id]->lowerbound)) )
      {
         if( SCIPsetIsRelGE(set, reopt->reopttree->reoptnodes[id]->lowerbound, SCIPnodeGetLowerbound(childnodes[0])) )
            SCIPnodeSetEstimate(childnodes[0], set, reopt->reopttree->reoptnodes[id]->lowerbound);
      }

      /* set the reopttype */
      assert(reoptnode->reopttype != (unsigned int)SCIP_REOPTTYPE_INFSUBTREE
          && reoptnode->reopttype != (unsigned int)SCIP_REOPTTYPE_INFSUBTREE);
      SCIPnodeSetReopttype(childnodes[0], (SCIP_REOPTTYPE)reoptnode->reopttype);

      /* set the unique id */
      SCIPnodeSetReoptID(childnodes[0], id);

      *success = TRUE;
   }

   return SCIP_OKAY;
}

#ifdef SCIP_DISABLED_CODE
/** Reoptimize the node stored at ID @p id in the fashion of interdiction branching,
 *  i.e. create and split the node in the current run, if necessary.
 *
 *  To reconstruct the pruned part we create @p nnodes nodes, whereby
 *  - nodes[0] corresponds to the original node
 *  - nodes[k] contains: var[0] = ... = var[k-1] = 0 and var[k] = 1
 *  where var are the (negated) variables fixed to 0 by dual reductions.
 */
SCIP_RETCODE SCIPreoptApplyInterdiction(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branching tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   SCIP_NODE**           nodes,              /**< array to store created nodes */
   int                   nnodes,             /**< size of the array */
   int                   id,                 /**< id of a stored node which should be reoptimized */
   int*                  permutation,        /**< permutation of the variable order (within the constraint) */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_REOPTNODE* reoptnode;
   int c;

   assert(reopt != NULL);
   assert(scip != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(branchcand != NULL);
   assert(eventqueue != NULL);
   assert(cliquetable != NULL);
   assert(nodes != NULL || nnodes == 0);
   assert(blkmem != NULL);

   SCIPdebugMessage("reoptimizing node at ID %d:\n", id);

   assert(reopt->reopttree->reoptnodes[id] != NULL);
   reoptnode = reopt->reopttree->reoptnodes[id];

   /* enough that the node need to split */
   assert(reoptnode->dualfixing);

   /* iterate over all nodes and change the necessary bounds (nodes[0] corresponds to the original one)
    * we need to do this in the reverse order because we want to transform the bound changes based on dual information
    * into branching decisions at nodes[0].
    */
   for( c = nnodes-1; c >= 0; c-- )
   {
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage(" change bounds at node %lld\n", SCIPnodeGetNumber(nodes[c]));
#endif

      /* change all bounds */
      SCIP_CALL( changeAncestorBranchings(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue,
            cliquetable, blkmem, nodes[c], NULL, id) );

      /* reconstruct the original node and the pruned part, respectively */
      if( c == 0 )
      {
         /* fix bound changes based on dual information and convert all these bound changes to normal bound changes */
         SCIP_CALL( fixBounds(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue, cliquetable,
               blkmem, nodes[c], id, FALSE) );
      }
      else
      {
         /* fix the first c bound changes and negate the (c+1)th */
         SCIP_CALL( fixInterdiction(reopt, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue, cliquetable,
               blkmem, nodes[c], id, permutation, c) );
      }

      /* add all local constraints to both nodes */
      SCIP_CALL( addLocalConss(scip, reopt, set, stat, blkmem, nodes[c], NULL, id) );

      /* set estimates */
      if( !SCIPsetIsInfinity(set, REALABS(reopt->reopttree->reoptnodes[id]->lowerbound)) )
      {
         if( SCIPsetIsRelGE(set, reopt->reopttree->reoptnodes[id]->lowerbound, SCIPnodeGetLowerbound(nodes[c])))
            SCIPnodeSetEstimate(nodes[c], set, reopt->reopttree->reoptnodes[id]->lowerbound);
      }
   }

   /* reset the stored dual constraints */
   SCIP_CALL( reoptnodeUpdateDualConss(reopt->reopttree->reoptnodes[id], blkmem) );

   return SCIP_OKAY;
}
#endif

/** returns the time needed to store the nodes for reoptimization */
SCIP_Real SCIPreoptGetSavingtime(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return SCIPclockGetTime(reopt->savingtime);
}

/** store a global constraint that should be added at the beginning of the next iteration */
SCIP_RETCODE SCIPreoptAddGlbCons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_VAR**            vars,               /**< array to store the variables of the constraint */
   SCIP_Real*            vals,               /**< array to store the coefficients of the variables */
   int                   nvars,              /**< pointer to store the size of the constraints */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(blkmem != NULL);

   if( nvars > 0 )
   {
      int pos;

      /* check the memory */
      SCIP_CALL( checkMemGlbCons(reopt, blkmem, reopt->nglbconss + 1) );
      assert(reopt->allocmemglbconss >= reopt->nglbconss+1);

      pos = reopt->nglbconss;

      /* allocate memory */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopt->glbconss[pos]) ); /*lint !e866*/
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->glbconss[pos]->vars, &vars, nvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->glbconss[pos]->vals, &vals, nvars) );
      reopt->glbconss[pos]->varssize = nvars;
      reopt->glbconss[pos]->nvars = nvars;

      ++reopt->nglbconss;
   }

   return SCIP_OKAY;
}

/** add the stored constraints globally to the problem */
SCIP_RETCODE SCIPreoptApplyGlbConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int c;

   assert(scip != NULL);
   assert(reopt != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);

   if( reopt->glbconss == NULL || reopt->nglbconss == 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("try to add %d glb constraints\n", reopt->nglbconss);

   for(c = 0; c < reopt->nglbconss; c++)
   {
      SCIP_CONS* cons;
      SCIP_VAR** consvars;
      int v;

      assert(reopt->glbconss[c]->nvars > 0);

      /* allocate buffer */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, reopt->glbconss[c]->nvars) );

      SCIPdebugMessage("-> add constraints with %d vars\n", reopt->glbconss[c]->nvars);

      for(v = 0; v < reopt->glbconss[c]->nvars; v++)
      {
         consvars[v] = SCIPvarGetTransVar(reopt->glbconss[c]->vars[v]);

         /* negate the variable if it was fixed to 1 */
         if( SCIPsetIsFeasEQ(set, reopt->glbconss[c]->vals[v], 1.0) )
         {
            SCIP_CALL( SCIPvarNegate(consvars[v], blkmem, set, stat, &consvars[v]) );
         }
      }

      /* create the logic-or constraint and add them to the problem */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, "glblogicor", reopt->glbconss[c]->nvars,
            consvars, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE) );

      SCIPdebugPrintCons(scip, cons, NULL);

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* delete the global constraints data */
      SCIPfreeBlockMemoryArrayNull(scip, &reopt->glbconss[c]->vals, reopt->glbconss[c]->nvars);
      SCIPfreeBlockMemoryArrayNull(scip, &reopt->glbconss[c]->vars, reopt->glbconss[c]->nvars);
      SCIPfreeBlockMemoryNull(scip, &reopt->glbconss[c]); /*lint !e866*/
      reopt->glbconss[c]->nvars = 0;

      /* free buffer */
      SCIPfreeBufferArray(scip, &consvars);
   }

   /* reset the number of global constraints */
#ifdef SCIP_DEBUG
   for(c = 0; c < reopt->nglbconss; c++)
   {
      assert(reopt->glbconss[c]->nvars == 0);
      assert(reopt->glbconss[c]->vars == NULL);
      assert(reopt->glbconss[c]->vals == NULL);
   }
#endif
   reopt->nglbconss = 0;

   return SCIP_OKAY;
}

/** check if the LP of the given node should be solved or not */
SCIP_Bool SCIPreoptGetSolveLP(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node of the current search tree */
   )
{
   unsigned int id;

   assert(reopt != NULL);
   assert(node != NULL);

   /* get the ID */
   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* return if the node is not part of the reoptimization tree */
   if( SCIPnodeGetDepth(node) > 0 && id == 0 )
      return TRUE;

   /* current node is the root */
   if( id == 0 )
   {
      if( reopt->reopttree->reoptnodes[0]->nchilds > 0 )
      {
         /* the objective function has changed only slightly */
         if( reopt->simtolastobj >= set->reopt_objsimrootlp )
            return FALSE;
      }
   }
   else
   {
      /* solve node LP if the node type is greater or equal to solvelp or there were too many bound changes at the current node */
      if( reopt->reopttree->reoptnodes[id]->nvars < set->reopt_solvelpdiff && (int) SCIPnodeGetReopttype(node) < set->reopt_solvelp )
      {
         assert(reopt->reopttree->reoptnodes[id]->nchilds > 0);
         return FALSE;
      }
   }

   return TRUE;
}

/** initialize an empty node */
void SCIPreoptnodeInit(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(reoptnode != NULL);
   assert(set != NULL);

   reoptnode->conss = NULL;
   reoptnode->nconss = 0;
   reoptnode->consssize = 0;
   reoptnode->childids = NULL;
   reoptnode->allocchildmem = 0;
   reoptnode->nchilds = 0;
   reoptnode->nvars = 0;
   reoptnode->nafterdualvars = 0;
   reoptnode->parentID = 0;
   reoptnode->dualfixing = FALSE;
   reoptnode->reopttype = (unsigned int)SCIP_REOPTTYPE_NONE;
   reoptnode->varssize = 0;
   reoptnode->afterdualvarssize = 0;
   reoptnode->vars = NULL;
   reoptnode->varbounds = NULL;
   reoptnode->varboundtypes = NULL;
   reoptnode->afterdualvars = NULL;
   reoptnode->afterdualvarbounds = NULL;
   reoptnode->afterdualvarboundtypes = NULL;
   reoptnode->dualconscur = NULL;
   reoptnode->dualconsnex = NULL;
   reoptnode->lowerbound = -SCIPsetInfinity(set);
}

/** reset the given reoptimization node */
SCIP_RETCODE SCIPreoptnodeReset(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPTNODE*       reoptnode           /**< reoptimization node */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(reoptnode != NULL);

   SCIP_CALL( reoptnodeReset(reoptnode, set, blkmem) );

   return SCIP_OKAY;
}

/** delete the given reoptimization node */
SCIP_RETCODE SCIPreoptnodeDelete(
   SCIP_REOPTNODE**      reoptnode,          /**< pointer of reoptnode */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(reoptnode != NULL);
   assert(blkmem != NULL);

   SCIP_CALL( reoptnodeDelete(reoptnode, blkmem) );

   return SCIP_OKAY;
}

/** add a variable to a given reoptnode */
SCIP_RETCODE SCIPreoptnodeAddBndchg(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             val,                /**< value of the variable */
   SCIP_BOUNDTYPE        boundtype           /**< boundtype of the variable */
   )
{
   int nvars;

   assert(reoptnode != NULL);
   assert(var != NULL);
   assert(blkmem != NULL);

   nvars = reoptnode->nvars;

   SCIP_CALL( reoptnodeCheckMemory(reoptnode, blkmem, nvars + 1, 0, 0) );

   reoptnode->vars[nvars] = var;
   reoptnode->varbounds[nvars] = val;
   reoptnode->varboundtypes[nvars] = boundtype;
   ++reoptnode->nvars;

   return SCIP_OKAY;
}

/** add a constraint to a given reoptnode */
SCIP_RETCODE SCIPreoptnodeAddCons(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            consvars,           /**< variables which are part of the constraint */
   SCIP_Real*            consvals,           /**< values of the variables */
   int                   nvars,              /**< number of variables */
   REOPT_CONSTYPE        constype            /**< type of the constraint */
   )
{
   int nconss;

   assert(reoptnode != NULL);
   assert(consvars != NULL);
   assert(consvals != NULL);
   assert(nvars > 0);
   assert(blkmem != NULL);

   /* the constraint can be interpreted as a normal bound change */
   if( nvars == 1 )
   {
      SCIPdebugMessage("-> constraint has size 1 -> save as normal bound change.\n");

      SCIP_CALL( SCIPreoptnodeAddBndchg(reoptnode, blkmem, consvars[0], 1-consvals[0],
            1-consvals[0] == 1 ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER) );
   }
   else
   {
      nconss = reoptnode->nconss;

      SCIP_CALL( reoptnodeCheckMemory(reoptnode, blkmem, 0, 0, nconss+1) );

      /* create the constraint */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reoptnode->conss[nconss]) ); /*lint !e866*/
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptnode->conss[nconss]->vars, consvars, nvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptnode->conss[nconss]->vals, consvals, nvars) );

      reoptnode->conss[nconss]->varssize = nvars;
      reoptnode->conss[nconss]->nvars = nvars;
      reoptnode->conss[nconss]->constype = constype;
      ++reoptnode->nconss;
   }
   return SCIP_OKAY;
}
