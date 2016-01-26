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
#include "scip/lp.h"
#include "scip/misc.h"
#include "scip/reopt.h"
#include "scip/tree.h"
#include "scip/primal.h"
#include "scip/sepastore.h"
#include "scip/prob.h"
#include "scip/cons.h"
#include "scip/cons_logicor.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/clock.h"
#include "scip/heur_reoptsols.h"
#include "scip/history.h"
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
   assert(SCIPvarGetType(SCIPeventGetVar(event)) != SCIP_VARTYPE_CONTINUOUS);

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
      if( SCIPvarGetType(vars[varnr]) != SCIP_VARTYPE_CONTINUOUS )
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
   SCIP_SOLNODE* sibling;
   int nsols;

   assert(solnode != NULL);

   if( solnode->child == NULL && solnode->sol == NULL )
      return 0;
   if( solnode->child == NULL && solnode->sol != NULL )
      return 1;

   nsols = 0;
   sibling = solnode->child;

   /* travers through the list */
   while( sibling != NULL )
   {
      nsols += soltreeNInducedSols(sibling);
      sibling = sibling->sibling;
   }

   return nsols;
}

/** returns the similarity of the objective functions of two given iterations */
static
SCIP_Real reoptSimilarity(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   obj1_id,            /**< id of one objective function */
   int                   obj2_id,            /**< id of the other objective function */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars               /**< number of problem variables */
   )
{
   SCIP_Real similarity;
   SCIP_Real norm_obj1;
   SCIP_Real norm_obj2;
   int v;

   assert(reopt != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);

   similarity = 0.0;
   norm_obj1 = 0.0;
   norm_obj2 = 0.0;

   /* calc similarity */
   for(v = 0; v < nvars; v++)
   {
      SCIP_VAR* origvar;
      SCIP_VAR* transvar;
      SCIP_Real c1;
      SCIP_Real c2;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real constant;
      SCIP_Real scalar;

      origvar = vars[v];

      /* get the original variable */
      if( !SCIPvarIsOriginal(origvar) )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );
      }
      assert(origvar != NULL && SCIPvarIsOriginal(origvar));

      /* get the transformed variable, we skip globally fixed variables */
      transvar = SCIPvarGetTransVar(origvar);
      assert(transvar != NULL);

      lb = SCIPvarGetLbLocal(transvar);
      ub = SCIPvarGetUbLocal(transvar);

      if( SCIPsetIsFeasLT(set, lb, ub) )
      {
         int probidx;

         probidx = SCIPvarGetProbindex(origvar);
         assert(0 <= probidx && probidx <= nvars);

         c1 = reopt->objs[obj1_id][probidx];
         c2 = reopt->objs[obj2_id][probidx];

         /* vector product */
         similarity += c1*c2;
         norm_obj1 += SQR(c1);
         norm_obj2 += SQR(c2);
      }
   }

   /* divide similarity by norms of the objective vectors */
   norm_obj1 = SQRT(norm_obj1);
   norm_obj2 = SQRT(norm_obj2);

   if( !SCIPsetIsZero(set, norm_obj1) && !SCIPsetIsZero(set, norm_obj2) )
      similarity /= (norm_obj1 * norm_obj2);
   else
      similarity = 0.0;

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
         BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->conss[c]->boundtypes, (*reoptnode)->conss[c]->varssize);
         BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->conss[c]->bounds, (*reoptnode)->conss[c]->varssize);
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
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconscur->boundtypes, (*reoptnode)->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconscur->bounds, (*reoptnode)->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconscur->vars, (*reoptnode)->dualconscur->varssize);
      BMSfreeBlockMemory(blkmem, &(*reoptnode)->dualconscur);
      (*reoptnode)->dualconscur = NULL;
   }

   if( (*reoptnode)->dualconsnex != NULL )
   {
      assert((*reoptnode)->dualconsnex->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconsnex->boundtypes, (*reoptnode)->dualconsnex->varssize);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconsnex->bounds, (*reoptnode)->dualconsnex->varssize);
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
         BMSfreeBlockMemoryArray(blkmem, &reoptnode->conss[c]->bounds, reoptnode->conss[c]->varssize);
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
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->boundtypes, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->bounds, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem ,&reoptnode->dualconscur->vars, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualconscur);
      reoptnode->dualconscur = NULL;
   }

   if( reoptnode->dualconsnex != NULL )
   {
      assert(reoptnode->dualconsnex->varssize > 0);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconsnex->boundtypes, reoptnode->dualconsnex->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconsnex->bounds, reoptnode->dualconsnex->varssize);
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
   soltree->root->value = SCIP_INVALID;
   soltree->root->updated = FALSE;
   soltree->root->father = NULL;
   soltree->root->child = NULL;
   soltree->root->sibling = NULL;

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
   SCIP_SOLNODE* child;
   SCIP_SOLNODE* sibling;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(primal != NULL || set->stage == SCIP_STAGE_INIT);
   assert(solnode != NULL);
   assert(blkmem != NULL);

   child = (*solnode)->child;

   /* travers through the list and free recursive all subtree */
   while( child != NULL )
   {
      SCIP_CALL( soltreefreeNode(reopt, set, primal, blkmem, &child) );
      assert(child != NULL);

      sibling = child->sibling;
      BMSfreeBlockMemoryNull(blkmem, &child);
      child = sibling;
   }

   if( (*solnode)->sol != NULL )
   {
      assert(set->stage == SCIP_STAGE_PROBLEM);

      SCIP_CALL( SCIPsolFree(&(*solnode)->sol, blkmem, primal) );
   }

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
   BMSfreeBlockMemoryNull(blkmem, &reopt->soltree->root);

   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->sols, reopt->runsize);
   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->nsols, reopt->runsize);
   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->solssize, reopt->runsize);

   BMSfreeMemory(&reopt->soltree);

   return SCIP_OKAY;
}

/** creates and adds a solution node to the solution tree */
static
SCIP_RETCODE solnodeAddChild(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SOLNODE*         curnode,            /**< current node in the solution tree */
   SCIP_SOLNODE**        child,              /**< pointer to store the node representing the solution value */
   SCIP_Real             val,                /**< value the child shell represent */
   SCIP_Bool*            added               /**< TRUE iff we created a new node, i.e, we have not seen this solution so far */
   )
{
   SCIP_SOLNODE* solnode;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(curnode != NULL);
   assert(child != NULL && *child == NULL);
   assert(!SCIPsetIsInfinity(set, -val) && !SCIPsetIsInfinity(set, val));

   /* get the first node of the child node list */
   *child = curnode->child;

   /* this is the first solution in the subtree induced by the current node */
   if( *child == NULL )
   {
      assert(soltreeNInducedSols(curnode) == 0);

      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &solnode) );
      solnode->sol = NULL;
      solnode->updated = FALSE;
      solnode->father = curnode;
      solnode->child = NULL;
      solnode->sibling = NULL;
      solnode->value = val;

      *added = TRUE;
      *child = solnode;

      curnode->child = *child;

      SCIPdebugMessage("-> create new node %p: value=%g, sibling=%p\n", (void*) solnode, solnode->value,
            (void*) solnode->sibling);
   }
   else
   {
      /* we traverse through all children */
      while( *child != NULL )
      {
         SCIPdebugMessage("-> check %p: father=%p, value=%g, sibling=%p\n", (void*) *child, (void*) (*child)->father,
               (*child)->value, (void*) (*child)->sibling);

         /* we found a node repesenting this solution value */
         if( SCIPsetIsEQ(set, val, (*child)->value) )
            break;

         /* we are at the end of the list */
         if( (*child)->sibling == NULL )
         {
            /* create a new solnode */
            SCIP_ALLOC( BMSallocBlockMemory(blkmem, &solnode) );
            solnode->sol = NULL;
            solnode->updated = FALSE;
            solnode->father = curnode;
            solnode->child = NULL;
            solnode->value = val;
            *added = TRUE;

            /* we have to append the new node at the end of the list. but we have to check whether the insertion before
             * the current node would be correct. in that case, we switch the values, the child pointer, and the solution */
            solnode->sibling = NULL;
            (*child)->sibling = solnode;

            SCIPdebugMessage("-> create new node %p: value=%g, sibling=%p\n", (void*) solnode, solnode->value,
                  (void*) solnode->sibling);

            /* the given value is lower than the current, insertion before the current node would be correct
             * in this case we do not have to change the child pointer
             */
            if( SCIPsetIsLT(set, val, (*child)->value) )
            {
               SCIPdebugMessage("-> need to switch:");
               SCIPdebugMessage("   before switching: node %p witch child=%p, sibling=%p, sol=%p, value=%g\n",
                     (void*) (*child), (void*) (*child)->child, (void*) (*child)->sibling, (void*) (*child)->sol, (*child)->value);
               SCIPdebugMessage("                     node %p witch child=%p, sibling=%p, sol=%p, value=%g\n",
                     (void*) solnode, (void*) solnode->child, (void*) solnode->sibling, (void*) solnode->sol, solnode->value);

               /* switch child pointer */
               solnode->child = (*child)->child;
               (*child)->child = NULL;

               /* switch solution values */
               solnode->value = (*child)->value;
               (*child)->value = val;
               assert(SCIPsetIsLT(set, (*child)->value, solnode->value));

               /* switch solution pointer */
               solnode->sol = (*child)->sol;
               (*child)->sol = NULL;

               SCIPdebugMessage("    after switching: node %p witch child=%p, sibling=%p, sol=%p, value=%g\n",
                     (void*) (*child), (void*) (*child)->child, (void*) (*child)->sibling, (void*) (*child)->sol, (*child)->value);
               SCIPdebugMessage("                     node %p witch child=%p, sibling=%p, sol=%p, value=%g\n",
                     (void*) solnode, (void*) solnode->child, (void*) solnode->sibling, (void*) solnode->sol, solnode->value);
            }
            /* set the child pointer to the new created solnode */
            else
               (*child) = solnode;

            break;
         }

         /* the next sibling represents a solution value of larger size.
          * we insert a new node between the current child and the next sibling.
          */
         if( SCIPsetIsLT(set, val, (*child)->sibling->value) )
         {
            /* create a new solnode that points to the sibling of the current child */
            SCIP_ALLOC( BMSallocBlockMemory(blkmem, &solnode) );
            solnode->sol = NULL;
            solnode->updated = FALSE;
            solnode->father = curnode;
            solnode->child = NULL;
            solnode->sibling = (*child)->sibling;
            solnode->value = val;
            *added = TRUE;

            /* change the poiter of the next sibling to the new node */
            (*child)->sibling = solnode;

            *child = solnode;

            SCIPdebugMessage("-> create new node %p: value=%g, sibling=%p\n", (void*) solnode, solnode->value,
                  (void*) solnode->sibling);

            break;
         }

         /* go to the next sibling */
         *child = (*child)->sibling;
      }

#ifdef SCIP_DEBUG
      /* check whether the insert was correct and the list is increasing */
      solnode = curnode->child;
      assert(solnode != NULL);

      while( solnode->sibling != NULL )
      {
         assert(SCIPsetIsLT(set, solnode->value, solnode->sibling->value));
         solnode = solnode->sibling;
      }
#endif
   }
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
   *added = FALSE;

   if( set->reopt_savesols > 0 )
   {
      SCIPdebugMessage("try to add solution found by <%s>\n", SCIPsolGetHeur(sol) == NULL ? "relaxation" : SCIPheurGetName(SCIPsolGetHeur(sol)));

      for( varid = 0; varid < nvars; varid++ )
      {
         if( SCIPvarGetType(vars[varid]) != SCIP_VARTYPE_CONTINUOUS )
         {
            SCIP_SOLNODE* child;

            child = NULL;
            SCIP_CALL( solnodeAddChild(set, blkmem, cursolnode, &child, SCIPsolGetVal(sol, set, stat, vars[varid]), added) );
            assert(child != NULL);
            cursolnode = child;
         }
      }

      /* the solution was added or is an optimal solution */
      if( *added || bestsol )
      {
         SCIP_SOL* copysol;

         assert(cursolnode->child == NULL);

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

   if( node->child != NULL )
   {
      SCIP_SOLNODE* child;

      /* the node is no leaf */
      assert(node->sol == NULL);
      assert(!node->updated);

      child = node->child;

      /* travers through the list of siblings */
      while( child != NULL )
      {
         soltreeResetMarks(child);
         child = child->sibling;
      }
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

   assert(reopttree->nreoptnodes + SCIPqueueNElems(reopttree->openids) == (int)reopttree->reoptnodessize);

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
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->bounds, size) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->boundtypes, size) );
      reopt->dualcons->varssize = size;
      reopt->dualcons->nvars = 0;
   }
   else
   {
      if( reopt->dualcons->varssize > 0 )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->dualcons->vars, reopt->dualcons->varssize, size) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->dualcons->bounds, reopt->dualcons->varssize, size) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->dualcons->boundtypes, reopt->dualcons->varssize, size) );
      }
      else
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->vars, size) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->bounds, size) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->boundtypes, size) );
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

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->boundtypes, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->bounds, reoptnode->dualconscur->varssize);
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

      if( id == 0 )
         reopt->nlocrestarts = 0;

      sim = reoptSimilarity(reopt, set, reopt->run, reopt->run-1, transvars, ntransvars);

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
         SCIP_CONSHDLR* conshdlr;
         SCIP_Real* bounds;
         SCIP_BOUNDTYPE* boundtypes;
         SCIP_Bool success;

         conshdlr = SCIPconsGetHdlr(addedcons[consnr]);
         assert(strcmp(SCIPconshdlrGetName(conshdlr), "bounddisjunction") == 0
             || strcmp(SCIPconshdlrGetName(conshdlr), "logicor") == 0);

         SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopttree->reoptnodes[id]->conss[nconss]) ); /*lint !e866*/

         success = FALSE;
         SCIP_CALL( SCIPconsGetNVars(addedcons[consnr], set, &reopttree->reoptnodes[id]->conss[nconss]->nvars, &success) );
         assert(success);
         reopttree->reoptnodes[id]->conss[nconss]->varssize = reopttree->reoptnodes[id]->conss[nconss]->nvars;

         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->conss[nconss]->vars,
               reopttree->reoptnodes[id]->conss[nconss]->nvars) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->conss[nconss]->bounds,
               reopttree->reoptnodes[id]->conss[nconss]->nvars) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->conss[nconss]->boundtypes,
               reopttree->reoptnodes[id]->conss[nconss]->nvars) );

         success = FALSE;
         SCIP_CALL( SCIPconsGetVars(addedcons[consnr], set, reopttree->reoptnodes[id]->conss[nconss]->vars,
         reopttree->reoptnodes[id]->conss[nconss]->nvars, &success) );
         assert(success);

         /* only needed for bounddisjuction constraints, thus we set them to NULL to avoid compiler warnings */
         bounds = NULL;
         boundtypes = NULL;

         /* get bounds and bound types for bound disjunction constraints */
         if( strcmp(SCIPconshdlrGetName(conshdlr), "bounddisjunction") == 0 )
         {
            bounds = SCIPgetBoundsBounddisjunction(NULL, addedcons[consnr]);
            boundtypes = SCIPgetBoundtypesBounddisjunction(NULL, addedcons[consnr]);
         }

         if( strcmp("infsubtree", SCIPconsGetName(addedcons[consnr])) == 0 )
            reopttree->reoptnodes[id]->conss[nconss]->constype = REOPT_CONSTYPE_INFSUBTREE;
         else if( strcmp("splitcons", SCIPconsGetName(addedcons[consnr])) == 0 )
            reopttree->reoptnodes[id]->conss[nconss]->constype = REOPT_CONSTYPE_DUALREDS;

         assert(reopttree->reoptnodes[id]->conss[nconss]->constype == REOPT_CONSTYPE_INFSUBTREE
             || reopttree->reoptnodes[id]->conss[nconss]->constype == REOPT_CONSTYPE_DUALREDS);

         for(var = 0; var < reopttree->reoptnodes[id]->conss[nconss]->nvars; var++)
         {
            constant = 0;
            scalar = 1;

            if(!SCIPvarIsOriginal(reopttree->reoptnodes[id]->conss[nconss]->vars[var]))
            {
               if(SCIPvarIsNegated(reopttree->reoptnodes[id]->conss[nconss]->vars[var]))
               {
                  SCIP_CALL(SCIPvarGetOrigvarSum(&reopttree->reoptnodes[id]->conss[nconss]->vars[var], &scalar, &constant));

                  if( strcmp(SCIPconshdlrGetName(conshdlr), "logicor") == 0 )
                  {
                     reopttree->reoptnodes[id]->conss[nconss]->bounds[var] = 1;
                     reopttree->reoptnodes[id]->conss[nconss]->boundtypes[var] = SCIP_BOUNDTYPE_LOWER;
                  }
                  else
                  {
                     reopttree->reoptnodes[id]->conss[nconss]->bounds[var] = (bounds[var] - constant) / scalar;
                     reopttree->reoptnodes[id]->conss[nconss]->boundtypes[var] = boundtypes[var];
                  }
               }
               else
               {
                  SCIP_CALL(SCIPvarGetOrigvarSum(&reopttree->reoptnodes[id]->conss[nconss]->vars[var], &scalar, &constant));
                  if( strcmp(SCIPconshdlrGetName(conshdlr), "logicor") == 0 )
                  {
                     reopttree->reoptnodes[id]->conss[nconss]->bounds[var] = 0;
                     reopttree->reoptnodes[id]->conss[nconss]->boundtypes[var] = SCIP_BOUNDTYPE_UPPER;
                  }
                  else
                  {
                     reopttree->reoptnodes[id]->conss[nconss]->bounds[var] = (bounds[var] - constant) / scalar;
                     reopttree->reoptnodes[id]->conss[nconss]->boundtypes[var] = boundtypes[var];
                  }
               }
               assert(reopttree->reoptnodes[id]->conss[nconss]->vars[var] != NULL );
            }
            assert(SCIPvarIsOriginal(reopttree->reoptnodes[id]->conss[nconss]->vars[var]));
         }

         reopttree->reoptnodes[id]->conss[nconss]->lhs = 1.0;
         reopttree->reoptnodes[id]->conss[nconss]->rhs = SCIPsetInfinity(set);

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
   SCIP_SET*             set,                /**< global SCIP settings */
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
            reopt->dualcons->bounds,
            reopt->dualcons->boundtypes,
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
         reopt->dualcons->bounds[v] = (reopt->dualcons->bounds[v] - constant) / scalar;

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
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconscur->bounds,
            reopt->dualcons->bounds, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconscur->boundtypes,
            reopt->dualcons->boundtypes, nbndchgs) );

      reopt->reopttree->reoptnodes[id]->dualconscur->nvars = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconscur->varssize = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconscur->lhs = 1.0;
      reopt->reopttree->reoptnodes[id]->dualconscur->rhs = SCIPsetInfinity(set);
      reopt->reopttree->reoptnodes[id]->dualconscur->constype = reopttype == SCIP_REOPTTYPE_STRBRANCHED ? REOPT_CONSTYPE_DUALREDS : REOPT_CONSTYPE_INFSUBTREE;

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

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconsnex->vars,
            reopt->dualcons->vars, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconsnex->bounds,
            reopt->dualcons->bounds, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconsnex->boundtypes,
            reopt->dualcons->boundtypes, nbndchgs) );
      reopt->reopttree->reoptnodes[id]->dualconsnex->nvars = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconsnex->varssize = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconsnex->lhs = 1.0;
      reopt->reopttree->reoptnodes[id]->dualconsnex->rhs = SCIPsetInfinity(set);
      reopt->reopttree->reoptnodes[id]->dualconsnex->constype = reopttype == SCIP_REOPTTYPE_STRBRANCHED ? REOPT_CONSTYPE_DUALREDS : REOPT_CONSTYPE_INFSUBTREE;

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
         SCIPdebugMessage(" -> nvars: %d, ncons: %d, parentID: %d, reopttype: %d\n",
               reopt->reopttree->reoptnodes[id]->nvars,
               reopt->reopttree->reoptnodes[id]->nconss,
               reopt->reopttree->reoptnodes[id]->parentID, reopttype);
#ifdef SCIP_MORE_DEBUG
         int varnr;
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
      SCIPdebugMessage("save node #%lld successful\n", SCIPnodeGetNumber(node));
      SCIPdebugMessage(" -> ID %d, nvars %d, ncons %d, reopttype %d\n",
            id, reopt->reopttree->reoptnodes[id]->nvars + reopt->reopttree->reoptnodes[id]->nafterdualvars,
            reopt->reopttree->reoptnodes[id]->nconss,
            reopttype);
#ifdef SCIP_MORE_DEBUG
      int varnr;
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
         SCIP_CALL( collectDualInformation(reopt, set, blkmem, node, id, reopttype) );

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

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->boundtypes, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->bounds, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconscur->vars, reoptnode->dualconscur->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualconscur);
      reoptnode->dualconscur = NULL;
   }

   if( reoptnode->dualconsnex != NULL )
   {
      SCIPdebugMessage("reset dual (2) information\n");

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconsnex->boundtypes, reoptnode->dualconsnex->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconsnex->bounds, reoptnode->dualconsnex->varssize);
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->dualconsnex->vars, reoptnode->dualconsnex->varssize);
      BMSfreeBlockMemory(blkmem, &reoptnode->dualconsnex);
      reoptnode->dualconsnex = NULL;
   }

   reoptnode->dualfixing = FALSE;

   return SCIP_OKAY;
}


/** transform given set of variables, bounds and boundtypes into a global cut.
 *  note: boundtypes can be NULL if all variables are binary or a MIP solution should be separated.
 *  note: continouse variables will be skiped if boundtypes is NULL
 */
static
SCIP_RETCODE addGlobalCut(
   SCIP_REOPT*           reopt,
   BMS_BLKMEM*           blkmem,
   SCIP_SET*             set,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   SCIP_BOUNDTYPE*       boundtypes,
   int                   nvars,
   int                   nbinvars,
   int                   nintvars
   )
{
   int nglbconss;
   int nvarsadded;
   int v;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(nbinvars + nintvars == nvars);

   nglbconss = reopt->nglbconss;
   nvarsadded = 0;

   /* check whether we have enough memory allocated */
   SCIP_CALL( checkMemGlbCons(reopt, blkmem, nglbconss+1) );

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reopt->glbconss[nglbconss]) ); /*lint !e866*/
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vars, (int)(nbinvars+2*nintvars)) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->bounds, (int)(nbinvars+2*nintvars)) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->boundtypes, (int)(nbinvars+2*nintvars)) );
   reopt->glbconss[nglbconss]->constype = REOPT_CONSTYPE_CUT;
   reopt->glbconss[nglbconss]->varssize = (int)(nbinvars+2*nintvars);
   reopt->glbconss[nglbconss]->lhs = 1.0;
   reopt->glbconss[nglbconss]->rhs = SCIPsetInfinity(set);
   reopt->glbconss[nglbconss]->nvars = 0;

   for( v = 0; v < nvars; v++ )
   {
      assert(nvarsadded < reopt->glbconss[nglbconss]->varssize);
      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS || SCIPsetIsIntegral(set, vals[v]));

      /* if no boundtypes are given we skip continuous variables, otherwise we would add trivial clauses:
       * a)       x <= ub
       * b) lb <= x
       * c) (x <= val) or (x >= val)
       */
      if( boundtypes == NULL && SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY )
      {
         assert(SCIPvarIsOriginal(vars[v]));

         reopt->glbconss[nglbconss]->vars[nvarsadded] = vars[v];

         if( SCIPsetIsEQ(set, vals[v], 1.0) )
         {
            assert(boundtypes == NULL || boundtypes[v] == SCIP_BOUNDTYPE_LOWER);
            reopt->glbconss[nglbconss]->bounds[nvarsadded] = 0.0;
            reopt->glbconss[nglbconss]->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_UPPER;
         }
         else
         {
            assert(SCIPsetIsEQ(set, vals[v], 0.0));
            assert(boundtypes == NULL || boundtypes[v] == SCIP_BOUNDTYPE_UPPER);
            reopt->glbconss[nglbconss]->bounds[nvarsadded] = 1.0;
            reopt->glbconss[nglbconss]->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_LOWER;
         }
         ++nvarsadded;
      }
      else if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS)
      {
         assert(boundtypes != NULL);

         reopt->glbconss[nglbconss]->bounds[nvarsadded] = vals[v];
         reopt->glbconss[nglbconss]->boundtypes[nvarsadded] = (SCIP_BOUNDTYPE)(1-boundtypes[v]);
         ++nvarsadded;
      }
      else
      {
         SCIP_Real roundedval;
         SCIP_Real ubglb;
         SCIP_Real lbglb;

         assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_IMPLINT);

         reopt->glbconss[nglbconss]->vars[nvarsadded] = vars[v];

         ubglb = SCIPvarGetUbGlobal(vars[v]);
         lbglb = SCIPvarGetLbGlobal(vars[v]);

         /* case 1  :      x == val == ub -> x <= ub-1
          * case 2  :      x == val == lb -> x >= lb+1
          * case 3.1:      x <= val <  ub -> x >= y+1
          * case 3.2:      x >= val >  lb -> x <= y-1
          * case 4  : lb < x == val <  ub -> (x <= y-1) or (x >= y+1)
          */

         /* case 1 */
         if( SCIPsetIsEQ(set, vals[v], ubglb) )
         {
            assert(boundtypes == NULL || boundtypes[v] == SCIP_BOUNDTYPE_LOWER);
            reopt->glbconss[nglbconss]->bounds[nvarsadded] = ubglb - 1.0;
            reopt->glbconss[nglbconss]->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_UPPER;
            ++nvarsadded;
         }
         /* case 2 */
         else if( SCIPsetIsEQ(set, vals[v], lbglb) )
         {
            assert(boundtypes == NULL || boundtypes[v] == SCIP_BOUNDTYPE_UPPER);
            reopt->glbconss[nglbconss]->bounds[nvarsadded] = lbglb + 1.0;
            reopt->glbconss[nglbconss]->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_LOWER;
            ++nvarsadded;
         }
         else if( boundtypes != NULL )
         {
            /* we round the solution value to get a 'clean' bound */
            assert(SCIPsetIsIntegral(set, vals[v]));
            roundedval = SCIPsetRound(set, vals[v]);

            /* case 3.1 */
            if( boundtypes[v] == SCIP_BOUNDTYPE_UPPER )
            {
               reopt->glbconss[nglbconss]->bounds[nvarsadded] = roundedval + 1.0;
               reopt->glbconss[nglbconss]->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_LOWER;
               ++nvarsadded;
            }
            /* case 3.2 */
            else
            {
               assert(boundtypes[v] == SCIP_BOUNDTYPE_LOWER);
               reopt->glbconss[nglbconss]->bounds[nvarsadded] = roundedval - 1.0;
               reopt->glbconss[nglbconss]->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_UPPER;
               ++nvarsadded;
            }
         }
         /* case 4: in this case we have to add two clauses: (x <= val-1) and (x >= val+1) */
         else
         {
            /* we round the solution value to get a 'clean' bound */
            assert(SCIPsetIsIntegral(set, vals[v]));
            roundedval = SCIPsetRound(set, vals[v]);

            /* first clause: x <= val-1 */
            reopt->glbconss[nglbconss]->bounds[nvarsadded] = roundedval - 1.0;
            reopt->glbconss[nglbconss]->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_UPPER;
            ++nvarsadded;

            /* second clause:  x >= val+1 */
            reopt->glbconss[nglbconss]->vars[nvarsadded] = vars[v];
            reopt->glbconss[nglbconss]->bounds[nvarsadded] = roundedval + 1.0;
            reopt->glbconss[nglbconss]->boundtypes[nvarsadded] = SCIP_BOUNDTYPE_LOWER;
            ++nvarsadded;
         }
      }
   }
   assert(nvars <= nvarsadded);
   assert(nvarsadded == nbinvars + 2*nintvars);

   reopt->glbconss[nglbconss]->nvars = nvarsadded;
   ++reopt->nglbconss;

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
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_BOUNDTYPE* boundtypes;
      int allocmem;
      int nbranchvars;
      int nbinvars;
      int nintvars;
      int v;

      /* allocate memory to store the infeasible path */
      allocmem = SCIPnodeGetDepth(node);
      SCIP_CALL( SCIPsetAllocBufferArray(set, &vars, allocmem) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, allocmem) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &boundtypes, allocmem) );

      /* get the branching path */
      SCIPnodeGetAncestorBranchings(node, vars, vals, boundtypes, &nbranchvars, allocmem);

      if( allocmem < nbranchvars )
      {
         SCIP_CALL( SCIPsetReallocBufferArray(set, &vars, nbranchvars) );
         SCIP_CALL( SCIPsetReallocBufferArray(set, &vals, nbranchvars) );
         SCIP_CALL( SCIPsetReallocBufferArray(set, &boundtypes, nbranchvars) );
         allocmem = nbranchvars;

         SCIPnodeGetAncestorBranchings(node, vars, vals, boundtypes, &nbranchvars, allocmem);
      }

      /* we count the number of binary and (impl) integer variables */
      nbinvars = 0;
      nintvars = 0;
      for( v = 0; v < nbranchvars; v++ )
      {
         if( SCIPvarIsBinary(vars[v]) == SCIP_VARTYPE_BINARY )
            ++nbinvars;
         if( SCIPvarIsBinary(vars[v]) == SCIP_VARTYPE_INTEGER || SCIPvarIsBinary(vars[v]) == SCIP_VARTYPE_IMPLINT )
            ++nintvars;
      }
      assert(nbinvars + nintvars == nbranchvars);

      SCIP_CALL( addGlobalCut(reopt, blkmem, set, vars, vals, boundtypes, nbranchvars, nbinvars, nintvars) );

      /* free buffer */
      SCIPsetFreeBufferArray(set, &boundtypes);
      SCIPsetFreeBufferArray(set, &vals);
      SCIPsetFreeBufferArray(set, &vars);
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

   assert(reopt->reopttree->reoptnodes[id]->dualconscur->constype == REOPT_CONSTYPE_DUALREDS
         || reopt->reopttree->reoptnodes[id]->dualconscur->constype == REOPT_CONSTYPE_INFSUBTREE);

   if( reopt->reopttree->reoptnodes[id]->dualconscur->constype == REOPT_CONSTYPE_DUALREDS )
   {
      SCIPdebugMessage(" create a split-node #%lld\n", SCIPnodeGetNumber(node));
   }
   else
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
      newbound = reopt->reopttree->reoptnodes[id]->dualconscur->bounds[0];
      boundtype = reopt->reopttree->reoptnodes[id]->dualconscur->boundtypes[0];

      assert(SCIPvarIsOriginal(var));
      SCIP_CALL( SCIPvarGetProbvarBound(&var, &newbound, &boundtype) );
      assert(SCIPvarIsTransformed(var));

      oldlb = SCIPvarGetLbLocal(var);
      oldub = SCIPvarGetUbLocal(var);

      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         newbound = reopt->reopttree->reoptnodes[id]->dualconscur->bounds[0] - 1.0;
         assert(SCIPisLE(scip, newbound, oldub));
      }
      else
      {
         newbound = reopt->reopttree->reoptnodes[id]->dualconscur->bounds[0] + 1.0;
         assert(SCIPisGE(scip, newbound, oldlb));
      }
      boundtype = (SCIP_BOUNDTYPE) (1 - (int)boundtype);
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

      SCIPdebugMessage("  -> constraint consists of only one variable: <%s> %s %g\n", SCIPvarGetName(var),
            boundtype == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);
   }
   else
   {
      SCIP_REOPTCONSDATA* consdata;
      SCIP_VAR** consvars;
      SCIP_Real consval;
      SCIP_BOUNDTYPE consboundtype;
      int nbinvars = 0;
      int nintvars = 0;
      int ncontvars = 0;

      consdata = reopt->reopttree->reoptnodes[id]->dualconscur;

      /* allocate buffer */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, consdata->nvars) );

      /* count number of binary, integer, and continuous variables */
      for( v = 0; v < consdata->nvars; v++ )
      {
         switch ( SCIPvarGetType(consdata->vars[v]) ) {
         case SCIP_VARTYPE_BINARY:
            ++nbinvars;
            break;
         case SCIP_VARTYPE_IMPLINT:
         case SCIP_VARTYPE_INTEGER:
            if( SCIPisEQ(scip, SCIPvarGetLbLazy(consdata->vars[v]), 0.0)
             && SCIPisEQ(scip, SCIPvarGetUbLazy(consdata->vars[v]), 1.0) )
               ++nbinvars;
            else
               ++nintvars;
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            ++ncontvars;
            break;
         default:
            SCIPerrorMessage("Variable <%s> has to be either binary, (implied) integer, or continuous.\n",
               SCIPvarGetName(consdata->vars[v]));
            return SCIP_INVALIDDATA;
         }
      }

      if( consdata->constype == REOPT_CONSTYPE_INFSUBTREE )
         name = "infsubtree";
      else
      {
         assert(consdata->constype == REOPT_CONSTYPE_DUALREDS);
         name = "splitcons";
      }

      /* case 1: all variables are binary. we use a logic-or constraint. */
      if( consdata->nvars == nbinvars )
      {
         for(v = 0; v < consdata->nvars; v++)
         {
            consvars[v] = consdata->vars[v];
            consval = consdata->bounds[v];
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

         SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, consdata->nvars, consvars,
               FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
      }
      /* case 2: at least one variable is integer or continuous. we use a bounddisjunction constraint. */
      else
      {
         SCIP_Real* consvals;
         SCIP_BOUNDTYPE* consboundtypes;

         /* alloc buffer memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consboundtypes, consdata->nvars) );

         /* iterate over all variable and transform them */
         for( v = 0; v < consdata->nvars; v++ )
         {
            consvars[v] = consdata->vars[v];
            consvals[v] = consdata->bounds[v];
            consboundtypes[v] = consdata->boundtypes[v];

            /* we have to switch the bounds.
             * case 1: integer variable with bound x <= u is transformed to u+1 <= x
             *                                 and l <= x is transformed to   x <= l-1
             * case 2: continuous variable with bound x <= u is transformed to u <= x
             *                                    and l <= x is transformed to x <= l
             */
            if( SCIPvarGetType(consvars[v]) == SCIP_VARTYPE_BINARY
             || SCIPvarGetType(consvars[v]) == SCIP_VARTYPE_INTEGER
             || SCIPvarGetType(consvars[v]) == SCIP_VARTYPE_IMPLINT )
            {
               if( consboundtypes[v] == SCIP_BOUNDTYPE_UPPER )
               {
                  consvals[v] += 1.0;
                  assert(SCIPsetIsLE(set, consvals[v], SCIPvarGetUbLocal(consvars[v])));
               }
               else
               {
                  consvals[v] -= 1.0;
                  assert(SCIPsetIsGE(set, consvals[v], SCIPvarGetLbLocal(consvars[v])));
               }
            }

            consboundtypes[v] = (SCIP_BOUNDTYPE)(1 - consboundtypes[v]);

            assert(SCIPvarIsOriginal(consvars[v]));
            SCIP_CALL( SCIPvarGetProbvarBound(&consvars[v], &consvals[v], &consboundtypes[v]) );
            assert(SCIPvarIsTransformed(consvars[v]));
            assert(SCIPvarGetStatus(consvars[v]) != SCIP_VARSTATUS_MULTAGGR);
         }

         /* create the constraints and add them to the corresponding nodes */
         SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, consdata->nvars, consvars, consboundtypes,
               consvals, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

         /* free buffer memory */
         SCIPfreeBufferArray(scip, &consboundtypes);
         SCIPfreeBufferArray(scip, &consvals);
      }

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
      val = reoptnode->dualconscur->bounds[v];
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

         reoptnode->varbounds[pos] = reoptnode->dualconscur->bounds[v];
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
      SCIP_REOPTCONSDATA* consdata;
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      SCIP_BOUNDTYPE* consboundtypes;
      int v;

      consdata = reopt->reopttree->reoptnodes[id]->conss[c];
      assert(consdata != NULL);
      assert(consdata->nvars > 0);
      assert(consdata->varssize >= consdata->nvars);

      if( consdata->constype == REOPT_CONSTYPE_INFSUBTREE || consdata->constype == REOPT_CONSTYPE_DUALREDS )
      {
         int nbinvars = 0;
         int nintvars = 0;
         int ncontvars = 0;

         /* count number of binary, integer, and continuous variables */
         for( v = 0; v < consdata->nvars; v++ )
         {
            switch ( SCIPvarGetType(consdata->vars[v]) ) {
            case SCIP_VARTYPE_BINARY:
               ++nbinvars;
               break;
            case SCIP_VARTYPE_IMPLINT:
            case SCIP_VARTYPE_INTEGER:
               if( SCIPisEQ(scip, SCIPvarGetLbLazy(consdata->vars[v]), 0.0)
                && SCIPisEQ(scip, SCIPvarGetUbLazy(consdata->vars[v]), 1.0) )
                  ++nbinvars;
               else
                  ++nintvars;
               break;
            case SCIP_VARTYPE_CONTINUOUS:
               ++ncontvars;
               break;
            default:
               SCIPerrorMessage("Variable <%s> has to be either binary, (implied) integer, or continuous.\n",
                  SCIPvarGetName(consdata->vars[v]));
               return SCIP_INVALIDDATA;
            }
         }

         if( consdata->constype == REOPT_CONSTYPE_INFSUBTREE )
            name = "infsubtree";
         else
            name = "splitcons";


         /* allocate buffer */
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consboundtypes, consdata->nvars) );

         /* case 1: all variables are binary. we use a logic-or constraint. */
         if( consdata->nvars == nbinvars )
         {
            /* iterate over all variable and transform them */
            for( v = 0; v < consdata->nvars; v++ )
            {
               consvars[v] = consdata->vars[v];
               consvals[v] = consdata->bounds[v];
               consboundtypes[v] = consdata->boundtypes[v];
               assert((SCIPsetIsFeasEQ(set, consvals[v], 0.0) && consboundtypes[v] == SCIP_BOUNDTYPE_UPPER)
                  || (SCIPsetIsFeasEQ(set, consvals[v], 1.0) && consboundtypes[v] == SCIP_BOUNDTYPE_LOWER));

               assert(SCIPvarIsOriginal(consvars[v]));
               SCIP_CALL( SCIPvarGetProbvarBound(&consvars[v], &consvals[v], &consboundtypes[v]) );
               assert(SCIPvarIsTransformed(consvars[v]));
               assert(SCIPvarGetStatus(consvars[v]) != SCIP_VARSTATUS_MULTAGGR);

               if( SCIPsetIsFeasEQ(set, consvals[v], 1.0) )
               {
                  SCIP_CALL( SCIPvarNegate(consvars[v], blkmem, set, stat, &consvars[v]) );
                  assert(SCIPvarIsNegated(consvars[v]));
               }
            }

            /* create the constraints and add them to the corresponding nodes */
            SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, consdata->nvars, consvars,
                  FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
         }
         /* case 2: at least one variable is integer or continuous. we use a bounddisjunction constraint. */
         else
         {
            /* iterate over all variable and transform them */
            for( v = 0; v < consdata->nvars; v++ )
            {
               consvars[v] = consdata->vars[v];
               consvals[v] = consdata->bounds[v];
               consboundtypes[v] = consdata->boundtypes[v];

               /* we have to switch the bounds.
                * case 1: integer variable with bound x <= u is transformed to u+1 <= x
                *                                 and l <= x is transformed to   x <= l-1
                * case 2: continuous variable with bound x <= u is transformed to u <= x
                *                                    and l <= x is transformed to x <= l
                */
               if( SCIPvarGetType(consvars[v]) == SCIP_VARTYPE_BINARY
                || SCIPvarGetType(consvars[v]) == SCIP_VARTYPE_INTEGER
                || SCIPvarGetType(consvars[v]) == SCIP_VARTYPE_IMPLINT )
               {
                  if( consboundtypes[v] == SCIP_BOUNDTYPE_UPPER )
                  {
                     consvals[v] += 1.0;
                     assert(SCIPsetIsLE(set, consvals[v], SCIPvarGetUbLocal(consvars[v])));
                  }
                  else
                  {
                     consvals[v] -= 1.0;
                     assert(SCIPsetIsGE(set, consvals[v], SCIPvarGetLbLocal(consvars[v])));
                  }
               }

               consboundtypes[v] = (SCIP_BOUNDTYPE)(1 - consboundtypes[v]);

               assert(SCIPvarIsOriginal(consvars[v]));
               SCIP_CALL( SCIPvarGetProbvarBound(&consvars[v], &consvals[v], &consboundtypes[v]) );
               assert(SCIPvarIsTransformed(consvars[v]));
               assert(SCIPvarGetStatus(consvars[v]) != SCIP_VARSTATUS_MULTAGGR);
            }

            /* create the constraints and add them to the corresponding nodes */
            SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, consdata->nvars, consvars, consboundtypes,
                  consvals, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
         }

         SCIPdebugPrintCons(scip, cons, NULL);

         SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }

      /* free buffer */
      SCIPfreeBufferArray(scip, &consboundtypes);
      SCIPfreeBufferArray(scip, &consvals);
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
         SCIP_CALL( SCIPqueueInsert(reopt->reopttree->openids, (void*) (size_t) redchilds[nredchilds-1]) );

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
   SCIP_VAR**            origvars,           /**< original problem variables */
   int                   norigvars           /**< number of original problem variables */
   )
{
   int probidx;
   int v;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(origvars != NULL);
   assert(norigvars >= 0);

   /* check memory */
   SCIP_CALL( ensureRunSize(reopt, reopt->run, blkmem) );

   /* get memory */
   SCIP_ALLOC( BMSallocClearMemoryArray(&reopt->objs[reopt->run-1], norigvars) ); /*lint !e866*/

   /* save coefficients */
   for( v = 0; v < norigvars; v++ )
   {
      assert(SCIPvarIsOriginal(origvars[v]));

      probidx = SCIPvarGetProbindex(origvars[v]);
      assert(0 <= probidx && probidx <= norigvars);

      reopt->objs[reopt->run-1][probidx] = SCIPvarGetObj(origvars[v]);

      /* update flag to remember if the objective function has changed */
      if( !reopt->objhaschanged && reopt->run >= 2
          && SCIPsetIsEQ(set, reopt->objs[reopt->run-2][probidx], reopt->objs[reopt->run-1][probidx]) )
         reopt->objhaschanged = TRUE;

      /* mark this objective as the first non empty */
      if( reopt->firstobj == -1 && reopt->objs[reopt->run-1][probidx] != 0 )
         reopt->firstobj = reopt->run-1;
   }

   /* calculate similarity to last objective */
   if( reopt->run-1 >= 1 )
   {
      /* calculate similarity to last objective */
      reopt->simtolastobj = reoptSimilarity(reopt, set, reopt->run-1, reopt->run-2, origvars, norigvars);

      SCIPdebugMessage("new objective has similarity of %g.\n", reopt->simtolastobj);
      printf("new objective has similarity of %g.\n", reopt->simtolastobj);
   }

   SCIPdebugMessage("saved obj for run %d.\n", reopt->run);

   return SCIP_OKAY;
}

/** orders the variable by infernce score */
static
SCIP_RETCODE getInferenceOrder(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR**            vars,               /**< variable array to permute */
   SCIP_Real*            bounds,             /**< bound array to permute in the same order */
   SCIP_BOUNDTYPE*       boundtypes,         /**< boundtype array to permute in the same order */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             isroot              /**< is the current node the root node */
   )
{
   SCIP_Real* infscore;
   int v;

   assert(set != NULL);
   assert(vars != NULL);
   assert(bounds != NULL);
   assert(boundtypes != NULL);
   assert(nvars >= 0);

   /* allocate buffer for the scores */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &infscore, nvars) );

   for( v = 0; v < nvars; v++ )
   {
      if( SCIPsetIsEQ(set, bounds[v], 1.0) )
      {
         infscore[v] = 0.75 * SCIPvarGetAvgInferences(vars[v], stat, SCIP_BRANCHDIR_UPWARDS)
            + 0.25 * SCIPvarGetAvgInferences(vars[v], stat, SCIP_BRANCHDIR_DOWNWARDS);
      }
      else
      {
         infscore[v] = 0.25 * SCIPvarGetAvgInferences(vars[v], stat, SCIP_BRANCHDIR_UPWARDS)
               + 0.75 * SCIPvarGetAvgInferences(vars[v], stat, SCIP_BRANCHDIR_DOWNWARDS);
      }
   }

   /* sort vars and vals by score */
   SCIPsortDownRealRealPtrPtr(infscore, bounds, (void*) vars, (void*) boundtypes, nvars);

   /* free buffer */
   SCIPsetFreeBufferArray(set, &infscore);

   return SCIP_OKAY;
}

/** permute the variable and bound array randomly */
static
void permuteRandom(
   SCIP_VAR**            vars,               /**< variable array to permute */
   SCIP_Real*            bounds,             /**< bound array to permute in the same order */
   SCIP_BOUNDTYPE*       boundtypes,         /**< boundtype array to permute in the same order */
   int                   nvars,              /**< number of variables */
   unsigned int*         randseed            /**< seed value for the random generator */
   )
{
   SCIP_VAR* tmpvar;
   SCIP_Real tmpbound;
   SCIP_BOUNDTYPE tmpboundtype;
   int end;
   int i;

   assert(vars != NULL);
   assert(bounds != NULL);
   assert(boundtypes != NULL);

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

      /* swap the last bound and the random bound */
      tmpbound = bounds[i];
      bounds[i] = bounds[end];
      bounds[end] = tmpbound;

      /* swap the last boundtype and the random boundtype */
      tmpboundtype = boundtypes[i];
      boundtypes[i] = boundtypes[end];
      boundtypes[end] = tmpboundtype;
   }
}

static
SCIP_RETCODE separateSolution(
   SCIP_REOPT*           reopt,
   BMS_BLKMEM*           blkmem,
   SCIP_SET*             set,
   SCIP_STAT*            stat,
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_VAR**            vars,               /**< array of original problem variables */
   int                   nvars               /**< number of original problem variables */
   )
{
   SCIP_VAR** origvars;
   SCIP_Real* vals;
   int nintvars;
   int nbinvars;
   int v;
   int w;

   assert(reopt != NULL);
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(vars != NULL);
   assert(nvars != 0);
   assert(SCIPsolIsOriginal(sol));

   /* allocate buffer memory */
   SCIP_ALLOC( BMSallocMemoryArray(&origvars, nvars) );
   SCIP_ALLOC( BMSallocMemoryArray(&vals, nvars) );

   nbinvars = 0;
   nintvars = 0;

   /* get the solution values of the variables */
   for( v = 0, w = 0; v < nvars; v++ )
   {
      assert(SCIPvarIsOriginal(vars[v]));
      assert(nbinvars + nintvars == w);

      /* we do not want to create cuts for continous variables */
      if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY )
         ++nbinvars;
      if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_IMPLINT )
         ++nintvars;

      origvars[v] = vars[v];
      assert(origvars[v] != NULL);
      assert(SCIPvarIsOriginal(origvars[v]));

      vals[w] = SCIPsolGetVal(sol, set, stat, origvars[v]);
      ++w;
   }

   SCIP_CALL( addGlobalCut(reopt, blkmem, set, origvars, vals, NULL, w, nbinvars, nintvars) );

   /* free buffer memory */
   BMSfreeMemoryArray(&vals);
   BMSfreeMemoryArray(&origvars);

   return SCIP_OKAY;
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
   (*reopt)->varhistory = NULL;

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

         /* we have to free all optimal solution separatly, because those solutions are not stored in the
          * solution reopt_sepabestsol = TRUE
          */
         if( set->reopt_sepabestsol && (*reopt)->prevbestsols[p] != NULL )
         {
            SCIP_CALL( SCIPsolFree(&(*reopt)->prevbestsols[p], blkmem, origprimal) );
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
         BMSfreeBlockMemoryArray(blkmem, &(*reopt)->dualcons->boundtypes, (*reopt)->dualcons->varssize);
         BMSfreeBlockMemoryArray(blkmem, &(*reopt)->dualcons->bounds, (*reopt)->dualcons->varssize);
         BMSfreeBlockMemoryArray(blkmem, &(*reopt)->dualcons->vars, (*reopt)->dualcons->varssize);
         BMSfreeBlockMemory(blkmem, &(*reopt)->dualcons);
         (*reopt)->dualcons = NULL;
      }
   }

   if( (*reopt)->glbconss != NULL && (*reopt)->allocmemglbconss > 0 )
   {
      /* free all constraint */
      while( (*reopt)->nglbconss > 0 )
      {
         int c;
         c = (*reopt)->nglbconss;

         if( (*reopt)->glbconss[c] != NULL )
         {
            if( (*reopt)->glbconss[c]->varssize > 0 )
            {
               BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss[c]->boundtypes, (*reopt)->glbconss[c]->varssize);
               BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss[c]->bounds, (*reopt)->glbconss[c]->varssize);
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
   SCIP_PRIMAL*          origprimal,         /**< original primal */
   SCIP_VAR**            vars,               /**< original problem variables */
   int                   nvars               /**< number of original problem variables */
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

   /* store a global constraint that cutsoff the solution */
   if( set->reopt_sepabestsol )
   {
      SCIP_CALL( separateSolution(reopt, blkmem, set, stat, sol, vars, nvars) );
   }

   return SCIP_OKAY;
}

/** add a new iteration after changing the objective function */
SCIP_RETCODE SCIPreoptAddRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data sturcture */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            origvars,           /**< original problem variables */
   int                   norigvars,          /**< number of original variables */
   int                   size                /**< number of expected solutions */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem !=  NULL);
   assert(origvars != NULL);

   /* increase number of runs */
   ++reopt->run;

   /* check memory */
   SCIP_CALL( ensureRunSize(reopt, reopt->run, blkmem) );

   /* allocate memory */
   reopt->soltree->solssize[reopt->run-1] = size;
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->soltree->sols[reopt->run-1], size) ); /*lint !e866*/

   /* reset flag */
   reopt->objhaschanged = FALSE;

   /* save the objective function */
   SCIP_CALL( reoptSaveNewObj(reopt, set, blkmem, origvars, norigvars) );

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

   if( reopt->soltree->root->child != NULL )
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
   SCIP_VAR**            origvars,           /**< original problem variables */
   int                   norigvars           /**< number of original problem variables */
   )
{
   assert(reopt != NULL);
   assert(run1 > 0 && run1 <= reopt->run);
   assert(run2 > 0 && run2 <= reopt->run);
   assert(origvars != NULL);
   assert(norigvars >= 0);

   return reoptSimilarity(reopt, set, run1-1, run2-1, origvars, norigvars);
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
      return reopt->prevbestsols[reopt->run-2];
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

/** reset solving specific paramters */
SCIP_RETCODE SCIPreoptReset(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int c;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* clean addedconss array */
   for( c = 0; c < reopt->naddedconss; c++)
   {
      SCIP_CONS* cons;

      cons = reopt->addedconss[c];
      assert(cons != NULL);

      SCIP_CALL( SCIPconsRelease(&cons, blkmem, set) );
      reopt->addedconss[c] = 0;
   }

   reopt->naddedconss = 0;
   reopt->consadded = FALSE;
   reopt->objhaschanged = FALSE;

   return SCIP_OKAY;
}

/** reset marks of stored solutions to not updated */
void SCIPreoptResetSolMarks(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   )
{
   SCIP_SOLNODE* child;

   assert(reopt != NULL);
   assert(reopt->soltree != NULL);
   assert(reopt->soltree->root != NULL);

   child = reopt->soltree->root->child;

   /* travers through the list */
   while( child != NULL )
   {
      soltreeResetMarks(child);
      child = child->sibling;
   }
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

   assert(reoptnode->dualconscur->constype == REOPT_CONSTYPE_DUALREDS);

   /* copy the variable information */
   for(v = 0; v < *nvars; v++)
   {
      vars[v] = reoptnode->dualconscur->vars[v];
      vals[v] = reoptnode->dualconscur->bounds[v];
   }

   *constype = reoptnode->dualconscur->constype;

   return;
}

/** returns all added constraints at ID id */
void SCIPreoptnodeGetConss(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR***           vars,               /**< 2-dim array of variables */
   SCIP_Real**           bounds,             /**< 2-dim array of bounds */
   SCIP_BOUNDTYPE**      boundtypes,         /**< 2-dim array of boundtypes */
   int                   mem,                /**< allocated memory for constraints */
   int*                  nconss,             /**< pointer to store the number of constraints */
   int*                  nvars               /**< pointer to store the number of variables */
   )
{
   int c;

   assert(reoptnode != NULL);
   assert(vars != NULL);
   assert(bounds != NULL);
   assert(boundtypes != NULL);
   assert(nvars != NULL);

   (*nconss) = reoptnode->nconss;

   if( mem < *nconss )
      return;

   for(c = 0; c < *nconss; c++)
   {
      assert(vars[c] != NULL);
      assert(bounds[c] != NULL);

      vars[c] = reoptnode->conss[c]->vars;
      bounds[c] = reoptnode->conss[c]->bounds;
      boundtypes[c] = reoptnode->conss[c]->boundtypes;
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
      SCIP_CALL( saveGlobalCons(reopt, set, blkmem, node, REOPT_CONSTYPE_CUT) );
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
   SCIP_LP*              lp,
   SCIP_LPSOLSTAT        lpsolstat,          /**< solution status of the LP */
   SCIP_Bool             isrootnode,         /**< the node is the root */
   SCIP_Bool             isfocusnode,        /**< the node is the current focus node */
   SCIP_Real             lowerbound,         /**< lower bound of the node */
   int                   effectiverootdepth  /**< effective root depth */
   )
{
   SCIP_Bool strongbranched;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(lp != NULL);
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
      /* store active cuts */
      if( set->reopt_usecuts )
      {
         SCIP_ROW** lprows;
         int nlprows;
         int r;

         lprows = SCIPlpGetRows(lp);
         nlprows = SCIPlpGetNRows(lp);

         for( r = 0; r < nlprows; r++ )
         {
            /* we can break if we reach teh first row that is not part of the current LP */
            if( SCIProwGetLPPos(lprows[r]) == -1 )
               break;

            /* currently we only want to store cuts generated by a seperator */
            if( SCIProwGetOrigintype(lprows[r]) == SCIP_ROWORIGINTYPE_SEPA && SCIProwGetAge(lprows[r]) <= set->reopt_maxcutage )
            {
               SCIP_VAR** cutvars;
               SCIP_COL** cols;
               SCIP_Real* cutvals;
               SCIP_Real lhs;
               SCIP_Real rhs;
               int ncutvars;
               int c;

               ncutvars = SCIProwGetNLPNonz(lprows[r]);
               lhs = SCIProwGetLhs(lprows[r]);
               rhs = SCIProwGetRhs(lprows[r]);

               /* subtract row constant */
               if( !SCIPsetIsInfinity(set, -lhs) )
                  lhs -= SCIProwGetConstant(lprows[r]);
               if( !SCIPsetIsInfinity(set, rhs) )
                  rhs -= SCIProwGetConstant(lprows[r]);

               cutvals = SCIProwGetVals(lprows[r]);
               cols = SCIProwGetCols(lprows[r]);

               SCIP_CALL( SCIPsetAllocBufferArray(set, &cutvars, ncutvars) );

               for( c = 0; c < ncutvars; c++ )
               {
                  SCIP_Real constant;
                  SCIP_Real scalar;

                  cutvars[c] = SCIPcolGetVar(cols[c]);
                  assert(cutvars[c] != NULL);

                  constant = 0.0;
                  scalar = 1.0;

                  SCIP_CALL( SCIPvarGetOrigvarSum(&cutvars[c], &scalar, &constant) );
                  assert(cutvars[c] != NULL);
                  cutvals[c] = (cutvals[c] - constant)/scalar;
               }

               SCIP_CALL( SCIPreoptnodeAddCons(reopt->reopttree->reoptnodes[0], set, blkmem, cutvars, cutvals, NULL,
                     lhs, rhs, ncutvars, REOPT_CONSTYPE_CUT) );

               SCIPsetFreeBufferArray(set, &cutvars);
            }
         }

      }

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
               SCIPdebugMessage(" -> new constype    : %d\n", REOPT_CONSTYPE_DUALREDS);
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
               SCIPdebugMessage(" -> new constype    : %d\n", REOPT_CONSTYPE_DUALREDS);
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
      SCIP_BOUNDTYPE boundtype;
      int allocmem;

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

      if( SCIPsetIsEQ(set, oldval, newval) )
      {
         SCIPerrorMessage("cannot store equal bounds: old = %g, new = %g\n", oldval, newval);
         return SCIP_INVALIDDATA;
      }

      if( SCIPsetIsLT(set, newval, oldval) )
         boundtype = SCIP_BOUNDTYPE_UPPER;
      else
         boundtype = SCIP_BOUNDTYPE_LOWER;

      reopt->dualcons->vars[reopt->dualcons->nvars] = var;
      reopt->dualcons->bounds[reopt->dualcons->nvars] = newval;
      reopt->dualcons->boundtypes[reopt->dualcons->nvars] = boundtype;
      ++reopt->dualcons->nvars;

      SCIPdebugMessage(">> store %s bound change of <%s>: %g -> %g\n", boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper",
            SCIPvarGetName(var), oldval, newval);
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

/** merges the variable history of the current run with the stored history */
SCIP_RETCODE SCIPreoptMergeVarHistory(
   SCIP_REOPT*           reopt,
   SCIP_STAT*            stat,
   SCIP_VAR**            vars,
   int                   nvars
   )
{
   int idx;
   int v;

   assert(reopt != NULL);
   assert(stat != NULL);
   assert(nvars >= 0);

   if( reopt->simtolastobj < 0.975 || reopt->varhistory == NULL )
   {
      if( reopt->varhistory != NULL )
      {
         for( v = 0; v < nvars; v++ )
         {
            idx = SCIPvarGetProbindex(vars[v]);
            assert(idx >= 0 && idx < nvars);

            SCIPhistoryReset(reopt->varhistory[idx]);
         }
      }

      return SCIP_OKAY;
   }

   for( v = 0; v < nvars; v++ )
   {
      idx = SCIPvarGetProbindex(vars[v]);
      assert(0 <= idx && idx <= nvars);

      SCIPvarSetHistory(vars[v], reopt->varhistory[idx], stat);
   }

   return SCIP_OKAY;
}

/** updates the variable history */
SCIP_RETCODE SCIPreoptUpdateVarHistory(
   SCIP_REOPT*           reopt,
   SCIP_SET*             set,
   SCIP_STAT*            stat,
   SCIP_PROB*            origprob,
   SCIP_PROB*            transprob,
   BMS_BLKMEM*           blkmem,
   SCIP_VAR**            vars,
   int                   nvars
   )
{
   SCIP_SOL* lastoptsol;
   SCIP_Real objdelta;
   SCIP_Real oldobj;
   int v;

   assert(reopt != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(nvars >= 0);

   SCIPdebugMessage("updating variable history\n");
   printf("updating variable history\n");

   if( reopt->varhistory == NULL )
   {
      /* allocate memory */
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->varhistory, nvars) );

      for( v = 0; v < nvars; v++ )
      {
         SCIP_CALL( SCIPhistoryCreate(&(reopt->varhistory[v]), blkmem) );
      }
   }

   /* optimal solution of the current run */
   lastoptsol = reopt->prevbestsols[reopt->run-1];

   /* best objective value of the previous run */
   if( reopt->run == 1 )
      oldobj = 0;
   else
   {
      assert(reopt->run >= 2);

      oldobj = 0;
      for( v = 0; v < nvars; v++ )
      {
         SCIP_VAR* origvar;
         SCIP_Real solval;
         int probidx;

         origvar = vars[v];
         if( !SCIPvarIsOriginal(origvar) )
         {
            SCIP_Real constant;
            SCIP_Real scalar;

            constant = 0.0;
            scalar = 1.0;

            SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );
         }
         assert(origvar != NULL && SCIPvarIsOriginal(origvar));

         probidx = SCIPvarGetProbindex(origvar);
         solval = SCIPsolGetVal(lastoptsol, set, stat, origvar);
         oldobj += reopt->objs[reopt->run-2][probidx] * solval;
      }
   }

   /* the objective delta is defined as: new LP obj - old LP obj and has to be non-negative. for the reoptimization we
    * are interested in the obj value of the best solution of the current iteration wrt the objective function of the
    * previous iteration instead of the old LP obj.
    * hence, c_{i}(x_{i}) <= c_{i-1}(x_{i}) <==> c_{i}(x_{i}) - c_{i-1}(x_{i}) <= 0 we multiply with -1.
    */
   objdelta = SCIPsolGetObj(lastoptsol, set, transprob, origprob) - oldobj;
   objdelta *= -1.0;
   assert(!SCIPsetIsNegative(set, objdelta));

   /* update the history and scale them */
   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real solvedelta;
      SCIP_Real avginference[2];
      SCIP_Real avgcutoff[2];
      int idx;
      int d;

      idx = SCIPvarGetProbindex(vars[v]);
      assert(idx >= 0 && idx < nvars);

      solvedelta = SCIPvarGetAvgSol(vars[v]);

      /* get and scale stored data and data from the current iteration for both directions */
      for( d = 0; d <= 1; d++ )
      {
         /* inference score */
         avginference[d] = SCIPvarGetAvgInferences(vars[v], stat, (SCIP_BRANCHDIR)d);
         avginference[d] += SCIPhistoryGetAvgInferences(reopt->varhistory[idx], (SCIP_BRANCHDIR)d);
         avginference[d] /= 2;

         /* cutoff score */
         avgcutoff[d] = SCIPvarGetAvgCutoffs(vars[v], stat, (SCIP_BRANCHDIR)d);
         avgcutoff[d] += SCIPhistoryGetAvgCutoffs(reopt->varhistory[idx], (SCIP_BRANCHDIR)d);
         avgcutoff[d] /= 2;
      }

      /* reset the history */
      SCIPhistoryReset(reopt->varhistory[idx]);

      /* set the updated history for both directions */
      for( d = 0; d <= 1; d++ )
      {
         SCIPhistoryIncNBranchings(reopt->varhistory[idx], (SCIP_BRANCHDIR)d, 1);

         if( d == 1 )
            SCIPhistoryUpdatePseudocost(reopt->varhistory[idx], set, solvedelta, 0, 1);

         SCIPhistoryIncInferenceSum(reopt->varhistory[idx], (SCIP_BRANCHDIR)d, avginference[d]);
         SCIPhistoryIncCutoffSum(reopt->varhistory[idx], (SCIP_BRANCHDIR)d, avgcutoff[d]);
      }
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
         int v;

         assert(representatives[r]->nvars <= representatives[r]->varssize);

         for( v = 0; v < representatives[r]->nvars; v++ )
         {
            SCIP_CALL( SCIPreoptnodeAddBndchg(reopttree->reoptnodes[id], blkmem, representatives[r]->vars[v],
                  representatives[r]->varbounds[v], representatives[r]->varboundtypes[v]) );
         }
      }

      if( representatives[r]->nconss > 0 )
      {
         int c;

         assert(representatives[r]->nconss <= representatives[r]->consssize);

         for( c = 0; c < representatives[r]->nconss; c++ )
         {
            SCIP_CALL( SCIPreoptnodeAddCons(reopttree->reoptnodes[id], set, blkmem, representatives[r]->conss[c]->vars,
                  representatives[r]->conss[c]->bounds, representatives[r]->conss[c]->boundtypes,
                  representatives[r]->conss[c]->lhs, representatives[r]->conss[c]->rhs,
                  representatives[r]->conss[c]->nvars, representatives[r]->conss[c]->constype) );
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
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          randseed,           /**< seed value for random generator */
   int*                  ncreatedchilds,     /**< pointer to store the number of created nodes */
   int*                  naddedconss         /**< pointer to store the number added constraints */
   )
{
   SCIP_REOPTTREE* reopttree;
   SCIP_REOPTCONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* bounds;
   SCIP_BOUNDTYPE* boundtypes;
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
   assert(stat != NULL);
   assert(blkmem != NULL);

   reopttree = reopt->reopttree;

   nchilds = reopttree->reoptnodes[0]->nchilds;

   assert(reopttree->reoptnodes[0]->dualconscur != NULL);
   nbndchgs = reopttree->reoptnodes[0]->dualconscur->nvars;

   vars = NULL;
   bounds = NULL;
   boundtypes = NULL;
   nvars = 0;

   (*ncreatedchilds) = 0;
   (*naddedconss) = 0;

   if( !set->reopt_usesplitcons )
   {
      vars = reopttree->reoptnodes[0]->dualconscur->vars;
      bounds = reopttree->reoptnodes[0]->dualconscur->bounds;
      boundtypes = reopttree->reoptnodes[0]->dualconscur->boundtypes;
      nvars = reopttree->reoptnodes[0]->dualconscur->nvars;

      /* calculate the order of the variables */
      switch (set->reopt_varorderinterdiction) {
         /* default order */
         case 'd':
            break;

         /* inference order */
         case 'i':
            SCIP_CALL( getInferenceOrder(set, stat, vars, bounds, boundtypes, nvars, TRUE) );
            break;

         /* random order */
         case 'r':
            permuteRandom(vars, bounds, boundtypes, nvars, &randseed);
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
      reopttree->reoptnodes[id]->varbounds[v] = reopttree->reoptnodes[0]->dualconscur->bounds[v];
      reopttree->reoptnodes[id]->varboundtypes[v] = SCIPsetIsFeasEQ(set, reopttree->reoptnodes[0]->dualconscur->bounds[v], 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER; // ???????????????????????
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
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &consdata->bounds, reopttree->reoptnodes[0]->dualconscur->bounds, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &consdata->boundtypes, reopttree->reoptnodes[0]->dualconscur->boundtypes, nbndchgs) );

      consdata->varssize = nbndchgs;
      consdata->nvars = nbndchgs;
      consdata->lhs = reopttree->reoptnodes[0]->dualconscur->lhs;
      consdata->rhs = reopttree->reoptnodes[0]->dualconscur->rhs;
      consdata->constype = REOPT_CONSTYPE_DUALREDS;

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
      assert(bounds != NULL);
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
            reopttree->reoptnodes[id]->varbounds[v] = bounds[v];
            reopttree->reoptnodes[id]->varboundtypes[v] = SCIPsetIsFeasEQ(set, bounds[v], 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
            ++reopttree->reoptnodes[id]->nvars;
         }

         /* set bound change v+1 (= c) to 1-vals[c] */
         assert(v == c);
         reopttree->reoptnodes[id]->vars[c] = vars[c];
         reopttree->reoptnodes[id]->varbounds[c] = 1-bounds[c];
         reopttree->reoptnodes[id]->varboundtypes[c] = SCIPsetIsFeasEQ(set, 1-bounds[c], 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
         ++reopttree->reoptnodes[id]->nvars;

         /* add dummy1 as a child of the root node */
         SCIP_CALL( reoptAddChild(reopttree, blkmem, 0, id) );

         ++(*ncreatedchilds);
      }

      assert(*ncreatedchilds == nvars+1);
   }

   /* free the current dualconscur and assign dualconsnex */
   assert(reopttree->reoptnodes[0]->dualconscur->vars != NULL);
   assert(reopttree->reoptnodes[0]->dualconscur->bounds != NULL);
   assert(reopttree->reoptnodes[0]->dualconscur->boundtypes != NULL);

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
   SCIP_CALL( SCIPqueueInsert(reopt->reopttree->openids, (void*) (size_t) id) );

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
            assert(reoptnode->dualconscur->constype == REOPT_CONSTYPE_DUALREDS);
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
               assert(reoptnode->dualconscur->constype == REOPT_CONSTYPE_DUALREDS);
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
         SCIP_Real* bounds;
         SCIP_BOUNDTYPE* boundtypes;
         int nvars;

         vars = reoptnode->dualconscur->vars;
         bounds = reoptnode->dualconscur->bounds;
         boundtypes = reoptnode->dualconscur->boundtypes;
         nvars = reoptnode->dualconscur->nvars;

         /* calculate the order of the variables */
         switch (set->reopt_varorderinterdiction)
         {
         /* default order */
         case 'd':
            break;

         /* inference order */
         case 'i':
            SCIP_CALL( getInferenceOrder(set, stat, vars, bounds, boundtypes, nvars, FALSE) );
            break;

         /* random order */
         case 'r':
            permuteRandom(vars, bounds, boundtypes, nvars, &randseed);
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
                     blkmem, childnodes[c], id, vars, bounds, nvars, c) );
            }

            /* add all local constraints */
            SCIP_CALL( addLocalConss(scip, reopt, set, stat, blkmem, childnodes[c], id) );

            /* set estimates */
            if( !SCIPsetIsInfinity(set, REALABS(reopt->reopttree->reoptnodes[id]->lowerbound)) )
            {
               if( SCIPsetIsRelGE(set, reoptnode->lowerbound, SCIPnodeGetLowerbound(childnodes[c])) )
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

#ifdef SCIP_DISABLED_CODE
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
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->glbconss[pos]->bounds, &vals, nvars) );
      reopt->glbconss[pos]->varssize = nvars;
      reopt->glbconss[pos]->nvars = nvars;

      ++reopt->nglbconss;
   }

   return SCIP_OKAY;
}
#endif

/** add the stored constraints globally to the problem */
SCIP_RETCODE SCIPreoptApplyGlbConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   char name[SCIP_MAXSTRLEN];
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
      int nbinvars;
      int nintvars;
      int v;

      assert(reopt->glbconss[c]->nvars > 0);

      cons = NULL;
      consvars = NULL;
      nbinvars = 0;
      nintvars = 0;

      /* check if we can use a logic-or or if we have to use a bounddisjuction constraint */
      for(v = 0; v < reopt->glbconss[c]->nvars; v++)
      {
         if( SCIPvarGetType(reopt->glbconss[c]->vars[v]) == SCIP_VARTYPE_BINARY )
            ++nbinvars;
         else if( SCIPvarGetType(reopt->glbconss[c]->vars[v]) == SCIP_VARTYPE_INTEGER
               || SCIPvarGetType(reopt->glbconss[c]->vars[v]) == SCIP_VARTYPE_IMPLINT )
            ++nintvars;
         else
         {
            SCIPerrorMessage("Expected variable type binary or (impl.) integer for variable <%s> in global constraint at pos. %d.\n",
                  SCIPvarGetName(reopt->glbconss[c]->vars[v]), c);
            return SCIP_INVALIDDATA;
         }
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "glb_%s_%d", reopt->glbconss[c]->constype == REOPT_CONSTYPE_CUT ? "cut" : "inf", reopt->run);

      /* all variables are binary, we can create a logic-or constraint */
      if( nbinvars == reopt->glbconss[c]->nvars )
      {
         SCIPdebugMessage("-> add logic-or constraints with %d binvars\n", nbinvars);

         /* allocate buffer */
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, reopt->glbconss[c]->nvars) );

         for(v = 0; v < reopt->glbconss[c]->nvars; v++)
         {
            consvars[v] = reopt->glbconss[c]->vars[v];
            assert(SCIPvarIsOriginal(consvars[v]));

            /* negate the variable if it was fixed to 1 */
            if( SCIPsetIsFeasEQ(set, reopt->glbconss[c]->bounds[v], 0.0) )
            {
               assert(reopt->glbconss[c]->boundtypes[v] == SCIP_BOUNDTYPE_UPPER);
               SCIP_CALL( SCIPvarNegate(reopt->glbconss[c]->vars[v], blkmem, set, stat, &consvars[v]) );
            }
         }

         /* create the logic-or constraint */
         SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, reopt->glbconss[c]->nvars,
               consvars, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

         /* free buffer */
         SCIPfreeBufferArray(scip, &consvars);
      }
      /* not all variables are binary, we need a bounddisjunction constraint */
      else
      {
         assert(reopt->glbconss[c]->nvars = nbinvars + 2*nintvars);

         SCIPdebugMessage("-> add bounddisjuction constraints with %d binvars, %d intvars\n", nbinvars, (int) (2*nintvars));

         /* create the bounddisjuction constraint */
         SCIP_CALL( SCIPcreateConsBasicBounddisjunction(scip, &cons, name, reopt->glbconss[c]->nvars, reopt->glbconss[c]->vars,
               reopt->glbconss[c]->boundtypes, reopt->glbconss[c]->bounds) );
      }

#ifdef SCIP_MORE_DEBUG
      SCIPdebugPrintCons(scip, cons, NULL);
#endif

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* delete the global constraints data */
      SCIPfreeBlockMemoryArrayNull(scip, &reopt->glbconss[c]->boundtypes, reopt->glbconss[c]->nvars);
      SCIPfreeBlockMemoryArrayNull(scip, &reopt->glbconss[c]->bounds, reopt->glbconss[c]->nvars);
      SCIPfreeBlockMemoryArrayNull(scip, &reopt->glbconss[c]->vars, reopt->glbconss[c]->nvars);
      SCIPfreeBlockMemoryNull(scip, &reopt->glbconss[c]); /*lint !e866*/
   }

   /* reset the number of global constraints */
#ifdef SCIP_DEBUG
   for(c = 0; c < reopt->nglbconss; c++)
   {
      assert(reopt->glbconss[c]->nvars == 0);
      assert(reopt->glbconss[c]->vars == NULL);
      assert(reopt->glbconss[c]->bounds == NULL);
   }
#endif
   reopt->nglbconss = 0;

   return SCIP_OKAY;
}

/** add the stored cuts to the separation storage */
SCIP_RETCODE SCIPreoptApplyCuts(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node,               /**< current focus node */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_Bool             root                /**< bool whether the current node is the root */
   )
{
   SCIP_REOPTNODE* reoptnode;
   SCIP_Bool infeasible;
   unsigned int id;
   int ncuts;
   int c;

   assert(reopt != NULL);
   assert(node != NULL);
   assert(sepastore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(eventqueue != NULL);
   assert(eventfilter != NULL);
   assert(lp != NULL);

   id = SCIPnodeGetReoptID(node);
   assert(id < reopt->reopttree->reoptnodessize);

   /* skip nodes that are node part of the reoptimization tree */
   if( id == 0 && SCIPnodeGetDepth(node) > 0 )
      return SCIP_OKAY;

   reoptnode = reopt->reopttree->reoptnodes[id];
   assert(reoptnode != NULL);

   ncuts = 0;
   for( c = 0; c < reoptnode->nconss; c++ )
   {
      SCIP_REOPTCONSDATA* cons;

      cons = reoptnode->conss[c];
      assert(cons != NULL);

      if( cons->constype == REOPT_CONSTYPE_CUT )
      {
         SCIP_ROW* cut;
         SCIP_COL** cols;
         SCIP_Real* vals;
         char cutname[SCIP_MAXSTRLEN];
         int ncols;
         int v;

         SCIP_CALL( SCIPsetAllocBufferArray(set, &cols, cons->nvars) );
         SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, cons->nvars) );

         ncols = 0;
         for( v = 0; v < cons->nvars; v++ )
         {
            SCIP_VAR* transvar;

            assert(SCIPvarIsOriginal(cons->vars[v]));

            transvar = SCIPvarGetTransVar(cons->vars[v]);
            assert(transvar != NULL);
            assert(SCIPvarGetStatus(transvar) == SCIP_VARSTATUS_COLUMN);

            vals[ncols] = cons->bounds[v];
            cols[ncols] = SCIPvarGetCol(transvar);
            assert(cols[ncols] != NULL);

            ++ncols;
         }
         assert(ncols == cons->nvars);

         (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "reoptcut_%d_%d", id, ncuts);
         SCIP_CALL( SCIProwCreate(&cut, blkmem, set, stat, lp, cutname, ncols, cols, vals, cons->lhs, cons->rhs,
               SCIP_ROWORIGINTYPE_REOPT, NULL, TRUE, TRUE, TRUE) );

         SCIPdebugMessage("add cut <%s> of size %d, [lhs, rhs] = [%g,%g] to node %lld\n", cutname, ncols, cons->lhs,
               cons->rhs, SCIPnodeGetNumber(node));

         SCIPsetFreeBufferArray(set, &vals);
         SCIPsetFreeBufferArray(set, &cols);

         infeasible = FALSE;

         SCIP_CALL( SCIPsepastoreAddCut(sepastore, blkmem, set, stat, eventqueue, eventfilter, lp, NULL, cut, FALSE, root,
               &infeasible) );

         if( infeasible )
         {
            SCIPdebugMessage("constraint %d stored at node %llu (id: %u) is infeasible.\n", c, SCIPnodeGetNumber(node), id);
         }
         else
            ++ncuts;
      }
   }
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
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            vars,               /**< variables which are part of the constraint */
   SCIP_Real*            bounds,             /**< bounds of the variables */
   SCIP_BOUNDTYPE*       boundtypes,         /**< boundtypes of the variables (or NULL is the constraint is a cut) */
   SCIP_Real             lhs,                /**< lhs of the constraint */
   SCIP_Real             rhs,                /**< rhs of the constraint */
   int                   nvars,              /**< number of variables */
   REOPT_CONSTYPE        constype            /**< type of the constraint */
   )
{
   int nconss;

   assert(reoptnode != NULL);
   assert(set != NULL);
   assert(vars != NULL);
   assert(bounds != NULL);
   assert(REOPT_CONSTYPE_CUT || boundtypes != NULL);
   assert(nvars > 0);
   assert(blkmem != NULL);

   /* the constraint can be interpreted as a normal bound change */
   if( nvars == 1 )
   {
      assert(constype == REOPT_CONSTYPE_DUALREDS || constype == REOPT_CONSTYPE_INFSUBTREE);

      SCIPdebugMessage("-> constraint has size 1 -> save as normal bound change.\n");

      if( SCIPvarGetType(vars[0]) == SCIP_VARTYPE_BINARY )
      {
         SCIP_CALL( SCIPreoptnodeAddBndchg(reoptnode, blkmem, vars[0], 1-bounds[0],
               1-bounds[0] == 1 ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER) );
      }
      else
      {
         SCIP_Real newbound;
         SCIP_BOUNDTYPE newboundtype;

         assert(SCIPvarGetType(vars[0]) == SCIP_VARTYPE_INTEGER);

         if( boundtypes[0] == SCIP_BOUNDTYPE_UPPER )
         {
            newbound = bounds[0] + 1.0;
            assert(SCIPsetIsLE(set, newbound, SCIPvarGetUbLocal(vars[0])));

            newboundtype = SCIP_BOUNDTYPE_LOWER;
         }
         else
         {
            newbound = bounds[0] - 1.0;
            assert(SCIPsetIsGE(set, newbound, SCIPvarGetLbLocal(vars[0])));

            newboundtype = SCIP_BOUNDTYPE_UPPER;
         }

         SCIP_CALL( SCIPreoptnodeAddBndchg(reoptnode, blkmem, vars[0], newbound, newboundtype) );
      }
   }
   else
   {
      nconss = reoptnode->nconss;

      SCIP_CALL( reoptnodeCheckMemory(reoptnode, blkmem, 0, 0, nconss+1) );

      /* create the constraint */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &reoptnode->conss[nconss]) ); /*lint !e866*/
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptnode->conss[nconss]->vars, vars, nvars) );

      if( constype != REOPT_CONSTYPE_CUT )
      {
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptnode->conss[nconss]->bounds, bounds, nvars) );
      }
      else
         reoptnode->conss[nconss]->bounds = NULL;

      reoptnode->conss[nconss]->varssize = nvars;
      reoptnode->conss[nconss]->nvars = nvars;
      reoptnode->conss[nconss]->lhs = lhs;
      reoptnode->conss[nconss]->rhs = rhs;
      reoptnode->conss[nconss]->constype = constype;
      ++reoptnode->nconss;
   }
   return SCIP_OKAY;
}

/** add a consraint to the reoptimization data structure */
SCIP_RETCODE SCIPreoptAddCons(
   SCIP_REOPT*           reopt,
   SCIP_SET*             set,
   BMS_BLKMEM*           blkmem,
   SCIP_CONS*            cons
   )
{
   int i;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(cons != NULL);

   /* check memory */
   if( reopt->addedconsssize == 0 )
   {
      assert(reopt->addedconss == NULL);

      reopt->addedconsssize = 10;
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->addedconss, reopt->addedconsssize) );

      /* clear the array */
      for( i = 0; i < reopt->addedconsssize; i++ )
         reopt->addedconss[i] = NULL;
   }
   else if( reopt->naddedconss == reopt->addedconsssize )
   {
      int oldsize;

      oldsize = reopt->addedconsssize;
      reopt->addedconsssize *= 2;
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->addedconss, oldsize, reopt->addedconsssize) );

      /* clear the array */
      for( i = oldsize; i < reopt->addedconsssize; i++ )
         reopt->addedconss[i] = NULL;
   }
   assert(reopt->naddedconss < reopt->addedconsssize);
   assert(reopt->addedconss[reopt->naddedconss] == NULL);

   reopt->addedconss[reopt->naddedconss] = cons;
   reopt->consadded = TRUE;
   ++reopt->naddedconss;

   /* capture the constraint */
   SCIPconsCapture(cons);

   return SCIP_OKAY;
}
