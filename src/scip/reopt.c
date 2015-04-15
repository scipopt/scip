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

/**@file   reopt.c
 * @brief  data structures and methods for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/mem.h"
#include "scip/scip.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/var.h"
#include "scip/misc.h"
#include "scip/reopt.h"
#include "scip/tree.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/cons_logicor.h"
#include "scip/clock.h"
#include "scip/heur_reoptsols.h"
#include "scip/cons.h"

#define DEFAULT_MEM_VARAFTERDUAL    10
#define DEFAULT_MEM_VAR             10
#define DEFAULT_MEM_NODES         1000
#define DEFAULT_MEM_RUN            200
#define DEFAULT_MEM_DUALCONS        10

/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that sols[pos] array can store at least num entries */
static
SCIP_RETCODE ensureSolsSize(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   num,                     /**< minimum number of entries to store */
   int                   runidx                   /**< run index for which the memory should checked */
   )
{
   assert(runidx >= 0);
   assert(runidx <= reopt->runsize);

   if( num > reopt->soltree->solssize[runidx] )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_CALL( SCIPreallocMemoryArray(set->scip, &reopt->soltree->sols[runidx], newsize) );
      reopt->soltree->solssize[runidx] = newsize;
   }
   assert(num <= reopt->soltree->solssize[runidx]);

   return SCIP_OKAY;
}

/** ensures, that sols array can store at least num entries */
static
SCIP_RETCODE ensureRunSize(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   int                   num,                     /**< minimum number of entries to store */
   BMS_BLKMEM*           blkmem                   /**< block memory */
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
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->objs, reopt->runsize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->lastbestsol, reopt->runsize, newsize) );

      for(s = reopt->runsize; s < newsize; s++)
      {
         reopt->lastbestsol[s] = NULL;
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

/*
 * check the memory of the reopttree and of necessary reallocate
 */
static
SCIP_RETCODE reopttreeCheckMemory(
   SCIP_REOPTTREE*       reopttree,               /**< reoptimization tree */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   assert(reopttree != NULL);
   assert(blkmem != NULL);

   if( SCIPqueueIsEmpty(reopttree->openids) )
   {
      int id;

      assert(reopttree->nsavednodes == reopttree->allocmemnodes-1);

      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopttree->reoptnodes, reopttree->allocmemnodes, 2*reopttree->allocmemnodes) );

      for(id = reopttree->allocmemnodes; id < 2*reopttree->allocmemnodes; id++)
      {
         SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void*) (size_t) id) );
         reopttree->reoptnodes[id] = NULL;
      }

      reopttree->allocmemnodes *= 2;
   }

   return SCIP_OKAY;
}

/*
 * check allocated memory of a node within the reopttree and
 * if necessary reallocate
 */
static
SCIP_RETCODE reoptnodeCheckMemory(
   SCIP_REOPTNODE*       reoptnode,               /**< node of the reopttree */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   var_mem,                 /**< memory for variables */
   int                   child_mem,               /**< memory for child nodes */
   int                   conss_mem                /**< memory for constraints */
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
      if( reoptnode->allocvarmem == 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->vars, var_mem) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->varbounds, var_mem) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->varboundtypes, var_mem) );
         reoptnode->allocvarmem = var_mem;
      }
      else if( reoptnode->allocvarmem < var_mem )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->vars, reoptnode->allocvarmem, var_mem) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->varbounds, reoptnode->allocvarmem, var_mem) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->varboundtypes, reoptnode->allocvarmem, var_mem) );
         reoptnode->allocvarmem = var_mem;
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
      if( reoptnode->allocmemconss == 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reoptnode->conss, conss_mem) );
         reoptnode->nconss = 0;
         reoptnode->allocmemconss = conss_mem;
      }
      else if( reoptnode->allocmemconss < conss_mem )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reoptnode->conss, reoptnode->allocmemconss, conss_mem) );
         reoptnode->allocmemconss = conss_mem;
      }
   }

   return SCIP_OKAY;
}

/*
 * local methods
 */

static
int soltreeNInducedtSols(
   SCIP_SOLNODE*         node                      /**< node within the solution tree */
   )
{
   assert(node != NULL);

   if( node->father == NULL && node->rchild == NULL && node->lchild == NULL )
      return 0;
   else if( node->rchild == NULL && node->lchild == NULL )
      return 1;
   else
   {
      if( node->rchild == NULL )
         return soltreeNInducedtSols(node->lchild);
      else if( node->lchild == NULL )
         return soltreeNInducedtSols(node->rchild);
      else
         return soltreeNInducedtSols(node->rchild) + soltreeNInducedtSols(node->lchild);
   }
}

/* returns similarity of two objective functions */
static
SCIP_Real reoptSimilarity(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   int                   obj1_id,                 /**< id of one objective function */
   int                   obj2_id                  /**< id of the other objective fucntion */
   )
{
   SCIP_Real similarity;
   SCIP_Bool onediffertozero;
   int id;

   onediffertozero = FALSE;

   /* calc similarity */
   similarity = 0.0;
   for(id = 0; id < reopt->nobjvars; id++)
   {
      SCIP_Real c1;
      SCIP_Real c2;

      c1 = reopt->objs[obj1_id][id];
      c2 = reopt->objs[obj2_id][id];

      if( c1 != 0 || c2 != 0 )
         onediffertozero = TRUE;

      /** vector product */
      similarity += c1*c2;
   }

   if( !onediffertozero )
      return -2.0;
   else
      return similarity;
}

/**
 * delete the data for node nodeID in nodedata
 */
static
SCIP_RETCODE reoptnodeDelete(
   SCIP_REOPTNODE**      reoptnode,               /**< node of the reoptimization tree */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   assert((*reoptnode) != NULL );
   assert(blkmem != NULL );

   /** delete data for constraints */
   if((*reoptnode)->allocmemconss > 0 )
   {
      int c;

      assert((*reoptnode)->conss != NULL);

      for(c = 0; c < (*reoptnode)->nconss; c++)
      {
         BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->conss[c]->vals, (*reoptnode)->conss[c]->allocmem);
         BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->conss[c]->vars, (*reoptnode)->conss[c]->allocmem);
         BMSfreeMemory(&(*reoptnode)->conss[c]);
      }
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->conss, (*reoptnode)->allocmemconss);
      (*reoptnode)->nconss = 0;
      (*reoptnode)->allocmemconss = 0;
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
      assert((*reoptnode)->dualconscur->allocmem > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconscur->vals, (*reoptnode)->dualconscur->allocmem);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconscur->vars, (*reoptnode)->dualconscur->allocmem);
      BMSfreeMemory(&(*reoptnode)->dualconscur);
      (*reoptnode)->dualconscur = NULL;
   }

   if( (*reoptnode)->dualconsnex != NULL )
   {
      assert((*reoptnode)->dualconsnex->allocmem > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconsnex->vals, (*reoptnode)->dualconsnex->allocmem);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->dualconsnex->vars, (*reoptnode)->dualconsnex->allocmem);
      BMSfreeMemory(&(*reoptnode)->dualconsnex);
      (*reoptnode)->dualconsnex = NULL;
   }

   /* free boundtypes */
   if ((*reoptnode)->varboundtypes != NULL )
   {
      assert((*reoptnode)->allocvarmem > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->varboundtypes, (*reoptnode)->allocvarmem);
      (*reoptnode)->varboundtypes = NULL;
   }

   /* free bounds */
   if ((*reoptnode)->varbounds != NULL )
   {
      assert((*reoptnode)->allocvarmem > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->varbounds, (*reoptnode)->allocvarmem);
      (*reoptnode)->varbounds = NULL;
   }

   /* free variables */
   if ((*reoptnode)->vars != NULL )
   {
      assert((*reoptnode)->allocvarmem > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->vars, (*reoptnode)->allocvarmem);
      (*reoptnode)->vars = NULL;
   }

   (*reoptnode)->allocvarmem = 0;

   /* free afterdual-boundtypes */
   if ((*reoptnode)->afterdualvarboundtypes != NULL )
   {
      assert((*reoptnode)->allocafterdualvarmem > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->afterdualvarboundtypes, (*reoptnode)->allocafterdualvarmem);
      (*reoptnode)->afterdualvarboundtypes = NULL;
   }

   /* free afterdual-bounds */
   if ((*reoptnode)->afterdualvarbounds != NULL )
   {
      assert((*reoptnode)->allocafterdualvarmem > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->afterdualvarbounds, (*reoptnode)->allocafterdualvarmem);
      (*reoptnode)->afterdualvarbounds = NULL;
   }

   /* free afterdual-variables */
   if ((*reoptnode)->afterdualvars != NULL )
   {
      assert((*reoptnode)->allocafterdualvarmem > 0);
      BMSfreeBlockMemoryArray(blkmem, &(*reoptnode)->afterdualvars, (*reoptnode)->allocafterdualvarmem);
      (*reoptnode)->afterdualvars = NULL;
   }

   (*reoptnode)->allocafterdualvarmem = 0;

   BMSfreeMemory(&(*reoptnode));
   (*reoptnode) = NULL;

   return SCIP_OKAY;
}

/**
 * delete the data for node nodeID in nodedata
 */
static
SCIP_RETCODE reopttreeDeleteNode(
   SCIP_REOPTTREE*       reopttree,               /**< reoptimization tree */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   nodeID,                  /**< id of a node */
   SCIP_Bool             exitsolve                /**< delete at the end of the solving process (hard delete) */
   )
{
   assert(reopttree != NULL );
   assert(reopttree->reoptnodes[nodeID] != NULL );

   if( exitsolve )
   {
      SCIP_CALL( reoptnodeDelete(&reopttree->reoptnodes[nodeID], blkmem) );
   }
   else
   {
      /** remove and delete all constraints */
      if( reopttree->reoptnodes[nodeID]->nconss > 0 )
      {
         int c;

         assert(reopttree->reoptnodes[nodeID]->conss != NULL);
         assert(reopttree->reoptnodes[nodeID]->allocmemconss > 0);

         for(c = 0; c < reopttree->reoptnodes[nodeID]->nconss; c++)
         {
            BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[nodeID]->conss[c]->vals, reopttree->reoptnodes[nodeID]->conss[c]->allocmem);
            BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[nodeID]->conss[c]->vars, reopttree->reoptnodes[nodeID]->conss[c]->allocmem);
            BMSfreeMemory(&reopttree->reoptnodes[nodeID]->conss[c]);
         }
         reopttree->reoptnodes[nodeID]->nconss = 0;
      }

      /* remove all children */
      if (reopttree->reoptnodes[nodeID]->childids != NULL )
      {
         reopttree->reoptnodes[nodeID]->nchilds = 0;
      }

      /* delete dual constraint */
      if( reopttree->reoptnodes[nodeID]->dualconscur != NULL )
      {
         assert(reopttree->reoptnodes[nodeID]->dualconscur->allocmem > 0);
         BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[nodeID]->dualconscur->vals, reopttree->reoptnodes[nodeID]->dualconscur->allocmem);
         BMSfreeBlockMemoryArray(blkmem ,&reopttree->reoptnodes[nodeID]->dualconscur->vars, reopttree->reoptnodes[nodeID]->dualconscur->allocmem);
         BMSfreeMemory(&reopttree->reoptnodes[nodeID]->dualconscur);
         reopttree->reoptnodes[nodeID]->dualconscur = NULL;
      }

      if( reopttree->reoptnodes[nodeID]->dualconsnex != NULL )
      {
         assert(reopttree->reoptnodes[nodeID]->dualconsnex->allocmem > 0);
         BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[nodeID]->dualconsnex->vals, reopttree->reoptnodes[nodeID]->dualconsnex->allocmem);
         BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[nodeID]->dualconsnex->vars, reopttree->reoptnodes[nodeID]->dualconsnex->allocmem);
         BMSfreeMemory(&reopttree->reoptnodes[nodeID]->dualconsnex);
         reopttree->reoptnodes[nodeID]->dualconsnex = NULL;
      }

      reopttree->reoptnodes[nodeID]->nvars = 0;
      reopttree->reoptnodes[nodeID]->nafterdualvars = 0;
      reopttree->reoptnodes[nodeID]->dualfixing = FALSE;
      reopttree->reoptnodes[nodeID]->reopttype = SCIP_REOPTTYPE_NONE;
      reopttree->reoptnodes[nodeID]->lowerbound = -1e+20;
   }

   assert(reopttree->reoptnodes[nodeID] == NULL
        || reopttree->reoptnodes[nodeID]->conss == NULL
        || reopttree->reoptnodes[nodeID]->nconss == 0);
   assert(reopttree->reoptnodes[nodeID] == NULL
       || reopttree->reoptnodes[nodeID]->childids == NULL
       || reopttree->reoptnodes[nodeID]->nchilds == 0);

   reopttree->nsavednodes--;

   return SCIP_OKAY;
}


static
SCIP_RETCODE createSolTree(
   SCIP_SOLTREE*         soltree,                 /**< solution tree */
   BMS_BLKMEM*           blkmem                   /**< block memory */
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
   SCIP_ALLOC( BMSallocMemory(&soltree->root) );
   soltree->root->sol = NULL;
   soltree->root->updated = FALSE;
   soltree->root->father = NULL;
   soltree->root->rchild = NULL;
   soltree->root->lchild = NULL;

   return SCIP_OKAY;
}

static
SCIP_RETCODE soltreefreeNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_PRIMAL*          primal,                  /**< the primal */
   SCIP_SOLNODE*         node,                    /**< node within the solution tree */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(primal != NULL);
   assert(node != NULL);
   assert(blkmem != NULL);

   /* free recursive right subtree */
   if( node->rchild != NULL )
   {
      SCIP_CALL( soltreefreeNode(reopt, set, primal, node->rchild, blkmem) );
   }

   /* free recursive left subtree */
   if( node->lchild != NULL )
   {
      SCIP_CALL( soltreefreeNode(reopt, set, primal, node->lchild, blkmem) );
   }

   if( node->sol != NULL )
   {
      assert(set->stage == SCIP_STAGE_PROBLEM);

      SCIP_CALL( SCIPsolFree(&node->sol, blkmem, primal) );
   }

   /* free this nodes */
   BMSfreeMemory(&node);

   return SCIP_OKAY;
}

/* free solution tree */
static
SCIP_RETCODE freeSolTree(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_PRIMAL*          origprimal,              /**< the origprimal */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(reopt->soltree != NULL);
   assert(reopt->soltree->root != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   /* free all nodes recursive */
   SCIP_CALL( soltreefreeNode(reopt, set, origprimal, reopt->soltree->root, blkmem) );

   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->sols, reopt->runsize);
   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->nsols, reopt->runsize);
   BMSfreeBlockMemoryArray(blkmem, &reopt->soltree->solssize, reopt->runsize);

   BMSfreeMemory(&reopt->soltree);

   return SCIP_OKAY;
}

/* add a node to the solution tree */
static
SCIP_RETCODE soltreeAddNode(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SOLNODE*         father,                  /**< father of the node to add */
   SCIP_Bool             rchild,                  /**< 0-branch? */
   SCIP_Bool             lchild                   /**< 1-branch? */
   )
{
   SCIP_SOLNODE* newnode;

   assert(reopt != NULL);
   assert(father != NULL);
   assert(rchild == !lchild);
   assert((rchild && father->rchild == NULL) || (lchild && father->lchild == NULL));

   SCIP_ALLOC( BMSallocMemory(&newnode) );
   newnode->sol = NULL;
   newnode->updated = FALSE;
   newnode->father = father;
   newnode->rchild = NULL;
   newnode->lchild = NULL;

   if( rchild )
      father->rchild = newnode;
   else
      father->lchild = newnode;

   return SCIP_OKAY;
}

/* add a solution */
static
SCIP_RETCODE soltreeAddSol(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,              /**< origprimal */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_VAR**            vars,                    /**< array of orig vars */
   SCIP_SOL*             sol,                     /**< solution to add */
   SCIP_SOLNODE**        solnode,                 /**< current solution node */
   int                   nvars,                   /**< number of variables */
   SCIP_Bool             bestsol,                 /**< is the solution an optimal (best found) solution */
   SCIP_Bool*            added                    /**< pointer to store the result */
   )
{
   SCIP_SOLNODE* cursolnode;
   int varid;

   assert(reopt != NULL);
   assert(sol != NULL);

   cursolnode = reopt->soltree->root;
   (*added) = FALSE;

   if( set->reopt_savesols > 0 )
   {
      for(varid = 0; varid < nvars; varid++)
      {
         if( SCIPvarGetType(vars[varid]) == SCIP_VARTYPE_BINARY
          || SCIPvarGetType(vars[varid]) == SCIP_VARTYPE_INTEGER
          || SCIPvarGetType(vars[varid]) == SCIP_VARTYPE_IMPLINT )
         {
            SCIP_Real objval;

            objval = SCIPsolGetVal(sol, set, stat, vars[varid]);
            if( SCIPsetIsFeasEQ(set, objval, 0) )
            {
               if( cursolnode->rchild == NULL )
               {
                  SCIP_CALL( soltreeAddNode(reopt, cursolnode, TRUE, FALSE) );
                  assert(cursolnode->rchild != NULL);
                  (*added) = TRUE;
               }
               cursolnode = cursolnode->rchild;
            }
            else
            {
               assert(SCIPsetIsFeasEQ(set, objval, 1));
               if( cursolnode->lchild == NULL )
               {
                  SCIP_CALL( soltreeAddNode(reopt, cursolnode, FALSE, TRUE) );
                  assert(cursolnode->lchild != NULL);
                  (*added) = TRUE;
               }
               cursolnode = cursolnode->lchild;
            }
         }
      }

      /* the solution was added */
      if( *added )
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
            assert(reopt->lastbestsol != NULL);
            assert(cursolnode->sol != NULL);

            reopt->lastbestsol[reopt->run-1] = cursolnode->sol;
            printf(">> add best sol at pos %d\n", reopt->run-1);
         }

         (*solnode) = cursolnode;
      }
      else if( bestsol )
      {
         SCIP_SOL* copysol;
         SCIP_CALL( SCIPsolCopy(&copysol, blkmem, set, stat, origprimal, sol) );
         reopt->lastbestsol[reopt->run-1] = copysol;
         printf(">> add best sol at pos %d\n", reopt->run-1);
      }
   }
   else if( bestsol )
   {
      SCIP_SOL* copysol;
      SCIP_CALL( SCIPsolCopy(&copysol, blkmem, set, stat, origprimal, sol) );
      reopt->lastbestsol[reopt->run-1] = copysol;
      printf(">> add best sol at pos %d\n", reopt->run-1);
   }

   return SCIP_OKAY;
}

/* set all marks updated to FALSE */
static
void soltreeResetMarks(
   SCIP_SOLNODE*         node                     /**< node within the solution tree */
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

/* return the number of used solutions */
static
int soltreeGetNUsedSols(
   SCIP_SOLNODE*         node                     /**< node within the solution tree */
   )
{
   int nusedsols;

   assert(node != NULL);

   nusedsols = 0;

   if(node->lchild != NULL)
      nusedsols += soltreeGetNUsedSols(node->lchild);
   if(node->rchild != NULL)
      nusedsols += soltreeGetNUsedSols(node->rchild);
   if(node->rchild == NULL && node->lchild == NULL)
      nusedsols = 1;

   return nusedsols;
}

/*
 * allocate memory for a node within the reopttree
 */
static
SCIP_RETCODE createReoptnode(
   SCIP_REOPTTREE*       reopttree,               /**< reoptimization tree */
   int                   nodeID                   /**< id of the node to create */
   )
{
   assert(reopttree != NULL );
   assert(0 <= nodeID && nodeID < reopttree->allocmemnodes);

   SCIPdebugMessage("create a reoptnode at ID %d\n", nodeID);

   if(reopttree->reoptnodes[nodeID] == NULL )
   {
      SCIP_ALLOC( BMSallocMemory(&reopttree->reoptnodes[nodeID]) );
      reopttree->reoptnodes[nodeID]->conss = NULL;
      reopttree->reoptnodes[nodeID]->nconss = 0;
      reopttree->reoptnodes[nodeID]->allocmemconss = 0;
      reopttree->reoptnodes[nodeID]->lpistate = NULL;
      reopttree->reoptnodes[nodeID]->childids = NULL;
      reopttree->reoptnodes[nodeID]->allocchildmem = 0;
      reopttree->reoptnodes[nodeID]->nchilds = 0;
      reopttree->reoptnodes[nodeID]->nvars = 0;
      reopttree->reoptnodes[nodeID]->nafterdualvars = 0;
      reopttree->reoptnodes[nodeID]->parentID = -1;
      reopttree->reoptnodes[nodeID]->dualfixing = FALSE;
      reopttree->reoptnodes[nodeID]->reopttype = SCIP_REOPTTYPE_NONE;
      reopttree->reoptnodes[nodeID]->allocvarmem = 0;
      reopttree->reoptnodes[nodeID]->allocafterdualvarmem = 0;
      reopttree->reoptnodes[nodeID]->vars = NULL;
      reopttree->reoptnodes[nodeID]->varbounds = NULL;
      reopttree->reoptnodes[nodeID]->varboundtypes = NULL;
      reopttree->reoptnodes[nodeID]->afterdualvars = NULL;
      reopttree->reoptnodes[nodeID]->afterdualvarbounds = NULL;
      reopttree->reoptnodes[nodeID]->afterdualvarboundtypes = NULL;
      reopttree->reoptnodes[nodeID]->dualconscur = NULL;
      reopttree->reoptnodes[nodeID]->dualconsnex = NULL;
      reopttree->reoptnodes[nodeID]->lowerbound = -1e+20;
   }
   else
   {
      assert(reopttree->reoptnodes[nodeID]->nvars == 0);
      assert(reopttree->reoptnodes[nodeID]->nafterdualvars == 0);
      reopttree->reoptnodes[nodeID]->reopttype = SCIP_REOPTTYPE_NONE;
      reopttree->reoptnodes[nodeID]->lowerbound = -1e+20;
   }

   /* increase the counter */
   reopttree->nsavednodes++;

   return SCIP_OKAY;
}

/* create the reopttree */
static
SCIP_RETCODE createReopttree(
   SCIP_REOPTTREE*       reopttree,               /**< pointer to reopttree */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   int id;

   assert(reopttree != NULL);

   /* allocate memory */
   reopttree->allocmemnodes = DEFAULT_MEM_NODES;
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopttree->reoptnodes, reopttree->allocmemnodes) );

   /* initialize the queue of open IDs */
   SCIP_CALL( SCIPqueueCreate(&reopttree->openids, reopttree->allocmemnodes, 2) );

   /* fill the queue, but reserve the 0 for the root */
   for(id = 1; id < reopttree->allocmemnodes; id++)
   {
      reopttree->reoptnodes[id] = NULL;
      SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void*) (size_t) id) );
   }
   assert(SCIPqueueNElems(reopttree->openids) == reopttree->allocmemnodes-1);

   /* initialize the root node */
   reopttree->reoptnodes[0] = NULL;
   SCIP_CALL( createReoptnode(reopttree, 0) );

   reopttree->nsavednodes = 0;
   reopttree->nbranchednodes = 0;
   reopttree->nbranchednodesround = 0;
   reopttree->nfeasnodes = 0;
   reopttree->nfeasnodesround = 0;
   reopttree->ninfeasnodes = 0;
   reopttree->ninfeasnodesround = 0;
   reopttree->nprunednodes = 0;
   reopttree->nprunednodesround = 0;
   reopttree->nrediednodes = 0;
   reopttree->nrediednodesround = 0;

   return SCIP_OKAY;
}

/*
 * Clear the reopttree, e.g., to restart and solve the next problem from scratch
 */
static
SCIP_RETCODE clearReoptnodes(
   SCIP_REOPTTREE*       reopttree,               /**< reoptimization tree */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_Bool             exitsolve                /**< delete nodes before exit the solving process */
   )
{
   int id;

   assert(reopttree != NULL );

   /** clear queue with open IDs */
   SCIPqueueClear(reopttree->openids);
   assert(SCIPqueueNElems(reopttree->openids) == 0);

   /** delete all data about nodes */
   for(id = 0; id < reopttree->allocmemnodes; id++)
   {
      if(reopttree->reoptnodes[id] != NULL )
      {
         SCIP_CALL( reopttreeDeleteNode(reopttree, blkmem, id, exitsolve) );
         assert(reopttree->reoptnodes[id] == NULL || reopttree->reoptnodes[id]->nvars == 0);
      }

      if(id > 0 && !exitsolve)
      {
         SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void* ) (size_t ) id) );
      }
   }
   assert(exitsolve || SCIPqueueNElems(reopttree->openids) == reopttree->allocmemnodes-1);

   reopttree->nsavednodes = 0;

   return SCIP_OKAY;
}

/* free the reopttree */
static
SCIP_RETCODE freeReoptTree(
   SCIP_REOPTTREE*       reopttree,               /**< tree data */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   assert(reopttree != NULL);
   assert(blkmem != NULL);

   /* free nodes */
   SCIP_CALL( clearReoptnodes(reopttree, blkmem, TRUE) );

   /* free the data */
   BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes, reopttree->allocmemnodes);
   SCIPqueueFree(&reopttree->openids);

   /* free the tree itself */
   BMSfreeMemory(&reopttree);

   return SCIP_OKAY;
}

/* check memory for the constraint to handle bound changes based on dual information */
static
SCIP_RETCODE checkMemDualCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   size                     /**< size which need to be allocated */
   )
{
   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(size > 0);

   if( reopt->dualcons == NULL )
   {
      SCIP_ALLOC( BMSallocMemory(&reopt->dualcons) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->vars, size) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->vals, size) );
      reopt->dualcons->allocmem = size;
      reopt->dualcons->nvars = 0;
   }
   else if( reopt->dualcons->allocmem < size )
   {
      if( reopt->dualcons->allocmem > 0 )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->dualcons->vars, reopt->dualcons->allocmem, size) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->dualcons->vals, reopt->dualcons->allocmem, size) );
      }
      else
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->vars, size) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->dualcons->vals, size) );
         reopt->dualcons->nvars = 0;
      }

      reopt->dualcons->allocmem = size;
   }

   return SCIP_OKAY;
}

/* check the memory to store global constraints */
static
SCIP_RETCODE checkMemGlbCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   mem                      /**< memory which has to be allocated */
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

/* update the constraint propagations made in the current iteration;
 * stop saving the bound changes if we reach a branching decision based on a
 * dual information
 * */
static
SCIP_RETCODE updateConstraintPropagation(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   int                   nodeID,                  /**< id of the node */
   SCIP_Bool*            transintoorig            /**< transform variables into originals */
   )
{
   int nvars;
   int nconsprops;
   int naddedbndchgs;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(0 < nodeID && nodeID < reopt->reopttree->allocmemnodes);
   assert(reopt->reopttree->reoptnodes[nodeID] != NULL );

   /* get the number of all stored constraint propagations */
   nconsprops = SCIPnodeGetNDomchg(node, FALSE, TRUE, FALSE);
   nvars = reopt->reopttree->reoptnodes[nodeID]->nvars;

   if( nconsprops > 0 )
   {
      /* check the memory */
      SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[nodeID], blkmem, nvars + nconsprops, 0, 0) );

      SCIPnodeGetConsProps(node,
            &reopt->reopttree->reoptnodes[nodeID]->vars[nvars],
            &reopt->reopttree->reoptnodes[nodeID]->varbounds[nvars],
            &reopt->reopttree->reoptnodes[nodeID]->varboundtypes[nvars],
            &naddedbndchgs,
            reopt->reopttree->reoptnodes[nodeID]->allocvarmem-nvars);

      assert(nvars + naddedbndchgs <= reopt->reopttree->reoptnodes[nodeID]->allocvarmem);

      reopt->reopttree->reoptnodes[nodeID]->nvars += naddedbndchgs;

      *transintoorig = TRUE;
   }

   return SCIP_OKAY;
}

/*
 * save bound changes made after dual methods, e.g., strong branching.
 */
static
SCIP_RETCODE saveAfterDualBranchings(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   int                   nodeID,                  /**< id of the node */
   SCIP_Bool*            transintoorig            /**< transform variables into originals */
   )
{
   int nbranchvars;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(0 < nodeID && nodeID < reopt->reopttree->allocmemnodes);
   assert(reopt->reopttree->reoptnodes[nodeID] != NULL );

   nbranchvars = 0;

   /* allocate memory */
   if (reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem == 0)
   {
      assert(reopt->reopttree->reoptnodes[nodeID]->afterdualvars == NULL );
      assert(reopt->reopttree->reoptnodes[nodeID]->afterdualvarbounds == NULL );
      assert(reopt->reopttree->reoptnodes[nodeID]->afterdualvarboundtypes == NULL );

      /** allocate block memory for node information */
      reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem = DEFAULT_MEM_VARAFTERDUAL;
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[nodeID]->afterdualvars), reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[nodeID]->afterdualvarbounds), reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[nodeID]->afterdualvarboundtypes), reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem) );
   }

   assert(reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem > 0);
   assert(reopt->reopttree->reoptnodes[nodeID]->nafterdualvars >= 0);

   SCIPnodeGetAfterDualBranchingsReopt(node,
         reopt->reopttree->reoptnodes[nodeID]->afterdualvars,
         reopt->reopttree->reoptnodes[nodeID]->afterdualvarbounds,
         reopt->reopttree->reoptnodes[nodeID]->afterdualvarboundtypes,
         reopt->reopttree->reoptnodes[nodeID]->nafterdualvars,
         &nbranchvars,
         reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem);

   if( nbranchvars > reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem )
   {
      int newsize;
      newsize = nbranchvars + 1;
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[nodeID]->afterdualvars), reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[nodeID]->afterdualvarbounds), reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &(reopt->reopttree->reoptnodes[nodeID]->afterdualvarboundtypes), reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem, newsize) );
      reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem = newsize;

      SCIPnodeGetAfterDualBranchingsReopt(node,
            reopt->reopttree->reoptnodes[nodeID]->afterdualvars,
            reopt->reopttree->reoptnodes[nodeID]->afterdualvarbounds,
            reopt->reopttree->reoptnodes[nodeID]->afterdualvarboundtypes,
            reopt->reopttree->reoptnodes[nodeID]->nafterdualvars,
            &nbranchvars,
            reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem);
   }

   /* the stored variables of this node need to be transformed into the original space */
   if( nbranchvars > 0 )
      *transintoorig = TRUE;

   SCIPdebugMessage(" -> save %d bound changes after dual reductions\n", nbranchvars);

   assert(nbranchvars <= reopt->reopttree->reoptnodes[nodeID]->allocafterdualvarmem); /* this should be the case */

   reopt->reopttree->reoptnodes[nodeID]->nafterdualvars = nbranchvars;

   return SCIP_OKAY;
}

/**
 * transform variable and bounds back to the originals
 */
static
SCIP_RETCODE transformIntoOrig(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   nodeID                   /**< id of the node */
   )
{
   int varnr;

   assert(reopt != NULL );
   assert(nodeID >= 1);
   assert(reopt->reopttree->reoptnodes[nodeID] != NULL );

   /* transform branching variables and bound changes applied before the first dual reduction */
   for(varnr = 0; varnr < reopt->reopttree->reoptnodes[nodeID]->nvars; varnr++)
   {
      SCIP_Real constant;
      SCIP_Real scalar;

      scalar = 1;
      constant = 0;

      if(!SCIPvarIsOriginal(reopt->reopttree->reoptnodes[nodeID]->vars[varnr]))
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->reopttree->reoptnodes[nodeID]->vars[varnr], &scalar, &constant)) ;
         reopt->reopttree->reoptnodes[nodeID]->varbounds[varnr] = (reopt->reopttree->reoptnodes[nodeID]->varbounds[varnr] - constant) / scalar;
      }
      assert(SCIPvarIsOriginal(reopt->reopttree->reoptnodes[nodeID]->vars[varnr]));
   }

   /* transform bound changes affected by dual reduction */
   for(varnr = 0; varnr < reopt->reopttree->reoptnodes[nodeID]->nafterdualvars; varnr++)
   {
      SCIP_Real constant;
      SCIP_Real scalar;

      scalar = 1;
      constant = 0;

      if(!SCIPvarIsOriginal(reopt->reopttree->reoptnodes[nodeID]->afterdualvars[varnr]))
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->reopttree->reoptnodes[nodeID]->afterdualvars[varnr], &scalar, &constant)) ;
         reopt->reopttree->reoptnodes[nodeID]->afterdualvarbounds[varnr] = (reopt->reopttree->reoptnodes[nodeID]->afterdualvarbounds[varnr] - constant) / scalar;
      }
      assert(SCIPvarIsOriginal(reopt->reopttree->reoptnodes[nodeID]->afterdualvars[varnr]));
   }

   return SCIP_OKAY;
}

/*
 * search the next node along the root path that
 * is saved by reoptimization
 */
static
SCIP_RETCODE getLastSavedNode(
   SCIP_REOPT*           reopt,                   /**< reotimization data structure */
   SCIP_NODE*            node,                    /**< node of the search tree */
   SCIP_NODE**           parent,                  /**< parent node within the search tree */
   int*                  parentID,                /**< id of the parent node */
   int*                  nbndchgs                 /**< number of bound changes */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(reopt->reopttree->reoptnodes != NULL);

   (*nbndchgs) = 0;
   (*parent) = node;

   /** look for a saved parent along the root-path */
   while( SCIPnodeGetDepth(*parent) != 0 )
   {
      (*nbndchgs) += SCIPnodeGetNDomchg((*parent), TRUE, TRUE, FALSE);
      (*parent) = SCIPnodeGetParent(*parent);

      if( SCIPnodeGetDepth(*parent) == 0)
      {
         (*parentID) = 0;
         break;
      }
      else if( SCIPnodeGetReopttype((*parent)) >= SCIP_REOPTTYPE_TRANSIT )
      {
         assert(SCIPnodeGetReoptID((*parent)) < reopt->reopttree->allocmemnodes);
         (*parentID) = SCIPnodeGetReoptID((*parent));
         break;
      }
   }

   return SCIP_OKAY;
}

/* returns the number of bound changes along the root path up to the
 * next stored node */
static
int lengthBranchPath(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node                     /**< node of the search tree */
   )
{
   int length;

   assert(reopt != NULL);
   assert(node != 0);

   if( SCIPnodeGetDepth(node) == 0 )
      return 0;
   else
   {
      SCIP_NODE* parent;

      assert(SCIPnodeGetReoptID(node) >= 0);
      assert(reopt->reopttree->reoptnodes[SCIPnodeGetReoptID(node)] != NULL);

      parent = SCIPnodeGetParent(node);
      length = reopt->reopttree->reoptnodes[SCIPnodeGetReoptID(node)]->nvars;
      while(SCIPnodeGetDepth(parent) != 0)
      {
         if( SCIPnodeGetReopttype(parent) >= SCIP_REOPTTYPE_TRANSIT )
         {
            assert(reopt->reopttree->reoptnodes[SCIPnodeGetReoptID(parent)] != NULL);
            length += reopt->reopttree->reoptnodes[SCIPnodeGetReoptID(parent)]->nvars;
         }
         parent = SCIPnodeGetParent(parent);
      }
   }

   return length;
}

/* adds the id @param childid to the array of child nodes of @param parentid */
static
SCIP_RETCODE reoptAddChild(
   SCIP_REOPTTREE*       reopttree,               /**< reoptimization tree */
   int                   parentid,                /**< id of the parent node */
   int                   childid,                 /**< id of the child node */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   int nchilds;

   assert(reopttree != NULL);
   assert(blkmem != NULL);
   assert(0 <= parentid && parentid < reopttree->allocmemnodes);
   assert(0 <= childid && childid < reopttree->allocmemnodes);
   assert(reopttree->reoptnodes[parentid] != NULL);

   nchilds = reopttree->reoptnodes[parentid]->nchilds;

   /* ensure that the array is large enough */
   SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[parentid], blkmem, 0, nchilds+1, 0) );
   assert(reopttree->reoptnodes[parentid]->allocchildmem > nchilds);

   /* add the child */
   reopttree->reoptnodes[parentid]->childids[nchilds] = childid;
   reopttree->reoptnodes[parentid]->nchilds++;

   SCIPdebugMessage("add ID %d as a child of ID %d.\n", childid, parentid);

   return SCIP_OKAY;
}

/*
 * move all children to the next node stored by reoptimization.
 */
static
SCIP_RETCODE moveChildrenUp(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   nodeID,                  /**< id of the node */
   int                   parentID                 /**< id of the parent node */
   )
{
   int childID;
   int varnr;
   int nvars;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(nodeID >= 1);
   assert(parentID >= 0);
   assert(reopt->reopttree->reoptnodes[nodeID]->childids != NULL);

   /* ensure that enough memory at the parentID is available */
   SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[parentID], blkmem, 0, reopt->reopttree->reoptnodes[parentID]->nchilds + reopt->reopttree->reoptnodes[nodeID]->nchilds, 0) );

   while( reopt->reopttree->reoptnodes[nodeID]->nchilds > 0 )
   {
      int nchilds;

      nchilds = reopt->reopttree->reoptnodes[nodeID]->nchilds;
      childID = reopt->reopttree->reoptnodes[nodeID]->childids[nchilds-1];

      /* check the memory */
      SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[childID], blkmem, reopt->reopttree->reoptnodes[childID]->nvars + reopt->reopttree->reoptnodes[nodeID]->nvars, 0, 0) );
      assert(reopt->reopttree->reoptnodes[childID]->allocvarmem >= reopt->reopttree->reoptnodes[childID]->nvars + reopt->reopttree->reoptnodes[nodeID]->nvars);

      /** save branching information */
      for(varnr = 0; varnr < reopt->reopttree->reoptnodes[nodeID]->nvars; varnr++)
      {
         nvars = reopt->reopttree->reoptnodes[childID]->nvars;
         reopt->reopttree->reoptnodes[childID]->vars[nvars] = reopt->reopttree->reoptnodes[nodeID]->vars[varnr];
         reopt->reopttree->reoptnodes[childID]->varbounds[nvars] = reopt->reopttree->reoptnodes[nodeID]->varbounds[varnr];
         reopt->reopttree->reoptnodes[childID]->varboundtypes[nvars] = reopt->reopttree->reoptnodes[nodeID]->varboundtypes[varnr];
         reopt->reopttree->reoptnodes[childID]->nvars++;
      }

      /* update the ID of the parent node */
      reopt->reopttree->reoptnodes[childID]->parentID = parentID;

      /* insert the node as a child */
      SCIP_CALL( reoptAddChild(reopt->reopttree, parentID, childID, blkmem) );

      /* reduce the number of child nodes by 1 */
      reopt->reopttree->reoptnodes[nodeID]->nchilds--;
   }

   return SCIP_OKAY;
}

/* delete all nodes in the subtree induced by nodeID */
static
SCIP_RETCODE deleteChildrenBelow(
   SCIP_REOPTTREE*       reopttree,               /**< reoptimization tree */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   nodeID,                  /**< id of the node */
   SCIP_Bool             delnodeitself,           /**< should the node deleted after deleting the induced subtree? */
   SCIP_Bool             exitsolve                /**< will the solving process end after deletion */
   )
{
   assert(reopttree != NULL );
   assert(blkmem != NULL);
   assert(nodeID >= 0);
   assert(reopttree->reoptnodes[nodeID] != NULL);

   /** delete all children below */
   if( reopttree->reoptnodes[nodeID]->childids != NULL && reopttree->reoptnodes[nodeID]->nchilds > 0 )
   {
      SCIPdebugMessage("-> delete subtree induced by ID %d (hard remove = %u)\n", nodeID, exitsolve);

      while( reopttree->reoptnodes[nodeID]->nchilds > 0 )
      {
         int nchilds;
         int childID;

         nchilds = reopttree->reoptnodes[nodeID]->nchilds;
         childID = reopttree->reoptnodes[nodeID]->childids[nchilds-1];

         SCIP_CALL( deleteChildrenBelow(reopttree, blkmem, childID, TRUE, exitsolve) );

         reopttree->reoptnodes[nodeID]->nchilds--;
      }
   }

   /** delete node data*/
   if( delnodeitself )
   {
      SCIP_CALL( reopttreeDeleteNode(reopttree, blkmem, nodeID, exitsolve) );
      SCIP_CALL( SCIPqueueInsert(reopttree->openids, (void*) (size_t) nodeID) );
   }

   return SCIP_OKAY;
}

/*
 * replace transit nodes by stored child nodes
 */
static
SCIP_RETCODE shrinkNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node,                    /**< node of the search tree */
   int                   nodeID,                  /**< id of the node */
   SCIP_Bool*            shrank,                  /**< pointer to store if the node was shrank */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   SCIP_NODE* parent;
   int ndomchgs;
   int parentID;

   assert(reopt != NULL);
   assert(node != NULL);
   assert(reopt->reopttree->reoptnodes[nodeID] != NULL);

   if( reopt->reopttree->reoptnodes[nodeID]->childids != NULL
    && reopt->reopttree->reoptnodes[nodeID]->nchilds > 0 )
   {
      ndomchgs = 0;
      parentID = -1;
      parent = NULL;

      SCIP_CALL( getLastSavedNode(reopt, node, &parent, &parentID, &ndomchgs) );

      assert(parentID != nodeID);
      assert(reopt->reopttree->reoptnodes[parentID] != NULL );
      assert(reopt->reopttree->reoptnodes[parentID]->childids != NULL && reopt->reopttree->reoptnodes[parentID]->nchilds);

      /* check if we want move all children to the next saved node above
       * we want to shrink the path if either
       * - the maximal number of bound changes fix and the number of bound changes is
       *   less than the given threshold reopt->maxdiffofnodes
       * or
       * - the number is calculated dynamically and the number of bound changes
       *   is less than log2(SCIPgetNBinVars - (#vars of parent))
       * */
      if( (!reopt->dynamicdiffofnodes && ndomchgs <= reopt->maxdiffofnodes) )
        /* ||(reopt->dynamicdiffofnodes && ndomchgs <= ceil(log10((SCIP_Real) (SCIPgetNOrigBinVars(scip) - MIN(SCIPgetNOrigBinVars(scip)-1, lengthBranchPath(reopt, parent))))/log10(2.0))) ) */
      {
         int c;

         SCIPdebugMessage(" -> shrink node %lld at ID %d, replaced by %d child nodes.\n", SCIPnodeGetNumber(node), nodeID, reopt->reopttree->reoptnodes[nodeID]->nchilds);

         /* copy the references of child nodes to the parent*/
         SCIP_CALL( moveChildrenUp(reopt, blkmem, nodeID, parentID) );

         /* delete the current node */
         c = 0;
         while( reopt->reopttree->reoptnodes[parentID]->childids[c] != nodeID
             && c < reopt->reopttree->reoptnodes[parentID]->nchilds )
         {
            c++;
         }
         assert(reopt->reopttree->reoptnodes[parentID]->childids[c] == nodeID);

         /* replace the childid at position c by the last one */
         reopt->reopttree->reoptnodes[parentID]->childids[c] = reopt->reopttree->reoptnodes[parentID]->childids[reopt->reopttree->reoptnodes[parentID]->nchilds-1];
         reopt->reopttree->reoptnodes[parentID]->nchilds--;

         SCIP_CALL( reopttreeDeleteNode(reopt->reopttree, blkmem, nodeID, TRUE) );
         SCIP_CALL( SCIPqueueInsert(reopt->reopttree->openids, (void*) (size_t) nodeID) );

         *shrank = TRUE;
      }
   }

   return SCIP_OKAY;
}

/*
 * change the reopttype of the subtree induced by nodeID
 */
static
SCIP_RETCODE changeReopttypeOfSubtree(
   SCIP_REOPTTREE*       reopttree,               /**< reopttree */
   int                   nodeID,                  /**< id of the node */
   SCIP_REOPTTYPE        reopttype                /**< reopttype */
   )
{
   assert(reopttree != NULL);
   assert(nodeID >= 0);
   assert(reopttree->reoptnodes[nodeID] != NULL);

   if( reopttree->reoptnodes[nodeID]->childids != NULL && reopttree->reoptnodes[nodeID]->nchilds > 0 )
   {
      int childID;
      int nchildIDs;
      int seenIDs;

      nchildIDs = reopttree->reoptnodes[nodeID]->nchilds;
      seenIDs = 0;

      while( seenIDs < nchildIDs )
      {
         /* get childID */
         childID = reopttree->reoptnodes[nodeID]->childids[seenIDs];
         assert(reopttree->reoptnodes[childID] != NULL);

         /* change the reopttype of the node iff the node is neither infeasible nor indices an
          * infeasible subtree and if the node contains no bound changes based on dual decisions */
         if( reopttree->reoptnodes[childID]->reopttype != SCIP_REOPTTYPE_STRBRANCHED
          && reopttree->reoptnodes[childID]->reopttype != SCIP_REOPTTYPE_INFSUBTREE
          && reopttree->reoptnodes[childID]->reopttype != SCIP_REOPTTYPE_INFEASIBLE )
            reopttree->reoptnodes[childID]->reopttype = reopttype;

         /* change reopttype of subtree */
         SCIP_CALL( changeReopttypeOfSubtree(reopttree, childID, reopttype) );

         seenIDs++;
      }
   }

   return SCIP_OKAY;
}

/**
 * save ancestor branching information up to the
 * next stored node
 */
static
SCIP_RETCODE saveAncestorBranchings(
   SCIP_REOPTTREE*       reopttree,               /**< reopttree */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node */
   SCIP_NODE*            parent,                  /**< parent node */
   int                   nodeID,                  /**< id of the node */
   int                   parentID                 /**< id of the parent node */
   )
{
   int nbranchvars;

   assert(reopttree != NULL );
   assert(node != NULL );
   assert(parent != NULL );
   assert(nodeID >= 1 && nodeID < reopttree->allocmemnodes);
   assert(reopttree->reoptnodes[nodeID] != NULL );
   assert(parentID == 0 || reopttree->reoptnodes[parentID] != NULL ); /* if the root is the next saved node, the nodedata can be NULL */

   SCIPdebugMessage(" -> save ancestor branchings\n");

   /* allocate memory */
   if (reopttree->reoptnodes[nodeID]->allocvarmem == 0)
   {
      assert(reopttree->reoptnodes[nodeID]->vars == NULL );
      assert(reopttree->reoptnodes[nodeID]->varbounds == NULL );
      assert(reopttree->reoptnodes[nodeID]->varboundtypes == NULL );

      /** allocate memory for node information */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[nodeID], blkmem, DEFAULT_MEM_VAR, 0, 0) );
   }

   assert(reopttree->reoptnodes[nodeID]->allocvarmem > 0);
   assert(reopttree->reoptnodes[nodeID]->nvars == 0);

   SCIPnodeGetAncestorBranchingsReopt(node, parent,
         reopttree->reoptnodes[nodeID]->vars,
         reopttree->reoptnodes[nodeID]->varbounds,
         reopttree->reoptnodes[nodeID]->varboundtypes,
         &nbranchvars,
         reopttree->reoptnodes[nodeID]->allocvarmem);

   if( nbranchvars >  reopttree->reoptnodes[nodeID]->allocvarmem )
   {
      /* reallocate memory */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[nodeID], blkmem, nbranchvars, 0, 0) );

      SCIPnodeGetAncestorBranchingsReopt(node, parent,
            reopttree->reoptnodes[nodeID]->vars,
            reopttree->reoptnodes[nodeID]->varbounds,
            reopttree->reoptnodes[nodeID]->varboundtypes,
            &nbranchvars,
            reopttree->reoptnodes[nodeID]->allocvarmem);
   }

   assert(nbranchvars <= reopttree->reoptnodes[nodeID]->allocvarmem); /* this should be the case */

   reopttree->reoptnodes[nodeID]->nvars = nbranchvars;

   assert(nbranchvars <= reopttree->reoptnodes[nodeID]->allocvarmem);
   assert(reopttree->reoptnodes[nodeID]->vars != NULL );

   return SCIP_OKAY;
}

static
SCIP_RETCODE saveLocalConssData(
   SCIP_REOPTTREE*       reopttree,               /**< reopttree */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_NODE*            node,                    /**< node of teh search tree */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   nodeID                   /**< id of the node*/
   )
{
   SCIP_CONS** addedcons;
   SCIP_Real constant;
   SCIP_Real scalar;
   int var;
   int consnr;
   int naddedcons;
   int nconss;

   assert(node != NULL );
   assert(reopttree != NULL);

   /** save the added pseudo-constraint */
   if(SCIPnodeGetNAddedcons(node) > 0)
   {
      naddedcons = SCIPnodeGetNAddedcons(node);

      SCIPdebugMessage(" -> save %d locally added constraints\n", naddedcons);

      /** get memory */
      SCIP_ALLOC( BMSallocMemoryArray(&addedcons, naddedcons) );
      SCIP_CALL( SCIPnodeGetAddedcons(node, addedcons) );

      nconss = reopttree->reoptnodes[nodeID]->nconss;

      /* check memory for added constraints */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[nodeID], blkmem, 0, 0, nconss+naddedcons) );

      for(consnr = 0; consnr < naddedcons; consnr++)
      {
         SCIP_Bool success;

         SCIP_CALL( SCIPallocMemory(scip, &reopttree->reoptnodes[nodeID]->conss[nconss]) );

         success = FALSE;
         SCIP_CALL( SCIPconsGetNVars(addedcons[consnr], set, &reopttree->reoptnodes[nodeID]->conss[nconss]->nvars, &success) );
         assert(success);
         reopttree->reoptnodes[nodeID]->conss[nconss]->allocmem = reopttree->reoptnodes[nodeID]->conss[nconss]->nvars;

         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopttree->reoptnodes[nodeID]->conss[nconss]->vars, reopttree->reoptnodes[nodeID]->conss[nconss]->nvars) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopttree->reoptnodes[nodeID]->conss[nconss]->vals, reopttree->reoptnodes[nodeID]->conss[nconss]->nvars) );

         success = FALSE;
         SCIPconsGetVars(addedcons[consnr], set, reopttree->reoptnodes[nodeID]->conss[nconss]->vars, reopttree->reoptnodes[nodeID]->conss[nconss]->nvars, &success);
         assert(success);

         if( strcmp("sepasol", SCIPconsGetName(addedcons[consnr])) == 0 )
            reopttree->reoptnodes[nodeID]->conss[nconss]->constype = REOPT_CONSTYPE_SEPASOLUTION;
         else if( strcmp("infsubtree", SCIPconsGetName(addedcons[consnr])) == 0 )
            reopttree->reoptnodes[nodeID]->conss[nconss]->constype = REOPT_CONSTYPE_INFSUBTREE;
         else if( strcmp("splitcons", SCIPconsGetName(addedcons[consnr])) == 0 )
            reopttree->reoptnodes[nodeID]->conss[nconss]->constype = REOPT_CONSTYPE_STRBRANCHED;

         assert(reopttree->reoptnodes[nodeID]->conss[nconss]->constype == REOPT_CONSTYPE_SEPASOLUTION
             || reopttree->reoptnodes[nodeID]->conss[nconss]->constype == REOPT_CONSTYPE_INFSUBTREE
             || reopttree->reoptnodes[nodeID]->conss[nconss]->constype == REOPT_CONSTYPE_STRBRANCHED);

         for(var = 0; var < reopttree->reoptnodes[nodeID]->conss[nconss]->nvars; var++)
         {
            constant = 0;
            scalar = 1;

            if(!SCIPvarIsOriginal(reopttree->reoptnodes[nodeID]->conss[nconss]->vars[var]))
            {
               if(SCIPvarIsNegated(reopttree->reoptnodes[nodeID]->conss[nconss]->vars[var]))
               {
                  SCIP_CALL(SCIPvarGetOrigvarSum(&reopttree->reoptnodes[nodeID]->conss[nconss]->vars[var], &scalar, &constant));
                  reopttree->reoptnodes[nodeID]->conss[nconss]->vals[var] = 1;
               }
               else
               {
                  SCIP_CALL(SCIPvarGetOrigvarSum(&reopttree->reoptnodes[nodeID]->conss[nconss]->vars[var], &scalar, &constant));
                  reopttree->reoptnodes[nodeID]->conss[nconss]->vals[var] = 0;
               }
               assert(reopttree->reoptnodes[nodeID]->conss[nconss]->vars[var] != NULL );
            }
            assert(SCIPvarIsOriginal(reopttree->reoptnodes[nodeID]->conss[nconss]->vars[var]));
         }

         /* increase the counter for added constraints */
         reopttree->reoptnodes[nodeID]->nconss++;
         nconss++;
      }

      assert(reopttree->reoptnodes[nodeID]->nconss == naddedcons);
      SCIPfreeMemoryArray(scip, &addedcons);
   }

   return SCIP_OKAY;
}

/**
 * save the LPI state
 */
static
SCIP_RETCODE saveLPIstate(
   SCIP_REOPTTREE*       reopttree,               /**< reopttree */
   SCIP_NODE*            node,                    /**< node of the search tree */
   int                   nodeID                   /**< id of the node */
   )
{
   printf("TODO: implement saveLPIstate\n");

   return SCIP_ERROR;
}

/* collect all bound changes based on dual information
 *
 * if the bound changes are global, all information are already stored because
 * they were caught by an event handler. otherwise, we need to use
 * SCIPnodeGetPseudoBranchings.
 *
 * afterwards, we check if the constraint will be added next or after
 * splitting the node.
 */
static
SCIP_RETCODE collectDualInformation(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node,                    /**< node of the search tree */
   int                   id,                      /**< id of the node */
   SCIP_REOPTTYPE        reopttype,               /**< reopttype */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   SCIP_Real constant;
   SCIP_Real scalar;
   SCIP_Bool cons_is_next;
   int nbndchgs;
   int v;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(0 <= id && id < reopt->reopttree->allocmemnodes);
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
      SCIPnodeGetPseudoBranchings(node,
            reopt->dualcons->vars,
            reopt->dualcons->vals,
            &nbndchgs,
            reopt->dualcons->allocmem);

      assert(nbndchgs <= reopt->dualcons->allocmem);

      reopt->dualcons->nvars = nbndchgs;

      /* transform the variables into the original space */
      for(v = 0; v < nbndchgs; v++)
      {
         constant = 0;
         scalar = 1;

         SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->dualcons->vars[v], &scalar, &constant) );
         reopt->dualcons->vals[v] = (reopt->dualcons->vals[v] - constant) / scalar;

         assert(SCIPvarIsOriginal(reopt->dualcons->vars[v]));
      }
   }

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
      SCIP_ALLOC( BMSallocMemory(&reopt->reopttree->reoptnodes[id]->dualconscur) );
      reopt->reopttree->reoptnodes[id]->dualconscur->nvars = -1;

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconscur->vars, reopt->dualcons->vars, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconscur->vals, reopt->dualcons->vals, nbndchgs) );
      reopt->reopttree->reoptnodes[id]->dualconscur->nvars = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconscur->allocmem = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconscur->constype = reopttype == SCIP_REOPTTYPE_STRBRANCHED ? REOPT_CONSTYPE_STRBRANCHED : REOPT_CONSTYPE_INFSUBTREE;

      SCIPdebugMessage(" -> save dual information: node %lld, nvars %d, constype %d\n",
            SCIPnodeGetNumber(node), reopt->reopttree->reoptnodes[id]->dualconscur->nvars,
            reopt->reopttree->reoptnodes[id]->dualconscur->constype);
   }
   /* the constraint will be added after next */
   else
   {
      assert(reopt->reopttree->reoptnodes[id]->dualconsnex == NULL);
      SCIP_ALLOC( BMSallocMemory(&reopt->reopttree->reoptnodes[id]->dualconsnex) );
      reopt->reopttree->reoptnodes[id]->dualconsnex->nvars = -1;

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconsnex->vars, reopt->dualcons->vars, nbndchgs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconsnex->vals, reopt->dualcons->vals, nbndchgs) );
      reopt->reopttree->reoptnodes[id]->dualconsnex->nvars = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconsnex->allocmem = nbndchgs;
      reopt->reopttree->reoptnodes[id]->dualconsnex->constype = reopttype == SCIP_REOPTTYPE_STRBRANCHED ? REOPT_CONSTYPE_STRBRANCHED : REOPT_CONSTYPE_INFSUBTREE;

      SCIPdebugMessage(" -> save dual information: node %lld, nvars %d, constype %d\n",
            SCIPnodeGetNumber(node), reopt->reopttree->reoptnodes[id]->dualconsnex->nvars,
            reopt->reopttree->reoptnodes[id]->dualconsnex->constype);
   }


   return SCIP_OKAY;
}

/*
 * Add a pruned node the data structure.
 */
static
SCIP_RETCODE addNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< current node */
   SCIP_REOPTTYPE        reopttype,               /**< reason for storing the node*/
   SCIP_Bool             saveafterdual,           /**< save branching decisions after the first dual? */
   SCIP_Bool             isrootnode,              /**< node is the root node */
   SCIP_Bool             isfocusnode,             /**< node is the focus node */
   SCIP_LPSOLSTAT        lpsolstat,               /**< LP solution status */
   SCIP_Real             lowerbound               /**< lower bound of the node */
   )
{
   SCIP_NODE* parent;
   SCIP_Bool shrank;
   int nodeID;
   int parentID;

   assert(reopt != NULL);
   assert(node != NULL);

   parentID = -1;
   parent = NULL;
   shrank = FALSE;

   if( reopt->maxsavednodes == 0 )
      return SCIP_OKAY;

   assert(reopttype == SCIP_REOPTTYPE_TRANSIT
       || reopttype == SCIP_REOPTTYPE_INFSUBTREE
       || reopttype == SCIP_REOPTTYPE_STRBRANCHED
       || reopttype == SCIP_REOPTTYPE_LOGICORNODE
       || reopttype == SCIP_REOPTTYPE_LEAF
       || reopttype == SCIP_REOPTTYPE_PRUNED
       || reopttype == SCIP_REOPTTYPE_FEASIBLE);

   /** start clock */
   SCIPclockStart(reopt->savingtime, set);

   /* the node was created by reoptimization, i.e., we need to update the
    * stored data */
   if( SCIPnodeGetReoptID(node) >= 1 )
   {
      SCIP_Bool transintoorig;

      assert(reopttype != SCIP_REOPTTYPE_LEAF);
      assert(!isrootnode);

      nodeID = SCIPnodeGetReoptID(node);
      assert(nodeID < reopt->reopttree->allocmemnodes);
      assert(reopt->reopttree->reoptnodes[nodeID] != NULL);

      SCIPdebugMessage("update node %lld at ID %u:\n", SCIPnodeGetNumber(node), nodeID);

      transintoorig = FALSE;

      /* store in*/
      if( saveafterdual )
      {
         SCIP_CALL( saveAfterDualBranchings(reopt, blkmem, node, nodeID, &transintoorig) );
      }

      /* update constraint propagations */
      SCIP_CALL( updateConstraintPropagation(reopt, blkmem, node, nodeID, &transintoorig) );

      /* ensure that all variables are original */
      if( transintoorig )
      {
         SCIP_CALL( transformIntoOrig(reopt, nodeID) );
      }

      /* update the lowerbound if it was not set */
      if( SCIPsetIsLE(set, reopt->reopttree->reoptnodes[nodeID]->lowerbound, -1e+20) )
         reopt->reopttree->reoptnodes[nodeID]->lowerbound = lowerbound;

#ifdef SCIP_DEBUG
         int varnr;

         SCIPdebugMessage(" -> nvars: %d, ncons: %d, parentID: %d, reopttype: %d\n",
               reopt->reopttree->reoptnodes[nodeID]->nvars,
               reopt->reopttree->reoptnodes[nodeID]->nconss,
               reopt->reopttree->reoptnodes[nodeID]->parentID, reopttype);
         SCIPdebugMessage(" -> saved variables:\n");

         for (varnr = 0; varnr < reopt->reopttree->reoptnodes[nodeID]->nvars; varnr++)
            SCIPdebugMessage("  <%s> %s %f\n", SCIPvarGetName(reopt->reopttree->reoptnodes[nodeID]->vars[varnr]),
                  reopt->reopttree->reoptnodes[nodeID]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                  "=>" : "<=", reopt->reopttree->reoptnodes[nodeID]->varbounds[varnr]);

         for (varnr = 0; varnr < reopt->reopttree->reoptnodes[nodeID]->nafterdualvars; varnr++)
            SCIPdebugMessage("  <%s> %s %f (after dual red.)\n", SCIPvarGetName(reopt->reopttree->reoptnodes[nodeID]->afterdualvars[varnr]),
                  reopt->reopttree->reoptnodes[nodeID]->afterdualvarboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                  "=>" : "<=", reopt->reopttree->reoptnodes[nodeID]->afterdualvarbounds[varnr]);
#endif

      /** update LPI state if node is pseudobranched or feasible */
      switch( reopttype ) {
         case SCIP_REOPTTYPE_TRANSIT:
            assert(reopt->reopttree->reoptnodes[nodeID]->nconss == 0);

            if( reopt->shrinknodepath )
            {
               SCIP_CALL( shrinkNode(reopt, node, nodeID, &shrank, blkmem) );
            }

            goto TRANSIT;

            break;

         case SCIP_REOPTTYPE_LOGICORNODE:
         case SCIP_REOPTTYPE_LEAF:
            goto TRANSIT;
            break;

         case SCIP_REOPTTYPE_INFSUBTREE:
            /* delete the whole subtree induced be the current node */
            SCIP_CALL( deleteChildrenBelow(reopt->reopttree, blkmem, nodeID, FALSE, FALSE) );
            goto PSEUDO;
            break;

         case SCIP_REOPTTYPE_STRBRANCHED:
            /* dive through all children and change the reopttype to LEAF */
            SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, nodeID, SCIP_REOPTTYPE_PRUNED) );
            goto PSEUDO;
            break;

         case SCIP_REOPTTYPE_FEASIBLE:
            /* delete the subtree */
            if( reopt->reducetofrontier )
            {
               SCIP_CALL( deleteChildrenBelow(reopt->reopttree, blkmem, nodeID, FALSE, FALSE) );
               SCIPreoptResetDualcons(reopt, node, blkmem);
            }
            /* dive through all children and change the reopttype to PRUNED */
            else
            {
               SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, nodeID, SCIP_REOPTTYPE_PRUNED) );
            }
            goto FEASIBLE;
            break;

         case SCIP_REOPTTYPE_PRUNED:
            /* delete the subtree */
            if( reopt->reducetofrontier )
            {
               SCIP_CALL( deleteChildrenBelow(reopt->reopttree, blkmem, nodeID, FALSE, FALSE) );
               SCIPreoptResetDualcons(reopt, node, blkmem);
            }
            /* dive through all children and change the reopttype to LEAF */
            else
            {
               SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, nodeID, SCIP_REOPTTYPE_PRUNED) );
            }

            reopt->reopttree->nrediednodes++;
            reopt->reopttree->nrediednodesround++;

            goto PRUNED;

         default:
            break;
      }

      /** stop clock */
      SCIPclockStart(reopt->savingtime, set);

      return SCIP_OKAY;
   }

   /* get new IDs */
   SCIP_CALL( reopttreeCheckMemory(reopt->reopttree, blkmem) );

   /** the current node is the root node */
   if( isrootnode )
   {
      nodeID = 0;

      switch( reopttype ) {
         case SCIP_REOPTTYPE_TRANSIT:
            /* ensure that no dual constraints are stored */
            SCIPreoptResetDualcons(reopt, node, blkmem);

            /* update the lowerbound if it was not set */
            if( SCIPsetIsLE(set, reopt->reopttree->reoptnodes[nodeID]->lowerbound, -1e+20) )
               reopt->reopttree->reoptnodes[nodeID]->lowerbound = lowerbound;

            goto TRANSIT;
            break;

         case SCIP_REOPTTYPE_INFSUBTREE:
         case SCIP_REOPTTYPE_STRBRANCHED:
            reopt->reopttree->reoptnodes[0]->reopttype = reopttype;
            reopt->reopttree->reoptnodes[0]->dualfixing = TRUE;
            reopt->reopttree->reoptnodes[0]->nvars = 0;

            if( reopttype == SCIP_REOPTTYPE_INFSUBTREE )
            {
               /* delete the whole subtree induced be the current node */
               SCIP_CALL( deleteChildrenBelow(reopt->reopttree, blkmem, 0, FALSE, FALSE) );
            }

            SCIPdebugMessage("update node %d at ID %u:\n", 1, 0);
            SCIPdebugMessage(" -> nvars: 0, ncons: 0, parentID: -, reopttype: %d\n", reopttype);

            /* update the lowerbound if it was not set */
            if( SCIPsetIsLE(set, reopt->reopttree->reoptnodes[nodeID]->lowerbound, -1e+20) )
               reopt->reopttree->reoptnodes[nodeID]->lowerbound = lowerbound;

            goto PSEUDO;
            break;

         case SCIP_REOPTTYPE_FEASIBLE:
            reopt->reopttree->nfeasnodes++;
            reopt->reopttree->nfeasnodesround++;
            reopt->reopttree->reoptnodes[0]->reopttype = SCIP_REOPTTYPE_FEASIBLE;
            reopt->reopttree->reoptnodes[0]->dualfixing = FALSE;

            if( reopt->reopttree->reoptnodes[0]->childids != NULL && reopt->reopttree->reoptnodes[0]->nchilds > 0 )
            {
              /* delete the subtree */
               if( reopt->reducetofrontier )
               {
                  SCIP_CALL( deleteChildrenBelow(reopt->reopttree, blkmem, 0, FALSE, FALSE) );
                  SCIPreoptResetDualcons(reopt, node, blkmem);
               }
               /* dive through all children and change the reopttype to LEAF */
               else
               {
                  SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, 0, SCIP_REOPTTYPE_PRUNED) );
               }
            }
            else
               SCIPreoptResetDualcons(reopt, node, blkmem);

            /* update the lowerbound */
            reopt->reopttree->reoptnodes[nodeID]->lowerbound = lowerbound;

            SCIPdebugMessage("update node %d at ID %u:\n", 1, 0);
            SCIPdebugMessage(" -> nvars: 0, ncons: 0, parentID: -, reopttype: %d\n", reopttype);

            break;

         case SCIP_REOPTTYPE_PRUNED:
            reopt->reopttree->nprunednodes++;
            reopt->reopttree->nprunednodesround++;
            reopt->reopttree->reoptnodes[0]->reopttype = SCIP_REOPTTYPE_PRUNED;
            reopt->reopttree->reoptnodes[0]->dualfixing = FALSE;

            if( reopt->reopttree->reoptnodes[0]->childids != NULL && reopt->reopttree->reoptnodes[0]->nchilds > 0 )
            {
               /* delete the subtree */
               if( reopt->reducetofrontier )
               {
                  SCIP_CALL( deleteChildrenBelow(reopt->reopttree, blkmem, 0, FALSE, FALSE) );
                  SCIPreoptResetDualcons(reopt, node, blkmem);
               }
               /* dive through all children and change the reopttype to LEAF */
               else
               {
                  SCIP_CALL( changeReopttypeOfSubtree(reopt->reopttree, 0, SCIP_REOPTTYPE_PRUNED) );
               }
            }
            else
               SCIPreoptResetDualcons(reopt, node, blkmem);

            /* update the lowerbound if it was not set */
            if( SCIPsetIsLE(set, reopt->reopttree->reoptnodes[nodeID]->lowerbound, -1e+20) )
               reopt->reopttree->reoptnodes[nodeID]->lowerbound = lowerbound;

            SCIPdebugMessage("update node %d at ID %u:\n", 1, 0);
            SCIPdebugMessage(" -> nvars: 0, ncons: 0, parentID: -, reopttype: %d\n", reopttype);

            break;

         default:
            assert(reopttype == SCIP_REOPTTYPE_TRANSIT
                || reopttype == SCIP_REOPTTYPE_INFSUBTREE
                || reopttype == SCIP_REOPTTYPE_STRBRANCHED
                || reopttype == SCIP_REOPTTYPE_PRUNED
                || reopttype == SCIP_REOPTTYPE_FEASIBLE);
            break;
      }

      /* reset the information of dual bound changes */
      reopt->currentnode = -1;
      if( reopt->dualcons != NULL )
         reopt->dualcons->nvars = 0;

      /** stop clock */
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
      SCIP_CALL( getLastSavedNode(reopt, node, &parent, &parentID, &nbndchgdiff) );

      if( reopttype < SCIP_REOPTTYPE_INFSUBTREE
        && ((!reopt->dynamicdiffofnodes && nbndchgdiff <= reopt->maxdiffofnodes)
              || TRUE))
           /* ||(reopt->dynamicdiffofnodes && nbndchgdiff <= ceil(log10((SCIP_Real)(SCIPgetNOrigBinVars(scip) - MIN(SCIPgetNOrigBinVars(scip)-1,lengthBranchPath(reopt, parent))))/log10(2.0))) ) ) */
      {
         SCIPdebugMessage(" -> skip saving\n");

         /** stop clock */
         SCIPclockStop(reopt->savingtime, set);

         return SCIP_OKAY;
      }

      /** check if there are free slots to store the node */
      SCIP_CALL( reopttreeCheckMemory(reopt->reopttree, blkmem) );

      nodeID = (int) (size_t) SCIPqueueRemove(reopt->reopttree->openids);

      SCIPdebugMessage(" -> save at ID %d\n", nodeID);

      assert(reopt->reopttree->reoptnodes[nodeID] == NULL
         || (reopt->reopttree->reoptnodes[nodeID]->nvars == 0 && reopt->reopttree->reoptnodes[nodeID]->nconss == 0));
      assert(nodeID >= 1 && nodeID < reopt->reopttree->allocmemnodes);
      assert(!isrootnode);

      /** get memory for nodedata */
      assert(reopt->reopttree->reoptnodes[nodeID] == NULL || reopt->reopttree->reoptnodes[nodeID]->nvars == 0);
      SCIP_CALL( createReoptnode(reopt->reopttree, nodeID) );
      reopt->reopttree->reoptnodes[nodeID]->parentID = parentID;

      assert(parent != NULL );
      assert((SCIPnodeGetDepth(parent) == 0 && parentID == 0) || (SCIPnodeGetDepth(parent) >= 1 && parentID > 0));
      assert(nodeID >= 1);

      /** create the array of "child nodes" if they not exist */
      if( reopt->reopttree->reoptnodes[parentID]->childids == NULL || reopt->reopttree->reoptnodes[parentID]->allocchildmem == 0 )
      {
         SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[parentID], blkmem, 0, 10, 0) );
      }

      /** add the "child node" */
      SCIP_CALL( reoptAddChild(reopt->reopttree, parentID, nodeID, blkmem) );

      /* save branching path */
      SCIP_CALL( saveAncestorBranchings(reopt->reopttree, blkmem, node, parent, nodeID, parentID) );

      /* save bound changes after some dual reduction */
      if( saveafterdual )
      {
         SCIP_CALL( saveAfterDualBranchings(reopt, blkmem, node, nodeID, &transintoorig) );
      }
      else
      {
         SCIPdebugMessage(" -> skip saving bound changes after dual reductions.\n");
      }

      /** transform all bounds of branched variables and ensure that they are original. */
      SCIP_CALL( transformIntoOrig(reopt, nodeID) );

      /** save pseudo-constraints (if one exists) */
      if (SCIPnodeGetNAddedcons(node) >= 1)
      {
         assert(reopt->reopttree->reoptnodes[nodeID]->nconss == 0);

         SCIP_CALL( saveLocalConssData(reopt->reopttree, set, node, blkmem, nodeID) );
      }

      /* update the lowerbound if it was not set */
      if( SCIPsetIsLE(set, reopt->reopttree->reoptnodes[nodeID]->lowerbound, -1e+20) )
         reopt->reopttree->reoptnodes[nodeID]->lowerbound = lowerbound;

      /* set ID */
      SCIPnodeSetReoptID(node, nodeID);

      /* set the REOPTTYPE */
      SCIPnodeSetReopttype(node, reopttype);

#ifdef SCIP_DEBUG
      int varnr;
      SCIPdebugMessage("save node #%lld successful\n", SCIPnodeGetNumber(node));
      SCIPdebugMessage(" -> ID %d, nvars %d, ncons %d, reopttype %d\n",
            nodeID, reopt->reopttree->reoptnodes[nodeID]->nvars + reopt->reopttree->reoptnodes[nodeID]->nafterdualvars,
            reopt->reopttree->reoptnodes[nodeID]->nconss,
            reopttype);
      for (varnr = 0; varnr < reopt->reopttree->reoptnodes[nodeID]->nvars; varnr++)
      {
         SCIPdebugMessage("  <%s> %s %f\n", SCIPvarGetName(reopt->reopttree->reoptnodes[nodeID]->vars[varnr]),
               reopt->reopttree->reoptnodes[nodeID]->varboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                     "=>" : "<=", reopt->reopttree->reoptnodes[nodeID]->varbounds[varnr]);
      }
      for (varnr = 0; varnr < reopt->reopttree->reoptnodes[nodeID]->nafterdualvars; varnr++)
      {
         SCIPdebugMessage("  <%s> %s %f (after dual red.)\n", SCIPvarGetName(reopt->reopttree->reoptnodes[nodeID]->afterdualvars[varnr]),
               reopt->reopttree->reoptnodes[nodeID]->afterdualvarboundtypes[varnr] == SCIP_BOUNDTYPE_LOWER ?
                     "=>" : "<=", reopt->reopttree->reoptnodes[nodeID]->afterdualvarbounds[varnr]);
      }
#endif
   }

   switch( reopttype ) {
      case SCIP_REOPTTYPE_TRANSIT:
      case SCIP_REOPTTYPE_LOGICORNODE:
      case SCIP_REOPTTYPE_LEAF:
         TRANSIT:

         if( !shrank )
         {
            reopt->reopttree->reoptnodes[nodeID]->reopttype = reopttype;

            if( reopt->savelpbasis
             && reopttype != SCIP_REOPTTYPE_LOGICORNODE
             && isfocusnode && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
            {
               SCIP_CALL( saveLPIstate(reopt->reopttree, node, nodeID) );
            }
         }
         else
         {
            SCIPnodeSetReoptID(node, -1);
            SCIPnodeSetReopttype(node, SCIP_REOPTTYPE_NONE);
         }
         break;

      case SCIP_REOPTTYPE_INFSUBTREE:
      case SCIP_REOPTTYPE_STRBRANCHED:
         PSEUDO:

         assert(reopt->currentnode == SCIPnodeGetNumber(node));

         reopt->reopttree->reoptnodes[nodeID]->reopttype = reopttype;
         reopt->reopttree->reoptnodes[nodeID]->dualfixing = TRUE;

         /* save the basis if the node */
         if( reopt->savelpbasis
          && reopttype == SCIP_REOPTTYPE_STRBRANCHED
          && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
         {
            SCIP_CALL( saveLPIstate(reopt->reopttree, node, nodeID) );
         }

         /* get all the dual information and decide if the constraint need
          * to be added next or after next */
         SCIP_CALL( collectDualInformation(reopt, node, nodeID, reopttype, blkmem) );

         break;

      case SCIP_REOPTTYPE_FEASIBLE:
         FEASIBLE:
         reopt->reopttree->reoptnodes[nodeID]->reopttype = SCIP_REOPTTYPE_FEASIBLE;
         reopt->reopttree->reoptnodes[nodeID]->dualfixing = FALSE;
         reopt->reopttree->nfeasnodes++;
         reopt->reopttree->nfeasnodesround++;

         /**
          * save all information of the current feasible solution to separate this
          * solution in a following round (but only if all variablea are binary)
          * TODO: Verbesserungswuerdig
          */
         if( reopt->sepasolsloc && nodeID > 0 )
         {
            printf("TODO: implement storing a solution separating constraint.");
         }

         /* save the basis if the node */
         if( reopt->savelpbasis )
         {
            SCIP_CALL( saveLPIstate(reopt->reopttree, node, nodeID) );
         }

         break;

      case SCIP_REOPTTYPE_PRUNED:
         PRUNED:

         reopt->reopttree->reoptnodes[nodeID]->reopttype = SCIP_REOPTTYPE_PRUNED;
         reopt->reopttree->reoptnodes[nodeID]->dualfixing = FALSE;
         reopt->reopttree->nprunednodes++;
         reopt->reopttree->nprunednodesround++;

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
   }

   /* stop clock */
   SCIPclockStop(reopt->savingtime, set);

   /* reset the information of dual bound changes */
   reopt->currentnode = -1;
   if( reopt->dualcons != NULL )
      reopt->dualcons->nvars = 0;

   return SCIP_OKAY;
}

/* delete the stored information about dual bound changes of
 * the last focused node */
static
void deleteLastDualBndchgs(
   SCIP_REOPT*           reopt
   )
{
   assert(reopt != NULL);

   if( reopt->dualcons != NULL && reopt->dualcons->nvars > 0 )
   {
      SCIPdebugMessage("delete %d dual variable information about node %lld\n", reopt->dualcons->nvars, reopt->currentnode);
      reopt->dualcons->nvars = 0;
      reopt->currentnode = -1;
   }
}

/* build a global constrain to separate an infeasible subtree */
static
SCIP_RETCODE saveGlobalCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   REOPT_CONSTYPE        consttype                /**< reopttype of the constraint */
   )
{
   assert(reopt != NULL);
   assert(node != NULL);

   if( consttype == REOPT_CONSTYPE_INFSUBTREE )
   {
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
      SCIP_ALLOC( BMSallocMemory(&reopt->glbconss[nglbconss]) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vars, nvars) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vals, nvars) );

      reopt->glbconss[nglbconss]->allocmem = nvars;
      reopt->glbconss[nglbconss]->constype = REOPT_CONSTYPE_INFSUBTREE;

      assert(nvars <= reopt->consallocmem);

      SCIPnodeGetAncestorBranchings(node,
            reopt->glbconss[nglbconss]->vars,
            reopt->glbconss[nglbconss]->vals,
            reopt->consbounds,
            &nbranchvars,
            nvars);

      assert(nbranchvars <= reopt->consallocmem);

      if( nvars < nbranchvars )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vars, nvars, nbranchvars) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vals, nvars, nbranchvars) );
         nvars = nbranchvars;
         reopt->glbconss[nglbconss]->allocmem = nvars;

         SCIPnodeGetAncestorBranchings(node,
               reopt->glbconss[nglbconss]->vars,
               reopt->glbconss[nglbconss]->vals,
               reopt->consbounds,
               &nbranchvars,
               nvars);
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

         assert(SCIPsetIsFeasEQ(set, reopt->glbconss[nglbconss]->vals[v], 0) || SCIPsetIsFeasEQ(set, reopt->glbconss[nglbconss]->vals[v], 1));
      }

      /* increase the number of global constraints */
      reopt->nglbconss++;
   }

   return SCIP_OKAY;
}


/* move all id of child nodes from id1 to id2 */
static
SCIP_RETCODE reoptMoveIDs(
   SCIP_REOPTTREE*       reopttree,               /**< reopttree */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   id1,                     /**< source id */
   int                   id2                      /**< target id*/
   )
{
   int c;
   int nchilds_id1;
   int nchilds_id2;

   assert(reopttree != NULL);
   assert(blkmem != NULL);
   assert(0 <= id1 && id1 < reopttree->allocmemnodes);
   assert(0 <= id2 && id2 < reopttree->allocmemnodes);
   assert(reopttree->reoptnodes[id1] != NULL);
   assert(reopttree->reoptnodes[id2] != NULL);

   nchilds_id1 = reopttree->reoptnodes[id1]->nchilds;
   nchilds_id2 = reopttree->reoptnodes[id2]->nchilds;

   /* ensure that the array storing the child id's is large enough */
   SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id2], blkmem, 0, nchilds_id1+nchilds_id2, 0) );
   assert(reopttree->reoptnodes[id2]->allocchildmem >= nchilds_id1+nchilds_id2);

   SCIPdebugMessage("move %d IDs: %d -> %d\n", nchilds_id1, id1, id2);

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

/* apply all bound changes along the root path */
static
SCIP_RETCODE changeAncestorBranchings(
   SCIP_REOPTTREE*       reopttree,               /**< reopttree */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PROB*            transprob,               /**< transformed problem */
   SCIP_PROB*            origprob,                /**< original problem */
   SCIP_TREE*            tree,                    /**< search tree */
   SCIP_LP*              lp,                      /**< current LP */
   SCIP_BRANCHCAND*      branchcand,              /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,              /**< event queue */
   SCIP_NODE*            node_fix,                /**< node corresponding to the fixed part */
   SCIP_NODE*            node_cons,               /**< node corresponding to the pruned part */
   int                   id,                      /**< id of stored node */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   SCIP_REOPTNODE* reoptnode;
   int v;

   assert(reopttree != NULL);
   assert(node_fix != NULL || node_cons != NULL);
   assert(blkmem != NULL);
   assert(0 <= id && id < reopttree->allocmemnodes);
   assert(reopttree->reoptnodes[id] != NULL);

   reoptnode = reopttree->reoptnodes[id];

   /** copy memory to ensure that only original variables are saved */
   if( reoptnode->nvars == 0 && reoptnode->nafterdualvars == 0)
      return SCIP_OKAY;

   /* change the bounds along the branching path */
   SCIPdebugMessage(" -> change bound along the branching path:\n");

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

         if( node_fix != NULL )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(node_fix, blkmem, set, stat, transprob, origprob,
                  tree, lp, branchcand, eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
         }

         if( node_cons != NULL )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(node_cons, blkmem, set, stat, transprob, origprob,
                  tree, lp, branchcand, eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
      }
      else if(boundtype == SCIP_BOUNDTYPE_UPPER
           && SCIPsetIsLT(set, newbound, oldub)
           && SCIPsetIsFeasGE(set, newbound, oldlb))
      {
         SCIPvarAdjustUb(var, set, &newbound);

         if( node_fix != NULL )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(node_fix, blkmem, set, stat, transprob, origprob,
                  tree, lp, branchcand, eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }

         if( node_cons != NULL )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(node_cons, blkmem, set, stat, transprob, origprob,
                  tree, lp, branchcand, eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
      }

      SCIPdebugMessage("    <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);
   }

   /* fix bound affected by dual information at node_fix only */
   if( node_fix != NULL && reoptnode->nafterdualvars > 0 )
   {
      /* check the memory to convert this bound changes into 'normal' */
      SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], blkmem, reoptnode->nvars + reoptnode->nafterdualvars, 0, 0) );

      /* change the bounds */
      SCIPdebugMessage(" -> change bounds affected by dual information:\n");

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
            SCIP_CALL( SCIPnodeAddBoundchg(node_fix, blkmem, set, stat, transprob, origprob,
                  tree, lp, branchcand, eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );

            bndchgd = TRUE;
         }
         else if(boundtype == SCIP_BOUNDTYPE_UPPER
              && SCIPsetIsLT(set, newbound, oldub)
              && SCIPsetIsFeasGE(set, newbound, oldlb))
         {
            SCIPvarAdjustUb(var, set, &newbound);
            SCIP_CALL( SCIPnodeAddBoundchg(node_fix, blkmem, set, stat, transprob, origprob,
                  tree, lp, branchcand, eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );

            bndchgd = TRUE;
         }

         assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

         SCIPdebugMessage("    <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);

         if( bndchgd )
         {
            int nvars;

            nvars = reoptnode->nvars;
            reoptnode->vars[nvars] = reoptnode->afterdualvars[v];
            reoptnode->varbounds[nvars] = reoptnode->afterdualvarbounds[v];
            reoptnode->varboundtypes[nvars] = reoptnode->afterdualvarboundtypes[v];

            reoptnode->nvars++;
         }
      }

      /* free the afterdualvars, -bounds, and -boundtypes */
      BMSfreeBlockMemoryArray(blkmem, &reoptnode->afterdualvarboundtypes, reoptnode->allocafterdualvarmem);
      reoptnode->afterdualvarboundtypes = NULL;

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->afterdualvarbounds, reoptnode->allocafterdualvarmem);
      reoptnode->afterdualvarbounds = NULL;

      BMSfreeBlockMemoryArray(blkmem, &reoptnode->afterdualvars, reoptnode->allocafterdualvarmem);
      reoptnode->afterdualvars = NULL;

      reoptnode->nafterdualvars = 0;
      reoptnode->allocafterdualvarmem = 0;
   }

   return SCIP_OKAY;
}

/* add a constraint to ensure that at least one variable bound gets different */
static
SCIP_RETCODE addSplitcons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_PROB*            transprob,               /**< transformed problem */
   SCIP_PROB*            origprob,                /**< original problem */
   SCIP_TREE*            tree,                    /**< search tree */
   SCIP_LP*              lp,                      /**< current LP */
   SCIP_BRANCHCAND*      branchcand,              /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,              /**< event queue */
   SCIP_NODE*            node_cons,               /**< node corresponding to the pruned part */
   int                   id                       /**< id of stored node */
   )
{
   SCIP_CONS* cons;
   const char* name;
   int v;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
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
   assert(node_cons != NULL);

   if( reopt->reopttree->reoptnodes[id]->dualconscur->constype == REOPT_CONSTYPE_STRBRANCHED )
   {
      SCIPdebugMessage(" create a split-node #%lld\n", SCIPnodeGetNumber(node_cons));
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
      boundtype = SCIPsetIsFeasEQ(set, newbound, 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;

      assert(SCIPvarIsOriginal(var));
      SCIP_CALL( SCIPvarGetProbvarBound(&var, &newbound, &boundtype) );
      assert(SCIPvarIsTransformed(var));

      oldlb = SCIPvarGetLbLocal(var);
      oldub = SCIPvarGetUbLocal(var);

      /* negate the bound */
      newbound = 1 - newbound;
      boundtype = (SCIP_BOUNDTYPE) (1 - boundtype);

      if(boundtype == SCIP_BOUNDTYPE_LOWER
      && SCIPsetIsGT(set, newbound, oldlb)
      && SCIPsetIsFeasLE(set, newbound, oldub))
      {
         SCIPvarAdjustLb(var, set, &newbound);
         SCIP_CALL( SCIPnodeAddBoundchg(node_cons, blkmem, set, stat, transprob, origprob,
               tree, lp, branchcand, eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      }
      else if(boundtype == SCIP_BOUNDTYPE_UPPER
           && SCIPsetIsLT(set, newbound, oldub)
           && SCIPsetIsFeasGE(set, newbound, oldlb))
      {
         SCIPvarAdjustUb(var, set, &newbound);
         SCIP_CALL( SCIPnodeAddBoundchg(node_cons, blkmem, set, stat, transprob, origprob,
               tree, lp, branchcand, eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      }

      assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

      SCIPdebugMessage("  -> constraint consists of only one variable: <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? "=>" : "<=", newbound);
   }
   else
   {
      assert(reopt->consallocmem >= reopt->reopttree->reoptnodes[id]->dualconscur->nvars);

      for(v = 0; v < reopt->reopttree->reoptnodes[id]->dualconscur->nvars; v++)
      {
         reopt->consvars[v] = reopt->reopttree->reoptnodes[id]->dualconscur->vars[v];
         reopt->consvals[v] = reopt->reopttree->reoptnodes[id]->dualconscur->vals[v];
         reopt->consbounds[v] = SCIPsetIsFeasEQ(set, reopt->consvals[v], 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;

        assert(SCIPvarIsOriginal(reopt->consvars[v]));
        SCIP_CALL( SCIPvarGetProbvarBound(&reopt->consvars[v], &reopt->consvals[v], &reopt->consbounds[v]) );
        assert(SCIPvarIsTransformed(reopt->consvars[v]));
        assert(SCIPvarGetStatus(reopt->consvars[v]) != SCIP_VARSTATUS_MULTAGGR);

        if ( SCIPsetIsFeasEQ(set, reopt->consvals[v], 1) )
        {
           SCIP_CALL( SCIPvarNegate(reopt->consvars[v], blkmem, set, stat, &reopt->consvars[v]) );
           assert(SCIPvarIsNegated(reopt->consvars[v]));
        }
      }

      if( reopt->reopttree->reoptnodes[id]->dualconscur->constype == REOPT_CONSTYPE_INFSUBTREE )
         name = "infsubtree";
      else
      {
         assert(reopt->reopttree->reoptnodes[id]->dualconscur->constype == REOPT_CONSTYPE_STRBRANCHED);
         name = "splitcons";
      }

      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name,
            reopt->reopttree->reoptnodes[id]->dualconscur->nvars, reopt->consvars,
            FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

      SCIPdebugMessage(" -> added constraint in node #%lld\n", SCIPnodeGetNumber(node_cons));
      SCIPdebugPrintCons(scip, cons, NULL);

      SCIP_CALL( SCIPaddConsNode(scip, node_cons, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}

/* fix all bounds stored in dualconscur in the given @param node_fix */
static
SCIP_RETCODE fixBounds(
   SCIP_REOPTTREE*       reopttree,               /**< reopttree */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PROB*            transprob,               /**< transformed problem */
   SCIP_PROB*            origprob,                /**< original problem */
   SCIP_TREE*            tree,                    /**< search tree */
   SCIP_LP*              lp,                      /**< current LP */
   SCIP_BRANCHCAND*      branchcand,              /**< branching candidates */
   SCIP_EVENTQUEUE*      eventqueue,              /**< event queue */
   SCIP_NODE*            node_fix,                /**< node corresponding to the fixed part */
   int                   id,                      /**< id of stored node */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   int v;

   assert(reopttree != NULL);
   assert(node_fix != NULL);
   assert(blkmem != NULL);
   assert(0 < id && id < reopttree->allocmemnodes);
   assert(reopttree->reoptnodes[id] != NULL);
   assert(reopttree->reoptnodes[id]->dualfixing);
   assert(reopttree->reoptnodes[id]->dualconscur != NULL);

   /* ensure that the arrays to store the bound changes are large enough */
   SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[id], blkmem, reopttree->reoptnodes[id]->nvars + reopttree->reoptnodes[id]->dualconscur->nvars, 0, 0) );

   SCIPdebugMessage(" -> reconstruct bound changes at node #%lld (save at ID %d):\n", SCIPnodeGetNumber(node_fix), id);

   for(v = 0; v < reopttree->reoptnodes[id]->dualconscur->nvars; v++)
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_BOUNDTYPE boundtype;
      SCIP_Bool bndchgd;

      var = reopttree->reoptnodes[id]->dualconscur->vars[v];
      val = reopttree->reoptnodes[id]->dualconscur->vals[v];
      boundtype = SCIPsetIsFeasEQ(set, val, 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;

      SCIP_CALL(SCIPvarGetProbvarBound(&var, &val, &boundtype));
      assert(SCIPvarIsTransformedOrigvar(var));

      bndchgd = FALSE;

      if(boundtype == SCIP_BOUNDTYPE_LOWER
      && SCIPsetIsGT(set, val, SCIPvarGetLbLocal(var))
      && SCIPsetIsFeasLE(set, val, SCIPvarGetUbLocal(var)))
      {
         SCIPvarAdjustLb(var, set, &val);
         SCIP_CALL( SCIPnodeAddBoundchg(node_fix, blkmem, set, stat, transprob, origprob,
               tree, lp, branchcand, eventqueue, var, val, SCIP_BOUNDTYPE_LOWER, FALSE) );

         bndchgd = TRUE;
      }
      else if(boundtype == SCIP_BOUNDTYPE_UPPER
           && SCIPsetIsLT(set, val, SCIPvarGetUbLocal(var))
           && SCIPsetIsFeasGE(set, val, SCIPvarGetLbLocal(var)))
      {
         SCIPvarAdjustUb(var, set, &val);
         SCIP_CALL( SCIPnodeAddBoundchg(node_fix, blkmem, set, stat, transprob, origprob,
               tree, lp, branchcand, eventqueue, var, val, SCIP_BOUNDTYPE_UPPER, FALSE) );

         bndchgd = TRUE;
      }
      else if(boundtype != SCIP_BOUNDTYPE_LOWER && boundtype != SCIP_BOUNDTYPE_UPPER)
      {
         printf("** Unknown boundtype: %d **\n", boundtype);
         assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);
      }

      SCIPdebugMessage("  <%s> %s %g\n", SCIPvarGetName(var), boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", val);

      /** add variable and bound to branching path information, because we don't want to delete this data */
      if( bndchgd )
      {
         int pos;
         SCIP_Real constant;
         SCIP_Real scalar;

         pos = reopttree->reoptnodes[id]->nvars;

         reopttree->reoptnodes[id]->vars[pos] = var;
         SCIP_CALL( SCIPvarGetOrigvarSum(&reopttree->reoptnodes[id]->vars[pos], &scalar, &constant) );
         assert(SCIPvarIsOriginal(reopttree->reoptnodes[id]->vars[pos]));

         reopttree->reoptnodes[id]->varbounds[pos] = reopttree->reoptnodes[id]->dualconscur->vals[v];
         reopttree->reoptnodes[id]->varboundtypes[pos] = (SCIPsetIsFeasEQ(set, reopttree->reoptnodes[id]->varbounds[pos], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
         reopttree->reoptnodes[id]->nvars++;
      }
   }

   /* delete dualconscur and move dualconsnex -> dualconscur */
   BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->dualconscur->vals, reopttree->reoptnodes[id]->dualconscur->allocmem);
   BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->dualconscur->vars, reopttree->reoptnodes[id]->dualconscur->allocmem);
   SCIPfreeMemory(scip, &reopttree->reoptnodes[id]->dualconscur);
   reopttree->reoptnodes[id]->dualconscur = NULL;

   if( reopttree->reoptnodes[id]->dualconsnex != NULL )
   {
      reopttree->reoptnodes[id]->dualconscur = reopttree->reoptnodes[id]->dualconsnex;
      reopttree->reoptnodes[id]->dualconsnex = NULL;
   }

   return SCIP_OKAY;
}

/* add all local constraints stored at ID id */
static
SCIP_RETCODE addLocalConss(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node_fix,                /**< node correspinding to the fixed part */
   SCIP_NODE*            node_cons,               /**< node corresponding to the pruned part */
   int                   id                       /**< id of stored node */
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
   assert(node_fix != NULL || node_cons != NULL);
   assert(0 < id && id < reopt->reopttree->allocmemnodes);

   if( reopt->reopttree->reoptnodes[id]->nconss == 0 )
      return SCIP_OKAY;

   for( c = 0; c < reopt->reopttree->reoptnodes[id]->nconss; c++ )
   {
      SCIP_CONS* cons_node_fix;
      SCIP_CONS* cons_node_cons;
      LOGICORDATA* consdata;
      int v;

      consdata = reopt->reopttree->reoptnodes[id]->conss[c];
      assert(consdata != NULL);
      assert(consdata->nvars > 0);
      assert(consdata->allocmem >= consdata->nvars);
      assert(reopt->consallocmem >= consdata->nvars);

      /* iterate over all variable and transform them */
      for(v = 0; v < consdata->nvars; v++)
      {
         reopt->consvars[v] = consdata->vars[v];
         reopt->consvals[v] = consdata->vals[v];
         reopt->consbounds[v] = SCIPsetIsFeasEQ(set, reopt->consvals[v], 0) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER;

         assert(SCIPvarIsOriginal(reopt->consvars[v]));
         SCIP_CALL( SCIPvarGetProbvarBound(&reopt->consvars[v], &reopt->consvals[v], &reopt->consbounds[v]) );
         assert(SCIPvarIsTransformed(reopt->consvars[v]));
         assert(SCIPvarGetStatus(reopt->consvars[v]) != SCIP_VARSTATUS_MULTAGGR);

         if( SCIPsetIsFeasEQ(set, reopt->consvals[v], 1) )
         {
            SCIP_CALL( SCIPvarNegate(reopt->consvars[v], blkmem, set, stat, &reopt->consvars[v]) );
            assert(SCIPvarIsNegated(reopt->consvars[v]));
         }
      }

      assert(consdata->constype == REOPT_CONSTYPE_INFSUBTREE || consdata->constype == REOPT_CONSTYPE_STRBRANCHED);

      if( consdata->constype == REOPT_CONSTYPE_INFSUBTREE )
         name = "infsubtree";
      else
         name = "splitcons";


      /* create the constraints and add them to the corresponding nodes */
      if( node_fix != NULL )
      {
         SCIP_CALL( SCIPcreateConsLogicor(scip, &cons_node_fix, name, consdata->nvars, reopt->consvars,
               FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

         SCIP_CALL( SCIPaddConsNode(scip, node_fix, cons_node_fix, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons_node_fix) );
      }

      if( node_cons != NULL )
      {
         SCIP_CALL( SCIPcreateConsLogicor(scip, &cons_node_cons, name, consdata->nvars, reopt->consvars,
               FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );

         SCIP_CALL( SCIPaddConsNode(scip, node_cons, cons_node_cons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons_node_cons) );
      }
   }

   if( node_fix != NULL && node_cons != NULL )
   {
      SCIPdebugMessage(" -> added %d constraint(s) at node #%lld and #%lld\n", c, SCIPnodeGetNumber(node_fix), SCIPnodeGetNumber(node_cons));
   }
   else
   {
      assert(node_fix != NULL || node_cons != NULL);
      SCIPdebugMessage(" -> added %d constraint(s) at node #%lld\n", c, node_fix != NULL ? SCIPnodeGetNumber(node_fix) : SCIPnodeGetNumber(node_cons));
   }

   return SCIP_OKAY;
}

static
void resetStats(
   SCIP_REOPT*           reopt
   )
{
   assert(reopt != NULL);

   reopt->lastbranched = -1;
   reopt->currentnode = -1;
   reopt->lastseennode = -1;
   reopt->reopttree->nbranchednodesround = 0;
   reopt->reopttree->nfeasnodesround = 0;
   reopt->reopttree->ninfeasnodesround= 0;
   reopt->reopttree->nprunednodesround = 0;
   reopt->reopttree->nrediednodesround = 0;

   return;
}

/*
 * check if the node is infeasible or redundant due to strong branching
 */
static
SCIP_RETCODE dryBranch(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_Bool*            runagain,                /**< pointer to store of this method should run again */
   int                   id                       /**< id of stored node */
   )
{
   SCIP_REOPTNODE* reoptnode;
   int* cutoffchilds;
   int ncutoffchilds;
   int* redchilds;
   int nredchilds;
   int c;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(0 <= id && id < reopt->reopttree->allocmemnodes);
   assert(reopt->reopttree->reoptnodes != NULL);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   reoptnode = reopt->reopttree->reoptnodes[id];

   *runagain = FALSE;
   ncutoffchilds = 0;
   nredchilds = 0;

   SCIPdebugMessage("start dry branching of node at ID %d:\n", id);

   /* allocate buffer arrays */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cutoffchilds, reoptnode->nchilds) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &redchilds, reoptnode->nchilds) );

   /* iterate over all child nodes and check each bound changes
    * for redundancy and conflict */
   for(c = 0; c < reoptnode->nchilds; c++)
   {
      SCIP_REOPTNODE* child;
      SCIP_Bool cutoff;
      SCIP_Bool redundant;
      int* redundantvars;
      int nredundantvars;
      int v;
      int childid;

      cutoff = FALSE;
      redundant = FALSE;
      nredundantvars = 0;

      childid = reoptnode->childids[c];
      child = reopt->reopttree->reoptnodes[childid];
      assert(child != NULL);

      SCIPdebugMessage("-> check child at ID %d (%d vars, %d conss):\n", childid, child->nvars, child->nconss);

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
               nredundantvars++;
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
               child->nvars--;
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

      /* the node is redundant because all bound changes were redundant */
      if( child->nvars > 0 && child->nvars == nredundantvars )
      {
         redundant = TRUE;
         SCIPdebugMessage(" -> redundant node found.\n");
      }

      if( cutoff )
      {
         cutoffchilds[ncutoffchilds] = childid;
         ncutoffchilds++;
      }
      else if( redundant )
      {
         redchilds[nredchilds] = childid;
         nredchilds++;
      }
   }

   SCIPdebugMessage("-> found %d redundant and %d infeasible nodes\n", nredchilds, ncutoffchilds);

   c = 0;

   /* delete all nodes that can be cutoff */
   while( ncutoffchilds > 0 )
   {
      /* delete the node and the induced subtree */
      SCIP_CALL( deleteChildrenBelow(reopt->reopttree, blkmem, cutoffchilds[ncutoffchilds-1], TRUE, TRUE) );

      /* find the position in the childid array */
      c = 0;
      while( reoptnode->childids[c] != cutoffchilds[ncutoffchilds-1] && c < reoptnode->nchilds )
         c++;
      assert(reoptnode->childids[c] == cutoffchilds[ncutoffchilds-1]);

      /* replace the ID at position c by the last ID */
      reoptnode->childids[c] = reoptnode->childids[reoptnode->nchilds-1];
      reoptnode->nchilds--;

      /* decrease the number of nodes to cutoff */
      ncutoffchilds--;
   }

   c = 0;

   /* replace all redundant nodes their child nodes or cutoff the node if it is a leaf */
   while( nredchilds > 0 )
   {
      /* find the position in the childid array */
      c = 0;
      while( reoptnode->childids[c] != redchilds[nredchilds-1] && c < reoptnode->nchilds )
         c++;
      assert(reoptnode->childids[c] == redchilds[nredchilds-1]);

      /* the node is a leaf and we can cutoff them  */
      if( reopt->reopttree->reoptnodes[redchilds[nredchilds-1]]->nchilds == 0 )
      {
         /* delete the node and the induced subtree */
         SCIP_CALL( deleteChildrenBelow(reopt->reopttree, blkmem, redchilds[nredchilds-1], TRUE, TRUE) );

         /* replace the ID at position c by the last ID */
         reoptnode->childids[c] = reoptnode->childids[reoptnode->nchilds-1];
         reoptnode->nchilds--;

         /* decrease the number of redundant nodes */
         nredchilds--;
      }
      else
      {
         int cc;
         int ncc;

         /* replace the ID at position c by the last ID */
         reoptnode->childids[c] = reoptnode->childids[reoptnode->nchilds-1];
         reoptnode->nchilds--;

         ncc = reopt->reopttree->reoptnodes[redchilds[nredchilds-1]]->nchilds;

         /* check the memory */
         SCIP_CALL( reoptnodeCheckMemory(reopt->reopttree->reoptnodes[id], blkmem, 0, reoptnode->nchilds+ncc, 0) );

         /* add all IDs of child nodes to the current node */
         for(cc = 0; cc < ncc; cc++)
         {
            reoptnode->childids[reoptnode->nchilds] = reopt->reopttree->reoptnodes[redchilds[nredchilds-1]]->childids[cc];
            reoptnode->nchilds++;
         }

         /* delete the redundant node */
         SCIP_CALL( reopttreeDeleteNode(reopt->reopttree, blkmem, redchilds[nredchilds-1], TRUE) );

         /* decrease the number of redundant nodes */
         nredchilds--;

         /* update the flag to rerun this method */
         *runagain = TRUE;
      }
   }

   /* free buffer arrays */
   SCIPsetFreeBufferArray(set, &cutoffchilds);
   SCIPsetFreeBufferArray(set, &redchilds);

   return SCIP_OKAY;
}

/* returns the number of leave nodes of the induced subtree */
static
int reoptGetNLeaves(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   id                       /**< id of stored node */
   )
{
   int i;
   int nleaves;

   assert(reopt != NULL);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   nleaves = 0;

   for(i = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++)
   {
      int childid;

      childid = reopt->reopttree->reoptnodes[id]->childids[i];

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
         nleaves++;
      else
         nleaves += reoptGetNLeaves(reopt, childid);
   }

   return nleaves;
}

/* returns all leave if of the subtree induced by @param node */
static
SCIP_RETCODE reoptGetLeaves(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure*/
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   id,                      /**< id of stored node */
   int*                  leaves,                  /**< array of leave nodes */
   int                   mem,                     /**< allocated memory */
   int*                  nleaves                  /**< pointer to store the number of leave nodes */
   )
{
   int i;
   int l;

   assert(reopt != NULL);
   assert(mem > 0 && leaves != NULL);
   assert((*nleaves) >= 0);
   assert(id >= 0);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   for(i = 0, l = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++)
   {
      int childid;

      assert(*nleaves + 1 <= mem);

      childid = reopt->reopttree->reoptnodes[id]->childids[i];

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
      {
         leaves[l] = reopt->reopttree->reoptnodes[id]->childids[i];
         l++;
         (*nleaves)++;
      }
      else
      {
         int nleaves2;

         nleaves2 = 0;
         SCIP_CALL( reoptGetLeaves(reopt, set, blkmem, childid, &leaves[l], mem, &nleaves2) );
         l += nleaves2;
         (*nleaves) += nleaves2;
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE reoptResetTree(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(blkmem != NULL);

   /* clear the tree */
   SCIP_CALL( clearReoptnodes(reopt->reopttree, blkmem, FALSE) );

   /* reset the dual constraint */
   if( reopt->dualcons != NULL )
      reopt->dualcons->nvars = 0;

   reopt->currentnode = -1;

   return SCIP_OKAY;
}

static
SCIP_RETCODE reoptRestart(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);

   /* clear the tree */
   SCIP_CALL( reoptResetTree(reopt, blkmem) );

   reopt->nrestarts += 1;

   return SCIP_OKAY;
}

static
int reopttreeGetNNodes(
   SCIP_REOPTTREE*       reopttree,               /**< reopttree */
   int                   id                       /**< id of stored node */
   )
{
   int nnodes;
   int i;

   assert(reopttree != NULL);
   assert(0 <= id && id < reopttree->allocmemnodes);

   nnodes = 0;

   for(i = 0; i < reopttree->reoptnodes[id]->nchilds; i++)
   {
      nnodes += reopttreeGetNNodes(reopttree, reopttree->reoptnodes[id]->childids[i]);
   }

   return nnodes + 1;
}

/*
 * public methods
 */

/** creates reopt data */
SCIP_RETCODE SCIPreoptCreate(
   SCIP_REOPT**          reopt,                   /**< pointer to reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   int s;

   assert(reopt != NULL);

   SCIP_ALLOC( BMSallocMemory(reopt) );
   (*reopt)->runsize = DEFAULT_MEM_RUN;
   (*reopt)->run = 0;
   (*reopt)->nobjvars = 0;
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
   (*reopt)->nrestarts = 0;
   (*reopt)->consallocmem = 0;
   (*reopt)->consvars = NULL;
   (*reopt)->consvals = NULL;
   (*reopt)->consbounds = NULL;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*reopt)->objs, (*reopt)->runsize) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*reopt)->lastbestsol, (*reopt)->runsize) );

   for(s = 0; s < (*reopt)->runsize; s++)
   {
      (*reopt)->objs[s] = NULL;
      (*reopt)->lastbestsol[s] = NULL;
   }

   /* clocks */
   SCIP_CALL( SCIPclockCreate(&(*reopt)->savingtime, SCIP_CLOCKTYPE_DEFAULT) );

   /* get parameters */
   SCIP_CALL( SCIPsetGetBoolParam(set, "reoptimization/globalcons/sepainfsubtrees", &(*reopt)->sepasubtreesglb) );
   SCIP_CALL( SCIPsetGetBoolParam(set, "reoptimization/globalcons/sepasols", &(*reopt)->sepasolsglb) );

   SCIP_CALL( SCIPsetGetBoolParam(set, "reoptimization/localcons/sepasols", &(*reopt)->sepasolsloc) );

   SCIP_CALL( SCIPsetGetBoolParam(set, "reoptimization/reducetofrontier", &(*reopt)->reducetofrontier) );
   SCIP_CALL( SCIPsetGetBoolParam(set, "reoptimization/savelpbasis", &(*reopt)->savelpbasis) );
   SCIP_CALL( SCIPsetGetBoolParam(set, "reoptimization/shrinktransit", &(*reopt)->shrinknodepath) );
   SCIP_CALL( SCIPsetGetBoolParam(set, "reoptimization/dynamicdiffofnodes", &(*reopt)->dynamicdiffofnodes));

   SCIP_CALL( SCIPsetGetRealParam(set, "reoptimization/delay", &(*reopt)->localdelay) );
   SCIP_CALL( SCIPsetGetRealParam(set, "reoptimization/objsimrootLP", &(*reopt)->objsimrootlp) );

   SCIP_CALL( SCIPsetGetIntParam(set, "reoptimization/maxsavednodes", &(*reopt)->maxsavednodes) );
   SCIP_CALL( SCIPsetGetIntParam(set, "reoptimization/maxdiffofnodes", &(*reopt)->maxdiffofnodes) );
   SCIP_CALL( SCIPsetGetIntParam(set, "reoptimization/solvelp", &(*reopt)->solvelp) );
   SCIP_CALL( SCIPsetGetIntParam(set, "reoptimization/solvelpdiff", &(*reopt)->solvelpdiff) );
   SCIP_CALL( SCIPsetGetIntParam(set, "reoptimization/forceheurrestart", &(*reopt)->forceheurrestart) );

   /* create and initialize SCIP_SOLTREE */
   SCIP_ALLOC( BMSallocMemory(&(*reopt)->soltree) );
   SCIP_CALL( createSolTree((*reopt)->soltree, blkmem) );

   /* create and initialize SCIP_REOPTTREE */
   SCIP_ALLOC( BMSallocMemory(&(*reopt)->reopttree) );
   SCIP_CALL( createReopttree((*reopt)->reopttree, blkmem) );

   return SCIP_OKAY;
}

/** frees reopt data */
SCIP_RETCODE SCIPreoptFree(
   SCIP_REOPT**          reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_PRIMAL*          origprimal,              /**< original primal */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   int p;

   assert(reopt != NULL);
   assert(*reopt != NULL);
   assert(set != NULL);
   assert(origprimal != NULL);
   assert(blkmem != NULL);

   /* free reopttree */
   SCIP_CALL( freeReoptTree((*reopt)->reopttree, blkmem) );

   /* free solutions */
   for( p = (*reopt)->run-1; p >= 0; p-- )
   {
      if( (*reopt)->soltree->sols[p] != NULL )
      {
         BMSfreeMemoryArray(&(*reopt)->soltree->sols[p]);
         (*reopt)->soltree->sols[p] = NULL;
      }

      if( (*reopt)->objs[p] != NULL )
      {
         BMSfreeMemoryArray(&(*reopt)->objs[p]);
      }
   }

   /* free solution tree */
   SCIP_CALL( freeSolTree((*reopt), set, origprimal, blkmem) );

   if( (*reopt)->dualcons != NULL )
   {
      if( (*reopt)->dualcons->allocmem > 0 )
      {
         BMSfreeBlockMemoryArray(blkmem, &(*reopt)->dualcons->vals, (*reopt)->dualcons->allocmem);
         BMSfreeBlockMemoryArray(blkmem, &(*reopt)->dualcons->vars, (*reopt)->dualcons->allocmem);
         BMSfreeMemoryArray(&(*reopt)->dualcons);
         (*reopt)->dualcons = NULL;
      }
   }

   if( (*reopt)->glbconss != NULL && (*reopt)->allocmemglbconss > 0 )
   {
      (*reopt)->nglbconss--;

      /* free all constraint */
      while( (*reopt)->nglbconss > 0 )
      {
         int c;
         c = (*reopt)->nglbconss;

         if( (*reopt)->glbconss[c] != NULL )
         {
            if( (*reopt)->glbconss[c]->allocmem > 0 )
            {
               BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss[c]->vals, (*reopt)->glbconss[c]->allocmem);
               BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss[c]->vars, (*reopt)->glbconss[c]->allocmem);
               (*reopt)->glbconss[c]->allocmem = 0;
            }
            BMSfreeMemory(&(*reopt)->glbconss[c]);
         }

         (*reopt)->nglbconss--;
      }
      assert((*reopt)->nglbconss == 0);

      BMSfreeBlockMemoryArray(blkmem, &(*reopt)->glbconss, (*reopt)->allocmemglbconss);
      (*reopt)->allocmemglbconss = 0;
   }


   /* clocks */
   SCIPclockFree(&(*reopt)->savingtime);

   if( (*reopt)->consallocmem > 0 )
   {
      BMSfreeBlockMemoryArray(blkmem, &(*reopt)->consbounds, (*reopt)->consallocmem);
      BMSfreeBlockMemoryArray(blkmem, &(*reopt)->consvals, (*reopt)->consallocmem);
      BMSfreeBlockMemoryArray(blkmem, &(*reopt)->consvars, (*reopt)->consallocmem);
   }

   BMSfreeBlockMemoryArray(blkmem, &(*reopt)->lastbestsol, (*reopt)->runsize);
   BMSfreeBlockMemoryArray(blkmem, &(*reopt)->objs, (*reopt)->runsize);
   BMSfreeMemory(reopt);

   return SCIP_OKAY;
}

/* returns the number constraints added by reoptimization plug-in */
int SCIPreoptGetNAddedConss(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node                     /**< node of the search tree */
)
{
   int id;

   assert(reopt != NULL);
   assert(node != NULL);

   id = SCIPnodeGetReoptID(node);

   if( id >= 1 && reopt->reopttree->reoptnodes[id]->nconss > 0 )
      return MAX(SCIPnodeGetNAddedcons(node), reopt->reopttree->reoptnodes[id]->nconss);
   else
      return SCIPnodeGetNAddedcons(node);
}

/* returns the number of stored nodes providing an improving feasible LP solution */
int SCIPreoptGetNFeasNodes(
   SCIP_REOPT*           reopt                    /* reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->nfeasnodes;
}

/* returns the number of stored nodes that exceeded the cutoff bound */
int SCIPreoptGetNPrunedNodes(
   SCIP_REOPT*           reopt                    /* reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->nprunednodes;
}

/* rerturns the number of reoptimized nodes that were cutoff in the same iteration */
int SCIPreoptGetNRediedNodes(
   SCIP_REOPT*           reopt                    /* reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->nrediednodes;
}

/* returns the number of stored nodes with an infeasible LP */
int SCIPreoptGetNInfNodes(
   SCIP_REOPT*           reopt                    /* reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->reopttree->ninfeasnodes;
}

/** add a solution to the last run */
SCIP_RETCODE SCIPreoptAddSol(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,              /**< original primal */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_SOL*             sol,                     /**< solution to add */
   SCIP_Bool             bestsol,                 /**< is the current solution an optimal solution? */
   SCIP_Bool*            added,                   /**< pointer to store the information if the soltion was added */
   SCIP_VAR**            vars,                    /**< variable array */
   int                   nvars,                   /**< number of variables */
   int                   run                      /**< number of the current run (1,2,...) */
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
   if( heur != NULL && strcmp(SCIPheurGetName(heur), "reoptsols") == 0 )
   {
      *added = FALSE;

      if( bestsol )
      {
         SCIP_SOL* copysol;

         reopt->noptsolsbyreoptsol++;

         /* we add the solution to the best solutions found so far */
         SCIP_CALL( SCIPsolCopy(&copysol, blkmem, set, stat, origprimal, sol) );
         reopt->lastbestsol[reopt->run-1] = copysol;
      }

      return SCIP_OKAY;
   }

   if( bestsol )
      reopt->noptsolsbyreoptsol = 0;

   /* check memory */
   SCIP_CALL( ensureSolsSize(reopt, set, blkmem, reopt->soltree->nsols[run-1], run-1) );

   solnode = NULL;

   /** add solution to solution tree */
   SCIP_CALL( soltreeAddSol(reopt, set, stat, origprimal, blkmem, vars, sol, &solnode, nvars, bestsol, added) );

   if( (*added) )
   {
      assert(solnode != NULL);

      /** add solution */
      insertpos = reopt->soltree->nsols[run-1];
      reopt->soltree->sols[run-1][insertpos] = solnode;
      reopt->soltree->nsols[run-1]++;
      assert(reopt->soltree->nsols[run-1] <= set->reopt_savesols);
   }

   return SCIP_OKAY;
}

/* add optimal solution */
SCIP_RETCODE SCIPreoptAddOptSol(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SOL*             sol,                     /**< solution to add */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal               /**< original primal */
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
   reopt->lastbestsol[reopt->run-1] = solcopy;
   printf(">> add best sol at pos %d\n", reopt->run-1);

   return SCIP_OKAY;
}

/* add a run */
SCIP_RETCODE SCIPreoptAddRun(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_VAR**            origvars,                /**< original variables */
   int                   norigvars,               /**< number of original variables */
   int                   size                     /**< number of expected solutions */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem !=  NULL);
   assert(origvars != NULL);

   /* increase number of runs */
   reopt->run++;

   /* check memory */
   SCIP_CALL( ensureRunSize(reopt, set, reopt->run, blkmem) );

   /* allocate memory */
   reopt->soltree->solssize[reopt->run-1] = size;
   SCIP_ALLOC( BMSallocMemoryArray(&reopt->soltree->sols[reopt->run-1], size) );

   if( reopt->nobjvars == 0 )
   {
      assert(reopt->run == 1);
      reopt->nobjvars = norigvars;

      /* allocate permanent memory */
      if( norigvars > 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->consvars, norigvars) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->consvals, norigvars) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->consbounds, norigvars) );
         reopt->consallocmem = norigvars;
      }
   }
   assert(reopt->nobjvars == norigvars);

   /* save the objective function */
   SCIP_CALL( SCIPreoptSaveNewObj(reopt, set, blkmem, origvars, norigvars) );

   resetStats(reopt);

   return SCIP_OKAY;
}

/* get the number of checked during the reoptimization process */
int SCIPreoptGetNCheckedsols(
   SCIP_REOPT*           reopt                    /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->ncheckedsols;
}

/* update the number of checked during the reoptimization process */
void SCIPreoptSetNCheckedsols(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   ncheckedsols             /**< number of updated solutions */
   )
{
   assert(reopt != NULL);

   reopt->ncheckedsols += ncheckedsols;
}

/* get the number of checked during the reoptimization process */
int SCIPreoptGetNImprovingsols(
   SCIP_REOPT*           reopt                    /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->nimprovingsols;
}

/* update the number of checked during the reoptimization process */
void SCIPreoptSetNImprovingsols(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   nimprovingsols           /**< number of improving solutions */
   )
{
   assert(reopt != NULL);

   reopt->nimprovingsols += nimprovingsols;
}

/* returns number of solution of a given run */
int SCIPreoptGetNSolsRun(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   int                   run                      /**< number of the run (1,2,..) */
   )
{
   assert(reopt != NULL);
   assert(0 < run && run <= reopt->runsize);

   if( reopt->soltree->sols[run-1] == NULL )
      return 0;
   else
      return reopt->soltree->nsols[run-1];
}

/* returns number of all solutions */
int SCIPreoptGetNSols(
   SCIP_REOPT*           reopt                    /**< reopt data */
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

/* return the stored solutions of a given run */
SCIP_RETCODE SCIPreoptGetSolsRun(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   int                   run,                     /**< number of the run (1,2,...) */
   SCIP_SOL**            sols,                    /**< array of solutions to fill */
   int                   allocmem,                /**< length of the array */
   int*                  nsols                    /**< pointer to store the number of added solutions */
   )
{
   int s;

   assert(reopt != NULL);
   assert(run > 0 && run <= reopt->run);
   assert(sols != NULL);

   assert(allocmem > 0);
   *nsols = 0;

   for(s = 0; s < MIN(reopt->soltree->nsols[run-1], allocmem); s++)
   {
      if( !reopt->soltree->sols[run-1][s]->updated )
         (*nsols)++;
   }

   if( allocmem < (*nsols) )
      return SCIP_OKAY;

   (*nsols) = 0;
   for(s = 0; s < reopt->soltree->nsols[run-1]; s++)
   {
      if( !reopt->soltree->sols[run-1][s]->updated )
      {
         sols[*nsols] = reopt->soltree->sols[run-1][s]->sol;
         reopt->soltree->sols[run-1][s]->updated = TRUE;
         (*nsols)++;
      }
   }

   return SCIP_OKAY;
}

/* returns the number of saved solutions overall runs */
int SCIPreoptNSavedSols(
    SCIP_REOPT*           reopt                   /**< reopt data */
   )
{
   int nsavedsols;

   assert(reopt != NULL);
   assert(reopt->soltree->root != NULL);

   nsavedsols = 0;

   if( reopt->soltree->root->lchild != NULL
    || reopt->soltree->root->rchild != NULL)
      nsavedsols = soltreeNInducedtSols(reopt->soltree->root);

   return nsavedsols;
}

/* returns the number of reused sols over all runs */
int SCIPreoptNUsedSols(
   SCIP_REOPT*           reopt                    /**< reopt data */
   )
{
   int nsolsused;

   assert(reopt != NULL);

   nsolsused = 0;

   if( reopt->soltree->root != NULL )
      nsolsused = soltreeGetNUsedSols(reopt->soltree->root);

   return nsolsused;
}

/* save objective function */
SCIP_RETCODE SCIPreoptSaveNewObj(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_VAR**            vars,                    /**< array od variables */
   int                   nvars                    /**< number of variables */
   )
{
   SCIP_Real norm;
   int v;
   int id;

   assert(reopt != NULL);

   /* check memory */
   SCIP_CALL( ensureRunSize(reopt, set, reopt->run, blkmem) );

   norm = 0;

   /* get memory */
   SCIP_ALLOC( BMSallocMemoryArray(&reopt->objs[reopt->run-1], reopt->nobjvars) );

   /* save coefficients */
   for(v = 0; v < reopt->nobjvars; v++)
   {
      id = SCIPvarGetIndex(vars[v]);
      reopt->objs[reopt->run-1][id] = SCIPvarGetObj(vars[id]);
      norm += (SCIPvarGetObj(vars[id]) * SCIPvarGetObj(vars[id]));

      /* mark this objective as the first non empty */
      if( reopt->firstobj == -1 && reopt->objs[reopt->run-1][id] != 0 )
         reopt->firstobj = reopt->run-1;
   }
   assert(norm >= 0);
   norm = sqrt(norm);

   /* normalize the coefficients */
   for(v = 0; v < reopt->nobjvars && norm > 0; v++)
   {
      id = SCIPvarGetIndex(vars[v]);
      reopt->objs[reopt->run-1][id] /= norm;
   }

   /* calculate similarity to last objective */
   if( reopt->run-1 > 1 )
   {
      /* calculate similarity to first objective */
      if( reopt->run-1 > 1 && reopt->firstobj < reopt->run-1 )
         reopt->simtofirstobj = reoptSimilarity(reopt, reopt->run-1, reopt->firstobj);

      /* calculate similarity to last objective */
      reopt->simtolastobj = reoptSimilarity(reopt, reopt->run-1, reopt->run-2);

      SCIPdebugMessage("new objective has similarity of %g/%g compared to first/previous.\n", reopt->simtofirstobj, reopt->simtolastobj);
      printf("new objective has similarity of %g/%g compared to first/previous.\n", reopt->simtofirstobj, reopt->simtolastobj);
   }

   SCIPdebugMessage("saved obj for run %d.\n", reopt->run);

   return SCIP_OKAY;
}

/* check if the current and the previous objective are similar enough
 * returns TRUE if we want to restart, otherwise FALSE */
SCIP_RETCODE SCIPreoptCheckRestart(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   SCIP_Bool restart;
   SCIP_Real sim;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);

   sim = 1.0;
   restart = FALSE;

   if( reopt->run > 0 && set->reopt_delay > -1.0 )
   {
      sim = reopt->simtolastobj;
   }

   if( SCIPsetIsFeasLT(set, sim, set->reopt_delay) )
   {
      SCIPdebugMessage("-> restart reoptimization (objective functions are not similar enough)\n");
      restart = TRUE;
   }
   else if( reopt->reopttree->nsavednodes > reopt->maxsavednodes )
   {
      SCIPdebugMessage("-> restart reoptimization (node limit reached)\n");
      restart = TRUE;
   }
   else if( reopt->noptsolsbyreoptsol >= reopt->forceheurrestart )
   {
      SCIPdebugMessage("-> restart reoptimization (found last %d optimal solutions by <reoptsols>)\n", reopt->noptsolsbyreoptsol);
      printf("-> restart reoptimization (found last %d optimal solutions by <reoptsols>)\n", reopt->noptsolsbyreoptsol);
      reopt->noptsolsbyreoptsol = 0;
      restart = TRUE;
   }

   if( restart )
   {
      SCIP_CALL( reoptRestart(reopt, blkmem) );
   }

   return SCIP_OKAY;
}

/*
 * returns the similarity to the previous objective function,
 * if no objective functions are saved the similarity is -2.0.
 */
SCIP_Real SCIPreoptGetSimToPrevious(
      SCIP_REOPT*        reopt                    /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);
   return reopt->simtolastobj;
}

/*
 * returns the similarity to the first objective function,
 * if no objective functions are saved the similarity is -2.0.
 */
SCIP_Real SCIPreoptGetSimToFirst(
      SCIP_REOPT*        reopt                    /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);
   return reopt->simtofirstobj;
}

/*
 * return the similarity between two of objective functions of two given runs
 */
SCIP_Real SCIPreoptGetSim(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   int                   run1,                    /**< number of the first run */
   int                   run2                     /**< number of the second run */
   )
{
   assert(reopt != NULL);
   assert(run1 > 0 && run1 <= reopt->run);
   assert(run2 > 0 && run2 <= reopt->run);

   return reoptSimilarity(reopt, run1-1, run2-1);
}

/*
 * returns the best solution of the last run
 */
SCIP_SOL* SCIPreoptGetLastBestSol(
   SCIP_REOPT*           reopt                    /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);
   assert(reopt->lastbestsol != NULL);

   if( reopt->run-2 < 0 )
      return NULL;
   else
   {
      assert(reopt->lastbestsol[reopt->run-2] != NULL);
      return reopt->lastbestsol[reopt->run-2];
   }
}

/*
 * returns the coefficient of variable with index @param idx in run @param run
 */
SCIP_Real SCIPreoptGetObjCoef(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   int                   run,                     /**< number of the run */
   int                   idx                      /**< index of variable */
   )
{
   assert(reopt != NULL);
   assert(0 <= run-1 && run < reopt->runsize);

   return reopt->objs[run-1][idx];
}

/* checks the changes of the objective coefficients */
void SCIPreoptGetVarCoefChg(
   SCIP_REOPT*           reopt,                   /**< reopt data */
   int                   varidx,                  /**< index of the variable */
   SCIP_Bool*            negated,                 /**< variable switches the impact, i.e., coefficient changed the sign */
   SCIP_Bool*            entering,                /**< variable enters the objective, i.e., coefficient gets different to 0 */
   SCIP_Bool*            leaving                  /**< variable leaves the objective, i.e., coefficient gets 0 */
   )
{
   assert(reopt != NULL);
   assert(varidx >= 0 && varidx < reopt->nobjvars);

   *negated = FALSE;
   *entering = FALSE;
   *leaving = FALSE;

   if( reopt->objs[reopt->run-2] == NULL || reopt->run-2 <= 0)
      return;

   /* variable has objective coefficients with opposed sign */
   if( reopt->objs[reopt->run-1] != NULL && reopt->run >= 1 )
   {
      *negated = (SCIP_Real)reopt->objs[reopt->run-1][varidx]/reopt->objs[reopt->run-2][varidx] < 0 ? TRUE : FALSE;
   }
   /* variable leaves the objective */
   else if( reopt->objs[reopt->run-2][varidx] == 0 && reopt->objs[reopt->run-3][varidx] != 0 )
   {
      *leaving = TRUE;
   }
   /* variable enters the objective function */
   else if( reopt->objs[reopt->run-2][varidx] != 0 && reopt->objs[reopt->run-3][varidx] == 0 )
   {
      *entering = TRUE;
   }

   return;
}

/*
 * print optimal solutions of all previous runs
 */
SCIP_RETCODE SCIPreoptPrintOptSols(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_REOPT*           reopt                    /**< reoptimization data structure */
   )
{
   int run;

   assert(scip != NULL);
   assert(reopt != NULL);
   assert(reopt->run > 0);
   assert(reopt->lastbestsol != NULL);

   printf(">> best %d solutions:\n", reopt->run-1);

   for(run = reopt->run-1; run >= 0; run--)
   {
      SCIP_SOL* tmp_sol;

      assert(reopt->lastbestsol != NULL);
      assert(reopt->lastbestsol[run] != NULL);

      SCIP_CALL( SCIPcreateSolCopy(scip, &tmp_sol, reopt->lastbestsol[run]) );

      printf(">> optimal solution of run %d:\n", run);
      SCIP_CALL( SCIPprintSol(scip, tmp_sol, NULL, FALSE) );
      printf("\n");

      SCIP_CALL( SCIPfreeSol(scip, &tmp_sol) );
   }

   return SCIP_OKAY;
}

/*
 * return all optimal solutions of the previous runs
 * depending on the current stage the method copies the solutions into
 * the origprimal or primal space. That means, all solutions need to be
 * freed before starting a new iteration!!!
 */
SCIP_RETCODE SCIPreoptGetOptSols(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_REOPT*           reopt,                   /**< reopt data */
   SCIP_SOL**            sols,                    /**< array to store solution pointers */
   int*                  nsols                    /**< pointer to store number of solutions */
   )
{
   SCIP_SOL* sol;
   int run;

   assert(scip != NULL);
   assert(reopt != NULL);
   assert(sols != NULL);

   for(run = 0; run < reopt->run; run++)
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &sol, reopt->lastbestsol[run]) );
      sols[run] = sol;
   }

   return SCIP_OKAY;
}

/* reset marks of stored solutions to not updated */
void SCIPreoptResetSolMarks(
   SCIP_REOPT*           reopt                    /**< reoptimization data structure */
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

/* returns the number of stored nodes in the subtree induced by @param node */
int SCIPreoptGetNNodes(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node                     /**< node of the search tree */
   )
{
   assert(reopt != NULL);

   if( node == NULL || SCIPnodeGetDepth(node) == 0 )
      return reopt->reopttree->nsavednodes + 1;

   if( SCIPnodeGetReoptID(node) == -1 )
      return 0;

   assert(0 < SCIPnodeGetReoptID(node) && SCIPnodeGetReoptID(node) < reopt->reopttree->allocmemnodes);

   return reopttreeGetNNodes(reopt->reopttree, SCIPnodeGetReoptID(node));
}

/* returns the number of child nodes of @param node */
int SCIPreoptGetNChildren(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node                     /**< node of the search tree */
   )
{
   assert(reopt != NULL);

   if( node == NULL || SCIPnodeGetDepth(node) == 0 )
      return reopt->reopttree->reoptnodes[0]->nchilds;

   if( SCIPnodeGetReoptID(node) == -1 )
      return 0;

   assert(0 < SCIPnodeGetReoptID(node) && SCIPnodeGetReoptID(node) < reopt->reopttree->allocmemnodes);

   return reopt->reopttree->reoptnodes[SCIPnodeGetReoptID(node)]->nchilds;
}

/* returns the number of leave nodes of the subtree induced
 * by @param node (of the whole tree of node == NULL) */
int SCIPreoptGetNLeaves(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node                     /**< node of the search tree */
   )
{
   int nleaves;
   int id;
   int i;

   assert(reopt != NULL);

   nleaves = 0;
   id = (node == NULL) ? 0 : SCIPnodeGetReoptID(node);

   /* the node is a leaf node */
   if( id == -1 )
      return nleaves;

   for(i = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++)
   {
      int childid;

      childid = reopt->reopttree->reoptnodes[id]->childids[i];

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
         nleaves++;
      else
         nleaves += reoptGetNLeaves(reopt, childid);
   }

   return nleaves;
}

/*
 *  Save information if infeasible nodes
 */
SCIP_RETCODE SCIPreoptAddInfNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node                     /**< node of the search tree */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);

   if( reopt->sepasubtreesglb )
   {
      SCIP_CALL( saveGlobalCons(reopt, set, blkmem, node, REOPT_CONSTYPE_INFSUBTREE) );
   }

   reopt->reopttree->ninfeasnodesround++;
   reopt->reopttree->ninfeasnodes++;

   return SCIP_OKAY;
}

/**
 * check the reason for cut off a node and if necessary store the node
 */
SCIP_RETCODE SCIPreoptCheckCutoff(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memery */
   SCIP_NODE*            node,                    /**< node of the search tree */
   SCIP_EVENT*           event,                   /**< event */
   SCIP_LPSOLSTAT        lpsolstat,               /**< solution status of the LP */
   SCIP_Bool             isrootnode,              /**< the node is the root */
   SCIP_Bool             isfocusnode,             /**< the node is the current focus node */
   SCIP_Real             lowerbound,              /**< lower bound of the node */
   int                   effectiverootdepth       /**< effective root depth */
   )
{
   SCIP_EVENTTYPE eventtype;
   SCIP_Bool strongbranched;

   assert(reopt != NULL);
   assert(node != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEBRANCHED
       || SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEFEASIBLE
       || SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEINFEASIBLE);
   assert(SCIPeventGetNode(event) == node);

   if( reopt->lastseennode == SCIPnodeGetNumber(node) )
      return SCIP_OKAY;

   reopt->lastseennode = SCIPnodeGetNumber(node);

   eventtype = SCIPeventGetType(event);

   SCIPdebugMessage("catch event %x for node %lld\n", eventtype, SCIPnodeGetNumber(node));

   /* case 1: the current node is the root node
    * we can skip if the root is (in)feasible or branched w/o bound
    * changes based on dual information.
    *
    * case 2: we need to store the current node if it contains
    * bound changes based on dual information or is a leave node */

   if( isrootnode )
   {
      if( SCIPreoptGetNDualBndchs(reopt, node) > 0 )
      {
         goto CHECK;
      }
      else if( eventtype == SCIP_EVENTTYPE_NODEBRANCHED )
      {
         /* store or update the information */
         SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_TRANSIT, TRUE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
      }
      else if( eventtype == SCIP_EVENTTYPE_NODEFEASIBLE )
      {
         /* delete saved dual information which would lead to split the node in a further iteration */
         SCIPreoptResetDualcons(reopt, node, blkmem);

         /* store or update the information */
         SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_FEASIBLE, FALSE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
      }
      else if( eventtype == SCIP_EVENTTYPE_NODEINFEASIBLE )
      {
         /* delete saved dual information which would lead to split the node in a further iteration */
         SCIPreoptResetDualcons(reopt, node, blkmem);

         /* store or update the information */
         SCIP_CALL( addNode(reopt, set, blkmem, node, reopt->currentnode == 1 ? SCIP_REOPTTYPE_INFSUBTREE : SCIP_REOPTTYPE_PRUNED, FALSE,
               isrootnode, isfocusnode, lpsolstat, lowerbound) );
      }

      assert(reopt->currentnode == -1);
      assert(reopt->dualcons == NULL || reopt->dualcons->nvars == 0);

      return SCIP_OKAY;
   }

   CHECK:

   if( effectiverootdepth == SCIPnodeGetDepth(node) )
   {
      strongbranched = SCIPreoptGetNDualBndchs(reopt, node) > 0 ? TRUE : FALSE;
   }
   else
   {
      strongbranched = SCIPnodeGetNDualBndchgs(node) > 0 ? TRUE : FALSE;
   }

   SCIPdebugMessage("check the reason of cutoff for node %lld:\n", SCIPnodeGetNumber(node));
   SCIPdebugMessage(" -> focusnode: %s\n", isfocusnode ? "yes" : "no");
   SCIPdebugMessage(" -> depth: %d, eff. root depth: %d\n", SCIPnodeGetDepth(node), effectiverootdepth);
   SCIPdebugMessage(" -> strong branched: %s\n", strongbranched ? "yes" : "no");
   SCIPdebugMessage(" -> LP lpsolstat     : %d\n", lpsolstat);

   switch( SCIPeventGetType(event) ) {
      case SCIP_EVENTTYPE_NODEFEASIBLE:
         /** current node has to be the eventnode */
         assert(isfocusnode);

         SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_FEASIBLE);

         /* delete strong branching information of some exists */
         deleteLastDualBndchgs(reopt);

         SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_FEASIBLE, FALSE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
         break;

      case SCIP_EVENTTYPE_NODEINFEASIBLE:
         /**
          * we have to check if the current node is the event node.
          * if the current node is not the event node, we have to save this node, else we have to
          * look at LP lpsolstat and decide.
          */
         if( isfocusnode )
         {
            /**
             * an after-branch heuristic says NODEINFEASIBLE, maybe the cutoff bound is reached.
             * because the node is already branched we have all children and can delete this node.
             */
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
                   * sense of effective root depth */
                  if( SCIPnodeGetDepth(node) > effectiverootdepth )
                  {
                     SCIP_CALL( SCIPreoptAddDualBndchg(reopt, set, blkmem, node, NULL, 0, 1) );
                  }

                  SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_INFSUBTREE);
                  SCIPdebugMessage(" -> new constraint of type: %d\n", REOPT_CONSTYPE_INFSUBTREE);

                  /* save the node as a strong branched node */
                  SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_INFSUBTREE, FALSE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
               }
               else
               {
                  assert(SCIP_LPSOLSTAT_OBJLIMIT || SCIP_LPSOLSTAT_OPTIMAL || SCIP_LPSOLSTAT_NOTSOLVED);

                  SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_PRUNED);

                  /* delete strong branching information of some exists */
                  deleteLastDualBndchgs(reopt);

                  SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_PRUNED, FALSE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
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
                  SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_INFEASIBLE);
                  SCIP_CALL( SCIPreoptAddInfNode(reopt, set, blkmem, node) );
               }
               else
               {
                  SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_PRUNED);

                  /* store the node */
                  SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_PRUNED, TRUE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
               }
            }
         }
         else
         {
            SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_PRUNED);

            /* if the node was created by branch_nodereopt, nothing happens */
            SCIP_CALL(addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_PRUNED, TRUE, isrootnode, isfocusnode, lpsolstat, lowerbound) );

         }
         break;

      case SCIP_EVENTTYPE_NODEBRANCHED:
         /** current node has to be the eventnode */
         assert(isfocusnode);

         reopt->lastbranched = SCIPnodeGetNumber(node);

         /**
          * we have to check the depth of the current node. if the depth is equal to the effective
          * root depth, then all information about bound changes based on dual information already exists, \
          * else we have to look at the domchg-data-structure.
          */
         if (SCIPnodeGetDepth(node) == effectiverootdepth)
         {
            /*
             * Save the node if there are added constraints, because this means the node is a copy create by the
             * reoptimization plug-in and contains at least one logic-or-constraint
             */
            if( strongbranched )
            {
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_STRBRANCHED);
               SCIPdebugMessage(" -> new constraint of type: %d\n", REOPT_CONSTYPE_STRBRANCHED);
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_STRBRANCHED, TRUE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
            }
            else if( SCIPreoptGetNAddedConss(reopt, node) > 0 )
            {
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_LOGICORNODE);
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_LOGICORNODE, TRUE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
            }
            else
            {
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_TRANSIT);
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_TRANSIT, TRUE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
            }
         }
         else
         {
            /**
             * we only branch on binary variables and var == NULL indicates memory allocation w/o saving information.
             *
             * we have to do this in the following order:
             * 1) all bound-changes are local, thats way we have to mark the node to include bound changes based
             *    on dual information.
             * 2) save or update the node.
             */
            if( strongbranched )
            {
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_STRBRANCHED);
               SCIPdebugMessage(" -> new constraint of type: %d\n", REOPT_CONSTYPE_STRBRANCHED);
               SCIP_CALL( SCIPreoptAddDualBndchg(reopt, set, blkmem, node, NULL, 0, 1) );
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_STRBRANCHED, TRUE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
            }
            else if( SCIPreoptGetNAddedConss(reopt, node) > 0 )
            {
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_LOGICORNODE);
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_LOGICORNODE, TRUE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
            }
            else
            {
               SCIPdebugMessage(" -> new reopttype: %d\n", SCIP_REOPTTYPE_TRANSIT);
               SCIP_CALL( addNode(reopt, set, blkmem, node, SCIP_REOPTTYPE_TRANSIT, TRUE, isrootnode, isfocusnode, lpsolstat, lowerbound) );
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

/** store bound changes based on dual information */
SCIP_RETCODE SCIPreoptAddDualBndchg(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   SCIP_VAR*             var,                     /**< variables */
   SCIP_Real             newval,                  /**< new bound */
   SCIP_Real             oldval                   /**< old bound */
   )
{
   SCIP_Real constant;
   SCIP_Real scalar;

   assert(reopt != NULL);
   assert(node != NULL);

   constant = 0;
   scalar = 1;

   /**
    * If var == NULL, we save all information by calling SCIPreoptNodeFinished().
    * In that case, all bound changes were not global and we can find them within the
    * domchg data structure.
    * Otherwise, we allocate memory and store the information.
    */
   if( var != NULL )
   {
      int allocmem;

      assert(SCIPsetIsFeasEQ(set, newval, 0) || SCIPsetIsFeasEQ(set, newval, 1));

      allocmem = (reopt->dualcons == NULL || reopt->dualcons->allocmem == 0) ? DEFAULT_MEM_DUALCONS : reopt->dualcons->allocmem+2;

      /* allocate memory of necessary */
      SCIP_CALL( checkMemDualCons(reopt, blkmem, allocmem) );

      assert(reopt->dualcons->allocmem > 0);
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
      reopt->dualcons->nvars++;

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

/* returns the number of bound changes based on dual information */
int SCIPreoptGetNDualBndchs(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node                     /**< node of the search tree */
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

/* returns the child nodes of @param node that need to be
 * reoptimized next or NULL if @param node is a leaf */
SCIP_RETCODE SCIPreoptGetChildIDs(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   int*                  childs,                  /**< array to store the child ids */
   int                   mem,                     /**< allocated memory */
   int*                  nchilds                  /**< pointer to store the number of child nodes */
   )
{
   SCIP_Bool runagain;
   int id;

   assert(reopt != NULL);
   assert(node == NULL || SCIPnodeGetReoptID(node) >= 0);
   assert(mem > 0 && childs != NULL);

   if( node == NULL )
      id = 0;
   else
      id = SCIPnodeGetDepth(node) == 0 ? 0 : SCIPnodeGetReoptID(node);

   (*nchilds) = 0;

   assert(reopt->reopttree->reoptnodes[id] != NULL);

   /* check if there are redundant bound changes or infeasible nodes */
   runagain = TRUE;

   while( runagain && reopt->reopttree->reoptnodes[id]->nchilds > 0 )
   {
      SCIP_CALL( dryBranch(reopt, set, blkmem, &runagain, id) );
   }

   /* return the list of child nodes of some exists;
    * otherwise return NULL */
   if( reopt->reopttree->reoptnodes[id]->childids != NULL && reopt->reopttree->reoptnodes[id]->nchilds > 0 )
   {
      int c;

      (*nchilds) = reopt->reopttree->reoptnodes[id]->nchilds;

      if( mem < *nchilds )
         return SCIP_OKAY;

      for(c = 0; c < *nchilds; c++)
      {
         childs[c] = reopt->reopttree->reoptnodes[id]->childids[c];
      }
   }

   return SCIP_OKAY;
}

/* returns all leave if of the subtree induced by @param node */
SCIP_RETCODE SCIPreoptGetLeaves(
   SCIP_REOPT*           reopt,                   /**< reoptimization data */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   int*                  leaves,                  /**< array to the the ids */
   int                   mem,                     /**< allocated memory */
   int*                  nleaves                  /**< pointer to store the number of leav node */
   )
{
   int id;
   int i;

   assert(reopt != NULL);
   assert(node == NULL || SCIPnodeGetReoptID(node) >= 0);
   assert(mem > 0 && leaves != NULL);
   assert((*nleaves) >= 0);

   if( node == NULL )
      id = 0;
   else
      id = SCIPnodeGetDepth(node) == 0 ? 0 : SCIPnodeGetReoptID(node);

   assert(reopt->reopttree->reoptnodes[id] != NULL);

   for(i = 0; i < reopt->reopttree->reoptnodes[id]->nchilds; i++)
   {
      int childid;

      assert(*nleaves + 1 <= mem);

      childid = reopt->reopttree->reoptnodes[id]->childids[i];

      if( reopt->reopttree->reoptnodes[childid]->nchilds == 0 )
      {
         leaves[(*nleaves)] = reopt->reopttree->reoptnodes[id]->childids[i];
         (*nleaves)++;
      }
      else
      {
         SCIP_CALL( reoptGetLeaves(reopt, set, blkmem, childid, &leaves[*nleaves], mem, nleaves) );
      }
   }

   return SCIP_OKAY;
}

/* add the node @param node to the reopttree */
SCIP_RETCODE SCIPreoptAddNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   SCIP_REOPTTYPE        reopttype,               /**< reopttype of the node */
   SCIP_LPSOLSTAT        lpsolstat,               /**< solution status of the LP */
   SCIP_Bool             saveafterduals,          /**< store decisions directly after the first dual decision */
   SCIP_Bool             isrootnode,              /**< the node is the root */
   SCIP_Bool             isfocusnode,             /**< the node is the current focus node */
   SCIP_Real             lowerbound               /**< lower bound of the node */
   )
{
   assert(reopt != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);

   SCIP_CALL( addNode(reopt, set, blkmem, node, reopttype, saveafterduals, isrootnode, isfocusnode, lpsolstat, lowerbound) );

   return SCIP_OKAY;
}

/* add all unprocessed nodes to the reoptimization tree */
SCIP_RETCODE SCIPreoptSaveOpenNodes(
   SCIP_REOPT*           reopt,                   /* reoptimization data structure */
   SCIP_SET*             set,                     /* global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /* block memory */
   SCIP_NODE**           leaves,                  /* array of open leave nodes */
   int                   nleaves,                 /* number of open leave nodes */
   SCIP_NODE**           childs,                  /* array of open children nodes */
   int                   nchilds,                 /* number of open leave nodes */
   SCIP_NODE**           siblings,                /* array of open sibling nodes */
   int                   nsiblings                /* number of open leave nodes */
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
      SCIP_CALL( addNode(reopt, set, blkmem, leaves[n], SCIP_REOPTTYPE_PRUNED, FALSE,
            FALSE, FALSE, SCIP_LPSOLSTAT_NOTSOLVED, SCIPnodeGetLowerbound(leaves[n])) );
   }

   /* save open children */
   for( n = 0; n < nchilds; n++ )
   {
      SCIP_CALL( addNode(reopt, set, blkmem, childs[n], SCIP_REOPTTYPE_PRUNED, FALSE,
            FALSE, FALSE, SCIP_LPSOLSTAT_NOTSOLVED, SCIPnodeGetLowerbound(childs[n])) );
   }

   /* save open siblings */
   for( n = 0; n < nsiblings; n++ )
   {
      SCIP_CALL( addNode(reopt, set, blkmem, siblings[n], SCIP_REOPTTYPE_PRUNED, FALSE,
            FALSE, FALSE, SCIP_LPSOLSTAT_NOTSOLVED, SCIPnodeGetLowerbound(siblings[n])) );
   }

   return SCIP_OKAY;
}

/* calculates a local similarity of a given node and returns if the subproblem
 * should be solved from scratch */
SCIP_RETCODE SCIPreoptCheckLocalRestart(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_NODE*            node,                    /**< node of the search tree */
   SCIP_VAR**            origvars,                /**< original variables */
   int                   norigvars,               /**< number of original variables */
   SCIP_Bool*            localrestart             /**< pointer to store if we want to restart solving the (sub)problem */
   )
{
   int id;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(origvars != NULL);
   assert(reopt->nobjvars == norigvars);

   *localrestart = FALSE;
   id = SCIPnodeGetDepth(node) == 0 ? 0 : SCIPnodeGetReoptID(node);

   if( id > -1 && reopt->localdelay > -1 )
   {
      SCIP_Real sim;
      SCIP_Real scalar;
      SCIP_Real oldnorm;
      SCIP_Real newnorm;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real oldcoef;
      SCIP_Real newcoef;
      int v;
      int vid;

      if( id == 0 )
         reopt->nlocalrestarts = 0;

      sim = 0.0;
      scalar = 0.0;
      oldnorm = 0.0;
      newnorm = 0.0;

      /* dot-product and norm */
      for(v = 0; v < norigvars; v++)
      {
         lb = SCIPvarGetLbLocal(origvars[v]);
         ub = SCIPvarGetUbLocal(origvars[v]);

         if( SCIPsetIsFeasLT(set, lb, ub) )
         {
            vid = SCIPvarGetIndex(origvars[v]);
            oldcoef = SCIPreoptGetObjCoef(reopt, reopt->run-1, vid);
            newcoef = SCIPreoptGetObjCoef(reopt, reopt->run, vid);

            scalar += (oldcoef * newcoef);
            oldnorm += pow(oldcoef, 2);
            newnorm += pow(newcoef, 2);
         }
      }

      /* normalize the dot-product */
      if( newnorm == 0 || oldnorm == 0 || scalar == 0 )
         sim = 0.0;
      else
         sim = scalar/(sqrt(oldnorm)*sqrt(newnorm));

      /* delete the stored subtree and information about bound changes
       * based on dual information */
      if( SCIPsetIsLT(set, sim, reopt->localdelay) )
      {
         /* set the flag */
         *localrestart = TRUE;

         reopt->nlocalrestarts++;

         /* delete the stored subtree */
         SCIP_CALL( deleteChildrenBelow(reopt->reopttree, blkmem, id, FALSE, FALSE) );

         /* delete the stored constraints */
         if( reopt->reopttree->reoptnodes[id]->dualfixing )
         {
            if( reopt->reopttree->reoptnodes[id]->dualconscur != NULL )
            {
               BMSfreeBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconscur->vars, reopt->reopttree->reoptnodes[id]->dualconscur->allocmem);
               BMSfreeBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconscur->vals, reopt->reopttree->reoptnodes[id]->dualconscur->allocmem);
               BMSfreeMemory(&reopt->reopttree->reoptnodes[id]->dualconscur);
               reopt->reopttree->reoptnodes[id]->dualconscur = NULL;
            }

            if( reopt->reopttree->reoptnodes[id]->dualconsnex != NULL )
            {
               BMSfreeBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconsnex->vars, reopt->reopttree->reoptnodes[id]->dualconsnex->allocmem);
               BMSfreeBlockMemoryArray(blkmem, &reopt->reopttree->reoptnodes[id]->dualconsnex->vals, reopt->reopttree->reoptnodes[id]->dualconsnex->allocmem);
               BMSfreeMemory(&reopt->reopttree->reoptnodes[id]->dualconsnex);
               reopt->reopttree->reoptnodes[id]->dualconsnex = NULL;
            }

            reopt->reopttree->reoptnodes[id]->dualfixing = FALSE;
            reopt->reopttree->reoptnodes[id]->reopttype = SCIP_REOPTTYPE_LEAF;
         }
      }

      SCIPdebugMessage(" -> local similarity: %.4f%s\n", sim, *localrestart ? " (solve subproblem from scratch)" : "");
   }

   return SCIP_OKAY;
}

/* reset the complete tree an set a new search frontier */
SCIP_RETCODE SCIPreoptApplyCompression(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   SCIP_REOPTNODE**      representatives,         /**< array of representatives */
   int                   nrepresentatives,        /**< number of representatives */
   SCIP_Bool*            success                  /**< pointer to store if the method was successful */
   )
{
   SCIP_REOPTTREE* reopttree;
   int id;
   int r;

   assert(reopt != NULL);
   assert(blkmem != NULL);
   assert(representatives != NULL);
   assert(nrepresentatives > 0);

   reopttree = reopt->reopttree;

   /* reset the current search tree */
   SCIP_CALL( reoptResetTree(reopt, blkmem) );
   assert(reopttree->nsavednodes == 0);

   /* create a new root node */
   id = 0;
   SCIP_CALL( createReoptnode(reopttree, id) );

   /* set the reopttype */
   reopttree->reoptnodes[0]->reopttype = SCIP_REOPTTYPE_TRANSIT;

   /* add all representatives */
   for(r = 0; r < nrepresentatives; r++)
   {
      /* get an empty slot*/
      id = (int) (size_t) SCIPqueueRemove(reopttree->openids);
      assert(1 <= id && id < reopttree->allocmemnodes);

      /* ensure that the slot is empty */
      if( reopttree->reoptnodes[id] != NULL )
      {
         assert(reopttree->reoptnodes[id]->nvars == 0);
         assert(reopttree->reoptnodes[id]->nconss == 0);
         assert(reopttree->reoptnodes[id]->nchilds == 0);

         /* hard reset */
         SCIP_CALL( reopttreeDeleteNode(reopttree, blkmem, id, TRUE) );
      }

      SCIP_CALL( createReoptnode(reopttree, id) );
      assert(reopttree->reoptnodes[id] != NULL);

      /* set the new node
       * 1. copy all variables, bounds, and boundtypes
       * 2. copy all constraints
       * 3. set the parent relation */
      if( representatives[r]->nvars > 0 )
      {
         assert(representatives[r]->nvars <= representatives[r]->allocvarmem);
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->vars, representatives[r]->vars, representatives[r]->nvars) );
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->varbounds, representatives[r]->varbounds, representatives[r]->nvars) );
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->varboundtypes, representatives[r]->varboundtypes, representatives[r]->nvars) );
         reopttree->reoptnodes[id]->allocvarmem = representatives[r]->allocvarmem;
         reopttree->reoptnodes[id]->nvars = representatives[r]->nvars;
      }

      if( representatives[r]->nconss > 0 )
      {
         int c;

         assert(representatives[r]->nconss <= representatives[r]->allocmemconss);

         for( c = 0; c < representatives[r]->nconss; c++ )
         {
            SCIP_CALL( SCIPreoptnodeAddCons(reopttree->reoptnodes[id], representatives[r]->conss[c]->vars,
                  representatives[r]->conss[c]->vals, representatives[r]->conss[c]->nvars, representatives[r]->conss[c]->constype, blkmem) );
         }
      }

      reopttree->reoptnodes[id]->parentID = representatives[r]->parentID;

      assert(reopttree->reoptnodes[id]->parentID == 0);
      assert(reopttree->reoptnodes[id]->nvars >= 0);
      assert(reopttree->reoptnodes[id]->nvars <= reopttree->reoptnodes[id]->allocvarmem);
      assert(reopttree->reoptnodes[id]->nconss >= 0);

      /* set the reopttype */
      if( reopttree->reoptnodes[id]->nconss == 0 )
         reopttree->reoptnodes[id]->reopttype = SCIP_REOPTTYPE_LEAF;
      else
         reopttree->reoptnodes[id]->reopttype = SCIP_REOPTTYPE_LOGICORNODE;

      /* add the representative as a child of the root */
      SCIP_CALL( reoptAddChild(reopttree, 0, id, blkmem) );
   }

   SCIPdebugMessage("-> new tree consists of %d nodes, the root has %d child nodes.\n",
         reopttree->nsavednodes, reopttree->reoptnodes[0]->nchilds);

   (*success) = TRUE;

   return SCIP_OKAY;
}

/* returns if a node need to be split because some bound changes
 * were based on dual information */
SCIP_Bool SCIPreoptSplitNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node                     /**< node of the search tree */
   )
{
   int id;

   assert(reopt != NULL);
   assert(node != NULL);

   id = SCIPnodeGetReoptID(node);

   assert(0 <= id && id < reopt->reopttree->allocmemnodes);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   if( reopt->reopttree->reoptnodes[id]->dualfixing )
   {
      assert(reopt->reopttree->reoptnodes[id]->dualconscur != NULL);
      assert(reopt->reopttree->reoptnodes[id]->dualconscur->nvars > 0);

      return TRUE;
   }

   return FALSE;
}

/**
 * create the constraint which splits the node stored at ID id on the basis of
 * the stored dual information.
 */
void SCIPreoptCreateSplitCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   id,                      /**< id of stored node */
   LOGICORDATA*          consdata                 /**< logicor data*/
   )
{
   assert(reopt != NULL);
   assert(consdata != NULL);
   assert(consdata->allocmem > 0);
   assert(consdata->vars != NULL);
   assert(consdata->vals != NULL);
   assert(consdata->nvars == 0);
   assert(0 <= id && id < reopt->reopttree->allocmemnodes);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   /* copy the variable information */
   if( reopt->reopttree->reoptnodes[id]->dualconscur != NULL
    && consdata->allocmem >= reopt->reopttree->reoptnodes[id]->dualconscur->nvars )
   {
      int v;
      for(v = 0; v < reopt->reopttree->reoptnodes[id]->dualconscur->nvars; v++)
      {
         consdata->vars[v] = reopt->reopttree->reoptnodes[id]->dualconscur->vars[v];
         consdata->vals[v] = reopt->reopttree->reoptnodes[id]->dualconscur->vals[v];
      }
      consdata->nvars = reopt->reopttree->reoptnodes[id]->dualconscur->nvars;
      consdata->constype = reopt->reopttree->reoptnodes[id]->dualconscur->constype;
   }

   return;
}

/* split the root node and move all children to one of the two resulting nodes */
SCIP_RETCODE SCIPreoptSplitRoot(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   SCIP_REOPTTREE* reopttree;
   LOGICORDATA* consdata;
   int dummy1;
   int dummy2;
   int v;
   int nbndchgs;
   int nchilds;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(reopt->reopttree->reoptnodes[0] != NULL);
   assert(reopt->reopttree->reoptnodes[0]->dualfixing);
   assert(reopt->reopttree->reoptnodes[0]->reopttype == SCIP_REOPTTYPE_STRBRANCHED);
   assert(set != NULL);
   assert(blkmem != NULL);

   reopttree = reopt->reopttree;

   nchilds = reopttree->reoptnodes[0]->nchilds;

   assert(reopttree->reoptnodes[0]->dualconscur != NULL);
   nbndchgs = reopttree->reoptnodes[0]->dualconscur->nvars;

   /* ensure that two free slots are available  */
   SCIP_CALL( reopttreeCheckMemory(reopttree, blkmem) );
   dummy1 = (int) (size_t) SCIPqueueRemove(reopttree->openids);

   SCIP_CALL( reopttreeCheckMemory(reopttree, blkmem) );
   dummy2 = (int) (size_t) SCIPqueueRemove(reopttree->openids);

   assert(dummy1 > 0 && dummy2 > 0);
   assert(reopttree->reoptnodes[dummy1] == NULL || reopttree->reoptnodes[dummy1]->nvars == 0);
   assert(reopttree->reoptnodes[dummy2] == NULL || reopttree->reoptnodes[dummy2]->nvars == 0);

   SCIPdebugMessage("split the root into two dummy nodes.\n");
   SCIPdebugMessage(" -> store the node with identical bnd chgs at ID %d\n", dummy1);
   SCIPdebugMessage(" -> store the node with logic-or cons at ID %d\n", dummy2);

   /* dummy1:
    *   1. create the node
    *   2. add all bound changes
    *   3. convert all childids of the root to childids dummy1
    *   4. add the ID dummy1 as a child of the root node
    * */
   SCIP_CALL( createReoptnode(reopttree, dummy1) );
   reopttree->reoptnodes[dummy1]->parentID = 0;
   reopttree->reoptnodes[dummy1]->reopttype = SCIP_REOPTTYPE_TRANSIT;

   /* check memory */
   SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[dummy1], blkmem, nbndchgs, nchilds, 0) );
   assert(reopttree->reoptnodes[dummy1]->allocvarmem >= nbndchgs);
   assert(reopttree->reoptnodes[dummy1]->nvars == 0);
   assert(reopttree->reoptnodes[dummy1]->vars != NULL);
   assert(reopttree->reoptnodes[dummy1]->varbounds != NULL);
   assert(reopttree->reoptnodes[dummy1]->varboundtypes != NULL);

   /* copy bounds */
   for(v = 0; v < nbndchgs; v++)
   {
      reopttree->reoptnodes[dummy1]->vars[v] = reopttree->reoptnodes[0]->dualconscur->vars[v];
      reopttree->reoptnodes[dummy1]->varbounds[v] = reopttree->reoptnodes[0]->dualconscur->vals[v];
      reopttree->reoptnodes[dummy1]->varboundtypes[v] = SCIPsetIsFeasEQ(set, reopttree->reoptnodes[0]->dualconscur->vals[v], 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;
      reopttree->reoptnodes[dummy1]->nvars++;
   }

   /* move the children */
   SCIP_CALL( reoptMoveIDs(reopttree, blkmem, 0, dummy1) );
   assert(reopttree->reoptnodes[0]->nchilds == 0);

   /* add dummy1 as a child of the root node */
   SCIP_CALL( reoptAddChild(reopttree, 0, dummy1, blkmem) );

   /* dummy2:
    *   1. create the node
    *   2. add the constraint to ensure that at least one
    *      variable gets different
    *   3. add the ID dummy2 as a child of the root node
    * */
   SCIP_CALL( createReoptnode(reopttree, dummy2) );
   reopttree->reoptnodes[dummy2]->parentID = 0;
   reopttree->reoptnodes[dummy2]->reopttype = SCIP_REOPTTYPE_LOGICORNODE;

   /* create the constraint */
   SCIP_ALLOC( BMSallocMemory(&consdata) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &consdata->vars, nbndchgs) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &consdata->vals, nbndchgs) );
   consdata->allocmem = nbndchgs;
   consdata->nvars = nbndchgs;
   consdata->constype = REOPT_CONSTYPE_STRBRANCHED;

   for(v = 0; v < nbndchgs; v++)
   {
      consdata->vars[v] = reopttree->reoptnodes[0]->dualconscur->vars[v];
      consdata->vals[v] = reopttree->reoptnodes[0]->dualconscur->vals[v];
   }

   /* check memory for added constraints */
   SCIP_CALL( reoptnodeCheckMemory(reopttree->reoptnodes[dummy2], blkmem, 0, 0, 10) );

   /* add the constraint */
   reopttree->reoptnodes[dummy2]->conss[reopttree->reoptnodes[dummy2]->nconss] = consdata;
   reopttree->reoptnodes[dummy2]->nconss++;

   /* add dummy2 as a child of the root node */
   SCIP_CALL( reoptAddChild(reopttree, 0, dummy2, blkmem) );

   /* free the current dualconscur and assign dualconsnex */
   assert(reopttree->reoptnodes[0]->dualconscur->vars != NULL);
   assert(reopttree->reoptnodes[0]->dualconscur->vals != NULL);

   BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[0]->dualconscur->vals, reopttree->reoptnodes[0]->dualconscur->allocmem);
   BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[0]->dualconscur->vars, reopttree->reoptnodes[0]->dualconscur->allocmem);
   BMSfreeMemory(&reopttree->reoptnodes[0]->dualconscur);
   reopttree->reoptnodes[0]->dualconscur = NULL;

   if( reopttree->reoptnodes[0]->dualconsnex != NULL )
   {
      reopttree->reoptnodes[0]->dualconscur = reopttree->reoptnodes[0]->dualconsnex;
      reopttree->reoptnodes[0]->dualconsnex = NULL;
   }

   /* check if the flag dualfixing can be removed */
   reopttree->reoptnodes[0]->dualfixing = (reopttree->reoptnodes[0]->dualconscur != NULL);

   return SCIP_OKAY;
}

/* reset the stored information abound bound changes based on dual information */
void SCIPreoptResetDualcons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node,                    /**< node of the search tree */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   SCIP_REOPTTREE* reopttree;
   int id;

   assert(reopt != NULL);
   assert(node != NULL);

   reopttree = reopt->reopttree;
   assert(reopttree != NULL);

   id = SCIPnodeGetReoptID(node);

   /* the node was node create by reoptimization */
   if( id == -1 )
      return;

   if( reopttree->reoptnodes[id] != NULL && reopttree->reoptnodes[id]->dualconscur != NULL )
   {
      SCIPdebugMessage("reset dual (1) information at ID %d\n", id);

      BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->dualconscur->vals, reopttree->reoptnodes[id]->dualconscur->allocmem);
      BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->dualconscur->vars, reopttree->reoptnodes[id]->dualconscur->allocmem);
      BMSfreeMemory(&reopttree->reoptnodes[id]->dualconscur);
      reopttree->reoptnodes[id]->dualconscur = NULL;
   }

   if( reopttree->reoptnodes[id] != NULL && reopttree->reoptnodes[id]->dualconsnex != NULL )
   {
      SCIPdebugMessage("reset dual (2) information at ID %d\n", id);

      BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->dualconsnex->vals, reopttree->reoptnodes[id]->dualconsnex->allocmem);
      BMSfreeBlockMemoryArray(blkmem, &reopttree->reoptnodes[id]->dualconsnex->vars, reopttree->reoptnodes[id]->dualconsnex->allocmem);
      BMSfreeMemory(&reopttree->reoptnodes[id]->dualconsnex);
      reopttree->reoptnodes[id]->dualconsnex = NULL;
   }

   reopt->reopttree->reoptnodes[id]->dualfixing = FALSE;

   return;
}

/* returns the number of bound changes based on primal information including bound
 * changes directly after the first bound change based on dual information at the node
 * stored at ID id */
int SCIPreoptnodeGetIdNVars(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   id                       /**< id of a stored node */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   return reopt->reopttree->reoptnodes[id]->nvars + reopt->reopttree->reoptnodes[id]->nafterdualvars;
}

/* returns the number of bound changes at the node
 * stored at ID id */
int SCIPreoptnodeGetNConss(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   id                       /**< id of a stored node */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   return reopt->reopttree->reoptnodes[id]->nconss;
}

/* return the branching path stored at ID id */
void SCIPreoptnodeGetPath(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   id,                      /**< id of the stored node */
   SCIP_VAR**            vars,                    /**< array for variables */
   SCIP_Real*            vals,                    /**< array for values */
   SCIP_BOUNDTYPE*       boundtypes,              /**< array for bound types */
   int                   mem,                     /**< allocated memory */
   int*                  nvars,                   /**< pointer to store the number of variables */
   int*                  nafterdualvars           /**< pointer to store the number of variables stored after
                                                    *  the first dual decision */
   )
{
   SCIP_REOPTTREE* reopttree;
   int v;
   int nvars2;
   int nafterdualvars2;

   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(0 <= id && id <= reopt->reopttree->allocmemnodes);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(boundtypes != NULL);

   reopttree = reopt->reopttree;

   (*nvars) = reopttree->reoptnodes[id]->nvars;
   (*nafterdualvars) = reopttree->reoptnodes[id]->nafterdualvars;

   if( mem == 0 || mem < *nvars + *nafterdualvars )
      return;

   for(v = 0; v < *nvars; v++)
   {
      vars[v] = reopttree->reoptnodes[id]->vars[v];
      vals[v] = reopttree->reoptnodes[id]->varbounds[v];
      boundtypes[v] = reopttree->reoptnodes[id]->varboundtypes[v];
   }

   for(; v < *nvars + *nafterdualvars; v++)
   {
      vars[v] = reopttree->reoptnodes[id]->afterdualvars[v-(*nvars)];
      vals[v] = reopttree->reoptnodes[id]->afterdualvarbounds[v-(*nvars)];
      boundtypes[v] = reopttree->reoptnodes[id]->afterdualvarboundtypes[v-(*nvars)];
   }

   if( reopttree->reoptnodes[id]->parentID != 0 )
   {
      SCIPreoptnodeGetPath(reopt, reopttree->reoptnodes[id]->parentID, &vars[v],
            &vals[v], &boundtypes[v], mem, &nvars2, &nafterdualvars2);

      (*nvars) += nvars2;
      (*nafterdualvars) += nafterdualvars2;
   }

   return;
}

/* returns all added constraints at ID id */
void SCIPreoptnodeGetConss(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   id,                      /**< id of a stored node */
   SCIP_VAR***           vars,                    /**< 2-dim array of variables */
   SCIP_Real**           vals,                    /**< 2-dim array of values */
   int                   mem,                     /**< allocated memory for constraints */
   int*                  nconss,                  /**< pointer to store the number of constraints */
   int*                  nvars                    /**< pointer to store the number of variables */
   )
{
   int c;

   assert(reopt != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(nvars != NULL);
   assert(reopt->reopttree->reoptnodes[id] != NULL);


   (*nconss) = reopt->reopttree->reoptnodes[id]->nconss;

   if( mem < *nconss )
      return;

   for(c = 0; c < *nconss; c++)
   {
      assert(vars[c] != NULL);
      assert(vals[c] != NULL);

      vars[c] = reopt->reopttree->reoptnodes[id]->conss[c]->vars;
      vals[c] = reopt->reopttree->reoptnodes[id]->conss[c]->vals;
      nvars[c] = reopt->reopttree->reoptnodes[id]->conss[c]->nvars;
   }

   return;
}

/* replace the node stored at ID id by its child nodes */
SCIP_RETCODE SCIPreoptShrinkNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   id                       /**< id of a stored node */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(reopt->reopttree->reoptnodes[id] != NULL);
   assert(blkmem != NULL);

   SCIPdebugMessage(" -> shrink node at ID %d, replaced by %d child nodes.\n", id, reopt->reopttree->reoptnodes[id]->nchilds);

   /* move all children to the parent node */
   SCIP_CALL( moveChildrenUp(reopt, blkmem, id, reopt->reopttree->reoptnodes[id]->parentID) );

   /* delete the node */
   SCIP_CALL( reopttreeDeleteNode(reopt->reopttree , blkmem, id, TRUE) );

   /* add the ID the list of open IDs */
   SCIP_CALL( SCIPqueueInsert(reopt->reopttree->openids, (void*) (size_t) id) );

   return SCIP_OKAY;
}

/* delete a node stored in the reopttree */
SCIP_RETCODE SCIPreopttreeDeleteNode(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   id,                      /**< id of a stored node */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(reopt->reopttree->reoptnodes[id] != NULL);
   assert(blkmem != NULL);

   SCIP_CALL( reopttreeDeleteNode(reopt->reopttree, blkmem, id, TRUE) );

   return SCIP_OKAY;
}

/* reoptimize the node stored at ID id */
SCIP_RETCODE SCIPreoptApply(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   SCIP_PROB*            transprob,               /**< transformed problem */
   SCIP_PROB*            origprob,                /**< original problem */
   SCIP_TREE*            tree,                    /**< branching tree */
   SCIP_LP*              lp,                      /**< current LP */
   SCIP_BRANCHCAND*      branchcand,              /**< branching candidated */
   SCIP_EVENTQUEUE*      eventqueue,              /**< event queue */
   SCIP_NODE*            node_fix,                /**< node curresponding to the fixed part */
   SCIP_NODE*            node_cons,               /**< node corresponding to the pruned part */
   int                   id,                      /**< id of a stored node which should be reoptimized */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
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
   assert(node_fix != NULL || node_cons != NULL);
   assert(blkmem != NULL);

   SCIPdebugMessage("reoptimizing node at ID %d:\n", id);

   /* change all bounds */
   SCIP_CALL( changeAncestorBranchings(reopt->reopttree, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue, node_fix, node_cons, id, blkmem) );

   /* add the constraint to node_cons */
   if( node_cons != NULL && reopt->reopttree->reoptnodes[id]->dualconscur != NULL )
   {
      SCIP_CALL( addSplitcons(reopt, scip, set, stat, blkmem, transprob, origprob, tree, lp, branchcand, eventqueue, node_cons, id) );
   }

   /* fix all bound changes based on dual information in node and
    * convert all these bound changes to 'normal' bound changes */
   if( node_fix != NULL && reopt->reopttree->reoptnodes[id]->dualconscur != NULL )
   {
      SCIP_CALL( fixBounds(reopt->reopttree, set, stat, transprob, origprob, tree, lp, branchcand, eventqueue, node_fix, id, blkmem) );
   }

   /* add all local constraints to both nodes */
   SCIP_CALL( addLocalConss(scip, reopt, set, stat, blkmem, node_fix, node_cons, id) );

   /* set estimates */
   if( !SCIPsetIsInfinity(set, REALABS(reopt->reopttree->reoptnodes[id]->lowerbound)) )
   {
      if( node_fix != NULL && SCIPsetIsRelGE(set, reopt->reopttree->reoptnodes[id]->lowerbound, SCIPnodeGetLowerbound(node_fix)))
         SCIPnodeSetEstimate(node_fix, set, reopt->reopttree->reoptnodes[id]->lowerbound);

      if( node_cons != NULL && SCIPsetIsRelGE(set, reopt->reopttree->reoptnodes[id]->lowerbound, SCIPnodeGetLowerbound(node_cons)))
         SCIPnodeSetEstimate(node_cons, set, reopt->reopttree->reoptnodes[id]->lowerbound);
   }

   return SCIP_OKAY;
}

/* returns the reopttype of a node stored at ID id */
SCIP_REOPTTYPE SCIPreoptnodeGetType(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   int                   id                       /**< id a stored node */
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(0 <= id && id < reopt->reopttree->allocmemnodes);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   return reopt->reopttree->reoptnodes[id]->reopttype;
}


/* returns the time needed to store the nodes */
SCIP_Real SCIPreoptGetSavingtime(
   SCIP_REOPT*           reopt
   )
{
   assert(reopt != NULL);

   return SCIPclockGetTime(reopt->savingtime);
}

/* store a global constraint that should be added at the beginning of the next iteration */
SCIP_RETCODE SCIPreoptAddGlbCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   LOGICORDATA*          consdata,                /**< logicor data */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   assert(reopt != NULL);
   assert(consdata != NULL);
   assert(blkmem != NULL);

   if( consdata->nvars > 0 )
   {
      int pos;

      /* check the memory */
      SCIP_CALL( checkMemGlbCons(reopt, blkmem, reopt->nglbconss + 1) );
      assert(reopt->allocmemglbconss >= reopt->nglbconss+1);

      pos = reopt->nglbconss;

      /* allocate memory */
      SCIP_ALLOC( BMSallocMemory(&reopt->glbconss[pos]) );
      reopt->glbconss[pos]->allocmem = consdata->nvars;
      reopt->glbconss[pos]->nvars = consdata->nvars;
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->glbconss[pos]->vars, &consdata->vars, consdata->nvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reopt->glbconss[pos]->vals, &consdata->vals, consdata->nvars) );

      reopt->nglbconss++;
   }

   return SCIP_OKAY;
}

/* add the stored constraints globally to the problem */
SCIP_RETCODE SCIPreoptApplyGlbConss(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem                   /**< block memory */
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
      int v;

      assert(reopt->glbconss[c]->nvars > 0);
      assert(reopt->consallocmem >= reopt->glbconss[c]->nvars);

      SCIPdebugMessage("-> add constraints with %d vars\n", reopt->glbconss[c]->nvars);

      for(v = 0; v < reopt->glbconss[c]->nvars; v++)
      {
         reopt->consvars[v] = SCIPvarGetTransVar(reopt->glbconss[c]->vars[v]);

         /* negate the variable if it was fixed to 1 */
         if( SCIPsetIsFeasEQ(set, reopt->glbconss[c]->vals[v], 1) )
         {
            SCIP_CALL( SCIPvarNegate(reopt->consvars[v], blkmem, set, stat, &reopt->consvars[v]) );
         }
      }

      /* create the logic-or constraint and add them to the problem */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, "glblogicor", reopt->glbconss[c]->nvars,
            reopt->consvars, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE) );

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* delete the global constraints data */
      SCIPfreeBlockMemoryArrayNull(scip, &reopt->glbconss[c]->vals, reopt->glbconss[c]->nvars);
      SCIPfreeBlockMemoryArrayNull(scip, &reopt->glbconss[c]->vars, reopt->glbconss[c]->nvars);
      reopt->glbconss[c]->nvars = 0;
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

/**
 * construct a constraint to cut off a given solution.
 * add this constraint globally at the beginning of the next iteration.
 */
SCIP_RETCODE SCIPreoptAddGlbSolCons(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_SOL*             sol,                     /**< solution */
   SCIP_VAR**            vars,                    /**< array of original variables */
   SCIP_SET*             set,                     /**< global SCIP settings */
   SCIP_STAT*            stat,                    /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,                  /**< block memory */
   int                   nvars                    /**< number of original variables */
   )
{
   int nglbconss;
   int v;

   assert(reopt != NULL);
   assert(sol != NULL);
   assert(vars != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(nvars >= 0);

   nglbconss = reopt->nglbconss;

   /* allocate memory */
   SCIP_CALL( checkMemGlbCons(reopt, blkmem, nglbconss+1) );

   SCIP_ALLOC( BMSallocMemory(&reopt->glbconss[nglbconss]) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vars, nvars) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &reopt->glbconss[nglbconss]->vals, nvars) );

   reopt->glbconss[nglbconss]->allocmem = nvars;
   reopt->glbconss[nglbconss]->nvars = 0;
   reopt->glbconss[nglbconss]->constype = REOPT_CONSTYPE_SEPASOLUTION;

   /* save all variables */
   for(v = 0; v < nvars; v++)
   {
      SCIP_Real constant;
      SCIP_Real scalar;

      constant = 0;
      scalar = 1;

      reopt->glbconss[nglbconss]->vars[v] = vars[v];
      reopt->glbconss[nglbconss]->vals[v] = SCIPsolGetVal(sol, set, stat, vars[v]);

      /* transform into the original space */
      SCIP_CALL( SCIPvarGetOrigvarSum(&reopt->glbconss[nglbconss]->vars[v], &scalar, &constant) );
      reopt->glbconss[nglbconss]->vals[v] = (reopt->glbconss[nglbconss]->vals[v] - constant) / scalar;

      assert(SCIPsetIsFeasEQ(set, reopt->glbconss[nglbconss]->vals[v], 0) || SCIPsetIsFeasEQ(set, reopt->glbconss[nglbconss]->vals[v], 1));

      reopt->glbconss[nglbconss]->nvars++;
   }

   /* increase the number of global constraints */
   reopt->nglbconss++;

   return SCIP_OKAY;
}

/*
 * check if the LP of the given node should be solved or not.
 */
SCIP_RETCODE SCIPreoptGetSolveLP(
   SCIP_REOPT*           reopt,                   /**< reoptimization data structure */
   SCIP_NODE*            node,                    /**< node of the current search tree */
   SCIP_Bool*            solvelp                  /**< pointer to store whether the LP of the node should be solved */
   )
{
   int id;

   assert(reopt != NULL);
   assert(node != NULL);

   /* get the ID */
   id = SCIPnodeGetReoptID(node);

   (*solvelp) = TRUE;

   if( id == 0 )
   {
      if( reopt->reopttree->reoptnodes[0]->nchilds > 0 )
      {
         if( reopt->simtolastobj >= reopt->objsimrootlp )
            (*solvelp) = FALSE;
      }
   }
   else
      switch (reopt->solvelp) {
      /* solve all LPs */
      case 0:
         if( SCIPnodeGetReopttype(node) < SCIP_REOPTTYPE_LEAF )
         {
            if( reopt->reopttree->reoptnodes[id]->nvars < reopt->solvelpdiff)
               (*solvelp) = FALSE;
         }
         break;

      default:
         if( reopt->reopttree->reoptnodes[id]->nchilds > 0 )
         {
            if( reopt->reopttree->reoptnodes[id]->nvars < reopt->solvelpdiff && (int) SCIPnodeGetReopttype(node) < reopt->solvelp )
               (*solvelp) = FALSE;
         }
         break;
   }

   assert(*solvelp || reopt->reopttree->reoptnodes[id]->nchilds > 0);

   return SCIP_OKAY;
}

/* returns the number of restarts */
int SCIPreoptGetNRestarts(
   SCIP_REOPT*           reopt                    /**< reoptimization data structure */
   )
{
   assert(reopt != NULL);

   return reopt->nrestarts;
}

SCIP_Real SCIPreoptGetNodeLb(
   SCIP_REOPT*           reopt,
   int                   id
   )
{
   assert(reopt != NULL);
   assert(reopt->reopttree != NULL);
   assert(reopt->reopttree->reoptnodes[id] != NULL);

   return reopt->reopttree->reoptnodes[id]->lowerbound;
}

/*
 * initialize an empty node
 */
void SCIPreoptnodeCreateClear(
   SCIP_REOPTNODE*       reoptnode                /**< node of the reopttree */
   )
{
   assert(reoptnode != NULL);

   reoptnode->conss = NULL;
   reoptnode->nconss = 0;
   reoptnode->allocmemconss = 0;
   reoptnode->lpistate = NULL;
   reoptnode->childids = NULL;
   reoptnode->allocchildmem = 0;
   reoptnode->nchilds = 0;
   reoptnode->nvars = 0;
   reoptnode->nafterdualvars = 0;
   reoptnode->parentID = -1;
   reoptnode->dualfixing = FALSE;
   reoptnode->reopttype = SCIP_REOPTTYPE_NONE;
   reoptnode->allocvarmem = 0;
   reoptnode->allocafterdualvarmem = 0;
   reoptnode->vars = NULL;
   reoptnode->varbounds = NULL;
   reoptnode->varboundtypes = NULL;
   reoptnode->afterdualvars = NULL;
   reoptnode->afterdualvarbounds = NULL;
   reoptnode->afterdualvarboundtypes = NULL;
   reoptnode->dualconscur = NULL;
   reoptnode->dualconsnex = NULL;
   reoptnode->lowerbound = -1e+20;
}

/*
 * clear the node
 */
SCIP_RETCODE SCIPreoptnodeClear(
   SCIP_REOPTNODE**      reoptnode,               /**< pointer of reoptnode */
   BMS_BLKMEM*           blkmem                   /**< block memory */
   )
{
   assert(reoptnode != NULL);

   SCIP_CALL( reoptnodeDelete(reoptnode, blkmem) );

   return SCIP_OKAY;
}

/* add a variable to a given reoptnode */
SCIP_RETCODE SCIPreoptnodeAddVar(
   SCIP_REOPTNODE*       reoptnode,               /**< node of the reopttree */
   SCIP_VAR*             var,                     /**< variable to add */
   SCIP_Real             val,                     /**< value of the variable */
   SCIP_BOUNDTYPE        boundtype,               /**< boundtype of the variable */
   BMS_BLKMEM*           blkmem                   /**< block memory */
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
   reoptnode->nvars++;

   return SCIP_OKAY;
}

/* add a constraint to a given reoptnode */
SCIP_RETCODE SCIPreoptnodeAddCons(
   SCIP_REOPTNODE*       reoptnode,               /**< node of the reopttree */
   SCIP_VAR**            consvars,                /**< variables which are part of the constraint */
   SCIP_Real*            consvals,                /**< values of the variables */
   int                   nvars,                   /**< number of variables */
   REOPT_CONSTYPE        constype,                /**< type of the constraint */
   BMS_BLKMEM*           blkmem                   /**< block memory */
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

      SCIP_CALL( SCIPreoptnodeAddVar(reoptnode, consvars[0], 1-consvals[0],
            1-consvals[0] == 1 ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER, blkmem ) );
   }
   else
   {
      nconss = reoptnode->nconss;

      SCIP_CALL( reoptnodeCheckMemory(reoptnode, blkmem, 0, 0, nconss+1) );

      /* create the constraint */
      SCIP_ALLOC( BMSallocMemory(&reoptnode->conss[nconss]) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptnode->conss[nconss]->vars, consvars, nvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &reoptnode->conss[nconss]->vals, consvals, nvars) );
      reoptnode->conss[nconss]->allocmem = nvars;
      reoptnode->conss[nconss]->nvars = nvars;
      reoptnode->conss[nconss]->constype = constype;
      reoptnode->nconss++;
   }
   return SCIP_OKAY;
}

/* returns the number of bound changes stored in the reopttree at ID id*/
int SCIPreoptnodeGetNVars(
   SCIP_REOPTNODE*       reoptnode                /**< node of the roepttree */
   )
{
   assert(reoptnode != NULL);

   return reoptnode->nvars;
}

/* set the parent id */
void SCIPreoptnodeSetParentID(
   SCIP_REOPTNODE*       reoptnode,               /**< node of the reopttree */
   int                   parentid                 /**< id of the parent node */
   )
{
   assert(reoptnode != NULL);
   assert(0 <= parentid);

   reoptnode->parentID = parentid;
}