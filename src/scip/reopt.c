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
 * @brief  methods for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#define SCIP_DEBUG
#include "scip/def.h"
#include "scip/scip.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/misc.h"
#include "scip/reopt.h"
#include "scip/primal.h"
#include "scip/prob.h"


/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that sols[pos] array can store at least num entries */
static
SCIP_RETCODE ensureSolsSize(
   SCIP_REOPT*           reopt,             /**< primal data */
   SCIP_SET*             set,               /**< global SCIP settings */
   int                   num,               /**< minimum number of entries to store */
   int                   run                /**< run for which the memory should checked */
)
{
   assert(run >= 0);
   assert(run <= reopt->runsize);

   if( num > reopt->solssize[run] )
   {
      int newsize;

      newsize = num;
      SCIP_ALLOC( BMSreallocMemoryArray(&reopt->sols[run], newsize) );
      reopt->solssize[run] = newsize;
   }
   assert(num <= reopt->solssize[run]);

   return SCIP_OKAY;
}

/** ensures, that sols array can store at least num entries */
static
SCIP_RETCODE ensureRunSize(
   SCIP_REOPT*           reopt,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
)
{
   if( num >= reopt->runsize )
   {
      int newsize;
      int s;

      newsize = 2*reopt->runsize;
      SCIP_ALLOC( BMSreallocMemoryArray(&reopt->sols, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&reopt->nsols, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&reopt->solssize, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&reopt->objs, newsize) );

      for(s = reopt->runsize; s < newsize; s++)
      {
         reopt->sols[s] = NULL;
         reopt->objs[s] = NULL;
         reopt->nsols[s] = 0;
         reopt->solssize[s] = 0;
      }

      reopt->runsize = newsize;
   }
   assert(num <= reopt->runsize);

   return SCIP_OKAY;
}

/*
 * local methods
 */

static
int soltreeNInducedtSols(
   SCIP_SOLNODE*         node
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

/* subroutine to calculate Hamming-Distance */
static
SCIP_Real soltreeHammingDistSub(
   SCIP*                 scip,
   SCIP_SET*             set,
   SCIP_STAT*            stat,
   SCIP_SOL*             sol,
   SCIP_SOLNODE*         node,
   SCIP_VAR**            vars,
   int                   id
)
{
   SCIP_Real hamdist;
   SCIP_Real objval;

   assert(scip != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(sol != NULL);
   assert(node != NULL);
   assert(vars != NULL);

   if( node->lchild == NULL && node->rchild == NULL )
      return 0;

   assert(0 <= id && id < SCIPgetNVars(scip));

   hamdist = 0.0;
   objval = SCIPsolGetVal(sol, set, stat, vars[id]);

   if( SCIPsetIsFeasEQ(set, objval, 0) )
   {
      if( node->lchild != NULL )
      {
         hamdist += (soltreeNInducedtSols(node->lchild)
                   + soltreeHammingDistSub(scip, set, stat, sol, node->lchild, vars, id+1));
      }
      if( node->rchild != NULL )
      {
         hamdist += soltreeHammingDistSub(scip, set, stat, sol, node->rchild, vars, id+1);
      }
   }
   else
   {
      if( node->rchild != NULL )
      {
         hamdist += (soltreeNInducedtSols(node->rchild)
                   + soltreeHammingDistSub(scip, set, stat, sol, node->rchild, vars, id+1));
      }

      if( node->lchild != NULL )
      {
         hamdist += soltreeHammingDistSub(scip, set, stat, sol, node->lchild, vars, id+1);
      }
   }

   return hamdist;
}

/* return the average Hamming-Distance of a given solution
 * to all saved solutions */
static
SCIP_Real soltreeGetHammingDist(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   SCIP_SET*             set,
   SCIP_STAT*            stat,
   SCIP_SOL*             sol
)
{
   SCIP_VAR** vars;
   SCIP_Real hamdist;
   SCIP_Real objval;

   assert(scip != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(sol != NULL);

   hamdist = 0.0;

   /** there are no solutions */
   if( reopt->soltree->root->lchild == NULL && reopt->soltree->root->rchild == NULL )
      return 1.0;

   vars = SCIPgetVars(scip);
   objval = SCIPsolGetVal(sol, set, stat, vars[0]);

   if( SCIPsetIsFeasEQ(set, objval, 0) )
   {
      if( reopt->soltree->root->lchild != NULL )
      {
         hamdist += (soltreeNInducedtSols(reopt->soltree->root->lchild)
                   + soltreeHammingDistSub(scip, set, stat, sol, reopt->soltree->root->lchild, vars, 1));
      }

      if( reopt->soltree->root->rchild != NULL )
      {
         hamdist += soltreeHammingDistSub(scip, set, stat, sol, reopt->soltree->root->rchild, vars, 1);
      }
   }
   else
   {
      if( reopt->soltree->root->rchild != NULL )
      {
         hamdist += (soltreeNInducedtSols(reopt->soltree->root->rchild)
                   + soltreeHammingDistSub(scip, set, stat, sol, reopt->soltree->root->rchild, vars, 1));
      }

      if( reopt->soltree->root->lchild != NULL )
      {
         hamdist += soltreeHammingDistSub(scip, set, stat, sol, reopt->soltree->root->lchild, vars, 1);
      }
   }

   return hamdist/(reopt->soltree->nsols*SCIPgetNBinVars(scip));
}

/* add solutions to origprimal space */
static
SCIP_RETCODE soltreeUpdateSols(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   SCIP_Bool             sepabestsol,
   int                   run,
   int*                  naddedsols
)
{
   int s;
   SCIP_HEUR* heur;

   (*naddedsols) = 0;

   heur = SCIPfindHeur(scip, "reoptsols");
   assert(heur != NULL);

   for(s = 0; s < reopt->nsols[run]; s++)
   {
      SCIP_SOL* sol;
      SCIP_Real solobj;

      /* add a solution one once, if it occurs in more
       * than one run */
      if( !reopt->sols[run][s]->updated && !reopt->sols[run][s]->infeasible )
      {
         sol = reopt->sols[run][s]->sol;
         SCIPsolSetNodenum(sol, 0);
         SCIPsolSetHeur(sol, heur);
         SCIP_CALL( SCIPrecomputeSol(scip, sol) );
         solobj = SCIPgetSolTransObj(scip, sol);

         /* we do not want to add solutions with objective value +infinity */
         if( !SCIPisInfinity(scip, solobj) && !SCIPisInfinity(scip, -solobj)
         && SCIPisFeasLT(scip, solobj, SCIPgetCutoffbound(scip)) )
         {
            SCIP_SOL* bestsol = SCIPgetBestSol(scip);
            SCIP_Bool stored;
            SCIP_Bool feasible;

            if( sepabestsol )
            {
               SCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, FALSE, FALSE) );
            }
            else
               feasible = TRUE;

            if( feasible)
            {
               /* add primal solution to solution storage by copying it */
               SCIP_CALL( SCIPaddSol(scip, sol, &stored) );

               if( stored )
               {
                  if( bestsol != SCIPgetBestSol(scip) )
                     SCIPstoreSolutionGap(scip);

                  /* mark the solution as already added and used */
                  reopt->sols[run][s]->updated = TRUE;
                  reopt->sols[run][s]->used = TRUE;

                  (*naddedsols)++;
               }
            }
            else
               reopt->sols[run][s]->infeasible = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}

/* returns similariry of two objective functions */
static
SCIP_Real reoptSimilarity(
   SCIP_REOPT*           reopt,
   int                   nvars,
   int                   obj1_id,
   int                   obj2_id
)
{
   SCIP_Real norm_n1;
   SCIP_Real norm_n2;
   SCIP_Real scalarproduct;
   int id;

   /* calc similarity */
   norm_n1 = 0.0;
   norm_n2 = 0.0;
   scalarproduct = 0.0;
   for(id = 0; id < nvars; id++)
   {
      SCIP_Real c1;
      SCIP_Real c2;

      c1 = reopt->objs[obj1_id][id];
      c2 = reopt->objs[obj2_id][id];

      /** vector product */
      scalarproduct += c1*c2;

      /** norm of normalvector to obj1 */
      norm_n1 += c1*c1;

      /** norm of normalvector to obj1 */
      norm_n2 += c2*c2;
   }

   norm_n1 = sqrt(norm_n1);
   norm_n2 = sqrt(norm_n2);

   if( norm_n1*norm_n2 == 0 )
      return -2.0;
   else
      return scalarproduct/(norm_n1*norm_n2);
}

static
SCIP_RETCODE createSolTree(
   SCIP_REOPT*           reopt
)
{
   assert(reopt != NULL);

   SCIP_ALLOC( BMSallocMemory(&reopt->soltree) );
   reopt->soltree->nsols = 0;

   SCIP_ALLOC( BMSallocMemory(&reopt->soltree->root) );
   reopt->soltree->root->sol = NULL;
   reopt->soltree->root->updated = FALSE;
   reopt->soltree->root->father = NULL;
   reopt->soltree->root->rchild = NULL;
   reopt->soltree->root->lchild = NULL;
   reopt->soltree->root->used = FALSE;
   reopt->soltree->root->infeasible = FALSE;
   reopt->soltree->root->val = -1;

   return SCIP_OKAY;
}

static
SCIP_RETCODE soltreefreeNode(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   SCIP_SOLNODE*         node
)
{
   assert(reopt != NULL);
   assert(node != NULL);

   /* free recursive right subtree */
   if( node->rchild != NULL )
   {
      SCIP_CALL( soltreefreeNode(scip, reopt, node->rchild) );
   }

   /* free recursive left subtree */
   if( node->lchild != NULL )
   {
      SCIP_CALL( soltreefreeNode(scip, reopt, node->lchild) );
   }

   if( node->sol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &node->sol) );
   }

   /* free this nodes */
   BMSfreeMemory(&node);

   return SCIP_OKAY;
}

/* free solution tree */
static
SCIP_RETCODE freeSolTree(
   SCIP*                 scip,
   SCIP_REOPT*           reopt
)
{
   assert(reopt != NULL);
   assert(reopt->soltree != NULL);
   assert(reopt->soltree->root != NULL);

   /* free all nodes recursive */
   SCIP_CALL( soltreefreeNode(scip, reopt, reopt->soltree->root) );

   BMSfreeMemory(&reopt->soltree);

   return SCIP_OKAY;
}

/* add a node to the solution tree */
static
SCIP_RETCODE soltreeAddNode(
   SCIP_REOPT*           reopt,
   SCIP_SOLNODE*         father,
   SCIP_Bool             rchild,
   SCIP_Bool             lchild,
   SCIP_Real             val
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
   newnode->used = FALSE;
   newnode->infeasible = FALSE;
   newnode->father = father;
   newnode->rchild = NULL;
   newnode->lchild = NULL;
   newnode->val = val;

   if( rchild )
      father->rchild = newnode;
   else
      father->lchild = newnode;

   return SCIP_OKAY;
}

/* add a solution */
static
SCIP_RETCODE soltreeAddSol(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   SCIP_SET*             set,
   SCIP_STAT*            stat,
   SCIP_VAR**            vars,
   SCIP_SOL*             sol,
   SCIP_SOLNODE**        solnode,
   int                   nvars,
   SCIP_Bool             bestsol,
   SCIP_Bool*            added
)
{
   SCIP_SOLNODE* cursolnode;
   SCIP_Real hamdist;
   int varid;

   assert(reopt != NULL);
   assert(sol != NULL);

   cursolnode = reopt->soltree->root;
   (*added) = FALSE;
   hamdist = 0.0;

   /** either all variables are delete or the given minimal Hamming-Distance is at most
    *  1/SCIPgetNVars(scip), which is the minimum possible Hamming-Distance.
    */
   if( SCIPgetNVars(scip) == 0 || set->reopt_minavghamdist == 1 || set->reopt_minavghamdist < (SCIP_Real)1/(SCIPgetNBinVars(scip)*reopt->soltree->nsols) )
   {
      hamdist = 1;
   }
   else
   {
      hamdist = soltreeGetHammingDist(scip, reopt, set, stat, sol);
   }

   /* add the solution iff the solution differs in at least one variable
    * to all saved solution and the avarage Hamming-Distance is greater or
    * equal to reopt_minavghamdist */
   for(varid = 0; varid < nvars && hamdist >= set->reopt_minavghamdist; varid++)
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
               SCIP_CALL( soltreeAddNode(reopt, cursolnode, TRUE, FALSE, objval) );
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
               SCIP_CALL( soltreeAddNode(reopt, cursolnode, FALSE, TRUE, objval) );
               assert(cursolnode->lchild != NULL);
               (*added) = TRUE;
            }
            cursolnode = cursolnode->lchild;
         }
      }
   }

   /* the solution was added */
   if( *added || hamdist >= set->reopt_minavghamdist )
   {
      SCIP_SOL* copysol;

      assert(cursolnode->lchild == NULL && cursolnode->rchild == NULL);

      if( *added )
      {
         SCIP_CALL( SCIPcreateSolCopy(scip, &copysol, sol) );
         SCIPsolSetHeur(sol, NULL);
         SCIPsolSetNodenum(sol, 0);
         cursolnode->sol = copysol;
         reopt->soltree->nsols++;
      }
      else
         /* this is a pseudo add; we do not want to save this solution
          * more than once, but we will link this solution to the solution
          * storage of this round */
         (*added) = TRUE;

      if( bestsol )
         reopt->lastbestsol = cursolnode->sol;

      (*solnode) = cursolnode;
   }

#ifdef SCIP_DEBUG
   {
      if( set->reopt_minavghamdist <= 0.5 && set->reopt_minavghamdist > 0)
         SCIPdebugMessage("** reoptimization ** solution%s added (Hamming-Distance %.4f).\n", (*added) ? "" : " not", hamdist);
   }
#endif

   return SCIP_OKAY;
}

/* set all marks updated to FALSE */
static
void soltreeResetMarks(
   SCIP_SOLNODE*         node
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
   SCIP_SOLNODE*         node
)
{
   int nusedsols;

   assert(node != NULL);

   nusedsols = 0;

   if(node->lchild != NULL)
      nusedsols += soltreeGetNUsedSols(node->lchild);
   if(node->rchild != NULL)
      nusedsols += soltreeGetNUsedSols(node->rchild);
   if(node->rchild == NULL && node->lchild == NULL && node->used)
      nusedsols = 1;

   return nusedsols;
}

/*
 * public methods
 */

/** creates reopt data */
SCIP_RETCODE SCIPreoptCreate(
   SCIP_REOPT**          reopt                   /**< pointer to reopt data */
)
{
   int s;

   assert(reopt != NULL);

   SCIP_ALLOC( BMSallocMemory(reopt) );
   (*reopt)->sols = NULL;
   (*reopt)->nsols = NULL;
   (*reopt)->solssize = NULL;
   (*reopt)->runsize = 200;
   (*reopt)->run = -1;
   (*reopt)->simtolastobj = -2.0;
   (*reopt)->simtofirstobj = -2.0;
   (*reopt)->lastbestsol = NULL;
   (*reopt)->firstobj = -1;

   SCIP_ALLOC( BMSallocMemoryArray(&(*reopt)->sols, (*reopt)->runsize) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*reopt)->nsols, (*reopt)->runsize) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*reopt)->solssize, (*reopt)->runsize) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*reopt)->objs, (*reopt)->runsize) );

   for(s = 0; s < (*reopt)->runsize; s++)
   {
      (*reopt)->nsols[s] = 0;
      (*reopt)->solssize[s] = 0;
      (*reopt)->sols[s] = NULL;
      (*reopt)->objs[s] = NULL;
   }

   /* create SCIP_SOLTREE */
   SCIP_CALL( createSolTree((*reopt)) );

   return SCIP_OKAY;
}

/** frees reopt data */
SCIP_RETCODE SCIPreoptFree(
   SCIP*                 scip,
   SCIP_REOPT**          reopt,              /**< pointer to primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int p;

   assert(reopt != NULL);
   assert(*reopt != NULL);

   /* free solution tree */
   SCIP_CALL( freeSolTree(scip, (*reopt)) );

   /* free solutions */
   for( p = (*reopt)->runsize-1; p >= 0; --p )
   {
      if( (*reopt)->sols[p] != NULL )
      {
         BMSfreeMemoryArray(&(*reopt)->sols[p]);
      }

      if( (*reopt)->objs[p] != NULL )
      {
         BMSfreeMemoryArray(&(*reopt)->objs[p]);
      }
   }

   BMSfreeMemoryArray(&(*reopt)->sols);
   BMSfreeMemoryArray(&(*reopt)->nsols);
   BMSfreeMemoryArray(&(*reopt)->solssize);
   BMSfreeMemoryArray(&(*reopt)->objs);
   BMSfreeMemory(reopt);

   return SCIP_OKAY;
}

/** add a solution to the last run */
SCIP_RETCODE SCIPreoptAddSol(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   SCIP_SET*             set,
   SCIP_STAT*            stat,
   SCIP_SOL*             sol,
   SCIP_Bool             bestsol,
   SCIP_Bool*            added,
   int                   run
)
{
   SCIP_SOLNODE* solnode;
   int insertpos;

   assert(reopt != NULL);
   assert(set != NULL);
   assert(sol != NULL);
   assert(run >= 0);

   assert(reopt->sols[run] != NULL);

   /* check memory */
   SCIP_CALL( ensureSolsSize(reopt, set, reopt->nsols[run], run) );

   solnode = NULL;

   /** ad solution to solution tree */
   SCIP_CALL( soltreeAddSol(scip, reopt, set, stat, SCIPgetOrigVars(scip), sol, &solnode, SCIPgetNOrigVars(scip), bestsol, added) );

   if( (*added) )
   {
      assert(solnode != NULL);

      /** add solution */
      insertpos = reopt->nsols[run];
      reopt->sols[run][insertpos] = solnode;
      reopt->nsols[run]++;
      assert(reopt->nsols[run] <= set->reopt_savesols);
   }

   return SCIP_OKAY;
}

/* add a run */
SCIP_RETCODE SCIPreoptAddRun(
   SCIP_SET*             set,
   SCIP_REOPT*           reopt,
   int                   run,
   int                   size
)
{
   assert(reopt != NULL);

   /* check memory */
   SCIP_CALL( ensureRunSize(reopt, set, run) );

   /* set number of last run */
   reopt->run = run;

   /* allocate memory */
   reopt->solssize[run] = size;
   SCIP_ALLOC( BMSallocMemoryArray(&reopt->sols[run], size) );

   return SCIP_OKAY;
}

/* returns number of solution */
int SCIPreoptGetNSols(
   SCIP_REOPT*           reopt,
   int                   run
)
{
   assert(reopt != NULL);
   assert(run <= reopt->runsize);

   if( reopt->sols[run] == NULL )
      return 0;
   else
      return reopt->nsols[run];
}

/* add solutions to origprimal space */
SCIP_RETCODE SCIPreoptUpdateSols(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   int                   savesols,
   SCIP_Real             objsim,
   SCIP_Real             sepabestsol
)
{
   int naddedsols;
   int naddedsolsrun;
   int nsols;
   int run;

   assert(reopt != NULL);

   /** we can skip this method if we have no solutions saved */
   if( savesols == 0 )
      return SCIP_OKAY;

   naddedsols = 0;
   naddedsolsrun = -1;
   nsols = reopt->soltree->nsols;

   for(run = reopt->run; run >= 0; run--)
   {
      SCIP_Real sim;

      naddedsolsrun = 0;

      if( objsim == 0 )
         sim = 1;
      else
         sim = reoptSimilarity(reopt, SCIPgetNOrigVars(scip), run, reopt->run+1);

      if( sim >= objsim )
      {
         SCIP_CALL( soltreeUpdateSols(scip, reopt, sepabestsol, run, &naddedsolsrun) );

         naddedsols += naddedsolsrun;

#ifdef SCIP_DEBUG
         {
            if( naddedsolsrun > 0 )
            {
               if( objsim == 0  )
                  SCIPdebugMessage("** reoptimization ** add %d solutions from run %d.\n", naddedsolsrun, run);
               else
                  SCIPdebugMessage("** reoptimization ** add %d solutions from run %d (lambda %.4f).\n", naddedsolsrun, run, sim);
            }
         }
#endif
      }
   }

   printf("%u/%u reoptimized solutions, feasible solution found by reoptsols heuristic, best objective value %.6e\n", naddedsols, nsols, SCIPgetPrimalbound(scip));

   /* reset the marks for added solutions */
   if(  reopt->soltree ->nsols > 0 )
       soltreeResetMarks(reopt->soltree->root);

    return SCIP_OKAY;
 }

 /* returns the number of saved solutions overall runs */
 int SCIPreoptNSavedSols(
    SCIP_REOPT*           reopt
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
   SCIP_REOPT*           reopt
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
SCIP_RETCODE SCIPreoptSaveObj(
   SCIP*                 scip,
   SCIP_REOPT*           reopt,
   SCIP_SET*             set,
   int                   run
)
{
   SCIP_VAR** vars;
   int id;

   assert(reopt != NULL);

   /* check memory */
   SCIP_CALL( ensureRunSize(reopt, set, run) );

   vars = SCIPgetOrigVars(scip);

   /* get memory */
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &reopt->objs[run], SCIPgetNOrigVars(scip)) );

   /* save coefficients */
   vars = SCIPgetVars(scip);
   for(id = 0; id < SCIPgetNOrigVars(scip); id++)
   {
      reopt->objs[run][id] = SCIPvarGetObj(vars[id]);

      /* mark this objective as the first non empty */
      if( reopt->firstobj == -1 && reopt->objs[run][id] != 0 )
         reopt->firstobj = run;
   }

   /* calculate similarity to last objective */
   if( run >= 1 )
   {
      /* calculate similarity to first objective */
      if( run > 1 && reopt->firstobj < run )
         reopt->simtofirstobj = reoptSimilarity(reopt, SCIPgetNOrigVars(scip), run, reopt->firstobj);

      /* calculate similarity to last objective */
      reopt->simtolastobj = reoptSimilarity(reopt, SCIPgetNOrigVars(scip), run, run-1);

      SCIPdebugMessage("new objective has similarity of %.4f/%.4f compared to first/previous.\n", reopt->simtofirstobj, reopt->simtolastobj);
   }

   SCIPdebugMessage("saved obj for run %d.\n", run);

   return SCIP_OKAY;
}

/* check if the current and the previous objective are similar enough
 * returns TRUE if we want to restart, otherwise FALSE */
SCIP_Bool SCIPreoptCheckRestart(
   SCIP_REOPT*           reopt,
   SCIP_SET*             set,
   int                   nvars,
   SCIP_Real*            sim
)
{
   assert(reopt != NULL);
   assert(set != NULL);

   (*sim) = 0.0;

   if( reopt->run > 0 && set->reopt_delay > 0.0 )
   {
      (*sim) = reopt->simtolastobj;
   }

   if( (*sim) >= set->reopt_delay )
      return FALSE;
   else
      return TRUE;
}

/*
 * returns the similarity to the previous objective function,
 * if no objective functions are saved the similarity is -2.0.
 */
SCIP_Real SCIPreoptGetSimToPrevious(
      SCIP_REOPT*        reopt
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
      SCIP_REOPT*        reopt
)
{
   assert(reopt != NULL);
   return reopt->simtofirstobj;
}

/*
 * returns the best solution of the last run
 */
SCIP_SOL* SCIPreoptGetLastBestSol(
   SCIP_REOPT*           reopt
)
{
   assert(reopt != NULL);
   return reopt->lastbestsol;
}

/*
 * returns true if the sign of the objective coefficient has changed
 */
SCIP_Bool SCIPreoptIsObjCoefNegated(
   SCIP_REOPT*           reopt,
   int                   varnr
)
{
   assert(reopt != NULL);

   if( reopt->objs[reopt->run][varnr] == 0 || reopt->objs[reopt->run-1][varnr] == 0 )
      return FALSE;

   if( reopt->objs[reopt->run] != NULL && reopt->run >= 1 )
   {
      return ( (SCIP_Real)reopt->objs[reopt->run][varnr]/reopt->objs[reopt->run-1][varnr] < 0 ? TRUE : FALSE );
   }

   return FALSE;
}

/*
 * returns if a variable enters or leaves the objective function.
 */
void SCIPreoptEnterOrLeaveObj(
   SCIP_REOPT*           reopt,
   int                   varnr,
   SCIP_Bool*            entering,
   SCIP_Bool*            leaving
)
{
   assert(reopt != NULL);

   *entering = FALSE;
   *leaving = FALSE;

   if( reopt->objs[reopt->run] == NULL || reopt->run == 0)
      return;

   /* variable leaves the objective */
   if( reopt->objs[reopt->run][varnr] == 0 && reopt->objs[reopt->run-1][varnr] != 0 )
   {
      *leaving = TRUE;
   }
   /* variable enters the objective function */
   else if( reopt->objs[reopt->run][varnr] != 0 && reopt->objs[reopt->run-1][varnr] == 0 )
   {
      *entering = TRUE;
   }

   return;
}
