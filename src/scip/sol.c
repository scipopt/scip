/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sol.c
 * @brief  datastructures and methods for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <limits.h>

#include "sol.h"



/** primal CIP solution
 *  For reasons of efficiency, a working solution only stores values that have been accessed at least once,
 *  or that have been changed from the value in the solution's source.
 *  The user has to call SCIPsolUnlink() in order to retrieve all non-cached elements from the solution's source
 *  and to store the values in the solution's own array. This changes the solution's origin to SCIP_SOLORIGIN_ZERO.
 *  A linked solution with origin SCIP_SOLORIGIN_LPSOL or SCIP_SOLORIGIN_PSEUDOSOL becomes invalid after the
 *  next node is activated (i.e. the LP and pseudo solutions changed) and cannot be accessed anymore.
 */
struct Sol
{
   REALARRAY*       vals;               /**< solution values for variables */
   BOOLARRAY*       valid;              /**< for solutions originating from LPSOL or PSEUDOSOL: TRUE iff variable's val
                                         *   is valid; otherwise the value has to be retrieved from the origin */
   HEUR*            heur;               /**< heuristic that found the solution (or NULL if it's an LP solution) */
   Real             obj;                /**< objective value of solution */
   Longint          nodenum;            /**< last node number, where this solution was modified */
   unsigned int     solorigin:2;        /**< origin of solution: where to retrieve uncached elements */
};




/** creates primal CIP solution, initialized to zero */
RETCODE SCIPsolCreate(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, sol) );   
   CHECK_OKAY( SCIPrealarrayCreate(&(*sol)->vals, memhdr) );
   (*sol)->valid = NULL;
   (*sol)->heur = heur;
   (*sol)->obj = 0.0;
   (*sol)->nodenum = stat->nnodes;
   (*sol)->solorigin = SCIP_SOLORIGIN_ZERO;

   debugMessage("created empty solution %p\n", *sol);

   return SCIP_OKAY;
}

/** creates a copy of a primal CIP solution */
RETCODE SCIPsolCopy(
   SOL**            sol,                /**< pointer to store the copy of the primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SOL*             sourcesol           /**< primal CIP solution to copy */
   )
{
   assert(sol != NULL);
   assert(sourcesol != NULL);
   assert((sourcesol->solorigin == SCIP_SOLORIGIN_ZERO) ^ (sourcesol->valid != NULL));

   debugMessage("copying solution %p\n", sourcesol);

   ALLOC_OKAY( allocBlockMemory(memhdr, sol) );   
   CHECK_OKAY( SCIPrealarrayCopy(&(*sol)->vals, memhdr, sourcesol->vals) );
   if( sourcesol->solorigin == SCIP_SOLORIGIN_ZERO )
      (*sol)->valid = NULL;
   else
   {
      CHECK_OKAY( SCIPboolarrayCopy(&(*sol)->valid, memhdr, sourcesol->valid) );
   }
   (*sol)->heur = sourcesol->heur;
   (*sol)->obj = sourcesol->obj;
   (*sol)->nodenum = sourcesol->nodenum;
   (*sol)->solorigin = sourcesol->solorigin;

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the actual LP solution */
RETCODE SCIPsolCreateLPSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp,                 /**< actual LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   debugMessage("creating solution from LP\n");

   CHECK_OKAY( SCIPsolCreate(sol, memhdr, stat, heur) );
   CHECK_OKAY( SCIPsolLinkLPSol(*sol, memhdr, stat, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the actual pseudo solution */
RETCODE SCIPsolCreatePseudoSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(tree != NULL);

   debugMessage("creating solution from pseudo solution\n");

   CHECK_OKAY( SCIPsolCreate(sol, memhdr, stat, heur) );
   CHECK_OKAY( SCIPsolLinkPseudoSol(*sol, memhdr, set, stat, tree) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the actual solution */
RETCODE SCIPsolCreateActSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(tree != NULL);

   debugMessage("creating solution from actual solution\n");

   if( tree->actnodehaslp )
   {
      CHECK_OKAY( SCIPsolCreateLPSol(sol, memhdr, stat, lp, heur) );
   }
   else
   {
      CHECK_OKAY( SCIPsolCreatePseudoSol(sol, memhdr, set, stat, tree, heur) );
   }

   return SCIP_OKAY;
}

/** frees primal CIP solution */
RETCODE SCIPsolFree(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(sol != NULL);
   assert(*sol != NULL);
   assert(((*sol)->solorigin == SCIP_SOLORIGIN_ZERO) ^ ((*sol)->valid != NULL));

   CHECK_OKAY( SCIPrealarrayFree(&(*sol)->vals) );
   if( (*sol)->solorigin != SCIP_SOLORIGIN_ZERO )
   {
      CHECK_OKAY( SCIPboolarrayFree(&(*sol)->valid) );
   }
   freeBlockMemory(memhdr, sol);

   return SCIP_OKAY;
}

/** copies actual LP solution into CIP solution by linking */
RETCODE SCIPsolLinkLPSol(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   debugMessage("linking solution to LP\n");

   /* clear the old solution */
   SCIPrealarrayClear(sol->vals);
   if( sol->solorigin == SCIP_SOLORIGIN_ZERO )
   {
      assert(sol->valid == NULL);
      CHECK_OKAY( SCIPboolarrayCreate(&sol->valid, memhdr) );
   }
   else
   {
      assert(sol->valid != NULL);
      SCIPboolarrayClear(sol->valid);
   }

   /* link solution to LP solution */
   sol->obj = SCIPlpGetObjval(lp);
   sol->solorigin = SCIP_SOLORIGIN_LPSOL;
   sol->nodenum = stat->nnodes;
   
   debugMessage(" -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies actual pseudo solution into CIP solution by linking */
RETCODE SCIPsolLinkPseudoSol(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);

   debugMessage("linking solution to pseudo solution\n");

   /* clear the old solution */
   SCIPrealarrayClear(sol->vals);
   if( sol->solorigin == SCIP_SOLORIGIN_ZERO )
   {
      assert(sol->valid == NULL);
      CHECK_OKAY( SCIPboolarrayCreate(&sol->valid, memhdr) );
   }
   else
   {
      assert(sol->valid != NULL);
      SCIPboolarrayClear(sol->valid);
   }

   /* link solution to pseudo solution */
   sol->obj = SCIPtreeGetActPseudoobjval(tree, set);
   sol->solorigin = SCIP_SOLORIGIN_PSEUDOSOL;
   sol->nodenum = stat->nnodes;

   debugMessage(" -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies actual solution (LP or pseudo solution) into CIP solution by linking */
RETCODE SCIPsolLinkActSol(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(tree != NULL);

   debugMessage("linking solution to actual solution\n");

   if( tree->actnodehaslp )
   {
      CHECK_OKAY( SCIPsolLinkLPSol(sol, memhdr, stat, lp) );
   }
   else
   {
      CHECK_OKAY( SCIPsolLinkPseudoSol(sol, memhdr, set, stat, tree) );
   }

   return SCIP_OKAY;
}

/** clears primal CIP solution */
RETCODE SCIPsolClear(
   SOL*             sol,                /**< primal CIP solution */
   STAT*            stat                /**< problem statistics data */
   )
{
   assert(sol != NULL);
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) ^ (sol->valid != NULL));

   SCIPrealarrayClear(sol->vals);
   sol->obj = 0.0;
   if( sol->solorigin != SCIP_SOLORIGIN_ZERO )
   {
      CHECK_OKAY( SCIPboolarrayFree(&sol->valid) );
      sol->solorigin = SCIP_SOLORIGIN_ZERO;
   }
   sol->nodenum = stat->nnodes;

   return SCIP_OKAY;
}

/** stores solution value of variable in solution's own array */
static
RETCODE solUnlinkVar(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var                 /**< problem variable */
   )
{
   assert(sol != NULL);
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) ^ (sol->valid != NULL));
   assert(var != NULL);
   assert(var->varstatus == SCIP_VARSTATUS_COLUMN || var->varstatus == SCIP_VARSTATUS_LOOSE);

   switch( sol->solorigin )
   {
   case SCIP_SOLORIGIN_ZERO:
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_LPSOL:
      /*debugMessage("completing variable <%s> in LP solution %p\n", var->name, sol);*/
      if( !SCIPboolarrayGetVal(sol->valid, var->index) )
      {
         assert(SCIPrealarrayGetVal(sol->vals, var->index) == 0.0);
         CHECK_OKAY( SCIPrealarraySetVal(sol->vals, set, var->index, SCIPvarGetLPSol(var)) );
         CHECK_OKAY( SCIPboolarraySetVal(sol->valid, set, var->index, TRUE) );
      }
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_PSEUDOSOL:
      /*debugMessage("completing variable <%s> in pseudo solution %p\n", var->name, sol);*/
      if( !SCIPboolarrayGetVal(sol->valid, var->index) )
      {
         assert(SCIPrealarrayGetVal(sol->vals, var->index) == 0.0);
         CHECK_OKAY( SCIPrealarraySetVal(sol->vals, set, var->index, SCIPvarGetPseudoSol(var)) );
         CHECK_OKAY( SCIPboolarraySetVal(sol->valid, set, var->index, TRUE) );
      }
      return SCIP_OKAY;

   default:
      errorMessage("unknown variable status");
      abort();
   }
}

/** stores solution values of variables in solution's own array */
RETCODE SCIPsolUnlink(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   )
{
   int v;

   assert(sol != NULL);
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) ^ (sol->valid != NULL));

   if( sol->solorigin != SCIP_SOLORIGIN_ZERO )
   {
      debugMessage("completing solution %p\n", sol);

      for( v = 0; v < prob->nvars; ++v )
      {
         assert(prob->vars[v] != NULL);
         CHECK_OKAY( solUnlinkVar(sol, set, prob->vars[v]) );
      }
      
      CHECK_OKAY( SCIPboolarrayFree(&sol->valid) );
      sol->solorigin = SCIP_SOLORIGIN_ZERO;
   }

   return SCIP_OKAY;
}

/** sets value of variable in primal CIP solution */
RETCODE SCIPsolSetVal(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   )
{
   Real oldval;

   assert(sol != NULL);
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) ^ (sol->valid != NULL));
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO || sol->nodenum == stat->nnodes);
   assert(var != NULL);

   debugMessage("setting value of <%s> in solution %p to %g\n", var->name, sol, val);

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      return SCIPsolSetVal(sol, set, stat, var->data.transvar, val);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      CHECK_OKAY( SCIPsolGetVal(sol, set, stat, var, &oldval) );
      if( !SCIPsetIsEQ(set, val, oldval) )
      {
         CHECK_OKAY( SCIPrealarraySetVal(sol->vals, set, var->index, val) );
         sol->obj += var->obj * (val - oldval);
         sol->nodenum = stat->nnodes;
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot set solution value for fixed variable");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      return SCIPsolSetVal(sol, set, stat, var->data.aggregate.var,
         (val - var->data.aggregate.constant)/var->data.aggregate.scalar);

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot set solution value for multiple aggregated variable");
      return SCIP_INVALIDDATA;

   default:
      errorMessage("unknown variable status");
      abort();
   }
}

/** increases value of variable in primal CIP solution */
RETCODE SCIPsolIncVal(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   VAR*             var,                /**< variable to increase solution value for */
   Real             incval              /**< increment for solution value of variable */
   )
{
   assert(sol != NULL);
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) ^ (sol->valid != NULL));
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO || sol->nodenum == stat->nnodes);
   assert(var != NULL);

   debugMessage("increasing value of <%s> in solution %p by %g\n", var->name, sol, incval);

   if( SCIPsetIsZero(set, incval) )
      return SCIP_OKAY;

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      return SCIPsolIncVal(sol, set, stat, var->data.transvar, incval);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      CHECK_OKAY( solUnlinkVar(sol, set, var) );
      CHECK_OKAY( SCIPrealarrayIncVal(sol->vals, set, var->index, incval) );
      sol->obj += var->obj * incval;
      sol->nodenum = stat->nnodes;
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot set solution value for fixed variable");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPsetIsZero(set, var->data.aggregate.scalar));
      return SCIPsolIncVal(sol, set, stat, var->data.aggregate.var,
         (incval - var->data.aggregate.constant)/var->data.aggregate.scalar);

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot set solution value for multiple aggregated variable");
      return SCIP_INVALIDDATA;

   default:
      errorMessage("unknown variable status");
      abort();
   }
}

/** returns value of variable in primal CIP solution */
RETCODE SCIPsolGetVal(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   VAR*             var,                /**< variable to get value for */
   Real*            solval              /**< pointer to store the solution value */
   )
{
   Real val;
   int i;

   assert(sol != NULL);
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) ^ (sol->valid != NULL));
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO || sol->nodenum == stat->nnodes);
   assert(var != NULL);
   assert(solval != NULL);

   /*debugMessage("getting value of <%s> in solution %p\n", var->name, sol);*/

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( var->varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      CHECK_OKAY( SCIPsolGetVal(sol, set, stat, var->data.transvar, solval) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      CHECK_OKAY( solUnlinkVar(sol, set, var) );
      *solval = SCIPrealarrayGetVal(sol->vals, var->index);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      assert(var->glbdom.lb == var->glbdom.ub);
      assert(var->actdom.lb == var->actdom.ub);
      assert(var->glbdom.lb == var->actdom.lb);
      *solval = var->glbdom.lb;
      return SCIP_OKAY;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      CHECK_OKAY( SCIPsolGetVal(sol, set, stat, var->data.aggregate.var, solval) );
      *solval *= var->data.aggregate.scalar;
      *solval += var->data.aggregate.constant;
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      *solval = var->data.multaggr.constant;
      for( i = 0; i < var->data.multaggr.nvars; ++i )
      {
         CHECK_OKAY( SCIPsolGetVal(sol, set, stat, var->data.multaggr.vars[i], &val) );
         *solval += var->data.multaggr.scalars[i] * val;
      }
      return SCIP_OKAY;

   default:
      errorMessage("unknown variable status");
      abort();
   }
}

/** checks primal CIP solution for feasibility */
RETCODE SCIPsolCheck(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   Bool*            feasible            /**< stores whether solution is feasible */
   )
{
   RESULT result;
   int h;

   assert(set != NULL);
   assert(feasible != NULL);

   *feasible = TRUE;
   for( h = 0; h < set->nconshdlrs && *feasible; ++h )
   {
      CHECK_OKAY( SCIPconshdlrCheck(set->conshdlrs[h], memhdr, set, prob, sol, chckintegrality, chcklprows, &result) );
      *feasible &= (result == SCIP_FEASIBLE);
   }

   return SCIP_OKAY;
}

/** gets objective value of primal CIP solution */
Real SCIPsolGetObj(
   SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->obj;
}

/** gets node number, where this solution was found */
Longint SCIPsolGetNodenum(
   SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->nodenum;
}

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
HEUR* SCIPsolGetHeur(
   SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->heur;
}

/** outputs non-zero elements of solution to file stream */
RETCODE SCIPsolPrint(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< problem data */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   Real solval;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);

   if( file == NULL )
      file = stdout;

   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);
      CHECK_OKAY( SCIPsolGetVal(sol, set, stat, prob->vars[v], &solval) );
      if( !SCIPsetIsZero(set, solval) )
      {
         fprintf(file, "%-32s %f\n", prob->vars[v]->name, solval);
      }
   }

   return SCIP_OKAY;
}

