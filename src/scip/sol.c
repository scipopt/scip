/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sol.c,v 1.29 2004/01/26 15:10:18 bzfpfend Exp $"

/**@file   sol.c
 * @brief  methods and datastructures for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "stat.h"
#include "clock.h"
#include "misc.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
#include "sol.h"
#include "tree.h"
#include "cons.h"

#include "struct_sol.h"



/** creates primal CIP solution, initialized to zero */
RETCODE SCIPsolCreate(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);
   assert(tree != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, sol) );   
   CHECK_OKAY( SCIPrealarrayCreate(&(*sol)->vals, memhdr) );
   (*sol)->valid = NULL;
   (*sol)->heur = heur;
   (*sol)->obj = 0.0;
   (*sol)->time = SCIPclockGetTime(stat->solvingtime);
   (*sol)->nodenum = stat->nnodes;
   (*sol)->solorigin = SCIP_SOLORIGIN_ZERO;
   (*sol)->depth = (tree->actnode != NULL ? SCIPnodeGetDepth(tree->actnode) : -1);

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
   assert((sourcesol->solorigin == SCIP_SOLORIGIN_ZERO) == (sourcesol->valid == NULL));

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
   (*sol)->time = sourcesol->time;
   (*sol)->nodenum = sourcesol->nodenum;
   (*sol)->solorigin = sourcesol->solorigin;
   (*sol)->depth = sourcesol->depth;

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the actual LP solution */
RETCODE SCIPsolCreateLPSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< actual LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   debugMessage("creating solution from LP\n");

   CHECK_OKAY( SCIPsolCreate(sol, memhdr, stat, tree, heur) );
   CHECK_OKAY( SCIPsolLinkLPSol(*sol, memhdr, set, stat, tree, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the actual pseudo solution */
RETCODE SCIPsolCreatePseudoSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);

   debugMessage("creating solution from pseudo solution\n");

   CHECK_OKAY( SCIPsolCreate(sol, memhdr, stat, tree, heur) );
   CHECK_OKAY( SCIPsolLinkPseudoSol(*sol, memhdr, set, stat, tree, lp) );

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
      CHECK_OKAY( SCIPsolCreateLPSol(sol, memhdr, set, stat, tree, lp, heur) );
   }
   else
   {
      CHECK_OKAY( SCIPsolCreatePseudoSol(sol, memhdr, set, stat, tree, lp, heur) );
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
   assert(((*sol)->solorigin == SCIP_SOLORIGIN_ZERO) == ((*sol)->valid == NULL));

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
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->diving || !lp->divingobjchg);

   debugMessage("linking solution to LP\n");

   /* clear the old solution */
   CHECK_OKAY( SCIPrealarrayClear(sol->vals) );
   if( sol->solorigin == SCIP_SOLORIGIN_ZERO )
   {
      assert(sol->valid == NULL);
      CHECK_OKAY( SCIPboolarrayCreate(&sol->valid, memhdr) );
   }
   else
   {
      assert(sol->valid != NULL);
      CHECK_OKAY( SCIPboolarrayClear(sol->valid) );
   }

   /* link solution to LP solution */
   if( lp->divingobjchg )
   {
      /* the objective value has to be calculated manually, because the LP's value is invalid;
       * use objective values of variables, because columns objective values are changed to dive values
       */
      sol->obj = SCIPlpGetLooseObjval(lp, set);
      if( !SCIPsetIsInfinity(set, -sol->obj) )
      {
         VAR* var;
         COL** cols;
         int ncols;
         int c;
         
         cols = SCIPlpGetCols(lp);
         ncols = SCIPlpGetNCols(lp);
         for( c = 0; c < ncols; ++c )
         {
            var = SCIPcolGetVar(cols[c]);
            sol->obj += SCIPvarGetObj(var) * cols[c]->primsol;
         }
      }
   }
   else
   {
      /* the objective value in the columns is correct, s.t. the LP's objective value is also correct */
      sol->obj = SCIPlpGetObjval(lp, set);
   }
   sol->solorigin = SCIP_SOLORIGIN_LPSOL;
   sol->time = SCIPclockGetTime(stat->solvingtime);
   sol->nodenum = stat->nnodes;
   sol->depth = (tree->actnode != NULL ? SCIPnodeGetDepth(tree->actnode) : -1);

   debugMessage(" -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies actual pseudo solution into CIP solution by linking */
RETCODE SCIPsolLinkPseudoSol(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);

   debugMessage("linking solution to pseudo solution\n");

   /* clear the old solution */
   CHECK_OKAY( SCIPrealarrayClear(sol->vals) );
   if( sol->solorigin == SCIP_SOLORIGIN_ZERO )
   {
      assert(sol->valid == NULL);
      CHECK_OKAY( SCIPboolarrayCreate(&sol->valid, memhdr) );
   }
   else
   {
      assert(sol->valid != NULL);
      CHECK_OKAY( SCIPboolarrayClear(sol->valid) );
   }

   /* link solution to pseudo solution */
   sol->obj = SCIPlpGetPseudoObjval(lp, set);
   sol->solorigin = SCIP_SOLORIGIN_PSEUDOSOL;
   sol->time = SCIPclockGetTime(stat->solvingtime);
   sol->nodenum = stat->nnodes;
   sol->depth = (tree->actnode != NULL ? SCIPnodeGetDepth(tree->actnode) : -1);

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
      CHECK_OKAY( SCIPsolLinkLPSol(sol, memhdr, set, stat, tree, lp) );
   }
   else
   {
      CHECK_OKAY( SCIPsolLinkPseudoSol(sol, memhdr, set, stat, tree, lp) );
   }

   return SCIP_OKAY;
}

/** clears primal CIP solution */
RETCODE SCIPsolClear(
   SOL*             sol,                /**< primal CIP solution */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree                /**< branch-and-bound tree */
   )
{
   assert(sol != NULL);
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) == (sol->valid == NULL));
   assert(stat != NULL);
   assert(tree != NULL);

   CHECK_OKAY( SCIPrealarrayClear(sol->vals) );
   sol->obj = 0.0;
   if( sol->solorigin != SCIP_SOLORIGIN_ZERO )
   {
      CHECK_OKAY( SCIPboolarrayFree(&sol->valid) );
      sol->solorigin = SCIP_SOLORIGIN_ZERO;
   }
   sol->time = SCIPclockGetTime(stat->solvingtime);
   sol->nodenum = stat->nnodes;
   sol->depth = (tree->actnode != NULL ? SCIPnodeGetDepth(tree->actnode) : -1);

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
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) == (sol->valid == NULL));
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

   switch( sol->solorigin )
   {
   case SCIP_SOLORIGIN_ZERO:
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_LPSOL:
      /*debugMessage("completing variable <%s> in LP solution %p\n", SCIPvarGetName(var), sol);*/
      if( !SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var)) )
      {
         assert(SCIPrealarrayGetVal(sol->vals, SCIPvarGetIndex(var)) == 0.0);
         CHECK_OKAY( SCIPrealarraySetVal(sol->vals, set, SCIPvarGetIndex(var), SCIPvarGetLPSol(var)) );
         CHECK_OKAY( SCIPboolarraySetVal(sol->valid, set, SCIPvarGetIndex(var), TRUE) );
      }
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_PSEUDOSOL:
      /*debugMessage("completing variable <%s> in pseudo solution %p\n", SCIPvarGetName(var), sol);*/
      if( !SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var)) )
      {
         assert(SCIPrealarrayGetVal(sol->vals, SCIPvarGetIndex(var)) == 0.0);
         CHECK_OKAY( SCIPrealarraySetVal(sol->vals, set, SCIPvarGetIndex(var), SCIPvarGetPseudoSol(var)) );
         CHECK_OKAY( SCIPboolarraySetVal(sol->valid, set, SCIPvarGetIndex(var), TRUE) );
      }
      return SCIP_OKAY;

   default:
      errorMessage("unknown variable status\n");
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
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) == (sol->valid == NULL));

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
   TREE*            tree,               /**< branch-and-bound tree */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   )
{
   Real oldval;

   assert(sol != NULL);
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) == (sol->valid == NULL));
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO || sol->nodenum == stat->nnodes);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(var != NULL);

   debugMessage("setting value of <%s> in solution %p to %g\n", SCIPvarGetName(var), sol, val);

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetTransVar(var), val);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      CHECK_OKAY( SCIPsolGetVal(sol, set, stat, var, &oldval) );
      if( !SCIPsetIsEQ(set, val, oldval) )
      {
         CHECK_OKAY( SCIPrealarraySetVal(sol->vals, set, SCIPvarGetIndex(var), val) );
         sol->obj += SCIPvarGetObj(var) * (val - oldval);
         sol->time = SCIPclockGetTime(stat->solvingtime);
         sol->nodenum = stat->nnodes;
         sol->depth = (tree->actnode != NULL ? SCIPnodeGetDepth(tree->actnode) : -1);
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot set solution value for fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPsetIsZero(set, SCIPvarGetAggrScalar(var)));
      return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetAggrVar(var),
         (val - SCIPvarGetAggrConstant(var))/SCIPvarGetAggrScalar(var));

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot set solution value for multiple aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), SCIPvarGetNegationConstant(var) - val);
      
   default:
      errorMessage("unknown variable status\n");
      abort();
   }
}

/** increases value of variable in primal CIP solution */
RETCODE SCIPsolIncVal(
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   VAR*             var,                /**< variable to increase solution value for */
   Real             incval              /**< increment for solution value of variable */
   )
{
   assert(sol != NULL);
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) == (sol->valid == NULL));
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO || sol->nodenum == stat->nnodes);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(var != NULL);

   debugMessage("increasing value of <%s> in solution %p by %g\n", SCIPvarGetName(var), sol, incval);

   if( SCIPsetIsZero(set, incval) )
      return SCIP_OKAY;

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      return SCIPsolIncVal(sol, set, stat, tree, SCIPvarGetTransVar(var), incval);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      CHECK_OKAY( solUnlinkVar(sol, set, var) );
      CHECK_OKAY( SCIPrealarrayIncVal(sol->vals, set, SCIPvarGetIndex(var), incval) );
      sol->obj += SCIPvarGetObj(var) * incval;
      sol->time = SCIPclockGetTime(stat->solvingtime);
      sol->nodenum = stat->nnodes;
      sol->depth = (tree->actnode != NULL ? SCIPnodeGetDepth(tree->actnode) : -1);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      errorMessage("cannot set solution value for fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPsetIsZero(set, SCIPvarGetAggrScalar(var)));
      return SCIPsolIncVal(sol, set, stat, tree, SCIPvarGetAggrVar(var), incval/SCIPvarGetAggrScalar(var));

   case SCIP_VARSTATUS_MULTAGGR:
      errorMessage("cannot set solution value for multiple aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPsolIncVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), -incval);

   default:
      errorMessage("unknown variable status\n");
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
   VAR** vars;
   Real* scalars;
   Real val;
   int nvars;
   int i;

   assert(sol != NULL);
   assert((sol->solorigin == SCIP_SOLORIGIN_ZERO) == (sol->valid == NULL));
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO || sol->nodenum == stat->nnodes);
   assert(var != NULL);
   assert(solval != NULL);

   /*debugMessage("getting value of <%s> in solution %p\n", SCIPvarGetName(var), sol);*/

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      CHECK_OKAY( SCIPsolGetVal(sol, set, stat, SCIPvarGetTransVar(var), solval) );
      return SCIP_OKAY;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      CHECK_OKAY( solUnlinkVar(sol, set, var) );
      *solval = SCIPrealarrayGetVal(sol->vals, SCIPvarGetIndex(var));
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetUbGlobal(var)); /*lint !e777*/
      assert(SCIPvarGetLbLocal(var) == SCIPvarGetUbLocal(var)); /*lint !e777*/
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetLbLocal(var)); /*lint !e777*/
      *solval = SCIPvarGetLbGlobal(var);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      CHECK_OKAY( SCIPsolGetVal(sol, set, stat, SCIPvarGetAggrVar(var), solval) );
      (*solval) *= SCIPvarGetAggrScalar(var);
      (*solval) += SCIPvarGetAggrConstant(var);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_MULTAGGR:
      nvars = SCIPvarGetMultaggrNVars(var);
      vars = SCIPvarGetMultaggrVars(var);
      scalars = SCIPvarGetMultaggrScalars(var);
      *solval = SCIPvarGetMultaggrConstant(var);
      for( i = 0; i < nvars; ++i )
      {
         CHECK_OKAY( SCIPsolGetVal(sol, set, stat, vars[i], &val) );
         (*solval) += scalars[i] * val;
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_NEGATED:
      CHECK_OKAY( SCIPsolGetVal(sol, set, stat, SCIPvarGetNegationVar(var), solval) );
      (*solval) = SCIPvarGetNegationConstant(var) - (*solval);
      return SCIP_OKAY;

   default:
      errorMessage("unknown variable status\n");
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

   debugMessage("checking solution with objective value %g (nodenum=%lld, origin=%d)\n", 
      sol->obj, sol->nodenum, sol->solorigin);

   *feasible = TRUE;
   for( h = 0; h < set->nconshdlrs && *feasible; ++h )
   {
      CHECK_OKAY( SCIPconshdlrCheck(set->conshdlrs[h], memhdr, set, prob, sol, chckintegrality, chcklprows, &result) );
      *feasible = *feasible && (result == SCIP_FEASIBLE);
   }

   return SCIP_OKAY;
}

/** try to round given solution */
RETCODE SCIPsolRound(
   SOL*             sol,                /**< primal solution */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch-and-bound tree */
   Bool*            success             /**< pointer to store whether rounding was successful */
   )
{
   int nvars;
   int v;

   assert(sol != NULL);
   assert(success != NULL);

   /* store solution in its own arrays */
   CHECK_OKAY( SCIPsolUnlink(sol, set, prob) );

   /* round all roundable fractional variables in the corresponding direction as long as no unroundable var was found */
   nvars = prob->nbinvars + prob->nintvars;
   for( v = 0; v < nvars; ++v )
   {
      VAR* var;
      Real solval;
      Bool mayrounddown;
      Bool mayroundup;

      var = prob->vars[v];
      CHECK_OKAY( SCIPsolGetVal(sol, set, stat, var, &solval) );

      /* if solution value is already integral, there is nothing to do */
      if( SCIPsetIsIntegral(set, solval) )
         continue;

      /* get rounding possibilities */
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);

      /* choose rounding direction */
      if( mayrounddown && mayroundup )
      {
         /* we can round in both directions: round in objective function direction */
         if( SCIPvarGetObj(var) >= 0.0 )
            solval = SCIPsetFloor(set, solval);
         else
            solval = SCIPsetCeil(set, solval);
      }
      else if( mayrounddown )
         solval = SCIPsetFloor(set, solval);
      else if( mayroundup )
         solval = SCIPsetCeil(set, solval);
      else
         break;

      /* store new solution value */
      CHECK_OKAY( SCIPsolSetVal(sol, set, stat, tree, var, solval) );
   }

   /* check, if rounding was successful */
   *success = (v == nvars);

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

/** gets clock time, when this solution was found */
Real SCIPsolGetTime(
   SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->time;
}

/** gets node number, where this solution was found */
Longint SCIPsolGetNodenum(
   SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->nodenum;
}

/** gets node's depth, where this solution was found */
int SCIPsolGetDepth(
   SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->depth;
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
         if( SCIPsetIsInfinity(set, solval) )
            fprintf(file, "%-32s +infinity\n", SCIPvarGetName(prob->vars[v]));
         else if( SCIPsetIsInfinity(set, -solval) )
            fprintf(file, "%-32s -infinity\n", SCIPvarGetName(prob->vars[v]));
         else
            fprintf(file, "%-32s %f \t(obj:%g)\n", SCIPvarGetName(prob->vars[v]), solval, SCIPvarGetObj(prob->vars[v]));
      }
   }

   return SCIP_OKAY;
}

