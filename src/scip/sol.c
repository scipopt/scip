/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sol.c,v 1.49 2005/01/13 16:20:49 bzfpfend Exp $"

/**@file   sol.c
 * @brief  methods for storing primal CIP solutions
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
#include "primal.h"
#include "tree.h"
#include "cons.h"

#ifndef NDEBUG
#include "struct_sol.h"
#endif



/** clears solution arrays of primal CIP solution */
static
RETCODE solClearArrays(
   SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   CHECK_OKAY( SCIPrealarrayClear(sol->vals) );
   CHECK_OKAY( SCIPboolarrayClear(sol->valid) );

   return SCIP_OKAY;
}

/** sets value of variable in the solution's array */
static
RETCODE solSetArrayVal(
   SOL*             sol,                /**< primal CIP solution */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable */
   Real             val                 /**< value to set variable to */
   )
{
   int index;

   assert(sol != NULL);

   index = SCIPvarGetIndex(var);

   /* mark the variable valid */
   CHECK_OKAY( SCIPboolarraySetVal(sol->valid, set, index, TRUE) );

   /* set the value in the solution array */
   CHECK_OKAY( SCIPrealarraySetVal(sol->vals, set, index, val) );

   return SCIP_OKAY;
}

/** increases value of variable in the solution's array */
static
RETCODE solIncArrayVal(
   SOL*             sol,                /**< primal CIP solution */
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable */
   Real             incval              /**< increase of variable's solution value */
   )
{
   int index;

   assert(sol != NULL);

   index = SCIPvarGetIndex(var);

   /* mark the variable valid */
   CHECK_OKAY( SCIPboolarraySetVal(sol->valid, set, index, TRUE) );

   /* increase the value in the solution array */
   CHECK_OKAY( SCIPrealarrayIncVal(sol->vals, set, index, incval) );

   return SCIP_OKAY;
}

/** returns the value of the variable in the given solution */
static
Real solGetArrayVal(
   SOL*             sol,                /**< primal CIP solution */
   VAR*             var                 /**< problem variable */
   )
{
   int index;

   assert(sol != NULL);

   index = SCIPvarGetIndex(var);

   /* check, if the variable's value is valid */
   if( SCIPboolarrayGetVal(sol->valid, index) )
   {
      return SCIPrealarrayGetVal(sol->vals, index);
   }
   else
   {
      assert(SCIPrealarrayGetVal(sol->vals, index) == 0.0);

      /* return the variable's value corresponding to the origin */
      switch( sol->solorigin )
      {
      case SCIP_SOLORIGIN_ZERO:
         return 0.0;

      case SCIP_SOLORIGIN_LPSOL:
         return SCIPvarGetLPSol(var);

      case SCIP_SOLORIGIN_PSEUDOSOL:
         return SCIPvarGetPseudoSol(var);

      default:
         errorMessage("unknown solution origin <%d>\n", sol->solorigin);
         abort();
      }
   }
}

/** stores solution value of variable in solution's own array */
static
RETCODE solUnlinkVar(
   SOL*             sol,                /**< primal CIP solution */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< problem variable */
   )
{
   Real solval;

   assert(sol != NULL);
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

   /* if variable is already valid, nothing has to be done */
   if( SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var)) )
      return SCIP_OKAY;

   assert(SCIPrealarrayGetVal(sol->vals, SCIPvarGetIndex(var)) == 0.0);

   debugMessage("unlinking solution value of variable <%s>\n", SCIPvarGetName(var));

   /* store the correct solution value into the solution array */
   switch( sol->solorigin )
   {
   case SCIP_SOLORIGIN_ZERO:
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_LPSOL:
      solval = SCIPvarGetLPSol(var);
      if( !SCIPsetIsZero(set, solval) )
      {
         CHECK_OKAY( solSetArrayVal(sol, set, var, solval) );
      }
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_PSEUDOSOL:
      solval = SCIPvarGetPseudoSol(var);
      if( !SCIPsetIsZero(set, solval) )
      {
         CHECK_OKAY( solSetArrayVal(sol, set, var, solval) );
      }
      return SCIP_OKAY;

   default:
      errorMessage("unknown solution origin <%d>\n", sol->solorigin);
      return SCIP_INVALIDDATA;
   }
}

/** sets the solution time, nodenum, runnum, and depth stamp to the current values */
static
void solStamp(
   SOL*             sol,                /**< primal CIP solution */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);

   sol->time = SCIPclockGetTime(stat->solvingtime);
   sol->nodenum = stat->nnodes;
   sol->runnum = stat->nruns;
   sol->depth = SCIPtreeGetCurrentDepth(tree);
}

/** creates primal CIP solution, initialized to zero */
RETCODE SCIPsolCreate(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(memhdr != NULL);
   assert(stat != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, sol) );   
   CHECK_OKAY( SCIPrealarrayCreate(&(*sol)->vals, memhdr) );
   CHECK_OKAY( SCIPboolarrayCreate(&(*sol)->valid, memhdr) );
   (*sol)->heur = heur;
   (*sol)->solorigin = SCIP_SOLORIGIN_ZERO;
   (*sol)->obj = 0.0;
   (*sol)->primalindex = -1;
   solStamp(*sol, stat, tree);

   CHECK_OKAY( SCIPprimalSolCreated(primal, memhdr, set, *sol) );

   return SCIP_OKAY;
}

/** creates a copy of a primal CIP solution */
RETCODE SCIPsolCopy(
   SOL**            sol,                /**< pointer to store the copy of the primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   PRIMAL*          primal,             /**< primal data */
   SOL*             sourcesol           /**< primal CIP solution to copy */
   )
{
   assert(sol != NULL);
   assert(sourcesol != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, sol) );   
   CHECK_OKAY( SCIPrealarrayCopy(&(*sol)->vals, memhdr, sourcesol->vals) );
   CHECK_OKAY( SCIPboolarrayCopy(&(*sol)->valid, memhdr, sourcesol->valid) );
   (*sol)->heur = sourcesol->heur;
   (*sol)->obj = sourcesol->obj;
   (*sol)->primalindex = -1;
   (*sol)->time = sourcesol->time;
   (*sol)->nodenum = sourcesol->nodenum;
   (*sol)->solorigin = sourcesol->solorigin;
   (*sol)->runnum = sourcesol->runnum;
   (*sol)->depth = sourcesol->depth;

   CHECK_OKAY( SCIPprimalSolCreated(primal, memhdr, set, *sol) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current LP solution */
RETCODE SCIPsolCreateLPSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(lp != NULL);
   assert(lp->solved);

   CHECK_OKAY( SCIPsolCreate(sol, memhdr, set, stat, primal, tree, heur) );
   CHECK_OKAY( SCIPsolLinkLPSol(*sol, memhdr, set, stat, tree, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current pseudo solution */
RETCODE SCIPsolCreatePseudoSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);

   CHECK_OKAY( SCIPsolCreate(sol, memhdr, set, stat, primal, tree, heur) );
   CHECK_OKAY( SCIPsolLinkPseudoSol(*sol, memhdr, set, stat, tree, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current solution */
RETCODE SCIPsolCreateCurrentSol(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PRIMAL*          primal,             /**< primal data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(tree != NULL);

   if( SCIPtreeHasCurrentNodeLP(tree) )
   {
      CHECK_OKAY( SCIPsolCreateLPSol(sol, memhdr, set, stat, primal, tree, lp, heur) );
   }
   else
   {
      CHECK_OKAY( SCIPsolCreatePseudoSol(sol, memhdr, set, stat, primal, tree, lp, heur) );
   }

   return SCIP_OKAY;
}

/** frees primal CIP solution */
RETCODE SCIPsolFree(
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   PRIMAL*          primal              /**< primal data */
   )
{
   assert(sol != NULL);
   assert(*sol != NULL);

   SCIPprimalSolFreed(primal, *sol);

   CHECK_OKAY( solClearArrays(*sol) );
   CHECK_OKAY( SCIPrealarrayFree(&(*sol)->vals) );
   CHECK_OKAY( SCIPboolarrayFree(&(*sol)->valid) );
   freeBlockMemory(memhdr, sol);

   return SCIP_OKAY;
}

/** informs the solution that it now belongs to the given primal heuristic */
void SCIPsolSetHeur(
   SOL*             sol,                /**< primal CIP solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);

   sol->heur = heur;
}

/** copies current LP solution into CIP solution by linking */
RETCODE SCIPsolLinkLPSol(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(SCIPlpDiving(lp) || !SCIPlpDivingObjChanged(lp));

   debugMessage("linking solution to LP\n");

   /* clear the old solution arrays */
   CHECK_OKAY( solClearArrays(sol) );

   /* link solution to LP solution */
   if( SCIPlpDivingObjChanged(lp) )
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
   solStamp(sol, stat, tree);

   debugMessage(" -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current pseudo solution into CIP solution by linking */
RETCODE SCIPsolLinkPseudoSol(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);

   debugMessage("linking solution to pseudo solution\n");

   /* clear the old solution arrays */
   CHECK_OKAY( solClearArrays(sol) );

   /* link solution to pseudo solution */
   sol->obj = SCIPlpGetPseudoObjval(lp, set);
   sol->solorigin = SCIP_SOLORIGIN_PSEUDOSOL;
   solStamp(sol, stat, tree);

   debugMessage(" -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current solution (LP or pseudo solution) into CIP solution by linking */
RETCODE SCIPsolLinkCurrentSol(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);

   debugMessage("linking solution to current solution\n");

   if( SCIPtreeHasCurrentNodeLP(tree) )
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
   TREE*            tree                /**< branch and bound tree */
   )
{
   assert(sol != NULL);

   CHECK_OKAY( solClearArrays(sol) );
   sol->solorigin = SCIP_SOLORIGIN_ZERO;
   sol->obj = 0.0;
   solStamp(sol, stat, tree);

   return SCIP_OKAY;
}

/** stores solution values of variables in solution's own array */
RETCODE SCIPsolUnlink(
   SOL*             sol,                /**< primal CIP solution */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   )
{
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(prob->nvars == 0 || prob->vars != NULL);

   if( sol->solorigin != SCIP_SOLORIGIN_ZERO )
   {
      debugMessage("completing solution %p\n", sol);

      for( v = 0; v < prob->nvars; ++v )
      {
         CHECK_OKAY( solUnlinkVar(sol, set, prob->vars[v]) );
      }
      
      sol->solorigin = SCIP_SOLORIGIN_ZERO;
   }

   return SCIP_OKAY;
}

/** sets value of variable in primal CIP solution */
RETCODE SCIPsolSetVal(
   SOL*             sol,                /**< primal CIP solution */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   )
{
   Real oldval;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO || (sol->nodenum == stat->nnodes && sol->runnum == stat->nruns));
   assert(stat != NULL);
   assert(var != NULL);

   debugMessage("setting value of <%s> in solution %p to %g\n", SCIPvarGetName(var), sol, val);

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetTransVar(var), val);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      oldval = solGetArrayVal(sol, var);
      if( !SCIPsetIsEQ(set, val, oldval) )
      {
         CHECK_OKAY( solSetArrayVal(sol, set, var, val) );
         sol->obj += SCIPvarGetObj(var) * (val - oldval);
         solStamp(sol, stat, tree);
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      oldval = SCIPvarGetLbGlobal(var);
      if( !SCIPsetIsEQ(set, val, oldval) )
      {
         errorMessage("cannot set solution value for variable <%s> fixed to %g to different value %g\n",
            SCIPvarGetName(var), oldval, val);
         return SCIP_INVALIDDATA;
      }
      return SCIP_OKAY;

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
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   VAR*             var,                /**< variable to increase solution value for */
   Real             incval              /**< increment for solution value of variable */
   )
{
   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO || (sol->nodenum == stat->nnodes && sol->runnum == stat->nruns));
   assert(stat != NULL);
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
      if( !SCIPsetIsZero(set, incval) )
      {
         CHECK_OKAY( solIncArrayVal(sol, set, var, incval) );
         sol->obj += SCIPvarGetObj(var) * incval;
         solStamp(sol, stat, tree);
      }
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
Real SCIPsolGetVal(
   SOL*             sol,                /**< primal CIP solution */
   STAT*            stat,               /**< problem statistics data */
   VAR*             var                 /**< variable to get value for */
   )
{
   VAR** vars;
   Real* scalars;
   Real solval;
   int nvars;
   int i;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO || (sol->nodenum == stat->nnodes && sol->runnum == stat->nruns));
   assert(var != NULL);

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      return SCIPsolGetVal(sol, stat, SCIPvarGetTransVar(var));

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return solGetArrayVal(sol, var);

   case SCIP_VARSTATUS_FIXED:
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetUbGlobal(var)); /*lint !e777*/
      assert(SCIPvarGetLbLocal(var) == SCIPvarGetUbLocal(var)); /*lint !e777*/
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetLbLocal(var)); /*lint !e777*/
      return SCIPvarGetLbGlobal(var);

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      return SCIPvarGetAggrScalar(var) * SCIPsolGetVal(sol, stat, SCIPvarGetAggrVar(var)) + SCIPvarGetAggrConstant(var);

   case SCIP_VARSTATUS_MULTAGGR:
      nvars = SCIPvarGetMultaggrNVars(var);
      vars = SCIPvarGetMultaggrVars(var);
      scalars = SCIPvarGetMultaggrScalars(var);
      solval = SCIPvarGetMultaggrConstant(var);
      for( i = 0; i < nvars; ++i )
         solval += scalars[i] * SCIPsolGetVal(sol, stat, vars[i]);
      return solval;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPvarGetNegationConstant(var) - SCIPsolGetVal(sol, stat, SCIPvarGetNegationVar(var));

   default:
      errorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** updates primal solutions after a change in a variable's objective value */
void SCIPsolUpdateVarObj(
   SOL*             sol,                /**< primal CIP solution */
   VAR*             var,                /**< problem variable */
   Real             oldobj,             /**< old objective value */
   Real             newobj              /**< new objective value */
   )
{
   Real solval;

   assert(sol != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   solval = solGetArrayVal(sol, var);
   sol->obj += (newobj - oldobj) * solval;
}

/** checks primal CIP solution for feasibility */
RETCODE SCIPsolCheck(
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
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
      CHECK_OKAY( SCIPconshdlrCheck(set->conshdlrs[h], memhdr, set, stat, prob, sol, 
            checkintegrality, checklprows, &result) );
      *feasible = *feasible && (result == SCIP_FEASIBLE);
   }

#ifdef DEBUG
   if( !(*feasible) )
      printf("  -> infeasibility detected in constraint handler <%s>\n",SCIPconshdlrGetName(set->conshdlrs[h-1])); 
#endif

   return SCIP_OKAY;
}

/** try to round given solution */
RETCODE SCIPsolRound(
   SOL*             sol,                /**< primal solution */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   Bool*            success             /**< pointer to store whether rounding was successful */
   )
{
   int nvars;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(prob->transformed);
   assert(success != NULL);

   /* round all roundable fractional variables in the corresponding direction as long as no unroundable var was found */
   nvars = prob->nbinvars + prob->nintvars;
   for( v = 0; v < nvars; ++v )
   {
      VAR* var;
      Real solval;
      Bool mayrounddown;
      Bool mayroundup;

      var = prob->vars[v];
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      solval = solGetArrayVal(sol, var);

      /* if solution value is already integral, there is nothing to do */
      if( SCIPsetIsFeasIntegral(set, solval) )
         continue;

      /* get rounding possibilities */
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);

      /* choose rounding direction */
      if( mayrounddown && mayroundup )
      {
         /* we can round in both directions: round in objective function direction */
         if( SCIPvarGetObj(var) >= 0.0 )
            solval = SCIPsetFeasFloor(set, solval);
         else
            solval = SCIPsetFeasCeil(set, solval);
      }
      else if( mayrounddown )
         solval = SCIPsetFeasFloor(set, solval);
      else if( mayroundup )
         solval = SCIPsetFeasCeil(set, solval);
      else
         break;

      /* store new solution value */
      CHECK_OKAY( SCIPsolSetVal(sol, set, stat, tree, var, solval) );
   }

   /* check, if rounding was successful */
   *success = (v == nvars);

   return SCIP_OKAY;
}

/** updates the solution value sums in variables by adding the value in the given solution */
void SCIPsolUpdateVarsum(
   SOL*             sol,                /**< primal CIP solution */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem data */
   Real             weight              /**< weight of solution in weighted average */
   )
{
   Real solval;
   int v;

   assert(sol != NULL);
   assert(0.0 <= weight && weight <= 1.0);

   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);
      solval = SCIPsolGetVal(sol, stat, prob->vars[v]);
      prob->vars[v]->primsolavg *= (1.0-weight);
      prob->vars[v]->primsolavg += weight*solval;
   }
}

/** outputs non-zero elements of solution to file stream */
RETCODE SCIPsolPrint(
   SOL*             sol,                /**< primal CIP solution */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< problem data (original or transformed) */
   PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   Real solval;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(prob->transformed || transprob != NULL);

   if( file == NULL )
      file = stdout;

   /* display variables of problem data */
   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);
      solval = SCIPsolGetVal(sol, stat, prob->vars[v]);
      if( !SCIPsetIsZero(set, solval) )
      {
         fprintf(file, "%-32s", SCIPvarGetName(prob->vars[v]));
         if( SCIPsetIsInfinity(set, solval) )
            fprintf(file, " +infinity");
         else if( SCIPsetIsInfinity(set, -solval) )
            fprintf(file, " -infinity");
         else
            fprintf(file, " %f", solval);
         fprintf(file, " \t(obj:%g)\n", SCIPvarGetObj(prob->vars[v]));
      }
   }

   /* display additional priced variables (if given problem data is original problem) */
   if( !prob->transformed )
   {
      for( v = transprob->startnvars; v < transprob->nvars; ++v )
      {
         assert(transprob->vars[v] != NULL);
         solval = SCIPsolGetVal(sol, stat, transprob->vars[v]);
         if( !SCIPsetIsZero(set, solval) )
         {
            fprintf(file, "%-32s", SCIPvarGetName(transprob->vars[v]));
            if( SCIPsetIsInfinity(set, solval) )
               fprintf(file, " +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               fprintf(file, " -infinity");
            else
               fprintf(file, " %f", solval);
            fprintf(file, " \t(obj:%g)\n", SCIPvarGetObj(transprob->vars[v]));
         }
      }
   }

   return SCIP_OKAY;
}




/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPsolGetObj
#undef SCIPsolGetTime
#undef SCIPsolGetNodenum
#undef SCIPsolGetRunnum
#undef SCIPsolGetDepth
#undef SCIPsolGetHeur
#undef SCIPsolGetPrimalIndex
#undef SCIPsolSetPrimalIndex

/** gets objective value of primal CIP solution in transformed problem */
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

/** gets branch and bound run number, where this solution was found */
int SCIPsolGetRunnum(
   SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->runnum;
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

/** gets current position of solution in array of existing solutions of primal data */
int SCIPsolGetPrimalIndex(
   SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->primalindex;
}

/** sets current position of solution in array of existing solutions of primal data */
void SCIPsolSetPrimalIndex(
   SOL*             sol,                /**< primal CIP solution */
   int              primalindex         /**< new primal index of solution */
   )
{
   assert(sol != NULL);

   sol->primalindex = primalindex;
}
