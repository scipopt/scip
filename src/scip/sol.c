/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sol.c
 * @ingroup OTHER_CFILES
 * @brief  methods for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "scip/clock.h"
#include "scip/cons.h"
#include "scip/lp.h"
#include "scip/lpexact.h"
#include "scip/misc.h"
#include "scip/nlp.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/relax.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/stat.h"
#include "scip/struct_lp.h"
#include "scip/struct_lpexact.h"
#include "scip/struct_prob.h"
#include "scip/struct_set.h"
#include "scip/struct_sol.h"
#include "scip/struct_stat.h"
#include "scip/struct_var.h"
#include "scip/tree.h"
#include "scip/var.h"


/** clears solution arrays of primal CIP solution */
static
SCIP_RETCODE solClearArrays(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   SCIP_CALL( SCIPboolarrayClear(sol->valid) );
   sol->hasinfval = FALSE;

   if( SCIPsolIsExact(sol) )
   {
      SCIP_CALL( SCIPboolarrayClear(sol->valsexact->valid) );
   }

   return SCIP_OKAY;
}

/** sets value of variable in the solution's array */
static
SCIP_RETCODE solSetArrayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             val                 /**< value to set variable to */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* from now on, variable must not be deleted */
   SCIPvarMarkNotDeletable(var);

   /* mark the variable valid */
   SCIP_CALL( SCIPboolarraySetVal(sol->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, TRUE) );

   /* set the value in the solution array */
   SCIP_CALL( SCIPrealarraySetVal(sol->vals, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, val) );

   /* store whether the solution has infinite values assigned to variables */
   if( val != SCIP_UNKNOWN ) /*lint !e777*/
      sol->hasinfval = (sol->hasinfval || SCIPsetIsInfinity(set, val) || SCIPsetIsInfinity(set, -val));

   return SCIP_OKAY;
}

/** sets value of variable in the exact solution's array */
static
SCIP_RETCODE solSetArrayValExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_RATIONAL*        val                 /**< value to set variable to */
   )
{
   int idx;

   assert(sol != NULL);
   assert(SCIPsolIsExact(sol));

   idx = SCIPvarGetIndex(var);

   /* from now on, variable must not be deleted */
   SCIPvarMarkNotDeletable(var);

   /* mark the variable valid */
   SCIP_CALL( SCIPboolarraySetVal(sol->valsexact->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, TRUE) );

   /* set the value in the solution array */
   SCIP_CALL( SCIPrationalarraySetVal(sol->valsexact->vals, idx, val) );

   /* store whether the solution has infinite values assigned to variables */
   if( SCIPrationalIsAbsInfinity(val) ) /*lint !e777*/
      sol->hasinfval = TRUE;

   return SCIP_OKAY;
}

/** increases value of variable in the solution's array */
static
SCIP_RETCODE solIncArrayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             incval              /**< increase of variable's solution value */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* from now on, variable must not be deleted */
   SCIPvarMarkNotDeletable(var);

   /* if the variable was not valid, mark it to be valid and set the value to the incval (it is 0.0 if not valid) */
   if( !SCIPboolarrayGetVal(sol->valid, idx) )
   {
      /* mark the variable valid */
      SCIP_CALL( SCIPboolarraySetVal(sol->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, TRUE) );

      /* set the value in the solution array */
      SCIP_CALL( SCIPrealarraySetVal(sol->vals, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, incval) );
   }
   else
   {
      /* increase the value in the solution array */
      SCIP_CALL( SCIPrealarrayIncVal(sol->vals, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, incval) );
   }

   /* store whether the solution has infinite values assigned to variables */
   incval = SCIPrealarrayGetVal(sol->vals, idx);
   if( incval != SCIP_UNKNOWN ) /*lint !e777*/
      sol->hasinfval = (sol->hasinfval || SCIPsetIsInfinity(set, incval) || SCIPsetIsInfinity(set, -incval));

   return SCIP_OKAY;
}

/** returns the value of the variable in the given solution */
static
SCIP_Real solGetArrayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* check, if the variable's value is valid */
   if( SCIPboolarrayGetVal(sol->valid, idx) )
   {
      return SCIPrealarrayGetVal(sol->vals, idx);
   }
   else
   {
      /* return the variable's value corresponding to the origin */
      switch( sol->solorigin )
      {
      case SCIP_SOLORIGIN_ORIGINAL:
      case SCIP_SOLORIGIN_ZERO:
         return 0.0;

      case SCIP_SOLORIGIN_LPSOL:
         return SCIPvarGetLPSol(var);

      case SCIP_SOLORIGIN_NLPSOL:
         return SCIPvarGetNLPSol(var);

      case SCIP_SOLORIGIN_RELAXSOL:
         return SCIPvarGetRelaxSolTransVar(var);

      case SCIP_SOLORIGIN_PSEUDOSOL:
         return SCIPvarGetPseudoSol(var);

      case SCIP_SOLORIGIN_PARTIAL:
      case SCIP_SOLORIGIN_UNKNOWN:
         return SCIP_UNKNOWN;

      default:
         SCIPerrorMessage("unknown solution origin <%d>\n", sol->solorigin);
         SCIPABORT();
         return 0.0; /*lint !e527*/
      }
   }
}

/** returns the value of the variable in the given exact solution */
static
void solGetArrayValExact(
   SCIP_RATIONAL*        res,                /**< buffer to store result */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   int idx;

   assert(sol != NULL);
   assert(SCIPsolIsExact(sol));

   idx = SCIPvarGetIndex(var);

   /* check, if the variable's value is valid */
   if( SCIPboolarrayGetVal(sol->valsexact->valid, idx) )
   {
      SCIPrationalarrayGetVal(sol->valsexact->vals, idx, res);
   }
   else
   {
      /* return the variable's value corresponding to the origin */
      switch( sol->solorigin )
      {
      case SCIP_SOLORIGIN_ORIGINAL:
      case SCIP_SOLORIGIN_ZERO:
         SCIPrationalSetReal(res, 0.0);
         break;

      case SCIP_SOLORIGIN_LPSOL:
         SCIPvarGetLPSolExact(var, res);
         break;

      case SCIP_SOLORIGIN_PSEUDOSOL:
         SCIPrationalSetRational(res, SCIPvarGetPseudoSolExact(var));
         break;

      case SCIP_SOLORIGIN_PARTIAL:
      case SCIP_SOLORIGIN_UNKNOWN:
      case SCIP_SOLORIGIN_RELAXSOL:
      case SCIP_SOLORIGIN_NLPSOL:
      default:
         SCIPerrorMessage("unknown solution origin <%d>\n", sol->solorigin);
         SCIPABORT();
         return; /*lint !e527*/
      }
   }
}

/** stores solution value of variable in solution's own array */
static
SCIP_RETCODE solUnlinkVar(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real solval;

   assert(sol != NULL);
   assert(var != NULL);
   assert(SCIPvarIsTransformed(var));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

   /* if variable is already valid, nothing has to be done */
   if( SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var)) )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "unlinking solution value of variable <%s>\n", SCIPvarGetName(var));

   /* store the correct solution value into the solution array */
   switch( sol->solorigin )
   {
   case SCIP_SOLORIGIN_ORIGINAL:
      SCIPerrorMessage("cannot unlink solutions of original problem space\n");
      return SCIP_INVALIDDATA;

   case SCIP_SOLORIGIN_ZERO:
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_LPSOL:
      solval = SCIPvarGetLPSol(var);
      SCIP_CALL( solSetArrayVal(sol, set, var, solval) );
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_NLPSOL:
      solval = SCIPvarGetNLPSol(var);
      SCIP_CALL( solSetArrayVal(sol, set, var, solval) );
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_RELAXSOL:
      solval = SCIPvarGetRelaxSolTransVar(var);
      SCIP_CALL( solSetArrayVal(sol, set, var, solval) );
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_PSEUDOSOL:
      solval = SCIPvarGetPseudoSol(var);
      SCIP_CALL( solSetArrayVal(sol, set, var, solval) );
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_PARTIAL:
   case SCIP_SOLORIGIN_UNKNOWN:
      SCIP_CALL( solSetArrayVal(sol, set, var, SCIP_UNKNOWN) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown solution origin <%d>\n", sol->solorigin);
      return SCIP_INVALIDDATA;
   }
}

/** stores solution value of variable in exact solution's own array */
static
SCIP_RETCODE solUnlinkVarExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_RATIONAL* solval;

   assert(sol != NULL);
   assert(var != NULL);
   assert(SCIPsolIsExact(sol));
   assert(SCIPvarIsTransformed(var));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

   /* if variable is already valid, nothing has to be done */
   if( SCIPboolarrayGetVal(sol->valsexact->valid, SCIPvarGetIndex(var)) )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "unlinking solution value of variable <%s>\n", SCIPvarGetName(var));

   /* store the correct solution value into the solution array */
   switch( sol->solorigin )
   {
   case SCIP_SOLORIGIN_ORIGINAL:
      SCIPerrorMessage("cannot unlink solutions of original problem space\n");
      return SCIP_INVALIDDATA;

   case SCIP_SOLORIGIN_ZERO:
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_LPSOL:
      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &solval) );
      SCIPvarGetLPSolExact(var, solval);
      SCIP_CALL( solSetArrayValExact(sol, set, var, solval) );
      SCIPrationalFreeBuffer(set->buffer, &solval);
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_PSEUDOSOL:
      SCIP_CALL( solSetArrayValExact(sol, set, var, SCIPvarGetPseudoSolExact(var)) );
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_RELAXSOL:
   case SCIP_SOLORIGIN_NLPSOL:
   case SCIP_SOLORIGIN_PARTIAL:
   case SCIP_SOLORIGIN_UNKNOWN:
   default:
      SCIPerrorMessage("unknown solution origin <%d>\n", sol->solorigin);
      return SCIP_INVALIDDATA;
   }
}

/** sets the solution time, nodenum, runnum, and depth stamp to the current values */
static
void solStamp(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_Bool             checktime           /**< should the time be updated? */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);

   if( checktime )
   {
      sol->time = SCIPclockGetTime(stat->solvingtime);
#ifndef NDEBUG
      sol->lpcount = stat->lpcount;
#endif
   }
   else
      sol->time = SCIPclockGetLastTime(stat->solvingtime);
   sol->nodenum = stat->nnodes;
   sol->runnum = stat->nruns;
   if( tree == NULL )
      sol->depth = -1;
   else
      sol->depth = SCIPtreeGetCurrentDepth(tree);
}

/** creates primal CIP solution, initialized to zero */
SCIP_RETCODE SCIPsolCreate(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPrealarrayCreate(&(*sol)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valid, blkmem) );

   (*sol)->solorigin = SCIP_SOLORIGIN_ZERO;
   (*sol)->obj = 0.0;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   (*sol)->valsexact = NULL;
#ifndef NDEBUG
   (*sol)->scip = set->scip;
#endif

   SCIPsolResetViolations(*sol);
   stat->solindex++;
   solStamp(*sol, stat, tree, TRUE);
   SCIPsolResetViolations(*sol);

   /* set solution type and creator depending on whether a heuristic or NULL is passed */
   SCIPsolSetHeur(*sol, heur);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** creates primal CIP solution with exact rational values, initialized to zero */
SCIP_RETCODE SCIPsolCreateExact(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   SCIP_CALL( SCIPsolCreate(sol, blkmem, set, stat, primal, tree, heur) );

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &(*sol)->valsexact) );
   SCIP_CALL( SCIPrationalarrayCreate(&(*sol)->valsexact->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valsexact->valid, blkmem) );
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(*sol)->valsexact->obj) );

   assert(SCIPsolIsExact(*sol));

   return SCIP_OKAY;
}

/** creates a copy of exact solution data */
SCIP_RETCODE SCIPvalsExactCopy(
   SCIP_VALSEXACT**      valsexact,          /**< pointer to store the copy of the primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VALSEXACT*       sourcevals          /**< primal CIP solution to copy */
   )
{
   assert(valsexact != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, valsexact) );
   SCIP_CALL( SCIPrationalarrayCopy(&(*valsexact)->vals, blkmem, sourcevals->vals) );
   SCIP_CALL( SCIPboolarrayCopy(&(*valsexact)->valid, blkmem, sourcevals->valid) );
   SCIP_CALL( SCIPrationalCopyBlock(blkmem, &(*valsexact)->obj, sourcevals->obj) );

   return SCIP_OKAY;
}

/** creates primal CIP solution in original problem space, initialized to the offset in the original problem */
SCIP_RETCODE SCIPsolCreateOriginal(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPrealarrayCreate(&(*sol)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valid, blkmem) );
   (*sol)->solorigin = SCIP_SOLORIGIN_ORIGINAL;
   (*sol)->obj = origprob->objoffset;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   (*sol)->valsexact = NULL;
#ifndef NDEBUG
   (*sol)->scip = set->scip;
#endif
   stat->solindex++;
   solStamp(*sol, stat, tree, TRUE);

   /* set solution type and creator depending on whether a heuristic or NULL is passed */
   SCIPsolSetHeur(*sol, heur);

   SCIPsolResetViolations(*sol);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** creates exact primal CIP solution in original problem space, initialized to the offset in the original problem */
SCIP_RETCODE SCIPsolCreateOriginalExact(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   SCIP_CALL( SCIPsolCreateOriginal(sol, blkmem, set, stat, origprob, primal, tree, heur) );

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &(*sol)->valsexact ) );
   SCIP_CALL( SCIPrationalarrayCreate(&(*sol)->valsexact->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valsexact->valid, blkmem) );
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &(*sol)->valsexact->obj) );

   assert(SCIPsolIsExact(*sol));

   return SCIP_OKAY;
}

/** creates a copy of a primal CIP solution */
SCIP_RETCODE SCIPsolCopy(
   SCIP_SOL**            sol,                /**< pointer to store the copy of the primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SOL*             sourcesol           /**< primal CIP solution to copy */
   )
{
   assert(sol != NULL);
   assert(sourcesol != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPrealarrayCopy(&(*sol)->vals, blkmem, sourcesol->vals) );
   SCIP_CALL( SCIPboolarrayCopy(&(*sol)->valid, blkmem, sourcesol->valid) );

   /* copy solution type and creator information */
   switch( sourcesol->type )
   {
   case SCIP_SOLTYPE_UNKNOWN:
   case SCIP_SOLTYPE_LPRELAX:
   case SCIP_SOLTYPE_STRONGBRANCH:
   case SCIP_SOLTYPE_PSEUDO:
      (*sol)->type = sourcesol->type;
      break;
   case SCIP_SOLTYPE_HEUR:
      SCIPsolSetHeur((*sol), SCIPsolGetHeur(sourcesol));
      break;
   case SCIP_SOLTYPE_RELAX:
      SCIPsolSetRelax((*sol), SCIPsolGetRelax(sourcesol));
      break;
   default:
      SCIPerrorMessage("Unknown source solution type %d!\n", sourcesol->type);
      return SCIP_INVALIDDATA;
   }
   (*sol)->obj = sourcesol->obj;
   (*sol)->primalindex = -1;
   (*sol)->time = sourcesol->time;
#ifndef NDEBUG
   (*sol)->lpcount = sourcesol->lpcount;
#endif
   (*sol)->nodenum = sourcesol->nodenum;
   (*sol)->solorigin = sourcesol->solorigin;
   (*sol)->runnum = sourcesol->runnum;
   (*sol)->depth = sourcesol->depth;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = sourcesol->hasinfval;
   stat->solindex++;
   (*sol)->viol.absviolbounds = sourcesol->viol.absviolbounds;
   (*sol)->viol.absviolcons = sourcesol->viol.absviolcons;
   (*sol)->viol.absviolintegrality = sourcesol->viol.absviolintegrality;
   (*sol)->viol.absviollprows = sourcesol->viol.absviollprows;
   (*sol)->viol.relviolbounds = sourcesol->viol.relviolbounds;
   (*sol)->viol.relviolcons = sourcesol->viol.relviolcons;
   (*sol)->viol.relviollprows = sourcesol->viol.relviollprows;
#ifndef NDEBUG
   (*sol)->scip = set->scip;
#endif

   /* copy rational values if solution is exact */
   if( SCIPsolIsExact(sourcesol) )
   {
      SCIP_CALL( SCIPvalsExactCopy( &(*sol)->valsexact, blkmem, sourcesol->valsexact) );
   }
   else
      (*sol)->valsexact = NULL;

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** transformes given original solution to the transformed space; a corresponding transformed solution has to be given
 *  which is copied into the existing solution and freed afterwards
 */
SCIP_RETCODE SCIPsolTransform(
   SCIP_SOL*             sol,                /**< primal CIP solution to change, living in original space */
   SCIP_SOL**            transsol,           /**< pointer to corresponding transformed primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal              /**< primal data */
   )
{  /*lint --e{715}*/
   SCIP_REALARRAY* tmpvals;
   SCIP_BOOLARRAY* tmpvalid;
   SCIP_SOL* tsol;

   assert(sol != NULL);
   assert(transsol != NULL);
   assert(set != NULL);
   assert(SCIPsolIsOriginal(sol));
   assert(sol->primalindex > -1);

   tsol = *transsol;
   assert(tsol != NULL);
   assert(!SCIPsolIsOriginal(tsol));

   /* switch vals and valid arrays; the exisiting solution gets the arrays of the transformed solution;
    * the transformed one gets the original arrays, because they have to be freed anyway and freeing the transsol
    * automatically frees its arrays
    */
   tmpvals = sol->vals;
   tmpvalid = sol->valid;
   sol->vals = tsol->vals;
   sol->valid = tsol->valid;
   tsol->vals = tmpvals;
   tsol->valid = tmpvalid;
   if( SCIPsolIsExact(sol) )
   {
      SCIP_VALSEXACT* tmpvalsexact;
      tmpvalsexact = sol->valsexact;
      sol->valsexact = tsol->valsexact;
      tsol->valsexact = tmpvalsexact;
   }

   /* copy solorigin and objective (should be the same, only to avoid numerical issues);
    * we keep the other statistics of the original solution, since that was the first time that this solution as found
    */
   sol->solorigin = tsol->solorigin;
   sol->obj = tsol->obj;

   SCIP_CALL( SCIPsolFree(transsol, blkmem, primal) );

   return SCIP_OKAY;
}

/** adjusts solution values of implied integral variables in handed solution, solution objective value is not
 *  deteriorated by this method
 */
SCIP_RETCODE SCIPsolAdjustImplicitSolVals(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< either original or transformed problem, depending on sol origin */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             uselprows           /**< should LP row information be considered for none-objective variables */
   )
{
   SCIP_VAR** vars;
   int nimplvarsbegin;
   int nimplvarsend;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);

   /* get number of implied integral variables */
   nimplvarsend = SCIPprobGetNImplVars(prob);

   if( nimplvarsend == 0 )
      return SCIP_OKAY;

   /* get range of implied integral variables */
   vars = SCIPprobGetVars(prob);
   nimplvarsbegin = SCIPprobGetNBinVars(prob) + SCIPprobGetNIntVars(prob);
   nimplvarsend += nimplvarsbegin;

   /* loop over implied integral variables and round them up or down */
   for( v = nimplvarsbegin; v < nimplvarsend; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real solval;
      SCIP_Real obj;
      SCIP_Real newsolval;
      SCIP_Bool roundup;
      SCIP_Bool rounddown;
      int nuplocks;
      int ndownlocks;

      var = vars[v];

      assert( SCIPvarIsImpliedIntegral(var) );
      solval = SCIPsolGetVal(sol, set, stat, var);

      /* we do not need to round integral solution values or those of variables which are not column variables */
      if( SCIPsetIsFeasIntegral(set, solval) || SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
         continue;

      nuplocks = SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL);
      ndownlocks = SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL);
      obj = SCIPvarGetUnchangedObj(var);

      roundup = FALSE;
      rounddown = FALSE;

      /* in case of a non-zero objective coefficient, there is only one possible rounding direction */
      if( SCIPsetIsFeasNegative(set, obj) )
         roundup = TRUE;
      else if( SCIPsetIsFeasPositive(set, obj) )
         rounddown = TRUE;
      else if( uselprows )
      {
         /* determine rounding direction based on row violations */
         SCIP_COL* col;
         SCIP_ROW** rows;
         SCIP_Real* vals;
         int nrows;
         int r;

         col = SCIPvarGetCol(var);
         vals = SCIPcolGetVals(col);
         rows = SCIPcolGetRows(col);
         nrows = SCIPcolGetNNonz(col);

         /* loop over rows and search for equations whose violation can be decreased by rounding */
         for( r = 0; r < nrows && !(roundup && rounddown); ++r )
         {
            SCIP_ROW* row;
            SCIP_Real activity;
            SCIP_Real rhs;
            SCIP_Real lhs;

            row = rows[r];

            if( SCIProwIsLocal(row) || !SCIProwIsInLP(row) )
               continue;

            rhs = SCIProwGetRhs(row);
            lhs = SCIProwGetLhs(row);

            if( SCIPsetIsInfinity(set, rhs) || SCIPsetIsInfinity(set, -lhs) )
               continue;

            activity = SCIProwGetSolActivity(row, set, stat, sol);
            if( SCIPsetIsFeasLE(set, activity, rhs) && SCIPsetIsFeasLE(set, lhs, activity) )
               continue;

            assert(! SCIPsetIsZero(set, vals[r]));
            if( (SCIPsetIsFeasGT(set, activity, rhs) && SCIPsetIsPositive(set, vals[r]))
                  || (SCIPsetIsFeasLT(set, activity, lhs) && SCIPsetIsNegative(set, vals[r])) )
               rounddown = TRUE;
            else
               roundup = TRUE;
         }
      }

      /* in case of a tie, we select the rounding step based on the number of variable locks */
      if( roundup == rounddown )
      {
         rounddown = ndownlocks <= nuplocks;
         roundup = !rounddown;
      }

      /* round the variable up or down */
      if( roundup )
      {
         newsolval = SCIPsetCeil(set, solval);
         assert(SCIPsetIsFeasLE(set, newsolval, SCIPvarGetUbGlobal(var)));
      }
      else
      {
         assert( rounddown ); /* should be true because of the code above */
         newsolval = SCIPsetFloor(set, solval);
         assert(SCIPsetIsFeasGE(set, newsolval, SCIPvarGetLbGlobal(var)));
      }

      SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, var, newsolval) );
   }

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current LP solution */
SCIP_RETCODE SCIPsolCreateLPSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(lp != NULL);
   assert(SCIPlpIsSolved(lp));

   SCIP_CALL( SCIPsolCreate(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolLinkLPSol(*sol, set, stat, prob, tree, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution with exact rational values, initialized to the current exact LP solution
 *  (will use exact safe dual solution if lp was not solved exactly)
 */
SCIP_RETCODE SCIPsolCreateLPSolExact(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(lp != NULL);
   assert(lp->solved);

   SCIP_CALL( SCIPsolCreateExact(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolLinkLPSolExact(*sol, set, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current NLP solution */
SCIP_RETCODE SCIPsolCreateNLPSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(nlp != NULL);

   SCIP_CALL( SCIPsolCreate(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolLinkNLPSol(*sol, stat, tree, nlp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current relaxation solution */
SCIP_RETCODE SCIPsolCreateRelaxSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(relaxation != NULL);
   assert(SCIPrelaxationIsSolValid(relaxation));

   SCIP_CALL( SCIPsolCreate(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolLinkRelaxSol(*sol, set, stat, tree, relaxation) );

   /* update solution type and store relaxator as creator only if no heuristic is specified as creator */
   if( heur == NULL )
      SCIPsolSetRelax(*sol, SCIPrelaxationGetSolRelax(relaxation));

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current pseudo solution */
SCIP_RETCODE SCIPsolCreatePseudoSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);

   SCIP_CALL( SCIPsolCreate(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolLinkPseudoSol(*sol, set, stat, prob, tree, lp) );

   /* update solution type to pseudo solution */
   if( heur == NULL )
      SCIPsolSetPseudo(*sol);

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current exact pseudo solution */
SCIP_RETCODE SCIPsolCreatePseudoSolExact(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);

   SCIP_CALL( SCIPsolCreateExact(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolLinkPseudoSolExact(*sol, set, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current solution */
SCIP_RETCODE SCIPsolCreateCurrentSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(tree != NULL);

   if( SCIPtreeHasCurrentNodeLP(tree) )
   {
      SCIP_CALL( SCIPsolCreateLPSol(sol, blkmem, set, stat, prob, primal, tree, lp, heur) );
   }
   else
   {
      SCIP_CALL( SCIPsolCreatePseudoSol(sol, blkmem, set, stat, prob, primal, tree, lp, heur) );
   }

   return SCIP_OKAY;
}

/** creates primal CIP solution with exact rational values, initialized to the current solution */
SCIP_RETCODE SCIPsolCreateCurrentSolExact(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(tree != NULL);

   if( SCIPtreeHasCurrentNodeLP(tree) )
   {
      assert(lp->solved);
      SCIP_CALL( SCIPsolCreateLPSolExact(sol, blkmem, set, stat, primal, tree, lp, heur) );
   }
   else
   {
      SCIP_CALL( SCIPsolCreatePseudoSolExact(sol, blkmem, set, stat, primal, tree, lp, heur) );
   }

   return SCIP_OKAY;
}

/** creates partial primal CIP solution, initialized to unknown values */
SCIP_RETCODE SCIPsolCreatePartial(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(primal != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPrealarrayCreate(&(*sol)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valid, blkmem) );
   (*sol)->solorigin = SCIP_SOLORIGIN_PARTIAL;
   (*sol)->obj = SCIP_UNKNOWN;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   (*sol)->valsexact = NULL;
#ifndef NDEBUG
   (*sol)->scip = set->scip;
#endif
   stat->solindex++;
   solStamp(*sol, stat, NULL, TRUE);
   SCIPsolResetViolations(*sol);

   /* set solution type and creator depending on whether a heuristic or NULL is passed */
   SCIPsolSetHeur(*sol, heur);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to unknown values */
SCIP_RETCODE SCIPsolCreateUnknown(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sol) );
   SCIP_CALL( SCIPrealarrayCreate(&(*sol)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valid, blkmem) );
   (*sol)->solorigin = SCIP_SOLORIGIN_UNKNOWN;
   (*sol)->obj = 0.0;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   (*sol)->valsexact = NULL;
   stat->solindex++;
   solStamp(*sol, stat, tree, TRUE);
   SCIPsolResetViolations(*sol);

   /* set solution type and creator depending on whether a heuristic or NULL is passed */
   SCIPsolSetHeur(*sol, heur);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** frees exact solution values */
static
SCIP_RETCODE valsExactFree(
   SCIP_VALSEXACT**      valsexact,          /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(valsexact != NULL);
   assert(*valsexact != NULL);

   SCIPrationalFreeBlock(blkmem, &(*valsexact)->obj);
   SCIP_CALL( SCIPrationalarrayFree(&(*valsexact)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayFree(&(*valsexact)->valid) );
   BMSfreeBlockMemory(blkmem, valsexact);

   return SCIP_OKAY;
}

/** frees primal CIP solution */
SCIP_RETCODE SCIPsolFree(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PRIMAL*          primal              /**< primal data */
   )
{
   assert(sol != NULL);
   assert(*sol != NULL);

   SCIPprimalSolFreed(primal, *sol);

   SCIP_CALL( SCIPrealarrayFree(&(*sol)->vals) );
   SCIP_CALL( SCIPboolarrayFree(&(*sol)->valid) );
   if( SCIPsolIsExact(*sol) )
   {
      SCIP_CALL( valsExactFree(&((*sol)->valsexact), blkmem) );
   }
   BMSfreeBlockMemory(blkmem, sol);

   return SCIP_OKAY;
}

/** copies current LP solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkLPSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(SCIPlpDiving(lp) || SCIPtreeProbing(tree) || !SCIPlpDivingObjChanged(lp));

   SCIPsetDebugMsg(set, "linking solution to LP\n");

   /* clear the old solution arrays */
   SCIP_CALL( solClearArrays(sol) );

   /* link solution to LP solution */
   if( SCIPlpDivingObjChanged(lp) )
   {
      /* the objective value has to be calculated manually, because the LP's value is invalid;
       * use objective values of variables, because columns objective values are changed to dive values
       */
      sol->obj = SCIPlpGetLooseObjval(lp, set, prob);
      if( !SCIPsetIsInfinity(set, -sol->obj) )
      {
         SCIP_VAR* var;
         SCIP_COL** cols;
         int ncols;
         int c;

         cols = SCIPlpGetCols(lp);
         ncols = SCIPlpGetNCols(lp);
         for( c = 0; c < ncols; ++c )
         {
            var = SCIPcolGetVar(cols[c]);
            sol->obj += SCIPvarGetUnchangedObj(var) * cols[c]->primsol;
         }
      }
   }
   else
   {
      /* the objective value in the columns is correct, s.t. the LP's objective value is also correct */
      sol->obj = SCIPlpGetObjval(lp, set, prob);
   }
   sol->solorigin = SCIP_SOLORIGIN_LPSOL;
   solStamp(sol, stat, tree, TRUE);

   SCIPsetDebugMsg(set, " -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current exact LP solution into exact CIP solution by linking */
SCIP_RETCODE SCIPsolLinkLPSolExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   )
{
   assert(sol != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(SCIPsolIsExact(sol));

   /* clear the old solution arrays */
   SCIP_CALL( solClearArrays(sol) );

   /* the objective value in the columns is correct, s.t. the LP's objective value is also correct */
   SCIPlpExactGetObjval(lp, set, sol->valsexact->obj);
   sol->obj = SCIPrationalRoundReal(sol->valsexact->obj, SCIP_R_ROUND_UPWARDS);
   sol->solorigin = SCIP_SOLORIGIN_LPSOL;

   return SCIP_OKAY;
}

/** copies current NLP solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkNLPSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(nlp != NULL);
   assert(SCIPnlpHasSolution(nlp));

   SCIPstatDebugMsg(stat, "linking solution to NLP\n");

   /* clear the old solution arrays */
   SCIP_CALL( solClearArrays(sol) );

   /* get objective value of NLP solution */
   if( SCIPnlpIsDivingObjChanged(nlp) )
   {
      /* the objective value has to be calculated manually, because the NLP's value is invalid */

      SCIP_VAR** vars;
      int nvars;
      int v;

      sol->obj = 0.0;

      vars = SCIPnlpGetVars(nlp);
      nvars = SCIPnlpGetNVars(nlp);
      for( v = 0; v < nvars; ++v )
      {
         assert(SCIPvarIsActive(vars[v]));
         sol->obj += SCIPvarGetUnchangedObj(vars[v]) * SCIPvarGetNLPSol(vars[v]);
      }
   }
   else
   {
      sol->obj = SCIPnlpGetObjval(nlp);
   }

   sol->solorigin = SCIP_SOLORIGIN_NLPSOL;
   solStamp(sol, stat, tree, TRUE);

   SCIPstatDebugMsg(stat, " -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current relaxation solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkRelaxSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   )
{  /*lint --e{715}*/
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(relaxation != NULL);
   assert(SCIPrelaxationIsSolValid(relaxation));

   SCIPsetDebugMsg(set, "linking solution to relaxation\n");

   /* clear the old solution arrays */
   SCIP_CALL( solClearArrays(sol) );

   /* the objective value in the columns is correct, s.t. the LP's objective value is also correct */
   sol->obj = SCIPrelaxationGetSolObj(relaxation);
   sol->solorigin = SCIP_SOLORIGIN_RELAXSOL;
   solStamp(sol, stat, tree, TRUE);

   SCIPsetDebugMsg(set, " -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current pseudo solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkPseudoSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(sol != NULL);
   assert(stat != NULL);
   assert(tree != NULL);

   SCIPsetDebugMsg(set, "linking solution to pseudo solution\n");

   /* clear the old solution arrays */
   SCIP_CALL( solClearArrays(sol) );

   /* link solution to pseudo solution */
   sol->obj = SCIPlpGetPseudoObjval(lp, set, prob);
   sol->solorigin = SCIP_SOLORIGIN_PSEUDOSOL;
   solStamp(sol, stat, tree, TRUE);

   SCIPsetDebugMsg(set, " -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current exact pseudo solution into exact CIP solution by linking */
SCIP_RETCODE SCIPsolLinkPseudoSolExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   )
{
   assert(sol != NULL);
   assert(SCIPsolIsExact(sol));

   /* clear the old solution arrays */
   SCIP_CALL( solClearArrays(sol) );

   /* link solution to pseudo solution */
   SCIPlpExactGetPseudoObjval(lp, set, sol->valsexact->obj);

   SCIPsetDebugMsg(set, " -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** copies current solution (LP or pseudo solution) into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkCurrentSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(tree != NULL);

   SCIPsetDebugMsg(set, "linking solution to current solution\n");

   if( SCIPtreeHasCurrentNodeLP(tree) && SCIPlpIsSolved(lp) )
   {
      SCIP_CALL( SCIPsolLinkLPSol(sol, set, stat, prob, tree, lp) );
   }
   else
   {
      SCIP_CALL( SCIPsolLinkPseudoSol(sol, set, stat, prob, tree, lp) );
   }

   return SCIP_OKAY;
}

/** clears primal CIP solution */
SCIP_RETCODE SCIPsolClear(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(sol != NULL);

   SCIP_CALL( solClearArrays(sol) );
   sol->solorigin = SCIP_SOLORIGIN_ZERO;
   sol->obj = 0.0;
   solStamp(sol, stat, tree, TRUE);

   if( SCIPsolIsExact(sol) )
      SCIPrationalSetReal(sol->valsexact->obj, 0.0);

   return SCIP_OKAY;
}

/** declares all entries in the primal CIP solution to be unknown */
SCIP_RETCODE SCIPsolSetUnknown(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(sol != NULL);

   SCIP_CALL( solClearArrays(sol) );
   sol->solorigin = SCIP_SOLORIGIN_UNKNOWN;
   sol->obj = 0.0;
   solStamp(sol, stat, tree, TRUE);

   return SCIP_OKAY;
}

/** stores solution values of variables in solution's own array */
SCIP_RETCODE SCIPsolUnlink(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< transformed problem data */
   )
{
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(prob->nvars == 0 || prob->vars != NULL);

   if( !SCIPsolIsOriginal(sol) && sol->solorigin != SCIP_SOLORIGIN_ZERO
      && sol->solorigin != SCIP_SOLORIGIN_UNKNOWN )
   {
      SCIPsetDebugMsg(set, "completing solution %p\n", (void*)sol);

      for( v = 0; v < prob->nvars; ++v )
      {
         SCIP_CALL( solUnlinkVar(sol, set, prob->vars[v]) );
      }

      sol->solorigin = SCIP_SOLORIGIN_ZERO;
   }

   return SCIP_OKAY;
}

/** stores solution values of variables in exact solution's own array */
SCIP_RETCODE SCIPsolUnlinkExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< transformed problem data */
   )
{
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(prob->nvars == 0 || prob->vars != NULL);
   assert(SCIPsolIsExact(sol));

   if( sol->solorigin != SCIP_SOLORIGIN_ORIGINAL && sol->solorigin != SCIP_SOLORIGIN_ZERO
      && sol->solorigin != SCIP_SOLORIGIN_UNKNOWN )
   {
      SCIPsetDebugMsg(set, "completing solution %p\n", (void*)sol);

      for( v = 0; v < prob->nvars; ++v )
      {
         SCIP_CALL( solUnlinkVarExact(sol, set, prob->vars[v]) );
      }
   }

   SCIP_CALL( SCIPsolUnlink(sol, set, prob) );

   return SCIP_OKAY;
}

/** sets value of variable in primal CIP solution */
SCIP_RETCODE SCIPsolSetVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_Real             val                 /**< solution value of variable */
   )
{
   SCIP_Real oldval;

   assert(sol != NULL);
   assert(stat != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL
      || sol->solorigin == SCIP_SOLORIGIN_ZERO
      || sol->solorigin == SCIP_SOLORIGIN_PARTIAL
      || sol->solorigin == SCIP_SOLORIGIN_UNKNOWN
      || (sol->nodenum == stat->nnodes && sol->runnum == stat->nruns));
   assert(var != NULL);
   assert(var->scip == sol->scip);
   assert(SCIPisFinite(val));

   SCIPsetDebugMsg(set, "setting value of <%s> in solution %p to %g\n", SCIPvarGetName(var), (void*)sol, val);

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( SCIPsolIsOriginal(sol) )
      {
         oldval = solGetArrayVal(sol, var);

         if( val != oldval )  /*lint !e777*/
         {
            SCIP_CALL( solSetArrayVal(sol, set, var, val) );

            /* update the objective value; we do not need to do this for invalid objectives or partial solutions */
            if( sol->obj != SCIP_INVALID && !SCIPsolIsPartial(sol) ) /*lint !e777*/
            {
               SCIP_Real obj;
               SCIP_Real oldobjcont;
               SCIP_Real newobjcont;

               /* an unknown solution value does not count towards the objective */
               obj = SCIPvarGetUnchangedObj(var);
               oldobjcont = (oldval == SCIP_UNKNOWN ? 0.0 : obj * oldval); /*lint !e777*/
               newobjcont = (val == SCIP_UNKNOWN ? 0.0 : obj * val); /*lint !e777*/

               /* we want to use a safe invalid if the contribution exchange contradicts the infinity status of the objective value */
               if( SCIPsetIsInfinity(set, sol->obj) )
               {
                  if( ( SCIPsetIsInfinity(set, oldobjcont) && !SCIPsetIsInfinity(set, newobjcont) )
                     || ( !SCIPsetIsInfinity(set, -oldobjcont) && SCIPsetIsInfinity(set, -newobjcont) ) )
                     sol->obj = SCIP_INVALID;
               }
               else if( SCIPsetIsInfinity(set, -sol->obj) )
               {
                  if( ( SCIPsetIsInfinity(set, -oldobjcont) && !SCIPsetIsInfinity(set, -newobjcont) )
                     || ( !SCIPsetIsInfinity(set, oldobjcont) && SCIPsetIsInfinity(set, newobjcont) ) )
                     sol->obj = SCIP_INVALID;
               }
               /* we want to use a clean infinity if the contribution exchange or the resulting objective hits the infinity bound */
               else
               {
                  if( !SCIPsetIsInfinity(set, MAX(ABS(oldobjcont), ABS(newobjcont))) )
                  {
                     sol->obj -= oldobjcont;
                     sol->obj += newobjcont;

                     if( SCIPsetIsInfinity(set, sol->obj) )
                        sol->obj = SCIPsetInfinity(set);
                     else if( SCIPsetIsInfinity(set, -sol->obj) )
                        sol->obj = -SCIPsetInfinity(set);
                  }
                  else if( !SCIPsetIsInfinity(set, MAX(oldobjcont, -newobjcont)) )
                     sol->obj = SCIPsetInfinity(set);
                  else if( !SCIPsetIsInfinity(set, MAX(-oldobjcont, newobjcont)) )
                     sol->obj = -SCIPsetInfinity(set);
               }
            }

            solStamp(sol, stat, tree, FALSE);
         }
         return SCIP_OKAY;
      }
      else
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetTransVar(var), val);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(!SCIPsolIsOriginal(sol));
      assert(sol->solorigin != SCIP_SOLORIGIN_LPSOL || SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var))
         || sol->lpcount == stat->lpcount);
      oldval = solGetArrayVal(sol, var);

      if( val != oldval )  /*lint !e777*/
      {
         SCIP_CALL( solSetArrayVal(sol, set, var, val) );

         /* update the objective value; we do not need to do this for invalid objectives */
         if( sol->obj != SCIP_INVALID ) /*lint !e777*/
         {
            SCIP_Real obj;
            SCIP_Real oldobjcont;
            SCIP_Real newobjcont;

            /* an unknown solution value does not count towards the objective */
            obj = SCIPvarGetUnchangedObj(var);
            oldobjcont = (oldval == SCIP_UNKNOWN ? 0.0 : obj * oldval); /*lint !e777*/
            newobjcont = (val == SCIP_UNKNOWN ? 0.0 : obj * val); /*lint !e777*/

            /* we want to use a safe invalid if the contribution exchange contradicts the infinity status of the objective value */
            if( SCIPsetIsInfinity(set, sol->obj) )
            {
               if( ( SCIPsetIsInfinity(set, oldobjcont) && !SCIPsetIsInfinity(set, newobjcont) )
                  || ( !SCIPsetIsInfinity(set, -oldobjcont) && SCIPsetIsInfinity(set, -newobjcont) ) )
                  sol->obj = SCIP_INVALID;
            }
            else if( SCIPsetIsInfinity(set, -sol->obj) )
            {
               if( ( SCIPsetIsInfinity(set, -oldobjcont) && !SCIPsetIsInfinity(set, -newobjcont) )
                  || ( !SCIPsetIsInfinity(set, oldobjcont) && SCIPsetIsInfinity(set, newobjcont) ) )
                  sol->obj = SCIP_INVALID;
            }
            /* we want to use a clean infinity if the contribution exchange or the resulting objective hits the infinity bound */
            else
            {
               if( !SCIPsetIsInfinity(set, MAX(ABS(oldobjcont), ABS(newobjcont))) )
               {
                  sol->obj -= oldobjcont;
                  sol->obj += newobjcont;

                  if( SCIPsetIsInfinity(set, sol->obj) )
                     sol->obj = SCIPsetInfinity(set);
                  else if( SCIPsetIsInfinity(set, -sol->obj) )
                     sol->obj = -SCIPsetInfinity(set);
               }
               else if( !SCIPsetIsInfinity(set, MAX(oldobjcont, -newobjcont)) )
                  sol->obj = SCIPsetInfinity(set);
               else if( !SCIPsetIsInfinity(set, MAX(-oldobjcont, newobjcont)) )
                  sol->obj = -SCIPsetInfinity(set);
            }
         }

         solStamp(sol, stat, tree, FALSE);
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      assert(!SCIPsolIsOriginal(sol));
      oldval = SCIPvarGetLbGlobal(var);
      if( val != oldval )  /*lint !e777*/
      {
         SCIPerrorMessage("cannot set solution value for variable <%s> fixed to %.15g to different value %.15g\n",
            SCIPvarGetName(var), oldval, val);
         return SCIP_INVALIDDATA;
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPsetIsZero(set, SCIPvarGetAggrScalar(var)));
      assert(!SCIPsetIsInfinity(set, SCIPvarGetAggrConstant(var)) && !SCIPsetIsInfinity(set, -SCIPvarGetAggrConstant(var)));
      assert(!SCIPsetIsInfinity(set, SCIPvarGetAggrScalar(var)) && !SCIPsetIsInfinity(set, -SCIPvarGetAggrScalar(var)));

      if( val == SCIP_UNKNOWN )/*lint !e777*/
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetAggrVar(var), val);
      if( SCIPsetIsInfinity(set, val) || SCIPsetIsInfinity(set, -val) )
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetAggrVar(var),  SCIPvarGetAggrScalar(var) > 0 ? val : -val);
      else
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetAggrVar(var), (val - SCIPvarGetAggrConstant(var))/SCIPvarGetAggrScalar(var));

   case SCIP_VARSTATUS_MULTAGGR:
      if ( SCIPvarGetMultaggrNVars(var) == 1 )
      {
         SCIP_VAR** multaggrvars;
         SCIP_Real* multaggrscalars;
         SCIP_Real multaggrconstant;

         multaggrvars = SCIPvarGetMultaggrVars(var);
         multaggrscalars = SCIPvarGetMultaggrScalars(var);
         multaggrconstant = SCIPvarGetMultaggrConstant(var);

         if( SCIPsetIsInfinity(set, multaggrconstant) || SCIPsetIsInfinity(set, -multaggrconstant) )
         {
            if( (SCIPsetIsInfinity(set, multaggrconstant) && !SCIPsetIsInfinity(set, val))
               || (SCIPsetIsInfinity(set, -multaggrconstant) && !SCIPsetIsInfinity(set, -val)) )
            {
               SCIPerrorMessage("cannot set solution value for variable <%s> fixed to %.15g to different value %.15g\n",
                  SCIPvarGetName(var), multaggrconstant, val);
               return SCIP_INVALIDDATA;
            }
            return SCIP_OKAY;
         }
         else
         {
            if( SCIPsetIsInfinity(set, val) || SCIPsetIsInfinity(set, -val) )
               return SCIPsolSetVal(sol, set, stat, tree, multaggrvars[0],  multaggrscalars[0] > 0 ? val : -val);
            else
               return SCIPsolSetVal(sol, set, stat, tree, multaggrvars[0], (val - multaggrconstant)/multaggrscalars[0]);
         }
      }
      SCIPerrorMessage("cannot set solution value for multiple aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      assert(!SCIPsetIsInfinity(set, SCIPvarGetNegationConstant(var)) && !SCIPsetIsInfinity(set, -SCIPvarGetNegationConstant(var)));

      if( val == SCIP_UNKNOWN )/*lint !e777*/
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), val);
      else if( SCIPsetIsInfinity(set, val) || SCIPsetIsInfinity(set, -val) )
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), -val);
      else
         return SCIPsolSetVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), SCIPvarGetNegationConstant(var) - val);

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** sets value of variable in exact primal CIP solution */
SCIP_RETCODE SCIPsolSetValExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_RATIONAL*        val                 /**< solution value of variable */
   )
{
   SCIP_RATIONAL* oldval;
   SCIP_RATIONAL* tmp;
   SCIP_RETCODE retcode;

   assert(sol != NULL);
   assert(stat != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL
      || sol->solorigin == SCIP_SOLORIGIN_ZERO
      || sol->solorigin == SCIP_SOLORIGIN_PARTIAL
      || sol->solorigin == SCIP_SOLORIGIN_UNKNOWN);
   assert(var != NULL);
   assert(!SCIPrationalIsAbsInfinity(val));
   assert(SCIPsolIsExact(sol));

   SCIPsetDebugMsg(set, "setting value of <%s> in exact solution %p to %g\n", SCIPvarGetName(var), (void*)sol, SCIPrationalGetReal(val));

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( SCIPvarGetStatusExact(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( sol->solorigin == SCIP_SOLORIGIN_ORIGINAL )
      {
         SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &oldval) );

         solGetArrayValExact(oldval, sol, var);

         if( !SCIPrationalIsEQ(val, oldval) )
         {
            SCIP_RATIONAL* obj;

            SCIP_CALL( solSetArrayValExact(sol, set, var, val) );
            obj = SCIPvarGetObjExact(var);
            SCIPrationalDiffProd(sol->valsexact->obj, obj, oldval);

            SCIPrationalAddProd(sol->valsexact->obj, obj, val);
         }

         SCIPrationalFreeBuffer(set->buffer, &oldval);
         return SCIP_OKAY;
      }
      else
         return SCIPsolSetValExact(sol, set, stat, tree, SCIPvarGetTransVar(var), val);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);
      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &oldval) );

      solGetArrayValExact(oldval, sol, var);

      if( !SCIPrationalIsEQ(val, oldval) )
      {
         SCIP_RATIONAL* obj;
         SCIP_CALL( solSetArrayValExact(sol, set, var, val) );
         obj = SCIPvarGetObjExact(var);
         SCIPrationalDiffProd(sol->valsexact->obj, obj, oldval);
         SCIPrationalAddProd(sol->valsexact->obj, obj, val);
      }

      SCIPrationalFreeBuffer(set->buffer, &oldval);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);
      if( !SCIPrationalIsEQ(val, SCIPvarGetLbGlobalExact(var)) )
      {
         SCIPerrorMessage("cannot set solution value for variable <%s> fixed to %.15g to different value %.15g\n",
            SCIPvarGetName(var), SCIPrationalGetReal(SCIPvarGetLbGlobalExact(var)), SCIPrationalGetReal(val));
         return SCIP_INVALIDDATA;
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPrationalIsZero(SCIPvarGetAggrScalarExact(var)));
      assert(!SCIPrationalIsAbsInfinity(SCIPvarGetAggrConstantExact(var)));
      assert(!SCIPrationalIsAbsInfinity(SCIPvarGetAggrScalarExact(var)));

      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );

      if( SCIPrationalIsAbsInfinity(val) )
      {
         if( !SCIPrationalIsPositive(SCIPvarGetAggrScalarExact(var)) )
            SCIPrationalNegate(tmp, val);
         retcode = SCIPsolSetValExact(sol, set, stat, tree, SCIPvarGetAggrVar(var),  tmp);
      }
      else
      {
         SCIPrationalDiff(tmp, val, SCIPvarGetAggrConstantExact(var));
         SCIPrationalDiv(tmp, tmp, SCIPvarGetAggrScalarExact(var));
         retcode = SCIPsolSetValExact(sol, set, stat, tree, SCIPvarGetAggrVar(var), tmp);
      }

      SCIPrationalFreeBuffer(set->buffer, &tmp);
      return retcode;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot set solution value for multiple aggregated variable\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      assert(!SCIPsetIsInfinity(set, SCIPvarGetNegationConstant(var)) && !SCIPsetIsInfinity(set, -SCIPvarGetNegationConstant(var)));

      return SCIPsolSetValExact(sol, set, stat, tree, SCIPvarGetNegationVar(var), val);

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** increases value of variable in primal CIP solution */
SCIP_RETCODE SCIPsolIncVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR*             var,                /**< variable to increase solution value for */
   SCIP_Real             incval              /**< increment for solution value of variable */
   )
{
   SCIP_Real oldval;

   assert(sol != NULL);
   assert(stat != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL
      || sol->solorigin == SCIP_SOLORIGIN_ZERO
      || (sol->nodenum == stat->nnodes && sol->runnum == stat->nruns));
   assert(var != NULL);
   assert(!SCIPsetIsInfinity(set, incval) && !SCIPsetIsInfinity(set, -incval));

   SCIPsetDebugMsg(set, "increasing value of <%s> in solution %p by %g\n", SCIPvarGetName(var), (void*)sol, incval);

   if( incval == 0.0 )
      return SCIP_OKAY;

   assert(sol->solorigin != SCIP_SOLORIGIN_LPSOL || SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var))
      || sol->lpcount == stat->lpcount);

   oldval = solGetArrayVal(sol, var);
   if( SCIPsetIsInfinity(set, oldval) || SCIPsetIsInfinity(set, -oldval) )
      return SCIP_OKAY;

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   /* @todo: handle strange cases, such as sums that yield infinite values */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( SCIPsolIsOriginal(sol) )
      {
         SCIP_CALL( solIncArrayVal(sol, set, var, incval) );
         sol->obj += SCIPvarGetUnchangedObj(var) * incval;
         solStamp(sol, stat, tree, FALSE);
         return SCIP_OKAY;
      }
      else
         return SCIPsolIncVal(sol, set, stat, tree, SCIPvarGetTransVar(var), incval);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(!SCIPsolIsOriginal(sol));
      SCIP_CALL( solIncArrayVal(sol, set, var, incval) );
      sol->obj += SCIPvarGetUnchangedObj(var) * incval;
      solStamp(sol, stat, tree, FALSE);
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      SCIPerrorMessage("cannot increase solution value for fixed variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      assert(!SCIPsetIsZero(set, SCIPvarGetAggrScalar(var)));
      return SCIPsolIncVal(sol, set, stat, tree, SCIPvarGetAggrVar(var), incval/SCIPvarGetAggrScalar(var));

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot increase solution value for multiple aggregated variable\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      return SCIPsolIncVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), -incval);

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** returns value of variable in primal CIP solution */
SCIP_Real SCIPsolGetVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* scalars;
   SCIP_Real solval;
   SCIP_Real solvalsum;
   int nvars;
   int i;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL
      || sol->solorigin == SCIP_SOLORIGIN_ZERO
      || sol->solorigin == SCIP_SOLORIGIN_PARTIAL
      || sol->solorigin == SCIP_SOLORIGIN_UNKNOWN
      || (sol->nodenum == stat->nnodes && sol->runnum == stat->nruns));
   assert(var != NULL);
   assert(var->scip == sol->scip);

   /* if the value of a transformed variable in an original solution is requested, we need to project the variable back
    * to the original space, the opposite case is handled below
    */
   if( SCIPsolIsOriginal(sol) && SCIPvarIsTransformed(var) )
   {
      SCIP_RETCODE retcode;
      SCIP_VAR* origvar;
      SCIP_Real scalar;
      SCIP_Real constant;

      /* we cannot get the value of a transformed variable for a solution that lives in the original problem space
       * -> get the corresponding original variable first
       */
      origvar = var;
      scalar = 1.0;
      constant = 0.0;
      retcode = SCIPvarGetOrigvarSum(&origvar, &scalar, &constant);
      if ( retcode != SCIP_OKAY )
         return SCIP_INVALID;
      if( origvar == NULL )
      {
         /* the variable has no original counterpart: in the original solution, it has a value of zero */
         return 0.0;
      }
      assert(!SCIPvarIsTransformed(origvar));

      solval = SCIPsolGetVal(sol, set, stat, origvar);
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         return SCIP_UNKNOWN;
      else
         return scalar * solval + constant;
   }

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( SCIPsolIsOriginal(sol) )
         return solGetArrayVal(sol, var);
      else
         return SCIPsolGetVal(sol, set, stat, SCIPvarGetTransVar(var));

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(!SCIPsolIsOriginal(sol));
      assert(sol->solorigin != SCIP_SOLORIGIN_LPSOL || SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var))
         || sol->lpcount == stat->lpcount);
      return solGetArrayVal(sol, var);

   case SCIP_VARSTATUS_FIXED:
      assert(!SCIPsolIsOriginal(sol));
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetUbGlobal(var) || (set->exact_enable && SCIPrationalIsEQ(SCIPvarGetLbGlobalExact(var), SCIPvarGetUbGlobalExact(var)))) ; /*lint !e777*/
      assert(SCIPvarGetLbLocal(var) == SCIPvarGetUbLocal(var) || (set->exact_enable && SCIPrationalIsEQ(SCIPvarGetLbLocalExact(var), SCIPvarGetUbLocalExact(var)))); /*lint !e777*/
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetLbLocal(var) || (set->exact_enable && SCIPrationalIsEQ(SCIPvarGetLbGlobalExact(var), SCIPvarGetLbLocalExact(var)))); /*lint !e777*/
      return SCIPvarGetLbGlobal(var);

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      solval = SCIPsolGetVal(sol, set, stat, SCIPvarGetAggrVar(var));
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         return SCIP_UNKNOWN;
      if( SCIPsetIsInfinity(set, solval) || SCIPsetIsInfinity(set, -solval) )
      {
         if( SCIPvarGetAggrScalar(var) * solval > 0.0 )
            return SCIPsetInfinity(set);
         if( SCIPvarGetAggrScalar(var) * solval < 0.0 )
            return -SCIPsetInfinity(set);
      }
      return SCIPvarGetAggrScalar(var) * solval + SCIPvarGetAggrConstant(var);

   case SCIP_VARSTATUS_MULTAGGR:
      nvars = SCIPvarGetMultaggrNVars(var);
      vars = SCIPvarGetMultaggrVars(var);
      scalars = SCIPvarGetMultaggrScalars(var);
      solvalsum = SCIPvarGetMultaggrConstant(var);
      for( i = 0; i < nvars; ++i )
      {
         solval = SCIPsolGetVal(sol, set, stat, vars[i]);
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            return SCIP_UNKNOWN;
         if( SCIPsetIsInfinity(set, solval) || SCIPsetIsInfinity(set, -solval) )
         {
            if( scalars[i] * solval > 0.0 )
               return SCIPsetInfinity(set);
            if( scalars[i] * solval < 0.0 )
               return -SCIPsetInfinity(set);
         }
         solvalsum += scalars[i] * solval;
      }
      return solvalsum;

   case SCIP_VARSTATUS_NEGATED:
      solval = SCIPsolGetVal(sol, set, stat, SCIPvarGetNegationVar(var));
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         return SCIP_UNKNOWN;
      if( SCIPsetIsInfinity(set, solval) )
         return -SCIPsetInfinity(set);
      if( SCIPsetIsInfinity(set, -solval) )
         return SCIPsetInfinity(set);
      return SCIPvarGetNegationConstant(var) - solval;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns value of variable in exact primal CIP solution */
void SCIPsolGetValExact(
   SCIP_RATIONAL*        res,                /**< resulting rational */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   SCIP_VAR** vars;
   SCIP_RATIONAL** scalars;
   SCIP_RATIONAL* solval;
   int nvars;
   int i;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL
      || sol->solorigin == SCIP_SOLORIGIN_ZERO
      || sol->solorigin == SCIP_SOLORIGIN_PARTIAL
      || sol->solorigin == SCIP_SOLORIGIN_LPSOL
      || sol->solorigin == SCIP_SOLORIGIN_UNKNOWN);
   assert(var != NULL);
   assert(SCIPsolIsExact(sol));

   /* if the value of a transformed variable in an original solution is requested, we need to project the variable back
    * to the original space, the opposite case is handled below
    */
   if( sol->solorigin == SCIP_SOLORIGIN_ORIGINAL && SCIPvarIsTransformed(var) )
   {
      SCIP_RETCODE retcode;
      SCIP_VAR* origvar;
      SCIP_RATIONAL* scalar;
      SCIP_RATIONAL* constant;

      (void) SCIPrationalCreateBuffer(set->buffer, &scalar);
      (void) SCIPrationalCreateBuffer(set->buffer, &constant);

      /* we cannot get the value of a transformed variable for a solution that lives in the original problem space
       * -> get the corresponding original variable first
       */
      origvar = var;
      SCIPrationalSetFraction(scalar, 1LL, 1LL);
      SCIPrationalSetReal(constant, 0.0);
      retcode = SCIPvarGetOrigvarSumExact(&origvar, scalar, constant);
      if ( retcode != SCIP_OKAY )
      {
         SCIPABORT();
         return;
      }
      if( origvar == NULL )
      {
         /* the variable has no original counterpart: in the original solution, it has a value of zero */
         SCIPrationalSetReal(res, 0.0);
         return;
      }

      assert(!SCIPvarIsTransformed(origvar));

      SCIPsolGetValExact(res, sol, set, stat, origvar);
      SCIPrationalMult(res, res, scalar);
      SCIPrationalAdd(res, res, constant);

      SCIPrationalFreeBuffer(set->buffer, &constant);
      SCIPrationalFreeBuffer(set->buffer, &scalar);

      return;
   }

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatusExact(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( sol->solorigin == SCIP_SOLORIGIN_ORIGINAL )
         solGetArrayValExact(res, sol, var);
      else
         SCIPsolGetValExact(res, sol, set, stat, SCIPvarGetTransVar(var));
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);
      assert(sol->solorigin != SCIP_SOLORIGIN_LPSOL || SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var))
      || sol->lpcount == stat->lpcount );
      solGetArrayValExact(res, sol, var);
      break;

   case SCIP_VARSTATUS_FIXED:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);
      assert(SCIPrationalIsEQ(SCIPvarGetLbGlobalExact(var), SCIPvarGetUbGlobalExact(var))); /*lint !e777*/
      assert(SCIPrationalIsEQ(SCIPvarGetLbLocalExact(var), SCIPvarGetUbLocalExact(var))); /*lint !e777*/
      assert(SCIPrationalIsEQ(SCIPvarGetLbGlobalExact(var), SCIPvarGetLbLocalExact(var))); /*lint !e777*/
      SCIPrationalSetRational(res, SCIPvarGetLbGlobalExact(var));
      break;

    case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      SCIPsolGetValExact(res, sol, set, stat, SCIPvarGetAggrVar(var));
      if( SCIPrationalIsAbsInfinity(res) )
      {
         if( SCIPrationalGetSign(res) * SCIPrationalGetSign(SCIPvarGetAggrScalarExact(var)) > 0 )
         {
            SCIPrationalSetInfinity(res);
            return;
         }
         if( SCIPrationalGetSign(res) * SCIPrationalGetSign(SCIPvarGetAggrScalarExact(var)) < 0 )
         {
            SCIPrationalSetNegInfinity(res);
            return;
         }
      }
      SCIPrationalMult(res, res, SCIPvarGetAggrScalarExact(var));
      SCIPrationalAdd(res, res, SCIPvarGetAggrConstantExact(var));
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      (void) SCIPrationalCreateBuffer(set->buffer, &solval);

      nvars = SCIPvarGetMultaggrNVars(var);
      vars = SCIPvarGetMultaggrVars(var);
      scalars = SCIPvarGetMultaggrScalarsExact(var);
      SCIPrationalSetRational(res, SCIPvarGetMultaggrConstantExact(var));
      for( i = 0; i < nvars; ++i )
      {
         SCIPsolGetValExact(solval, sol, set, stat, vars[i]);
         if( SCIPrationalIsAbsInfinity(solval) )
         {
            if( SCIPrationalGetSign(scalars[i]) == SCIPrationalGetSign(solval) )
               SCIPrationalSetInfinity(res);
            if( SCIPrationalGetSign(scalars[i]) != SCIPrationalGetSign(solval) && !SCIPrationalIsZero(scalars[i]) )
               SCIPrationalSetNegInfinity(res);
            break;
         }
         SCIPrationalAddProd(res, scalars[i], solval);
      }
      SCIPrationalFreeBuffer(set->buffer, &solval);
      break;

   case SCIP_VARSTATUS_NEGATED:
      SCIPsolGetValExact(res, sol, set, stat, SCIPvarGetNegationVar(var));
      SCIPrationalDiffReal(res, res, SCIPvarGetNegationConstant(var));
      SCIPrationalNegate(res, res);
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      SCIPrationalSetReal(res, 0.0); /*lint !e527*/
   }
}

/** returns value of variable in primal ray represented by primal CIP solution */
SCIP_Real SCIPsolGetRayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution, representing a primal ray */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* scalars;
   SCIP_Real solval;
   SCIP_Real solvalsum;
   int nvars;
   int i;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO);
   assert(var != NULL);

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      return SCIPsolGetRayVal(sol, set, stat, SCIPvarGetTransVar(var));

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return solGetArrayVal(sol, var);

   case SCIP_VARSTATUS_FIXED:
      assert(!SCIPsolIsOriginal(sol));
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetUbGlobal(var)); /*lint !e777*/
      assert(SCIPvarGetLbLocal(var) == SCIPvarGetUbLocal(var)); /*lint !e777*/
      assert(SCIPvarGetLbGlobal(var) == SCIPvarGetLbLocal(var)); /*lint !e777*/
      return 0.0; /* constants are ignored for computing the ray direction */

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      solval = SCIPsolGetRayVal(sol, set, stat, SCIPvarGetAggrVar(var));
      assert(solval != SCIP_UNKNOWN); /*lint !e777*/
      assert(!SCIPsetIsInfinity(set, REALABS(solval)));
      return SCIPvarGetAggrScalar(var) * solval; /* constants are ignored for computing the ray direction */

   case SCIP_VARSTATUS_MULTAGGR:
      nvars = SCIPvarGetMultaggrNVars(var);
      vars = SCIPvarGetMultaggrVars(var);
      scalars = SCIPvarGetMultaggrScalars(var);
      solvalsum = 0.0; /* constants are ignored for computing the ray direction */
      for( i = 0; i < nvars; ++i )
      {
         solval = SCIPsolGetRayVal(sol, set, stat, vars[i]);
         assert(solval != SCIP_UNKNOWN ); /*lint !e777*/
         assert(!SCIPsetIsInfinity(set, REALABS(solval)));
         solvalsum += scalars[i] * solval;
      }
      return solvalsum;

   case SCIP_VARSTATUS_NEGATED:
      solval = SCIPsolGetRayVal(sol, set, stat, SCIPvarGetNegationVar(var));
      assert(solval != SCIP_UNKNOWN); /*lint !e777*/
      assert(!SCIPsetIsInfinity(set, REALABS(solval)));
      return -solval; /* constants are ignored for computing the ray direction */

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** gets objective value of primal CIP solution in transformed problem */
SCIP_Real SCIPsolGetObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob            /**< original problem data */
   )
{
   assert(sol != NULL);

   /* for original solutions, sol->obj contains the external objective value */
   if( SCIPsolIsOriginal(sol) )
      return SCIPprobInternObjval(transprob, origprob, set, sol->obj);
   else
      return sol->obj;
}

/** gets objective value of exact primal CIP solution in transformed problem */
void SCIPsolGetObjExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_RATIONAL*        objval              /**< store the result here */
    )
{
   assert(sol != NULL);
   assert(SCIPsolIsExact(sol));

   /* for original solutions, sol->obj contains the external objective value */
   if( SCIPsolIsOriginal(sol) )
      SCIPprobInternObjvalExact(transprob, origprob, set, sol->valsexact->obj, objval);
   else
      SCIPrationalSetRational(objval, sol->valsexact->obj);
}

/** updates primal solutions after a change in a variable's objective value */
void SCIPsolUpdateVarObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldobj,             /**< old objective value */
   SCIP_Real             newobj              /**< new objective value */
   )
{
   SCIP_Real solval;

   assert(sol != NULL);
   assert(!SCIPsolIsOriginal(sol));
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   solval = solGetArrayVal(sol, var);
   if( solval != SCIP_UNKNOWN ) /*lint !e777*/
      sol->obj += (newobj - oldobj) * solval;
}

/* mark the given solution as partial solution */
SCIP_RETCODE SCIPsolMarkPartial(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars               /**< number of problem variables */
   )
{
   SCIP_Real* vals;
   int v;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL);
   assert(nvars == 0 || vars != NULL);

   if( nvars == 0 )
      return SCIP_OKAY;;

   SCIP_CALL( SCIPsetAllocBufferArray(set, &vals, nvars) );

   /* get values */
   for( v = 0; v < nvars; v++ )
   {
      assert(!SCIPvarIsTransformed(vars[v]));
      vals[v] = SCIPsolGetVal(sol, set, stat, vars[v]);
   }

   /* change origin to partial */
   sol->solorigin = SCIP_SOLORIGIN_PARTIAL;

   /* set values */
   for( v = 0; v < nvars; v++ )
   {
      int idx = SCIPvarGetIndex(vars[v]);

      if( vals[v] != SCIP_UNKNOWN ) /*lint !e777*/
      {
         /* from now on, variable must not be deleted */
         SCIPvarMarkNotDeletable(vars[v]);

         /* mark the variable valid */
         SCIP_CALL( SCIPboolarraySetVal(sol->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, TRUE) );

         /* set the value in the solution array */
         SCIP_CALL( SCIPrealarraySetVal(sol->vals, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, vals[v]) );
      }
      else
      {
         /* mark the variable invalid */
         SCIP_CALL( SCIPboolarraySetVal(sol->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, FALSE) );
      }
   }

   /* free buffer */
   SCIPsetFreeBufferArray(set, &vals);

   return SCIP_OKAY;
}

/** checks primal CIP solution for exact feasibility
 *  (either checks fp values exactly or rational values if it is a rational solution)
 */
static
SCIP_RETCODE solCheckExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            feasible            /**< stores whether solution is feasible */
   )
{
   SCIP_RESULT result;
   SCIP_RATIONAL* solval;
   int h;

   assert(sol != NULL);
   assert(!SCIPsolIsOriginal(sol));
   assert(set != NULL);
   assert(prob != NULL);
   assert(feasible != NULL);

   SCIPsetDebugMsg(set, "checking solution with objective value %g (nodenum=%" SCIP_LONGINT_FORMAT ", origin=%u)\n",
      sol->obj, sol->nodenum, sol->solorigin);

   *feasible = TRUE;

   if( !printreason )
      completely = FALSE;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &solval) );

   /* check whether the solution respects the global bounds of the variables */
   {
      int v;

      for( v = 0; v < prob->nvars && (*feasible || completely); ++v )
      {
         SCIP_VAR* var;

         var = prob->vars[v];
         if( SCIPsolIsExact(sol) )
            SCIPsolGetValExact(solval, sol, set, stat, var);
         else
            SCIPrationalSetReal(solval, SCIPsolGetVal(sol, set, stat, var));

         if( !SCIPrationalIsAbsInfinity(solval) ) /*lint !e777*/
         {
            SCIP_RATIONAL* lb;
            SCIP_RATIONAL* ub;

            lb = SCIPvarGetLbGlobalExact(var);
            ub = SCIPvarGetUbGlobalExact(var);

            /* if we have to check bound and one of the current bounds is violated */
            if( (!SCIPrationalIsNegInfinity(lb) && SCIPrationalIsLT(solval, lb)) || (!SCIPrationalIsInfinity(ub) && SCIPrationalIsGT(solval, ub)) )
            {
               *feasible = FALSE;

               if( printreason )
               {
                  SCIPmessagePrintInfo(messagehdlr, "solution value %g violates bounds of <%s>[%g,%g] by %g\n", SCIPrationalGetReal(solval), SCIPvarGetName(var),
                     SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), SCIPrationalIsGT(solval, ub) ? SCIPrationalGetReal(lb) - SCIPrationalGetReal(solval) : SCIPrationalGetReal(solval) - SCIPrationalGetReal(ub));
               }
#ifdef SCIP_DEBUG
               else
               {
                  SCIPrationalDebugMessage("  -> solution value %q violates bounds of <%s>[%g,%g]\n", solval, SCIPvarGetName(var),
                     SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
               }
#endif
            }

            /* check whether there are infinite variable values that lead to an objective value of +infinity */
            if( *feasible && sol->hasinfval )
            {
               *feasible = *feasible && (!SCIPrationalIsInfinity(solval) || !SCIPrationalIsPositive(SCIPvarGetObjExact(var)) );
               *feasible = *feasible && (!SCIPrationalIsNegInfinity(solval) || !SCIPrationalIsNegative(SCIPvarGetObjExact(var)) );

               if( ((SCIPrationalIsInfinity(solval) && SCIPrationalIsPositive(SCIPvarGetObjExact(var))))
                  || (SCIPrationalIsNegInfinity(solval) && SCIPrationalIsNegative(SCIPvarGetObjExact(var))) )
               {
                  if( printreason )
                  {
                     SCIPrationalDebugMessage("infinite solution value %q for variable  <%s> with obj %q implies objective value +infinity\n",
                        SCIPrationalGetReal(solval), SCIPvarGetName(var), SCIPvarGetUnchangedObj(var));
                  }
#ifdef SCIP_DEBUG
                  else
                  {
                     SCIPrationalDebugMessage("infinite solution value %q for variable  <%s> with obj %g implies objective value +infinity\n",
                        solval, SCIPvarGetName(var), SCIPvarGetUnchangedObj(var));
                  }
#endif
               }
            }
         }
      }
   }

   /* check whether the solution fulfills all constraints */
   for( h = 0; h < set->nconshdlrs && (*feasible || completely); ++h )
   {
      SCIP_CALL( SCIPconshdlrCheck(set->conshdlrs[h], blkmem, set, stat, sol,
            TRUE, checklprows, printreason, completely, &result) );
      *feasible = *feasible && (result == SCIP_FEASIBLE);

#ifdef SCIP_DEBUG
      if( !(*feasible) )
      {
         SCIPdebugPrintf("  -> infeasibility detected in constraint handler <%s>\n",
            SCIPconshdlrGetName(set->conshdlrs[h]));
      }
#endif
   }

   SCIPrationalFreeBuffer(set->buffer, &solval);

   return SCIP_OKAY;
}

/** checks solution for feasibility in original problem without adding it to the solution store
 *
 *  We first check the variable bounds. Then we loop over all constraint handlers and constraints, checking each in the
 *  order of their check priority.
 */
SCIP_RETCODE SCIPsolCheckOrig(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked if printreason is true? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             checkmodifiable,    /**< have modifiable constraint to be checked? */
   SCIP_Bool*            feasible            /**< stores whether given solution is feasible */
   )
{
   SCIP_RESULT result;
#ifndef NDEBUG
   int oldpriority;
#endif
   int v;
   int c;
   int h;

   assert(sol != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert(!prob->transformed);
   assert(feasible != NULL);

   *feasible = TRUE;

   SCIPsolResetViolations(sol);

   if( !printreason )
      completely = FALSE;

   /* check bounds */
   if( checkbounds )
   {
      for( v = 0; v < prob->nvars; ++v )
      {
         SCIP_VAR* var;
         SCIP_Real solval;
         SCIP_Real lb;
         SCIP_Real ub;

         var = prob->vars[v];
         solval = SCIPsolGetVal(sol, set, stat, var);

         lb = SCIPvarGetLbOriginal(var);
         ub = SCIPvarGetUbOriginal(var);

         if( SCIPprimalUpdateViolations(primal) )
         {
            SCIPsolUpdateBoundViolation(sol, lb - solval, SCIPrelDiff(lb, solval));
            SCIPsolUpdateBoundViolation(sol, solval - ub, SCIPrelDiff(solval, ub));
         }

         if( SCIPsetIsFeasLT(set, solval, lb) || SCIPsetIsFeasGT(set, solval, ub) )
         {
            *feasible = FALSE;

            if( printreason )
            {
               SCIPmessagePrintInfo(messagehdlr, "solution violates original bounds of variable <%s> [%g,%g] solution value <%g>\n",
                  SCIPvarGetName(var), lb, ub, solval);
            }

            if( !completely )
               return SCIP_OKAY;
         }
      }
   }

   /* sort original constraint according to check priority */
   SCIP_CALL( SCIPprobSortConssCheck(prob) );

   /* check original constraints
    *
    * in general modifiable constraints can not be checked, because the variables to fulfill them might be missing in
    * the original problem; however, if the solution comes from a heuristic during presolving, modifiable constraints
    * have to be checked;
    */
#ifndef NDEBUG
   oldpriority = INT_MAX;
#endif
   h = 0;
   for( c = 0; c < prob->nconss; ++c )
   {
      SCIP_CONS* cons;
      int priority;

      cons = prob->origcheckconss[c];
      assert( SCIPconsGetHdlr(cons) != NULL );
      priority = SCIPconshdlrGetCheckPriority(SCIPconsGetHdlr(cons));

#ifndef NDEBUG
      assert( priority <= oldpriority );
      oldpriority = priority;
#endif

      /* check constraints handlers without constraints that have a check priority at least as high as current
       * constraint */
      while( h < set->nconshdlrs && SCIPconshdlrGetCheckPriority(set->conshdlrs[h]) >= priority )
      {
         if( !SCIPconshdlrNeedsCons(set->conshdlrs[h]) )
         {
            SCIP_CALL( SCIPconshdlrCheck(set->conshdlrs[h], blkmem, set, stat, sol,
                  checkintegrality, checklprows, printreason, completely, &result) );

            if( result != SCIP_FEASIBLE )
            {
               *feasible = FALSE;

               if( !completely )
                  return SCIP_OKAY;
            }
         }
         ++h;
      }

      /* now check constraint */
      if( SCIPconsIsChecked(cons) && (checkmodifiable || !SCIPconsIsModifiable(cons)) )
      {
         /* check solution */
         SCIP_CALL( SCIPconsCheck(cons, set, sol, checkintegrality, checklprows, printreason, &result) );

         if( result != SCIP_FEASIBLE )
         {
            *feasible = FALSE;

            if( !completely )
               return SCIP_OKAY;
         }
      }
   }

   /* one final loop over the remaining constraints handlers without constraints */
   while( h < set->nconshdlrs )
   {
      if( !SCIPconshdlrNeedsCons(set->conshdlrs[h]) )
      {
         SCIP_CALL( SCIPconshdlrCheck(set->conshdlrs[h], blkmem, set, stat, sol,
               checkintegrality, checklprows, printreason, completely, &result) );

         if( result != SCIP_FEASIBLE )
         {
            *feasible = FALSE;

            if( !completely )
               return SCIP_OKAY;
         }
      }
      ++h;
   }

   return SCIP_OKAY;
}

/** checks primal CIP solution for feasibility
 *
 *  @note The difference between SCIPsolCheck() and SCIPcheckSolOrig() is that modifiable constraints are handled
 *        differently. There might be some variables which do not have an original counter part (e.g. in
 *        branch-and-price). Therefore, modifiable constraints can not be double-checked in the original space.
 */
SCIP_RETCODE SCIPsolCheck(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            feasible            /**< stores whether solution is feasible */
   )
{
   SCIP_RESULT result;
   int h;

   assert(sol != NULL);
   assert(!SCIPsolIsOriginal(sol));
   assert(set != NULL);
   assert(prob != NULL);
   assert(feasible != NULL);

   SCIPsetDebugMsg(set, "checking solution with objective value %g (nodenum=%" SCIP_LONGINT_FORMAT ", origin=%d)\n",
      sol->obj, sol->nodenum, sol->solorigin);

   *feasible = TRUE;

   /* have to check bounds without tolerances in exact solving mode */
   if( set->exact_enable )
   {
      SCIP_CALL( solCheckExact(sol, set, messagehdlr, blkmem, stat, prob, printreason,
            completely, checklprows, feasible) );
   }

   SCIPsolResetViolations(sol);

   if( !printreason )
      completely = FALSE;

   /* check whether the solution respects the global bounds of the variables */
   if( checkbounds || sol->hasinfval )
   {
      int v;

      for( v = 0; v < prob->nvars && (*feasible || completely); ++v )
      {
         SCIP_VAR* var;
         SCIP_Real solval;

         var = prob->vars[v];
         solval = SCIPsolGetVal(sol, set, stat, var);

         if( solval != SCIP_UNKNOWN ) /*lint !e777*/
         {
            SCIP_Real lb;
            SCIP_Real ub;

            lb = SCIPvarGetLbGlobal(var);
            ub = SCIPvarGetUbGlobal(var);

            /* if we have to check bound and one of the current bounds is violated */
            if( checkbounds && ((!SCIPsetIsInfinity(set, -lb) && SCIPsetIsFeasLT(set, solval, lb))
                  || (!SCIPsetIsInfinity(set, ub) && SCIPsetIsFeasGT(set, solval, ub))) )
            {
               *feasible = FALSE;

               if( printreason )
               {
                  SCIPmessagePrintInfo(messagehdlr, "solution value %g violates bounds of <%s>[%g,%g] by %g\n", solval, SCIPvarGetName(var),
                     SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), MAX(lb - solval, 0.0) + MAX(solval - ub, 0.0));
               }
#ifdef SCIP_DEBUG
               else
               {
                  SCIPsetDebugMsgPrint(set, "  -> solution value %g violates bounds of <%s>[%g,%g]\n", solval, SCIPvarGetName(var),
                     SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
               }
#endif
            }

            /* check whether there are infinite variable values that lead to an objective value of +infinity */
            if( *feasible && sol->hasinfval )
            {
               *feasible = *feasible && (!SCIPsetIsInfinity(set, solval) || SCIPsetIsLE(set, SCIPvarGetUnchangedObj(var), 0.0) );
               *feasible = *feasible && (!SCIPsetIsInfinity(set, -solval) || SCIPsetIsGE(set, SCIPvarGetUnchangedObj(var), 0.0) );

               if( ((SCIPsetIsInfinity(set, solval) && SCIPsetIsGT(set, SCIPvarGetUnchangedObj(var), 0.0)) || (SCIPsetIsInfinity(set, -solval) && SCIPsetIsLT(set, SCIPvarGetUnchangedObj(var), 0.0))) )
               {
                  if( printreason )
                  {
                     SCIPmessagePrintInfo(messagehdlr, "infinite solution value %g for variable  <%s> with obj %g implies objective value +infinity\n",
                        solval, SCIPvarGetName(var), SCIPvarGetUnchangedObj(var));
                  }
#ifdef SCIP_DEBUG
                  else
                  {
                     SCIPsetDebugMsgPrint(set, "infinite solution value %g for variable  <%s> with obj %g implies objective value +infinity\n",
                        solval, SCIPvarGetName(var), SCIPvarGetUnchangedObj(var));
                  }
#endif
               }
            }
         }
      }
   }

   /* check whether the solution fulfills all constraints */
   for( h = 0; h < set->nconshdlrs && (*feasible || completely); ++h )
   {
      SCIP_CALL( SCIPconshdlrCheck(set->conshdlrs[h], blkmem, set, stat, sol,
            checkintegrality, checklprows, printreason, completely, &result) );
      *feasible = *feasible && (result == SCIP_FEASIBLE);

#ifdef SCIP_DEBUG
      if( !(*feasible) )
      {
         SCIPdebugPrintf("  -> infeasibility detected in constraint handler <%s>\n",
            SCIPconshdlrGetName(set->conshdlrs[h]));
      }
#endif
   }

   return SCIP_OKAY;
}

/** try to round given solution */
SCIP_RETCODE SCIPsolRound(
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool*            success             /**< pointer to store whether rounding was successful */
   )
{
   int nvars;
   int v;

   assert(sol != NULL);
   assert(!SCIPsolIsOriginal(sol));
   assert(prob != NULL);
   assert(prob->transformed);
   assert(success != NULL);

   /* round all roundable fractional enforced integral variables as long as no unroundable var was found */
   nvars = prob->nvars - prob->ncontvars - prob->ncontimplvars;
   assert(nvars >= 0);
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real solval;
      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;

      var = prob->vars[v];
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      assert(sol->solorigin != SCIP_SOLORIGIN_LPSOL || SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var))
         || sol->lpcount == stat->lpcount);
      solval = solGetArrayVal(sol, var);

      /* solutions with unknown entries cannot be rounded */
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         break;

      /* if solution value is already integral with feastol, continue */
      if( SCIPsetIsFeasIntegral(set, solval) )
         continue;

      /* get rounding possibilities */
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);

      /* choose rounding direction */
      if( mayrounddown && mayroundup )
      {
         /* we can round in both directions: round in objective function direction */
         if( SCIPvarGetUnchangedObj(var) >= 0.0 )
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
      SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, var, solval) );
   }

   /* check, if rounding was successful */
   *success = (v == nvars);

   return SCIP_OKAY;
}

/** copies the real values to the exact arrays of the solution */
SCIP_RETCODE SCIPsolMakeExact(
   SCIP_SOL*             sol,                /**< primal solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob                /**< transformed problem data */
   )
{
   SCIP_RATIONAL* tmp;
   int v;

   if( SCIPsolIsExact(sol) )
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &sol->valsexact) );
   SCIP_CALL( SCIPrationalarrayCreate(&sol->valsexact->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&sol->valsexact->valid, blkmem) );
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &sol->valsexact->obj) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );

   SCIP_CALL( SCIPsolUnlink(sol, set, prob) );

   for( v = 0; v < prob->nvars; ++v )
   {
      SCIPrationalSetReal(tmp, solGetArrayVal(sol, prob->vars[v]));
      SCIP_CALL( solSetArrayValExact(sol, set, prob->vars[v], tmp) );
   }

   SCIPsolRecomputeInternObjExact(sol, set, stat, prob);

   SCIPrationalFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** approximates and copies the exact values to the real arrays of the solution and frees the exact data */
SCIP_RETCODE SCIPsolMakeReal(
   SCIP_SOL*             sol,                /**< primal solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob                /**< transformed problem data */
   )
{
   SCIP_RATIONAL* tmp;
   int v;

   if( !SCIPsolIsExact(sol) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );

   SCIP_CALL( SCIPsolUnlinkExact(sol, set, prob) );

   for( v = 0; v < prob->nvars; ++v )
   {
      solGetArrayValExact(tmp, sol, prob->vars[v]);
      SCIP_CALL( solSetArrayVal(sol, set, prob->vars[v], SCIPrationalGetReal(tmp)) );
   }

   SCIPsolRecomputeObj(sol, set, stat, prob);

   SCIPrationalFreeBuffer(set->buffer, &tmp);
   SCIPrationalFreeBlock(blkmem, &sol->valsexact->obj);
   SCIP_CALL( SCIPboolarrayFree(&sol->valsexact->valid) );
   SCIP_CALL( SCIPrationalarrayFree(&sol->valsexact->vals, blkmem) );
   BMSfreeBlockMemory(blkmem, &sol->valsexact);

   return SCIP_OKAY;
}

/** updates the solution value sums in variables by adding the value in the given solution */
void SCIPsolUpdateVarsum(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Real             weight              /**< weight of solution in weighted average */
   )
{
   SCIP_Real solval;
   int v;

   assert(sol != NULL);
   assert(!SCIPsolIsOriginal(sol));
   assert(0.0 <= weight && weight <= 1.0);

   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);
      solval = SCIPsolGetVal(sol, set, stat, prob->vars[v]);
      if( solval != SCIP_UNKNOWN ) /*lint !e777*/
      {
         prob->vars[v]->primsolavg *= (1.0-weight);
         prob->vars[v]->primsolavg += weight*solval;
      }
   }
}

/** retransforms solution to original problem space */
SCIP_RETCODE SCIPsolRetransform(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_Bool*            hasinfval           /**< pointer to store whether the solution has infinite values */
   )
{
   SCIP_VAR** transvars;
   SCIP_VAR** vars;
   SCIP_VAR** activevars;
   SCIP_Real* solvals;
   SCIP_Real* activevals;
   SCIP_Real* transsolvals;
   SCIP_Real constant;
   int requiredsize;
   int ntransvars;
   int nactivevars;
   int nvars;
   int v;
   int i;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO);
   assert(origprob != NULL);
   assert(transprob != NULL);
   assert(hasinfval != NULL);
   assert(!origprob->transformed);
   assert(transprob->transformed);

   *hasinfval = FALSE;

   /* transform exact values first (needs unchanged solorigin) */
   if( SCIPsolIsExact(sol) )
   {
      SCIP_CALL( SCIPsolRetransformExact(sol, set, stat, origprob, transprob, hasinfval) );
      SCIP_CALL( SCIPsolOverwriteFPSolWithExact(sol, set, stat, origprob, transprob, NULL) );
      return SCIP_OKAY;
   }

   /* This method was a performance bottleneck when retransforming a solution during presolving, before flattening the
    * aggregation graph. In that case, calling SCIPsolGetVal() on the original variable consumed too much
    * time. Therefore, we now first compute the active representation of each original variable using
    * SCIPvarGetActiveRepresentatives(), which is much faster, and sum up the solution values of the active variables by
    * hand for each original variable.
    */
   vars = origprob->vars;
   nvars = origprob->nvars;
   transvars = transprob->vars;
   ntransvars = transprob->nvars;

   /* allocate temporary memory for getting the active representation of the original variables, buffering the solution
    * values of all active variables and storing the original solution values
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &transsolvals, ntransvars + 1) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activevars, ntransvars + 1) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activevals, ntransvars + 1) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &solvals, nvars) );
   assert(transsolvals != NULL); /* for flexelint */
   assert(solvals != NULL); /* for flexelint */

   /* get the solution values of all active variables */
   for( v = 0; v < ntransvars; ++v )
   {
      transsolvals[v] = SCIPsolGetVal(sol, set, stat, transvars[v]);
   }

   /* get the solution in original problem variables */
   for( v = 0; v < nvars; ++v )
   {
      activevars[0] = vars[v];
      activevals[0] = 1.0;
      nactivevars = 1;
      constant = 0.0;

      /* get active representation of the original variable */
      SCIP_CALL( SCIPvarGetActiveRepresentatives(set, activevars, activevals, &nactivevars, ntransvars + 1, &constant,
            &requiredsize) );
      assert(requiredsize <= ntransvars);

      /* compute solution value of the original variable */
      solvals[v] = constant;
      for( i = 0; i < nactivevars; ++i )
      {
         assert(0 <= SCIPvarGetProbindex(activevars[i]) && SCIPvarGetProbindex(activevars[i]) < ntransvars);
         assert(!SCIPsetIsInfinity(set, -solvals[v]) || !SCIPsetIsInfinity(set, activevals[i] * transsolvals[SCIPvarGetProbindex(activevars[i])]));
         assert(!SCIPsetIsInfinity(set, solvals[v]) || !SCIPsetIsInfinity(set, -activevals[i] * transsolvals[SCIPvarGetProbindex(activevars[i])]));
         solvals[v] += activevals[i] * transsolvals[SCIPvarGetProbindex(activevars[i])];
      }

      if( SCIPsetIsInfinity(set, solvals[v]) )
      {
         solvals[v] = SCIPsetInfinity(set);
         *hasinfval = TRUE;
      }
      else if( SCIPsetIsInfinity(set, -solvals[v]) )
      {
         solvals[v] = -SCIPsetInfinity(set);
         *hasinfval = TRUE;
      }
   }

   /* clear the solution and convert it into original space */
   SCIP_CALL( solClearArrays(sol) );
   sol->solorigin = SCIP_SOLORIGIN_ORIGINAL;
   sol->obj = origprob->objoffset;

   /* reinsert the values of the original variables */
   for( v = 0; v < nvars; ++v )
   {
      assert(SCIPvarGetUnchangedObj(vars[v]) == SCIPvarGetObj(vars[v])); /*lint !e777*/

      if( solvals[v] != 0.0 )
      {
         SCIP_CALL( solSetArrayVal(sol, set, vars[v], solvals[v]) );
         if( solvals[v] != SCIP_UNKNOWN ) /*lint !e777*/
            sol->obj += SCIPvarGetUnchangedObj(vars[v]) * solvals[v];
      }
   }

   /**@todo remember the variables without original counterpart (priced variables) in the solution */

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &solvals);
   SCIPsetFreeBufferArray(set, &activevals);
   SCIPsetFreeBufferArray(set, &activevars);
   SCIPsetFreeBufferArray(set, &transsolvals);

   return SCIP_OKAY;
}

/** retransforms exact solution to original problem space */
SCIP_RETCODE SCIPsolRetransformExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_Bool*            hasinfval           /**< pointer to store whether the solution has infinite values */
   )
{
   SCIP_VAR** transvars;
   SCIP_VAR** vars;
   SCIP_VAR** activevars;
   SCIP_RATIONAL** solvals;
   SCIP_RATIONAL** activevals;
   SCIP_RATIONAL** transsolvals;
   SCIP_RATIONAL* constant;
   int requiredsize;
   int ntransvars;
   int nactivevars;
   int nvars;
   int v;
   int i;

   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ZERO);
   assert(origprob != NULL);
   assert(transprob != NULL);
   assert(hasinfval != NULL);
   assert(!origprob->transformed);
   assert(transprob->transformed);
   assert(SCIPsolIsExact(sol));

   *hasinfval = FALSE;

   /* This method was a performance bottleneck when retransforming a solution during presolving, before flattening the
    * aggregation graph. In that case, calling SCIPsolGetVal() on the original variable consumed too much
    * time. Therefore, we now first compute the active representation of each original variable using
    * SCIPvarGetActiveRepresentatives(), which is much faster, and sum up the solution values of the active variables by
    * hand for each original variable.
    */
   vars = origprob->vars;
   nvars = origprob->nvars;
   transvars = transprob->vars;
   ntransvars = transprob->nvars;

   /* allocate temporary memory for getting the active representation of the original variables, buffering the solution
    * values of all active variables and storing the original solution values
    */
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &transsolvals, ntransvars + 1) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activevars, ntransvars + 1) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &activevals, ntransvars + 1) );
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &solvals, nvars) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &constant) );

   assert(transsolvals != NULL); /* for flexelint */

   /* get the solution values of all active variables */
   for( v = 0; v < ntransvars; ++v )
   {
      SCIPsolGetValExact(transsolvals[v], sol, set, stat, transvars[v]);
   }

   /* get the solution in original problem variables */
   for( v = 0; v < nvars; ++v )
   {
      activevars[0] = vars[v];
      SCIPrationalSetReal(activevals[0], 1.0);
      nactivevars = 1;
      SCIPrationalSetReal(constant, 0.0);

      /* get active representation of the original variable */
      SCIP_CALL( SCIPvarGetActiveRepresentativesExact(set, activevars, activevals, &nactivevars, ntransvars + 1, constant,
            &requiredsize, TRUE) );
      assert(requiredsize <= ntransvars);

      /* compute solution value of the original variable */
      SCIPrationalSetRational(solvals[v], constant);
      for( i = 0; i < nactivevars; ++i )
      {
         assert(0 <= SCIPvarGetProbindex(activevars[i]) && SCIPvarGetProbindex(activevars[i]) < ntransvars);
         SCIPrationalAddProd(solvals[v], activevals[i], transsolvals[SCIPvarGetProbindex(activevars[i])]);
      }

      if( SCIPrationalIsAbsInfinity(solvals[v]) )
         *hasinfval = TRUE;
   }

   /* clear the solution and convert it into original space */
   SCIP_CALL( solClearArrays(sol) );
   SCIPrationalSetReal(sol->valsexact->obj, origprob->objoffset);
   sol->solorigin = SCIP_SOLORIGIN_ORIGINAL;

   /* reinsert the values of the original variables */
   for( v = 0; v < nvars; ++v )
   {
      /* we might require unchangedObjexact for this assert if exact probing mode is implemented */
      assert(SCIPvarGetUnchangedObj(vars[v]) == SCIPvarGetObj(vars[v])); /*lint !e777*/

      if( !SCIPrationalIsZero(solvals[v]) )
      {
         SCIP_CALL( solSetArrayValExact(sol, set, vars[v], solvals[v]) );
         SCIPrationalAddProd(sol->valsexact->obj, SCIPvarGetObjExact(vars[v]), solvals[v]);
      }
   }

   /* free temporary memory */
   SCIPrationalFreeBuffer(set->buffer, &constant);
   SCIPrationalFreeBufferArray(set->buffer, &solvals, nvars);
   SCIPrationalFreeBufferArray(set->buffer, &activevals, ntransvars + 1);
   SCIPsetFreeBufferArray(set, &activevars);
   SCIPrationalFreeBufferArray(set->buffer, &transsolvals, ntransvars + 1);

   return SCIP_OKAY;
}

/** recomputes the objective value of an original solution, e.g., when transferring solutions
 *  from the solution pool (objective coefficients might have changed in the meantime)
 */
void SCIPsolRecomputeObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob            /**< original problem */
   )
{
   SCIP_VAR** vars;
   SCIP_Real solval;
   int nvars;
   int v;

   assert(sol != NULL);
   assert(SCIPsolIsOriginal(sol));
   assert(origprob != NULL);

   vars = origprob->vars;
   nvars = origprob->nvars;

   /* recompute the objective value */
   sol->obj = SCIPprobGetObjoffset(origprob);
   for( v = 0; v < nvars; ++v )
   {
      solval = SCIPsolGetVal(sol, set, stat, vars[v]);
      if( solval != 0.0 && solval != SCIP_UNKNOWN ) /*lint !e777*/
      {
         sol->obj += SCIPvarGetUnchangedObj(vars[v]) * solval;
      }
   }

   if( SCIPsetIsInfinity(set, -sol->obj) )
      sol->obj = -SCIPsetInfinity(set);
}

/** recomputes the objective value of an exact solution, e.g., when initialized from a real solution */
void SCIPsolRecomputeInternObjExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob                /**< scip problem */
   )
{
   SCIP_VAR** vars;
   SCIP_RATIONAL* solval;
   int nvars;
   int v;

   assert(sol != NULL);
   assert(SCIPsolIsExact(sol));
   assert(prob != NULL);

   vars = prob->vars;
   nvars = prob->nvars;
   (void) SCIPrationalCreateBuffer(set->buffer, &solval);

   SCIPrationalSetFraction(sol->valsexact->obj, 0LL, 1LL);

   /* recompute the objective value */
   for( v = 0; v < nvars; ++v )
   {
      SCIPsolGetValExact(solval, sol, set, stat, vars[v]);
      if( !SCIPrationalIsZero(solval) ) /*lint !e777*/
      {
         SCIPrationalAddProd(sol->valsexact->obj, SCIPvarGetObjExact(vars[v]), solval);
      }
   }

   SCIPrationalFreeBuffer(set->buffer, &solval);
}

/** returns whether the given solutions (exact or floating point) are exactly equal */
static
SCIP_Bool solsAreEqualExact(
   SCIP_SOL*             sol1,               /**< first primal CIP solution */
   SCIP_SOL*             sol2,               /**< second primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob           /**< transformed problem after presolve, or NULL if both solution are
                                              *   defined in the original problem space */
   )
{
   SCIP_PROB* prob;
   SCIP_RATIONAL* tmp1;
   SCIP_RATIONAL* tmp2;
   int v;
   SCIP_Bool result = TRUE;

   assert(sol1 != NULL);
   assert(sol2 != NULL);
   assert(((SCIPsolGetOrigin(sol1) == SCIP_SOLORIGIN_ORIGINAL) && (SCIPsolGetOrigin(sol2) == SCIP_SOLORIGIN_ORIGINAL)) || transprob != NULL);

   (void) SCIPrationalCreateBuffer(set->buffer, &tmp1);
   (void) SCIPrationalCreateBuffer(set->buffer, &tmp2);

   /* if both solutions are original or both are transformed, take the objective values stored in the solutions */
   if( (SCIPsolGetOrigin(sol1) == SCIP_SOLORIGIN_ORIGINAL) == (SCIPsolGetOrigin(sol2) == SCIP_SOLORIGIN_ORIGINAL) )
   {
      SCIPsolIsExact(sol1) ? SCIPrationalSetRational(tmp1, sol1->valsexact->obj) : SCIPrationalSetReal(tmp1, sol1->obj);
      SCIPsolIsExact(sol2) ? SCIPrationalSetRational(tmp2, sol2->valsexact->obj) : SCIPrationalSetReal(tmp2, sol2->obj);
   }
   /* one solution is original and the other not, so we have to get for both the objective in the transformed problem */
   else
   {
      if( SCIPsolIsExact(sol1) )
         SCIPsolGetObjExact(sol1, set, transprob, origprob, tmp1);
      else
         SCIPrationalSetReal(tmp1, SCIPsolGetObj(sol1, set, transprob, origprob));
      if( SCIPsolIsExact(sol2) )
         SCIPsolGetObjExact(sol2, set, transprob, origprob, tmp2);
      else
         SCIPrationalSetReal(tmp2, SCIPsolGetObj(sol2, set, transprob, origprob));
   }

   /* solutions with different objective values cannot be the same */
   if( !SCIPrationalIsEQ(tmp1, tmp2) )
      result = FALSE;

   /* if one of the solutions is defined in the original space, the comparison has to be performed in the original
    * space
    */
   prob = transprob;
   if( sol1->solorigin == SCIP_SOLORIGIN_ORIGINAL || sol2->solorigin == SCIP_SOLORIGIN_ORIGINAL )
      prob = origprob;
   assert(prob != NULL);

   /* compare each variable value */
   for( v = 0; v < prob->nvars; ++v )
   {
      if( SCIPsolIsExact(sol1) )
         SCIPsolGetValExact(tmp1, sol1, set, stat, prob->vars[v]);
      else
         SCIPrationalSetReal(tmp1, SCIPsolGetVal(sol1, set, stat, prob->vars[v]));
      if( SCIPsolIsExact(sol2) )
         SCIPsolGetValExact(tmp2, sol2, set, stat, prob->vars[v]);
      else
         SCIPrationalSetReal(tmp2, SCIPsolGetVal(sol2, set, stat, prob->vars[v]));

      if( !SCIPrationalIsEQ(tmp1, tmp2) )
         result = FALSE;
   }

   SCIPrationalFreeBuffer(set->buffer, &tmp2);
   SCIPrationalFreeBuffer(set->buffer, &tmp1);

   return result;
}

/** returns whether the given solutions are equal */
SCIP_Bool SCIPsolsAreEqual(
   SCIP_SOL*             sol1,               /**< first primal CIP solution */
   SCIP_SOL*             sol2,               /**< second primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob           /**< transformed problem after presolve, or NULL if both solution are
                                              *   defined in the original problem space */
   )
{
   SCIP_PROB* prob;
   SCIP_Bool infobjs;
   SCIP_Real obj1;
   SCIP_Real obj2;
   int v;

   assert(sol1 != NULL);
   assert(sol2 != NULL);
   assert((SCIPsolIsOriginal(sol1) && SCIPsolIsOriginal(sol2)) || transprob != NULL);

   /* exact solutions should be checked exactly */
   if( set->exact_enable )
   {
      return solsAreEqualExact(sol1, sol2, set, stat, origprob, transprob);
   }

   /* if both solutions are original or both are transformed, take the objective values stored in the solutions */
   if( SCIPsolIsOriginal(sol1) == SCIPsolIsOriginal(sol2) )
   {
      obj1 = sol1->obj;
      obj2 = sol2->obj;
   }
   /* one solution is original and the other not, so we have to get for both the objective in the transformed problem */
   else
   {
      obj1 = SCIPsolGetObj(sol1, set, transprob, origprob);
      obj2 = SCIPsolGetObj(sol2, set, transprob, origprob);
   }

   /* solutions with different objective values cannot be the same; we consider two infinite objective values with the
    * same sign always to be different
    */
   infobjs = (SCIPsetIsInfinity(set, obj1) && SCIPsetIsInfinity(set, obj2))
      || (SCIPsetIsInfinity(set, -obj1) && SCIPsetIsInfinity(set, -obj2));
   if( !infobjs && !SCIPsetIsEQ(set, obj1, obj2) )
      return FALSE;

   /* if one of the solutions is defined in the original space, the comparison has to be performed in the original
    * space
    */
   prob = transprob;
   if( SCIPsolIsOriginal(sol1) || SCIPsolIsOriginal(sol2) )
      prob = origprob;
   assert(prob != NULL);

   /* compare each variable value */
   for( v = 0; v < prob->nvars; ++v )
   {
      SCIP_Real val1;
      SCIP_Real val2;

      val1 = SCIPsolGetVal(sol1, set, stat, prob->vars[v]);
      val2 = SCIPsolGetVal(sol2, set, stat, prob->vars[v]);
      if( !SCIPsetIsEQ(set, val1, val2) )
         return FALSE;
   }

   return TRUE;
}

/** outputs non-zero elements of solution to file stream */
SCIP_RETCODE SCIPsolPrint(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             mipstart,           /**< should only discrete variables be printed? */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_Real solval;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(SCIPsolIsOriginal(sol) || prob->transformed || transprob != NULL);
   assert(!mipstart || !SCIPsolIsPartial(sol));

   /* display variables of problem data */
   for( v = 0; v < prob->nfixedvars; ++v )
   {
      assert(prob->fixedvars[v] != NULL);

      /* skip non-discrete variables in a mip start */
      if( mipstart && !SCIPvarIsIntegral(prob->fixedvars[v]) )
         continue;

      solval = SCIPsolGetVal(sol, set, stat, prob->fixedvars[v]);
      if( printzeros || mipstart
         || (sol->solorigin != SCIP_SOLORIGIN_PARTIAL && !SCIPsetIsZero(set, solval))
         || (sol->solorigin == SCIP_SOLORIGIN_PARTIAL && solval != SCIP_UNKNOWN) ) /*lint !e777*/
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->fixedvars[v]));
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
         else if( SCIPsetIsInfinity(set, solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
         else if( SCIPsetIsInfinity(set, -solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g", solval);
         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->fixedvars[v]));
      }
   }

   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);

      /* skip non-discrete variables in a mip start */
      if( mipstart && !SCIPvarIsIntegral(prob->vars[v]) )
         continue;

      solval = SCIPsolGetVal(sol, set, stat, prob->vars[v]);
      if( printzeros || mipstart
         || (sol->solorigin != SCIP_SOLORIGIN_PARTIAL && !SCIPsetIsZero(set, solval))
         || (sol->solorigin == SCIP_SOLORIGIN_PARTIAL && solval != SCIP_UNKNOWN) ) /*lint !e777*/
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->vars[v]));
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
         else if( SCIPsetIsInfinity(set, solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
         else if( SCIPsetIsInfinity(set, -solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g", solval);
         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->vars[v]));
      }
   }

   /* display additional priced variables (if given problem data is original problem); consider these variables only
    * if there is at least one active pricer, otherwise we might print variables that have been added by, e.g., the
    * dual sparsify presolver (see #2946)
    */
   if( !prob->transformed && !SCIPsolIsOriginal(sol) && set->nactivepricers > 0 )
   {
      assert(transprob != NULL);
      for( v = 0; v < transprob->nfixedvars; ++v )
      {
         assert(transprob->fixedvars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->fixedvars[v]) )
            continue;

         /* skip non-discrete variables in a mip start */
         if( mipstart && !SCIPvarIsIntegral(transprob->fixedvars[v]) )
            continue;

         solval = SCIPsolGetVal(sol, set, stat, transprob->fixedvars[v]);
         if( printzeros || mipstart || !SCIPsetIsZero(set, solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->fixedvars[v]));
            if( solval == SCIP_UNKNOWN ) /*lint !e777*/
               SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
            else if( SCIPsetIsInfinity(set, solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g", solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->fixedvars[v]));
         }
      }
      for( v = 0; v < transprob->nvars; ++v )
      {
         assert(transprob->vars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->vars[v]) )
            continue;

         /* skip non-discrete variables in a mip start */
         if( mipstart && !SCIPvarIsIntegral(transprob->vars[v]) )
            continue;

         solval = SCIPsolGetVal(sol, set, stat, transprob->vars[v]);
         if( printzeros || !SCIPsetIsZero(set, solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->vars[v]));
            if( solval == SCIP_UNKNOWN ) /*lint !e777*/
               SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
            else if( SCIPsetIsInfinity(set, solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g", solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->vars[v]));
         }
      }
   }

   return SCIP_OKAY;
}

/** outputs non-zero elements of exact solution to file stream */
SCIP_RETCODE SCIPsolPrintExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             mipstart,           /**< should only discrete variables be printed? */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_RATIONAL* solval;
   char* solvalstr;
   int solvallen;
   int solvalsize = SCIP_MAXSTRLEN;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL || prob->transformed || transprob != NULL);
   assert(SCIPsolIsExact(sol));

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &solval) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &solvalstr, solvalsize) );

   /* display variables of problem data */
   for( v = 0; v < prob->nfixedvars; ++v )
   {
      assert(prob->fixedvars[v] != NULL);

      /* skip non-discrete variables in a mip start */
      if( mipstart && !SCIPvarIsIntegral(prob->fixedvars[v]) )
         continue;

      SCIPsolGetValExact(solval, sol, set, stat, prob->fixedvars[v]);
      if( printzeros || mipstart
         || (sol->solorigin != SCIP_SOLORIGIN_PARTIAL && !SCIPrationalIsZero(solval))
         || (sol->solorigin == SCIP_SOLORIGIN_PARTIAL) ) /*lint !e777*/
      {
         solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
         if( solvallen >= solvalsize )
         {
            solvalsize = SCIPrationalStrLen(solval) + 1;
            SCIP_CALL( SCIPsetReallocBufferArray(set, &solvalstr, solvalsize) );
            solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
            assert(solvallen < solvalsize);
         }

         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->fixedvars[v]));
         SCIPmessageFPrintInfo(messagehdlr, file, " %20s", solvalstr);

         solvallen = SCIPrationalToString(SCIPvarGetObjExact(prob->fixedvars[v]), solvalstr, solvalsize);
         if( solvallen >= solvalsize )
         {
            solvalsize = SCIPrationalStrLen(solval) + 1;
            SCIP_CALL( SCIPsetReallocBufferArray(set, &solvalstr, solvalsize) );
            solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
            assert(solvallen < solvalsize);
         }

         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%s)\n", solvalstr);
      }
   }

   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);

      /* skip non-discrete variables in a mip start */
      if( mipstart && !SCIPvarIsIntegral(prob->vars[v]) )
         continue;

      SCIPsolGetValExact(solval, sol, set, stat, prob->vars[v]);
      if( printzeros || mipstart
         || (sol->solorigin != SCIP_SOLORIGIN_PARTIAL && !SCIPrationalIsZero(solval)) ) /*lint !e777*/
      {
         solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
         if( solvallen >= solvalsize )
         {
            solvalsize = SCIPrationalStrLen(solval) + 1;
            SCIP_CALL( SCIPsetReallocBufferArray(set, &solvalstr, solvalsize) );
            solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
            assert(solvallen < solvalsize);
         }

         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->vars[v]));
         SCIPmessageFPrintInfo(messagehdlr, file, " %20s", solvalstr);

         solvallen = SCIPrationalToString(SCIPvarGetObjExact(prob->vars[v]), solvalstr, solvalsize);
         if( solvallen >= solvalsize )
         {
            solvalsize = SCIPrationalStrLen(solval) + 1;
            SCIP_CALL( SCIPsetReallocBufferArray(set, &solvalstr, solvalsize) );
            solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
            assert(solvallen < solvalsize);
         }

         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%s)\n", solvalstr);
      }
   }

   /* display additional priced variables (if given problem data is original problem) */
   if( !prob->transformed && sol->solorigin != SCIP_SOLORIGIN_ORIGINAL )
   {
      assert(transprob != NULL);
      for( v = 0; v < transprob->nfixedvars; ++v )
      {
         assert(transprob->fixedvars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->fixedvars[v]) )
            continue;

         /* skip non-discrete variables in a mip start */
         if( mipstart && !SCIPvarIsIntegral(transprob->fixedvars[v]) )
            continue;

         SCIPsolGetValExact(solval, sol, set, stat, transprob->fixedvars[v]);
         if( printzeros || mipstart || !SCIPrationalIsZero(solval) )
         {
            solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
            if( solvallen >= solvalsize )
            {
               solvalsize = SCIPrationalStrLen(solval) + 1;
               SCIP_CALL( SCIPsetReallocBufferArray(set, &solvalstr, solvalsize) );
               solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
               assert(solvallen < solvalsize);
            }

            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->fixedvars[v]));
            SCIPmessageFPrintInfo(messagehdlr, file, " %20s", solvalstr);

            solvallen = SCIPrationalToString(SCIPvarGetObjExact(transprob->fixedvars[v]), solvalstr, solvalsize);
            if( solvallen >= solvalsize )
            {
               solvalsize = SCIPrationalStrLen(solval) + 1;
               SCIP_CALL( SCIPsetReallocBufferArray(set, &solvalstr, solvalsize) );
               solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
               assert(solvallen < solvalsize);
            }

            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%s)\n", solvalstr);
         }
      }
      for( v = 0; v < transprob->nvars; ++v )
      {
         assert(transprob->vars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->vars[v]) )
            continue;

         /* skip non-discrete variables in a mip start */
         if( mipstart && !SCIPvarIsIntegral(transprob->vars[v]) )
            continue;

         SCIPsolGetValExact(solval, sol, set, stat, transprob->vars[v]);
         if( printzeros || !SCIPrationalIsZero(solval) )
         {
            solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
            if( solvallen >= solvalsize )
            {
               solvalsize = SCIPrationalStrLen(solval) + 1;
               SCIP_CALL( SCIPsetReallocBufferArray(set, &solvalstr, solvalsize) );
               solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
               assert(solvallen < solvalsize);
            }

            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->vars[v]));
            SCIPmessageFPrintInfo(messagehdlr, file, " %20s", solvalstr);

            solvallen = SCIPrationalToString(SCIPvarGetObjExact(transprob->vars[v]), solvalstr, solvalsize);
            if( solvallen >= solvalsize )
            {
               solvalsize = SCIPrationalStrLen(solval) + 1;
               SCIP_CALL( SCIPsetReallocBufferArray(set, &solvalstr, solvalsize) );
               solvallen = SCIPrationalToString(solval, solvalstr, solvalsize);
               assert(solvallen < solvalsize);
            }

            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%s)\n", solvalstr);
         }
      }
   }

   SCIPsetFreeBufferArray(set, &solvalstr);
   SCIPrationalFreeBuffer(set->buffer, &solval);

   return SCIP_OKAY;
}

/** outputs non-zero elements of solution representing a ray to file stream */
SCIP_RETCODE SCIPsolPrintRay(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_Real solval;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(SCIPsolIsOriginal(sol) || prob->transformed || transprob != NULL);

   /* display variables of problem data */
   for( v = 0; v < prob->nfixedvars; ++v )
   {
      assert(prob->fixedvars[v] != NULL);
      solval = SCIPsolGetRayVal(sol, set, stat, prob->fixedvars[v]);
      if( printzeros || !SCIPsetIsZero(set, solval) )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->fixedvars[v]));
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
         else if( SCIPsetIsInfinity(set, solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
         else if( SCIPsetIsInfinity(set, -solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g", solval);
         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->fixedvars[v]));
      }
   }
   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);
      solval = SCIPsolGetRayVal(sol, set, stat, prob->vars[v]);
      if( printzeros || !SCIPsetIsZero(set, solval) )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->vars[v]));
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
         else if( SCIPsetIsInfinity(set, solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
         else if( SCIPsetIsInfinity(set, -solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g", solval);
         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->vars[v]));
      }
   }

   /* display additional priced variables (if given problem data is original problem) */
   if( !prob->transformed && !SCIPsolIsOriginal(sol) )
   {
      assert(transprob != NULL);
      for( v = 0; v < transprob->nfixedvars; ++v )
      {
         assert(transprob->fixedvars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->fixedvars[v]) )
            continue;

         solval = SCIPsolGetRayVal(sol, set, stat, transprob->fixedvars[v]);
         if( printzeros || !SCIPsetIsZero(set, solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->fixedvars[v]));
            if( solval == SCIP_UNKNOWN ) /*lint !e777*/
               SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
            else if( SCIPsetIsInfinity(set, solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g", solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->fixedvars[v]));
         }
      }
      for( v = 0; v < transprob->nvars; ++v )
      {
         assert(transprob->vars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->vars[v]) )
            continue;

         solval = SCIPsolGetRayVal(sol, set, stat, transprob->vars[v]);
         if( printzeros || !SCIPsetIsZero(set, solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->vars[v]));
            if( solval == SCIP_UNKNOWN ) /*lint !e777*/
               SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
            else if( SCIPsetIsInfinity(set, solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g", solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->vars[v]));
         }
      }
   }

   return SCIP_OKAY;
}

/** set new origin type for a solution */
void SCIPsolSetOrigin(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SOLORIGIN        origin              /**< new origin type of the solution */
   )
{
   assert( sol != NULL );
   sol->solorigin = origin;
}

/*
 * methods for accumulated numerical violations of a solution
 */

/** reset violations of a solution */
void SCIPsolResetViolations(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   sol->viol.absviolbounds = 0.0;
   sol->viol.relviolbounds = 0.0;
   sol->viol.absviolintegrality = 0.0;
   sol->viol.absviollprows = 0.0;
   sol->viol.relviollprows = 0.0;
   sol->viol.absviolcons = 0.0;
   sol->viol.relviolcons = 0.0;
}

/** update integrality violation of a solution */
void SCIPsolUpdateIntegralityViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviolintegrality  /**< absolute violation of integrality */
   )
{
   assert(sol != NULL);

   sol->viol.absviolintegrality = MAX(sol->viol.absviolintegrality, absviolintegrality);
}

/** update bound violation of a solution */
void SCIPsolUpdateBoundViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviolbounds,      /**< absolute violation of bounds */
   SCIP_Real             relviolbounds       /**< relative violation of bounds */
   )
{
   assert(sol != NULL);

   sol->viol.absviolbounds = MAX(sol->viol.absviolbounds, absviolbounds);
   sol->viol.relviolbounds = MAX(sol->viol.relviolbounds, relviolbounds);
}

/** update LP row violation of a solution */
void SCIPsolUpdateLPRowViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviollprows,      /**< absolute violation of LP rows */
   SCIP_Real             relviollprows       /**< relative violation of LP rows */
   )
{
   assert(sol != NULL);

   sol->viol.absviollprows = MAX(sol->viol.absviollprows, absviollprows);
   sol->viol.relviollprows = MAX(sol->viol.relviollprows, relviollprows);
}

/** update constraint violation of a solution */
void SCIPsolUpdateConsViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviolcons,        /**< absolute violation of constraint */
   SCIP_Real             relviolcons         /**< relative violation of constraint */
   )
{
   assert(sol != NULL);

   sol->viol.absviolcons = MAX(sol->viol.absviolcons, absviolcons);
   sol->viol.relviolcons = MAX(sol->viol.relviolcons, relviolcons);
}

/** update violation of a constraint that is represented in the LP */
void SCIPsolUpdateLPConsViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation of constraint */
   SCIP_Real             relviol             /**< relative violation of constraint */
   )
{
   assert(sol != NULL);

   SCIPsolUpdateConsViolation(sol, absviol, relviol);
   SCIPsolUpdateLPRowViolation(sol, absviol, relviol);
}

/** get maximum absolute bound violation of solution */
SCIP_Real SCIPsolGetAbsBoundViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.absviolbounds;
}

/** get maximum relative bound violation of solution */
SCIP_Real SCIPsolGetRelBoundViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.relviolbounds;
}

/** get maximum absolute integrality violation of solution */
SCIP_Real SCIPsolGetAbsIntegralityViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.absviolintegrality;
}

/** get maximum absolute LP row violation of solution */
SCIP_Real SCIPsolGetAbsLPRowViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.absviollprows;
}

/** get maximum relative LP row violation of solution */
SCIP_Real SCIPsolGetRelLPRowViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.relviollprows;
}

/** get maximum absolute constraint violation of solution */
SCIP_Real SCIPsolGetAbsConsViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.absviolcons;
}

/** get maximum relative constraint violation of solution */
SCIP_Real SCIPsolGetRelConsViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->viol.relviolcons;
}

/** overwrite FP solution with exact values */
SCIP_RETCODE SCIPsolOverwriteFPSolWithExact(
   SCIP_SOL*             sol,                /**< exact primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< problem data */
   SCIP_PROB*            transprob,          /**< problem data */
   SCIP_TREE*            tree                /**< branch and bound tree, or NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_RATIONAL* solval;
   int nvars;
   int i;

   assert(sol != NULL);
   assert(SCIPsolIsExact(sol));

   vars = SCIPsolIsOriginal(sol) ? SCIPprobGetVars(origprob) : SCIPprobGetVars(transprob);
   nvars = SCIPsolIsOriginal(sol) ? SCIPprobGetNVars(origprob) : SCIPprobGetNVars(transprob);

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &solval) );

   /* overwrite all the variables */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_ROUNDMODE_RAT roundmode;
      SCIPsolGetValExact(solval, sol, set, stat, vars[i]);
      roundmode = vars[i]->obj > 0 ? SCIP_R_ROUND_UPWARDS : SCIP_R_ROUND_DOWNWARDS;

      SCIPrationalDebugMessage("overwriting value %g of var %s with value %g (%q) \n", SCIPsolGetVal(sol, set, stat, vars[i]),
           vars[i]->name, SCIPrationalRoundReal(solval, roundmode), solval);

      SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[i],
         SCIPrationalRoundReal(solval, roundmode)) );
   }

   if( SCIPsolIsOriginal(sol) )
   {
      SCIPrationalSetRational(solval, SCIPsolGetOrigObjExact(sol));
   }
   else
   {
      SCIPsolGetObjExact(sol, set, transprob, origprob, solval);
   }

   /* hard-set the obj value of the solution  */
   sol->obj = SCIPrationalRoundReal(solval, SCIP_R_ROUND_UPWARDS);

   SCIPrationalFreeBuffer(set->buffer, &solval);

   return SCIP_OKAY;
}

/** comparison method for sorting solution by decreasing objective value (best solution will be sorted to the end) */
SCIP_DECL_SORTPTRCOMP(SCIPsolComp)
{
   if( ((SCIP_SOL*)elem1)->obj - ((SCIP_SOL*)elem2)->obj < 0 )
      return 1;
   else if( ((SCIP_SOL*)elem1)->obj - ((SCIP_SOL*)elem2)->obj > 0 )
      return -1;
   else
      return 0;
}

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPsolGetOrigin
#undef SCIPsolIsOriginal
#undef SCIPsolGetOrigObj
#undef SCIPsolGetTime
#undef SCIPsolGetNodenum
#undef SCIPsolGetRunnum
#undef SCIPsolGetDepth
#undef SCIPsolGetHeur
#undef SCIPsolGetRelax
#undef SCIPsolOrigAddObjval
#undef SCIPsolGetPrimalIndex
#undef SCIPsolSetPrimalIndex
#undef SCIPsolGetIndex
#undef SCIPsolGetType
#undef SCIPsolSetLPRelaxation
#undef SCIPsolSetStrongbranching
#undef SCIPsolSetPseudo

/** gets origin of solution */
SCIP_SOLORIGIN SCIPsolGetOrigin(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->solorigin;
}

/** returns whether the given solution is defined on original variables */
SCIP_Bool SCIPsolIsOriginal(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return (sol->solorigin == SCIP_SOLORIGIN_ORIGINAL || sol->solorigin == SCIP_SOLORIGIN_PARTIAL);
}

/** returns whether a solution is an exact rational solution */
SCIP_Bool SCIPsolIsExact(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->valsexact != NULL;
}

/** returns whether the given solution is defined on original variables and containes unknown solution values */
SCIP_Bool SCIPsolIsPartial(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return (sol->solorigin == SCIP_SOLORIGIN_PARTIAL);
}

/** gets objective value of primal CIP solution which lives in the original problem space */
SCIP_Real SCIPsolGetOrigObj(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);
   assert(SCIPsolIsOriginal(sol));

   return sol->obj;
}

/** gets objective value of primal CIP solution which lives in the original problem space */
SCIP_RATIONAL* SCIPsolGetOrigObjExact(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);
   assert(SCIPsolIsOriginal(sol));
   assert(SCIPsolIsExact(sol));

   return sol->valsexact->obj;
}

/** adds value to the objective value of a given original primal CIP solution */
void SCIPsolOrigAddObjval(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             addval              /**< offset value to add */
   )
{
   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL);

   sol->obj += addval;
}

/** adds value to the objective value of a given original primal CIP solution */
void SCIPsolOrigAddObjvalExact(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RATIONAL*        addval              /**< offset value to add */
   )
{
   assert(sol != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL);
   assert(SCIPsolIsExact(sol));

   SCIPrationalAdd(sol->valsexact->obj, sol->valsexact->obj, addval);
   sol->obj = SCIPrationalGetReal(sol->valsexact->obj);
}

/** gets clock time, when this solution was found */
SCIP_Real SCIPsolGetTime(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->time;
}

/** gets branch and bound run number, where this solution was found */
int SCIPsolGetRunnum(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->runnum;
}

/** gets node number, where this solution was found */
SCIP_Longint SCIPsolGetNodenum(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->nodenum;
}

/** gets node's depth, where this solution was found */
int SCIPsolGetDepth(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->depth;
}

/** gets heuristic, that found this solution or NULL if solution has type different than SCIP_SOLTYPE_HEUR */
SCIP_HEUR* SCIPsolGetHeur(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->type == SCIP_SOLTYPE_HEUR ? sol->creator.heur : NULL;
}

/** gets current position of solution in array of existing solutions of primal data */
int SCIPsolGetPrimalIndex(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->primalindex;
}

/** sets current position of solution in array of existing solutions of primal data */
void SCIPsolSetPrimalIndex(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   primalindex         /**< new primal index of solution */
   )
{
   assert(sol != NULL);

   sol->primalindex = primalindex;
}

/** returns unique index of given solution */
int SCIPsolGetIndex(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->index;
}

/** informs the solution that it now belongs to the given primal heuristic. For convenience and backwards compatibility,
 *  the method accepts NULL as input for \p heur, in which case the solution type is set to SCIP_SOLTYPE_LPRELAX.
 *
 *  @note Relaxation handlers should use SCIPsolSetRelax() instead.
 */
void SCIPsolSetHeur(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_HEUR*            heur                /**< primal heuristic that found the solution, or NULL for LP solutions */
   )
{
   assert(sol != NULL);

   if( heur == NULL )
      SCIPsolSetLPRelaxation(sol);
   else
   {
      sol->type = SCIP_SOLTYPE_HEUR;
      sol->creator.heur = heur;
   }
}

/** gets information if solution was found by the LP, a primal heuristic, or a custom relaxator */
SCIP_SOLTYPE SCIPsolGetType(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->type;
}

/** gets relaxation handler that found this solution, or NULL if solution has different type than SCIP_SOLTYPE_RELAX */
SCIP_RELAX* SCIPsolGetRelax(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->type == SCIP_SOLTYPE_RELAX ? sol->creator.relax : NULL;
}

/** informs the solution that it now belongs to the given relaxation handler */
void SCIPsolSetRelax(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RELAX*           relax               /**< relaxator that found the solution */
   )
{
   assert(sol != NULL);
   assert(relax != NULL);

   sol->type = SCIP_SOLTYPE_RELAX;
   sol->creator.relax = relax;
}

/** informs the solution that it is an LP relaxation solution */
void SCIPsolSetLPRelaxation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   sol->type = SCIP_SOLTYPE_LPRELAX;
}

/** informs the solution that it is a solution found during strong branching */
void SCIPsolSetStrongbranching(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   sol->type = SCIP_SOLTYPE_STRONGBRANCH;
}

/** informs the solution that it originates from a pseudo solution */
void SCIPsolSetPseudo(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   sol->type = SCIP_SOLTYPE_PSEUDO;
}

