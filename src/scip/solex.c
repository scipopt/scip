/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   solex.c
 * @brief  methods for storing primal CIP solutions
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/clock.h"
#include "scip/cons.h"
#include "scip/lp.h"
#include "scip/lpex.h"
#include "scip/misc.h"
#include "scip/nlp.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/pub_varex.h"
#include "scip/relax.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solex.h"
#include "scip/stat.h"
#include "scip.h"
#include "scip/struct_lp.h"
#include "scip/struct_lpex.h"
#include "scip/struct_prob.h"
#include "scip/struct_set.h"
#include "scip/struct_sol.h"
#include "scip/struct_stat.h"
#include "scip/struct_var.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/varex.h"
#include "scip/rational.h"


/** clears solution arrays of primal CIP solution */
static
SCIP_RETCODE solexClearArrays(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   SCIP_CALL( SCIPboolarrayClear(sol->valsex->valid) );
   sol->hasinfval = FALSE;

   return SCIP_OKAY;
}

/** sets value of variable in the solution's array */
static
SCIP_RETCODE solexSetArrayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Rational*        val                 /**< value to set variable to */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* from now on, variable must not be deleted */
   SCIPvarMarkNotDeletable(var);

   /* mark the variable valid */
   SCIP_CALL( SCIPboolarraySetVal(sol->valsex->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, TRUE) );

   /* set the value in the solution array */
   SCIP_CALL( SCIPrationalarraySetVal(sol->valsex->vals, idx, val) );

   /* store whether the solution has infinite values assigned to variables */
   if( RatIsAbsInfinity(val) ) /*lint !e777*/
      sol->hasinfval = TRUE;

   return SCIP_OKAY;
}

/** increases value of variable in the solution's array */
static
SCIP_RETCODE solexIncArrayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Rational*        incval              /**< increase of variable's solution value */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* from now on, variable must not be deleted */
   SCIPvarMarkNotDeletable(var);

   /* if the variable was not valid, mark it to be valid and set the value to the incval (it is 0.0 if not valid) */
   if( !SCIPboolarrayGetVal(sol->valsex->valid, idx) )
   {
      /* mark the variable valid */
      SCIP_CALL( SCIPboolarraySetVal(sol->valsex->valid, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, TRUE) );

      /* set the value in the solution array */
      SCIP_CALL( SCIPrationalarraySetVal(sol->valsex->vals, idx, incval) );
   }
   else
   {
      /* increase the value in the solution array */
      SCIP_CALL( SCIPrationalarrayIncVal(sol->valsex->vals, idx, incval) );
   }

   /* store whether the solution has infinite values assigned to variables */
   if( RatIsAbsInfinity(incval) )
      sol->hasinfval = TRUE;

   return SCIP_OKAY;
}

/** returns the value of the variable in the given solution */
static
void solexGetArrayVal(
   SCIP_Rational*        res,
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   int idx;

   assert(sol != NULL);

   idx = SCIPvarGetIndex(var);

   /* check, if the variable's value is valid */
   if( SCIPboolarrayGetVal(sol->valsex->valid, idx) )
   {
      SCIPrationalarrayGetVal(sol->valsex->vals, idx, res);
   }
   else
   {
      /* return the variable's value corresponding to the origin */
      switch( sol->solorigin )
      {
      case SCIP_SOLORIGIN_ORIGINAL:
      case SCIP_SOLORIGIN_ZERO:
         RatSetReal(res, 0.0);
         break;

      case SCIP_SOLORIGIN_LPSOL:
         RatSet(res, SCIPvarGetLPexSolex(var));
         break;

      case SCIP_SOLORIGIN_PSEUDOSOL:
         RatSet(res, SCIPvarGetPseudoSolex(var));
         break;

      case SCIP_SOLORIGIN_PARTIAL:
      case SCIP_SOLORIGIN_UNKNOWN:
      default:
         SCIPerrorMessage("unknown solution origin <%d>\n", sol->solorigin);
         SCIPABORT();
         return; /*lint !e527*/
      }
   }
}

/** stores solution value of variable in solution's own array */
static
SCIP_RETCODE solexUnlinkVar(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Rational* solval;

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
      SCIP_CALL( solexSetArrayVal(sol, set, var, SCIPvarGetLPexSolex(var) ) );
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_PSEUDOSOL:
      SCIP_CALL( solexSetArrayVal(sol, set, var, SCIPvarGetPseudoSolex(var)) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("unknown solution origin <%d>\n", sol->solorigin);
      return SCIP_INVALIDDATA;
   }
}

/** creates primal CIP solution with exact rational values, initialized to zero */
SCIP_RETCODE SCIPsolexCreate(
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

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &(*sol)->valsex ) );
   SCIP_CALL( SCIPrationalarrayCreate(&(*sol)->valsex->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayCreate(&(*sol)->valsex->valid, blkmem) );
   SCIP_CALL( RatCreateBlock(blkmem, &(*sol)->valsex->obj) );

   return SCIP_OKAY;
}

/** creates a copy of a primal CIP solution */
SCIP_RETCODE SCIPvalsexCopy(
   SCIP_VALSEX**         valsex,             /**< pointer to store the copy of the primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VALSEX*          sourcevals          /**< primal CIP solution to copy */
   )
{
   assert(valsex != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, valsex) );
   SCIP_CALL( SCIPrationalarrayCopy(&(*valsex)->vals, blkmem, sourcevals->vals) );
   SCIP_CALL( SCIPboolarrayCopy(&(*valsex)->valid, blkmem, sourcevals->valid) );
   SCIP_CALL( RatCopy(blkmem, &(*valsex)->obj, sourcevals->obj) );

   return SCIP_OKAY;
}

/** creates primal CIP solution with exact rational values, initialized to the current LP solution */
SCIP_RETCODE SCIPsolexCreateLPexSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);
   assert(lp != NULL);
   assert(lp->solved);

   SCIP_CALL( SCIPsolexCreate(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolexLinkLPexSol(*sol, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to the current pseudo solution */
SCIP_RETCODE SCIPsolexCreatePseudoSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(sol != NULL);

   SCIP_CALL( SCIPsolexCreate(sol, blkmem, set, stat, primal, tree, heur) );
   SCIP_CALL( SCIPsolexLinkPseudoSol(*sol, set, prob, lp) );

   return SCIP_OKAY;
}

/** creates primal CIP solution with exact rational values, initialized to the current solution */
SCIP_RETCODE SCIPsolexCreateCurrentSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(tree != NULL);

   if( SCIPtreeHasCurrentNodeLP(tree) )
   {
      assert(lp->solved);
      SCIP_CALL( SCIPsolexCreateLPexSol(sol, blkmem, set, stat, prob, primal, tree, lp, heur) );
   }
   else
   {
      SCIP_CALL( SCIPsolexCreatePseudoSol(sol, blkmem, set, stat, prob, primal, tree, lp, heur) );
   }

   return SCIP_OKAY;
}

/** frees primal CIP solution */
SCIP_RETCODE SCIPvalsexFree(
   SCIP_VALSEX**         valsex,             /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_SOL* fpsol;

   assert(valsex != NULL);
   assert(*valsex != NULL);

   RatFreeBlock(blkmem, &(*valsex)->obj);
   SCIP_CALL( SCIPrationalarrayFree(&(*valsex)->vals, blkmem) );
   SCIP_CALL( SCIPboolarrayFree(&(*valsex)->valid) );
   BMSfreeBlockMemory(blkmem, valsex);

   return SCIP_OKAY;
}

/** copies current exact LP solution into CIP solution by linking */
SCIP_RETCODE SCIPsolexLinkLPexSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_LPEX*            lp                  /**< current LP data */
   )
{
   assert(sol != NULL);
   assert(lp != NULL);
   assert(lp->solved);

   /* clear the old solution arrays */
   SCIP_CALL( solexClearArrays(sol) );

   /* the objective value in the columns is correct, s.t. the LP's objective value is also correct */
   RatSet(sol->valsex->obj, lp->lpobjval);
   sol->obj = RatRoundReal(sol->valsex->obj, SCIP_ROUND_UPWARDS);
   sol->solorigin = SCIP_SOLORIGIN_LPSOL;

   return SCIP_OKAY;
}

/** copies current pseudo solution into CIP solution by linking */
SCIP_RETCODE SCIPsolexLinkPseudoSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_LPEX*            lp                  /**< current LP data */
   )
{
   assert(sol != NULL);

   /* clear the old solution arrays */
   SCIP_CALL( solexClearArrays(sol) );

   /* link solution to pseudo solution */
   SCIPlpexGetPseudoObjval(lp, set, prob, sol->valsex->obj);

   SCIPsetDebugMsg(set, " -> objective value: %g\n", sol->obj);

   return SCIP_OKAY;
}

/** clears primal CIP solution */
SCIP_RETCODE SCIPsolexClear(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   SCIP_CALL( solexClearArrays(sol) );
   RatSetReal(sol->valsex->obj, 0.0);

   return SCIP_OKAY;
}


/** stores solution values of variables in solution's own array */
SCIP_RETCODE SCIPsolexUnlink(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< transformed problem data */
   )
{
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(prob->nvars == 0 || prob->vars != NULL);

   if( sol->solorigin != SCIP_SOLORIGIN_ORIGINAL && sol->solorigin != SCIP_SOLORIGIN_ZERO
      && sol->solorigin != SCIP_SOLORIGIN_UNKNOWN )
   {
      SCIPsetDebugMsg(set, "completing solution %p\n", (void*)sol);

      for( v = 0; v < prob->nvars; ++v )
      {
         SCIP_CALL( solexUnlinkVar(sol, set, prob->vars[v]) );
      }

      sol->solorigin = SCIP_SOLORIGIN_ZERO;
   }

   return SCIP_OKAY;
}

/** sets value of variable in primal CIP solution */
SCIP_RETCODE SCIPsolexSetVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_Rational*        val                 /**< solution value of variable */
   )
{
   SCIP_Rational* oldval;
   SCIP_Rational* tmp;

   assert(sol != NULL);
   assert(stat != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL
      || sol->solorigin == SCIP_SOLORIGIN_ZERO
      || sol->solorigin == SCIP_SOLORIGIN_PARTIAL
      || sol->solorigin == SCIP_SOLORIGIN_UNKNOWN);
   assert(var != NULL);
   assert(!RatIsAbsInfinity(val));

   SCIPsetDebugMsg(set, "setting value of <%s> in exact solution %p to %g\n", SCIPvarGetName(var), (void*)sol, RatApproxReal(val));

   SCIP_CALL( RatCreateBuffer(set->buffer, &oldval) );

   /* we want to store only values for non fixed variables (LOOSE or COLUMN); others have to be transformed */
   switch( SCIPvarGetStatusExact(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( sol->solorigin == SCIP_SOLORIGIN_ORIGINAL )
      {
         solexGetArrayVal(oldval, sol, var);

         if( !RatIsEqual(val, oldval) )
         {
            SCIP_Rational* obj;

            SCIP_CALL( solexSetArrayVal(sol, set, var, val) );
            obj = SCIPvarGetObjExact(var);
            RatDiffProd(sol->valsex->obj, obj, oldval);

            RatAddProd(sol->valsex->obj, obj, val);
         }
      }
      else
         return SCIPsolexSetVal(sol, set, stat, tree, SCIPvarGetTransVar(var), val);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);

      solexGetArrayVal(oldval, sol, var);

      if( !RatIsEqual(val, oldval) )
      {
         SCIP_Rational* obj;
         SCIP_CALL( solexSetArrayVal(sol, set, var, val) );
         obj = SCIPvarGetObjExact(var);
         RatDiffProd(sol->valsex->obj, obj, oldval);
         RatAddProd(sol->valsex->obj, obj, val);
      }
      return SCIP_OKAY;

   case SCIP_VARSTATUS_FIXED:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);
      if( !RatIsEqual(val, SCIPvarGetLbGlobalExact(var)) )
      {
         SCIPerrorMessage("cannot set solution value for variable <%s> fixed to %.15g to different value %.15g\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobalExact(var), val);
         return SCIP_INVALIDDATA;
      }
      return SCIP_OKAY;

      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      /** @todo exip: presolving extension
       *  - implement this if exact version of SCIP supports aggregated variables
       */
      SCIPerrorMessage("cannot set solution value for aggregated variable\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("cannot set solution value for multiple aggregated variable\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_NEGATED:
      assert(!SCIPsetIsInfinity(set, SCIPvarGetNegationConstant(var)) && !SCIPsetIsInfinity(set, -SCIPvarGetNegationConstant(var)));

      return SCIPsolexSetVal(sol, set, stat, tree, SCIPvarGetNegationVar(var), val);

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
}

/** overwrite FP solution with exact values */
SCIP_RETCODE SCIPsolexOverwriteFPSol(
   SCIP_SOL*             sol,                /**< exact primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< problem data */
   SCIP_PROB*            transprob,          /**< problem data */
   SCIP_TREE*            tree                /**< branch and bound tree, or NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_Rational* solval;
   SCIP_Rational* obj;
   int nvars;
   int i;

   assert(sol != NULL);

   vars = SCIPprobGetVars(transprob);
   nvars = SCIPprobGetNVars(transprob);

   SCIP_CALL( RatCreateBuffer(set->buffer, &solval) );

   /* overwrite all the variables */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_ROUNDMODE roundmode;
      SCIPsolexGetVal(solval, sol, set, stat, vars[i]);
      roundmode = vars[i]->obj > 0 ? SCIP_ROUND_UPWARDS : SCIP_ROUND_DOWNWARDS;
      SCIP_CALL( SCIPsolSetVal(sol, set, stat, tree, vars[i],
         RatRoundReal(solval, roundmode)) );
   }

   obj = SCIPsolexGetObj(sol, set, transprob, origprob);
   /* hard-set the obj value of the solution  */
   sol->obj = RatRoundReal(obj, SCIP_ROUND_UPWARDS);

   RatFreeBuffer(set->buffer, &solval);

   return SCIP_OKAY;
}

/** returns value of variable in primal CIP solution */
void SCIPsolexGetVal(
   SCIP_Rational*        res,                /**< resulting rational */
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
      || sol->solorigin == SCIP_SOLORIGIN_LPSOL
      || sol->solorigin == SCIP_SOLORIGIN_UNKNOWN);
   assert(var != NULL);

   /* if the value of a transformed variable in an original solution is requested, we need to project the variable back
    * to the original space, the opposite case is handled below
    */
   if( sol->solorigin == SCIP_SOLORIGIN_ORIGINAL && SCIPvarIsTransformed(var) )
   {
      SCIP_RETCODE retcode;
      SCIP_VAR* origvar;
      SCIP_Rational* scalar;
      SCIP_Rational* constant;

      RatCreateBuffer(set->buffer, &scalar);
      RatCreateBuffer(set->buffer, &constant);
      /* we cannot get the value of a transformed variable for a solution that lives in the original problem space
       * -> get the corresponding original variable first
       */
      origvar = var;
      RatSetInt(scalar, 1, 1);
      RatSetReal(constant, 0.0);
      retcode = SCIPvarGetOrigvarSumExact(&origvar, scalar, constant);
      if ( retcode != SCIP_OKAY )
         return;
      if( origvar == NULL )
      {
         /* the variable has no original counterpart: in the original solution, it has a value of zero */
         RatSetReal(res, 0.0);
         return;
      }

      assert(!SCIPvarIsTransformed(origvar));

      SCIPsolexGetVal(res, sol, set, stat, origvar);
      RatMult(res, res, scalar);
      RatAdd(res, res, constant);

      RatFreeBuffer(set->buffer, &constant);
      RatFreeBuffer(set->buffer, &scalar);

      return;
   }

   /* only values for non fixed variables (LOOSE or COLUMN) are stored; others have to be transformed */
   switch( SCIPvarGetStatusExact(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( sol->solorigin == SCIP_SOLORIGIN_ORIGINAL )
         solexGetArrayVal(res, sol, var);
      else
         SCIPsolexGetVal(res, sol, set, stat, SCIPvarGetTransVar(var));
      break;

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);
      assert(sol->solorigin != SCIP_SOLORIGIN_LPSOL || SCIPboolarrayGetVal(sol->valid, SCIPvarGetIndex(var))
      || sol->lpcount == stat->lpcount );
      solexGetArrayVal(res, sol, var);
      break;

   case SCIP_VARSTATUS_FIXED:
      assert(sol->solorigin != SCIP_SOLORIGIN_ORIGINAL);
      assert(RatIsEqual(SCIPvarGetLbGlobalExact(var), SCIPvarGetUbGlobalExact(var))); /*lint !e777*/
      assert(RatIsEqual(SCIPvarGetLbLocalExact(var), SCIPvarGetUbLocalExact(var))); /*lint !e777*/
      assert(RatIsEqual(SCIPvarGetLbGlobalExact(var), SCIPvarGetLbLocalExact(var))); /*lint !e777*/
      SCIPvarGetLbGlobalExact(var);
      break;

    case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  =>  y = (x-c)/a */
      /** @todo exip: presolving extension
       *  - implement this if exact version of SCIP supports aggregated variables
       */
      SCIPerrorMessage("cannot get solution value of aggregated variable\n");
      SCIPABORT();
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      /** @todo exip: presolving extension
       *  - implement this if exact version of SCIP supports multiaggregated variables
       */
      SCIPerrorMessage("cannot get solution value of multiaggregated variable\n");
      SCIPABORT();
      break;

   case SCIP_VARSTATUS_NEGATED:
      SCIPsolexGetVal(res, sol, set, stat, SCIPvarGetNegationVar(var));
      RatDiffReal(res, res, SCIPvarGetNegationConstant(var));
      RatNegate(res, res);
      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      RatSetReal(res, 0.0); /*lint !e527*/
   }
}

/** gets objective value of primal CIP solution in transformed problem */
SCIP_Rational* SCIPsolexGetObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob            /**< original problem data */
   )
{
   assert(sol != NULL);

   /* for original solutions, sol->obj contains the external objective value */
   if( sol->solorigin == SCIP_SOLORIGIN_ORIGINAL )
   {
      /** @todo exip: heuristics extension
       *  - implement this if exact version of SCIP supports getting objective value of original solutions */
      SCIPerrorMessage("cannot get objectiv value of original solution\n");
      SCIPABORT();
      return NULL;
    //  return SCIPprobInternObjval(transprob, origprob, set, sol->obj);
   }
   else
      return sol->valsex->obj;
}

/** gets objective value of primal CIP solution which lives in the original problem space */
SCIP_Rational* SCIPsolexGetOrigObj(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);
   assert(SCIPsolIsOriginal(sol));

   return sol->valsex->obj;
}

/** retransforms solution to original problem space */
SCIP_RETCODE SCIPsolexRetransform(
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
   SCIP_Rational** transsolvals;
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
   SCIP_CALL( RatCreateBufferArray(set->buffer, &transsolvals, ntransvars) );
   assert(transsolvals != NULL); /* for flexelint */

   /* get the solution values of all active variables */
   for( v = 0; v < ntransvars; ++v )
   {
      SCIPsolexGetVal(transsolvals[v], sol, set, stat, transvars[v]);
   }

   /** @todo exip: presolving extension (once we have exact presolving, we need to do here what we do in the fp case */

   /* clear the solution and convert it into original space */
   SCIP_CALL( solexClearArrays(sol) );
   RatSetReal(sol->valsex->obj, origprob->objoffset);

   /* reinsert the values of the original variables */
   for( v = 0; v < nvars; ++v )
   {
      assert(SCIPvarGetUnchangedObj(vars[v]) == SCIPvarGetObj(vars[v])); /*lint !e777*/

      if( !RatIsZero(transsolvals[v]) )
      {
         SCIP_CALL( solexSetArrayVal(sol, set, vars[v], transsolvals[v]) );
         RatAddProd(sol->valsex->obj, SCIPvarGetObjExact(vars[v]), transsolvals[v]);
      }
   }

   /* free temporary memory */
   RatFreeBufferArray(set->buffer, &transsolvals, ntransvars);

   return SCIP_OKAY;
}

/** outputs non-zero elements of solution to file stream */
SCIP_RETCODE SCIPsolexPrint(
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
   SCIP_Rational* solval;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);
   assert(sol->solorigin == SCIP_SOLORIGIN_ORIGINAL || prob->transformed || transprob != NULL);

   SCIP_CALL( RatCreateBuffer(set->buffer, &solval) );

   /* display variables of problem data */
   for( v = 0; v < prob->nfixedvars; ++v )
   {
      assert(prob->fixedvars[v] != NULL);

      /* skip non-discrete variables in a mip start */
      if( mipstart && !SCIPvarIsIntegral(prob->fixedvars[v]) )
         continue;

      SCIPsolexGetVal(solval, sol, set, stat, prob->fixedvars[v]);
      if( printzeros || mipstart
         || (sol->solorigin != SCIP_SOLORIGIN_PARTIAL && !RatIsZero(solval))
         || (sol->solorigin == SCIP_SOLORIGIN_PARTIAL) ) /*lint !e777*/
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->fixedvars[v]));
         RatMessage(messagehdlr, file, solval);

         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->fixedvars[v]));
      }
   }

   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);

      /* skip non-discrete variables in a mip start */
      if( mipstart && !SCIPvarIsIntegral(prob->vars[v]) )
         continue;

      SCIPsolexGetVal(solval, sol, set, stat, prob->vars[v]);
      if( printzeros || mipstart
         || (sol->solorigin != SCIP_SOLORIGIN_PARTIAL && !RatIsZero(solval)) ) /*lint !e777*/
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->vars[v]));
         RatMessage(messagehdlr, file, solval);
         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->vars[v]));
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

         SCIPsolexGetVal(solval, sol, set, stat, transprob->fixedvars[v]);
         if( printzeros || mipstart || !RatIsZero(solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->fixedvars[v]));
            RatMessage(messagehdlr, file, solval);
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

         SCIPsolexGetVal(solval, sol, set, stat, transprob->vars[v]);
         if( printzeros || !RatIsZero(solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->vars[v]));
            RatMessage(messagehdlr, file, solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->vars[v]));
         }
      }
   }

   RatFreeBuffer(set->buffer, &solval);

   return SCIP_OKAY;
}

/** returns whether a solution is an exact rational solution */
SCIP_Bool SCIPsolIsExact(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->valsex != NULL;
}
