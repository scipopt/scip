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
#include "scip/struct_valsex.h"
#include "scip/struct_stat.h"
#include "scip/struct_var.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/varex.h"
#include "scip/rational.h"


/** clears solution arrays of primal CIP solution */
static
SCIP_RETCODE solexClearArrays(
   SCIP_SOL*             sol                /**< primal CIP solution */
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

#if 0
/** sets the solution time, nodenum, runnum, and depth stamp to the current values */
static
void solexStamp(
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
#endif

/** creates primal CIP solution, initialized to zero */
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

#if 0
/** creates primal CIP solution in original problem space, initialized to the offset in the original problem */
SCIP_RETCODE SCIPsolCreateOriginal(
   SCIP_SOLEX**          sol,                /**< pointer to primal CIP solution */
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
   (*sol)->heur = heur;
   (*sol)->solorigin = SCIP_SOLORIGIN_ORIGINAL;
   (*sol)->obj = origprob->objoffset;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   stat->solindex++;
   solStamp(*sol, stat, tree, TRUE);
   SCIPsolResetViolations(*sol);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}
#endif

/** creates a copy of a primal CIP solution */
SCIP_RETCODE SCIPvalsexCopy(
   SCIP_VALSEX**         valsex,             /**< pointer to store the copy of the primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VALSEX*          sourcevals           /**< primal CIP solution to copy */
   )
{
   assert(valsex != NULL);
   assert(*valsex != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, valsex) );
   SCIP_CALL( SCIPrationalarrayCopy(&(*valsex)->vals, blkmem, sourcevals->vals) );
   SCIP_CALL( SCIPboolarrayCopy(&(*valsex)->valid, blkmem, sourcevals->valid) );
   SCIP_CALL( RatCopy(blkmem, &(*valsex)->obj, sourcevals->obj) );

   return SCIP_OKAY;
}

#if 0
/** transformes given original solution to the transformed space; a corresponding transformed solution has to be given
 *  which is copied into the existing solution and freed afterwards
 */
SCIP_RETCODE SCIPsolTransform(
   SCIP_SOL*             sol,                /**< primal CIP solution to change, living in original space */
   SCIP_SOLEX**          transsol,           /**< pointer to corresponding transformed primal CIP solution */
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
   assert(sol->solorigin = SCIP_SOLORIGIN_ORIGINAL);
   assert(sol->primalindex > -1);

   tsol = *transsol;
   assert(tsol != NULL);
   assert(!SCIPsolIsOriginal(tsol));

   /* switch vals and valid arrays; the exisiting solution gets the arrays of the transformed solution;
    * the transformed one gets the original arrays, because they have to be freed anyway and freeing the transsol
    * automatically frees its arrays
    */
   tmpvals = sol->valsex->vals;
   tmpvalid = sol->valid;
   sol->valsex->vals = tsol->valsex->vals;
   sol->valid = tsol->valid;
   tsol->valsex->vals = tmpvals;
   tsol->valid = tmpvalid;

   /* copy solorigin and objective (should be the same, only to avoid numerical issues);
    * we keep the other statistics of the original solution, since that was the first time that this solution as found
    */
   sol->solorigin = tsol->solorigin;
   sol->obj = tsol->obj;

   SCIP_CALL( SCIPsolFree(transsol, blkmem, primal) );

   return SCIP_OKAY;
}


/** adjusts solution values of implicit integer variables in handed solution. Solution objective value is not
 *  deteriorated by this method.
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
   int nimplvars;
   int nbinvars;
   int nintvars;
   int v;

   assert(sol != NULL);
   assert(prob != NULL);

   /* get variable data */
   vars = SCIPprobGetVars(prob);
   nbinvars = SCIPprobGetNBinVars(prob);
   nintvars = SCIPprobGetNIntVars(prob);
   nimplvars = SCIPprobGetNImplVars(prob);

   if( nimplvars == 0 )
      return SCIP_OKAY;

   /* calculate the last array position of implicit integer variables */
   nimplvars = nbinvars + nintvars + nimplvars;

   /* loop over implicit integer variables and round them up or down */
   for( v = nbinvars + nintvars; v < nimplvars; ++v )
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

      assert( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT );
      SCIPsolexGetVal(solval, sol, set, stat, var);

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
#endif

/** creates primal CIP solution, initialized to the current LP solution */
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

#if 0

/** creates primal CIP solution, initialized to the current NLP solution */
SCIP_RETCODE SCIPsolCreateNLPSol(
   SCIP_SOLEX**          sol,                /**< pointer to primal CIP solution */
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
   SCIP_SOLEX**          sol,                /**< pointer to primal CIP solution */
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

   return SCIP_OKAY;
}
#endif

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

/** creates primal CIP solution, initialized to the current solution */
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

#if 0
/** creates partial primal CIP solution, initialized to unknown values */
SCIP_RETCODE SCIPsolCreatePartial(
   SCIP_SOLEX**          sol,                /**< pointer to primal CIP solution */
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
   (*sol)->heur = heur;
   (*sol)->solorigin = SCIP_SOLORIGIN_PARTIAL;
   (*sol)->obj = SCIP_UNKNOWN;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   stat->solindex++;
   solStamp(*sol, stat, NULL, TRUE);
   SCIPsolResetViolations(*sol);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}

/** creates primal CIP solution, initialized to unknown values */
SCIP_RETCODE SCIPsolCreateUnknown(
   SCIP_SOLEX**          sol,                /**< pointer to primal CIP solution */
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
   (*sol)->heur = heur;
   (*sol)->solorigin = SCIP_SOLORIGIN_UNKNOWN;
   (*sol)->obj = 0.0;
   (*sol)->primalindex = -1;
   (*sol)->index = stat->solindex;
   (*sol)->hasinfval = FALSE;
   stat->solindex++;
   solStamp(*sol, stat, tree, TRUE);
   SCIPsolResetViolations(*sol);

   SCIP_CALL( SCIPprimalSolCreated(primal, set, *sol) );

   return SCIP_OKAY;
}
#endif

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

#if 0
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
#endif

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
      /** @todo exiptodo: presolving extension
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
   SCIPsolSetObjVal(sol, RatRoundReal(obj, SCIP_ROUND_UPWARDS));

   RatFreeBuffer(set->buffer, &solval);

   return SCIP_OKAY;
}

#if 0
/** increases value of variable in primal CIP solution */
SCIP_RETCODE SCIPsolexIncVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR*             var,                /**< variable to increase solution value for */
   SCIP_Rational*        incval              /**< increment for solution value of variable */
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

   if( SCIPsetIsZero(set, incval) )
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
      if( sol->solorigin = SCIP_SOLORIGIN_ORIGINAL )
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
      assert(!sol->solorigin = SCIP_SOLORIGIN_ORIGINAL);
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
#endif

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
   /** @todo exip: maybe add check if lpsol */
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
   /** @todo exip: this check with lpcount has to be adapted */
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
      /** @todo exiptodo: presolving extension
       *  - implement this if exact version of SCIP supports aggregated variables
       */
      SCIPerrorMessage("cannot get solution value of aggregated variable\n");
      SCIPABORT();
      break;

   case SCIP_VARSTATUS_MULTAGGR:
      /** @todo exiptodo: presolving extension
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

#if 0
/** returns value of variable in primal ray represented by primal CIP solution */
SCIP_Real SCIPsolexGetRayVal(
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
      assert(!sol->solorigin = SCIP_SOLORIGIN_ORIGINAL);
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
#endif

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
      /** @todo exiptodo: heuristics extension
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

#if 0
/** updates primal solutions after a change in a variable's objective value */
void SCIPsolUpdateVarObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Rational*        oldobj,             /**< old objective value */
   SCIP_Rational*        newobj              /**< new objective value */
   )
{
   SCIP_Real solval;

   assert(sol != NULL);
   assert(!sol->solorigin = SCIP_SOLORIGIN_ORIGINAL);
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
         SCIP_CALL( SCIPrationalarraySetVal(sol->valsex->vals, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, vals[v]) );
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
#endif

/** checks primal CIP solution for feasibility
 *
 *  @note The difference between SCIPsolCheck() and SCIPcheckSolOrig() is that modifiable constraints are handled
 *        differently. There might be some variables which do not have an original counter part (e.g. in
 *        branch-and-price). Therefore, modifiable constraints can not be double-checked in the original space.
 */
SCIP_RETCODE SCIPsolexCheck(
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
   SCIP_Rational* solval;
   int h;

   assert(sol != NULL);
   assert(!SCIPsolIsOriginal(sol));
   assert(set != NULL);
   assert(prob != NULL);
   assert(feasible != NULL);

   SCIPsetDebugMsg(set, "checking solution with objective value %g (nodenum=%" SCIP_LONGINT_FORMAT ", origin=%u)\n",
      RatApproxReal(sol->valsex->obj), sol->nodenum, sol->solorigin);

   *feasible = TRUE;

   if( !printreason )
      completely = FALSE;

   SCIP_CALL( RatCreateBuffer(set->buffer, &solval) );

   /* check whether the solution respects the global bounds of the variables */
   if( checkbounds || sol->hasinfval )
   {
      int v;

      for( v = 0; v < prob->nvars && (*feasible || completely); ++v )
      {
         SCIP_VAR* var;

         var = prob->vars[v];
         SCIPsolexGetVal(solval, sol, set, stat, var);

         if( !RatIsAbsInfinity(solval) ) /*lint !e777*/
         {
            SCIP_Rational* lb;
            SCIP_Rational* ub;

            lb = SCIPvarGetLbGlobalExact(var);
            ub = SCIPvarGetUbGlobalExact(var);

            /* if we have to check bound and one of the current bounds is violated */
            if( checkbounds && ((!RatIsNegInfinity(lb) && RatIsLT(solval, lb))
                     || (!RatIsInfinity(ub) && RatIsGT(solval, ub))) )
            {
               *feasible = FALSE;

               if( printreason )
               {
                  SCIPmessagePrintInfo(messagehdlr, "solution value %g violates bounds of <%s>[%g,%g] by %g\n", solval, SCIPvarGetName(var),
                        SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), MAX(RatApproxReal(lb) - RatApproxReal(solval), 0.0) + MAX(RatApproxReal(solval) - RatApproxReal(ub), 0.0));
               }
#ifdef SCIP_DEBUG
               else
               {
                  SCIPsetDebugMsgPrint(set, "  -> solution value %g violates bounds of <%s>[%g,%g]\n", RgetRealApprox(solval), SCIPvarGetName(var),
                        SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
               }
#endif
            }

            /* check whether there are infinite variable values that lead to an objective value of +infinity */
            if( *feasible && sol->hasinfval )
            {
               *feasible = *feasible && (!RatIsInfinity(solval) || !RatIsPositive(SCIPvarGetObjExact(var)) );
               *feasible = *feasible && (!RatIsNegInfinity(solval) || !RatIsNegative(SCIPvarGetObjExact(var)) );

               if( ((RatIsInfinity(solval) && RatIsPositive(SCIPvarGetObjExact(var)))) 
                     || (RatIsNegInfinity(solval) && RatIsNegative(SCIPvarGetObjExact(var))) )
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
                        RgetRealApprox(solval), SCIPvarGetName(var), SCIPvarGetUnchangedObj(var));
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
      /** @todo: exip turn this into checkExact new callback */
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

   RatFreeBuffer(set->buffer, &solval);

   return SCIP_OKAY;
}

#if 0

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
   assert(!sol->solorigin = SCIP_SOLORIGIN_ORIGINAL);
   assert(prob != NULL);
   assert(prob->transformed);
   assert(success != NULL);

   /* round all roundable fractional variables in the corresponding direction as long as no unroundable var was found */
   nvars = prob->nbinvars + prob->nintvars;
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

/** updates the solution value sums in variables by adding the value in the given solution */
void SCIPsolUpdateVarsum(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Rational*        weight              /**< weight of solution in weighted average */
   )
{
   SCIP_Real solval;
   int v;

   assert(sol != NULL);
   assert(!sol->solorigin = SCIP_SOLORIGIN_ORIGINAL);
   assert(0.0 <= weight && weight <= 1.0);

   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);
      SCIPsolexGetVal(solval, sol, set, stat, prob->vars[v]);
      if( solval != SCIP_UNKNOWN ) /*lint !e777*/
      {
         prob->vars[v]->primsolavg *= (1.0-weight);
         prob->vars[v]->primsolavg += weight*solval;
      }
   }
}
#endif

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

   /** @todo exip: once we have exact presolving, we need to do here what we do in the fp case */

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

#if 0
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
   assert(sol->solorigin = SCIP_SOLORIGIN_ORIGINAL);
   assert(origprob != NULL);

   vars = origprob->vars;
   nvars = origprob->nvars;

   /* recompute the objective value */
   sol->obj = SCIPprobGetObjoffset(origprob);
   for( v = 0; v < nvars; ++v )
   {
      SCIPsolexGetVal(solval, sol, set, stat, vars[v]);
      if( !RatIsZero(solval) && solval != SCIP_UNKNOWN ) /*lint !e777*/
      {
         sol->obj += SCIPvarGetUnchangedObj(vars[v]) * solval;
      }
   }

   if( SCIPsetIsInfinity(set, -sol->obj) )
      sol->obj = -SCIPsetInfinity(set);
}
#endif

/** returns whether the given solutions are equal */
SCIP_Bool SCIPsolexsAreEqual(
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
   SCIP_Rational* tmp1;
   SCIP_Rational* tmp2;
   int v;
   SCIP_Bool result = TRUE;

   assert(sol1 != NULL);
   assert(sol2 != NULL);
   assert((sol1->solorigin == SCIP_SOLORIGIN_ORIGINAL) && (sol2->solorigin == SCIP_SOLORIGIN_ORIGINAL) || transprob != NULL);

   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp1) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp2) );
   /* if both solutions are original or both are transformed, take the objective values stored in the solutions */
   if( (sol1->solorigin == SCIP_SOLORIGIN_ORIGINAL) == (sol2->solorigin == SCIP_SOLORIGIN_ORIGINAL) )
   {
      RatSet(tmp1, sol1->valsex->obj);
      RatSet(tmp2, sol2->valsex->obj);
   }
   /* one solution is original and the other not, so we have to get for both the objective in the transformed problem */
   else
   {
      RatSet(tmp1, SCIPsolexGetObj(sol1, set, transprob, origprob));
      RatSet(tmp2, SCIPsolexGetObj(sol2, set, transprob, origprob));
   }

   /* solutions with different objective values cannot be the same */
   if( !RatIsEqual(tmp1, tmp2) )
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
      SCIPsolexGetVal(tmp1, sol1, set, stat, prob->vars[v]);
      SCIPsolexGetVal(tmp2, sol2, set, stat, prob->vars[v]);

      if( !RatIsEqual(tmp1, tmp2) )
         result = FALSE;
   }

   RatFreeBuffer(set->buffer, &tmp2);
   RatFreeBuffer(set->buffer, &tmp1);

   return result;
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

/** hard-set the obj value of a solution  */
void SCIPsolSetObjVal(
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Real             val                 /**< objective value */
   )
{
   assert(sol != NULL);

   sol->obj = val;
}

#if 0
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
   assert(sol->solorigin = SCIP_SOLORIGIN_ORIGINAL || prob->transformed || transprob != NULL);

   /* display variables of problem data */
   for( v = 0; v < prob->nfixedvars; ++v )
   {
      assert(prob->fixedvars[v] != NULL);
      solval = SCIPsolGetRayVal(sol, set, stat, prob->fixedvars[v]);
      if( printzeros || !RatIsZero(solval) )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(prob->fixedvars[v]));
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
         else if( SCIPsetIsInfinity(set, solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
         else if( SCIPsetIsInfinity(set, -solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " % 20.15g", solval);
         SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(prob->fixedvars[v]));
      }
   }
   for( v = 0; v < prob->nvars; ++v )
   {
      assert(prob->vars[v] != NULL);
      solval = SCIPsolGetRayVal(sol, set, stat, prob->vars[v]);
      if( printzeros || !RatIsZero(solval) )
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
   if( !prob->transformed && !sol->solorigin = SCIP_SOLORIGIN_ORIGINAL )
   {
      assert(transprob != NULL);
      for( v = 0; v < transprob->nfixedvars; ++v )
      {
         assert(transprob->fixedvars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->fixedvars[v]) )
            continue;

         solval = SCIPsolGetRayVal(sol, set, stat, transprob->fixedvars[v]);
         if( printzeros || !RatIsZero(solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->fixedvars[v]));
            if( solval == SCIP_UNKNOWN ) /*lint !e777*/
               SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
            else if( SCIPsetIsInfinity(set, solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " % 20.15g", solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->fixedvars[v]));
         }
      }
      for( v = 0; v < transprob->nvars; ++v )
      {
         assert(transprob->vars[v] != NULL);
         if( SCIPvarIsTransformedOrigvar(transprob->vars[v]) )
            continue;

         solval = SCIPsolGetRayVal(sol, set, stat, transprob->vars[v]);
         if( printzeros || !RatIsZero(solval) )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPvarGetName(transprob->vars[v]));
            if( solval == SCIP_UNKNOWN ) /*lint !e777*/
               SCIPmessageFPrintInfo(messagehdlr, file, "              unknown");
            else if( SCIPsetIsInfinity(set, solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity");
            else if( SCIPsetIsInfinity(set, -solval) )
               SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity");
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " % 20.15g", solval);
            SCIPmessageFPrintInfo(messagehdlr, file, " \t(obj:%.15g)\n", SCIPvarGetUnchangedObj(transprob->vars[v]));
         }
      }
   }

   return SCIP_OKAY;
}
#endif

SCIP_Bool SCIPsolIsExactSol(
   SCIP_SOL*             sol                /**< primal CIP solution */
   )
{
   assert(sol != NULL);

   return sol->valsex != NULL;
}