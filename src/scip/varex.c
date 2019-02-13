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

/**@file   varex.c
 * @brief  methods for problem variables
 * @author Leon Eifler
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#include "scip/cons.h"
#include "scip/event.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/lp.h"
#include "scip/lpex.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/pub_cons.h"
#include "scip/pub_history.h"
#include "scip/pub_implics.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_prop.h"
#include "scip/pub_var.h"
#include "scip/pub_varex.h"
#include "scip/relax.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/stat.h"
#include "scip/struct_event.h"
#include "scip/struct_lpex.h"
#include "scip/struct_prob.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_var.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/varex.h"
#include <string.h>

/** gets objective function value of variable */
SCIP_Rational* SCIPvarGetObjExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);

   return var->exactdata->obj;
}

/** converts loose transformed variable into column variable, creates LP column */
SCIP_RETCODE SCIPvarColumnExact(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LPEX*            lp                /**< current LP data */
   )
{
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(var->scip == set->scip);
   assert(var->exactdata != NULL);

   SCIPsetDebugMsg(set, "creating column for variable <%s>\n", var->name);

   /* switch variable status */
   var->varstatus = SCIP_VARSTATUS_COLUMN; /*lint !e641*/

   /* create column of variable */
   SCIP_CALL( SCIPcolCreate(&var->data.col, blkmem, set, stat, var, 0, NULL, NULL, var->removable) );
   SCIP_CALL( SCIPcolexCreate(&var->exactdata->excol, var->data.col, blkmem, set, stat, var, 0, NULL, NULL, var->removable) );


   if( var->probindex != -1 )
   {
      /* inform problem about the variable's status change */
      SCIP_CALL( SCIPprobVarChangedStatus(prob, blkmem, set, NULL, NULL, var) );

      /* inform LP, that problem variable is now a column variable and no longer loose */
      SCIP_CALL( SCIPlpexUpdateVarColumn(lp, set, var) );
      /* todo exip: implement this */
      /* SCIP_CALL( SCIPlpUpdateVarColumn(lp, set, var) ); */
   }

   return SCIP_OKAY;
}


/** resolves variable to columns and adds them with the coefficient to the
 *
 */
SCIP_RETCODE SCIPvarAddToRowExact(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LPEX*            lpex,               /**< current LP data */
   SCIP_ROWEX*           rowex,              /**< LP row */
   SCIP_Rational*        val                 /**< value of coefficient */
   )
{
   SCIP_Rational* tmp;
   int i;

   assert(var != NULL);
   assert(set != NULL);
   assert(var->scip == set->scip);
   assert(rowex != NULL);
   assert(!RisAbsInfinity(val));

   //SCIPsetDebugMsg(set, "adding coefficient %g<%s> to row <%s>\n", val, var->name, row->name);

   if ( RisZero(val) )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
      {
         SCIPerrorMessage("cannot add untransformed original variable <%s> to excact LP row <%s>\n", var->name, rowex->fprow->name);
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPvarAddToRowExact(var->data.original.transvar, blkmem, set, stat, eventqueue, prob, lpex, rowex, val) );
      break;;

   case SCIP_VARSTATUS_LOOSE:
      /* add globally fixed variables as constant */
      if( RisEqual(var->exactdata->glbdom.lb, var->exactdata->glbdom.ub) )
      {
         tmp = RcreateTemp(set->buffer);
         Rmult(tmp, val, var->exactdata->glbdom.lb);
         SCIP_CALL( SCIProwexAddConstant(rowex, blkmem, set, stat, eventqueue, lpex, tmp) );
         RdeleteTemp(set->buffer, &tmp);
         break;;
      }
      /* convert loose variable into column */
      SCIP_CALL( SCIPvarColumnExact(var, blkmem, set, stat, prob, lpex) );
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      /*lint -fallthrough*/

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      assert(var->data.col->var == var);
      SCIP_CALL( SCIProwexIncCoef(rowex, blkmem, set, eventqueue, lpex, var->exactdata->excol, val) );
      break;

   case SCIP_VARSTATUS_FIXED:
      assert(RisEqual(var->exactdata->glbdom.lb, var->exactdata->glbdom.ub)); /*lint !e777*/
      assert(RisEqual(var->exactdata->locdom.lb, var->exactdata->locdom.ub)); /*lint !e777*/
      assert(RisEqual(var->exactdata->locdom.lb, var->exactdata->glbdom.lb)); /*lint !e777*/
      assert(!RisAbsInfinity(var->exactdata->locdom.lb));

      tmp = RcreateTemp(set->buffer);

      Rmult(tmp, val, var->exactdata->locdom.lb);
      SCIP_CALL( SCIProwexAddConstant(rowex, blkmem, set, stat, eventqueue, lpex, tmp) );

      RdeleteTemp(set->buffer, &tmp);

      break;

   case SCIP_VARSTATUS_AGGREGATED:
      SCIPerrorMessage("aggregated variables not implemented in exact mode yet \n");
      return SCIP_ERROR;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("aggregated variables not implemented in exact mode yet \n");
      return SCIP_ERROR;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);

      tmp = RcreateTemp(set->buffer);

      Rneg(tmp, val);
      SCIP_CALL( SCIPvarAddToRowExact(var->negatedvar, blkmem, set, stat, eventqueue, prob, lpex, rowex, tmp) );

      RmultReal(tmp, val, var->data.negate.constant);
      SCIP_CALL( SCIProwexAddConstant(rowex, blkmem, set, stat, eventqueue, lpex, tmp) );

      RdeleteTemp(set->buffer, &tmp);

      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }
   //SCIP_CALL( SCIPvarAddToRow(var, blkmem, set, stat, eventqueue, prob, lpex->fplp, rowex->fprow, RgetRealApprox(val)) );

   return SCIP_OKAY;
}


/*
 * getter methods for exact varbounds
 */

/** gets original lower bound of original problem variable (i.e. the bound set in problem creation) */
SCIP_Rational* SCIPvarGetLbOriginalExact(
   SCIP_VAR*             var                 /**< original problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsOriginal(var));
   assert(var->exactdata != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      return var->exactdata->origdom.lb;
   else
   {
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) == SCIP_VARSTATUS_ORIGINAL);

      SCIPerrorMessage("negated var not implemented yet for rational data \n");
      SCIPABORT();
   }
}

/** gets original upper bound of original problem variable (i.e. the bound set in problem creation) */
SCIP_Rational* SCIPvarGetUbOriginalExact(
   SCIP_VAR*             var                 /**< original problem variable */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsOriginal(var));
   assert(var->exactdata != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      return var->exactdata->origdom.ub;
   else
   {
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) == SCIP_VARSTATUS_ORIGINAL);

      SCIPerrorMessage("negated var not implemented yet for rational data \n");
      SCIPABORT();
   }
}

/** gets global lower bound of variable */
SCIP_Rational* SCIPvarGetLbGlobalExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(var->exactdata->glbdom.lb != NULL);

   return var->exactdata->glbdom.lb;
}

/** gets global upper bound of variable */
SCIP_Rational* SCIPvarGetUbGlobalExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(var->exactdata->glbdom.ub != NULL);

   return var->exactdata->glbdom.ub;
}

/** gets best global bound of variable with respect to the objective function */
SCIP_Rational* SCIPvarGetBestBoundGlobalExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(var->exactdata->glbdom.lb != NULL);
   assert(var->exactdata->glbdom.ub != NULL);
   assert(var->exactdata->obj != NULL);

   if( RisPositive(var->exactdata->obj) || RisZero(var->exactdata->obj) )
      return var->exactdata->glbdom.lb;
   else
      return var->exactdata->glbdom.ub;
}

/** gets worst global bound of variable with respect to the objective function */
SCIP_Rational* SCIPvarGetWorstBoundGlobalExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(var->exactdata->glbdom.lb != NULL);
   assert(var->exactdata->glbdom.ub != NULL);
   assert(var->exactdata->obj != NULL);

   if( RisPositive(var->exactdata->obj) || RisZero(var->exactdata->obj) )
      return var->exactdata->glbdom.ub;
   else
      return var->exactdata->glbdom.lb;
}

/** gets current lower bound of variable */
SCIP_Rational* SCIPvarGetLbLocalExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(var->exactdata->locdom.lb != NULL);

   return var->exactdata->locdom.lb;
}

/** gets current upper bound of variable */
SCIP_Rational* SCIPvarGetUbLocalExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(var->exactdata->locdom.ub != NULL);

   return var->exactdata->locdom.ub;
}

/** gets best local bound of variable with respect to the objective function */
SCIP_Rational* SCIPvarGetBestBoundLocalExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(var->exactdata->locdom.lb != NULL);
   assert(var->exactdata->locdom.ub != NULL);
   assert(var->exactdata->obj != NULL);

   if( RisPositive(var->exactdata->obj) || RisZero(var->exactdata->obj) )
      return var->exactdata->locdom.ub;
   else
      return var->exactdata->locdom.ub;
}

/** gets worst local bound of variable with respect to the objective function */
SCIP_Rational* SCIPvarGetWorstBoundLocalExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(var->exactdata->locdom.lb != NULL);
   assert(var->exactdata->locdom.ub != NULL);
   assert(var->exactdata->obj != NULL);

   if( RisPositive(var->exactdata->obj) || RisZero(var->exactdata->obj) )
      return var->exactdata->locdom.ub;
   else
      return var->exactdata->locdom.lb;
}

/** gets type (lower or upper) of best bound of variable with respect to the objective function */
SCIP_BOUNDTYPE SCIPvarGetBestBoundTypeExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(var->exactdata->obj != NULL);

   if( RisPositive(var->exactdata->obj) || RisZero(var->exactdata->obj) )
      return SCIP_BOUNDTYPE_LOWER;
   else
      return SCIP_BOUNDTYPE_UPPER;
}

/** gets type (lower or upper) of worst bound of variable with respect to the objective function */
SCIP_BOUNDTYPE SCIPvarGetWorstBoundTypeExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(var->exactdata->obj != NULL);

   if( RisPositive(var->exactdata->obj) || RisZero(var->exactdata->obj) )
      return SCIP_BOUNDTYPE_UPPER;
   else
      return SCIP_BOUNDTYPE_LOWER;
}

/** todo exip: this does not work exactly yet */
/** retransforms given variable, scalar anqd constant to the corresponding original variable, scalar
 *  and constant, if possible; if the retransformation is impossible, NULL is returned as variable
 */
SCIP_RETCODE SCIPvarGetOrigvarSumExact(
   SCIP_VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   SCIP_Rational*        scalar,             /**< pointer to scalar a in sum a*x + c */
   SCIP_Rational*        constant            /**< pointer to constant c in sum a*x + c */
   )
{
   SCIP_VAR* parentvar;
   SCIP_Real tmp;

   assert(var != NULL);
   assert(*var != NULL);
   assert(scalar != NULL);
   assert(constant != NULL);

   while( !SCIPvarIsOriginal(*var) )
   {
      /* if the variable has no parent variables, it was generated during solving and has no corresponding original
       * var
       */
      if( (*var)->nparentvars == 0 )
      {
         /* negated variables do not need to have a parent variables, and negated variables can exist in original
          * space
          */
         if( SCIPvarGetStatus(*var) == SCIP_VARSTATUS_NEGATED &&
            ((*var)->negatedvar->nparentvars == 0 || (*var)->negatedvar->parentvars[0] != *var) )
         {
            Rneg(scalar, scalar);
            tmp = RgetRealApprox(scalar) * (*var)->data.negate.constant;
            RdiffReal(constant, constant, tmp);
            *var = (*var)->negatedvar;

            continue;
         }
         /* if the variables does not have any parent the variables was created during solving and has no original
          * counterpart
          */
         else
         {
            *var = NULL;

            return SCIP_OKAY;
         }
      }

      /* follow the link to the first parent variable */
      parentvar = (*var)->parentvars[0];
      assert(parentvar != NULL);

      switch( SCIPvarGetStatus(parentvar) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         break;

      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_MULTAGGR:
         SCIPerrorMessage("column, loose, fixed or multi-aggregated variable cannot be the parent of a variable\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + b  ->  y = (x-b)/a,  s*y + c = (s/a)*x + c-b*s/a */
         assert(parentvar->data.aggregate.var == *var);
         assert(parentvar->data.aggregate.scalar != 0.0);
         RdiffReal(scalar, scalar, parentvar->data.aggregate.scalar);
         tmp = RgetRealApprox(scalar) * parentvar->data.aggregate.constant;
         RdiffReal(constant, constant, tmp);
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = b - y  ->  y = b - x,  s*y + c = -s*x + c+b*s */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatus(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         Rneg(scalar, scalar);
         tmp = RgetRealApprox(scalar) * parentvar->data.negate.constant;
         RdiffReal(constant, constant, tmp);
         break;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }

      assert( parentvar != NULL );
      *var = parentvar;
   }

   return SCIP_OKAY;
}

/** create and set the exact variable bounds and objective value */
SCIP_RETCODE SCIPvarAddExactData(
   SCIP_VAR*             var,                /**< pointer to variable data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Rational*        lb,                 /**< lower bound of variable */
   SCIP_Rational*        ub,                 /**< upper bound of variable */
   SCIP_Rational*        obj                 /**< objective function value */
   )
{
   assert(var != NULL);
   assert(blkmem != NULL);

   assert(var->glbdom.ub == RgetRealApprox(ub));
   assert(var->glbdom.lb == RgetRealApprox(lb));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &(var->exactdata)) );

   var->exactdata->glbdom.lb = Rcopy(blkmem, lb);
   var->exactdata->glbdom.ub = Rcopy(blkmem, ub);

   var->exactdata->locdom.lb = Rcopy(blkmem, lb);
   var->exactdata->locdom.ub = Rcopy(blkmem, ub);

   if( obj != NULL )
      var->exactdata->obj = Rcopy(blkmem, obj);
   else
      var->exactdata->obj = Rcreate(blkmem);

   return SCIP_OKAY;
}

/** copy exact variable data from one variable to another */
SCIP_RETCODE SCIPvarCopyExactData(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             targetvar,          /**< variable that gets the exact data */
   SCIP_VAR*             sourcevar           /**< variable the data gets copied from */
   )
{
   assert(blkmem != NULL);
   assert(targetvar != NULL);
   assert(sourcevar != NULL);

   /* todo exip: what should happen here? error or just no copy? */
   if( sourcevar->exactdata == NULL )
   {
      return SCIP_OKAY;
   }
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &(targetvar->exactdata)) );

   targetvar->exactdata->glbdom.lb = Rcopy(blkmem, sourcevar->exactdata->glbdom.lb);
   targetvar->exactdata->glbdom.ub = Rcopy(blkmem, sourcevar->exactdata->glbdom.ub);
   targetvar->exactdata->locdom.lb = Rcopy(blkmem, sourcevar->exactdata->locdom.lb);
   targetvar->exactdata->locdom.ub = Rcopy(blkmem, sourcevar->exactdata->locdom.ub);
   targetvar->exactdata->obj = Rcopy(blkmem, sourcevar->exactdata->obj);
   targetvar->exactdata->excol = NULL;

   return SCIP_OKAY;
}

/** appends OBJCHANGED event to the event queue */
static
SCIP_RETCODE varEventObjChanged(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             oldobj,             /**< old objective value for variable */
   SCIP_Real             newobj              /**< new objective value for variable */
   )
{
   SCIP_EVENT* event;

   assert(var != NULL);
   assert(var->scip == set->scip);
   assert(var->eventfilter != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarIsTransformed(var));
   assert(!SCIPsetIsEQ(set, oldobj, newobj));

   SCIP_CALL( SCIPeventCreateObjChanged(&event, blkmem, var, oldobj, newobj) );
   SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, primal, lp, NULL, NULL, &event) );

   return SCIP_OKAY;
}

/** changes objective value of variable */
SCIP_RETCODE SCIPvarChgExactObj(
   SCIP_VAR*             var,                /**< variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Rational*        newobj              /**< new objective value for variable */
   )
{
   SCIP_Real newobjreal;
   SCIP_Real oldobj;
   SCIP_Rational* tmp = Rcreate(blkmem);

   assert(var != NULL);
   assert(set != NULL);
   assert(var->exactdata != NULL);
   assert(var->scip == set->scip);

   newobjreal = RgetRealApprox(newobj);

   SCIPsetDebugMsg(set, "changing exact objective value of <%s> from %g to %g\n", var->name, var->obj, newobjreal);

   if( !RisEqual(var->exactdata->obj, newobj) )
   {
      switch( SCIPvarGetStatus(var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.original.transvar != NULL )
         {
            assert(SCIPprobIsTransformed(prob));

            RmultReal(tmp, newobj, (SCIP_Real) prob->objsense/prob->objscale);

            SCIP_CALL( SCIPvarChgExactObj(var->data.original.transvar, blkmem, set, prob, primal, lp, eventqueue,
                  tmp) );
         }
         else
            assert(set->stage == SCIP_STAGE_PROBLEM);

         Rset(var->exactdata->obj, newobj);
         var->unchangedobj = newobjreal;
         var->obj = newobjreal;

         break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         oldobj = var->obj;
         Rset(var->exactdata->obj, newobj);
         var->obj = newobjreal;

         /* update unchanged objective value of variable */
         if( !lp->fplp->divingobjchg )
            var->unchangedobj = newobjreal;

         /* update the number of variables with non-zero objective coefficient;
          * we only want to do the update, if the variable is added to the problem;
          * since the objective of inactive variables cannot be changed, this corresponds to probindex != -1
          */
         if( SCIPvarIsActive(var) )
            SCIPprobUpdateNObjVars(prob, set, oldobj, var->obj);

         SCIP_CALL( varEventObjChanged(var, blkmem, set, primal, lp->fplp, eventqueue, oldobj, var->obj) );
         break;

      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_MULTAGGR:
      case SCIP_VARSTATUS_NEGATED:
         SCIPerrorMessage("cannot change objective value of a fixed, aggregated, multi-aggregated, or negated variable\n");
         return SCIP_INVALIDDATA;

      default:
         SCIPerrorMessage("unknown variable status\n");
         return SCIP_INVALIDDATA;
      }
   }

   Rdelete(blkmem, &tmp);

   return SCIP_OKAY;
}