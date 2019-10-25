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
   if( !set->misc_exactsolve )
      return SCIP_OKAY;

   assert(var != NULL);
   assert(var->exactdata->excol == NULL);
   assert(var->scip == set->scip);
   assert(var->exactdata != NULL);

   SCIPsetDebugMsg(set, "creating exact column for variable <%s>\n", var->name);


   /* switch variable status */
   var->exactdata->exvarstatus = SCIP_VARSTATUS_COLUMN; /*lint !e641*/

   /* create column of variable */
   SCIP_CALL( SCIPcolexCreate(&(var->exactdata->excol), SCIPvarGetCol(var), blkmem, set, stat, var, 0, NULL, NULL, var->removable) );

   if( var->probindex != -1 )
   {
      /* inform LP, that problem variable is now a column variable and no longer loose */
      SCIP_CALL( SCIPlpexUpdateVarColumn(lp, set, var) );
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
   assert(!RatIsAbsInfinity(val));

   //SCIPsetDebugMsg(set, "adding coefficient %g<%s> to row <%s>\n", val, var->name, row->name);

   if ( RatIsZero(val) )
      return SCIP_OKAY;

   switch( SCIPvarGetStatusExact(var) )
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
      if( RatIsEqual(var->exactdata->glbdom.lb, var->exactdata->glbdom.ub) )
      {
         SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );
         RatMult(tmp, val, var->exactdata->glbdom.lb);
         SCIP_CALL( SCIProwexAddConstant(rowex, blkmem, set, stat, eventqueue, lpex, tmp) );
         RatFreeBuffer(set->buffer, &tmp);
         break;;
      }
      /* convert loose variable into column */
      SCIP_CALL( SCIPvarColumnExact(var, blkmem, set, stat, prob, lpex) );
      assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);
      /*lint -fallthrough*/

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      assert(var->data.col->var == var);
      SCIP_CALL( SCIProwexIncCoef(rowex, blkmem, set, eventqueue, lpex, var->exactdata->excol, val) );
      break;

   case SCIP_VARSTATUS_FIXED:
      assert(RatIsEqual(var->exactdata->glbdom.lb, var->exactdata->glbdom.ub)); /*lint !e777*/
      assert(RatIsEqual(var->exactdata->locdom.lb, var->exactdata->locdom.ub)); /*lint !e777*/
      assert(RatIsEqual(var->exactdata->locdom.lb, var->exactdata->glbdom.lb)); /*lint !e777*/
      assert(!RatIsAbsInfinity(var->exactdata->locdom.lb));

      SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );

      RatMult(tmp, val, var->exactdata->locdom.lb);
      SCIP_CALL( SCIProwexAddConstant(rowex, blkmem, set, stat, eventqueue, lpex, tmp) );

      RatFreeBuffer(set->buffer, &tmp);

      break;

   case SCIP_VARSTATUS_AGGREGATED:
      SCIPerrorMessage("aggregated variables not implemented in exact mode yet \n");
      return SCIP_ERROR;

   case SCIP_VARSTATUS_MULTAGGR:
      SCIPerrorMessage("aggregated variables not implemented in exact mode yet \n");
      return SCIP_ERROR;

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatusExact(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);

      SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );

      RatNegate(tmp, val);
      SCIP_CALL( SCIPvarAddToRowExact(var->negatedvar, blkmem, set, stat, eventqueue, prob, lpex, rowex, tmp) );

      RatMultReal(tmp, val, var->data.negate.constant);
      SCIP_CALL( SCIProwexAddConstant(rowex, blkmem, set, stat, eventqueue, lpex, tmp) );

      RatFreeBuffer(set->buffer, &tmp);

      break;

   default:
      SCIPerrorMessage("unknown variable status\n");
      return SCIP_INVALIDDATA;
   }

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

   if( SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_ORIGINAL )
      return var->exactdata->origdom.lb;
   else
   {
      assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatusExact(var->negatedvar) == SCIP_VARSTATUS_ORIGINAL);

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

   if( SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_ORIGINAL )
      return var->exactdata->origdom.ub;
   else
   {
      assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatusExact(var->negatedvar) == SCIP_VARSTATUS_ORIGINAL);

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

   if( RatIsPositive(var->exactdata->obj) || RatIsZero(var->exactdata->obj) )
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

   if( RatIsPositive(var->exactdata->obj) || RatIsZero(var->exactdata->obj) )
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

   if( RatIsPositive(var->exactdata->obj) || RatIsZero(var->exactdata->obj) )
      return var->exactdata->locdom.lb;
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

   if( RatIsPositive(var->exactdata->obj) || RatIsZero(var->exactdata->obj) )
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

   if( RatIsPositive(var->exactdata->obj) || RatIsZero(var->exactdata->obj) )
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

   if( RatIsPositive(var->exactdata->obj) || RatIsZero(var->exactdata->obj) )
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
         if( SCIPvarGetStatusExact(*var) == SCIP_VARSTATUS_NEGATED &&
            ((*var)->negatedvar->nparentvars == 0 || (*var)->negatedvar->parentvars[0] != *var) )
         {
            RatNegate(scalar, scalar);
            tmp = RatApproxReal(scalar) * (*var)->data.negate.constant;
            RatDiffReal(constant, constant, tmp);
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

      switch( SCIPvarGetStatusExact(parentvar) )
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
         RatDiffReal(scalar, scalar, parentvar->data.aggregate.scalar);
         tmp = RatApproxReal(scalar) * parentvar->data.aggregate.constant;
         RatDiffReal(constant, constant, tmp);
         break;

      case SCIP_VARSTATUS_NEGATED: /* x = b - y  ->  y = b - x,  s*y + c = -s*x + c+b*s */
         assert(parentvar->negatedvar != NULL);
         assert(SCIPvarGetStatusExact(parentvar->negatedvar) != SCIP_VARSTATUS_NEGATED);
         assert(parentvar->negatedvar->negatedvar == parentvar);
         RatNegate(scalar, scalar);
         tmp = RatApproxReal(scalar) * parentvar->data.negate.constant;
         RatDiffReal(constant, constant, tmp);
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

/** gets column of COLUMN variable */
SCIP_COLEX* SCIPvarGetColExact(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);
   assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN);

   return var->exactdata->excol;
}

/** gets primal LP solution value of variable */
SCIP_Rational* SCIPvarGetLPexSolex_rec(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   switch( SCIPvarGetStatusExact(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return NULL;
      return SCIPvarGetLPexSolex(var->data.original.transvar);

   case SCIP_VARSTATUS_LOOSE:
      return SCIPvarGetBestBoundLocalExact(var);

   case SCIP_VARSTATUS_COLUMN:
      assert(var->data.col != NULL);
      return SCIPcolexGetPrimsol(var->exactdata->excol);

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->exactdata->locdom.lb;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }
}

/** gets pseudo solution value of variable at current node */
static
SCIP_Rational* SCIPvarGetPseudoSolex_rec(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real pseudosol;
   int i;

   assert(var != NULL);

   switch( SCIPvarGetStatusExact(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return NULL;
      return SCIPvarGetPseudoSolex(var->data.original.transvar);

   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_COLUMN:
      return SCIPvarGetBestBoundLocalExact(var);

   case SCIP_VARSTATUS_FIXED:
      assert(var->locdom.lb == var->locdom.ub); /*lint !e777*/
      return var->exactdata->locdom.lb;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }
}


/** gets primal LP solution value of variable */
SCIP_Rational* SCIPvarGetLPexSolex(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      return SCIPcolexGetPrimsol(var->exactdata->excol);
   else
      return SCIPvarGetLPexSolex_rec(var);
}

/** gets pseudo solution value of variable */
SCIP_Rational* SCIPvarGetPseudoSolex(
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(var != NULL);

   if( SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN )
      return SCIPvarGetBestBoundLocalExact(var);
   else
      return SCIPvarGetPseudoSolex_rec(var);
}

/** gets current LP or pseudo solution value of variable */
SCIP_Rational* SCIPvarGetSolex(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             getlpval            /**< should the LP solution value be returned? */
   )
{
   if( getlpval )
      return SCIPvarGetLPexSolex(var);
   else
      return SCIPvarGetPseudoSolex(var);
}

/** adds correct bound-data to negated variable */
SCIP_RETCODE SCIPvarNegateExactData(
   SCIP_VAR*             negvar,             /**< the negated variable */
   SCIP_VAR*             origvar,            /**< the original variable */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_Real constant;

   if( !set->misc_exactsolve )
      return SCIP_OKAY;

   assert(negvar != NULL);
   assert(origvar != NULL);
   assert(origvar->exactdata != NULL);
   assert(negvar->exactdata == NULL);

   constant = negvar->data.negate.constant;

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &(negvar->exactdata)) );

   SCIP_CALL( RatCreateBlock(blkmem, &negvar->exactdata->glbdom.ub) );
   SCIP_CALL( RatCreateBlock(blkmem, &negvar->exactdata->glbdom.lb) );

   SCIP_CALL( RatCreateBlock(blkmem, &negvar->exactdata->glbdom.lb) );
   SCIP_CALL( RatCreateBlock(blkmem, &negvar->exactdata->glbdom.ub) );

   SCIP_CALL( RatCreateBlock(blkmem, &negvar->exactdata->obj) );

   RatDiffReal(negvar->exactdata->glbdom.ub, origvar->exactdata->glbdom.lb, constant);
   RatNegate(negvar->exactdata->glbdom.ub, negvar->exactdata->glbdom.ub);

   RatDiffReal(negvar->exactdata->glbdom.lb, origvar->exactdata->glbdom.ub, constant);
   RatNegate(negvar->exactdata->glbdom.lb, negvar->exactdata->glbdom.lb);

   RatDiffReal(negvar->exactdata->locdom.ub, origvar->exactdata->locdom.lb, constant);
   RatNegate(negvar->exactdata->locdom.ub, negvar->exactdata->locdom.ub);

   RatDiffReal(negvar->exactdata->locdom.lb, origvar->exactdata->locdom.ub, constant);
   RatNegate(negvar->exactdata->locdom.lb, negvar->exactdata->locdom.lb);

   negvar->exactdata->exvarstatus = SCIP_VARSTATUS_NEGATED;

   assert(RatIsEqualReal(negvar->exactdata->glbdom.ub, negvar->glbdom.ub));
   assert(RatIsEqualReal(negvar->exactdata->locdom.ub, negvar->locdom.ub));
   assert(RatIsEqualReal(negvar->exactdata->glbdom.lb, negvar->glbdom.lb));
   assert(RatIsEqualReal(negvar->exactdata->locdom.lb, negvar->locdom.lb));

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

   assert(var->glbdom.ub == RatApproxReal(ub));
   assert(var->glbdom.lb == RatApproxReal(lb));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &(var->exactdata)) );

   SCIP_CALL( RatCopy(blkmem, &var->exactdata->glbdom.lb, lb) );
   SCIP_CALL( RatCopy(blkmem, &var->exactdata->glbdom.ub, ub) );

   SCIP_CALL( RatCopy(blkmem, &var->exactdata->locdom.lb, lb) );
   SCIP_CALL( RatCopy(blkmem, &var->exactdata->locdom.ub, ub) );

   var->exactdata->excol = NULL;
   var->exactdata->exvarstatus = var->varstatus;

   if( obj != NULL )
      SCIP_CALL( RatCopy(blkmem, &var->exactdata->obj, obj) );
   else
      SCIP_CALL( RatCreateBlock(blkmem, &var->exactdata->obj) );

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

   SCIP_CALL( RatCopy(blkmem, &targetvar->exactdata->glbdom.lb, sourcevar->exactdata->glbdom.lb) );
   SCIP_CALL( RatCopy(blkmem, &targetvar->exactdata->glbdom.ub, sourcevar->exactdata->glbdom.ub) );
   SCIP_CALL( RatCopy(blkmem, &targetvar->exactdata->locdom.lb, sourcevar->exactdata->locdom.lb) );
   SCIP_CALL( RatCopy(blkmem, &targetvar->exactdata->locdom.ub, sourcevar->exactdata->locdom.ub) );
   SCIP_CALL( RatCopy(blkmem, &targetvar->exactdata->obj, sourcevar->exactdata->obj) );
   targetvar->exactdata->excol = NULL;
   targetvar->exactdata->exvarstatus = SCIP_VARSTATUS_LOOSE;

   return SCIP_OKAY;
}

/** free exact variable data, if it exists */
SCIP_RETCODE SCIPvarFreeExactData(
   SCIP_VAR*             var,                /**< variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue (may be NULL, if it's not a column variable) */
   SCIP_LPEX*            lp                  /**< current LP data (may be NULL, if it's not a column variable) */
   )
{
   if( set->misc_exactsolve )
   {
      if( SCIPvarGetStatusExact(var) ==  SCIP_VARSTATUS_COLUMN )
      {
         SCIP_CALL( SCIPcolexFree(&(var->exactdata->excol), blkmem, set, eventqueue, lp) );
      }
      /* free exact variable data if it was created */
      if( var->exactdata != NULL )
      {
         RatFreeBlock(blkmem, &(var)->exactdata->glbdom.lb);
         RatFreeBlock(blkmem, &(var)->exactdata->glbdom.ub);
         RatFreeBlock(blkmem, &(var)->exactdata->locdom.lb);
         RatFreeBlock(blkmem, &(var)->exactdata->locdom.ub);
         RatFreeBlock(blkmem, &(var)->exactdata->obj );

         BMSfreeBlockMemory(blkmem, &(var)->exactdata);
         assert((var)->exactdata == NULL);
      }
   }
   else
   {
      assert(var->exactdata == NULL);
   }
   
   return SCIP_OKAY;
}

/* @todo exip: this is currently a blank */
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
   assert(SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatusExact(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarIsTransformed(var));
   assert(!SCIPsetIsEQ(set, oldobj, newobj));

   //SCIP_CALL( SCIPeventCreateObjChanged(&event, blkmem, var, oldobj, newobj) );
   //SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, primal, lp, NULL, NULL, &event) );

   return SCIP_OKAY;
}

/** changes objective value of variable */
SCIP_RETCODE SCIPvarScaleObjExact(
   SCIP_VAR*             var,                /**< variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             scale               /**< new objective value for variable */
   )
{
   SCIP_Rational* tmp;

   assert(var != NULL);
   assert(set != NULL);

   if( !set->misc_exactsolve )
      return SCIP_OKAY;

   assert(var->exactdata != NULL);
   assert(var->scip == set->scip);

   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );

   RatMultReal(tmp, SCIPvarGetObjExact(var), scale);

   SCIP_CALL( SCIPvarChgObjExact(var, blkmem, set, prob, primal, lp, eventqueue, tmp) );
   RatFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** changes objective value of variable */
SCIP_RETCODE SCIPvarChgObjExact(
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
   SCIP_Rational* tmp;

   assert(var != NULL);
   assert(set != NULL);

   if( !set->misc_exactsolve )
      return SCIP_OKAY;

   assert(var->exactdata != NULL);
   assert(var->scip == set->scip);

   SCIP_CALL( RatCreateBlock(blkmem, &tmp) );
   newobjreal = RatApproxReal(newobj);

   SCIPsetDebugMsg(set, "changing exact objective value of <%s> from %g to %g\n", var->name, var->obj, newobjreal);

   if( !RatIsEqual(var->exactdata->obj, newobj) )
   {
      switch( SCIPvarGetStatusExact(var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         if( var->data.original.transvar != NULL )
         {
            assert(SCIPprobIsTransformed(prob));

            RatMultReal(tmp, newobj, (SCIP_Real) prob->objsense/prob->objscale);

            SCIP_CALL( SCIPvarChgObjExact(var->data.original.transvar, blkmem, set, prob, primal, lp, eventqueue,
                  tmp) );
         }
         else
            assert(set->stage == SCIP_STAGE_PROBLEM);

         RatSet(var->exactdata->obj, newobj);

         break;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
         oldobj = var->obj;
         RatSet(var->exactdata->obj, newobj);

         /* @todo exip SCIP_CALL( varEventObjChanged(var, blkmem, set, primal, lp->fplp, eventqueue, oldobj, var->obj) ); */
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

   RatFreeBlock(blkmem, &tmp);

   return SCIP_OKAY;
}

/** return the status of the exact variable data */
SCIP_VARSTATUS SCIPvarGetStatusExact(
   SCIP_VAR*             var                /**< scip variabel */
   )
{
   assert(var != NULL);
   assert(var->exactdata != NULL);

   return var->exactdata->exvarstatus;
}