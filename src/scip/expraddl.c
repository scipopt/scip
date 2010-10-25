/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: expraddl.c,v 1.1 2010/10/25 04:27:33 bzfviger Exp $"

/**@file   scip/expraddl.c
 * @brief  additional methods for expressions and expression trees
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 *
 * This file contains methods for handling and manipulating expressions and expression trees
 * that are SCIP specific and thus not included in nlpi/expression.*
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"

#include "nlpi/struct_expr.h"

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) (val >= infty1 ? infty2 : val)

/** returns variables of expression tree */
SCIP_VAR** SCIPexprtreeGetVars(
   SCIP_EXPRTREE*        tree                /**< expression tree */
)
{
   assert(tree != NULL);

   return tree->vars;
}

/** stores array of variables in expression tree */
SCIP_RETCODE SCIPexprtreeSetVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars                /**< variables */
)
{
   assert(tree != NULL);
   assert(vars != NULL || nvars == 0);

   if( nvars == 0 )
   {
      BMSfreeBlockMemoryArray(tree->blkmem, &tree->vars, tree->nvars);
      tree->nvars = 0;
      return SCIP_OKAY;
   }

   if( tree->vars != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(tree->blkmem, &tree->vars, tree->nvars, nvars) );
      BMScopyMemoryArray(tree->vars, vars, nvars);
   }
   else
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(tree->blkmem, &tree->vars, vars, nvars) );
   }

   tree->nvars = nvars;

   return SCIP_OKAY;
}

/** adds variables to the expression tree variables array */
SCIP_RETCODE SCIPexprtreeAddVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars                /**< variables */
)
{
   assert(tree != NULL);
   assert(vars != NULL || nvars == 0);
   assert(tree->vars != NULL || tree->nvars == 0);

   if( nvars == 0 )
      return SCIP_OKAY;

   if( tree->nvars == 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(tree->blkmem, &tree->vars, vars, nvars) );
      tree->nvars = nvars;
      return SCIP_OKAY;
   }

   /* append vars to tree->vars array */
   SCIP_ALLOC( BMSreallocBlockMemoryArray(tree->blkmem, &tree->vars, tree->nvars, tree->nvars + nvars) );
   BMScopyMemoryArray(&tree->vars[tree->nvars], vars, nvars);  /*lint !e866*/
   tree->nvars += nvars;

   return SCIP_OKAY;
}

/** evaluates an expression tree for a primal solution or LP solution */
SCIP_RETCODE SCIPexprtreeEvalSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_SOL*             sol,                /**< a solution, or NULL for current LP solution */
   SCIP_Real*            val                 /**< buffer to store value */
)
{
   SCIP_Real* varvals;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   if( tree->nvars == 0 )
   {
      SCIP_CALL( SCIPexprEval(tree->root, NULL, tree->params, val) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, tree->nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, tree->nvars, tree->vars, varvals) );

   SCIP_CALL( SCIPexprEval(tree->root, varvals, tree->params, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

/** evaluates an expression tree w.r.t. current global bounds */
SCIP_RETCODE SCIPexprtreeEvalIntGlobalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        val                 /**< buffer to store result */
)
{
   SCIP_INTERVAL* varvals;
   int i;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   if( tree->nvars == 0 )
   {
      SCIP_CALL( SCIPexprEvalInt(tree->root, infinity, NULL, tree->params, val) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, tree->nvars) );
   for( i = 0; i < tree->nvars; ++i )
   {
      SCIPintervalSetBounds(&varvals[i],
         -infty2infty(SCIPinfinity(scip), infinity, -SCIPvarGetLbGlobal(tree->vars[i])),  /*lint !e666*/
          infty2infty(SCIPinfinity(scip), infinity,  SCIPvarGetUbGlobal(tree->vars[i]))); /*lint !e666*/
   }

   SCIP_CALL( SCIPexprEvalInt(tree->root, infinity, varvals, tree->params, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

/** evaluates an expression tree w.r.t. current local bounds */
SCIP_RETCODE SCIPexprtreeEvalIntLocalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        val                 /**< buffer to store result */
)
{
   SCIP_INTERVAL* varvals;
   int i;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   if( tree->nvars == 0 )
   {
      SCIP_CALL( SCIPexprEvalInt(tree->root, infinity, NULL, tree->params, val) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, tree->nvars) );
   for( i = 0; i < tree->nvars; ++i )
   {
      SCIPintervalSetBounds(&varvals[i],
         -infty2infty(SCIPinfinity(scip), infinity, -SCIPvarGetLbLocal(tree->vars[i])),  /*lint !e666*/
          infty2infty(SCIPinfinity(scip), infinity,  SCIPvarGetUbLocal(tree->vars[i]))); /*lint !e666*/
   }

   SCIP_CALL( SCIPexprEvalInt(tree->root, infinity, varvals, tree->params, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

/** prints an expression tree using variable names from variables array */
SCIP_RETCODE SCIPexprtreePrintWithNames(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   FILE*                 file                /**< file for printing, or NULL for stdout */
)
{
   const char** varnames;
   int i;

   assert(tree != NULL);

   if( tree->nvars == 0 )
   {
      SCIPexprtreePrint(tree, file, NULL, NULL);
      return SCIP_OKAY;
   }

   assert(tree->vars != NULL);

   SCIP_ALLOC( BMSallocMemoryArray(&varnames, tree->nvars) );
   for( i = 0; i < tree->nvars; ++i )
      varnames[i] = SCIPvarGetName(tree->vars[i]);

   SCIPexprtreePrint(tree, file, varnames, NULL);

   BMSfreeMemoryArray(&varnames);

   return SCIP_OKAY;
}

/** searches the variables array of an expression tree for a variable and returns its position, or -1 if not found
 * Note that this is an O(n) operation!
 */
int SCIPexprtreeFindVar(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_VAR*             var                 /**< variable to search for */
)
{
   int i;

   assert(tree != NULL);
   assert(var  != NULL);

   for( i = 0; i < tree->nvars; ++i )
      if( tree->vars[i] == var )
         return i;

   return -1;
}

/** removes fixed variables from an expression tree, so that at exit all variables are active */
SCIP_RETCODE SCIPexprtreeRemoveFixedVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Bool*            changed             /**< buffer to store whether the tree was changed, i.e., whether there was a fixed variable */
)
{
   SCIP_HASHMAP* varhash;
   int i, j;
   int nvarsold;
   SCIP_Bool havefixedvar;
   SCIP_VAR* var;
   SCIP_Real scalar;
   SCIP_Real constant;
   SCIP_EXPR** replaceexprs;
   int idx;
   SCIP_EXPR* tmp;
   SCIP_VAR* mvar;
   SCIP_Real mscalar;
   int* newpos;
   int offset;

   assert(tree != NULL);
   assert(tree->vars != NULL || tree->nvars == 0);
   assert(changed != NULL);

   *changed = FALSE;

   if( tree->nvars == 0 )
      return SCIP_OKAY;

   /* create hash map from variable to indices in tree->vars and check if there is a nonfixed variable */
   havefixedvar = FALSE;
   SCIP_CALL( SCIPhashmapCreate(&varhash, tree->blkmem, SCIPcalcHashtableSize(5 * tree->nvars)) );
   for( i = 0; i < tree->nvars; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(varhash, (void*)tree->vars[i], (void*)(size_t)i) );
      if( !SCIPvarIsActive(tree->vars[i]) )
         havefixedvar = TRUE;
   }

   if( !havefixedvar )
   {
      /* nothing to do */
      SCIPhashmapFree(&varhash);
      return SCIP_OKAY;
   }

   /* we will do something */
   *changed = TRUE;

   nvarsold = tree->nvars;

   /* array to store expressions that replace a variable expression in the tree */
   SCIP_ALLOC( BMSallocBlockMemoryArray(tree->blkmem, &replaceexprs, nvarsold) );
   BMSclearMemoryArray(replaceexprs, nvarsold);

   /* construct for each nonactive variable an expression that replaces this variable in the tree */
   for( i = 0; i < nvarsold; ++i )
   {
      var = tree->vars[i];

      if( SCIPvarIsActive(tree->vars[i]) )
         continue;

      scalar   = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&var, &scalar, &constant) );

      if( scalar == 0.0 )
      {
         /* variable is fixed, thus replace by constant expression in tree */
         SCIP_CALL( SCIPexprCreate(tree->blkmem, &replaceexprs[i], SCIP_EXPR_CONST, constant) );
         continue;
      }

      if( SCIPvarIsActive(var) )
      {
         /* variable was aggregated or negated, thus replace by scalar * var + constant */
         if( !SCIPhashmapExists(varhash, var) )
         {
            /* var not in tree yet, so add it */
            SCIP_CALL( SCIPexprtreeAddVars(tree, 1, &var) );
            idx = tree->nvars - 1;
            SCIP_CALL( SCIPhashmapInsert(varhash, (void*)var, (void*)(size_t)idx) );
         }
         else
         {
            idx = (int)(size_t) SCIPhashmapGetImage(varhash, (void*)var);
         }
         assert(idx >= 0 && idx < tree->nvars);
         assert(tree->vars[idx] == var);

         SCIP_CALL( SCIPexprCreate(tree->blkmem, &replaceexprs[i], SCIP_EXPR_VARIDX, idx) );
         if( scalar != 1.0 )
         {
            /* multiply by scalar */
            SCIP_CALL( SCIPexprCreate(tree->blkmem, &tmp, SCIP_EXPR_CONST, scalar) );
            SCIP_CALL( SCIPexprCreate(tree->blkmem, &replaceexprs[i], SCIP_EXPR_MUL, replaceexprs[i], tmp) );
         }
         if( constant != 0.0 )
         {
            /* add constant */
            SCIP_CALL( SCIPexprCreate(tree->blkmem, &tmp, SCIP_EXPR_CONST, constant) );
            SCIP_CALL( SCIPexprCreate(tree->blkmem, &replaceexprs[i], SCIP_EXPR_PLUS, replaceexprs[i], tmp) );
         }
         continue;
      }

      {
         SCIP_EXPR** summands;
         int         nsummands;
         int         summandssize;

         assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR );
         /* var is now multiaggregated, thus replace by scalar * (multaggrconst + sum_j multaggrscalar_j*multaggrvar_j) + constant
          * and remember if any of the variables in multaggrvar_j are multiaggregated again */

         /* allocate array for summands (number of aggr. variables, +1 if there is a constant) */
         summandssize = SCIPvarGetMultaggrNVars(var);
         if( constant != 0.0 || SCIPvarGetMultaggrConstant(var) != 0.0 )
            ++summandssize;
         SCIP_ALLOC( BMSallocBlockMemoryArray(tree->blkmem, &summands, summandssize) );
         nsummands = 0;

         /* @todo use SCIP_EXPR_LINEAR instead of SCIP_EXPR_SUM */

         /* linear part
          * turn each variable in SCIPvarGetMultaggrVars(var) into an active or multiaggregated one and add corresponding term to summands */
         for( j = 0; j < SCIPvarGetMultaggrNVars(var); ++j )
         {
            mvar      = SCIPvarGetMultaggrVars(var)[j];
            mscalar   = scalar * SCIPvarGetMultaggrScalars(var)[j];
            SCIP_CALL( SCIPvarGetProbvarSum(&mvar, &mscalar, &constant) );

            /* if variable mvar is fixed, constant has been added to constant and we can continue */
            if( mscalar == 0.0 )
               continue;

            /* add mvar to tree, if not in tree yet */
            if( !SCIPhashmapExists(varhash, mvar) )
            {
               /* var not in tree yet, so add it */
               SCIP_CALL( SCIPexprtreeAddVars(tree, 1, &mvar) );
               idx = tree->nvars - 1;
               SCIP_CALL( SCIPhashmapInsert(varhash, (void*)mvar, (void*)(size_t)idx) );
            }
            else
            {
               idx = (int)(size_t) SCIPhashmapGetImage(varhash, (void*)mvar);
            }
            assert(idx >= 0 && idx < tree->nvars);
            assert(tree->vars[idx] == mvar);

            assert(nsummands < summandssize);
            SCIP_CALL( SCIPexprCreate(tree->blkmem, &summands[nsummands], SCIP_EXPR_VARIDX, idx) );
            if( mscalar != 1.0 )
            {
               /* multiply by scalar */
               SCIP_CALL( SCIPexprCreate(tree->blkmem, &tmp, SCIP_EXPR_CONST, mscalar) );
               SCIP_CALL( SCIPexprCreate(tree->blkmem, &summands[nsummands], SCIP_EXPR_MUL, summands[nsummands], tmp) );
            }
            ++nsummands;
         }

         /* constant part */
         if( constant != 0.0 || SCIPvarGetMultaggrConstant(var) != 0.0 )
         {
            assert(nsummands < summandssize);
            SCIP_CALL( SCIPexprCreate(tree->blkmem, &summands[nsummands], SCIP_EXPR_CONST, scalar * SCIPvarGetMultaggrConstant(var) + constant) );
            ++nsummands;
         }

         if( nsummands == 0 )
         {
            /* somehow everything collapsed to zero */
            SCIP_CALL( SCIPexprCreate(tree->blkmem, &replaceexprs[i], SCIP_EXPR_CONST, 0.0) );
            BMSfreeBlockMemoryArray(tree->blkmem, &summands, summandssize);
            continue;
         }
         if( nsummands == 1 )
         {
            /* somehow everything collapsed to one summand -> use that one for replaceexprs[i]*/
            replaceexprs[i] = summands[0];
            BMSfreeBlockMemoryArray(tree->blkmem, &summands, summandssize);
            continue;
         }

         /* set replaceexprs[i] to summation over summands array */
         SCIP_CALL( SCIPexprCreate(tree->blkmem, &replaceexprs[i], SCIP_EXPR_SUM, nsummands, summands) );

         BMSfreeBlockMemoryArray(tree->blkmem, &summands, summandssize);
      }
   }

   /* replace variables in tree by assembled expressions */
   SCIP_CALL( SCIPexprtreeSubstituteVars(tree, replaceexprs) );
   /* free replaceexprs */
   for( i = 0; i < nvarsold; ++i )
      if( replaceexprs[i] != NULL )
         SCIPexprFreeDeep(tree->blkmem, &replaceexprs[i]);
   BMSfreeBlockMemoryArray(tree->blkmem, &replaceexprs, nvarsold);

   /* the varhash is not needed anymore */
   SCIPhashmapFree(&varhash);

   /* remove inactive variables from vars array and recompute variable indices */
   SCIP_ALLOC( BMSallocBlockMemoryArray(tree->blkmem, &newpos, tree->nvars) );
   offset = 0;
   for( i = 0; i < tree->nvars; ++i )
   {
      if( SCIPvarIsActive(tree->vars[i]) || i >= nvarsold )
      {
         /* a new variable need to be either active or multiaggregated */
         assert(i < nvarsold || SCIPvarIsActive(tree->vars[i]) || SCIPvarGetStatus(tree->vars[i]) == SCIP_VARSTATUS_MULTAGGR);
         newpos[i] = i - offset;
      }
      else
      {
         /* non-active variable are removed */
         newpos[i] = -1;
         ++offset;
      }
   }

   /* update indices in tree */
   SCIPexprReindexVars(tree->root, newpos);

   /* move variable in expression tree vars array
    * check if there is a fixed variable left */
   havefixedvar = FALSE;
   for( i = 0; i < tree->nvars; ++i )
   {
      if( newpos[i] == -1 )
      {
         /* variable was removed */
         assert(!SCIPvarIsActive(tree->vars[i]));
         continue;
      }
      /* variable is moved */
      tree->vars[newpos[i]] = tree->vars[i];
      if( !SCIPvarIsActive(tree->vars[i]) )
         havefixedvar = TRUE;
   }

   /* free newpos array; resize vars array */
   BMSfreeBlockMemoryArray(tree->blkmem, &newpos, tree->nvars);
   if( offset < tree->nvars )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(tree->blkmem, &tree->vars, tree->nvars, tree->nvars - offset) );
      tree->nvars -= offset;
   }
   else
   {
      /* all variables were removed */
      BMSfreeBlockMemoryArray(tree->blkmem, &tree->vars, tree->nvars);
      tree->nvars = 0;
   }

   if( havefixedvar )
   {
      SCIP_Bool dummy;
      /* call recursively if still unfixed variables are left */
      SCIP_CALL( SCIPexprtreeRemoveFixedVars(tree, &dummy) );
      assert(dummy == TRUE);
   }

   return SCIP_OKAY;
}
