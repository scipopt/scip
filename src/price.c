/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   price.c
 * @brief  methods and datastructures for pricing variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "set.h"
#include "mem.h"
#include "prob.h"
#include "lp.h"
#include "stat.h"
#include "var.h"
#include "scip.h"
#include "price.h"


/** storage for priced variables */
struct Price
{
   VAR**            vars;               /**< array with priced variables with neg. reduced costs sorted by score */
   Real*            score;              /**< score for each priced variable (e.g. -redcost/#nonzeros) */
   VAR**            bdviolvars;         /**< variables where zero violates the bounds */
   Real*            bdviolvarslb;       /**< lower bounds of bdviolvars */
   Real*            bdviolvarsub;       /**< upper bounds of bdbiolvars */
   int              varssize;           /**< size of vars and score arrays */
   int              nvars;              /**< number of priced variables (max. is set->maxpricevars) */
   int              bdviolvarssize;     /**< size of bdviolvars, bdviolvarslb, and bdviolvarsub arrays */
   int              nbdviolvars;        /**< number of variables, where zero violates the bounds */
};


/*
 * dymanic memory arrays
 */

static
RETCODE priceEnsureVarsMem(             /**< resizes vars and score arrays to be able to store at least num entries */
   PRICE*           price,              /**< pricing storage */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(price != NULL);
   assert(set != NULL);

   if( num > price->varssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(price->vars, newsize) );
      ALLOC_OKAY( reallocMemoryArray(price->score, newsize) );
      price->varssize = newsize;
   }
   assert(num <= price->varssize);

   return SCIP_OKAY;
}

static
RETCODE priceEnsureBdviolvarsMem(       /**< resizes bdviolvars arrays to be able to store at least num entries */
   PRICE*           price,              /**< pricing storage */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(price != NULL);
   assert(set != NULL);

   if( num > price->bdviolvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(price->bdviolvars, newsize) );
      ALLOC_OKAY( reallocMemoryArray(price->bdviolvarslb, newsize) );
      ALLOC_OKAY( reallocMemoryArray(price->bdviolvarsub, newsize) );
      price->bdviolvarssize = newsize;
   }
   assert(num <= price->bdviolvarssize);

   return SCIP_OKAY;
}


RETCODE SCIPpriceCreate(                /**< creates pricing storage */
   PRICE**          price               /**< pointer to store pricing storage */
   )
{
   assert(price != NULL);
   
   ALLOC_OKAY( allocMemory(*price) );
   
   (*price)->vars = NULL;
   (*price)->score = NULL;
   (*price)->bdviolvars = NULL;
   (*price)->bdviolvarslb = NULL;
   (*price)->bdviolvarsub = NULL;
   (*price)->varssize = 0;
   (*price)->nvars = 0;
   (*price)->bdviolvarssize = 0;
   (*price)->nbdviolvars = 0;

   return SCIP_OKAY;
}

RETCODE SCIPpriceFree(                  /**< frees pricing storage */
   PRICE**          price               /**< pointer to store pricing storage */
   )
{
   assert(price != NULL);

   freeMemoryArrayNull((*price)->vars);
   freeMemoryArrayNull((*price)->score);
   freeMemoryArrayNull((*price)->bdviolvars);
   freeMemoryArrayNull((*price)->bdviolvarslb);
   freeMemoryArrayNull((*price)->bdviolvarsub);
   freeMemory(*price);

   return SCIP_OKAY;
}

RETCODE SCIPpriceAddVar(                /**< adds variable to pricing storage */
   PRICE*           price,              /**< pricing storage */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< priced variable */
   Real             score               /**< pricing score of variable (the larger, the better the variable) */
   )
{
   int v;
   int i;

   assert(price != NULL);
   assert(set != NULL);
   assert(var != NULL);

   /* get enough memory to store "maxpricevars" variables */
   CHECK_OKAY( priceEnsureVarsMem(price, set, set->maxpricevars) );
   assert(price->nvars <= price->varssize);

   /* check, if variable belongs to the best "maxpricevars" pricing variables */
   if( price->nvars < set->maxpricevars || score > price->score[set->maxpricevars-1] )
   {
      /* insert variable in the sorted arrays */
      for( v = 0; v < price->nvars && score <= price->score[v]; ++v )
      {
      }
      assert(v < set->maxpricevars);
      price->nvars = MIN(price->nvars+1, set->maxpricevars);
      for( i = price->nvars-1; i > v; --i )
      {
         price->vars[i] = price->vars[i-1];
         price->score[i] = price->score[i-1];
      }
      price->vars[v] = var;
      price->score[v] = score;
   }

   return SCIP_OKAY;
}

RETCODE SCIPpriceAddBdviolvar(          /**< adds variable where zero violates the bounds to array of bound violations */
   PRICE*           price,              /**< pricing storage */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   VAR*             var                 /**< variable, where zero violates the bounds */
   )
{
   assert(price != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPsetIsPos(set, var->dom.lb) || SCIPsetIsNeg(set, var->dom.ub));

   debugMessage("zero violates bounds of <%s> (lb=%g, ub=%g)\n", var->name, var->dom.lb, var->dom.ub);

   /* get enough memory to store "maxpricevars" variables */
   CHECK_OKAY( priceEnsureBdviolvarsMem(price, set, price->nbdviolvars+1) );
   assert(price->nbdviolvars <= price->bdviolvarssize);

   /* insert variable in bdviolvars arrays */
   price->bdviolvars[price->nbdviolvars] = var;
   price->bdviolvarslb[price->nbdviolvars] = var->dom.lb;
   price->bdviolvarsub[price->nbdviolvars] = var->dom.ub;
   price->nbdviolvars++;

   /* Temporarily set bounds, such that zero is feasible, because we don't want to destroy
    * dual feasibility (by adding columns) and primal feasibility (by introducing violated bounds)
    * at the same time.
    * The correct bounds must be reset with a call to SCIPpriceResetBounds().
    */
   if( SCIPsetIsPos(set, var->dom.lb) )
   {
      CHECK_OKAY( SCIPvarChgLb(var, set, lp, 0.0) );
   }
   else
   {
      CHECK_OKAY( SCIPvarChgUb(var, set, lp, 0.0) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPpriceVars(                  /**< calls all external pricer, prices problem variables, and adds some columns
                                           with negative reduced costs to the LP */
   PRICE*           price,              /**< pricing storage */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            transprob,          /**< transformed problem after presolve */
   LP*              lp                  /**< LP data */
   )
{
   VAR* var;
   COL* col;
   ROW* row;
   int v;
   int r;

   assert(price != NULL);
   assert(set != NULL);
   assert(transprob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(price->nvars == 0);

   /* call external pricer algorithms */
   todoMessage("external pricing");

   /* price problem variables */
   for( v = 0; v < transprob->nvars; ++v )
   {
      var = transprob->vars[v];
      switch( var->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         errorMessage("Found original variable in transformed problem");
         abort();
      case SCIP_VARSTATUS_LOOSE:
         /* A loose variable is a pricing candidate, if it can contribute negatively to the objective function.
          * In addition, we have to add all variables, where zero violates the bounds.
          */
         if( SCIPsetIsNeg(set, var->dom.lb) )
         {
            if( SCIPsetIsNeg(set, var->dom.ub) )
            {
               CHECK_OKAY( SCIPpriceAddBdviolvar(price, set, lp, var) );
            }
            else if( SCIPsetIsPos(set, var->obj) )
            {
               CHECK_OKAY( SCIPpriceAddVar(price, set, var, -SCIP_PRICE_SCALELOOSE * var->obj * var->dom.lb) );
            }
         }
         else if( SCIPsetIsPos(set, var->dom.ub) )
         {
            if( SCIPsetIsPos(set, var->dom.lb) )
            {
               CHECK_OKAY( SCIPpriceAddBdviolvar(price, set, lp, var) );
            }
            else if( SCIPsetIsNeg(set, var->obj) )
            {
               CHECK_OKAY( SCIPpriceAddVar(price, set, var, -SCIP_PRICE_SCALELOOSE * var->obj * var->dom.ub) );
            }
         }
         break;
      case SCIP_VARSTATUS_COLUMN:
         col = var->data.col;
         assert(col != NULL);
         assert(col->var == var);
         assert(col->lpipos >= -1);
         assert(col->inLP || col->lpipos == -1);
         assert(!col->inLP || col->lpipos >= 0);
         assert(col->len >= 0);

         if( !col->inLP )
         {
            /* column is not in LP -> calculate reduced costs */
            SCIPcolCalcRedcost(col);
            assert(col->redcost < SCIP_INVALID);

            /* price in variable, if reduced costs are negative, 
             * or if variable can be negative and reduced costs are nonzero;
             * the score is |redcost| / (#nonzeros in column + 1) to prefer short columns
             */
            if( SCIPsetIsNeg(set, col->redcost)
               || (SCIPsetIsNeg(set, var->dom.lb) && !SCIPsetIsZero(set, col->redcost)) )
            {
               CHECK_OKAY( SCIPpriceAddVar(price, set, var, ABS(col->redcost) / (col->len+1)) );
            }
         }
         break;
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_AGGREGATED:
         /* we don't have to price fixed or aggregated variables */
         break;
      }
   }

   /* add the selected pricing variables to LP */
   for( v = 0; v < price->nvars; ++v )
   {
      var = price->vars[v];
      if( var->varstatus == SCIP_VARSTATUS_LOOSE )
      {
         /* transform loose variable into column variable */
         CHECK_OKAY( SCIPvarColumn(var, memhdr, set, lp, stat) );
      }
      assert(var->varstatus == SCIP_VARSTATUS_COLUMN);

      col = var->data.col;
      assert(col != NULL);
      assert(!col->inLP);
      assert(col->lpipos == -1);
      debugMessage("pricing in variable <%s> (score=%g)\n", var->name, price->score[v]);
      CHECK_OKAY( SCIPlpAddCol(lp, set, col) );
   }
   price->nvars = 0;

   return SCIP_OKAY;
}

RETCODE SCIPpriceResetBounds(           /**< reset variables' bounds violated by zero to its original value */
   PRICE*           price,              /**< pricing storage */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< LP data */
   )
{
   int v;

   assert(price != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);

   for( v = 0; v < price->nbdviolvars; ++v )
   {
      CHECK_OKAY( SCIPvarChgLb(price->bdviolvars[v], set, lp, price->bdviolvarslb[v]) );
      CHECK_OKAY( SCIPvarChgUb(price->bdviolvars[v], set, lp, price->bdviolvarsub[v]) );
   }
   price->nbdviolvars = 0;

   return SCIP_OKAY;
}
