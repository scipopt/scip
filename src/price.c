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
   VAR**            vars;               /**< array with priced variables with violated reduced costs sorted by score */
   Real*            score;              /**< score for each priced variable (e.g. |redcost|/#nonzeros) */
   VAR**            bdviolvars;         /**< variables where zero violates the bounds */
   Real*            bdviolvarslb;       /**< lower bounds of bdviolvars */
   Real*            bdviolvarsub;       /**< upper bounds of bdbiolvars */
   int              varssize;           /**< size of vars and score arrays */
   int              nvars;              /**< number of priced variables (max. is set->maxpricevars) */
   int              bdviolvarssize;     /**< size of bdviolvars, bdviolvarslb, and bdviolvarsub arrays */
   int              nbdviolvars;        /**< number of variables, where zero violates the bounds */
   int              naddedbdviolvars;   /**< number of bound violated variables already added to the LP */
   int              nfoundvars;         /**< total number of variables, that were added (and evtl. thrown away) */
};


/*
 * dynamic memory arrays
 */

/** resizes vars and score arrays to be able to store at least num entries */
static
RETCODE priceEnsureVarsMem(
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
      ALLOC_OKAY( reallocMemoryArray(&price->vars, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&price->score, newsize) );
      price->varssize = newsize;
   }
   assert(num <= price->varssize);

   return SCIP_OKAY;
}

/** resizes bdviolvars arrays to be able to store at least num entries */
static
RETCODE priceEnsureBdviolvarsMem(
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
      ALLOC_OKAY( reallocMemoryArray(&price->bdviolvars, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&price->bdviolvarslb, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&price->bdviolvarsub, newsize) );
      price->bdviolvarssize = newsize;
   }
   assert(num <= price->bdviolvarssize);

   return SCIP_OKAY;
}


/** creates pricing storage */
RETCODE SCIPpriceCreate(
   PRICE**          price               /**< pointer to store pricing storage */
   )
{
   assert(price != NULL);
   
   ALLOC_OKAY( allocMemory(price) );
   
   (*price)->vars = NULL;
   (*price)->score = NULL;
   (*price)->bdviolvars = NULL;
   (*price)->bdviolvarslb = NULL;
   (*price)->bdviolvarsub = NULL;
   (*price)->varssize = 0;
   (*price)->nvars = 0;
   (*price)->bdviolvarssize = 0;
   (*price)->nbdviolvars = 0;
   (*price)->naddedbdviolvars = 0;
   (*price)->nfoundvars = 0;

   return SCIP_OKAY;
}

/** frees pricing storage */
RETCODE SCIPpriceFree(
   PRICE**          price               /**< pointer to store pricing storage */
   )
{
   assert(price != NULL);

   freeMemoryArrayNull(&(*price)->vars);
   freeMemoryArrayNull(&(*price)->score);
   freeMemoryArrayNull(&(*price)->bdviolvars);
   freeMemoryArrayNull(&(*price)->bdviolvarslb);
   freeMemoryArrayNull(&(*price)->bdviolvarsub);
   freeMemory(price);

   return SCIP_OKAY;
}

/** adds variable to pricing storage and capture it */
RETCODE SCIPpriceAddVar(
   PRICE*           price,              /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   VAR*             var,                /**< priced variable */
   Real             score,              /**< pricing score of variable (the larger, the better the variable) */
   Bool             root                /**< are we at the root node? */
   )
{
   int maxpricevars;
   int v;
   int i;

   assert(price != NULL);
   assert(set != NULL);
   assert(var != NULL);

   price->nfoundvars++;

   /* get enough memory to store variables */
   CHECK_OKAY( priceEnsureVarsMem(price, set, price->nvars+1) );
   assert(price->nvars <= price->varssize);

   maxpricevars = root ? set->maxpricevarsroot : set->maxpricevars;
   assert(maxpricevars >= 1);
   assert(price->nvars <= maxpricevars);

   /* check, if variable belongs to the best "maxpricevars" pricing variables */
   if( price->nvars < maxpricevars || score > price->score[maxpricevars-1] )
   {
      /* capture variable */
      SCIPvarCapture(var);

      /* search the variable's position in the sorted arrays */
      for( v = 0; v < price->nvars && score <= price->score[v]; ++v )
      {
      }
      assert(v <= price->nvars);
      assert(v < maxpricevars);

      /* if the array consists of "maxpricevars" variables, release the worst variables */
      if( price->nvars == maxpricevars )
      {
         CHECK_OKAY( SCIPvarRelease(&price->vars[price->nvars-1], memhdr, set, lp) );
         price->nvars--;
      }
      assert(price->nvars < maxpricevars);

      /* insert variable in the sorted arrays */
      for( i = price->nvars; i > v; --i )
      {
         price->vars[i] = price->vars[i-1];
         price->score[i] = price->score[i-1];
      }
      price->vars[v] = var;
      price->score[v] = score;
      price->nvars++;
   }

   return SCIP_OKAY;
}

/** adds variable where zero violates the bounds to pricing storage, capture it */
RETCODE SCIPpriceAddBdviolvar(
   PRICE*           price,              /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var                 /**< variable, where zero violates the bounds */
   )
{
   assert(price != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPsetIsPositive(set, var->actdom.lb) || SCIPsetIsNegative(set, var->actdom.ub));
   assert(price->naddedbdviolvars <= price->nbdviolvars);

   debugMessage("zero violates bounds of <%s> (lb=%g, ub=%g)\n", var->name, var->actdom.lb, var->actdom.ub);

   price->nfoundvars++;

   /* get enough memory to store "maxpricevars" variables */
   CHECK_OKAY( priceEnsureBdviolvarsMem(price, set, price->nbdviolvars+1) );
   assert(price->nbdviolvars <= price->bdviolvarssize);

   /* capture variable */
   SCIPvarCapture(var);

   /* insert variable in bdviolvars arrays */
   price->bdviolvars[price->nbdviolvars] = var;
   price->bdviolvarslb[price->nbdviolvars] = var->actdom.lb;
   price->bdviolvarsub[price->nbdviolvars] = var->actdom.ub;
   price->nbdviolvars++;

   /* Temporarily set bounds, such that zero is feasible, because we don't want to destroy
    * dual feasibility (by adding columns) and primal feasibility (by introducing violated bounds)
    * at the same time.
    * The correct bounds must be reset with a call to SCIPpriceResetBounds().
    */
   if( SCIPsetIsPositive(set, var->actdom.lb) )
   {
      CHECK_OKAY( SCIPvarChgLbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, 0.0) );
   }
   else
   {
      CHECK_OKAY( SCIPvarChgUbLocal(var, memhdr, set, stat, tree, lp, branchcand, eventqueue, 0.0) );
   }

   return SCIP_OKAY;
}

/** prices problem variables */
static
RETCODE priceProbVars(
   PRICE*           price,              /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   VAR* var;
   COL* col;
   Bool root;
   int v;
   int abortpricevars;
   int maxpricevars;
   
   assert(price != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnodehaslp);
   assert(prob->nvars > lp->ncols);

   root = (tree->actnode->depth == 0);
   maxpricevars = root ? set->maxpricevarsroot : set->maxpricevars;
   assert(maxpricevars >= 1);
   abortpricevars = (int)(set->abortpricevarsfac * maxpricevars);
   assert(abortpricevars >= maxpricevars);
   
   todoMessage("test pricing: is abortpricevars a good idea? -> like strong branching, lookahead, ...");

   stat->nlppricings++;

   /* start timing */
   SCIPclockStart(stat->lppricingtime, set);
   
   /* price already existing problem variables */
   for( v = 0; v < prob->nvars && price->nfoundvars < abortpricevars; ++v )
   {
      var = prob->vars[v];
      switch( var->varstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         errorMessage("Found original variable in transformed problem");
         abort();
      case SCIP_VARSTATUS_LOOSE:
         /* A loose variable is a pricing candidate, if it can contribute negatively to the objective function.
          * In addition, we have to add all variables, where zero violates the bounds.
          */
         /*debugMessage("price loose variable <%s> in bounds [%g,%g]\n", var->name, var->actdom.lb, var->actdom.ub);*/
         if( SCIPsetIsNegative(set, var->actdom.lb) )
         {
            if( SCIPsetIsNegative(set, var->actdom.ub) )
            {
               CHECK_OKAY( SCIPpriceAddBdviolvar(price, memhdr, set, stat, tree, lp, branchcand, eventqueue, var) );
               stat->nlppricingvars++;
            }
            else if( SCIPsetIsPositive(set, var->obj) )
            {
               CHECK_OKAY( SCIPpriceAddVar(price, memhdr, set, lp, var, -var->obj * var->actdom.lb, root) );
               stat->nlppricingvars++;
            }
         }
         else if( SCIPsetIsPositive(set, var->actdom.ub) )
         {
            if( SCIPsetIsPositive(set, var->actdom.lb) )
            {
               CHECK_OKAY( SCIPpriceAddBdviolvar(price, memhdr, set, stat, tree, lp, branchcand, eventqueue, var) );
               stat->nlppricingvars++;
            }
            else if( SCIPsetIsNegative(set, var->obj) )
            {
               CHECK_OKAY( SCIPpriceAddVar(price, memhdr, set, lp, var, -var->obj * var->actdom.ub, root) );
               stat->nlppricingvars++;
            }
         }
         break;
      case SCIP_VARSTATUS_COLUMN:
         col = var->data.col;
         assert(col != NULL);
         assert(col->var == var);
         assert(col->lppos >= -1);
         assert(col->lpipos >= -1);
         assert((col->lppos >= 0) ^ (col->lpipos == -1));
         assert(col->len >= 0);
            
         /*debugMessage("price column variable <%s> in bounds [%g,%g], col->lppos=%d\n", 
           var->name, var->actdom.lb, var->actdom.ub, col->lppos);*/
         if( col->lppos == -1 )
         {
            Real feasibility;
            Bool added;
   
            added = FALSE;

            /* add variable, if zero is not feasible within the bounds */
            if( SCIPsetIsPositive(set, var->actdom.lb) || SCIPsetIsNegative(set, var->actdom.ub) )
            {
               CHECK_OKAY( SCIPpriceAddBdviolvar(price, memhdr, set, stat, tree, lp, branchcand, eventqueue, var) );
               stat->nlppricingvars++;
               added = TRUE;
            }
            else
            {
               Real bestbound;

               /* if zero is not best bound w.r.t. objective function */
               bestbound = SCIPvarGetBestBound(var);
               if( !SCIPsetIsZero(set, bestbound) )
               {
                  CHECK_OKAY( SCIPpriceAddVar(price, memhdr, set, lp, var, var->obj * var->actdom.lb, root) );
                  stat->nlppricingvars++;
                  added = TRUE;
               }
            }
            
            if( !added )
            {
               /* a column not in LP that doesn't have zero in its bounds was added by bound checking above */
               assert(!SCIPsetIsPositive(set, col->var->actdom.lb));
               assert(!SCIPsetIsNegative(set, col->var->actdom.ub));
               
               if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE )
               {
                  /* The LP was proven infeasible, so we have an infeasibility proof by the dual farkas values y.
                   * The valid inequality  y^T A x >= y^T b  is violated by all x, especially by the (for this
                   * inequality most feasible solution) x' defined by 
                   *    x'_i = ub_i, if y^T A_i > 0
                   *    x'_i = 0   , if y^T A_i = 0
                   *    x'_i = lb_i, if y^T A_i < 0.
                   * Pricing in this case means to add variables i with positive farkas value, i.e. y^T A_i x'_i > 0
                   */
                  feasibility = -SCIPcolGetFarkas(col, stat);
                  debugMessage("  <%s> farkas feasibility: %e\n", col->var->name, feasibility);
               }
               else
               {
                  /* The dual LP is feasible, and we have a feasible dual solution. Pricing in this case means to
                   * add variables with negative feasibility, that is negative reduced costs for non-negative
                   * variables, and non-zero reduced costs for variables that can be negative.
                   */
                  feasibility = SCIPcolGetFeasibility(col, stat);
                  debugMessage("  <%s> reduced cost feasibility: %e\n", col->var->name, feasibility);
               }
               
               /* the score is -feasibility / (#nonzeros in column + 1) to prefer short columns
                * we must add variables with negative feasibility, but in order to not get a too large lower bound
                * due to missing columns, we better also add variables, that have a very small feasibility
                */
               if( !SCIPsetIsPositive(set, feasibility) )
               {
                  CHECK_OKAY( SCIPpriceAddVar(price, memhdr, set, lp, var, -feasibility / (col->len+1), root) );
                  stat->nlppricingvars++;
               }
            }
         }
         break;
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_AGGREGATED:
      case SCIP_VARSTATUS_MULTAGGR:
      case SCIP_VARSTATUS_NEGATED:
         /* we don't have to price fixed or aggregated variables */
         break;
      }
   }

   /* stop timing */
   SCIPclockStop(stat->lppricingtime, set);

   return SCIP_OKAY;
}

/** calls all external pricer, prices problem variables, and adds some columns with negative reduced costs to the LP */
RETCODE SCIPpriceVars(
   PRICE*           price,              /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   VAR* var;
   COL* col;
   int v;
   
   assert(price != NULL);
   assert(price->nvars == 0);
   assert(price->nfoundvars == 0);
   assert(price->naddedbdviolvars == price->nbdviolvars);
   assert(set != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnodehaslp);

   /* if there are problem variables not in the LP, price them */
   if( prob->nvars > lp->ncols )
   {
      CHECK_OKAY( priceProbVars(price, memhdr, set, stat, prob, tree, lp, branchcand, eventqueue) );
   }

   /* call external pricer algorithms */
   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE )
   {
      todoMessage("external farkas pricing");
   }
   else
   {
      todoMessage("external reduced cost pricing");
   }

   debugMessage("adding %d variables (%d bound violated and %d priced vars) to %d LP columns\n",
      price->nvars + price->nbdviolvars - price->naddedbdviolvars, price->nbdviolvars - price->naddedbdviolvars,
      price->nvars, lp->ncols);

   /* add the variables with violated bounds to LP */
   for( v = price->naddedbdviolvars; v < price->nbdviolvars; ++v )
   {
      var = price->bdviolvars[v];
      assert(var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN);

      /* add variable to problem, if needed */
      if( var->probindex == -1 )
      {
         CHECK_OKAY( SCIPprobAddVar(prob, memhdr, set, tree, branchcand, var) );
      }
      assert(var->probindex >= 0);
      assert(var->nuses >= 2); /* at least used in pricing storage and in problem */

      if( var->varstatus == SCIP_VARSTATUS_LOOSE )
      {
         /* transform loose variable into column variable */
         CHECK_OKAY( SCIPvarColumn(var, memhdr, set, stat) );
      }
      assert(var->varstatus == SCIP_VARSTATUS_COLUMN);

      col = var->data.col;
      assert(col != NULL);
      assert(col->lppos == -1);
      assert(col->lpipos == -1);
      debugMessage("adding bound violated variable <%s> (lb=%g, ub=%g)\n", var->name, 
         price->bdviolvarslb[v], price->bdviolvarsub[v]);
      CHECK_OKAY( SCIPlpAddCol(lp, set, col) );
   }
   price->naddedbdviolvars = price->nbdviolvars;

   /* add the selected pricing variables to LP */
   for( v = 0; v < price->nvars; ++v )
   {
      var = price->vars[v];
      assert(var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN);

      /* add variable to problem, if needed */
      if( var->probindex == -1 )
      {
         CHECK_OKAY( SCIPprobAddVar(prob, memhdr, set, tree, branchcand, var) );
      }
      assert(var->probindex >= 0);
      assert(var->nuses >= 2); /* at least used in pricing storage and in problem */

      /* transform variable into column variable, if needed */
      if( var->varstatus == SCIP_VARSTATUS_LOOSE )
      {
         CHECK_OKAY( SCIPvarColumn(var, memhdr, set, stat) );
      }
      assert(var->varstatus == SCIP_VARSTATUS_COLUMN);

      col = var->data.col;
      assert(col != NULL);
      assert(col->lppos == -1);
      assert(col->lpipos == -1);
      debugMessage("adding priced variable <%s> (score=%g)\n", var->name, price->score[v]);
      CHECK_OKAY( SCIPlpAddCol(lp, set, col) );

      /* release the variable */
      CHECK_OKAY( SCIPvarRelease(&price->vars[v], memhdr, set, lp) );
   }

   /* clear the pricing storage */
   price->nvars = 0;
   price->nfoundvars = 0;

   return SCIP_OKAY;
}

/** reset variables' bounds violated by zero to its original value */
RETCODE SCIPpriceResetBounds(
   PRICE*           price,              /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   int v;

   assert(price != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(price->nvars == 0);
   assert(price->naddedbdviolvars == price->nbdviolvars);

   /* reset variables' bounds, release them, and clear the boundviolation storage */
   for( v = 0; v < price->nbdviolvars; ++v )
   {
      debugMessage("resetting bounds of <%s> from [%g,%g] to [%g,%g]\n", price->bdviolvars[v]->name, 
         price->bdviolvars[v]->actdom.lb, price->bdviolvars[v]->actdom.ub,
         price->bdviolvarslb[v], price->bdviolvarsub[v]);
      CHECK_OKAY( SCIPvarChgLbLocal(price->bdviolvars[v], memhdr, set, stat, tree, lp, branchcand, eventqueue,
                     price->bdviolvarslb[v]) );
      CHECK_OKAY( SCIPvarChgUbLocal(price->bdviolvars[v], memhdr, set, stat, tree, lp, branchcand, eventqueue,
                     price->bdviolvarsub[v]) );
      CHECK_OKAY( SCIPvarRelease(&price->bdviolvars[v], memhdr, set, lp) );
   }
   price->naddedbdviolvars = 0;
   price->nbdviolvars = 0;

   return SCIP_OKAY;
}
