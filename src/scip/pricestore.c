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
#pragma ident "@(#) $Id: pricestore.c,v 1.15 2004/07/09 08:11:34 bzfpfend Exp $"

/**@file   pricestore.c
 * @brief  methods for storing priced variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "clock.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
#include "tree.h"
#include "pricestore.h"

#include "struct_pricestore.h"



/*
 * dynamic memory arrays
 */

/** resizes vars and score arrays to be able to store at least num entries */
static
RETCODE pricestoreEnsureVarsMem(
   PRICESTORE*      pricestore,         /**< pricing storage */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(pricestore != NULL);
   assert(set != NULL);

   if( num > pricestore->varssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&pricestore->vars, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&pricestore->scores, newsize) );
      pricestore->varssize = newsize;
   }
   assert(num <= pricestore->varssize);

   return SCIP_OKAY;
}

/** resizes bdviolvars arrays to be able to store at least num entries */
static
RETCODE pricestoreEnsureBdviolvarsMem(
   PRICESTORE*      pricestore,         /**< pricing storage */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(pricestore != NULL);
   assert(set != NULL);

   if( num > pricestore->bdviolvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&pricestore->bdviolvars, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&pricestore->bdviolvarslb, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&pricestore->bdviolvarsub, newsize) );
      pricestore->bdviolvarssize = newsize;
   }
   assert(num <= pricestore->bdviolvarssize);

   return SCIP_OKAY;
}


/** creates pricing storage */
RETCODE SCIPpricestoreCreate(
   PRICESTORE**     pricestore          /**< pointer to store pricing storage */
   )
{
   assert(pricestore != NULL);
   
   ALLOC_OKAY( allocMemory(pricestore) );
   
   CHECK_OKAY( SCIPclockCreate(&(*pricestore)->probpricingtime, SCIP_CLOCKTYPE_DEFAULT) );
   (*pricestore)->vars = NULL;
   (*pricestore)->scores = NULL;
   (*pricestore)->bdviolvars = NULL;
   (*pricestore)->bdviolvarslb = NULL;
   (*pricestore)->bdviolvarsub = NULL;
   (*pricestore)->varssize = 0;
   (*pricestore)->nvars = 0;
   (*pricestore)->bdviolvarssize = 0;
   (*pricestore)->nbdviolvars = 0;
   (*pricestore)->naddedbdviolvars = 0;
   (*pricestore)->nprobpricings = 0;
   (*pricestore)->nprobvarsfound = 0;
   (*pricestore)->nvarsfound = 0;
   (*pricestore)->nvarsapplied = 0;
   (*pricestore)->initiallp = FALSE;

   return SCIP_OKAY;
}

/** frees pricing storage */
RETCODE SCIPpricestoreFree(
   PRICESTORE**     pricestore          /**< pointer to store pricing storage */
   )
{
   assert(pricestore != NULL);

   SCIPclockFree(&(*pricestore)->probpricingtime);
   freeMemoryArrayNull(&(*pricestore)->vars);
   freeMemoryArrayNull(&(*pricestore)->scores);
   freeMemoryArrayNull(&(*pricestore)->bdviolvars);
   freeMemoryArrayNull(&(*pricestore)->bdviolvarslb);
   freeMemoryArrayNull(&(*pricestore)->bdviolvarsub);
   freeMemory(pricestore);

   return SCIP_OKAY;
}

/** informs pricing storage, that the setup of the initial LP starts now */
void SCIPpricestoreStartInitialLP(
   PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);
   assert(!pricestore->initiallp);
   assert(pricestore->nvars == 0);

   pricestore->initiallp = TRUE;
}

/** informs pricing storage, that the setup of the initial LP is now finished */
void SCIPpricestoreEndInitialLP(
   PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);
   assert(pricestore->initiallp);
   assert(pricestore->nvars == 0);

   pricestore->initiallp = FALSE;
}

/** adds variable to pricing storage and capture it */
RETCODE SCIPpricestoreAddVar(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   VAR*             var,                /**< priced variable */
   Real             score,              /**< pricing score of variable (the larger, the better the variable) */
   Bool             root                /**< are we at the root node? */
   )
{
   int maxpricevars;
   int v;

   assert(pricestore != NULL);
   assert(set != NULL);
   assert(var != NULL);

   debugMessage("adding variable <%s> (lb=%g, ub=%g) to pricing storage\n", 
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

   if( pricestore->initiallp )
      maxpricevars = INT_MAX;
   else
   {
      pricestore->nvarsfound++;
      maxpricevars = root ? set->maxpricevarsroot : set->maxpricevars;
   }
   assert(maxpricevars >= 1);
   assert(pricestore->nvars <= maxpricevars);

   /* check, if variable belongs to the best "maxpricevars" pricing variables */
   if( pricestore->nvars < maxpricevars || score > pricestore->scores[maxpricevars-1] )
   {
      /* capture variable */
      SCIPvarCapture(var);

      /* if the array consists of "maxpricevars" variables, release the worst variables */
      if( pricestore->nvars == maxpricevars )
      {
         CHECK_OKAY( SCIPvarRelease(&pricestore->vars[pricestore->nvars-1], memhdr, set, lp) );
         pricestore->nvars--;
      }
      assert(pricestore->nvars < maxpricevars);

      /* get enough memory to store additional variable */
      CHECK_OKAY( pricestoreEnsureVarsMem(pricestore, set, pricestore->nvars+1) );
      assert(pricestore->nvars <= pricestore->varssize);

      /* insert the variable at the correct position in sorted arrays */
      for( v = pricestore->nvars; v > 0 && score > pricestore->scores[v-1]; --v )
      {
         pricestore->vars[v] = pricestore->vars[v-1];
         pricestore->scores[v] = pricestore->scores[v-1];
      }
      pricestore->vars[v] = var;
      pricestore->scores[v] = score;
      pricestore->nvars++;
   }

   return SCIP_OKAY;
}

/** adds variable where zero violates the bounds to pricing storage, capture it */
RETCODE SCIPpricestoreAddBdviolvar(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var                 /**< variable, where zero violates the bounds */
   )
{
   assert(pricestore != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPsetIsPositive(set, SCIPvarGetLbLocal(var)) || SCIPsetIsNegative(set, SCIPvarGetUbLocal(var)));
   assert(pricestore->naddedbdviolvars <= pricestore->nbdviolvars);

   debugMessage("zero violates bounds of <%s> (lb=%g, ub=%g)\n", 
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

   if( !pricestore->initiallp )
      pricestore->nvarsfound++;

   /* get enough memory to store additional variable */
   CHECK_OKAY( pricestoreEnsureBdviolvarsMem(pricestore, set, pricestore->nbdviolvars+1) );
   assert(pricestore->nbdviolvars <= pricestore->bdviolvarssize);

   /* capture variable */
   SCIPvarCapture(var);

   /* insert variable in bdviolvars arrays */
   pricestore->bdviolvars[pricestore->nbdviolvars] = var;
   pricestore->bdviolvarslb[pricestore->nbdviolvars] = SCIPvarGetLbLocal(var);
   pricestore->bdviolvarsub[pricestore->nbdviolvars] = SCIPvarGetUbLocal(var);
   pricestore->nbdviolvars++;

   /* Temporarily set bounds, such that zero is feasible, because we don't want to destroy
    * dual feasibility (by adding columns) and primal feasibility (by introducing violated bounds)
    * at the same time.
    * The correct bounds must be reset with a call to SCIPpricestoreResetBounds().
    * The inference information is unimportant for this temporary bound change.
    */
   if( SCIPsetIsPositive(set, SCIPvarGetLbLocal(var)) )
   {
      CHECK_OKAY( SCIPvarSetLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, 0.0) );
   }
   else
   {
      CHECK_OKAY( SCIPvarSetUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, 0.0) );
   }

   return SCIP_OKAY;
}

/** adds given problem variable to pricing storage, if zero is not best bound w.r.t. objective function */
static
RETCODE addBoundViolated(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var,                /**< problem variable */
   Bool*            added               /**< pointer to store whether variable was added to pricing storage */
   )
{
   assert(tree != NULL);
   assert(added != NULL);

   *added = FALSE;

   /* add variable, if zero is not feasible within the bounds */
   if( SCIPsetIsPositive(set, SCIPvarGetLbLocal(var)) || SCIPsetIsNegative(set, SCIPvarGetUbLocal(var)) )
   {
      debugMessage(" -> zero violates bounds of <%s> [%g,%g]\n",
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      CHECK_OKAY( SCIPpricestoreAddBdviolvar(pricestore, memhdr, set, stat, lp, branchcand, eventqueue, var) );
      *added = TRUE;
   }
   else
   {
      Real bestbound;

      /* add variable, if zero is not best bound w.r.t. objective function */
      bestbound = SCIPvarGetBestBound(var);
      if( !SCIPsetIsZero(set, bestbound) )
      {
         debugMessage(" -> best bound of <%s> [%g,%g] is not zero but %g\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), bestbound);
         CHECK_OKAY( SCIPpricestoreAddVar(pricestore, memhdr, set, lp, var, 
                        -SCIPvarGetObj(var) * bestbound, (SCIPnodeGetDepth(tree->actnode) == 0)) );
         *added = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** adds problem variables with negative reduced costs to pricing storage */
RETCODE SCIPpricestoreAddProbVars(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   VAR* var;
   COL* col;
   Bool root;
   Bool added;
   int v;
   int abortpricevars;
   int maxpricevars;
   int nfoundvars;

   assert(pricestore != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnodehaslp);
   assert(prob->nvars >= SCIPlpGetNCols(lp));

   /* if all problem variables of status COLUMN are already in the LP, nothing has to be done */
   if( prob->ncolvars == SCIPlpGetNCols(lp) )
      return SCIP_OKAY;

   root = (SCIPnodeGetDepth(tree->actnode) == 0);
   maxpricevars = SCIPsetGetMaxpricevars(set, root);
   assert(maxpricevars >= 1);
   abortpricevars = (int)(set->abortpricevarsfac * maxpricevars);
   assert(abortpricevars >= maxpricevars);
   
   /**@todo test pricing: is abortpricevars a good idea? -> like strong branching, lookahead, ... */

   pricestore->nprobpricings++;

   /* start timing */
   SCIPclockStart(pricestore->probpricingtime, set);
   
   /* price already existing problem variables */
   nfoundvars = 0;
   for( v = 0; v < prob->nvars && nfoundvars < abortpricevars; ++v )
   {
      var = prob->vars[v];
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      {
         col = SCIPvarGetCol(var);
         assert(col != NULL);
         assert(col->var == var);
         assert(col->len >= 0);
         assert(col->lppos >= -1);
         assert(col->lpipos >= -1);
         assert(SCIPcolIsInLP(col) == (col->lpipos >= 0));
            
         if( !SCIPcolIsInLP(col) )
         {
            debugMessage("price column variable <%s> in bounds [%g,%g]\n", 
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

            /** add variable to pricing storage, if zero is not best bound w.r.t. objective function */
            CHECK_OKAY( addBoundViolated(pricestore, memhdr, set, stat, tree, lp, branchcand, eventqueue, var, &added) );

            if( added )
            {
               pricestore->nprobvarsfound++;
               nfoundvars++;
            }
            else if( SCIPcolGetNNonz(col) > 0 )
            {
               Real feasibility;
   
               /* a column not in LP that doesn't have zero in its bounds was added by bound checking above */
               assert(!SCIPsetIsPositive(set, SCIPvarGetLbLocal(col->var)));
               assert(!SCIPsetIsNegative(set, SCIPvarGetUbLocal(col->var)));
               
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
                  feasibility = -SCIPcolGetFarkas(col, stat, lp);
                  debugMessage("  <%s> farkas feasibility: %e\n", SCIPvarGetName(col->var), feasibility);
               }
               else
               {
                  /* The dual LP is feasible, and we have a feasible dual solution. Pricing in this case means to
                   * add variables with negative feasibility, that is
                   *  - positive reduced costs for variables with negative lower bound
                   *  - negative reduced costs for variables with positive upper bound
                   */
                  feasibility = SCIPcolGetFeasibility(col, set, stat, lp);
                  debugMessage("  <%s> reduced cost feasibility: %e\n", SCIPvarGetName(col->var), feasibility);
               }
               
               /* the score is -feasibility / (#nonzeros in column + 1) to prefer short columns
                * we must add variables with negative feasibility, but in order to not get a too large lower bound
                * due to missing columns, we better also add variables, that have a very small feasibility
                */
               if( !SCIPsetIsPositive(set, feasibility) )
               {
                  CHECK_OKAY( SCIPpricestoreAddVar(pricestore, memhdr, set, lp, var, -feasibility / (col->len+1), root) );
                  pricestore->nprobvarsfound++;
                  nfoundvars++;
               }
            }
         }
      }
   }

   /* stop timing */
   SCIPclockStop(pricestore->probpricingtime, set);

   return SCIP_OKAY;
}

/** adds priced variables to the LP */
RETCODE SCIPpricestoreApplyVars(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< LP data */
   )
{
   VAR* var;
   COL* col;
   int v;

   assert(pricestore != NULL);
   assert(pricestore->naddedbdviolvars <= pricestore->nbdviolvars);
   assert(set != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->solved || pricestore->initiallp);
   assert(tree != NULL);
   assert(tree->actnode != NULL);
   assert(tree->actnodehaslp);

   debugMessage("adding %d variables (%d bound violated and %d priced vars) to %d LP columns\n",
      SCIPpricestoreGetNVars(pricestore), pricestore->nbdviolvars - pricestore->naddedbdviolvars,
      pricestore->nvars, SCIPlpGetNCols(lp));

   /* add the variables with violated bounds to LP */
   for( v = pricestore->naddedbdviolvars; v < pricestore->nbdviolvars; ++v )
   {
      var = pricestore->bdviolvars[v];
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetProbindex(var) >= 0);
      assert(var->nuses >= 2); /* at least used in pricing storage and in problem */

      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         /* transform loose variable into column variable */
         CHECK_OKAY( SCIPvarColumn(var, memhdr, set, stat, prob, lp) );
      }
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

      col = SCIPvarGetCol(var);
      assert(col != NULL);
      assert(col->lppos == -1);
      assert(col->lpipos == -1);
      debugMessage("adding bound violated variable <%s> (lb=%g, ub=%g)\n", SCIPvarGetName(var), 
         pricestore->bdviolvarslb[v], pricestore->bdviolvarsub[v]);
      CHECK_OKAY( SCIPlpAddCol(lp, set, col) );

      if( !pricestore->initiallp )
         pricestore->nvarsapplied++;
   }
   pricestore->naddedbdviolvars = pricestore->nbdviolvars;

   /* add the selected pricing variables to LP */
   for( v = 0; v < pricestore->nvars; ++v )
   {
      var = pricestore->vars[v];
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetProbindex(var) >= 0);
      assert(var->nuses >= 2); /* at least used in pricing storage and in problem */

      /* transform variable into column variable, if needed */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         CHECK_OKAY( SCIPvarColumn(var, memhdr, set, stat, prob, lp) );
      }
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

      col = SCIPvarGetCol(var);
      assert(col != NULL);
      assert(col->lppos == -1);
      assert(col->lpipos == -1);
      debugMessage("adding priced variable <%s> (score=%g)\n", SCIPvarGetName(var), pricestore->scores[v]);
      CHECK_OKAY( SCIPlpAddCol(lp, set, col) );

      /* release the variable */
      CHECK_OKAY( SCIPvarRelease(&pricestore->vars[v], memhdr, set, lp) );

      if( !pricestore->initiallp )
         pricestore->nvarsapplied++;
   }

   /* clear the pricing storage */
   pricestore->nvars = 0;

   return SCIP_OKAY;
}

/** reset variables' bounds violated by zero to its original value */
RETCODE SCIPpricestoreResetBounds(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   VAR* var;
   int v;

   assert(pricestore != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(lp->solved);
   assert(pricestore->nvars == 0);
   assert(pricestore->naddedbdviolvars == pricestore->nbdviolvars);

   /* reset variables' bounds, release them, and clear the boundviolation storage;
    * the inference information is unimportant in these removals of temporary bound changes
    */
   for( v = 0; v < pricestore->nbdviolvars; ++v )
   {
      var = pricestore->bdviolvars[v];
      assert(var != NULL);

      debugMessage("resetting bounds of <%s> from [%g,%g] to [%g,%g]\n", var->name, 
         SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), pricestore->bdviolvarslb[v], pricestore->bdviolvarsub[v]);
      CHECK_OKAY( SCIPvarSetLbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, pricestore->bdviolvarslb[v]) );
      CHECK_OKAY( SCIPvarSetUbLocal(var, memhdr, set, stat, lp, branchcand, eventqueue, pricestore->bdviolvarsub[v]) );
      CHECK_OKAY( SCIPvarRelease(&pricestore->bdviolvars[v], memhdr, set, lp) );
   }
   pricestore->naddedbdviolvars = 0;
   pricestore->nbdviolvars = 0;

   return SCIP_OKAY;
}

/** gets number of variables in pricing storage */
int SCIPpricestoreGetNVars(
   PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);
   assert(pricestore->nbdviolvars >= pricestore->naddedbdviolvars);

   return pricestore->nvars + pricestore->nbdviolvars - pricestore->naddedbdviolvars;
}

/** gets number of variables in pricing storage whose bounds must be reset */
int SCIPpricestoreGetNBoundResets(
   PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);
   assert(pricestore->nbdviolvars >= pricestore->naddedbdviolvars);

   return pricestore->nbdviolvars - pricestore->naddedbdviolvars;
}

/** gets time needed to price existing problem variables */
Real SCIPpricestoreGetProbPricingTime(
   PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);

   return SCIPclockGetTime(pricestore->probpricingtime);
}

/** gets total number of calls to problem variable pricing */
int SCIPpricestoreGetNProbPricings(
   PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->nprobpricings;
}

/** gets total number of times, a problem variable was priced in */
int SCIPpricestoreGetNProbvarsFound(
   PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->nprobvarsfound;
}

/** get total number of variables found so far in pricing */
int SCIPpricestoreGetNVarsFound(
   PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->nvarsfound;
}

/** get total number of variables priced into the LP so far */
int SCIPpricestoreGetNVarsApplied(
   PRICESTORE*      pricestore          /**< pricing storage */
   )
{
   assert(pricestore != NULL);

   return pricestore->nvarsapplied;
}

