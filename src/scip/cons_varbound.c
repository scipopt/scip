/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_varbound.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for variable bound constraints
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_varbound.h"
#include "scip/cons_linear.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "varbound"
#define CONSHDLR_DESC          "variable bounds  lhs <= x + c*y <= rhs, x non-binary, y non-continuous"
#define CONSHDLR_SEPAPRIORITY   +900000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -500000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -500000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "varbound"
#define EVENTHDLR_DESC         "bound change event handler for variable bound constraints"

#define LINCONSUPGD_PRIORITY     +50000 /**< priority of the constraint handler for upgrading of linear constraints */


/** variable bound constraint data */
struct SCIP_ConsData
{
   SCIP_Real             vbdcoef;            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs;                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs;                /**< right hand side of variable bound inequality */
   SCIP_VAR*             var;                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar;             /**< binary, integer or implicit integer bounding variable y */
   SCIP_ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   unsigned int          propagated:1;       /**< is the variable bound constraint already propagated? */
   unsigned int          presolved:1;        /**< is the variable bound constraint already presolved? */
   unsigned int          addvarbounds:1;     /**< are the globally valid variable bound are added? */
};



/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_1,                          /**< left hand side and bounds on y -> lower bound on x */
   PROPRULE_2,                          /**< left hand side and upper bound on x -> bound on y */
   PROPRULE_3,                          /**< right hand side and bounds on y -> upper bound on x */
   PROPRULE_4,                          /**< right hand side and lower bound on x -> bound on y */
   PROPRULE_INVALID                     /**< propagation was applied without a specific propagation rule */
};
typedef enum Proprule PROPRULE;




/*
 * Local methods
 */

/** catches events for variables */
static
SCIP_RETCODE catchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< variable bound constraint data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   
   assert(consdata != NULL);

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* catch bound change events on variables */
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->var, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, 
         (SCIP_EVENTDATA*)consdata, NULL) );
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vbdvar, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, 
         (SCIP_EVENTDATA*)consdata, NULL) );
   
   return SCIP_OKAY;
}

/** drops events for variables */
static
SCIP_RETCODE dropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< variable bound constraint data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   
   assert(consdata != NULL);

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* drop events on variables */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->var, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, 
         (SCIP_EVENTDATA*)consdata, -1) );
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vbdvar, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr,
         (SCIP_EVENTDATA*)consdata, -1) );

   return SCIP_OKAY;
}

/** creates a variable bound constraint data object */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the variable bound constraint data */
   SCIP_VAR*             var,                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   SCIP_Real             vbdcoef,            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs,                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs                 /**< right hand side of variable bound inequality */
   )
{
   assert(consdata != NULL);
   assert(SCIPvarGetType(vbdvar) != SCIP_VARTYPE_CONTINUOUS);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->var = var;
   (*consdata)->vbdvar = vbdvar;
   (*consdata)->vbdcoef = vbdcoef;
   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;
   (*consdata)->row = NULL;
   (*consdata)->propagated = FALSE;
   (*consdata)->presolved = FALSE;
   (*consdata)->addvarbounds = FALSE;

   /* capture variables */
   SCIP_CALL( SCIPcaptureVar(scip, var) );
   SCIP_CALL( SCIPcaptureVar(scip, vbdvar) );

   /* if we are in the transformed problem, get transformed variables, add variable bound information, and catch events */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->var, &(*consdata)->var) );
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->vbdvar, &(*consdata)->vbdvar) );

      /* catch events for variables */
      SCIP_CALL( catchEvents(scip, *consdata) );
   }

   return SCIP_OKAY;
}   

/** frees a variable bound constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to the variable bound constraint */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* drop events */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( dropEvents(scip, *consdata) );
   }

   /* release variables */
   SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->var) );
   SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vbdvar) );


   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** creates LP row corresponding to variable bound constraint */
static 
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< variable bound constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), consdata->lhs, consdata->rhs,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, consdata->var, 1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, consdata->vbdvar, consdata->vbdcoef) );

   return SCIP_OKAY;
}  

/** adds linear relaxation of variable bound constraint to the LP */
static 
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< variable bound constraint */
   )
{
   SCIP_CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->row != NULL);

   if( !SCIProwIsInLP(consdata->row) )
   {
      SCIPdebugMessage("adding relaxation of variable bound constraint <%s>: ", SCIPconsGetName(cons));
      SCIPdebug( SCIProwPrint(consdata->row, NULL) );
      SCIP_CALL( SCIPaddCut(scip, NULL, consdata->row, FALSE) );
   }

   return SCIP_OKAY;
}

/** separates the given variable bound constraint */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            separated           /**< pointer to store whether a cut was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real feasibility;

   assert(separated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("separating variable bound constraint <%s>\n", SCIPconsGetName(cons));

   *separated = FALSE;

   /* create LP relaxation if not yet existing */
   if( consdata->row == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->row != NULL);

   /* check non-LP rows for feasibility and add them as cut, if violated */
   if( sol != NULL || !SCIProwIsInLP(consdata->row) )
   {
      feasibility = SCIPgetRowSolFeasibility(scip, consdata->row, sol);
      if( SCIPisFeasNegative(scip, feasibility) )
      {
         SCIP_CALL( SCIPaddCut(scip, sol, consdata->row, FALSE) );
         *separated = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** returns whether the given solution is feasible for the given variable bound constraint */
static
SCIP_Bool checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool             checklprows         /**< should LP rows be checked? */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("checking variable bound constraint <%s> for feasibility of solution %p (lprows=%u)\n",
      SCIPconsGetName(cons), (void*)sol, checklprows);

   if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
   {
      SCIP_Real sum;

      sum = SCIPgetSolVal(scip, sol, consdata->var);
      sum += consdata->vbdcoef * SCIPgetSolVal(scip, sol, consdata->vbdvar);

      return SCIPisFeasGE(scip, sum, consdata->lhs) && SCIPisFeasLE(scip, sum, consdata->rhs);
   }
   else
      return TRUE;
}


/** resolves a propagation on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) left hand side and bounds on y -> lower bound on x
 *   (2) left hand side and upper bound on x -> bound on y
 *   (3) right hand side and bounds on y -> upper bound on x
 *   (4) right hand side and lower bound on x -> bound on y
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   PROPRULE              proprule,           /**< propagation rule that deduced the bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, consdata->vbdcoef));

   switch( proprule )
   {
   case PROPRULE_1:
      /* lhs <= x + c*y: left hand side and bounds on y -> lower bound on x */
      assert(infervar == consdata->var);
      assert(boundtype == SCIP_BOUNDTYPE_LOWER);
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      if( consdata->vbdcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->vbdvar, bdchgidx) );
      }
      else
      {
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->vbdvar, bdchgidx) );
      }
      break;

   case PROPRULE_2:
      /* lhs <= x + c*y: left hand side and upper bound on x -> bound on y */
      assert(infervar == consdata->vbdvar);
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      SCIP_CALL( SCIPaddConflictUb(scip, consdata->var, bdchgidx) );
      break;

   case PROPRULE_3:
      /* x + c*y <= rhs: right hand side and bounds on y -> upper bound on x */
      assert(infervar == consdata->var);
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      assert(!SCIPisInfinity(scip, consdata->rhs));
      if( consdata->vbdcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->vbdvar, bdchgidx) );
      }
      else
      {
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->vbdvar, bdchgidx) );
      }
      break;

   case PROPRULE_4:
      /* x + c*y <= rhs: right hand side and lower bound on x -> bound on y */
      assert(infervar == consdata->vbdvar);
      assert(!SCIPisInfinity(scip, consdata->rhs));
      SCIP_CALL( SCIPaddConflictLb(scip, consdata->var, bdchgidx) );
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in variable bound constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** analyze infeasibility */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   PROPRULE              proprule,           /**< propagation rule that deduced the bound change */
   SCIP_BOUNDTYPE        boundtype           /**< the type of the changed bound (lower or upper bound) */
   )
{
   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );

   /* add the bound which got violated */
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( SCIPaddConflictUb(scip, infervar, NULL) );
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      SCIP_CALL( SCIPaddConflictLb(scip, infervar, NULL) );
   }

   /* add the reason for the violated of the bound */
   SCIP_CALL( resolvePropagation(scip, cons, infervar, proprule, boundtype, NULL) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}


/** propagation method for variable bound constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds,            /**< pointer to count number of bound changes */
   int*                  ndelconss           /**< pointer to count number of deleted constraints, or NULL */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Real ylb;
   SCIP_Real yub;
   SCIP_Real newlb;
   SCIP_Real newub;
   SCIP_Bool tightened;
   SCIP_Bool tightenedround;

   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("propagating variable bound constraint <%s>: %.15g <= <%s> + %.15g<%s> <= %.15g\n",
      SCIPconsGetName(cons), consdata->lhs, SCIPvarGetName(consdata->var), consdata->vbdcoef,
      SCIPvarGetName(consdata->vbdvar), consdata->rhs);

   *cutoff = FALSE;

   /* check, if constraint is already propagated */
   if( consdata->propagated )
      return SCIP_OKAY;

   /* get current bounds of variables */
   xlb = SCIPvarGetLbLocal(consdata->var);
   xub = SCIPvarGetUbLocal(consdata->var);
   ylb = SCIPvarGetLbLocal(consdata->vbdvar);
   yub = SCIPvarGetUbLocal(consdata->vbdvar);

   /* tighten bounds of variables as long as possible */
   do
   {
      tightenedround = FALSE;

      /* propagate left hand side inequality: lhs <= x + c*y */
      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         assert(!(*cutoff));

         /* propagate bounds on x:
          *  (1) left hand side and bounds on y -> lower bound on x
          */
         if( SCIPvarGetStatus(consdata->var) != SCIP_VARSTATUS_MULTAGGR ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               if( !SCIPisInfinity(scip, yub) )
                  newlb = SCIPadjustedVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * yub);
               else
                  newlb = -SCIPinfinity(scip);
            }
            else
            {
               if( !SCIPisInfinity(scip, -ylb) )
                  newlb = SCIPadjustedVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * ylb);
               else
                  newlb = -SCIPinfinity(scip);
            }
            
            if( SCIPisLbBetter(scip, newlb, xlb, xub) || ylb > yub - 0.5 )
            {
               SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                  SCIPvarGetName(consdata->var), xlb, xub, newlb, xub);
               SCIP_CALL( SCIPinferVarLbCons(scip, consdata->var, newlb, cons, (int)PROPRULE_1, FALSE,
                     cutoff, &tightened) );
               
               if( *cutoff )
               {
                  assert(SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->var)));
                 
                  /* analyze infeasibility */
                  SCIP_CALL( analyzeConflict(scip, cons, consdata->var, PROPRULE_1, SCIP_BOUNDTYPE_LOWER) );
                  break;
               }
               
               if( tightened )
               {
                  tightenedround = TRUE;
                  (*nchgbds)++;
               }
               xlb = SCIPvarGetLbLocal(consdata->var);
            }
         }

         assert(!*cutoff);

         /* propagate bounds on y:
          *  (2) left hand side and upper bound on x -> bound on y
          */
         if( SCIPvarGetStatus(consdata->vbdvar) != SCIP_VARSTATUS_MULTAGGR && !SCIPisInfinity(scip, xub) ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               newlb = SCIPadjustedVarLb(scip, consdata->vbdvar, (consdata->lhs - xub)/consdata->vbdcoef);
               if( newlb > ylb + 0.5 )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", 
                     SCIPvarGetName(consdata->vbdvar), ylb, yub, newlb, yub);
                  SCIP_CALL( SCIPinferVarLbCons(scip, consdata->vbdvar, newlb, cons, (int)PROPRULE_2, FALSE,
                        cutoff, &tightened) );
                  
                  if( *cutoff )
                  {
                     assert(SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->vbdvar)));
                     
                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, PROPRULE_2, SCIP_BOUNDTYPE_LOWER) );
                     break;
                  }
                  
                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  ylb = SCIPvarGetLbLocal(consdata->vbdvar);
               }
            }
            else
            {
               newub = SCIPadjustedVarUb(scip, consdata->vbdvar, (consdata->lhs - xub)/consdata->vbdcoef);
               if( newub < yub - 0.5 )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", 
                     SCIPvarGetName(consdata->vbdvar), ylb, yub, ylb, newub);
                  SCIP_CALL( SCIPinferVarUbCons(scip, consdata->vbdvar, newub, cons, (int)PROPRULE_2, FALSE,
                        cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->vbdvar)));
                     
                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, PROPRULE_2, SCIP_BOUNDTYPE_UPPER) );
                     break;
                  }
                  
                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  yub = SCIPvarGetUbLocal(consdata->vbdvar);
               }
            }
         }
      }

      assert(!*cutoff);

      /* propagate right hand side inequality: x + c*y <= rhs */
      if( !SCIPisInfinity(scip, consdata->rhs) )
      {
         /* propagate bounds on x:
          *  (3) right hand side and bounds on y -> upper bound on x
          */
         if( SCIPvarGetStatus(consdata->var) != SCIP_VARSTATUS_MULTAGGR ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               if( !SCIPisInfinity(scip, -ylb) )
                  newub = SCIPadjustedVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * ylb);
               else 
                  newub = SCIPinfinity(scip);
            }
            else
            {
               if( !SCIPisInfinity(scip, yub) )
                  newub = SCIPadjustedVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * yub);
               else 
                  newub = SCIPinfinity(scip);
            }
            
            if( SCIPisUbBetter(scip, newub, xlb, xub) || ylb > yub - 0.5 )
            {
               SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                  SCIPvarGetName(consdata->var), xlb, xub, xlb, newub);
               SCIP_CALL( SCIPinferVarUbCons(scip, consdata->var, newub, cons, (int)PROPRULE_3, FALSE,
                     cutoff, &tightened) );

               if( *cutoff )
               {
                  assert(SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->var)));
                  
                  /* analyze infeasibility */
                  SCIP_CALL( analyzeConflict(scip, cons, consdata->var, PROPRULE_3, SCIP_BOUNDTYPE_UPPER) );
                  break;
               }
               
               if( tightened )
               {
                  tightenedround = TRUE;
                  (*nchgbds)++;
               }
               xub = SCIPvarGetUbLocal(consdata->var);
            }
         }

         assert(!*cutoff);

         /* propagate bounds on y:
          *  (4) right hand side and lower bound on x -> bound on y
          */
         if( SCIPvarGetStatus(consdata->vbdvar) != SCIP_VARSTATUS_MULTAGGR && !SCIPisInfinity(scip, -xlb) ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               newub = SCIPadjustedVarUb(scip, consdata->vbdvar, (consdata->rhs - xlb)/consdata->vbdcoef);
               if( newub < yub - 0.5 )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", 
                     SCIPvarGetName(consdata->vbdvar), ylb, yub, ylb, newub);
                  SCIP_CALL( SCIPinferVarUbCons(scip, consdata->vbdvar, newub, cons, (int)PROPRULE_4, FALSE,
                        cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->vbdvar)));
                  
                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, PROPRULE_4, SCIP_BOUNDTYPE_UPPER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  yub = SCIPvarGetUbLocal(consdata->vbdvar);
               }
            }
            else
            {
               newlb = SCIPadjustedVarLb(scip, consdata->vbdvar, (consdata->rhs - xlb)/consdata->vbdcoef);
               if( newlb > ylb + 0.5 )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", 
                     SCIPvarGetName(consdata->vbdvar), ylb, yub, newlb, yub);
                  SCIP_CALL( SCIPinferVarLbCons(scip, consdata->vbdvar, newlb, cons, (int)PROPRULE_4, FALSE,
                        cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->vbdvar)));
                     
                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, PROPRULE_4, SCIP_BOUNDTYPE_LOWER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  ylb = SCIPvarGetLbLocal(consdata->vbdvar);
               }
            }
         }
      }
      assert(!(*cutoff));
   }
   while( tightenedround );
   
   /* check for redundancy */
   if( !(*cutoff) && (SCIPisInfinity(scip, -consdata->lhs)
         || (consdata->vbdcoef > 0.0 && SCIPisFeasGE(scip, xlb + consdata->vbdcoef * ylb, consdata->lhs))
         || (consdata->vbdcoef < 0.0 && SCIPisFeasGE(scip, xlb + consdata->vbdcoef * yub, consdata->lhs)))
      && (SCIPisInfinity(scip, consdata->rhs)
         || (consdata->vbdcoef > 0.0 && SCIPisFeasLE(scip, xub + consdata->vbdcoef * yub, consdata->rhs))
         || (consdata->vbdcoef < 0.0 && SCIPisFeasLE(scip, xub + consdata->vbdcoef * ylb, consdata->rhs))) )
   {
      SCIPdebugMessage("variable bound constraint <%s> is redundant: <%s>[%.15g,%.15g], <%s>[%.15g,%.15g]\n",
         SCIPconsGetName(cons), 
         SCIPvarGetName(consdata->var), SCIPvarGetLbLocal(consdata->var), SCIPvarGetUbLocal(consdata->var),
         SCIPvarGetName(consdata->vbdvar), SCIPvarGetLbLocal(consdata->vbdvar), SCIPvarGetUbLocal(consdata->vbdvar));
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      if( ndelconss != NULL )
         (*ndelconss)++;
   }

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** replaces fixed and aggregated variables in variable bound constraint by active problem variables */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store whether an infeasibility was detected */
   int*                  nchgbds,            /**< pointer to count number of bound changes */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real varscalar;
   SCIP_Real varconstant;
   SCIP_VAR* vbdvar;
   SCIP_Real vbdvarscalar;
   SCIP_Real vbdvarconstant;
   SCIP_Bool varschanged;
   SCIP_Bool redundant;

   assert(cutoff != NULL);
   assert(nchgbds != NULL);
   assert(ndelconss != NULL);

   *cutoff = FALSE;
   redundant = FALSE;

   /* the variable bound constraint is: lhs <= x + c*y <= rhs */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get active problem variables of x and y */
   var = consdata->var;
   varscalar = 1.0;
   varconstant = 0.0;
   SCIP_CALL( SCIPvarGetProbvarSum(&var, &varscalar, &varconstant) );
   vbdvar = consdata->vbdvar;
   vbdvarscalar = 1.0;
   vbdvarconstant = 0.0;
   SCIP_CALL( SCIPvarGetProbvarSum(&vbdvar, &vbdvarscalar, &vbdvarconstant) );
   varschanged = (var != consdata->var || vbdvar != consdata->vbdvar);

   /**@todo fix bug: active variables might be MULTAGGR -> no bound changes possible! */

   /* if the variables are equal, the variable bound constraint reduces to standard bounds on the single variable */
   if( var == vbdvar && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_Real scalar;
      SCIP_Real constant;

      SCIPdebugMessage("variable bound constraint <%s> has equal variable and vbd variable <%s>\n",
         SCIPconsGetName(cons), SCIPvarGetName(var));

      /*      lhs <= a1*z + b1 + c(a2*z + b2) <= rhs
       * <=>  lhs <= (a1 + c*a2)z + (b1 + c*b2) <= rhs
       */
      scalar = varscalar + consdata->vbdcoef * vbdvarscalar;
      constant = varconstant + consdata->vbdcoef * vbdvarconstant;
      if( SCIPisZero(scip, scalar) )
      {
         /* no variable is left: the constraint is redundant or infeasible */
         if( SCIPisFeasLT(scip, constant, consdata->lhs) || SCIPisFeasGT(scip, constant, consdata->rhs) )
            *cutoff = TRUE;
      }
      else if( scalar > 0.0 )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;
            
            SCIP_CALL( SCIPtightenVarLb(scip, var, (consdata->lhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage(" -> tightened lower bound: <%s> >= %.15g\n", 
                  SCIPvarGetName(var), SCIPvarGetLbGlobal(var));
               (*nchgbds)++;
            }
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;
            
            SCIP_CALL( SCIPtightenVarUb(scip, var, (consdata->rhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage(" -> tightened upper bound: <%s> <= %.15g\n", 
                  SCIPvarGetName(var), SCIPvarGetUbGlobal(var));
               (*nchgbds)++;
            }
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;
            
            SCIP_CALL( SCIPtightenVarUb(scip, var, (consdata->lhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage(" -> tightened upper bound: <%s> <= %.15g\n", 
                  SCIPvarGetName(var), SCIPvarGetUbGlobal(var));
               (*nchgbds)++;
            }
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;
            
            SCIP_CALL( SCIPtightenVarLb(scip, var, (consdata->rhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage(" -> tightened lower bound: <%s> >= %.15g\n", 
                  SCIPvarGetName(var), SCIPvarGetLbGlobal(var));
               (*nchgbds)++;
            }
         }
      }
      redundant = TRUE;
   }
   else
   {
      /* if the variables should be replaced, drop the events and catch the events on the new variables afterwards */
      if( varschanged )
      {
         SCIP_CALL( dropEvents(scip, consdata) );
      }

      /* apply aggregation on x */
      if( SCIPisZero(scip, varscalar) )
      {
         SCIPdebugMessage("variable bound constraint <%s>: variable <%s> is fixed to %.15g\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->var), varconstant);

         /* cannot change bounds on multi-aggregated variables */
         if( SCIPvarGetStatus(vbdvar) != SCIP_VARSTATUS_MULTAGGR )
         {
            /* x is fixed to varconstant: update bounds of y and delete the variable bound constraint */
            if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
            {
               if( consdata->vbdcoef > 0.0 )
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarLb(scip, consdata->vbdvar, (consdata->lhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage(" -> tightened lower bound: <%s> >= %.15g\n", 
                        SCIPvarGetName(consdata->vbdvar), SCIPvarGetLbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
               else
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarUb(scip, consdata->vbdvar, (consdata->lhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage(" -> tightened upper bound: <%s> <= %.15g\n", 
                        SCIPvarGetName(consdata->vbdvar), SCIPvarGetUbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
            }
            if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
            {
               if( consdata->vbdcoef > 0.0 )
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarUb(scip, consdata->vbdvar, (consdata->rhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage(" -> tightened upper bound: <%s> <= %.15g\n", 
                        SCIPvarGetName(consdata->vbdvar), SCIPvarGetUbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
               else
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarLb(scip, consdata->vbdvar, (consdata->rhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage(" -> tightened lower bound: <%s> >= %.15g\n", 
                        SCIPvarGetName(consdata->vbdvar), SCIPvarGetLbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
            }
            redundant = TRUE;
         }
      }
      else if( var != consdata->var )
      {
         /* replace aggregated variable x in the constraint by its aggregation */
         if( varscalar > 0.0 )
         {
            /* lhs := (lhs - varconstant) / varscalar
             * rhs := (rhs - varconstant) / varscalar
             * c   := c / varscalar
             */
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs = (consdata->lhs - varconstant)/varscalar;
            if( !SCIPisInfinity(scip, consdata->rhs) )
               consdata->rhs = (consdata->rhs - varconstant)/varscalar;
            consdata->vbdcoef /= varscalar;
         }
         else
         {
            SCIP_Real lhs;
            
            assert(varscalar != 0.0);

            /* lhs := (rhs - varconstant) / varscalar
             * rhs := (lhs - varconstant) / varscalar
             * c   := c / varscalar
             */
            lhs = consdata->lhs;
            consdata->lhs = -consdata->rhs;
            consdata->rhs = -lhs;
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs = (consdata->lhs + varconstant)/(-varscalar);
            if( !SCIPisInfinity(scip, consdata->rhs) )
               consdata->rhs = (consdata->rhs + varconstant)/(-varscalar);
            consdata->vbdcoef /= varscalar;
         }
         consdata->var = var;
      }

      /* apply aggregation on y */
      if( SCIPisZero(scip, vbdvarscalar) )
      {
         SCIPdebugMessage("variable bound constraint <%s>: vbd variable <%s> is fixed to %.15g\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->vbdvar), vbdvarconstant);

         /* cannot change bounds on multi-aggregated variables */
         if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
         {
            /* y is fixed to vbdvarconstant: update bounds of x and delete the variable bound constraint */
            if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
            {
               SCIP_Bool tightened;

               SCIP_CALL( SCIPtightenVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * vbdvarconstant,
                     TRUE, cutoff, &tightened) );
               if( tightened )
               {
                  SCIPdebugMessage(" -> tightened lower bound: <%s> >= %.15g\n", 
                     SCIPvarGetName(consdata->var), SCIPvarGetLbGlobal(consdata->var));
                  (*nchgbds)++;
               }
            }
            if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
            {
               SCIP_Bool tightened;

               SCIP_CALL( SCIPtightenVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * vbdvarconstant,
                     TRUE, cutoff, &tightened) );
               if( tightened )
               {
                  SCIPdebugMessage(" -> tightened upper bound: <%s> <= %.15g\n", 
                     SCIPvarGetName(consdata->var), SCIPvarGetUbGlobal(consdata->var));
                  (*nchgbds)++;
               }
            }
            redundant = TRUE;
         }
      }
      else if( vbdvar != consdata->vbdvar )
      {
         /* replace aggregated variable y in the constraint by its aggregation:
          * lhs := lhs - c * vbdvarconstant
          * rhs := rhs - c * vbdvarconstant
          * c   := c * vbdvarscalar
          */
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhs -= consdata->vbdcoef * vbdvarconstant;
         if( !SCIPisInfinity(scip, consdata->rhs) )
            consdata->rhs -= consdata->vbdcoef * vbdvarconstant;
         consdata->vbdcoef *= vbdvarscalar;
         consdata->vbdvar = vbdvar;
      }

      /* catch the events again on the new variables */
      if( varschanged )
      {
         SCIP_CALL( catchEvents(scip, consdata) );
      }
   }

   /* delete a redundant constraint */
   if( !(*cutoff) && redundant )
   {
      SCIPdebugMessage(" -> variable bound constraint <%s> is redundant\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPdelCons(scip, cons) );
      (*ndelconss)++;
   }

   return SCIP_OKAY;
}

/** tightens variable bound coefficient by inspecting the global bounds of the involved variables;
 *  note: this is also performed by the linear constraint handler - only necessary if the user directly creates variable bound constraints
 */
static
SCIP_RETCODE tightenCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   int*                  nchgsides           /**< pointer to count the number of left and right hand sides */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real xlb;
   SCIP_Real xub;
   int oldnchgcoefs;
   int oldnchgsides;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* coefficient tightening only works for binary bound variable */
   if( !SCIPvarIsBinary(consdata->vbdvar) )
      return SCIP_OKAY;

   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   /* get bounds of variable x */
   xlb = SCIPvarGetLbGlobal(consdata->var);
   xub = SCIPvarGetUbGlobal(consdata->var);

   /* modification of coefficient can only be applied if only one side is finite */
   if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, consdata->rhs) )
   {
      /* lhs <= x + c*y  =>  x >= lhs - c*y */
      if( consdata->vbdcoef > 0.0 && SCIPisFeasGT(scip, xlb, consdata->lhs - consdata->vbdcoef) )
      {
         /* constraint has positive slack for the non-restricting case y = 1
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 1 and equivalent in the restricting case y = 0
          * -> c' = lhs - xlb
          */
         SCIPdebugMessage("tighten binary VLB <%s>[%.15g,%.15g] %+.15g<%s> >= %.15g to <%s> %+.15g<%s> >= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->lhs,
            SCIPvarGetName(consdata->var), consdata->lhs - xlb, SCIPvarGetName(consdata->vbdvar), consdata->lhs);
         consdata->vbdcoef = consdata->lhs - xlb;
         (*nchgcoefs)++;
      }
      else if( consdata->vbdcoef < 0.0 && SCIPisFeasGT(scip, xlb, consdata->lhs) )
      {
         /* constraint has positive slack for the non-restricting case y = 0
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 0 and equivalent in the restricting case y = 1
          * -> c' = c - lhs + xlb, lhs' = xlb
          */
         SCIPdebugMessage("tighten binary VLB <%s>[%.15g,%.15g] %+.15g<%s> >= %.15g to <%s> %+.15g<%s> >= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->lhs,
            SCIPvarGetName(consdata->var), consdata->vbdcoef - consdata->lhs + xlb, SCIPvarGetName(consdata->vbdvar), xlb);
         consdata->vbdcoef = consdata->vbdcoef - consdata->lhs + xlb;
         consdata->lhs = xlb;
         (*nchgcoefs)++;
         (*nchgsides)++;
      }
   }
   else if( SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* x + c*y <= rhs  =>  x <= rhs - c*y */
      if( consdata->vbdcoef < 0.0 && SCIPisFeasLT(scip, xub, consdata->rhs - consdata->vbdcoef) )
      {
         /* constraint has positive slack for the non-restricting case y = 1
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 1 and equivalent in the restricting case y = 0
          * -> c' = rhs - xub
          */
         SCIPdebugMessage("tighten binary VUB <%s>[%.15g,%.15g] %+.15g<%s> <= %.15g to <%s> %+.15g<%s> <= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs,
            SCIPvarGetName(consdata->var), consdata->rhs - xub, SCIPvarGetName(consdata->vbdvar), consdata->rhs);
         consdata->vbdcoef = consdata->rhs - xub;
         (*nchgcoefs)++;
      }
      else if( consdata->vbdcoef > 0.0 && SCIPisFeasLT(scip, xub, consdata->rhs) )
      {
         /* constraint has positive slack for the non-restricting case y = 0
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 0 and equivalent in the restricting case y = 1
          * -> c' = c - rhs + xub, rhs' = xub
          */
         SCIPdebugMessage("tighten binary VUB <%s>[%.15g,%.15g] %+.15g<%s> <= %.15g to <%s> %+.15g<%s> <= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs,
            SCIPvarGetName(consdata->var), consdata->vbdcoef - consdata->rhs + xub, SCIPvarGetName(consdata->vbdvar), xub);
         consdata->vbdcoef = consdata->vbdcoef - consdata->rhs + xub;
         consdata->rhs = xub;
         (*nchgcoefs)++;
         (*nchgsides)++;
      }
   }

   /* if something a coefficient or side of the varbound constraint was changed, ensure that the variable lower or
    * upper bounds of the variables are informed */
   if( *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
      consdata->addvarbounds = FALSE;

   return SCIP_OKAY;
}

/*
 * Linear constraint upgrading
 */

/** tries to upgrade a linear constraint into a variable bound constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdVarbound)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to a variable bound constraint  lhs <= x + a*y <= rhs
    * - there are exactly two variables
    * - one of the variables is non-binary (called the bounded variable x)
    * - one of the variables is non-continuous (called the bounding variable y)
    */
   upgrade = (nvars == 2) && (nposbin + nnegbin <= 1) && (nposcont + nnegcont <= 1);

   if( upgrade )
   {
      SCIP_VAR* var;
      SCIP_VAR* vbdvar;
      SCIP_Real vbdcoef;
      SCIP_Real vbdlhs;
      SCIP_Real vbdrhs;
      int vbdind;

      SCIPdebugMessage("upgrading constraint <%s> to variable bound constraint\n", SCIPconsGetName(cons));

      /* decide which variable we want to use as bounding variable y */
      if( SCIPvarGetType(vars[0]) < SCIPvarGetType(vars[1]) )
         vbdind = 0;
      else if( SCIPvarGetType(vars[0]) > SCIPvarGetType(vars[1]) )
         vbdind = 1;
      else if( SCIPisIntegral(scip, vals[0]) && !SCIPisIntegral(scip, vals[1]) )
         vbdind = 0;
      else if( !SCIPisIntegral(scip, vals[0]) && SCIPisIntegral(scip, vals[1]) )
         vbdind = 1;
      else if( REALABS(REALABS(vals[0]) - 1.0) < REALABS(REALABS(vals[1]) - 1.0) )
         vbdind = 0;
      else
         vbdind = 1;

      var = vars[1-vbdind];
      vbdvar = vars[vbdind];
      vbdcoef = vals[vbdind]/vals[1-vbdind];
      if( vals[1-vbdind] > 0.0 )
      {
         vbdlhs = SCIPisInfinity(scip, -lhs) ? -SCIPinfinity(scip) : lhs/vals[1-vbdind];
         vbdrhs = SCIPisInfinity(scip, rhs) ? SCIPinfinity(scip) : rhs/vals[1-vbdind];
      }
      else
      {
         vbdlhs = SCIPisInfinity(scip, rhs) ? -SCIPinfinity(scip) : rhs/vals[1-vbdind];
         vbdrhs = SCIPisInfinity(scip, -lhs) ? SCIPinfinity(scip) : lhs/vals[1-vbdind];
      }

      /* create the bin variable bound constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsVarbound(scip, upgdcons, SCIPconsGetName(cons), var, vbdvar, vbdcoef, vbdlhs, vbdrhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), 
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyVarbound)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrVarbound(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#define consFreeVarbound NULL


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitVarbound NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitVarbound NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreVarbound NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreVarbound NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolVarbound NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteVarbound)
{  /*lint --e{715}*/
   SCIP_CALL( consdataFree(scip, consdata) );
   
   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create target constraint data */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->var, sourcedata->vbdvar, sourcedata->vbdcoef, 
         sourcedata->lhs, sourcedata->rhs) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpVarbound)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      assert(SCIPconsIsInitial(conss[i]));
      SCIP_CALL( addRelaxation(scip, conss[i]) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpVarbound)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int i;

   *result = SCIP_DIDNOTFIND;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], NULL, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   /* separate remaining constraints */
   for( i = nusefulconss; i < nconss && *result == SCIP_DIDNOTFIND; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], NULL, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolVarbound)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int i;

   *result = SCIP_DIDNOTFIND;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], sol, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   /* separate remaining constraints */
   for( i = nusefulconss; i < nconss && *result == SCIP_DIDNOTFIND; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], sol, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpVarbound)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int i;

   *result = SCIP_FEASIBLE;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], NULL, FALSE) )
      {
         SCIP_CALL( separateCons(scip, conss[i], NULL, &separated) );
         if( separated )
         {
            *result = SCIP_SEPARATED;
            break;
         }
         else
            *result = SCIP_INFEASIBLE;
      }
   } 

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsVarbound)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], NULL, TRUE) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;  
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckVarbound)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], sol, checklprows) )
      {
         *result = SCIP_INFEASIBLE;

         if( printreason )
         {
            SCIP_Real sum;
            SCIP_CONSDATA* consdata;

            consdata = SCIPconsGetData(conss[i]);
            assert( consdata != NULL );
            
            sum = SCIPgetSolVal(scip, sol, consdata->var);
            sum += consdata->vbdcoef * SCIPgetSolVal(scip, sol, consdata->vbdvar);   
            
            SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
            if( !SCIPisFeasGE(scip, sum, consdata->lhs) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", consdata->lhs - sum);
            }
            if( !SCIPisFeasLE(scip, sum, consdata->rhs) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", sum - consdata->rhs);
            }
         }
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;
   
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropVarbound)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   int nchgbds;
   int i;

   cutoff = FALSE;
   nchgbds = 0;

   for( i = 0; i < nusefulconss && !cutoff; i++ )
   {
      SCIP_CALL( propagateCons(scip, conss[i], &cutoff, &nchgbds, NULL) );
   } 

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
   int oldnchgbds;
   int oldndelconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int i;

   cutoff = FALSE;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   for( i = 0; i < nconss && !cutoff && !SCIPisStopped(scip); i++ )
   {
      assert(!SCIPconsIsModifiable(conss[i]));

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->presolved = FALSE;

      if( consdata->presolved )
         continue;
      consdata->presolved = TRUE;

      /* make sure that the constraint is propagated */
      consdata->propagated = FALSE;

      /* incorporate fixings and aggregations in constraint */
      SCIP_CALL( applyFixings(scip, conss[i], &cutoff, nchgbds, ndelconss) );
      if( cutoff || !SCIPconsIsActive(conss[i]) )
         continue;

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, conss[i], &cutoff, nchgbds, ndelconss) );
      if( cutoff || !SCIPconsIsActive(conss[i]) )
         continue;

      /* tighten variable bound coefficient */
      SCIP_CALL( tightenCoefs(scip, conss[i], nchgcoefs, nchgsides) );

      /** informs once variable x about a globally valid variable lower or upper bound */
      if( !consdata->addvarbounds )
      {
         SCIP_Bool infeasible;
         int nlocalchgbds;
         
         nlocalchgbds = 0;

         /* if lhs is finite, we have a variable lower bound: lhs <= x + c*y  =>  x >= -c*y + lhs */
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIPdebugMessage("adding variable lower bound <%s> >= %g<%s> + %g\n", 
               SCIPvarGetName(consdata->var), -consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->lhs);

            SCIP_CALL( SCIPaddVarVlb(scip, consdata->var, consdata->vbdvar, -consdata->vbdcoef, consdata->lhs,
                  &infeasible, &nlocalchgbds) );
            assert(!infeasible);
            
            *nchgbds += nlocalchgbds;

            /* if lhs is finite, and x is not continuous we can add more variable bounds */
            if( SCIPvarGetType(consdata->var) != SCIP_VARTYPE_CONTINUOUS )
            {
               if( consdata->vbdcoef >= 0.0 )
               {
                  assert(consdata->vbdcoef != 0.0);

                  SCIPdebugMessage("adding variable lower bound <%s> >= %g<%s> + %g\n", 
                     SCIPvarGetName(consdata->vbdvar), -1.0/consdata->vbdcoef, SCIPvarGetName(consdata->var), 
                     consdata->lhs/consdata->vbdcoef);

                  /* if c > 0, we have a variable lower bound: lhs <= x + c*y  =>  y >= (lhs-x)/c */
                  SCIP_CALL( SCIPaddVarVlb(scip, consdata->vbdvar, consdata->var, 
                        -1.0/consdata->vbdcoef, consdata->lhs/consdata->vbdcoef, &infeasible, &nlocalchgbds) );
                  assert(!infeasible);

                  *nchgbds += nlocalchgbds;
               }
               else
               {
                  SCIPdebugMessage("adding variable upper bound <%s> <= %g<%s> + %g\n", 
                     SCIPvarGetName(consdata->vbdvar), -1.0/consdata->vbdcoef, SCIPvarGetName(consdata->var), 
                     consdata->lhs/consdata->vbdcoef);

                  /* if c < 0, we have a variable upper bound: lhs <= x + c*y  =>  y <= (lhs-x)/c */
                  SCIP_CALL( SCIPaddVarVub(scip, consdata->vbdvar, consdata->var, 
                        -1.0/consdata->vbdcoef, consdata->lhs/consdata->vbdcoef, &infeasible, &nlocalchgbds) );
                  assert(!infeasible);

                  *nchgbds += nlocalchgbds;
               }
            }
         }
         
         /* if rhs is finite, we have a variable upper bound: x + c*y <= rhs  =>  x <= -c*y + rhs */
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIPdebugMessage("adding variable upper bound <%s> <= %g<%s> + %g\n", 
               SCIPvarGetName(consdata->var), -consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs);

            SCIP_CALL( SCIPaddVarVub(scip, consdata->var, consdata->vbdvar, -consdata->vbdcoef, consdata->rhs,
                  &infeasible, &nlocalchgbds) );
            assert(!infeasible);

            *nchgbds += nlocalchgbds;

            /* if rhs is finite, and x is not continuous we can add more variable bounds */
            if( SCIPvarGetType(consdata->var) != SCIP_VARTYPE_CONTINUOUS )
            {
               if( consdata->vbdcoef > 0.0 )
               {
                  assert(consdata->vbdcoef != 0.0);

                  SCIPdebugMessage("adding variable upper bound <%s> <= %g<%s> + %g\n", 
                     SCIPvarGetName(consdata->vbdvar), -1.0/consdata->vbdcoef, SCIPvarGetName(consdata->var), 
                     consdata->rhs/consdata->vbdcoef);

                  /* if c > 0 we have a variable upper bound: x + c*y <= rhs  =>  y <= (rhs-x)/c */
                  SCIP_CALL( SCIPaddVarVub(scip, consdata->vbdvar, consdata->var, 
                        -1.0/consdata->vbdcoef, consdata->rhs/consdata->vbdcoef, &infeasible, &nlocalchgbds) );
                  assert(!infeasible);

                  *nchgbds += nlocalchgbds;
               }
               else
               {
                  SCIPdebugMessage("adding variable lower bound <%s> >= %g<%s> + %g\n", 
                     SCIPvarGetName(consdata->vbdvar), -1.0/consdata->vbdcoef, SCIPvarGetName(consdata->var), 
                     consdata->rhs/consdata->vbdcoef);
                  
                  /* if c < 0 we have a variable lower bound: x + c*y <= rhs  =>  y >= (rhs-x)/c */
                  SCIP_CALL( SCIPaddVarVlb(scip, consdata->vbdvar, consdata->var, 
                        -1.0/consdata->vbdcoef, consdata->rhs/consdata->vbdcoef, &infeasible, &nlocalchgbds) );
                  assert(!infeasible);

                  *nchgbds += nlocalchgbds;
               }
            }
         }
         consdata->addvarbounds = TRUE;
      }
   } 

   /**@todo preprocess pairs of variable bound constraints */

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nchgbds > oldnchgbds || *ndelconss > oldndelconss
      || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;
   
   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropVarbound)
{  /*lint --e{715}*/
   SCIP_CALL( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, boundtype, bdchgidx) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->var, nlockspos, nlocksneg) );
      if( consdata->vbdcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlockspos, nlocksneg) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlocksneg, nlockspos) );
      }
   }

   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->var, nlocksneg, nlockspos) );
      if( consdata->vbdcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlocksneg, nlockspos) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlockspos, nlocksneg) );
      }
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveVarbound NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveVarbound NULL


/** constraint enabling notification method of constraint handler */
#define consEnableVarbound NULL


/** constraint disabling notification method of constraint handler */
#define consDisableVarbound NULL


/** variable deletion method of constraint handler:
 *  varbound constraints are not modifiable and must have exactly two variables,
 *  so is is also not allowed to delete variables from them
 */
#define consDelVarsVarbound NULL


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);
   SCIPinfoMessage(scip, file, "<%s>[%c] %+.15g<%s>[%c]", SCIPvarGetName(consdata->var), 
      SCIPvarGetType(consdata->var) == SCIP_VARTYPE_BINARY ? 'B' :
      SCIPvarGetType(consdata->var) == SCIP_VARTYPE_INTEGER ? 'I' :
      SCIPvarGetType(consdata->var) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C',
      consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar),
      SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_BINARY ? 'B' :
      SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_INTEGER ? 'I' :
      SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');

   if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   
   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyVarbound)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   const char* consname;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, 2) );

   vars[0] = SCIPgetVarVarbound(sourcescip, sourcecons);
   vars[1] = SCIPgetVbdvarVarbound(sourcescip, sourcecons);

   coefs[0] = 1.0;
   coefs[1] = SCIPgetVbdcoefVarbound(sourcescip, sourcecons);

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   /* copy the varbound using the linear constraint copy method */
   SCIP_CALL( SCIPcopyConsLinear(scip, cons, sourcescip, consname, 2, vars, coefs,
         SCIPgetLhsVarbound(sourcescip, sourcecons), SCIPgetRhsVarbound(sourcescip, sourcecons), varmap, consmap, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
#define consParseVarbound NULL



/*
 * Event Handler
 */

/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   consdata->propagated = FALSE;
   consdata->presolved = FALSE;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for variable bound constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrVarbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   /* create variable bound constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyVarbound,
         consFreeVarbound, consInitVarbound, consExitVarbound, 
         consInitpreVarbound, consExitpreVarbound, consInitsolVarbound, consExitsolVarbound,
         consDeleteVarbound, consTransVarbound, consInitlpVarbound,
         consSepalpVarbound, consSepasolVarbound, consEnfolpVarbound, consEnfopsVarbound, consCheckVarbound, 
         consPropVarbound, consPresolVarbound, consRespropVarbound, consLockVarbound,
         consActiveVarbound, consDeactiveVarbound, 
         consEnableVarbound, consDisableVarbound,
         consDelVarsVarbound, consPrintVarbound, consCopyVarbound, consParseVarbound,
         conshdlrdata) );

   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint to varbound constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdVarbound, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }

   /* include event handler for bound change events */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecVarbound,
         eventhdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a variable bound constraint: lhs <= x + c*y <= rhs */
SCIP_RETCODE SCIPcreateConsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             var,                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   SCIP_Real             vbdcoef,            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs,                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs,                /**< right hand side of variable bound inequality */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the variable bound constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("variable bound constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, var, vbdvar, vbdcoef, lhs, rhs) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** gets left hand side of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_Real SCIPgetLhsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets right hand side of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_Real SCIPgetRhsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** gets bounded variable x of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_VAR* SCIPgetVarVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->var;
}

/** gets bounding variable y of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_VAR* SCIPgetVbdvarVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vbdvar;
}

/** gets bound coefficient c of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_Real SCIPgetVbdcoefVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vbdcoef;
}

/** gets the dual solution of the variable bound constraint in the current LP */
SCIP_Real SCIPgetDualsolVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets the dual Farkas value of the variable bound constraint in the current infeasible LP */
SCIP_Real SCIPgetDualfarkasVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}

/** returns the linear relaxation of the given variable bound constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->row;
}

