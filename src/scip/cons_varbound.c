/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_varbound.c,v 1.22 2005/01/31 12:20:57 bzfpfend Exp $"

/**@file   cons_varbound.c
 * @brief  constraint handler for varbound constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_varbound.h"
#include "cons_linear.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "varbound"
#define CONSHDLR_DESC          "variable bounds  lhs <= x + c*y <= rhs, x non-binary, y non-continuous"
#define CONSHDLR_SEPAPRIORITY   +900000
#define CONSHDLR_ENFOPRIORITY   -500000
#define CONSHDLR_CHECKPRIORITY  -500000
#define CONSHDLR_SEPAFREQ             1
#define CONSHDLR_PROPFREQ             1
#define CONSHDLR_EAGERFREQ          100
#define CONSHDLR_MAXPREROUNDS        -1
#define CONSHDLR_NEEDSCONS         TRUE

#define EVENTHDLR_NAME         "varbound"
#define EVENTHDLR_DESC         "bound change event handler for varbound constraints"

#define LINCONSUPGD_PRIORITY     +50000


/** varbound constraint data */
struct ConsData
{
   Real             vbdcoef;            /**< coefficient c of bounding variable y */
   Real             lhs;                /**< left hand side of variable bound inequality */
   Real             rhs;                /**< right hand side of variable bound inequality */
   VAR*             var;                /**< variable x that has variable bound */
   VAR*             vbdvar;             /**< binary, integer or implicit integer bounding variable y */
   ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   unsigned int     propagated:1;       /**< is the varbound constraint already propagated? */
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
RETCODE catchEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< varbound constraint data */
   )
{
   EVENTHDLR* eventhdlr;
   
   assert(consdata != NULL);

   eventhdlr = SCIPfindEventhdlr (scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* catch bound change events on variables */
   CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, 
         (EVENTDATA*)consdata, NULL) );
   CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vbdvar, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, 
         (EVENTDATA*)consdata, NULL) );
   
   return SCIP_OKAY;
}

/** drops events for variables */
static
RETCODE dropEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< varbound constraint data */
   )
{
   EVENTHDLR* eventhdlr;
   
   assert(consdata != NULL);

   eventhdlr = SCIPfindEventhdlr (scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* drop events on variables */
   CHECK_OKAY( SCIPdropVarEvent(scip, consdata->var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, 
         (EVENTDATA*)consdata, -1) );
   CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vbdvar, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
         (EVENTDATA*)consdata, -1) );

   return SCIP_OKAY;
}

/** creates a varbound constraint data object */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the varbound constraint data */
   VAR*             var,                /**< variable x that has variable bound */
   VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   Real             vbdcoef,            /**< coefficient c of bounding variable y */
   Real             lhs,                /**< left hand side of variable bound inequality */
   Real             rhs                 /**< right hand side of variable bound inequality */
   )
{
   assert(consdata != NULL);
   assert(SCIPvarGetType(vbdvar) != SCIP_VARTYPE_CONTINUOUS);

   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->var = var;
   (*consdata)->vbdvar = vbdvar;
   (*consdata)->vbdcoef = vbdcoef;
   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;
   (*consdata)->row = NULL;
   (*consdata)->propagated = FALSE;

   /* if we are in the transformed problem, get transformed variables, add variable bound information, and catch events */
   if( SCIPisTransformed(scip) )
   {
      CHECK_OKAY( SCIPgetTransformedVar(scip, (*consdata)->var, &(*consdata)->var) );
      CHECK_OKAY( SCIPgetTransformedVar(scip, (*consdata)->vbdvar, &(*consdata)->vbdvar) );

      /* if lhs is finite, we have a variable lower bound: lhs <= x + c*y  =>  x >= -c*y + lhs */
      if( !SCIPisInfinity(scip, -(*consdata)->lhs) )
      {
         CHECK_OKAY( SCIPaddVarVlb(scip, (*consdata)->var, (*consdata)->vbdvar, -(*consdata)->vbdcoef, (*consdata)->lhs) );
      }

      /* if rhs is finite, we have a variable upper bound: x + c*y <= rhs  =>  x <= -c*y + rhs */
      if( !SCIPisInfinity(scip, (*consdata)->rhs) )
      {
         CHECK_OKAY( SCIPaddVarVub(scip, (*consdata)->var, (*consdata)->vbdvar, -(*consdata)->vbdcoef, (*consdata)->rhs) );
      }

      /* catch events for variables */
      CHECK_OKAY( catchEvents(scip, *consdata) );
   }

   return SCIP_OKAY;
}   

/** frees a varbound constraint data */
static
RETCODE consdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata            /**< pointer to the varbound constraint */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* drop events */
   if( SCIPisTransformed(scip) )
   {
      CHECK_OKAY( dropEvents(scip, *consdata) );
   }

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** creates LP row corresponding to varbound constraint */
static 
RETCODE createRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< varbound constraint */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), consdata->lhs, consdata->rhs,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   CHECK_OKAY( SCIPaddVarToRow(scip, consdata->row, consdata->var, 1.0) );
   CHECK_OKAY( SCIPaddVarToRow(scip, consdata->row, consdata->vbdvar, consdata->vbdcoef) );

   return SCIP_OKAY;
}  

/** adds linear relaxation of varbound constraint to the LP */
static 
RETCODE addRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< varbound constraint */
   )
{
   CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      CHECK_OKAY( createRelaxation(scip, cons) );
   }
   assert(consdata->row != NULL);

   debugMessage("adding relaxation of varbound constraint <%s>: ", SCIPconsGetName(cons));
   debug( SCIProwPrint(consdata->row, NULL) );
   CHECK_OKAY( SCIPaddCut(scip, consdata->row, FALSE) );

   return SCIP_OKAY;
}

/** separates the given varbound constraint */
static
RETCODE separateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< varbound constraint */
   Bool*            separated           /**< pointer to store whether a cut was found */
   )
{
   CONSDATA* consdata;
   Real feasibility;

   assert(separated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("separating varbound constraint <%s>\n", SCIPconsGetName(cons));

   *separated = FALSE;

   /* create LP relaxation if not yet existing */
   if( consdata->row == NULL )
   {
      CHECK_OKAY( createRelaxation(scip, cons) );
   }

   /* check non-LP rows for feasibility and add them as cut, if violated */
   if( !SCIProwIsInLP(consdata->row) )
   {
      feasibility = SCIPgetRowLPFeasibility(scip, consdata->row);
      if( SCIPisFeasNegative(scip, feasibility) )
      {
         CHECK_OKAY( SCIPaddCut(scip, consdata->row, FALSE) );
         *separated = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** returns whether the given solution is feasible for the given varbound constraint */
static
Bool checkCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< varbound constraint */
   SOL*             sol,                /**< solution to check, NULL for current solution */
   Bool             checklprows         /**< should LP rows be checked? */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("checking varbound constraint <%s> for feasibility of solution %p (lprows=%d)\n",
      SCIPconsGetName(cons), sol, checklprows);

   if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
   {
      Real sum;

      sum = SCIPgetSolVal(scip, sol, consdata->var);
      sum += consdata->vbdcoef * SCIPgetSolVal(scip, sol, consdata->vbdvar);

      return SCIPisFeasGE(scip, sum, consdata->lhs) && SCIPisFeasLE(scip, sum, consdata->rhs);
   }
   else
      return TRUE;
}

/** propagation methode for varbound constraint */
static
RETCODE propagateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< varbound constraint */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   Bool*            redundant,          /**< pointer to store whether constraint is redundant */
   int*             nchgbds             /**< pointer to count number of bound changes */
   )
{
   CONSDATA* consdata;
   Real xlb;
   Real xub;
   Real ylb;
   Real yub;
   Real newlb;
   Real newub;
   Bool infeasible;
   Bool tightened;
   Bool tightenedround;

   assert(cutoff != NULL);
   assert(redundant != NULL);
   assert(nchgbds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("propagating varbound constraint <%s>: %g <= <%s> + %g<%s> <= %g\n",
      SCIPconsGetName(cons), consdata->lhs, SCIPvarGetName(consdata->var), consdata->vbdcoef,
      SCIPvarGetName(consdata->vbdvar), consdata->rhs);

   *cutoff = FALSE;
   *redundant = FALSE;

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
         /* propagate bounds on x:
          *  (1) left hand side and bounds on y -> lower bound on x
          */
         if( consdata->vbdcoef > 0.0 )
            newlb = SCIPadjustedVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * yub);
         else
            newlb = SCIPadjustedVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * ylb);
         if( newlb > xlb + 0.1 )
         {
            debugMessage(" -> tighten <%s>[%g,%g] -> [%g,%g]\n", SCIPvarGetName(consdata->var), xlb, xub, newlb, xub);
            CHECK_OKAY( SCIPinferVarLbCons(scip, consdata->var, newlb, cons, PROPRULE_1, &infeasible, &tightened) );
            *cutoff = *cutoff || infeasible;
            tightenedround = tightenedround || tightened;
            xlb = newlb;
            (*nchgbds)++;
         }

         /* propagate bounds on y:
          *  (2) left hand side and upper bound on x -> bound on y
          */
         if( consdata->vbdcoef > 0.0 )
         {
            newlb = SCIPadjustedVarLb(scip, consdata->vbdvar, (consdata->lhs - xub)/consdata->vbdcoef);
            if( newlb > ylb + 0.5 )
            {
               debugMessage(" -> tighten <%s>[%g,%g] -> [%g,%g]\n", 
                  SCIPvarGetName(consdata->vbdvar), ylb, yub, newlb, yub);
               CHECK_OKAY( SCIPinferVarLbCons(scip, consdata->vbdvar, newlb, cons, PROPRULE_2, &infeasible, &tightened) );
               *cutoff = *cutoff || infeasible;
               tightenedround = tightenedround || tightened;
               ylb = newlb;
               (*nchgbds)++;
            }
         }
         else
         {
            newub = SCIPadjustedVarUb(scip, consdata->vbdvar, (consdata->lhs - xub)/consdata->vbdcoef);
            if( newub < yub - 0.5 )
            {
               debugMessage(" -> tighten <%s>[%g,%g] -> [%g,%g]\n", 
                  SCIPvarGetName(consdata->vbdvar), ylb, yub, ylb, newub);
               CHECK_OKAY( SCIPinferVarUbCons(scip, consdata->vbdvar, newub, cons, PROPRULE_2, &infeasible, &tightened) );
               *cutoff = *cutoff || infeasible;
               tightenedround = tightenedround || tightened;
               yub = newub;
               (*nchgbds)++;
            }
         }
      }

      /* propagate right hand side inequality: x + c*y <= rhs */
      if( !SCIPisInfinity(scip, consdata->rhs) )
      {
         /* propagate bounds on x:
          *  (3) right hand side and bounds on y -> upper bound on x
          */
         if( consdata->vbdcoef > 0.0 )
            newub = SCIPadjustedVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * ylb);
         else
            newub = SCIPadjustedVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * yub);
         if( newub < xub - 0.1 )
         {
            debugMessage(" -> tighten <%s>[%g,%g] -> [%g,%g]\n", SCIPvarGetName(consdata->var), xlb, xub, xlb, newub);
            CHECK_OKAY( SCIPinferVarUbCons(scip, consdata->var, newub, cons, PROPRULE_3, &infeasible, &tightened) );
            *cutoff = *cutoff || infeasible;
            tightenedround = tightenedround || tightened;
            xub = newub;
            (*nchgbds)++;
         }

         /* propagate bounds on y:
          *  (4) right hand side and lower bound on x -> bound on y
          */
         if( consdata->vbdcoef > 0.0 )
         {
            newub = SCIPadjustedVarUb(scip, consdata->vbdvar, (consdata->rhs - xlb)/consdata->vbdcoef);
            if( newub < yub - 0.5 )
            {
               debugMessage(" -> tighten <%s>[%g,%g] -> [%g,%g]\n", 
                  SCIPvarGetName(consdata->vbdvar), ylb, yub, ylb, newub);
               CHECK_OKAY( SCIPinferVarUbCons(scip, consdata->vbdvar, newub, cons, PROPRULE_4, &infeasible, &tightened) );
               *cutoff = *cutoff || infeasible;
               tightenedround = tightenedround || tightened;
               yub = newub;
               (*nchgbds)++;
            }
         }
         else
         {
            newlb = SCIPadjustedVarLb(scip, consdata->vbdvar, (consdata->rhs - xlb)/consdata->vbdcoef);
            if( newlb > ylb + 0.5 )
            {
               debugMessage(" -> tighten <%s>[%g,%g] -> [%g,%g]\n", 
                  SCIPvarGetName(consdata->vbdvar), ylb, yub, newlb, yub);
               CHECK_OKAY( SCIPinferVarLbCons(scip, consdata->vbdvar, newlb, cons, PROPRULE_4, &infeasible, &tightened) );
               *cutoff = *cutoff || infeasible;
               tightenedround = tightenedround || tightened;
               ylb = newlb;
               (*nchgbds)++;
            }
         }
      }
   }
   while( !(*cutoff) && tightenedround );

   /* if one of the two variables is fixed, the constraint is redundant */
   if( SCIPisFeasEQ(scip, xlb, xub) || SCIPisFeasEQ(scip, ylb, yub) )
   {
      debugMessage("varbound constraint <%s> is redundant: <%s>[%g,%g], <%s>[%g,%g]\n",
         SCIPconsGetName(cons), 
         SCIPvarGetName(consdata->var), SCIPvarGetLbLocal(consdata->var), SCIPvarGetUbLocal(consdata->var),
         SCIPvarGetName(consdata->vbdvar), SCIPvarGetLbLocal(consdata->vbdvar), SCIPvarGetUbLocal(consdata->vbdvar));
      CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      *redundant = TRUE;
   }

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** resolves a propagation on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) left hand side and bounds on y -> lower bound on x
 *   (2) left hand side and upper bound on x -> bound on y
 *   (3) right hand side and bounds on y -> upper bound on x
 *   (4) right hand side and lower bound on x -> bound on y
 */
static
RETCODE resolvePropagation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint that inferred the bound change */
   VAR*             infervar,           /**< variable that was deduced */
   PROPRULE         proprule,           /**< propagation rule that deduced the bound change */
   BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   CONSDATA* consdata;

   assert(result != NULL);

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
         CHECK_OKAY( SCIPaddConflictUb(scip, consdata->vbdvar, bdchgidx) );
      }
      else
      {
         CHECK_OKAY( SCIPaddConflictLb(scip, consdata->vbdvar, bdchgidx) );
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_2:
      /* lhs <= x + c*y: left hand side and upper bound on x -> bound on y */
      assert(infervar == consdata->vbdvar);
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      CHECK_OKAY( SCIPaddConflictUb(scip, consdata->var, bdchgidx) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_3:
      /* x + c*y <= rhs: right hand side and bounds on y -> upper bound on x */
      assert(infervar == consdata->var);
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      assert(!SCIPisInfinity(scip, consdata->rhs));
      if( consdata->vbdcoef > 0.0 )
      {
         CHECK_OKAY( SCIPaddConflictLb(scip, consdata->vbdvar, bdchgidx) );
      }
      else
      {
         CHECK_OKAY( SCIPaddConflictUb(scip, consdata->vbdvar, bdchgidx) );
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_4:
      /* x + c*y <= rhs: right hand side and lower bound on x -> bound on y */
      assert(infervar == consdata->vbdvar);
      assert(!SCIPisInfinity(scip, consdata->rhs));
      CHECK_OKAY( SCIPaddConflictLb(scip, consdata->var, bdchgidx) );
      *result = SCIP_SUCCESS;
      break;

   default:
      errorMessage("invalid inference information %d in linear constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

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
DECL_CONSEXITSOL(consExitsolVarbound)
{
   CONSDATA* consdata;
   int c;

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         CHECK_OKAY( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteVarbound)
{  /*lint --e{715}*/
   CHECK_OKAY( consdataFree(scip, consdata) );
   
   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransVarbound)
{  /*lint --e{715}*/
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create target constraint data */
   CHECK_OKAY( consdataCreate(scip, &targetdata, sourcedata->var, sourcedata->vbdvar, sourcedata->vbdcoef, 
         sourcedata->lhs, sourcedata->rhs) );

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpVarbound)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( SCIPconsIsInitial(conss[i]) )
      {
         CHECK_OKAY( addRelaxation(scip, conss[i]) );
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler */
static
DECL_CONSSEPA(consSepaVarbound)
{  /*lint --e{715}*/
   Bool separated;
   int i;

   *result = SCIP_DIDNOTFIND;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss; ++i )
   {
      CHECK_OKAY( separateCons(scip, conss[i], &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   /* separate remaining constraints */
   for( i = nusefulconss; i < nconss && *result == SCIP_DIDNOTFIND; ++i )
   {
      CHECK_OKAY( separateCons(scip, conss[i], &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpVarbound)
{  /*lint --e{715}*/
   Bool separated;
   int i;

   *result = SCIP_FEASIBLE;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], NULL, FALSE) )
      {
         CHECK_OKAY( separateCons(scip, conss[i], &separated) );
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
DECL_CONSENFOPS(consEnfopsVarbound)
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
DECL_CONSCHECK(consCheckVarbound)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], sol, checklprows) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
DECL_CONSPROP(consPropVarbound)
{  /*lint --e{715}*/
   Bool cutoff;
   Bool redundant;
   int nchgbds;
   int i;

   cutoff = FALSE;
   nchgbds = 0;

   for( i = 0; i < nusefulconss && !cutoff; i++ )
   {
      CHECK_OKAY( propagateCons(scip, conss[i], &cutoff, &redundant, &nchgbds) );
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
DECL_CONSPRESOL(consPresolVarbound)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   Bool cutoff;
   Bool redundant;
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

   for( i = 0; i < nconss && !cutoff; i++ )
   {
      assert(!SCIPconsIsModifiable(conss[i]));

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->propagated = FALSE;

      if( consdata->propagated )
         continue;

      /* propagate constraint */
      CHECK_OKAY( propagateCons(scip, conss[i], &cutoff, &redundant, nchgbds) );
      if( redundant )
      {
         (*ndelconss)++;
         continue;
      }

      /**@todo tighten variable bound coefficient */
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
DECL_CONSRESPROP(consRespropVarbound)
{  /*lint --e{715}*/
   CHECK_OKAY( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, boundtype, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockVarbound)
{  /*lint --e{715}*/
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      CHECK_OKAY( SCIPaddVarLocks(scip, consdata->var, nlockspos, nlocksneg) );
      if( consdata->vbdcoef > 0.0 )
      {
         CHECK_OKAY( SCIPaddVarLocks(scip, consdata->vbdvar, nlockspos, nlocksneg) );
      }
      else
      {
         CHECK_OKAY( SCIPaddVarLocks(scip, consdata->vbdvar, nlocksneg, nlockspos) );
      }
   }

   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      CHECK_OKAY( SCIPaddVarLocks(scip, consdata->var, nlocksneg, nlockspos) );
      if( consdata->vbdcoef > 0.0 )
      {
         CHECK_OKAY( SCIPaddVarLocks(scip, consdata->vbdvar, nlocksneg, nlockspos) );
      }
      else
      {
         CHECK_OKAY( SCIPaddVarLocks(scip, consdata->vbdvar, nlockspos, nlocksneg) );
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

/** constraint display method of constraint handler */
static
DECL_CONSPRINT(consPrintVarbound)
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPisInfinity(scip, -consdata->lhs) )
      fprintf(file, "%g <= ", consdata->lhs);
   fprintf(file, "<%s> + %g<%s>", SCIPvarGetName(consdata->var), consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar));
   if( !SCIPisInfinity(scip, consdata->rhs) )
      fprintf(file, " <= %g", consdata->rhs);
   fprintf(file, "\n");

   return SCIP_OKAY;
}




/*
 * Linear constraint upgrading
 */

/** tries to upgrade a linear constraint into a varbound constraint */
static
DECL_LINCONSUPGD(linconsUpgdVarbound)
{  /*lint --e{715}*/
   Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to a variable bound constraint  lhs <= x + a*y <= rhs
    * - there are exactly two variables
    * - one of the variables is non-binary (called the bounded variable x)
    * - one of the variables is non-continuous (called the bounding variable y)
    */
   upgrade = (nvars == 2) && (nposbin + nnegbin <= 1) && (nposcont + nnegcont <= 1);

   if( upgrade )
   {
      VAR* var;
      VAR* vbdvar;
      Real vbdcoef;
      Real vbdlhs;
      Real vbdrhs;
      int vbdind;

      debugMessage("upgrading constraint <%s> to varbound constraint\n", SCIPconsGetName(cons));

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

      /* create the bin Varbound constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      CHECK_OKAY( SCIPcreateConsVarbound(scip, upgdcons, SCIPconsGetName(cons), var, vbdvar, vbdcoef, vbdlhs, vbdrhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   }

   return SCIP_OKAY;
}




/*
 * Event Handler
 */

/** execution methode of bound change event handler */
static
DECL_EVENTEXEC(eventExecVarbound)
{
   CONSDATA* consdata;

   consdata = (CONSDATA*)eventdata;
   assert(consdata != NULL);

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for varbound constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrVarbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;
   EVENTHDLRDATA* eventhdlrdata;

   /* create varbound constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, CONSHDLR_NEEDSCONS,
         consFreeVarbound, consInitVarbound, consExitVarbound, 
         consInitpreVarbound, consExitpreVarbound, consInitsolVarbound, consExitsolVarbound,
         consDeleteVarbound, consTransVarbound, consInitlpVarbound,
         consSepaVarbound, consEnfolpVarbound, consEnfopsVarbound, consCheckVarbound, 
         consPropVarbound, consPresolVarbound, consRespropVarbound, consLockVarbound,
         consActiveVarbound, consDeactiveVarbound, 
         consEnableVarbound, consDisableVarbound,
         consPrintVarbound,
         conshdlrdata) );

   /* include event handler for bound change events */
   eventhdlrdata = NULL;
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC, 
         NULL, NULL, NULL, NULL, eventExecVarbound,
         eventhdlrdata) );

   /* include the linear constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdVarbound, LINCONSUPGD_PRIORITY) );

   return SCIP_OKAY;
}

/** creates and captures a varbound constraint: lhs <= x + c*y <= rhs */
RETCODE SCIPcreateConsVarbound(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   VAR*             var,                /**< variable x that has variable bound */
   VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   Real             vbdcoef,            /**< coefficient c of bounding variable y */
   Real             lhs,                /**< left hand side of variable bound inequality */
   Real             rhs,                /**< right hand side of variable bound inequality */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   /* find the varbound constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("varbound constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, &consdata, var, vbdvar, vbdcoef, lhs, rhs) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, removeable) );

   return SCIP_OKAY;
}
