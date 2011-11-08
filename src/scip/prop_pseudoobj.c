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

/**@file   prop_pseudoobj.c
 * @brief  Pseudo objective propagator
 * @author Tobias Achterberg
 * @author Stefan Heinz
 *
 * This propagator propagates the objective function using the cutoff bound and the pseudo objective value. The pseudo
 * objective value can be seen as minimum activity of the linear objective function. Using this, this propagator checks
 * if variables with non-zero objective coefficients can exceed the cutoff bound. If this is the case the corresponding
 * bound can be tightened.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_pseudoobj.h"


#define PROP_NAME              "pseudoobj"
#define PROP_DESC              "pseudo objective function propagator"
#define PROP_TIMING             SCIP_PROPTIMING_ALWAYS
#define PROP_PRIORITY                 0 /**< propagator priority */ 
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_PRESOL_PRIORITY   +6000000 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOL_DELAY          TRUE /**< should presolving be delay, if other presolvers found reductions?  */
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */

#define EVENTHDLR_NAME         "pseudoobj"
#define EVENTHDLR_DESC         "bound change event handler for pseudo objective function propagator"

#define DEFAULT_MINUSELESS          100 /**< minimal number of successive none binary variable propagator whithout a
                                         *   bound reduction before aborted */
#define DEFAULT_MAXVARSFRAC         0.1 /**< maximal fraction of none binary variables with non-zero objective
                                         *   without a bound reduction before aborted */
#define DEFAULT_PROPFULLINROOT     TRUE /**< do we want to propagate full if we are propagating the root node, all variables */
#define DEFAULT_PROPCUTOFFBOUND    TRUE /**< propagate new cutoff bound directly globally */
#define DEFAULT_FORCE             FALSE /**< should the propagator be forced even active pricer are present? Note that
                                         *   can be done if it is known that the pseudo objective activity is given by
                                         *   the zero bound for all variables which are currently not present in the
                                         *   problem */
#define DEFAULT_MAXNEWVARS         1000 /**< number of variable added after the propgatore is reinitialized? */



/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for global bound change events */
   SCIP_VAR**            objvars;            /**< none binary variables with non-zero objective */
   SCIP_Real             lastlowerbound;     /**< last lower bound which propagated */
   SCIP_Real             cutoffbound;        /**< last cutoff bound used in presolving */
   SCIP_Real             glbpseudoobjval;    /**< last pseudo objective used in presolving */
   SCIP_Real             maxpseudoobjact;    /**< maximal global pseudo objective activity */
   SCIP_Real             maxvarsfrac;        /**< maximal fraction of none binary variables with non-zero objective
                                              *   without a bound reduction before aborted */
   int                   maxpseudoobjactinf; /**< number of coefficients contributing with infinite value to maxpseudoobjact */
   int                   nobjbinvars;        /**< number of binary variables with non-zero objective */
   int                   nobjvars;           /**< number of none binary variables with non-zero objective */
   int                   minuseless;         /**< minimal number of successive none binary variable propagator whithout
                                              *   a bound reduction before aborted */
   int                   lastvarnum;         /**< last none binary variable number that was looked at */
   int                   glbfirstnonfixed;   /**< index of first globally non-fixed binary variable */
   int                   firstnonfixed;      /**< index of first non-fixed binary variable */
   int                   nnewvars;           /**< counter for counting number of new variables added */
   int                   maxnewvars;         /**< number of variable added after the propgatore is reinitialized? */
   SCIP_Bool             glbpropagated;      /**< are global domains propagated */
   SCIP_Bool             propfullinroot;     /**< do we want to propagate full if we are propagating the root node, all variables */
   SCIP_Bool             propcutoffbound;    /**< propagate new cutoff bound directly globally */
   SCIP_Bool             force;              /**< should the propagator be forced even if active pricer are present? */
   SCIP_Bool             catchvaradded;      /**< do we catch the variable added event? */
};

/** compare variables w.r.t.
 *  (i)  the absolute value the objective coefficient;
 *  (ii) the locks which indicate most effect -- for the variables with a positive (negative) objective coefficient the
 *       down (up) lock is used since this lock indicates that tightened of the upper (lower) bound will triegger
 *       further domain propagations;
 * (iii) the other locks;
 * (iv)  variable problem index;
 */
static
SCIP_DECL_SORTPTRCOMP(varCompObj)
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int locks1;
   int locks2;

   var1 = (SCIP_VAR*)elem1;
   var2 = (SCIP_VAR*)elem2;

   assert(SCIPvarGetObj(var1) != 0.0);
   assert(SCIPvarGetObj(var2) != 0.0);

   /* first criteria is the objective coefficient */
   if( REALABS(SCIPvarGetObj(var1)) < REALABS(SCIPvarGetObj(var2)) )
      return -1;
   else if( REALABS(SCIPvarGetObj(var1)) > REALABS(SCIPvarGetObj(var2)) )
      return +1;

   /* second criteria the locks which indicate most effect */
   if( SCIPvarGetObj(var1) > 0.0 )
      locks1 = SCIPvarGetNLocksDown(var1);
   else
      locks1 = SCIPvarGetNLocksUp(var1);

   if( SCIPvarGetObj(var2) > 0.0 )
      locks2 = SCIPvarGetNLocksDown(var2);
   else
      locks2 = SCIPvarGetNLocksUp(var2);

   if( locks1 < locks2 )
      return -1;
   if( locks1 > locks2 )
      return 1;

   /* third criteria the other locks */
   if( SCIPvarGetObj(var1) > 0.0 )
      locks1 = SCIPvarGetNLocksUp(var1);
   else
      locks1 = SCIPvarGetNLocksDown(var1);

   if( SCIPvarGetObj(var2) >  0.0 )
      locks2 = SCIPvarGetNLocksUp(var2);
   else
      locks2 = SCIPvarGetNLocksDown(var2);

   if( locks1 < locks2 )
      return -1;
   if( locks1 > locks2 )
      return 1;

   /* forth criteria use the problem index */
   return SCIPvarCompare(var1, var2);
}

/*
 * Local methods
 */

/** drop events depending on the objective coefficient */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_VAR** vars;
   SCIP_EVENTHDLR* eventhdlr;
   int k;

   assert(scip != NULL);
   assert(propdata != NULL);

   eventhdlr = propdata->eventhdlr;
   assert(eventhdlr != NULL);

   vars = propdata->objvars;
   assert(vars != NULL || propdata->nobjvars == 0);

   /* drop global bound tighten events that decrease the maximal pseudo objective activity and release variables */
   for( k = 0; k < propdata->nobjbinvars; ++k )
   {
      SCIP_VAR* var;
      SCIP_Real obj;

      assert(vars != NULL);

      var = vars[k];
      assert(var != NULL);

      obj = SCIPvarGetObj(var);
      assert(!SCIPisZero(scip, obj));

      /* drop bound change event */
      if( obj > 0.0 )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_GUBCHANGED | SCIP_EVENTTYPE_BOUNDRELAXED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      }
      else
      {
         SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_BOUNDRELAXED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      }

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &vars[k]) );
   }
   assert(k == propdata->nobjbinvars);

   /* drop global bound tighten events that decrease the maximal pseudo objective activity and release variables */
   for( ; k < propdata->nobjvars; ++k )
   {
      SCIP_VAR* var;
      SCIP_Real obj;

      assert(vars != NULL);

      var = vars[k];
      assert(var != NULL);

      obj = SCIPvarGetObj(var);
      assert(!SCIPisZero(scip, obj));

      /* drop bound change event */
      if( obj > 0.0 )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_GUBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      }
      else
      {
         SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_GLBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      }

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &vars[k]) );
   }

   return SCIP_OKAY;
}

/** counts the number if variables with none zero objective coefficient */
static
int countNoneZeroObjVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars               /**< number of variables */
   )
{
   int count;
   int v;

   count = 0;

   for( v = 0; v < nvars; ++v )
   {
      SCIP_Real obj;

      assert(vars[v] != NULL);

      obj = SCIPvarGetObj(vars[v]);
      if( !SCIPisZero(scip, obj) )
         count++;
   }

   return count;
}

/** free propagator data */
static
SCIP_RETCODE propdataExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   /* drop events for the none binary variables */
   SCIP_CALL( dropVarEvents(scip, propdata) );

   /* free memory for non-zero objective variables */
   SCIPfreeMemoryArrayNull(scip, &propdata->objvars);

   propdata->nobjvars = 0;
   propdata->nobjbinvars = 0;
   propdata->maxpseudoobjact = SCIP_INVALID;
   propdata->maxpseudoobjactinf = 0;
   propdata->lastvarnum = -1;
   propdata->glbfirstnonfixed = 0;
   propdata->firstnonfixed = 0;
   propdata->nnewvars = 0;

   return SCIP_OKAY;
}

/** initializate the propgator */
static
SCIP_RETCODE propdataInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_VAR** vars;
   int nvars;
   int nbinvars;
   int nintvars;
   int ncontvars;
   int nobjvars;
   int nobjintvars;
   int nobjbinvars;
   int v;

   assert(scip != NULL);
   assert(propdata != NULL);

   eventhdlr = propdata->eventhdlr;

   /* get problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   nintvars = SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip);
   ncontvars = SCIPgetNContVars(scip);

   /* count binary variables with non-zero objective value */
   nobjbinvars = countNoneZeroObjVars(scip, vars, nbinvars);
   SCIPdebugMessage("There are %d binary variables in the objective function <%s>.\n", nobjbinvars, SCIPgetProbName(scip));

   /* count integer variables with non-zero objective value */
   nobjintvars = countNoneZeroObjVars(scip, &vars[nbinvars], nintvars);
   SCIPdebugMessage("There are %d interger variables (none binary) in the objective function.\n", nobjintvars);

   /* count continue variables with non-zero objective value */
   nobjvars = countNoneZeroObjVars(scip, &vars[nbinvars + nintvars], ncontvars);
   SCIPdebugMessage("There are %d continuous variables in the objective function.\n", nobjvars);

   nobjvars += nobjbinvars + nobjintvars;

   /* collect the variables with non-zero objective and catch global bound tighten events that decrease the
    * maximal pseudo objective activity
    */
   if( nobjvars > 0 )
   {
      int k;

      k = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->objvars, nobjvars) );

      SCIPdebugMessage("Store objective variables at address <%p>.\n", (void*)propdata->objvars);

      for( v = 0; v < nvars; v++ )
      {
         SCIP_Real obj;
         SCIP_VAR* var;

         var = vars[v];
         obj = SCIPvarGetObj(var);
         if( SCIPisZero(scip, obj) )
            continue;

         /* store and capture variable */
         assert(k < nobjvars);
         propdata->objvars[k] = var;
         SCIP_CALL( SCIPcaptureVar(scip, var) ) ;

         /* catch bound change events */
         if( obj > 0.0 )
         {
            if( k < nobjbinvars )
            {
               SCIP_CALL( SCIPcatchVarEvent(scip, var,
                     SCIP_EVENTTYPE_GUBCHANGED | SCIP_EVENTTYPE_BOUNDRELAXED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
            }
            else
            {
               SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_GUBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
            }
         }
         else
         {
            if( k < nobjbinvars )
            {
            SCIP_CALL( SCIPcatchVarEvent(scip, var,
                  SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_BOUNDRELAXED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
            }
            else
            {
               SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_GLBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
            }
         }
         k++;

         /* check if already visited all variable with non-zero objective coefficient */
         if( k == nobjvars )
            break;
      }

      /* sort binary variables with respect to their objective coefficient */
      SCIPsortDownPtr((void**)propdata->objvars, varCompObj, nobjbinvars);

      /* sort integer variables with respect to their objective coefficient */
      SCIPsortDownPtr((void**)(&propdata->objvars[nobjbinvars]), varCompObj, nobjintvars);

      /* sort continuous variables with respect to their objective coefficient */
      SCIPsortDownPtr((void**)(&propdata->objvars[nobjbinvars + nobjintvars]), varCompObj, nobjvars - nobjbinvars - nobjintvars);

      assert(k == nobjvars);

      SCIPdebugMessage("variables with non-zero objective coefficient: %5d binaries, %5d integers, %5d continuous\n", nobjbinvars, nobjintvars, nobjvars - nobjbinvars - nobjintvars);
   }

   propdata->nobjvars = nobjvars;
   propdata->nobjbinvars = nobjbinvars;
   propdata->maxpseudoobjact = SCIP_INVALID;
   propdata->maxpseudoobjactinf = 0;
   propdata->lastvarnum = nobjbinvars-1;
   propdata->glbfirstnonfixed = 0;
   propdata->firstnonfixed = 0;
   propdata->nnewvars = 0;

   return SCIP_OKAY;
}


/** resolves a propagation by supplying the variables whose bound changes increased the pseudo objective value */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Real             reqpseudoobjval,    /**< the required minactivity which has to be proven */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL for conflict analysis initialization */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real glbpseudoobjval;
   SCIP_Real obj;
   int nvars;
   int v;

   vars = propdata->objvars;
   nvars = propdata->nobjvars;
   assert(nvars > 0 && vars != NULL);

   /* we use the last stored global pseudo objective activity; Note that this is just a relaxation since global bounds
    * might be tighten since then.
    */
   glbpseudoobjval = propdata->glbpseudoobjval;
   assert(!SCIPisInfinity(scip, -glbpseudoobjval));

   /* the global pseudo objective activity should be smaller than required minactivity, otherwise SCIP should be
    * already stopped
    */
   assert(SCIPisLE(scip, glbpseudoobjval, reqpseudoobjval));

   SCIPdebugMessage("resolve propagation global pseudo objective <%g>, a required minactivity <%g>\n",
      glbpseudoobjval, reqpseudoobjval);

   reqpseudoobjval -= glbpseudoobjval;

   /* the variables responsible for the propagation are the ones with
    *  - obj > 0 and local lb > global lb
    *  - obj < 0 and local ub < global ub
    *
    *  (i) first collect all variables which contribute negatively to the pseudo objective activity (minactivity), hence
    *      adjusting the required minactivity
    * (ii) finally collect all variables which contribute positively to the pseudo objective activity (minactivity)
    *      until we reached the adjusted required minactivity
    */
   for( v = 0; v < nvars && SCIPisPositive(scip, reqpseudoobjval); ++v )
   {
      var = vars[v];
      if( var == infervar )
         continue;

      obj = SCIPvarGetObj(var);
      assert(!SCIPisZero(scip, obj));

      if( obj > 0.0 )
      {
         SCIP_Real loclb;
         SCIP_Real glblb;

         glblb = SCIPvarGetLbGlobal(var);
         loclb = SCIPvarGetLbAtIndex(var, bdchgidx, FALSE);
         assert(SCIPisFeasGE(scip, loclb, glblb));

         if( SCIPisGT(scip, loclb, glblb) )
         {
            SCIPdebugMessage("  add bound change <%s>[%g] >= <%g>\n", SCIPvarGetName(var), obj, loclb);
            SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );

            assert(SCIPisPositive(scip, (loclb - glblb) * obj));
            reqpseudoobjval -= (loclb - glblb) * obj;
         }
      }
      else
      {
         SCIP_Real locub;
         SCIP_Real glbub;

         glbub = SCIPvarGetUbGlobal(var);
         locub = SCIPvarGetUbAtIndex(var, bdchgidx, FALSE);
         assert(SCIPisFeasLE(scip, locub, glbub));

         if( SCIPisLT(scip, locub, glbub) )
         {
            SCIPdebugMessage("  add bound change <%s>[%g] <= <%g>\n", SCIPvarGetName(var), obj, locub);
            SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );

            assert(SCIPisPositive(scip, (locub - glbub) * obj));
            reqpseudoobjval -= (locub - glbub) * obj;
         }
      }
   }

   return SCIP_OKAY;
}

/** propagates the cutoff bound for the given variable */
static
SCIP_RETCODE propagateCutoffboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_VAR*             var,                /**< variable to propagate */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   int*                  nchgbds,            /**< pointer to store the number of changed bounds */
   SCIP_Bool             local               /**< propagate local bounds, otherwise global bounds */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(!SCIPisInfinity(scip, -pseudoobjval));
   assert(!SCIPisInfinity(scip, cutoffbound));
   assert(SCIPisLT(scip, pseudoobjval, cutoffbound) );

   if( local )
   {
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
   }
   else
   {
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
   }

   if( SCIPisFeasEQ(scip, lb, ub) )
      return SCIP_OKAY;

   obj = SCIPvarGetObj(var);
   assert(!SCIPisZero(scip, obj));

   if( obj > 0.0 )
   {
      SCIP_Real newub;

      newub = lb + (cutoffbound - pseudoobjval) / obj;

      if( local )
      {
         SCIP_CALL( SCIPinferVarUbProp(scip, var, newub, prop, 0, FALSE, &infeasible, &tightened) );
         assert(!infeasible);

         if( tightened ) /* might not be tightened due to numerical reasons */
         {
            SCIPdebugMessage(" -> new (local) upper bound of variable <%s>[%.10f,%.10f]: %.10f\n",
               SCIPvarGetName(var), lb, ub, newub);
            (*nchgbds)++;
         }
      }
      else
      {
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newub, FALSE, &infeasible, &tightened) );
         assert(!infeasible);

         if( tightened )
         {
            SCIPdebugMessage(" -> new (global) upper bound of variable <%s>[%.10f,%.10f]: %.10f\n",
               SCIPvarGetName(var), lb, ub, newub);

            (*nchgbds)++;
         }
      }
   }
   else
   {
      SCIP_Real newlb;

      newlb = ub + (cutoffbound - pseudoobjval) / obj;

      if( local )
      {
         SCIP_CALL( SCIPinferVarLbProp(scip, var, newlb, prop, 0, FALSE, &infeasible, &tightened) );
         assert(!infeasible);

         if( tightened ) /* might not be tightened due to numerical reasons */
         {
            SCIPdebugMessage(" -> new (local) lower bound of variable <%s>[%g,%g]: %g\n",
               SCIPvarGetName(var), lb, ub, newlb);

            (*nchgbds)++;
         }
      }
      else
      {
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newlb, FALSE, &infeasible, &tightened) );
         assert(!infeasible);

         if( tightened )
         {
            SCIPdebugMessage(" -> new (global) lower bound of variable <%s>[%g,%g]: %g\n",
               SCIPvarGetName(var), lb, ub, newlb);

            (*nchgbds)++;
         }
      }
   }

   return SCIP_OKAY;
}

/** globally propagate new cutoff bound or pseudo objective activity */
static
SCIP_RETCODE propagateGlobally(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop                /**< propagator */
   )
{
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* check if we have a new cutoff bound; in that case we globally propagate this new bound */
   if( propdata->propcutoffbound && !propdata->glbpropagated )
   {
      SCIP_VAR** objvars;
      SCIP_Real pseudoobjval;
      SCIP_Real cutoffbound;
      int nobjbinvars;
      int nobjvars;
      int nchgbds;
      int v;

      pseudoobjval = propdata->glbpseudoobjval;
      cutoffbound = propdata->cutoffbound;

      /* check if the global pseudo objective value (minimum activity of the objective function) is greater or equal to
       * the cutoff bound
       */
      if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
      {
         /* we are done with solving since a global pseudo activity is greater or equal to the cutoff bound */
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
         return SCIP_OKAY;
      }

      objvars = propdata->objvars;
      nobjbinvars = propdata->nobjbinvars;
      nobjvars = propdata->nobjvars;
      nchgbds = 0;

      if( !SCIPisInfinity(scip, -pseudoobjval) )
      {
         /* always propagate the binary variables completely */
         for( v = propdata->glbfirstnonfixed; v < nobjbinvars; ++v )
         {
            SCIP_VAR* var;

            var =  objvars[v];
            assert(var != NULL);

            /* check if the variables is already globally fixed; if so continue with the potential candidate */
            if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
               continue;

            /* try to tighten the variable bound */
            SCIP_CALL( propagateCutoffboundVar(scip, prop, var, cutoffbound, pseudoobjval, &nchgbds, FALSE) );

            /* the binary variables are sorted in non-increasing manner w.r.t. the absolute value of their objective
             * coefficient; These values are the increase in the pseudo objective activity we would get if we fix the
             * variable to its worse bound; hence we can stop if for a variable this potential increase is not enough
             * too exceed the cutoff bound; It is not enough to fix it.
             */
            if( SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5)
            {
               SCIPdebugMessage("interrupt global pseudo objective propagation w.r.t. cutoff bound <%.15g> for binary variables after %d from %d binary variables\n",
                  cutoffbound, v, nobjbinvars);
               break;
            }
         }
         propdata->glbfirstnonfixed = v;
         propdata->firstnonfixed = MAX(propdata->firstnonfixed, v);

         /* propagate the none binary variables completely */
         for( v = nobjbinvars; v < nobjvars; ++v )
         {
            SCIP_CALL( propagateCutoffboundVar(scip, prop, objvars[v], cutoffbound, pseudoobjval, &nchgbds, FALSE) );
         }

         propdata->glbpropagated = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** propagates the cutoff bound c*x <= cutoff */
static
SCIP_RETCODE propagateCutoffbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** objvars;
   SCIP_Real pseudoobjval;
   SCIP_Real cutoffbound;
   int nchgbds;
   int nobjbinvars;
   int nobjvars;
   int c;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   objvars = propdata->objvars;
   nobjvars = propdata->nobjvars;
   nobjbinvars = propdata->nobjbinvars;
   assert(nobjvars == 0 || objvars != NULL);
   assert(nobjbinvars <= nobjvars);

   /* nothing to do for empty objective */
   if( nobjvars == 0 )
      return SCIP_OKAY;

   /* get current pseudo objective value and cutoff bound */
   pseudoobjval = SCIPgetPseudoObjval(scip);
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;
   cutoffbound = SCIPgetCutoffbound(scip);
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;

   /* check pseudo objective value of the root node */
   if( SCIPgetDepth(scip) == 0 && pseudoobjval > propdata->glbpseudoobjval )
   {
      propdata->glbpropagated = FALSE;
      propdata->glbpseudoobjval = pseudoobjval;
   }

   /* check current cutoff bound */
   if( cutoffbound < propdata->cutoffbound )
   {
      propdata->glbpropagated = FALSE;
      propdata->cutoffbound = cutoffbound;
   }

   /* first globally propagate new cutoff bound or pseudo objective activity */
   SCIP_CALL( propagateGlobally(scip, prop) );

   /* check if the pseudo objective value (minimum activity of the objective function) is greater or equal to the cutoff
    * bound
    */
   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      SCIPdebugMessage("pseudo objective value %g exceeds cutoff bound %g\n", pseudoobjval, cutoffbound);

      /* check if we are not in the root node */
      if( SCIPgetDepth(scip) > 0 )
      {
         /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
         SCIP_CALL( SCIPinitConflictAnalysis(scip) );
         SCIP_CALL( resolvePropagation(scip, propdata, cutoffbound, NULL, NULL) );

         /* analyze the conflict */
         SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
      }
      *result = SCIP_CUTOFF;

      return SCIP_OKAY;
   }

   SCIPdebugMessage("propagating pseudo objective function (pseudoobj: %g, cutoffbound: %g)\n", pseudoobjval, cutoffbound);

   *result = SCIP_DIDNOTFIND;
   nchgbds = 0;

   /* always propagate the binary variables completely; note that the variables before the firstnonfixed indexed are
    * already locally fixed and those before glbfirstnonfixed are already globally fixed
    */
#ifndef NDEBUG
   /* check that the variables before glbfirstnonfixed are globally fixed */
   for( v = 0; v < propdata->glbfirstnonfixed; ++v )
   {
      SCIP_VAR* var;

      var =  objvars[v];
      assert(var != NULL);

      assert(SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5);
   }
   /* check that the variables before firstnonfixed are locally fixed */
   for( v = propdata->glbfirstnonfixed; v < propdata->firstnonfixed; ++v )
   {
      SCIP_VAR* var;

      var =  objvars[v];
      assert(var != NULL);

      assert(SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5);
   }
#endif
   for( v = propdata->firstnonfixed; v < nobjbinvars; ++v )
   {
      SCIP_VAR* var;

      var =  objvars[v];
      assert(var != NULL);

      /* check if the variable is already locally fixed; in that case we just continue with the next potential
       * candidate
       */
      if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5)
         continue;

      /* try to locally fix the variable bound */
      SCIP_CALL( propagateCutoffboundVar(scip, prop, var, cutoffbound, pseudoobjval, &nchgbds, TRUE) );

      /* the binary variables are sorted in non-increasing manner w.r.t. the absolute value of their objective
       * coefficient; These values are the increase in the pseudo objective activity we would get if we fix the variable
       * to its worse bound; hence we can stop if for a variable this potential increase is not enough too exceed the
       * cutoff bound; It is not enough to fix it.
       */
      if( SCIPvarGetLbLocal(var) < 0.5 && SCIPvarGetUbLocal(var) > 0.5)
      {
         SCIPdebugMessage("interrupt local pseudo objective propagation w.r.t. cutoff bound <%.15g> for binary variables after %d from %d binary variables\n",
            cutoffbound, v, nobjbinvars);
         break;
      }
   }
   propdata->firstnonfixed = v;

   /* tighten domains of none binary variables, if they would increase the pseudo objective value above the cutoff bound */
   if( propdata->propfullinroot && SCIPgetDepth(scip) == 0 )
   {
      for( v = nobjbinvars; v < nobjvars; ++v )
      {
         SCIP_CALL( propagateCutoffboundVar(scip, prop, objvars[v], cutoffbound, pseudoobjval, &nchgbds, FALSE) );
      }

      /* reset lastvarnum to first none binary variable */
      // propdata->lastvarnum = nobjbinvars-1;
   }
   else
   {
      int nmaxuseless;
      int nuseless;

      /* compute maximum number of useless propagations before aborting */
      nmaxuseless = MIN(MAX(propdata->minuseless, (int)propdata->maxvarsfrac*(nobjvars - nobjbinvars)), nobjvars - nobjbinvars);
      nuseless = 0;

      v = propdata->lastvarnum;
      for( c = 0; c < nobjvars - nobjbinvars && nuseless < nmaxuseless; ++c )
      {
         int oldnchgbds;

         v++;
         if( v >= nobjvars )
            v = nobjbinvars;

         oldnchgbds = nchgbds;
         SCIP_CALL( propagateCutoffboundVar(scip, prop, objvars[v], cutoffbound, pseudoobjval, &nchgbds, TRUE) );

         /* check if the domain of the variable was reduced */
         if( oldnchgbds == nchgbds )
            nuseless++;
      }
      propdata->lastvarnum = v;
   }

   /* check if we chanced bounds */
   if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}

/** recalculates the maximum objective pseudoactivity */
static
void calcMaxObjPseudoactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(propdata != NULL);

   /* get problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* calculate current max pseudo activity and largest contribution */
   propdata->maxpseudoobjact = 0.0;
   propdata->maxpseudoobjactinf = 0;
   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real obj;
      SCIP_Real contrib;

      obj = SCIPvarGetObj(vars[v]);
      if( SCIPisPositive(scip, obj) )
      {
         contrib = SCIPvarGetUbGlobal(vars[v]);
         if( !SCIPisInfinity(scip, contrib) )
            contrib *= obj;
      }
      else if( SCIPisNegative(scip, obj) )
      {
         contrib = SCIPvarGetLbGlobal(vars[v]);
         if( !SCIPisInfinity(scip, -contrib) )
            contrib *= obj;
         else
            contrib *= -1.0;
      }
      else
         continue;

      if( SCIPisInfinity(scip, contrib) )
         propdata->maxpseudoobjactinf++;
      else
         propdata->maxpseudoobjact += contrib;
   }
}

/** returns the pseudo objective activity */
static
SCIP_Real getMaxObjPseudoactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   /* if necessary, calculate the maximum pseudo objective activity */
   if( propdata->maxpseudoobjact == SCIP_INVALID ) /*lint !e777*/
      calcMaxObjPseudoactivity(scip, propdata);
   assert(propdata->maxpseudoobjact != SCIP_INVALID); /*lint !e777*/

   if( propdata->maxpseudoobjactinf > 0)
      return SCIPinfinity(scip);

   return propdata->maxpseudoobjact;
}

/** returns the residual pseudo objective activity without the given value */
static
SCIP_Real getMaxObjPseudoactivityResidualValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Real             contrib             /**< value to eliminate from pseudo objective activity */
   )
{
   SCIP_Real residual;

   assert(propdata != NULL);

   /* if necessary, calculate the maximum pseudo objective activity */
   if( propdata->maxpseudoobjact == SCIP_INVALID ) /*lint !e777*/
      calcMaxObjPseudoactivity(scip, propdata);
   assert(propdata->maxpseudoobjact != SCIP_INVALID); /*lint !e777*/

   if( SCIPisInfinity(scip, contrib) )
   {
      assert(propdata->maxpseudoobjactinf >= 1);
      /* check if this variable yields the only infinite contribution */
      if( propdata->maxpseudoobjactinf == 1 )
         residual = propdata->maxpseudoobjact;
      else
         residual = SCIPinfinity(scip);
   }
   else
   {
      /* check if there is an infinite contribution */
      if( propdata->maxpseudoobjactinf >= 1 )
         residual = SCIPinfinity(scip);
      else
         residual = propdata->maxpseudoobjact - contrib;
   }

   return residual;
}

/** returns the residual pseudo objective activity */
static
SCIP_Real getMaxObjPseudoactivityResidual(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var                 /**< variable to get residual activity for */
   )
{
   SCIP_Real obj;
   SCIP_Real contrib;

   assert(propdata != NULL);

   contrib = 0.0;
   obj = SCIPvarGetObj(var);
   if( SCIPisPositive(scip, obj) )
   {
      contrib = SCIPvarGetUbGlobal(var);
      if( !SCIPisInfinity(scip, contrib) )
         contrib *= obj;
   }
   else if( SCIPisNegative(scip, obj) )
   {
      contrib = SCIPvarGetLbGlobal(var);
      if( !SCIPisInfinity(scip, -contrib) )
         contrib *= obj;
      else
         contrib *= -1.0;
   }

   return getMaxObjPseudoactivityResidualValue(scip, propdata, contrib);
}

/** propagates the global domains of the given variable with non-negative objective coefficient against the lower bound
 * (c*x >= lowerbound)
 */
static
SCIP_RETCODE propagateLowerboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable to propagate */
   SCIP_Real             lowerbound,         /**< lower bound to use */
   SCIP_RESULT*          result              /**< pointer to store the result of the variable tightening */
   )
{
   SCIP_Real residual;
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   obj = SCIPvarGetObj(var);
   assert(!SCIPisZero(scip, obj));

   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);
   residual = getMaxObjPseudoactivityResidual(scip, propdata, var);

   if( SCIPisInfinity(scip, residual) )
      return SCIP_OKAY;

   if( obj > 0.0 )
   {
      SCIP_Real newlb;

      newlb = (lowerbound - residual) / obj;

      /* adjust variable bound w.r.t. integrality */
      newlb = SCIPadjustedVarLb(scip, var, newlb);

      SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newlb, FALSE, &infeasible, &tightened) );

      if( infeasible )
            *result = SCIP_CUTOFF;
      if( tightened ) /* might not be tightened due to numerical reasons */
      {
         SCIPdebugMessage(" -> new global lower bound of variable <%s>[%.10f,%.10f]: %.10f (obj=%g, residual=%g)\n",
            SCIPvarGetName(var), lb, ub, newlb, obj, residual);

         *result = SCIP_REDUCEDDOM;
      }
   }
   else
   {
      SCIP_Real newub;

      newub = (lowerbound - residual) / obj;

      SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newub, FALSE, &infeasible, &tightened) );

      if( infeasible )
         *result = SCIP_CUTOFF;
      if( tightened ) /* might not be tightened due to numerical reasons */
      {
         SCIPdebugMessage(" -> new global upper bound of variable <%s>[%.10f,%.10f]: %.10f (obj=%g, residual=%g)\n",
            SCIPvarGetName(var), lb, ub, newub, obj, residual);

         *result = SCIP_REDUCEDDOM;
      }
   }

   return SCIP_OKAY;
}

/** propagates the global lower (dual) bound c*x >= lowerbound */
static
SCIP_RETCODE propagateLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Real lowerbound;
   SCIP_Real maxactivity;
   SCIP_VAR** objvars;
   int nobjbinvars;
   int nobjvars;
   int k;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   objvars = propdata->objvars;
   nobjvars = propdata->nobjvars;
   nobjbinvars = propdata->nobjbinvars;
   assert(nobjvars == 0 || objvars != NULL);
   assert(nobjbinvars <= nobjvars);

   /* nothing to do for empty objective */
   if( nobjvars == 0 )
      return SCIP_OKAY;

   /* check if there is a chance to find a reduction */
   lowerbound = SCIPgetLowerbound(scip);

   if( SCIPisInfinity(scip, -lowerbound) )
      return SCIP_OKAY;

   /* if the lower bound did not change since the last propagation as well as the global bounds of the variables with a
    * non-zero objective coefficient we do nothing since there is no new information available
    */
   if( SCIPisLE(scip, lowerbound, propdata->lastlowerbound) && propdata->maxpseudoobjact < SCIP_INVALID)
      return SCIP_OKAY;

   maxactivity = getMaxObjPseudoactivity(scip, propdata);

   /* if more than one variable contribute an infinity to the maximal pseudo activity we can do nothing */
   assert(propdata->maxpseudoobjact < SCIP_INVALID);
   if( propdata->maxpseudoobjactinf > 1 )
      return SCIP_OKAY;

   for( k = 0; k < nobjbinvars && *result != SCIP_CUTOFF; ++k )
   {
      SCIP_VAR* var;
      SCIP_Real objval;

      var = objvars[k];
      assert(var != NULL);

      objval = SCIPvarGetObj(var);
      assert(!SCIPisZero(scip, objval));

      SCIP_CALL( propagateLowerboundVar(scip, propdata, objvars[k], lowerbound, result) );

      /* the binary variables are sorted in non-increasing manner w.r.t. the absolute value of their objective
       * coefficient; These values are the decrease we would get with the maximal pseudo objective activity if we fix
       * the variable to its best bound; hence we can stop if for a variable this potential decrease is not enough
       * anymore too fall below the lower bound.
       */
      if( propdata->maxpseudoobjactinf == 0 && SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5 )
      {
         assert(!SCIPisInfinity(scip, maxactivity));
         SCIPdebugMessage("interrupt pseudo objective propagation w.r.t. lower bound <%.15g> for binary variables after %d from %d binary variables\n", lowerbound, k, nobjbinvars);
         break;
      }
   }

   /* propagate c*x >= lowerbound */
   for( k = nobjbinvars; k < nobjvars && *result != SCIP_CUTOFF; ++k )
   {
      SCIP_CALL( propagateLowerboundVar(scip, propdata, objvars[k], lowerbound, result) );
   }

   if( *result == SCIP_CUTOFF )
   {
      /* we are done with solving since a global bound change is infeasible: cutoff root node */
      SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
   }

   /* remember the lower bound which we already propagated */
   propdata->lastlowerbound = lowerbound;

   return SCIP_OKAY;
}



/*
 * Callback methods of propagator
 */

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyPseudoobj)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of propagator */
   SCIP_CALL( SCIPincludePropPseudoobj(scip) );
 
   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreePseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   SCIPfreeMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}


/** initialization method of propagator (called after problem was transformed) */
#define propInitPseudoobj NULL


/** deinitialization method of propagator (called before transformed problem is freed) */
#define propExitPseudoobj NULL


/** presolving initialization method of propagator (called when presolving is about to begin) */
#define propInitprePseudoobj NULL


/** presolving deinitialization method of propagator (called after presolving has been finished) */
#define propExitprePseudoobj NULL


/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolPseudoobj)
{
   SCIP_PROPDATA* propdata;

   /* if a pricer is active we can do  nothing */
   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* do nothing if active pricer are present and force flag is not TRUE */
   if( !propdata->force && SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* initialize the propagator data structure */
   SCIP_CALL( propdataInit(scip, propdata) );

   /* if active pricer are present we want to catch the variable added event */
   if( SCIPgetNActivePricers(scip) > 0 )
   {
      assert(!propdata->catchvaradded);
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, propdata->eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
      propdata->catchvaradded = TRUE;
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolPseudoobj)
{
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   if( propdata->catchvaradded )
   {
      /* drop the variable added event */
      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, propdata->eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      propdata->catchvaradded = FALSE;
   }

   /* free propagator data */
   SCIP_CALL( propdataExit(scip, propdata) );

   return SCIP_OKAY;
}


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolPseudoobj)
{  /*lint --e{715}*/

   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_Real cutoffbound;
   SCIP_Real pseudoobjval;
   int oldnchgbds;
   int nvars;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* do nothing if active pricer are present and force flag is not TRUE */
   if( !propdata->force && SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   pseudoobjval = SCIPgetPseudoObjval(scip);
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;

   cutoffbound = SCIPgetCutoffbound(scip);
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;

   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if( cutoffbound < propdata->cutoffbound || pseudoobjval > propdata->glbpseudoobjval )
   {
      *result = SCIP_DIDNOTFIND;
      oldnchgbds = *nchgbds;

      /* get the problem variables */
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      /* scan the variables for pseudoobj bound reductions
       * (loop backwards, since a variable fixing can change the current and the subsequent slots in the vars array)
       */
      for( v = nvars - 1; v >= 0; --v )
      {
         if( SCIPisZero(scip, SCIPvarGetObj(vars[v])) )
            continue;

         SCIP_CALL( propagateCutoffboundVar(scip, prop, vars[v], cutoffbound, pseudoobjval, nchgbds, FALSE) );
      }

      if( *nchgbds > oldnchgbds )
         *result = SCIP_SUCCESS;

      /* store the old values */
      propdata->cutoffbound = cutoffbound;
      propdata->glbpseudoobjval = pseudoobjval;
      propdata->glbpropagated = TRUE;
   }

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   *result = SCIP_DIDNOTRUN;

   /* do nothing if active pricer are present and force flag is not TRUE */
   if( !propdata->force && SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* check if the enough new variable are added (due to cilumn generatition to reinitialized the propgator data */
   if( propdata->nnewvars > propdata->maxnewvars )
   {
      /* free current propdata data */
      SCIP_CALL( propdataExit(scip, propdata) );

      /* initialize propdata data from scratch */
      SCIP_CALL( propdataInit(scip, propdata) );
   }

   /* propagate c*x <= cutoff */
   SCIP_CALL( propagateCutoffbound(scip, prop, result) );

   if( *result != SCIP_CUTOFF && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_RESULT dualresult;

      /* propagate c*x >= lowerbound */
      SCIP_CALL( propagateLowerbound(scip, prop, &dualresult) );

      if( dualresult == SCIP_REDUCEDDOM || dualresult == SCIP_CUTOFF )
         *result = dualresult;
   }

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Real reqpseudoobjval;
   SCIP_Real glbbound;
   SCIP_Real newbound;
   SCIP_Real obj;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   reqpseudoobjval = SCIPgetCutoffbound(scip);
   assert(!SCIPisInfinity(scip, reqpseudoobjval));
   assert(infervar != NULL);

   obj = SCIPvarGetObj(infervar);
   assert(!SCIPisZero(scip, obj));

   if( obj > 0.0 )
   {
      newbound = SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE);
      glbbound = SCIPvarGetLbGlobal(infervar);
   }
   else
   {
      newbound = SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE);
      glbbound = SCIPvarGetUbGlobal(infervar);
   }

   /* in case the variable is integral we just need to prove the newbound plus/minus (1 - epsilon) since the this bound
    * would be rounded to newbound due to integrability which is global information
    */
   if( SCIPvarIsIntegral(infervar) )
   {
      if( obj > 0.0 )
         newbound += 1 - 10 * SCIPfeastol(scip);
      else
         newbound -= 1 - 10 * SCIPfeastol(scip);
   }

   /* adjust the cutoff bound by the portion the inference variable contributes to the presudo objective activitiy
    * (minactivity)
    */
   reqpseudoobjval -= obj * (newbound - glbbound);

   /* resolve the propagation of the inference variable w.r.t. required minactivity */
   SCIP_CALL( resolvePropagation(scip, propdata, reqpseudoobjval, infervar, bdchgidx) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}




/*
 * Event handler
 */

/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_EVENTTYPE eventtype;

   /* if a pricer is active we can do  nothing */
   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   propdata = (SCIP_PROPDATA*)eventdata;
   assert(propdata != NULL);

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   eventtype = SCIPeventGetType(event);

   switch( eventtype )
   {
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBRELAXED:
      /* if bound got relaxed we need to start up front for trial of bound tightening */
      propdata->firstnonfixed = 0;
      break;
   case SCIP_EVENTTYPE_VARADDED:
      propdata->nnewvars++;
      break;
   default:
      assert(eventtype == SCIP_EVENTTYPE_GLBCHANGED || eventtype == SCIP_EVENTTYPE_GUBCHANGED);

      /* invalidate the maximum pseudo objective activity */
      propdata->maxpseudoobjact = SCIP_INVALID;
      propdata->maxpseudoobjactinf = 0;
   }

   return SCIP_OKAY;
}




/*
 * propagator specific interface methods
 */

/** creates the pseudo objective function propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropPseudoobj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;

   /* include event handler for gloabl bound change events and variable added event (in case of pricing) */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, eventExecPseudoobj, NULL) );

   /* create pseudoobj propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );
   propdata->objvars = NULL;
   propdata->nobjbinvars = 0;
   propdata->nobjvars = 0;
   propdata->maxpseudoobjact = SCIP_INVALID;
   propdata->maxpseudoobjactinf = 0;
   propdata->lastvarnum = -1;
   propdata->glbpropagated = FALSE;
   propdata->cutoffbound = SCIPinfinity(scip);
   propdata->lastlowerbound = -SCIPinfinity(scip);
   propdata->glbpseudoobjval = -SCIPinfinity(scip);
   propdata->firstnonfixed = 0;
   propdata->nnewvars = 0;
   propdata->catchvaradded = FALSE;

   /* get event handler for bound change events */
   propdata->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( propdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for pseudo objective propagator not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOL_DELAY,
         propCopyPseudoobj, propFreePseudoobj, propInitPseudoobj, propExitPseudoobj, propInitprePseudoobj, propExitprePseudoobj, 
         propInitsolPseudoobj, propExitsolPseudoobj, propPresolPseudoobj, propExecPseudoobj, propRespropPseudoobj,
         propdata) );

   /* add pseudoobj propagator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/"PROP_NAME"/minuseless",
         "minimal number of successive none binary variable propagator whithout a bound reduction before aborted",
         &propdata->minuseless, TRUE, DEFAULT_MINUSELESS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "propagating/"PROP_NAME"/maxvarsfrac",
         "maximal fraction of none binary variables with non-zero objective without a bound reduction before aborted",
         &propdata->maxvarsfrac, TRUE, DEFAULT_MAXVARSFRAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/propfullinroot",
         "do we want to propagate full if we are propagating the root node, despite the number of maxcand",
         &propdata->propfullinroot, TRUE, DEFAULT_PROPFULLINROOT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/propcutoffbound",
         "propagate new cutoff bound directly globally",
         &propdata->propcutoffbound, TRUE, DEFAULT_PROPCUTOFFBOUND, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/force",
         "should the propagator be forced even active pricer are present?",
         &propdata->force, TRUE, DEFAULT_FORCE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/"PROP_NAME"/maxnewvars",
         "number of variable added after the propgatore is reinitialized?",
         &propdata->maxnewvars, TRUE, DEFAULT_MAXNEWVARS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** propagates the cutoff bound for the given variables */
SCIP_RETCODE SCIPpropagateCutoffboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator, or NULL */
   SCIP_VAR*             var,                /**< variables to propagate */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   int*                  nchgbds             /**< pointer to store the number of changed bounds */
   )
{
   SCIP_CALL( propagateCutoffboundVar(scip, prop, var, cutoffbound, pseudoobjval, nchgbds, FALSE) );

   return SCIP_OKAY;
}
