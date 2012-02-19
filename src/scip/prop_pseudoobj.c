/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
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
 * @todo use the complete implication to initialize the objective implication data structure
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_pseudoobj.h"


#define PROP_NAME              "pseudoobj"
#define PROP_DESC              "pseudo objective function propagator"
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP | SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY           3000000 /**< propagator priority */
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
#define DEFAULT_PROPFULLINROOT     TRUE /**< do we want to propagate all non-binary variables if we are propagating the root node? */
#define DEFAULT_PROPCUTOFFBOUND    TRUE /**< propagate new cutoff bound directly globally */
#define DEFAULT_FORCE             FALSE /**< should the propagator be forced even active pricer are present? Note that
                                         *   can be done if it is known that the pseudo objective activity is given by
                                         *   the zero bound for all variables which are currently not present in the
                                         *   problem */
#define DEFAULT_MAXNEWVARS         1000 /**< number of variable added after the propagator is reinitialized? */
#define DEFAULT_PROPUSEIMPLICS     TRUE /**< use implications to strengthen the propagation of binary variable (increasing the objective change)? */
#define DEFAULT_RESPROPUSEIMPLICS  TRUE /**< use implications to strengthen the resolve propagation of binary variable (increasing the objective change)? */
#define DEFAULT_MAXIMPLVARS       50000 /**< maximum number of binary variables the implications are used if turned on (-1: unlimited)? */


/*
 * Data structures
 */

/** implication data structure for objective contributions of a binary variable */
struct SCIP_ObjImplics
{
   SCIP_VAR**            objvars;            /**< variables y in implications y == 0 or y == 1, first we store the
                                              *   implications by x == 0 and second the implications x == 1 */
   SCIP_Real             maxobjchg;          /**< maximum objective contribution if variables x is fixed to zero or one */
   int                   nlbimpls;           /**< number of all implications result through for x == 0 */
   int                   nubimpls;           /**< number of all implications result through for x == 1 */
   int                   size;               /**< size of the objvars array */
};
typedef struct SCIP_ObjImplics SCIP_OBJIMPLICS; /**< implications in the form x == 0 or x == 1 ==> y == 0 or y == 1 for (x and y binary) */


/** propagator data */
struct SCIP_PropData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for global bound change events */
   SCIP_VAR**            minactvars;         /**< binary variables with non-zero objective contribution w.r.t. minimum activity of the objective function */
   SCIP_OBJIMPLICS**     minactimpls;        /**< implication data structure for the binary variables w.r.t. minimum activity */
   SCIP_VAR**            maxactvars;         /**< binary variables with non-zero objective contribution w.r.t. maximum activity of the objective function */
   SCIP_Real*            maxactchgs;         /**< the maximal potential change of the objective if the binary variable
                                              *   is fixed to its best bound w.r.t. maximum activity of the objective function */
   SCIP_VAR**            objintvars;         /**< none binary variable with non-zero objective coefficient */
   SCIP_HASHTABLE*       addedvars;          /**< hash table used during resolving of a bound change (conflict analysis) */
   SCIP_Real             lastlowerbound;     /**< last lower bound which was propagated */
   SCIP_Real             cutoffbound;        /**< last cutoff bound used for propagation */
   SCIP_Real             glbpseudoobjval;    /**< last global pseudo objective used in presolving */
   SCIP_Real             maxvarsfrac;        /**< maximal fraction of none binary variables with non-zero objective
                                              *   without a bound reduction before aborted */
   SCIP_Real             maxpseudoobjact;    /**< maximal global pseudo objective activity */
   int                   maxpseudoobjactinf; /**< number of coefficients contributing with infinite value to maxpseudoobjact */
   int                   nminactvars;        /**< number of binary variables with non-zero objective contribution w.r.t. minimum activity of the objective function */
   int                   nmaxactvars;        /**< number of binary variables with non-zero objective contribution w.r.t. maximum activity of the objective function */
   int                   nobjintvars;        /**< number of none binary variables with non-zero objective */
   int                   minuseless;         /**< minimal number of successive none binary variable propagator whithout
                                              *   a bound reduction before aborted */
   int                   lastvarnum;         /**< last none binary variable number that was looked at */
   int                   glbfirstnonfixed;   /**< index of first globally non-fixed binary variable in minactvars array */
   int                   maxactfirstnonfixed;/**< index of first globally non-fixed binary variable in maxctvars array */
   int                   firstnonfixed;      /**< index of first locally non-fixed binary variable in minactvars array */
   int                   nnewvars;           /**< counter for counting number of new variables added */
   int                   maxnewvars;         /**< number of variable added after the propagator is reinitialized? */
   int                   maximplvars;        /**< maximum number of binary variables the implications are used if turned on (-1: unlimited)? */
   SCIP_Bool             glbpropagated;      /**< are global domains propagated */
   SCIP_Bool             propfullinroot;     /**< do we want to propagate all non-binary variables if we are propagating the root node */
   SCIP_Bool             propcutoffbound;    /**< propagate new cutoff bound directly globally */
   SCIP_Bool             force;              /**< should the propagator be forced even if active pricer are present? */
   SCIP_Bool             catchvaradded;      /**< do we catch the variable added event? */
   SCIP_Bool             propuseimplics;     /**< use implications to strengthen the propagation of binary variable (increasing the objective change)? */
   SCIP_Bool             respropuseimplics;  /**< use implications to strengthen the resolve propagation of binary variable (increasing the objective change)? */
};

/*
 * Debug methods
 */

#ifndef NDEBUG
/** check that the implications are applied for a globally fixed variable */
static
void checkImplicsApplied(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to check the implications */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* bounds;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Bool varfixing;
   int nvars;
   int v;

   /* check that the given variable is locally or globally fixed */
   assert(SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5);

   /* get fixed value */
   varfixing = SCIPvarGetLbGlobal(var) > 0.5;

   vars = SCIPvarGetImplVars(var, varfixing);
   nvars = SCIPvarGetNImpls(var, varfixing);
   bounds = SCIPvarGetImplBounds(var, varfixing);
   boundtypes = SCIPvarGetImplTypes(var, varfixing);

   /* check that each implication was applied */
   for( v = 0; v < nvars; ++v )
   {
      if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_Real lb;

         lb = SCIPvarGetLbGlobal(vars[v]);
         assert(SCIPisLE(scip, lb, bounds[v]));
      }
      else
      {
         SCIP_Real ub;

         assert(boundtypes[v] == SCIP_BOUNDTYPE_UPPER);

         ub = SCIPvarGetLbGlobal(vars[v]);
         assert(SCIPisGE(scip, ub, bounds[v]));
      }
   }
}

/** check if the global fixed indices are correct */
static
void checkGlbfirstnonfixed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_VAR* var;
   int v;

   for( v = 0; v < propdata->glbfirstnonfixed; ++v )
   {
      var = propdata->minactvars[v];
      assert(var != NULL);

      assert(SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5);
   }

   for( v = 0; v < propdata->maxactfirstnonfixed; ++v )
   {
      var = propdata->maxactvars[v];
      assert(var != NULL);

      assert(SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5);
   }
}
#endif /* end of debug methods */

/*
 * Comparer
 */

/** compares objective implications w.r.t. their maximum contribution */
static
SCIP_DECL_SORTPTRCOMP(objimplicsComp)
{
   SCIP_OBJIMPLICS* objimplics1;
   SCIP_OBJIMPLICS* objimplics2;

   objimplics1 = (SCIP_OBJIMPLICS*)elem1;
   objimplics2 = (SCIP_OBJIMPLICS*)elem2;

   if( objimplics1->maxobjchg > objimplics2->maxobjchg )
      return +1;

   if( objimplics1->maxobjchg < objimplics2->maxobjchg )
      return -1;

   return 0;
}

/** compare variables w.r.t.
 *  (i)   the absolute value the objective coefficient;
 *  (ii)  the locks which indicate most effect -- for the variables with a positive (negative) objective coefficient the
 *        down (up) lock is used since this lock indicates that tightened of the upper (lower) bound will triegger
 *        further domain propagations;
 *  (iii) the other locks;
 *  (iv)  variable problem index;
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

   /* first criteria is the absolute value of objective coefficient */
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
 * methods for SCIP_OBJIMPLICS data structure
 */

/** creates an objective implication data structure */
static
SCIP_RETCODE objimplicsCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OBJIMPLICS**     objimplics,         /**< pointer to objective implication data structure */
   SCIP_VAR**            objvars,            /**< objective contributor variables, or NULL */
   SCIP_Real             maxlbobjchg,        /**< maximum objective contributor if variables id fixed to zero */
   SCIP_Real             maxubobjchg,        /**< maximum objective contributor if variables id fixed to one */
   int                   nlbimpls,           /**< number of variables contributing to to lower bound fix */
   int                   nubimpls            /**< number of variables contributing to to upper bound fix */
   )

{
   assert(scip != NULL);
   assert(objimplics != NULL);
   assert(!SCIPisNegative(scip, maxlbobjchg));
   assert(!SCIPisNegative(scip, maxubobjchg));

   /* allocate block memory for the implication data structure */
   SCIP_CALL( SCIPallocBlockMemory(scip, objimplics) );

   if( nlbimpls + nubimpls == 0 )
   {
      assert(nlbimpls == 0);
      assert(nubimpls == 0);
      (*objimplics)->objvars = NULL;
   }
   else
   {
      int v;

      assert(objvars != NULL);

      /* copy the binary objective implication variables */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*objimplics)->objvars, objvars, nlbimpls + nubimpls) );

      /* capture the variables */
      for( v = 0; v < nlbimpls + nubimpls; ++v )
      {
         assert(SCIPvarIsBinary(objvars[v]));
         assert(!SCIPisZero(scip, SCIPvarGetObj(objvars[v])));
         SCIP_CALL( SCIPcaptureVar(scip, objvars[v]) );
      }
   }

   (*objimplics)->maxobjchg = MAX(maxlbobjchg, maxubobjchg);
   (*objimplics)->nlbimpls = nlbimpls;
   (*objimplics)->nubimpls = nubimpls;
   (*objimplics)->size = nlbimpls + nubimpls;

   return SCIP_OKAY;
}

/** frees an objective implication data structure */
static
SCIP_RETCODE objimplicsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OBJIMPLICS**     objimplics          /**< pointer to objective implication data structure */
   )
{
   int v;

   assert(scip != NULL);
   assert(objimplics != NULL);
   assert(*objimplics != NULL);

   /* release all variables */
   for( v = 0; v < (*objimplics)->nlbimpls + (*objimplics)->nubimpls; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*objimplics)->objvars[v]) );
   }

   /* free objective variable array */
   SCIPfreeBlockMemoryArrayNull(scip, &(*objimplics)->objvars, (*objimplics)->size);

   /* frees block memory for the implication data structure */
   SCIPfreeBlockMemory(scip, objimplics);

   return SCIP_OKAY;
}

/** remove the given variable at the given pos from the objective implication data structure */
static
SCIP_RETCODE objimplicsDelPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OBJIMPLICS*      objimplics,         /**< objective implication data structure */
   int                   pos                 /**< position */
   )
{
   assert(0 <= pos);
   assert(pos < objimplics->nlbimpls + objimplics->nubimpls);

   SCIPdebugMessage("variable <%s> ", SCIPvarGetName(objimplics->objvars[pos]));

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &objimplics->objvars[pos]) );

   /* copy last lower bound variable to that position */
   if( pos < objimplics->nlbimpls )
   {
      objimplics->nlbimpls--;
      assert(objimplics->nlbimpls >= 0);

      /* copy last lower bound variable to that position */
      objimplics->objvars[pos] = objimplics->objvars[objimplics->nlbimpls];

      /* copy last upper bound variable to open slot */
      objimplics->objvars[objimplics->nlbimpls] = objimplics->objvars[objimplics->nlbimpls + objimplics->nubimpls];

      SCIPdebugPrintf("remove lower bound implication\n");
   }
   else
   {
      objimplics->nubimpls--;
      assert(objimplics->nubimpls >= 0);

      /* copy last upper bound variable to that position */
      objimplics->objvars[pos] = objimplics->objvars[objimplics->nlbimpls + objimplics->nubimpls];

      SCIPdebugPrintf("remove upper bound implication\n");
   }

   return SCIP_OKAY;
}

/*
 * Local methods
 */


/** catch bound change events if the variable has a non-zero objective coefficient to check if the maximum activity of
 *  the objective function changed
 */
static
SCIP_RETCODE catchObjEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for global bound change events */
   SCIP_VAR*             var                 /**< variable for which the event should be dropped */
   )
{
   SCIP_Real objval;

   assert(propdata != NULL);
   assert(eventhdlr != NULL);

   objval = SCIPvarGetObj(var);

   if( !SCIPisZero(scip, objval) )
   {
      if( objval > 0.0 )
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_GUBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
      }
      else
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_GLBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
      }
   }

   return SCIP_OKAY;
}

/** drop variable event w.r.t. objective coefficient */
static
SCIP_RETCODE dropObjEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for global bound change events */
   SCIP_VAR*             var                 /**< variable for which the event should be dropped */
   )
{
   SCIP_Real objval;

   assert(propdata != NULL);
   assert(eventhdlr != NULL);

   objval = SCIPvarGetObj(var);

   /* drop bound change event */
   if( !SCIPisZero(scip, objval) )
   {
      if( objval > 0.0 )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_GUBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      }
      else
      {
         SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_GLBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      }
   }
   return SCIP_OKAY;
}

/** drop all variable events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_VAR* var;
   int k;

   assert(scip != NULL);
   assert(propdata != NULL);

   eventhdlr = propdata->eventhdlr;
   assert(eventhdlr != NULL);

   /* drop all events and release variables */
   for( k = 0; k < propdata->nminactvars; ++k )
   {
      var =  propdata->minactvars[k];
      assert(var != NULL);
      assert(SCIPvarIsBinary(var));

      /* drop bound relax event which is caught for all binary variables which are used for propagation the objective
       * function via the minimum activity of the objective function
       */
      SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDRELAXED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   /* release variables */
   for( k = 0; k < propdata->nmaxactvars; ++k )
   {
      var = propdata->maxactvars[k];
      assert(var != NULL);
      assert(SCIPvarIsBinary(var));

      /* drop events which are needed for evaluating the maximum activity of the objective function */
      SCIP_CALL( dropObjEvent(scip, propdata, eventhdlr, var) );

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   /* drop all events and release variables */
   for( k = 0; k < propdata->nobjintvars; ++k )
   {
      var = propdata->objintvars[k];
      assert(var != NULL);

      /* drop events which are needed for evaluating the maximum activity of the objective function */
      SCIP_CALL( dropObjEvent(scip, propdata, eventhdlr, var) );

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   return SCIP_OKAY;
}

/** reset propagatore data structure */
static
void propdataReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   propdata->minactvars = NULL;
   propdata->minactimpls = NULL;
   propdata->maxactvars = NULL;
   propdata->maxactchgs = NULL;
   propdata->objintvars = NULL;
   propdata->nminactvars = 0;
   propdata->nmaxactvars = 0;
   propdata->nobjintvars = 0;
   propdata->maxpseudoobjact = SCIP_INVALID;
   propdata->maxpseudoobjactinf = 0;
   propdata->lastvarnum = -1;
   propdata->glbpropagated = FALSE;
   propdata->cutoffbound = SCIPinfinity(scip);
   propdata->lastlowerbound = -SCIPinfinity(scip);
   propdata->glbpseudoobjval = -SCIPinfinity(scip);
   propdata->glbfirstnonfixed = 0;
   propdata->maxactfirstnonfixed = 0;
   propdata->firstnonfixed = 0;
   propdata->nnewvars = 0;
   propdata->catchvaradded = FALSE;
}

/** free propagator data */
static
SCIP_RETCODE propdataExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int v;

   if( propdata->addedvars != NULL )
      SCIPhashtableFree(&propdata->addedvars);

   /* drop events for the variables */
   SCIP_CALL( dropVarEvents(scip, propdata) );

   for( v = 0; v < propdata->nminactvars; ++v )
   {
      SCIP_CALL( objimplicsFree(scip, &(propdata->minactimpls[v])) );
   }

   /* free memory for non-zero objective variables */
   SCIPfreeMemoryArrayNull(scip, &propdata->minactvars);
   SCIPfreeMemoryArrayNull(scip, &propdata->minactimpls);
   SCIPfreeMemoryArrayNull(scip, &propdata->maxactvars);
   SCIPfreeMemoryArrayNull(scip, &propdata->maxactchgs);
   SCIPfreeMemoryArrayNull(scip, &propdata->objintvars);

   /* reset propagator data structure */
   propdataReset(scip, propdata);

   return SCIP_OKAY;
}

/** returns the objective change for the given binary variable */
static
SCIP_Real getVarObjchg(
   SCIP_VAR*             var,                /**< variable to get objective change for */
   SCIP_BOUNDTYPE        boundtype,          /**< bound type to consider */
   SCIP_BOUNDTYPE        bound               /**< fixing bound */
   )
{
   assert(SCIPvarIsBinary(var));
   assert((int)SCIP_BOUNDTYPE_LOWER == 0);
   assert((int)SCIP_BOUNDTYPE_UPPER == 1);

   /* collect contribution of variable itself */
   return (SCIP_Real)((int)bound - (int)(boundtype == SCIP_BOUNDTYPE_UPPER)) * SCIPvarGetObj(var);
}

/** returns the objective change provided by the implication variable by fixing it to the given bound
 *  w.r.t. minimum activity of the objective function; additionally it collects all contributors for that objective
 *  change;
 */
static
SCIP_Real collectMinactImplicVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_Bool             bound,              /**< bound to which the variable gets fixed due to an implication r */
   int*                  binobjidxs,         /**< sorted array which contains all variable problem indices of variable which have a non-zero objective coefficient */
   SCIP_Bool*            collectedvars,      /**< temporary buffer to mark collected variables */
   int                   nbinobjvars,        /**< number of binary variables with non-zero objective coefficient */
   SCIP_VAR**            contributors,       /**< array to store the contributors */
   int*                  ncontributors       /**< pointer to store number of contributor to the objective contribution */
   )
{
   SCIP_Real objval;
   SCIP_Bool diff;
   int pos;

   assert(var != NULL);

   /* ignore global fixed variables */
   if( SCIPvarGetLbGlobal(var) > 0.5 && SCIPvarGetUbGlobal(var) < 0.5 )
      return 0.0;

   objval = SCIPvarGetObj(var);

   /* ignore variables with zero objective coefficient */
   if( SCIPisZero(scip, objval) )
      return 0.0;

   assert(SCIPsortedvecFindInt(binobjidxs, SCIPvarGetProbindex(var), nbinobjvars, &pos));
   (void) SCIPsortedvecFindInt(binobjidxs, SCIPvarGetProbindex(var), nbinobjvars, &pos);

   /* check if the variables was already collected through other cliques */
   if( collectedvars[pos] )
      return 0.0;

   diff = (bound != (SCIP_Bool)SCIPvarGetBestBoundType(var));

   if( diff )
   {
      /* collect variable */
      assert(*ncontributors < nbinobjvars);
      contributors[*ncontributors] = var;
      (*ncontributors)++;

      /* mark variable to be collected */
      collectedvars[pos] = TRUE;

      /* return the absolute value of the objective coefficient as constriction */
      return REALABS(objval);
   }

   return 0.0;
}

/** returns the objective change provided by the implications of the given variable by fixing it to the given bound
 *  w.r.t. minimum activity of the objective function; additionally it collects all contributors for that objective
 *  change;
 *
 *  Let I(0) and I(1) be all implications of the given variable which follow by fixing it to given bound and evaluate to
 *  fixing the implication variable to zero (I(0)) or one (I(1)), respectively. The objective change provided by the
 *  implications are:
 *
 *  \f[
 *  \displaystyle
 *  sum_{x\in I(1)} (1 - \mbox{bestbound}(x)) \cdot \mbox{objval}(x) - sum_{x\in I(1)} \mbox{bestbound}(x) \cdot \mbox{objval}(x)
 *  =
 *  sum_{x\in I(0) \cup I(1)} (\mbox{impliedbound}(x) - \mbox{bestbound}(x)) \cdot \mbox{objval}(x)
 *  \f]
 */
static
SCIP_RETCODE collectMinactImplicVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   int*                  binobjidxs,         /**< sorted array which contains all variable problem indices of variable which have a non-zero objective coefficient */
   SCIP_Bool*            collectedvars,      /**< temporary buffer to mark collected variables */
   int                   nbinobjvars,        /**< number of binary variables with non-zero objective coefficient */
   SCIP_VAR**            contributors,       /**< array to store the contributors */
   int*                  ncontributors,      /**< pointer to store number of contributor to the objective contribution */
   SCIP_Real*            objchg              /**< pointer to store the objective change */
   )
{
   SCIP_CLIQUE** cliques;
   SCIP_CLIQUE* clique;
   SCIP_VAR** vars;
   SCIP_VAR* implvar;
   SCIP_Real* bounds;
   SCIP_Bool* values;
   SCIP_Bool varfixing;
   int nbinvars;
   int ncliques;
   int c;
   int v;

   assert(SCIPvarIsBinary(var));
   assert(SCIPvarGetLbGlobal(var) < 0.5);
   assert(SCIPvarGetUbGlobal(var) > 0.5);
   assert(bound == SCIP_BOUNDTYPE_LOWER || bound == SCIP_BOUNDTYPE_UPPER);
   assert(objchg != NULL);
   assert(contributors != NULL);
   assert(ncontributors != NULL);
   assert(*ncontributors == 0);

   assert((SCIP_Bool)SCIP_BOUNDTYPE_LOWER == FALSE);
   assert((SCIP_Bool)SCIP_BOUNDTYPE_UPPER == TRUE);
   varfixing = (SCIP_Bool)bound;

   cliques = SCIPvarGetCliques(var, varfixing);
   ncliques = SCIPvarGetNCliques(var, varfixing);

#ifndef NDEBUG
   /* check that the marker array is reset */
   for( c = 0; c < nbinobjvars; ++c )
      assert(collectedvars[c] == FALSE);
#endif

   /* collect all implication which are given via cliques */
   for( c = 0; c < ncliques; ++c )
   {
      clique = cliques[c];
      assert(clique != NULL);

      nbinvars = SCIPcliqueGetNVars(clique);
      vars = SCIPcliqueGetVars(clique);
      values = SCIPcliqueGetValues(clique);

      for( v = 0; v < nbinvars; ++v )
      {
         implvar = vars[v];
         assert(implvar != NULL);

         if( implvar == var )
            continue;

         (*objchg) += collectMinactImplicVar(scip, implvar, !values[v], binobjidxs, collectedvars, nbinobjvars, contributors, ncontributors);
      }
   }

   /* collect implications */
   vars = SCIPvarGetImplVars(var, varfixing);
   nbinvars = SCIPvarGetNBinImpls(var, varfixing);
   bounds = SCIPvarGetImplBounds(var, varfixing);

   /* loop over all implications */
   for( v = 0; v < nbinvars; ++v )
   {
      assert(vars[v] != NULL);
      (*objchg) += collectMinactImplicVar(scip, vars[v], bounds[v] > 0.5, binobjidxs, collectedvars, nbinobjvars, contributors, ncontributors);
   }

   return SCIP_OKAY;
}

/** returns the objective change provided by the implications of the given variable by fixing it to the given bound
 *  w.r.t. minimum activity of the objective function
 *
 *  Let I(0) and I(1) be all implications of the given variable which follow by fixing it to given bound and evaluate to
 *  fixing the implication variable to zero (I(0)) or one (I(1)), respectively. The objective change provided by the
 *  implications are:
 *
 *  \f[
 *  \displaystyle
 *  sum_{x\in I(1)} (1 - \mbox{bestbound}(x)) \cdot \mbox{objval}(x) - sum_{x\in I(1)} \mbox{bestbound}(x) \cdot \mbox{objval}(x)
 *  =
 *  sum_{x\in I(0) \cup I(1)} (\mbox{impliedbound}(x) - \mbox{bestbound}(x)) \cdot \mbox{objval}(x)
 *  \f]
 *
 *  This can be done w.r.t. global variable bounds (local == FALSE), w.r.t. local variable bounds (local == TRUE &&
 *  bdchgidx == NULL), and w.r.t. given time stamp (local == TRUE && bdchgidx != NULL)
 */
static
SCIP_RETCODE getMinactImplicObjchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_OBJIMPLICS*      objimplics,         /**< objective implication data for the given variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, or NULL */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   SCIP_Bool             local,              /**< propagate local bounds, otherwise global bounds */
   SCIP_Real*            objchg              /**< pointer to store the objective change */
   )
{
   SCIP_VAR* implvar;
   SCIP_Bool lb;
   SCIP_Bool ub;
   int nbinvars;
   int v;

   assert(SCIPvarIsBinary(var));
   assert(!local || SCIPvarGetLbAtIndex(var, bdchgidx, FALSE) < 0.5);
   assert(!local || SCIPvarGetUbAtIndex(var, bdchgidx, FALSE) > 0.5);
   assert(SCIPvarGetLbGlobal(var) < 0.5);
   assert(SCIPvarGetUbGlobal(var) > 0.5);
   assert(bound == SCIP_BOUNDTYPE_LOWER || bound == SCIP_BOUNDTYPE_UPPER);

   if( bound == SCIP_BOUNDTYPE_LOWER )
   {
      v = 0;
      nbinvars = objimplics->nlbimpls;
   }
   else
   {
      assert(bound == SCIP_BOUNDTYPE_UPPER);
      v = objimplics->nlbimpls;
      nbinvars = objimplics->nlbimpls + objimplics->nubimpls;
   }

   /* loop over all implications */
   while( v < nbinvars )
   {
      implvar = objimplics->objvars[v];
      assert(implvar != NULL);
      assert(!SCIPisZero(scip, SCIPvarGetObj(implvar)));

      if( local )
      {
         lb = SCIPvarGetLbAtIndex(implvar, bdchgidx, FALSE) > 0.5;
         ub = SCIPvarGetUbAtIndex(implvar, bdchgidx, FALSE) > 0.5;

         /* check if variable is fixed */
         if( lb == TRUE || ub == FALSE )
         {
            v++;
            continue;
         }
      }
      else
      {
         lb = SCIPvarGetLbGlobal(implvar) > 0.5;
         ub = SCIPvarGetUbGlobal(implvar) > 0.5;

         /* check if variable is global fixed; if so remove it from the objective implication data structure and
          * continue with the next candidate
          */
         if( lb == TRUE || ub == FALSE )
         {
            SCIP_CALL( objimplicsDelPos(scip, objimplics, v) );
            nbinvars--;
            continue;
         }
      }

      /* add objective change */
      (*objchg) += REALABS(SCIPvarGetObj(implvar));
      v++;
   }

   return SCIP_OKAY;
}

/** computes for the given binary variable the objective contribution by fixing it to given bound w.r.t. minimum
 *  activity of the objective function; additionally it collects all contributors for that objective change;
 */
static
SCIP_RETCODE collectMinactObjchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   int*                  binobjidxs,         /**< sorted array which contains all variable problem indices of variable which have a non-zero objective coefficient */
   SCIP_Bool*            collectedvars,      /**< temporary buffer to mark collected variables */
   int                   nbinobjvars,        /**< number of binary variables with non-zero objective coefficient */
   SCIP_VAR**            contributors,       /**< array to store the contriboters */
   int*                  ncontributors,      /**< pointer to store number of contributor to the objective contribution */
   SCIP_Real*            objchg              /**< pointer to store the objective change */
   )
{
   assert(SCIPvarIsBinary(var));
   assert(contributors != NULL);
   assert(ncontributors != NULL);

   /* collects the contribution of the variable itself w.r.t. the best bound */
   (*objchg) = getVarObjchg(var, SCIPvarGetBestBoundType(var), bound);

   (*ncontributors) = 0;

   /* add the objective contribution from the implication variable */
   SCIP_CALL( collectMinactImplicVars(scip, var, bound, binobjidxs, collectedvars, nbinobjvars, contributors, ncontributors, objchg) );

   return SCIP_OKAY;
}

/** computes for the given binary variable the objective contribution by fixing it to given bound w.r.t. minimum
 *  activity of the objective function; this can be done w.r.t. global variable bounds (local == FALSE), w.r.t. local
 *  variable bounds (local == TRUE && bdchgidx == NULL), and w.r.t. given time stamp (local == TRUE && bdchgidx != NULL)
 */
static
SCIP_RETCODE getMinactObjchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_OBJIMPLICS*      objimplics,         /**< objective implication data for the given variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, or NULL */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   SCIP_Bool             local,              /**< propagate local bounds, otherwise global bounds */
   SCIP_Real*            objchg              /**< pointer to store the objective change */
   )
{
   assert(SCIPvarIsBinary(var));

   /* collects the contribution of the variable itself w.r.t. the best bound */
   (*objchg) = getVarObjchg(var, SCIPvarGetBestBoundType(var), bound);

   /* add the objective contribution from the implication variable */
   SCIP_CALL( getMinactImplicObjchg(scip, var, objimplics, bdchgidx, bound, local, objchg) );

   return SCIP_OKAY;
}

/** returns the global (that means w.r.t. global bounds of the variables) objective change provided by the implication of
 *  the given variable by fixing it to the given bound w.r.t. maximum activity of the objective function
 *
 *  Let I(0) and I(1) be all implications of the given variable which follow by fixing it to given bound and evaluate to
 *  fixing the implication variable to zero (I(0)) or one (I(1)), respectively. The objective change provided by the
 *  implications are:
 *
 *  \f[
 *  \displaystyle
 *  sum_{x\in I(1)} (1 - \mbox{worstbound}(x)) \cdot \mbox{objval}(x) - sum_{x\in I(1)} \mbox{worst}(x) \cdot \mbox{objval}(x)
 *  =
 *  sum_{x\in I(0) \cup I(1)} (\mbox{impliedbound}(x) - \mbox{worstbound}(x)) \cdot \mbox{objval}(x)
 *  \f]
 */
static
SCIP_Real getMaxactImplicObjchg(
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_BOUNDTYPE        bound               /**< bound to check for */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* implvar;
   SCIP_Real* bounds;
   SCIP_Real objchg;
   SCIP_Bool varfixing;
   int nbinvars;
   int v;

   assert(SCIPvarIsBinary(var));
   assert(SCIPvarGetLbGlobal(var) < 0.5);
   assert(SCIPvarGetUbGlobal(var) > 0.5);
   assert(bound == SCIP_BOUNDTYPE_LOWER || bound == SCIP_BOUNDTYPE_UPPER);

   varfixing = (SCIP_Bool)bound;
   assert((SCIP_Bool)SCIP_BOUNDTYPE_LOWER == FALSE);
   assert((SCIP_Bool)SCIP_BOUNDTYPE_UPPER == TRUE);

   vars = SCIPvarGetImplVars(var, varfixing);
   nbinvars = SCIPvarGetNBinImpls(var, varfixing);
   bounds = SCIPvarGetImplBounds(var, varfixing);

   objchg = 0.0;

   /* loop over all implications */
   for( v = 0; v < nbinvars; ++v )
   {
      implvar = vars[v];
      assert(implvar != NULL);

      /* ignore globally fixed variables */
      if( SCIPvarGetLbGlobal(implvar) < 0.5 && SCIPvarGetUbGlobal(implvar) > 0.5 )
         objchg += (bounds[v] - SCIPvarGetWorstBoundGlobal(implvar)) * SCIPvarGetObj(implvar);
   }

   return objchg;
}

/** computes for the given binary variable the gloabl (that means w.r.t. global bounds of the variables) objective
 *  contribution by fixing it to given bound w.r.t. maximum activity of the objective function
 */
static
SCIP_Real getMaxactObjchg(
   SCIP_VAR*             var,                /**< variable to computes the objective contribution */
   SCIP_BOUNDTYPE        bound,              /**< bound to check for */
   SCIP_Bool             useimplics          /**< should implications be used */
   )
{
   SCIP_Real objchg;

   assert(SCIPvarIsBinary(var));

   /* collects the contribution of the variable itself w.r.t. the worst bound */
   objchg = getVarObjchg(var, SCIPvarGetWorstBoundType(var), bound);

   /* check if the implications should be used to increase the objective contribution for given variable */
   if( useimplics )
   {
      /* add the objective contribution from the implication variable */
      objchg += getMaxactImplicObjchg(var, bound);
   }

   return objchg;
}

/** check if the given variable should be collected for the minimum activity propagation */
static
SCIP_RETCODE collectMinactVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   SCIP_OBJIMPLICS**     objimplics,         /**< pointer to store the objective implication data structure w.r.t. minimum activity */
   SCIP_Bool             useimplics,         /**< should implications be used */
   int*                  binobjidxs,         /**< sorted array which contains all variable problem indices of variable which have a non-zero objective coefficient */
   SCIP_Bool*            collectedvars,      /**< temporary buffer to mark collected variables */
   int                   nbinobjvars,        /**< number of binary variables with non-zero objective coefficient */
   SCIP_VAR**            contributors,       /**< temporary buffer to use for collecting contributors */
   SCIP_Bool*            collect             /**< pointer to store if the variable should be stored */
   )
{
   SCIP_Real lbobjchg;
   SCIP_Real ubobjchg;
   SCIP_Real objval;
   int nlbimplics;
   int nubimplics;
   int nlbcliques;
   int nubcliques;

   assert(objimplics != NULL);

   objval = SCIPvarGetObj(var);
   (*objimplics) = NULL;

   if( SCIPisZero(scip, objval) )
      (*collect) = FALSE;
   else
      (*collect) = TRUE;

   nlbimplics = SCIPvarGetNBinImpls(var, FALSE);
   nubimplics = SCIPvarGetNBinImpls(var, TRUE);

   nlbcliques = SCIPvarGetNCliques(var,  FALSE);
   nubcliques = SCIPvarGetNCliques(var,  TRUE);

   /* check if implications should be used and if implications are existing */
   if( useimplics && nlbimplics + nubimplics + nlbcliques + nubcliques > 0 )
   {
      int nlbcontributors;
      int nubcontributors;
      int pos;
      int c;

      assert((SCIP_Bool)SCIP_BOUNDTYPE_LOWER == FALSE);
      assert((SCIP_Bool)SCIP_BOUNDTYPE_UPPER == TRUE);

      /* get contribution of variable by fixing it to its lower bound w.r.t. minimum activity of the objective function */
      SCIP_CALL( collectMinactObjchg(scip, var, SCIP_BOUNDTYPE_LOWER, binobjidxs, collectedvars, nbinobjvars,
            contributors, &nlbcontributors, &lbobjchg) );
      assert(!SCIPisNegative(scip, lbobjchg));

      SCIPdebugMessage("variable <%s> fixed to bound=%d implies %d(%d)\n", SCIPvarGetName(var),
         SCIP_BOUNDTYPE_LOWER, SCIPvarGetNBinImpls(var, (SCIP_Bool)SCIP_BOUNDTYPE_LOWER), nlbcontributors);

      /* reset variables array which marks variables which are collected */
      for( c = 0; c < nlbcontributors; ++c )
      {
         assert(SCIPsortedvecFindInt(binobjidxs, SCIPvarGetProbindex(contributors[c]), nbinobjvars, &pos));
         (void) SCIPsortedvecFindInt(binobjidxs, SCIPvarGetProbindex(contributors[c]), nbinobjvars, &pos);
         collectedvars[pos] = FALSE;
      }

      /* ignore implications if the variable has a zero objective coefficient and implications only one variable, since
       * this covered by that implied variable
       */
      if( !(*collect) && nlbcontributors == 1 )
      {
         assert(SCIPisZero(scip, objval));
         nlbcontributors = 0;
      }

      /* get contribution of variable by fixing it to its upper bound w.r.t. minimum activity of the objective function */
      SCIP_CALL( collectMinactObjchg(scip, var, SCIP_BOUNDTYPE_UPPER, binobjidxs, collectedvars, nbinobjvars,
            &contributors[nlbcontributors], &nubcontributors, &ubobjchg) );
      assert(!SCIPisNegative(scip, ubobjchg));

      SCIPdebugMessage("variable <%s> fixed to bound=%d implies %d(%d)\n", SCIPvarGetName(var),
         SCIP_BOUNDTYPE_UPPER, SCIPvarGetNBinImpls(var, (SCIP_Bool)SCIP_BOUNDTYPE_UPPER), nubcontributors);

      /* reset variables array which marks variables which are collected */
      for( c = 0; c < nubcontributors; ++c )
      {
         assert(SCIPsortedvecFindInt(binobjidxs, SCIPvarGetProbindex(contributors[nlbcontributors + c]), nbinobjvars, &pos));
         (void) SCIPsortedvecFindInt(binobjidxs, SCIPvarGetProbindex(contributors[nlbcontributors + c]), nbinobjvars, &pos);
         collectedvars[pos] = FALSE;
      }

      /* ignore implications if the variable has a zero objective coefficient and implications only one variable, since
       * this covered by that implied variable
       */
      if( !(*collect) && nubcontributors == 1 )
      {
         assert(SCIPisZero(scip, objval));
         nubcontributors = 0;
      }

      if( (*collect) || nlbcontributors > 1 || nubcontributors > 1 )
      {
         SCIP_CALL( objimplicsCreate(scip, objimplics, contributors, lbobjchg, ubobjchg, nlbcontributors, nubcontributors) );
         (*collect) = TRUE;
      }
   }
   else if( (*collect) )
   {
      lbobjchg = getVarObjchg(var, SCIPvarGetBestBoundType(var), SCIP_BOUNDTYPE_LOWER);
      ubobjchg = getVarObjchg(var, SCIPvarGetBestBoundType(var), SCIP_BOUNDTYPE_UPPER);
      assert(!SCIPisZero(scip, lbobjchg) || !SCIPisZero(scip, ubobjchg));
      assert(!SCIPisNegative(scip, lbobjchg));
      assert(!SCIPisNegative(scip, ubobjchg));

      SCIP_CALL( objimplicsCreate(scip, objimplics, NULL, lbobjchg, ubobjchg, 0, 0) );
   }

   return SCIP_OKAY;
}

/** check if the given variable should be collected for the maximum activity propagation */
static
SCIP_Bool collectMaxactVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   SCIP_Bool             useimplics,         /**< should implications be used */
   SCIP_Real*            objchg              /**< pointer to store the objective change w.r.t. maximum activity */
   )
{
   SCIP_Real lbobjchg;
   SCIP_Real ubobjchg;

   /* get contribution of variable by fixing it to its lower bound w.r.t. maximum activity of the objective function */
   lbobjchg = getMaxactObjchg(var, SCIP_BOUNDTYPE_LOWER, useimplics);
   assert(!SCIPisPositive(scip, lbobjchg));

   /* get contribution of variable by fixing it to its upper bound w.r.t. maximum activity of the objective function */
   ubobjchg = getMaxactObjchg(var, SCIP_BOUNDTYPE_UPPER, useimplics);
   assert(!SCIPisPositive(scip, ubobjchg));

   (*objchg) = MIN(lbobjchg, ubobjchg);

   /* only consider variables with non-zero objective contribution */
   if( SCIPisZero(scip, (*objchg)) )
      return FALSE;

   return TRUE;
}

/** initializate the propagator */
static
SCIP_RETCODE propdataInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   int* binobjvars;
   int sizebinobjvars;
   int nvars;
   int nbinvars;
   int nintvars;
   int nminactvars;
   int nmaxactvars;
   int nobjintvars;
   int nobjcontvars;
   int nobjvars;
   int nbinobjvars;
   int v;

   assert(scip != NULL);
   assert(propdata != NULL);

   /* get problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nintvars = nvars - SCIPgetNContVars(scip);

   nbinvars = 0;
   nobjvars = 0;
   nbinobjvars = 0;
   sizebinobjvars = 100;

   SCIP_CALL( SCIPallocBufferArray(scip, &binobjvars, sizebinobjvars) );

   /* count and collect variable problem indices of variables with non-zero objective coefficient */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);

      if( !SCIPisZero(scip, SCIPvarGetObj(var)) )
      {
         nobjvars++;

         if( SCIPvarIsBinary(var) )
         {
            if( nbinobjvars == sizebinobjvars )
            {
               sizebinobjvars = SCIPcalcMemGrowSize(scip, nbinobjvars + 1);
               SCIP_CALL( SCIPreallocBufferArray(scip, &binobjvars, sizebinobjvars) );
            }

            binobjvars[nbinobjvars] = SCIPvarGetProbindex(var);
            nbinobjvars++;
         }
      }

      if( SCIPvarIsBinary(var) )
         nbinvars++;
   }

   nminactvars = 0;
   nmaxactvars = 0;
   nobjintvars = 0;
   nobjcontvars = 0;

   if( nobjvars > 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;
      SCIP_OBJIMPLICS* objimplics;
      SCIP_VAR** contributors;
      SCIP_Bool* collectedvars;
      SCIP_Bool collect;
      SCIP_Bool useimplics;
      SCIP_Real objval;
      SCIP_Real objchg;

      eventhdlr = propdata->eventhdlr;
      assert(eventhdlr != NULL);

      useimplics = (propdata->propuseimplics && nbinobjvars < propdata->maximplvars);

      /* allocate memory for all arrays */
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->minactvars, nbinvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->minactimpls, nbinvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->maxactvars, nbinvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->maxactchgs, nbinvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->objintvars, nobjvars - nbinobjvars) );

      if( useimplics )
      {
         /* create temporary buffer */
         SCIP_CALL( SCIPallocBufferArray(scip, &contributors, nbinobjvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &collectedvars, nbinobjvars) );
         BMSclearMemoryArray(collectedvars, nbinobjvars);
      }
      else
      {
         contributors = NULL;
         collectedvars = NULL;
      }

      /* collect the variables with non-zero objective contribution and catch global bound tighten events that decrease the
       * maximal pseudo objective activity
       */
      for( v = 0; v < nvars && (nobjintvars == 0 || nobjintvars < nobjvars - nbinobjvars); ++v )
      {
         var = vars[v];
         assert(var != NULL);

         objval = SCIPvarGetObj(var);

         if( SCIPvarIsBinary(var) )
         {
            /* ignore variables which are globally fixed */
            if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5 )
            {
#ifndef NDEBUG
               /* check that the binary implications are applied for binary variables which are globally fixed */
               checkImplicsApplied(scip, var);
#endif
               continue;
            }

            /* check if the variable should be collected for the minimum activity propagation */
            SCIP_CALL( collectMinactVar(scip, var, &objimplics, useimplics, binobjvars, collectedvars, nbinobjvars, contributors, &collect) );

            if( collect )
            {
               assert(nminactvars < nbinvars);
               assert(objimplics != NULL);
               assert(objimplics->nlbimpls + objimplics->nubimpls <= nbinobjvars);

               /* collect the binary variable with non-zero objective contribution */
               propdata->minactvars[nminactvars] = var;
               propdata->minactimpls[nminactvars] = objimplics;
               nminactvars++;

               /* catch bound relax event for the binary variable to handel the firstnonfixed index correctly */
               SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDRELAXED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );

               SCIPdebugMessage("variable <%s>[obj: <%g>] implicit objective change %g\n",
                  SCIPvarGetName(var), objval, objimplics->maxobjchg);

               /* captures the variable */
               SCIP_CALL( SCIPcaptureVar(scip, var) ) ;
            }

            /* check if the variable should be collected for the maximum activity propagation */
            collect = collectMaxactVar(scip, var, useimplics, &objchg);

            if( collect )
            {
               assert(nmaxactvars < nbinvars);

               /* collect the binary variable with non-zero objective contribution */
               propdata->maxactvars[nmaxactvars] = var;
               propdata->maxactchgs[nmaxactvars] = -objchg;
               nmaxactvars++;

               /* catch bound change events if the variable has a non-zero objective coefficient to check if the maximum
                * activity of the objective function changed
                */
               SCIP_CALL( catchObjEvent(scip, propdata, eventhdlr, var) );

               /* captures the variable */
               SCIP_CALL( SCIPcaptureVar(scip, var) ) ;
            }
         }
         else
         {
            /* only consider none binary variables with a non-zero objective */
            if( SCIPisZero(scip, objval) )
               continue;

            assert(nobjintvars < nobjvars - nbinobjvars);

            propdata->objintvars[nobjintvars] = var;
            nobjintvars++;

            if( v >= nintvars )
               nobjcontvars++;

            /* catch bound change events if the variable has a non-zero objective coefficient to check if the maximum
             * activity of the objective function changed
             */
            SCIP_CALL( catchObjEvent(scip, propdata, eventhdlr, var) );

            /* captures the variable */
            SCIP_CALL( SCIPcaptureVar(scip, var) );
         }
      }

      if( useimplics )
      {
         SCIPfreeBufferArray(scip, &collectedvars);
         SCIPfreeBufferArray(scip, &contributors);
      }

      if( nminactvars == 0 )
      {
         SCIPfreeMemoryArray(scip, &propdata->minactvars);
         SCIPfreeMemoryArray(scip, &propdata->minactimpls);
         propdata->minactvars = NULL;
         propdata->minactimpls = NULL;
      }
      else
      {
         /* sort binary variables with respect to the absolute value of their maximal potential objective contribution for
          * the minimum activity of the objective function
          */
         SCIPsortDownPtrPtr((void**)propdata->minactimpls, (void**)propdata->minactvars, objimplicsComp, nminactvars);

         SCIPdebugMessage("%d binary variables with non-zero objective contribution w.r.t. the minimum activity of the objective function\n", nminactvars);
      }

      if( nmaxactvars == 0 )
      {
         SCIPfreeMemoryArray(scip, &propdata->maxactvars);
         SCIPfreeMemoryArray(scip, &propdata->maxactchgs);
         propdata->maxactvars = NULL;
         propdata->maxactchgs = NULL;
      }
      else
      {
         /* sort binary variables with respect to the absolute value of their maximal potential objective contribution for
          * the maximum activity of the objective function
          */
         SCIPsortDownRealPtr(propdata->maxactchgs, (void**)propdata->maxactvars, nmaxactvars);

         SCIPdebugMessage("%d binary variables with non-zero objective contribution w.r.t. the maximum activity of the objective function\n", nmaxactvars);
      }

      if( nobjintvars == 0 )
      {
         SCIPfreeMemoryArray(scip, &propdata->objintvars);
         propdata->objintvars = NULL;
      }
      else
      {
         /* sort integer variables with respect to the absolute value of their objective coefficient */
         SCIPsortDownPtr((void**)propdata->objintvars, varCompObj, nobjintvars - nobjcontvars);

         /* sort continuous variables with respect to the absolute value of their objective coefficient */
         SCIPsortDownPtr((void**)(&propdata->objintvars[nobjintvars - nobjcontvars]), varCompObj, nobjcontvars);

         SCIPdebugMessage("%d integer variables and %d continuous variables with non-zero objective contribution\n",
            nobjintvars - nobjcontvars, nobjcontvars);
      }
   }

   SCIPfreeBufferArray(scip, &binobjvars);

   propdata->nminactvars = nminactvars;
   propdata->nmaxactvars = nmaxactvars;
   propdata->nobjintvars = nobjintvars;
   propdata->maxpseudoobjact = SCIP_INVALID;
   propdata->maxpseudoobjactinf = 0;
   propdata->lastvarnum = -1;
   propdata->glbfirstnonfixed = 0;
   propdata->maxactfirstnonfixed = 0;
   propdata->firstnonfixed = 0;
   propdata->nnewvars = 0;

   /* due to scaling after presolving we need to update the global pseudoactivity and the cutoffbound */
   propdata->glbpropagated = FALSE;
   propdata->glbpseudoobjval = SCIPgetGlobalPseudoObjval(scip);
   propdata->cutoffbound = SCIPgetCutoffbound(scip);
   assert(SCIPisFeasEQ(scip, propdata->glbpseudoobjval, SCIPgetPseudoObjval(scip)));

   /* create hash table which is used for resolving bound changes */
   if( nminactvars > 0 )
   {
      SCIP_CALL( SCIPhashtableCreate(&propdata->addedvars, SCIPblkmem(scip), SCIPcalcHashtableSize(5*nvars),
            SCIPvarGetHashkey, SCIPvarIsHashkeyEq, SCIPvarGetHashkeyVal, NULL) );
   }
   else
      propdata->addedvars = NULL;


   return SCIP_OKAY;
}

/** adds for the given none binary variable a conflict bound depending on its objective contribution */
static
SCIP_RETCODE addConflictBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check for objective contribution */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_Real*            reqpseudoobjval     /**< pointer to store the remaining minimum activity which has to be proven */
   )
{
   SCIP_Real objval;

   objval = SCIPvarGetObj(var);
   assert(!SCIPisZero(scip, objval));

   if( objval > 0.0 )
   {
      SCIP_Real loclb;
      SCIP_Real glblb;

      glblb = SCIPvarGetLbGlobal(var);
      loclb = SCIPvarGetLbAtIndex(var, bdchgidx, FALSE);
      assert(SCIPisFeasGE(scip, loclb, glblb));

      /* check if the local lower bound (at time stamp bdchgidx) is larger than the global lower bound */
      if( SCIPisGT(scip, loclb, glblb) )
      {
         SCIPdebugMessage("  add bound change <%s>[%g] >= <%g>\n", SCIPvarGetName(var), objval, loclb);
         SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );

         assert(SCIPisPositive(scip, (loclb - glblb) * objval));
         (*reqpseudoobjval) -= (loclb - glblb) * objval;
      }
   }
   else
   {
      SCIP_Real locub;
      SCIP_Real glbub;

      glbub = SCIPvarGetUbGlobal(var);
      locub = SCIPvarGetUbAtIndex(var, bdchgidx, FALSE);
      assert(SCIPisFeasLE(scip, locub, glbub));

      /* check if the local upper bound (at time stamp bdchgidx) is smaller than the global upper bound */
      if( SCIPisLT(scip, locub, glbub) )
      {
         SCIPdebugMessage("  add bound change <%s>[%g] <= <%g>\n", SCIPvarGetName(var), objval, locub);
         SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );

         assert(SCIPisPositive(scip, (locub - glbub) * objval));
         (*reqpseudoobjval) -= (locub - glbub) * objval;
      }
   }

   return SCIP_OKAY;
}

/** check for the given implication variables of they also contribute to the required minimum activity */
static
SCIP_RETCODE getConflictImplics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable to check for objective contribution */
   int                   start,              /**< start index */
   int                   end,                /**< end index */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_HASHTABLE*       addedvars,          /**< hash table containing variables which are already add directly or implicitly due to implications */
   SCIP_Real*            reqpseudoobjval,    /**< pointer to store the remaining minimum activity which has to be proven */
   SCIP_Bool*            foundimplics        /**< pointer to store if an implication is found */
   )
{
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   int v;

   assert(foundimplics != NULL);
   assert(*foundimplics == FALSE);

   for( v = start; v < end; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarIsBinary(var));

      lb = SCIPvarGetLbAtIndex(var, bdchgidx, FALSE);
      ub = SCIPvarGetUbAtIndex(var, bdchgidx, FALSE);

      if( lb < 0.5 && ub > 0.5 && !SCIPhashtableExists(addedvars, (void*)var) )
      {
         (*reqpseudoobjval) -= REALABS(SCIPvarGetObj(var));
         SCIPdebugMessage("  implicated variables <%s>[%g] bdchgidx [%g,%g] -> remaining <%g>\n", SCIPvarGetName(var), SCIPvarGetObj(var), lb, ub, *reqpseudoobjval);

         SCIP_CALL( SCIPhashtableInsert(addedvars, (void*)var) );
         assert(SCIPhashtableExists(addedvars, (void*)var));
         (*foundimplics) = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** adds for the given binary variable a conflict bound depending on its objective contribution */
static
SCIP_RETCODE addConflictBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check for objective contribution */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_OBJIMPLICS*      objimplics,         /**< objective implication data for the given variable */
   SCIP_HASHTABLE*       addedvars,          /**< hash table containing variables which are already add directly or implicitly due to implications */
   SCIP_Bool             respropuseimplics,  /**< should implications be used */
   SCIP_Real*            reqpseudoobjval     /**< pointer to store the remaining minimum activity which has to be proven */
   )
{
   SCIP_Real objval;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool foundimplics;

   assert(SCIPvarIsBinary(var));

   if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5 )
      return SCIP_OKAY;

   lb = SCIPvarGetLbAtIndex(var, bdchgidx, FALSE);
   ub = SCIPvarGetUbAtIndex(var, bdchgidx, FALSE);

   objval = SCIPvarGetObj(var);
   foundimplics = FALSE;

   /* only consider variables which are fixed */
   if( lb > 0.5 )
   {
      if( respropuseimplics )
      {
         SCIP_CALL( getConflictImplics(scip, objimplics->objvars, objimplics->nlbimpls, objimplics->nlbimpls + objimplics->nubimpls,
               bdchgidx, addedvars, reqpseudoobjval, &foundimplics) );
      }

      /* check if the binary variable has a positive contribution (positive objective coefficient since it is fixed to
       * one) or is needed due a positive contribution of an implied variable
       */
      if( foundimplics || SCIPisPositive(scip, objval) )
      {
         SCIPdebugMessage("  add bound change <%s>[%g] >= <%g> bdchgidx [%g,%g]\n", SCIPvarGetName(var), objval, lb, lb, ub);
         SCIP_CALL( SCIPaddConflictLb(scip, var, NULL) );

         (*reqpseudoobjval) -= MAX(0.0, objval);

         if( addedvars != NULL )
         {
            assert(!SCIPhashtableExists(addedvars, (void*)var));
            SCIP_CALL( SCIPhashtableInsert(addedvars, (void*)var) );
         }
      }
   }
   else if( ub < 0.5 )
   {
      if( respropuseimplics )
      {
         SCIP_CALL( getConflictImplics(scip, objimplics->objvars, 0, objimplics->nlbimpls,
               bdchgidx, addedvars, reqpseudoobjval, &foundimplics) );
      }

      /* check if the binary variable has a positive contribution (negative objective coefficient since it is fixed to
       * zero) or is needed due a positive contribution of an implied variable
       */
      if( foundimplics || SCIPisNegative(scip, objval) )
      {
         SCIPdebugMessage("  add bound change <%s>[%g] <= <%g> bdchgidx=[%g,%g]\n", SCIPvarGetName(var), objval, ub, lb, ub);
         SCIP_CALL( SCIPaddConflictUb(scip, var, NULL) );

         (*reqpseudoobjval) +=  MIN(0.0, objval);

         if( addedvars != NULL )
         {
            assert(!SCIPhashtableExists(addedvars, (void*)var));
            SCIP_CALL( SCIPhashtableInsert(addedvars, (void*)var) );
         }
      }
   }

   return SCIP_OKAY;
}


/** resolves a propagation by supplying the variables whose bound changes increased the pseudo objective value above the
 *  cutoff bound
 */
static
SCIP_RETCODE adjustCutoffbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable that was deduced */
   int                   inferinfo,          /**< inference information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_HASHTABLE*       addedvars,          /**< hash table which contains variable which are already added or implict given as reason for the resolve, or NULL */
   SCIP_Real*            cutoffbound         /**< pointer to store the adjusted cutoff bound */
   )
{
   if( inferinfo != -1 )
   {
      SCIP_OBJIMPLICS* objimplics;
      SCIP_Bool foundimplics;
      int start;
      int end;

      assert(var != NULL);
      assert(SCIPvarIsBinary(var));
      assert(bdchgidx != NULL);
      assert(SCIPisEQ(scip, SCIPvarGetLbAtIndex(var, bdchgidx, TRUE), SCIPvarGetUbAtIndex(var, bdchgidx, TRUE)));
      assert(inferinfo >= 0);
      assert(inferinfo < propdata->nminactvars);
      assert((SCIP_Bool)SCIP_BOUNDTYPE_LOWER == FALSE);
      assert((SCIP_Bool)SCIP_BOUNDTYPE_UPPER == TRUE);

      objimplics = propdata->minactimpls[inferinfo];
      assert(objimplics != NULL);

      /* get the objective contribution if we would fix the binray inference variable to its other bound */
      (*cutoffbound) -= getVarObjchg(var, SCIPvarGetBestBoundType(var), boundtype);
      foundimplics = FALSE;

      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         start = 0;
         end = objimplics->nlbimpls;
      }
      else
      {
         start = objimplics->nlbimpls;
         end = objimplics->nlbimpls + objimplics->nubimpls;
      }

      SCIP_CALL( getConflictImplics(scip, objimplics->objvars, start, end, bdchgidx, addedvars, cutoffbound, &foundimplics) );
   }
   else
   {
      SCIP_Real glbbound;
      SCIP_Real newbound;
      SCIP_Real objval;

      objval = SCIPvarGetObj(var);

      assert(!SCIPisZero(scip, objval));

      if( objval > 0.0 )
      {
         newbound = SCIPvarGetUbAtIndex(var, bdchgidx, TRUE);
         glbbound = SCIPvarGetLbGlobal(var);
      }
      else
      {
         newbound = SCIPvarGetLbAtIndex(var, bdchgidx, TRUE);
         glbbound = SCIPvarGetUbGlobal(var);
      }

      /* in case the variable is integral we just need to prove the newbound plus/minus (1 - epsilon) since the this bound
       * would be rounded to newbound due to integrability which is global information
       */
      if( SCIPvarIsIntegral(var) )
      {
         if( objval > 0.0 )
            newbound += 1 - 10 * SCIPfeastol(scip);
         else
            newbound -= 1 - 10 * SCIPfeastol(scip);
      }

      /* adjust the cutoff bound by the portion the inference variable contributes to the presudo objective activitiy
       * (minactivity)
       */
      assert(!SCIPisNegative(scip, objval * (newbound - glbbound)));
      (*cutoffbound) -= objval * (newbound - glbbound);
   }

   return SCIP_OKAY;
}


/** resolves a propagation by supplying the variables whose bound changes increased the pseudo objective value above the
 *  cutoff bound
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Real             cutoffbound,        /**< the global cutoff */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL for conflict analysis initialization */
   int                   inferinfo,          /**< inference information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   SCIP_HASHTABLE* addedvars;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real glbpseudoobjval;
   SCIP_Real reqpseudoobjval;
   SCIP_Bool infinity;
   int nvars;
   int v;

   infinity = FALSE;
   addedvars = NULL;
   nvars = propdata->nminactvars;

   /* the global pseudo value gives us a global valid minimal activity
    *
    * @note The global pseudo objective activity can be minus infinity. In that case all variable are part of the
    *       reason/explanation
    *
    * @note If the global pseudo objective activity is greater than the required minactivity, the local bound change
    *       which has to explained is actually (now) a global one. That means, the reason/explanation is empty
    */
   glbpseudoobjval = SCIPgetGlobalPseudoObjval(scip);

   if( SCIPisInfinity(scip, -glbpseudoobjval) )
   {
      infinity = TRUE;
      reqpseudoobjval = cutoffbound;
   }
   else
   {
      /* clear hash table for storing variables which are not needed to add the reason due to global implications or
       * already added
       */
      if( nvars > 0 )
      {
         addedvars = propdata->addedvars;
         SCIPhashtableClear(addedvars);
      }

      if( infervar != NULL )
      {
         SCIP_CALL( adjustCutoffbound(scip, propdata, infervar, inferinfo, boundtype, bdchgidx, addedvars, &cutoffbound) );
      }

      reqpseudoobjval = cutoffbound - glbpseudoobjval;
   }

   SCIPdebugMessage("resolve propagation global pseudo objective <%g>, cutoff bounda <%g>, required minactivity <%g>\n",
      glbpseudoobjval, cutoffbound, reqpseudoobjval);

   /* the variables responsible for the propagation are the ones with
    *  - obj > 0 and local lb > global lb
    *  - obj < 0 and local ub < global ub
    *
    * collect all variables which contribute positively to the pseudo objective value (minimum activity) until we
    * reached the (adjusted) required minimum activity for the inference bound chnage
    */

   /* first consider the binary variables */
   if( nvars > 0 )
   {
      SCIP_OBJIMPLICS** minactimpls;

      vars = propdata->minactvars;
      assert(vars != NULL);

      minactimpls = propdata->minactimpls;
      assert(minactimpls != NULL);

#ifndef NDEBUG
      checkGlbfirstnonfixed(scip, propdata);
#endif

      if( infinity )
      {
         /* if the required minimum activity is minus infinity, we have to add all variables which contribute the local
          * prseudo objective activity
          */

         for( v = propdata->glbfirstnonfixed; v < nvars; ++v )
         {
            var = vars[v];
            assert(var != NULL);

            /* @note binary variables can have a zero objective value */

            if( var == infervar )
               continue;

            SCIP_CALL( addConflictBinvar(scip, var, bdchgidx, NULL, NULL, FALSE, &reqpseudoobjval) );
         }
      }
      else
      {
         assert(addedvars != NULL);

         for( v = propdata->glbfirstnonfixed; v < nvars && SCIPisPositive(scip, reqpseudoobjval); ++v )
         {
            var = vars[v];
            assert(var != NULL);

            /* @note binary variables can have a zero objective value */

            if( var == infervar )
               continue;

            if( SCIPhashtableExists(addedvars, (void*)var) )
               continue;

            SCIP_CALL( addConflictBinvar(scip, var, bdchgidx, minactimpls[v], addedvars, propdata->respropuseimplics, &reqpseudoobjval) );
         }
      }
   }

   vars = propdata->objintvars;
   nvars = propdata->nobjintvars;
   assert(nvars == 0 || vars != NULL);

   /* second consider the none binary variables */
   for( v = 0; v < nvars && (infinity || SCIPisPositive(scip, reqpseudoobjval)); ++v )
   {
      var = vars[v];
      assert(var != NULL);
      assert(!SCIPisZero(scip, SCIPvarGetObj(var)));

      if( var == infervar )
         continue;

      SCIP_CALL( addConflictBounds(scip, var, bdchgidx, &reqpseudoobjval) );
   }

   return SCIP_OKAY;
}

/** propagates the given variable against the cutoff bound */
static
SCIP_RETCODE propagateCutoffboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator, or NULL */
   SCIP_VAR*             var,                /**< variable to propagate */
   int                   inferinfo,          /**< inference information to store with the bound change */
   SCIP_Real             objchg,             /**< objective change */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   SCIP_Bool             local,              /**< local or global propagation */
   SCIP_Bool*            tightened           /**< pointer to store if the variable domain was tightened */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real newbd;
   SCIP_Bool infeasible;

   assert(!SCIPisZero(scip, objchg));
   assert(!SCIPisInfinity(scip, -pseudoobjval));
   assert(!SCIPisInfinity(scip, cutoffbound));
   assert(SCIPisLT(scip, pseudoobjval, cutoffbound) );
   assert(tightened != NULL);

   *tightened = FALSE;

   /* collect bounds of the variable */
   if( local )
   {
      assert(prop != NULL);
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

   /* depending on the objective contribution we can try to tighten the lower or upper bound of the variable */
   if( objchg > 0.0 )
   {
      newbd = lb + (cutoffbound - pseudoobjval) / objchg;

      if( local )
      {
         SCIP_CALL( SCIPinferVarUbProp(scip, var, newbd, prop, inferinfo, FALSE, &infeasible, tightened) );
         assert(!infeasible);

         if( *tightened ) /* might not be tightened due to numerical reasons */
         {
            SCIPdebugMessage(" -> new (local) upper bound of variable <%s>[%g,%g]: %g, objective change <%g>, pseudo objective <%g>, cutoff bound <%g>\n",
               SCIPvarGetName(var), lb, ub, newbd, objchg, pseudoobjval, cutoffbound);
         }
      }
      else
      {
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newbd, FALSE, &infeasible, tightened) );
         assert(!infeasible);

         if( *tightened )
         {
            SCIPdebugMessage(" -> new (global) upper bound of variable <%s>[%g,%g]: %g, objective change <%g>, pseudo objective <%g>, cutoff bound <%g>\n",
               SCIPvarGetName(var), lb, ub, newbd, objchg, pseudoobjval, cutoffbound);
         }
      }
   }
   else
   {
      newbd = ub + (cutoffbound - pseudoobjval) / objchg;

      if( local )
      {
         SCIP_CALL( SCIPinferVarLbProp(scip, var, newbd, prop, inferinfo, FALSE, &infeasible, tightened) );
         assert(!infeasible);

         if( *tightened ) /* might not be tightened due to numerical reasons */
         {
            SCIPdebugMessage(" -> new (local) lower bound of variable <%s>[%g,%g]: %g, objective change <%g>, pseudo objective <%g>, cutoff bound <%g>\n",
               SCIPvarGetName(var), lb, ub, newbd, objchg, pseudoobjval, cutoffbound);
         }
      }
      else
      {
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newbd, FALSE, &infeasible, tightened) );
         assert(!infeasible);

         if( *tightened )
         {
            SCIPdebugMessage(" -> new (global) lower bound of variable <%s>[%g,%g]: %g, objective change <%g>, pseudo objective <%g>, cutoff bound <%g>\n",
               SCIPvarGetName(var), lb, ub, newbd, objchg, pseudoobjval, cutoffbound);
         }
      }
   }

   return SCIP_OKAY;
}

/** propagates the given binary variable against the cutoff bound */
static
SCIP_RETCODE propagateCutoffboundBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator, or NULL */
   SCIP_VAR*             var,                /**< variable to propagate */
   int                   pos,                /**< position of the variable in the corresponding propdata variable array */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   SCIP_Bool*            tightened,          /**< pointer to store if the variable domain was tightened */
   SCIP_Bool             local               /**< propagate local bounds, otherwise global bounds */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_OBJIMPLICS* objimplics;
   SCIP_Real lbobjchg;
   SCIP_Real ubobjchg;
   SCIP_Real objchg;

   assert(SCIPvarIsBinary(var));

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   objimplics = propdata->minactimpls[pos];
   assert(objimplics != NULL);

   /* get objective change in case of fixing the variable to its lower bound (that is zero) */
   SCIP_CALL( getMinactObjchg(scip, var, objimplics, NULL, SCIP_BOUNDTYPE_LOWER, local, &lbobjchg) );
   assert(!SCIPisNegative(scip, lbobjchg));

   /* get objective change in case of fixing the variable to its upper bound (that is one) */
   SCIP_CALL( getMinactObjchg(scip, var, objimplics, NULL, SCIP_BOUNDTYPE_UPPER, local, &ubobjchg) );
   assert(!SCIPisNegative(scip, ubobjchg));

   (*tightened) = FALSE;

   /* nothing can be done if the objective contribution is zero independently of the bound */
   if( SCIPisZero(scip, lbobjchg) && SCIPisZero(scip, ubobjchg) )
      return SCIP_OKAY;

   /* try to tighten the variable bound use the larger objective contribution */
   if( lbobjchg > ubobjchg )
      objchg = -lbobjchg;
   else
      objchg = ubobjchg;

   SCIP_CALL( propagateCutoffboundVar(scip, prop, var, pos, objchg, cutoffbound, pseudoobjval, local, tightened) );

   return SCIP_OKAY;
}

/** globally propagates if a new cutoff bound or global pseudo objective value (minimum activity of the objective
 *  function) is available
 */
static
SCIP_RETCODE propagateCutoffboundGlobally(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** minactvars;
   SCIP_VAR** objintvars;
   SCIP_VAR* var;
   SCIP_Bool tightened;
   SCIP_Real pseudoobjval;
   SCIP_Real cutoffbound;
   int nminactvars;
   int nobjintvars;
   int v;

   /* this method should not be called in the root node of the search tree since the standard propagation already does
    * the job
    */
   assert(SCIPgetDepth(scip) > 0);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   pseudoobjval = SCIPgetGlobalPseudoObjval(scip);
   cutoffbound = propdata->cutoffbound;

   /* nothing can be done if the global pseudo objective is minus infinity */
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;

   /* check if the global pseudo objective value (minimum activity of the objective function) is greater or equal to
    * the cutoff bound
    */
   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   minactvars = propdata->minactvars;
   objintvars = propdata->objintvars;
   nminactvars = propdata->nminactvars;
   nobjintvars = propdata->nobjintvars;

#ifndef NDEBUG
   checkGlbfirstnonfixed(scip, propdata);
#endif

   /* always propagate the binary variables completely */
   for( v = propdata->glbfirstnonfixed; v < nminactvars; ++v )
   {
      var = minactvars[v];
      assert(var != NULL);

      /* check if the variables is already globally fixed; if so continue with the potential candidate */
      if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, FALSE) );

      /* the binary variables are sorted in non-increasing manner w.r.t. the absolute value of their objective
       * contribution w.r.t. minimum activity (pseudo objective value) of the objective function; these values are the
       * increase in the pseudo objective activity we would get if we fix the variable to its worse bound; hence, we can
       * stop if for a variable this potential increase is not enough too exceed the cutoff bound;
       */
      if( !tightened )
      {
         SCIPdebugMessage("interrupt global pseudo objective propagation w.r.t. cutoff bound <%.15g> for binary variables after %d from %d binary variables\n",
            cutoffbound, v, nminactvars);
         break;
      }

      /* @note The variable might not be globally fixed right away since this would destroy the local internal
       *       data structure of a search node; the bound change is in that case pending; hence we cannot assert
       *       that the variable is globally fixed
       */
      (*nchgbds)++;
   }
   propdata->glbfirstnonfixed = v;
   propdata->firstnonfixed = MAX(propdata->firstnonfixed, v);

   /* check all binary variables which could potentially be fixed */
   for( ; v < nminactvars && cutoffbound - pseudoobjval <  propdata->minactimpls[v]->maxobjchg; ++v )
   {
      var = minactvars[v];
      assert(var != NULL);

      /* check if the variables is already globally fixed; if so continue with the potential candidate */
      if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, FALSE) );

      /* check if the domain of the variable was reduced */
      if( tightened )
         (*nchgbds)++;
   }

#ifndef NDEBUG
   /* check that the abort criteria for the binary variables works */
   for( ; v < nminactvars; ++v )
   {
      assert(cutoffbound - pseudoobjval >=  propdata->minactimpls[v]->maxobjchg);

      var = minactvars[v];
      assert(var != NULL);

      /* check if the variable is already locally fixed; in that case we just continue with the next potential
       * candidate
       */
      if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, FALSE) );
      assert(!tightened);
   }
#endif

   /* propagate the none binary variables completely */
   for( v = 0; v < nobjintvars; ++v )
   {
      var = objintvars[v];
      assert(var != NULL);
      assert(!SCIPisZero(scip, SCIPvarGetObj(var)));

      /* try to tighten the bound of the variable */
      SCIP_CALL( propagateCutoffboundVar(scip, NULL, var, -1, SCIPvarGetObj(var), cutoffbound, pseudoobjval, FALSE, &tightened) );

      /* check if the domain of the variable was reduced */
      if( tightened )
         (*nchgbds)++;
   }

   propdata->glbpropagated = TRUE;

   return SCIP_OKAY;
}

/** propagates the cutoff bound for binary variables (c*x <= cutoff) */
static
SCIP_RETCODE propagateCutoffboundBinvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   int*                  nfixedvars          /**< pointer to store the number of fixed variables */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** minactvars;
   SCIP_VAR* var;
   SCIP_Bool tightened;
   int nminactvars;
   int v;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   minactvars = propdata->minactvars;
   nminactvars = propdata->nminactvars;
   assert(nminactvars == 0 || minactvars != NULL);

   /* always propagate the binary variables completely; note that the variables before the firstnonfixed indexed are
    * already locally fixed and those before glbfirstnonfixed are already globally fixed
    */

#ifndef NDEBUG
   /* check that the variables before glbfirstnonfixed are globally fixed */
   checkGlbfirstnonfixed(scip, propdata);

   /* check that the variables before firstnonfixed are locally fixed */
   for( v = propdata->glbfirstnonfixed; v < propdata->firstnonfixed; ++v )
   {
      var =  minactvars[v];
      assert(var != NULL);

      assert(SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5);
   }
#endif

   for( v = propdata->firstnonfixed; v < nminactvars; ++v )
   {
      var =  minactvars[v];
      assert(var != NULL);

      /* check if the variable is already locally fixed; in that case we just continue with the next potential
       * candidate
       */
      if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, TRUE) );

      /* the binary variables are sorted in non-increasing manner w.r.t. the absolute value of their objective
       * contribution w.r.t. minimum activity of the objective function; These values are the increase in the pseudo
       * objective activity (minimum activity of the objective function) we would get if we fix the variable to its
       * worse bound; hence, we can stop if for a variable this potential increase is not enough too exceed the cutoff
       * bound;
       */
      if( !tightened )
      {
         SCIPdebugMessage("interrupt local pseudo objective propagation w.r.t. cutoff bound <%.15g> for binary variables after %d from %d binary variables\n",
            cutoffbound, v, nminactvars);
         break;
      }

      /* check that the binary variable is locally fixed */
      assert(SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5);
      (*nfixedvars)++;
   }
   propdata->firstnonfixed = v;

   /* check all binary variables which could potentially be fixed */
   for( ; v < nminactvars && cutoffbound - pseudoobjval < propdata->minactimpls[v]->maxobjchg; ++v )
   {
      var =  minactvars[v];
      assert(var != NULL);

      /* check if the variable is already locally fixed; in that case we just continue with the next potential
       * candidate
       */
      if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, TRUE) );

      if( tightened )
      {
         assert(SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5);
         (*nfixedvars)++;
      }
   }

#ifndef NDEBUG
   /* check that the abort criteria for the binary variables works */
   for( ; v < nminactvars; ++v )
   {
      var = minactvars[v];
      assert(var != NULL);

      assert(cutoffbound - pseudoobjval >= propdata->minactimpls[v]->maxobjchg);

      /* check if the variable is already locally fixed; in that case we just continue with the next potential
       * candidate
       */
      if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5)
         continue;

      /* propagates the cutoff bound for the given binary variable */
      SCIP_CALL( propagateCutoffboundBinvar(scip, prop, var, v, cutoffbound, pseudoobjval, &tightened, TRUE) );
      assert(!tightened);
   }
#endif

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
   SCIP_Real pseudoobjval;
   SCIP_Real cutoffbound;
   SCIP_Bool cutoff;
   SCIP_Bool tightened;
   int nchgbds;

   assert(result != NULL);

   (*result) = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* get current pseudo objective value (minimum activity of the objective function) and cutoff bound */
   pseudoobjval = SCIPgetPseudoObjval(scip);
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;
   cutoffbound = SCIPgetCutoffbound(scip);
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;

   /* @note A new global pseudo objective value could be used to retrive global fixings. There is, however, no need to
    *       check if a new global pseudo objective value is available. This is the case since a new (better) global
    *       pseudo activity implicis that a global bound change was performed. That causes that the root node of the
    *       search tree get marked for repropagation. That will result in propagation call of the pseudo objective
    *       propagator.
    */

   /* check current cutoff bound */
   if( cutoffbound < propdata->cutoffbound )
   {
      propdata->glbpropagated = FALSE;
      propdata->cutoffbound = cutoffbound;
   }

   nchgbds = 0;
   cutoff = FALSE;
   (*result) = SCIP_DIDNOTFIND;

   /* check if we have a new cutoff bound; in that case we globally propagate this new bound
    *
    * @note there is no need to propagate the cutoff pound if we are in the root node since this will be done by the
    *       standard local propagation
    */
   if( propdata->propcutoffbound && !propdata->glbpropagated && SCIPgetDepth(scip) > 0 )
   {
      /* first globally propagate new cutoff bound or pseudo objective activity */
      SCIP_CALL( propagateCutoffboundGlobally(scip, prop, &nchgbds, &cutoff) );

      if( cutoff )
      {
         /* we are done with solving since a global pseudo activity is greater or equal to the cutoff bound */
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );

         (*result) = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* check if the global propagation cut off the active/current node */
      if( SCIPgetCutoffdepth(scip) <= SCIPgetDepth(scip) )
      {
         (*result) = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   /* check if the pseudo objective value (minimum activity of the objective function) is greater or equal to the cutoff
    * bound
    */
   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      SCIPdebugMessage("pseudo objective value <%g> exceeds cutoff bound <%g>\n", pseudoobjval, cutoffbound);

      /* check if conflict analysis is applicable */
      if( SCIPisConflictAnalysisApplicable(scip) )
      {
         assert(SCIPgetDepth(scip) > 0);

         /* initialize conflict analysis */
         SCIP_CALL( SCIPinitConflictAnalysis(scip) );

         /* add all variable whose best bound changes increased the pseudo objective value above to cutoff bound */
         SCIP_CALL( resolvePropagation(scip, propdata, cutoffbound, NULL, -1, SCIP_BOUNDTYPE_UPPER, NULL) );

         /* analyze the conflict */
         SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
      }
      (*result) = SCIP_CUTOFF;

      return SCIP_OKAY;
   }

   SCIPdebugMessage("propagating pseudo objective function (pseudoobj: %g, cutoffbound: %g)\n", pseudoobjval, cutoffbound);

   /* propagate binary variables */
   SCIP_CALL( propagateCutoffboundBinvars(scip, prop, cutoffbound, pseudoobjval, &nchgbds) );

   /* tighten domains of none binary variables, if they would increase the pseudo objective value above the cutoff
    * bound
    */
   if( propdata->propfullinroot && SCIPgetDepth(scip) == 0 )
   {
      SCIP_VAR** objintvars;
      SCIP_VAR* var;
      SCIP_Real objval;
      int nobjintvars;
      int v;

      objintvars = propdata->objintvars;
      nobjintvars = propdata->nobjintvars;
      assert(nobjintvars == 0 || objintvars != NULL);

      /* propagate all none binary variables */
      for( v = 0; v < nobjintvars; ++v )
      {
         var = objintvars[v];
         assert(var != NULL);

         objval = SCIPvarGetObj(var);
         assert(!SCIPisZero(scip, objval));

         /* try to tighten the bound of the variable */
         SCIP_CALL( propagateCutoffboundVar(scip, NULL, var, -1, objval, cutoffbound, pseudoobjval, FALSE, &tightened) );

         /* check if the domain of the variable was reduced */
         if( tightened )
            nchgbds++;
      }
   }
   else
   {
      SCIP_VAR** objintvars;
      SCIP_VAR* var;
      SCIP_Real objval;
      int nobjintvars;
      int nmaxuseless;
      int nuseless;
      int c;
      int v;

      objintvars = propdata->objintvars;
      nobjintvars = propdata->nobjintvars;
      assert(nobjintvars == 0 || objintvars != NULL);

      /* compute maximum number of useless propagations before aborting */
      nmaxuseless = MAX(propdata->minuseless, (int)propdata->maxvarsfrac*(nobjintvars));

      nuseless = 0;
      v = propdata->lastvarnum;

      for( c = 0; c < nobjintvars && nuseless < nmaxuseless; ++c )
      {
         v++;
         if( v >= nobjintvars )
            v = 0;

         var = objintvars[v];
         assert(var != NULL);

         objval = SCIPvarGetObj(var);
         assert(!SCIPisZero(scip, objval));

         /* try to tighten the bound of the variable */
         SCIP_CALL( propagateCutoffboundVar(scip, prop, var, -1, objval, cutoffbound, pseudoobjval, TRUE, &tightened) );

         /* check if the domain of the variable was reduced */
         if( tightened )
            nchgbds++;
         else
            nuseless++;
      }
      propdata->lastvarnum = v;
   }

   /* check if we chanced bounds */
   if( nchgbds > 0 )
      (*result) = SCIP_REDUCEDDOM;

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
   SCIP_Real objval;
   SCIP_Real contrib;
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
      objval = SCIPvarGetObj(vars[v]);
      if( SCIPisPositive(scip, objval) )
      {
         contrib = SCIPvarGetUbGlobal(vars[v]);
         if( !SCIPisInfinity(scip, contrib) )
            contrib *= objval;
      }
      else if( SCIPisNegative(scip, objval) )
      {
         contrib = SCIPvarGetLbGlobal(vars[v]);
         if( !SCIPisInfinity(scip, -contrib) )
            contrib *= objval;
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

/** updates the pseudo objective activity if necessary */
static
void updateMaxObjPseudoactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   /* if necessary, calculate the maximum pseudo objective activity */
   if( propdata->maxpseudoobjact == SCIP_INVALID ) /*lint !e777*/
      calcMaxObjPseudoactivity(scip, propdata);
   assert(propdata->maxpseudoobjact != SCIP_INVALID); /*lint !e777*/
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
   SCIP_Real objval;
   SCIP_Real contrib;

   assert(propdata != NULL);

   contrib = 0.0;
   objval = SCIPvarGetObj(var);
   if( SCIPvarGetWorstBoundType(var) == SCIP_BOUNDTYPE_UPPER )
   {
      contrib = SCIPvarGetUbGlobal(var);
      if( !SCIPisInfinity(scip, contrib) )
         contrib *= objval;
   }
   else
   {
      assert(SCIPvarGetWorstBoundType(var) == SCIP_BOUNDTYPE_LOWER);
      contrib = SCIPvarGetLbGlobal(var);
      if( !SCIPisInfinity(scip, -contrib) )
         contrib *= objval;
      else
         contrib *= -1.0;
   }

   return getMaxObjPseudoactivityResidualValue(scip, propdata, contrib);
}

/** returns the maximum pseudo objective activity of the objective function */
static
SCIP_Real getMaxObjPseudoactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   return getMaxObjPseudoactivityResidualValue(scip, propdata, 0.0);
}

/** propagates the global domain of the given binary variable against the lower bound (c*x >= lowerbound) */
static
SCIP_RETCODE propagateLowerboundBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to propagate */
   SCIP_Real             lowerbound,         /**< lower bound to use */
   SCIP_Real             maxpseudoobjact,    /**< maximum pseudo objective activity */
   SCIP_Bool             useimplics,         /**< should implications be used */
   SCIP_Bool*            infeasible,         /**< pointer to store if the variable domain got empty, infeasible */
   SCIP_Bool*            tightened           /**< pointer to store if the variable domain was tightened */
   )
{
   SCIP_Real lbobjchg;
   SCIP_Real ubobjchg;

   assert(SCIPvarIsBinary(var));
   assert(SCIPisLE(scip, lowerbound, maxpseudoobjact));

   /* get contribution of variable by fixing it to its lower bound w.r.t. maximum activity of the objective function */
   lbobjchg = getMaxactObjchg(var, SCIP_BOUNDTYPE_LOWER, useimplics);
   assert(!SCIPisPositive(scip, lbobjchg));

   /* get contribution of variable by fixing it to its upper bound w.r.t. maximum activity of the objective function */
   ubobjchg = getMaxactObjchg(var, SCIP_BOUNDTYPE_UPPER, useimplics);
   assert(!SCIPisPositive(scip, ubobjchg));

   (*infeasible) = FALSE;
   (*tightened) = FALSE;

   /* if the maximum activity of the objective function without the contribution of the given variable shrinks below the
    * global lower bound, the contribution of the variable is need; hence, we can fix it to corresponding bound globally
    */
   if( maxpseudoobjact + lbobjchg < lowerbound )
   {
      SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, 1.0, FALSE, infeasible, tightened) );
   }
   else if( maxpseudoobjact + ubobjchg < lowerbound )
   {
      SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, 0.0, FALSE, infeasible, tightened) );
   }

   return SCIP_OKAY;
}

/** propagates the global domains of the given variable with non-zero objective coefficient against the lower bound
 *  (c*x >= lowerbound)
 */
static
SCIP_RETCODE propagateLowerboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable to propagate */
   SCIP_Real             lowerbound,         /**< lower bound to use */
   SCIP_Bool*            infeasible,         /**< pointer to store if the variable domain got empty, infeasible */
   SCIP_Bool*            tightened           /**< pointer to store if the variable domain was tightened */
   )
{
   SCIP_Real residual;
   SCIP_Real newbd;
   SCIP_Real objval;

   objval = SCIPvarGetObj(var);
   assert(!SCIPisZero(scip, objval));

   /* get residual pseudo objective activity, that is the pseudo objective activity without the given variable */
   residual = getMaxObjPseudoactivityResidual(scip, propdata, var);

   if( SCIPisInfinity(scip, residual) )
      return SCIP_OKAY;

   /* compute potential mew bound */
   newbd = (lowerbound - residual) / objval;

   if( objval > 0.0 )
   {
      SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newbd, FALSE, infeasible, tightened) );
   }
   else
   {
      SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newbd, FALSE, infeasible, tightened) );
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
   SCIP_Real maxpseudoobjact;
   SCIP_Bool cutoff;
   int nchgbds;

   assert(result != NULL);

   (*result) = SCIP_DIDNOTRUN;
   cutoff = FALSE;
   nchgbds = 0;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->nminactvars > 0 || propdata->nobjintvars > 0);

   /* check if there is a chance to find a reduction */
   lowerbound = SCIPgetLowerbound(scip);

   if( SCIPisInfinity(scip, -lowerbound) )
      return SCIP_OKAY;

   /* if the lower bound did not change since the last propagation as well as the global bounds of the variables with a
    * non-zero objective coefficient we do nothing since there is no new information available
    */
   if( SCIPisLE(scip, lowerbound, propdata->lastlowerbound) && propdata->maxpseudoobjact < SCIP_INVALID )
      return SCIP_OKAY;

   /* updates the pseudo objective activity if necessary */
   updateMaxObjPseudoactivity(scip, propdata);

   /* if more than one variable contributes an infinity to the maximal pseudo activity we can do nothing */
   assert(propdata->maxpseudoobjact < SCIP_INVALID);
   if( propdata->maxpseudoobjactinf > 1 )
      return SCIP_OKAY;

   maxpseudoobjact = getMaxObjPseudoactivity(scip, propdata);

#ifndef NDEBUG
   /* check that the global indices are correct */
   checkGlbfirstnonfixed(scip, propdata);
#endif

   /* if the maximum pseudo objective activity is smaller than the lower bound the problem is infeasible */
   if( SCIPisLT(scip, maxpseudoobjact, lowerbound) )
      cutoff = TRUE;
   else
   {
      SCIP_VAR** objintvars;
      SCIP_VAR* var;
      SCIP_Bool tightened;
      int nobjintvars;
      int v;

      if( propdata->maxpseudoobjactinf == 0 )
      {
         SCIP_VAR** maxactvars;
         int nmaxactvars;

         maxactvars = propdata->maxactvars;
         nmaxactvars = propdata->nmaxactvars;
         assert(nmaxactvars == 0 || maxactvars != NULL);

         for( v = propdata->maxactfirstnonfixed; v < nmaxactvars && !cutoff; ++v )
         {
            var = maxactvars[v];
            assert(var != NULL);

            /* check if the variables is already globally fixed; if so continue with the next potential candidate */
            if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
               continue;

            /* try to propagate variable domain globally */
            SCIP_CALL( propagateLowerboundBinvar(scip, var, lowerbound, maxpseudoobjact, propdata->propuseimplics, &cutoff, &tightened) );

            /* the binary variables are sorted in non-increasing manner w.r.t. the absolute value of their objective
             * contribution w.r.t. maximum activity of the objective function; These values are the decrease we would
             * get with the maximum pseudo objective activity if we fix the variable to its best bound; hence, we can
             * stop if for a variable this potential decrease is not enough anymore too fall below the lower bound.
             *
             * @note In case a fixing was performed. The variable might not be globally fixed right away since this would
             *       destroy the local internal data structure of a search node; the bound change is in that case pending;
             *       hence we cannot assert that the variable is globally fixed
             */
            if( !tightened )
            {
               assert(!SCIPisInfinity(scip, propdata->maxpseudoobjact));
               SCIPdebugMessage("interrupt pseudo objective propagation w.r.t. lower bound <%.15g> for binary variables after %d from %d binary variables\n",
                  lowerbound, v, nmaxactvars);
               break;
            }

            nchgbds++;
         }

         /* update globally fixed index if abort criteria was applied */
         propdata->maxactfirstnonfixed = v;

         /* check all binary variables which could potentially be fixed */
         for( ; v < nmaxactvars && maxpseudoobjact - lowerbound < propdata->maxactchgs[v] ; ++v )
         {
            var =  maxactvars[v];
            assert(var != NULL);

            /* check if the variables is already globally fixed; if so continue with the potential candidate */
            if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
               continue;

            /* propagates the cutoff bound for the given binary variable */
            SCIP_CALL( propagateLowerboundBinvar(scip, var, lowerbound, maxpseudoobjact, propdata->propuseimplics, &cutoff, &tightened) );

            if( tightened )
               nchgbds++;
         }

#ifndef NDEBUG
         /* check that the abort criteria for the binary variables works */
         for( ; v < nmaxactvars && !cutoff; ++v )
         {
            var = maxactvars[v];
            assert(var != NULL);

            /* check if the variables is already globally fixed; if so continue with the next potential candidate */
            if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
               continue;

            /* try to propagate variable domain globally */
            SCIP_CALL( propagateLowerboundBinvar(scip, var, lowerbound, maxpseudoobjact, propdata->propuseimplics, &cutoff, &tightened) );
            assert(!tightened);
            assert(!cutoff);
         }
#endif
      }

      objintvars = propdata->objintvars;
      nobjintvars = propdata->nobjintvars;
      assert(nobjintvars == 0 || objintvars != NULL);

      /* propagate c*x >= lowerbound for non-binary variables */
      for( v = 0; v < nobjintvars && !cutoff; ++v )
      {
         var = objintvars[v];
         assert(var != NULL);
         assert(!SCIPisZero(scip, SCIPvarGetObj(var)));

         /* try to propagate variable domain globally */
         SCIP_CALL( propagateLowerboundVar(scip, propdata, var, lowerbound, &cutoff, &tightened) );

         if( tightened )
            nchgbds++;
      }
   }

   /* evaluate propagation results */
   if( cutoff )
   {
      /* we are done with solving since a global bound change is infeasible: cutoff root node */
      SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
      (*result) = SCIP_CUTOFF;
   }
   else if( nchgbds > 0 )
      (*result) = SCIP_REDUCEDDOM;

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

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   (*result) = SCIP_DIDNOTRUN;

   /* do nothing if active pricer are present and force flag is not TRUE */
   if( !propdata->force && SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   pseudoobjval = SCIPgetGlobalPseudoObjval(scip);
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;

   cutoffbound = SCIPgetCutoffbound(scip);
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;

   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      (*result) = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* only propagate if a new cutoff bound or globale pseudo objective value is available */
   if( cutoffbound < propdata->cutoffbound || pseudoobjval > propdata->glbpseudoobjval )
   {
      SCIP_Real objval;
      SCIP_Bool tightened;

      (*result) = SCIP_DIDNOTFIND;
      oldnchgbds = *nchgbds;

      /* get the problem variables */
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      /* scan the variables for pseudoobj bound reductions
       * (loop backwards, since a variable fixing can change the current and the subsequent slots in the vars array)
       */
      for( v = nvars - 1; v >= 0; --v )
      {
         objval = SCIPvarGetObj(vars[v]);

         if( SCIPisZero(scip, objval) )
            continue;

         SCIP_CALL( propagateCutoffboundVar(scip, NULL, vars[v], -1, objval, cutoffbound, pseudoobjval, FALSE, &tightened) );

         if( tightened )
            (*nchgbds)++;
      }

      /* evaluate if bound change was detected */
      if( *nchgbds > oldnchgbds )
         (*result) = SCIP_SUCCESS;

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

   (*result) = SCIP_DIDNOTRUN;

   /* do nothing if active pricer are present and force flag is not TRUE */
   if( !propdata->force && SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* check if enough new variable are added (due to column generatition to reinitialized the propgator data) */
   if( propdata->nnewvars > propdata->maxnewvars )
   {
      /* free current propdata data */
      SCIP_CALL( propdataExit(scip, propdata) );

      /* initialize propdata data from scratch */
      SCIP_CALL( propdataInit(scip, propdata) );
   }

   /* nothing to do for zero objective */
   if( propdata->nminactvars == 0 && propdata->nmaxactvars == 0 && propdata->nobjintvars == 0 )
      return SCIP_OKAY;

   /* propagate c*x <= cutoff */
   SCIP_CALL( propagateCutoffbound(scip, prop, result) );

   if( (*result) != SCIP_CUTOFF && (propdata->nmaxactvars > 0 || propdata->nobjintvars > 0) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_RESULT dualresult;

      /* propagates the global lower (dual) bound c*x >= lowerbound */
      SCIP_CALL( propagateLowerbound(scip, prop, &dualresult) );

      if( dualresult == SCIP_REDUCEDDOM || dualresult == SCIP_CUTOFF )
         (*result) = dualresult;
   }

   return SCIP_OKAY;
}

/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Real cutoffbound;

   assert(!SCIPisEQ(scip, SCIPvarGetLbGlobal(infervar), SCIPvarGetUbGlobal(infervar)));

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   cutoffbound = SCIPgetCutoffbound(scip);
   assert(!SCIPisInfinity(scip, cutoffbound));
   assert(infervar != NULL);

   SCIPdebugMessage("resolve bound change <%s> %s <%g>(%g), cutoff bound <%g>\n", SCIPvarGetName(infervar),
      boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE),
      SCIPvarGetLbAtIndex(infervar, bdchgidx, FALSE), cutoffbound);

   /* resolve the propagation of the inference variable w.r.t. required minactivity */
   SCIP_CALL( resolvePropagation(scip, propdata, cutoffbound, infervar, inferinfo, boundtype, bdchgidx) );

   (*result) = SCIP_SUCCESS;

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

   /* reset propagator data structure */
   propdataReset(scip, propdata);

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
         "do we want to propagate all none binary variables if we are propagating the root node",
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

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/propuseimplics",
         "use implications to strengthen the propagation of binary variable (increasing the objective change)?",
         &propdata->propuseimplics, TRUE, DEFAULT_PROPUSEIMPLICS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/respropuseimplics",
         "use implications to strengthen the resolve propagation of binary variable (increasing the objective change)?",
         &propdata->respropuseimplics, TRUE, DEFAULT_RESPROPUSEIMPLICS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/"PROP_NAME"/maximplvars",
         "maximum number of binary variables the implications are used if turned on (-1: unlimited)?",
         &propdata->maximplvars, TRUE, DEFAULT_MAXIMPLVARS, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** propagates the cutoff bound for the given variables */
SCIP_RETCODE SCIPpropagateCutoffboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator, or NULL */
   SCIP_VAR*             var,                /**< variables to propagate */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   SCIP_Bool*            tightened           /**< pointer to if the domain was tightened */
   )
{
   SCIP_Real objval;

   objval = SCIPvarGetObj(var);

   SCIP_CALL( propagateCutoffboundVar(scip, prop, var, -1, objval, cutoffbound, pseudoobjval, TRUE, tightened) );

   return SCIP_OKAY;
}
