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

/**@file   cons_logicor.c
 * @brief  constraint handler for logic or constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_logicor.h"
#include "cons_linear.h"


#define CONSHDLR_NAME          "logicor"
#define CONSHDLR_DESC          "logic or constraints"
#define CONSHDLR_SEPAPRIORITY   +800000
#define CONSHDLR_ENFOPRIORITY   +800000
#define CONSHDLR_CHECKPRIORITY  -800000
#define CONSHDLR_SEPAFREQ             4
#define CONSHDLR_PROPFREQ             1
#define CONSHDLR_NEEDSCONS         TRUE

#define LINCONSUPGD_PRIORITY    +800000

#define EVENTHDLR_NAME         "logicor"
#define EVENTHDLR_DESC         "bound tighten event handler for logic or constraints"

#define CONFLICTHDLR_NAME      "logicor"
#define CONFLICTHDLR_DESC      "conflict handler for logic or constraints"
#define CONFLICTHDLR_PRIORITY  LINCONSUPGD_PRIORITY

#ifdef BRANCHLP
#define MINBRANCHWEIGHT             0.3  /**< minimum weight of both sets in binary set branching */
#define MAXBRANCHWEIGHT             0.7  /**< maximum weight of both sets in binary set branching */
#endif
#define DEFAULT_NPSEUDOBRANCHES       2  /**< number of children created in pseudo branching */
#define DEFAULT_MAXVARUSEFAC        1.0  /**< branching factor to weigh maximum of positive and negative variable uses */
#define DEFAULT_MINVARUSEFAC       -0.2  /**< branching factor to weigh minimum of positive and negative variable uses */


/** constraint handler data */
struct ConsHdlrData
{
   INTARRAY*        posvaruses;         /**< number of positive literal of a variable in the active logic or constraints */
   INTARRAY*        negvaruses;         /**< number of negative literal of a variable in the active logic or constraints */
   EVENTHDLR*       eventhdlr;          /**< event handler for bound tighten events on watched variables */
   Real             maxvarusefac;       /**< branching factor to weigh maximum of positive and negative variable uses */
   Real             minvarusefac;       /**< branching factor to weigh minimum of positive and negative variable uses */
   int              npseudobranches;    /**< number of children created in pseudo branching */
};

/** logic or constraint data */
struct ConsData
{
   ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   VAR**            vars;               /**< variables of the constraint */
   int              varssize;           /**< size of vars array */
   int              nvars;              /**< number of variables in the constraint */
   int              watchedvar1;        /**< position of the first watched variable */
   int              watchedvar2;        /**< position of the second watched variable */
   int              watchedfeasvar;     /**< position of the feasible watched variable (a variable fixed to one) */
   int              watchedsolvar;      /**< position of the variable making the last solution feasible */
   unsigned int     changed:1;          /**< was constraint changed since last preprocess/propagate call? */
};




/*
 * Local methods
 */

/** creates constaint handler data for logic or constraint handler */
static
RETCODE conshdlrdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, conshdlrdata) );
   CHECK_OKAY( SCIPcreateIntarray(scip, &(*conshdlrdata)->posvaruses) );
   CHECK_OKAY( SCIPcreateIntarray(scip, &(*conshdlrdata)->negvaruses) );
   (*conshdlrdata)->maxvarusefac = DEFAULT_MAXVARUSEFAC;
   (*conshdlrdata)->minvarusefac = DEFAULT_MINVARUSEFAC;
   (*conshdlrdata)->npseudobranches = DEFAULT_NPSEUDOBRANCHES;

   /* get event handler for catching bound tighten events on watched variables */
   (*conshdlrdata)->eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      errorMessage("event handler for logic or constraints not found");
      return SCIP_PLUGINNOTFOUND;
   }
   
   return SCIP_OKAY;
}

/** frees constraint handler data for logic or constraint handler */
static
RETCODE conshdlrdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   CHECK_OKAY( SCIPfreeIntarray(scip, &(*conshdlrdata)->posvaruses) );
   CHECK_OKAY( SCIPfreeIntarray(scip, &(*conshdlrdata)->negvaruses) );
   SCIPfreeMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** increases the usage counter of the given variable */
static
RETCODE conshdlrdataIncVaruses(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   VAR*             var                 /**< variable to increase usage counter for */
   )
{
   assert(conshdlrdata != NULL);
   assert(var != NULL);

   /* check, if the literal is positive or negative, and increase the corresponding varuses counter */
   if( SCIPvarIsNegated(var) )
   {
      CHECK_OKAY( SCIPgetNegatedVar(scip, var, &var) );
      CHECK_OKAY( SCIPincIntarrayVal(scip, conshdlrdata->negvaruses, SCIPvarGetIndex(var), +1) );
   }
   else
   {
      CHECK_OKAY( SCIPincIntarrayVal(scip, conshdlrdata->posvaruses, SCIPvarGetIndex(var), +1) );
   }
   /*debugMessage("increased varuses of <%s>: %d+/%d-\n", SCIPvarGetName(var),
     SCIPgetIntarrayVal(scip, conshdlrdata->posvaruses, SCIPvarGetIndex(var)),
     SCIPgetIntarrayVal(scip, conshdlrdata->negvaruses, SCIPvarGetIndex(var)));*/

   return SCIP_OKAY;
}

/** decreases the usage counter of the given variable */
static
RETCODE conshdlrdataDecVaruses(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   VAR*             var                 /**< variable to increase usage counter for */
   )
{
   assert(conshdlrdata != NULL);
   assert(var != NULL);

   /* check, if the literal is positive or negative, and decrease the corresponding varuses counter */
   if( SCIPvarIsNegated(var) )
   {
      CHECK_OKAY( SCIPgetNegatedVar(scip, var, &var) );
      CHECK_OKAY( SCIPincIntarrayVal(scip, conshdlrdata->negvaruses, SCIPvarGetIndex(var), -1) );
   }
   else
   {
      CHECK_OKAY( SCIPincIntarrayVal(scip, conshdlrdata->posvaruses, SCIPvarGetIndex(var), -1) );
   }
   /*debugMessage("decreased varuses of <%s>: %d+/%d-\n", SCIPvarGetName(var),
     SCIPgetIntarrayVal(scip, conshdlrdata->posvaruses, SCIPvarGetIndex(var)),
     SCIPgetIntarrayVal(scip, conshdlrdata->negvaruses, SCIPvarGetIndex(var)));*/

   assert(SCIPgetIntarrayVal(scip, conshdlrdata->posvaruses, SCIPvarGetIndex(var)) >= 0);
   assert(SCIPgetIntarrayVal(scip, conshdlrdata->negvaruses, SCIPvarGetIndex(var)) >= 0);

   return SCIP_OKAY;
}

/** catches events for variable at given position */
static
RETCODE consdataCatchEvent(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< logic or constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   /* catch bound tighten events on variable */
   CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
                  (EVENTDATA*)consdata) );
   
   return SCIP_OKAY;
}

/** drops events for variable at given position */
static
RETCODE consdataDropEvent(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< logic or constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   /* drop events on variable */
   CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
                  (EVENTDATA*)consdata) );

   return SCIP_OKAY;
}

/** locks rounding for variable at given position in transformed logic or constraint */
static
void lockRounding(
   VAR*             var,                /**< problem variable */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   /*debugMessage("locking coefficient <%s> in logic or constraint\n", SCIPvarGetName(var));*/

   /* forbid rounding of variable */
   SCIPvarLock(var, nlockspos, nlocksneg);
}

/** unlocks rounding for variable at given position in transformed logic or constraint */
static
void unlockRounding(
   VAR*             var,                /**< problem variable */
   int              nunlockspos,        /**< decrease in number of rounding locks for constraint */
   int              nunlocksneg         /**< decrease in number of rounding locks for constraint's negation */
   )
{
   /*debugMessage("unlocking coefficient <%s> in logic or constraint\n", SCIPvarGetName(var));*/

   /* allow rounding of variable */
   SCIPvarUnlock(var, nunlockspos, nunlocksneg);
}

/** locks rounding for all variables in transformed logic or constraint */
static
void consdataLockAllRoundings(
   CONSDATA*        consdata,           /**< logic or constraint data */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   int i;

   assert(consdata != NULL);

   /* lock every single coefficient */
   for( i = 0; i < consdata->nvars; ++i )
      lockRounding(consdata->vars[i], nlockspos, nlocksneg);
}

/** unlocks rounding for all variables in transformed logic or constraint */
static
void consdataUnlockAllRoundings(
   CONSDATA*        consdata,           /**< logic or constraint data */
   int              nunlockspos,        /**< decrease in number of rounding locks for constraint */
   int              nunlocksneg         /**< decrease in number of rounding locks for constraint's negation */
   )
{
   int i;

   assert(consdata != NULL);

   /* unlock every single coefficient */
   for( i = 0; i < consdata->nvars; ++i )
      unlockRounding(consdata->vars[i], nunlockspos, nunlocksneg);
}

/** creates a logic or constraint data object */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the logic or constraint data */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars                /**< variables of the constraint */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->row = NULL;
   if( nvars > 0 )
   {
      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      (*consdata)->varssize = nvars;
      (*consdata)->nvars = nvars;
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->varssize = 0;
      (*consdata)->nvars = 0;
   }
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->watchedfeasvar = -1;
   (*consdata)->watchedsolvar = -1;
   (*consdata)->changed = TRUE;

   return SCIP_OKAY;
}   

/** creates a transformed logic or constraint data object */
static
RETCODE consdataCreateTransformed(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the logic or constraint data */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars                /**< variables of the constraint */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, consdata, nvars, vars) );

   /* transform the variables */
   CHECK_OKAY( SCIPgetTransformedVars(scip, nvars, (*consdata)->vars, (*consdata)->vars) );

   return SCIP_OKAY;
}

/** frees a logic or constraint data */
static
RETCODE consdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the logic or constraint */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert(eventhdlr != NULL);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* drop bound tighten events for watched variables */
   if( (*consdata)->watchedvar1 != -1 )
   {
      assert(0 <= (*consdata)->watchedvar1 && (*consdata)->watchedvar1 < (*consdata)->nvars);
      assert(SCIPvarIsTransformed((*consdata)->vars[(*consdata)->watchedvar1]));
      CHECK_OKAY( consdataDropEvent(scip, *consdata, eventhdlr, (*consdata)->watchedvar1) );
   }
   if( (*consdata)->watchedvar2 != -1 )
   {
      assert(0 <= (*consdata)->watchedvar2 && (*consdata)->watchedvar2 < (*consdata)->nvars);
      assert(SCIPvarIsTransformed((*consdata)->vars[(*consdata)->watchedvar2]));
      CHECK_OKAY( consdataDropEvent(scip, *consdata, eventhdlr, (*consdata)->watchedvar2) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints logic or constraint to file stream */
static
void consdataPrint(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< logic or constraint data */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);

   if( file == NULL )
      file = stdout;

   /* print coefficients */
   fprintf(file, "logicor(");
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      if( v > 0 )
         fprintf(file, ", ");
      fprintf(file, "%s", SCIPvarGetName(consdata->vars[v]));
   }
   fprintf(file, ")\n");
}

/** deletes coefficient at given position from logic or constraint data */
static
RETCODE delCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< position of coefficient to delete */
   )
{
   CONSDATA* consdata;

   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(consdata->vars[pos]));

   /* if necessary, update the rounding locks of variable */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));
      unlockRounding(consdata->vars[pos], (int)SCIPconsIsLockedPos(cons), (int)SCIPconsIsLockedNeg(cons));
   }

   if( SCIPconsIsTransformed(cons) )
   {
      /* if the position is watched, drop bound tighten events and stop watching the position */
      if( pos == consdata->watchedvar1 )
      {
         assert(pos != consdata->watchedvar2);
         CHECK_OKAY( consdataDropEvent(scip, consdata, eventhdlr, pos) );
         consdata->watchedvar1 = -1;
      }
      if( pos == consdata->watchedvar2 )
      {
         assert(pos != consdata->watchedvar1);
         CHECK_OKAY( consdataDropEvent(scip, consdata, eventhdlr, pos) );
         consdata->watchedvar2 = -1;
      }
      if( pos == consdata->watchedfeasvar )
         consdata->watchedfeasvar = -1;
      if( pos == consdata->watchedsolvar )
         consdata->watchedsolvar = -1;
   }
   assert(pos != consdata->watchedvar1);
   assert(pos != consdata->watchedvar2);
   assert(pos != consdata->watchedfeasvar);
   assert(pos != consdata->watchedsolvar);

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->nvars--;

   consdata->changed = TRUE;

   return SCIP_OKAY;
}

/** deletes all zero-fixed variables, checks for variables fixed to one */
static
RETCODE applyFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   Bool*            redundant           /**< returns whether a variable fixed to one exists in the constraint */
   )
{
   CONSDATA* consdata;
   VAR* var;
   int v;

   assert(eventhdlr != NULL);
   assert(redundant != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->vars != NULL);

   *redundant = FALSE;
   v = 0;
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 1.0) )
      {
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         *redundant = TRUE;
         return SCIP_OKAY;
      }
      else if( SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 0.0) )
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         CHECK_OKAY( delCoefPos(scip, cons, eventhdlr, v) );
      }
      else
         ++v;
   }

   debugMessage("after fixings: ");
   debug(consdataPrint(scip, consdata, NULL));

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint, and adds conflict clause to problem */
static
RETCODE analyzeConflict(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< logic or constraint to be analyzed */
   )
{
   int v;

   assert(consdata != NULL);

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   CHECK_OKAY( SCIPinitConflictAnalysis(scip) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      CHECK_OKAY( SCIPaddConflictVar(scip, consdata->vars[v]) );
   }

   /* analyze the conflict */
   CHECK_OKAY( SCIPanalyzeConflict(scip, NULL) );

   return SCIP_OKAY;
}

/** checks constraint for violation only looking at the watched variables, applies fixings if possible */
static
RETCODE processWatchedVars(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint to be processed */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   Bool*            addcut,             /**< pointer to store whether this constraint must be added as a cut */
   Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   CONSDATA* consdata;
   VAR** vars;
   Real lb;
   Real ub;
   int nvars;
   int watchedvar1;
   int watchedvar2;
   int v;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(reduceddom != NULL);
   assert(addcut != NULL);
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->watchedvar1 == -1 || consdata->watchedvar1 != consdata->watchedvar2);

   *addcut = FALSE;
   *mustcheck = FALSE;

   /* don't process the constraint, if the watched variables were not changed since last processing */
   if( !consdata->changed )
   {
      assert(0 <= consdata->watchedvar1 && consdata->watchedvar1 < consdata->nvars);
      assert(0 <= consdata->watchedvar2 && consdata->watchedvar2 < consdata->nvars);
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[consdata->watchedvar1]), 0.0));
      assert(SCIPisEQ(scip, SCIPvarGetUbLocal(consdata->vars[consdata->watchedvar1]), 1.0));
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[consdata->watchedvar2]), 0.0));
      assert(SCIPisEQ(scip, SCIPvarGetUbLocal(consdata->vars[consdata->watchedvar2]), 1.0));
      *mustcheck = TRUE;
      return SCIP_OKAY;
   }

   vars = consdata->vars;
   nvars = consdata->nvars;
   assert(nvars == 0 || vars != NULL);
   watchedvar1 = consdata->watchedvar1;
   watchedvar2 = consdata->watchedvar2;

   /* check, if the old watched variables are still valid */
   if( consdata->watchedfeasvar >= 0 )
   {
      lb = SCIPvarGetLbLocal(vars[consdata->watchedfeasvar]);
      if( SCIPisEQ(scip, lb, 1.0) )
      {
         /* the watched feasible variable is fixed to one, making the constraint redundant;
          * we can disable it and don't have to check it for feasibility
          */
         debugMessage("disabling constraint <%s> (watched feasible variable still fixed to 1.0)\n", SCIPconsGetName(cons));
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
         return SCIP_OKAY;
      }
   }
   if( watchedvar1 >= 0 )
   {
      lb = SCIPvarGetLbLocal(vars[watchedvar1]);
      ub = SCIPvarGetUbLocal(vars[watchedvar1]);
      if( SCIPisEQ(scip, lb, ub) )
      {
         if( SCIPisEQ(scip, lb, 1.0) )
         {
            /* the variable is fixed to one, making the constraint redundant;
             * remember the variable and disable the constraint
             */
            debugMessage("disabling constraint <%s> (watchedvar1 fixed to 1.0)\n", SCIPconsGetName(cons));
            consdata->watchedfeasvar = watchedvar1;
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
            return SCIP_OKAY;
         }
         else
         {
            /* the variable is fixed to zero and can no longer be used as watched variable */
            assert(SCIPisEQ(scip, ub, 0.0));
            watchedvar1 = -1;
         }
      }
   }
   if( watchedvar2 >= 0 )
   {
      lb = SCIPvarGetLbLocal(vars[watchedvar2]);
      ub = SCIPvarGetUbLocal(vars[watchedvar2]);
      if( SCIPisEQ(scip, lb, ub) )
      {
         if( SCIPisEQ(scip, lb, 1.0) )
         {
            /* the variable is fixed to one, making the constraint redundant;
             * remember the variable and disable the constraint
             */
            debugMessage("disabling constraint <%s> (watchedvar2 fixed to 1.0)\n", SCIPconsGetName(cons));
            consdata->watchedfeasvar = watchedvar2;
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
            return SCIP_OKAY;
         }
         else
         {
            /* the variable is fixed to zero and can no longer be used as watched variable */
            assert(SCIPisEQ(scip, ub, 0.0));
            watchedvar2 = -1;
         }
      }
   }

   /* check, if both watched variables are still unfixed */
   if( watchedvar1 >= 0 && watchedvar2 >= 0 )
   {
      /* there are at least two unfixed variables -> the constraint must be checked manually */
      assert(0 <= watchedvar1 && watchedvar1 < nvars);
      assert(0 <= watchedvar2 && watchedvar2 < nvars);
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(vars[watchedvar1]), 0.0));
      assert(SCIPisEQ(scip, SCIPvarGetUbLocal(vars[watchedvar1]), 1.0));
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(vars[watchedvar2]), 0.0));
      assert(SCIPisEQ(scip, SCIPvarGetUbLocal(vars[watchedvar2]), 1.0));
      *mustcheck = TRUE;
      consdata->changed = FALSE;
      return SCIP_OKAY;
   }

   /* make sure, that if there is a valid watched variable, it is the first one */
   if( watchedvar1 == -1 )
   {
      watchedvar1 = watchedvar2;
      watchedvar2 = -1;
   }

   /* we have to search new watched variables: loop through variables until two unfixed variables are found */
   for( v = 0; v < nvars && watchedvar2 == -1; ++v )
   {
      /* don't use the same variable in both watch slots */
      if( v == watchedvar1 )
         continue;

      /* check, if the variable is fixed */
      lb = SCIPvarGetLbLocal(vars[v]);
      ub = SCIPvarGetUbLocal(vars[v]);
      if( SCIPisEQ(scip, lb, ub) )
      {
         if( SCIPisEQ(scip, lb, 1.0) )
         {
            /* the variable is fixed to one, making the constraint redundant;
             * remember the variable and disable the constraint
             */
            debugMessage("disabling constraint <%s> (found variable fixed to 1.0)\n", SCIPconsGetName(cons));
            consdata->watchedfeasvar = v;
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
            return SCIP_OKAY;
         }
      }
      else
      {
         /* the variable is unfixed and can be used as watched variable */
         if( watchedvar1 == -1 )
            watchedvar1 = v;
         else
            watchedvar2 = v;
      }
   }
   assert(watchedvar1 >= 0 || watchedvar2 == -1);

   if( watchedvar1 == -1 )
   {
      /* there is no unfixed variable left -> the constraint is infeasible
       *  - a modifiable constraint must be added as a cut and further pricing must be performed in the LP solving loop
       *  - an unmodifiable constraint is infeasible and the node can be cut off
       */
      assert(watchedvar2 == -1);
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      if( SCIPconsIsModifiable(cons) )
         *addcut = TRUE;
      else
      {
         if( SCIPconsIsGlobal(cons) )
         {
            /* use conflict analysis to get a conflict clause out of the conflicting assignment */
            CHECK_OKAY( analyzeConflict(scip, consdata) );
         }

         /* mark the node to be cut off */
         *cutoff = TRUE;
      }
   }
   else if( watchedvar2 == -1 )
   {
      /* there is only one unfixed variable:
       * - a modifiable constraint must be checked manually
       * - an unmodifiable constraint is feasible and can be disabled after the remaining variable is fixed to one
       */
      assert(0 <= watchedvar1 && watchedvar1 < nvars);
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(vars[watchedvar1]), 0.0));
      assert(SCIPisEQ(scip, SCIPvarGetUbLocal(vars[watchedvar1]), 1.0));
      if( SCIPconsIsModifiable(cons) )
         *mustcheck = TRUE;
      else
      {
         /* fixed remaining variable to one and disable constraint */
         debugMessage("single-literal constraint <%s> (fix <%s> to 1.0) at depth %d\n", 
            SCIPconsGetName(cons), SCIPvarGetName(vars[watchedvar1]), SCIPgetActDepth(scip));
         CHECK_OKAY( SCIPinferBinVar(scip, vars[watchedvar1], TRUE, cons) ); /* provide cons for conflict analysis */
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
         *reduceddom = TRUE;
      }
   }
   else
   {
      /* there are at least two unfixed variables -> the constraint must be checked manually */
      assert(watchedvar1 != watchedvar2);
      assert(0 <= watchedvar1 && watchedvar1 < nvars);
      assert(0 <= watchedvar2 && watchedvar2 < nvars);
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(vars[watchedvar1]), 0.0));
      assert(SCIPisEQ(scip, SCIPvarGetUbLocal(vars[watchedvar1]), 1.0));
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(vars[watchedvar2]), 0.0));
      assert(SCIPisEQ(scip, SCIPvarGetUbLocal(vars[watchedvar2]), 1.0));
      *mustcheck = TRUE;

      /* drop events on old watched variables */
      if( consdata->watchedvar1 != -1
         && consdata->watchedvar1 != watchedvar1 && consdata->watchedvar1 != watchedvar2 )
      {
         CHECK_OKAY( consdataDropEvent(scip, consdata, eventhdlr, consdata->watchedvar1) );
      }
      if( consdata->watchedvar2 != -1
         && consdata->watchedvar2 != watchedvar1 && consdata->watchedvar2 != watchedvar2 )
      {
         CHECK_OKAY( consdataDropEvent(scip, consdata, eventhdlr, consdata->watchedvar2) );
      }

      /* catch events on new watched variables */
      if( watchedvar1 != consdata->watchedvar1 && watchedvar1 != consdata->watchedvar2 )
      {
         CHECK_OKAY( consdataCatchEvent(scip, consdata, eventhdlr, watchedvar1) );
      }
      if( watchedvar2 != consdata->watchedvar1 && watchedvar2 != consdata->watchedvar2 )
      {
         CHECK_OKAY( consdataCatchEvent(scip, consdata, eventhdlr, watchedvar2) );
      }

      /* remember the two watched variables for next time */
      consdata->watchedvar1 = watchedvar1;
      consdata->watchedvar2 = watchedvar2;
      consdata->changed = FALSE;

      CHECK_OKAY( SCIPincConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, returns TRUE iff constraint is feasible */
static
Bool checkCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< logic or constraint to be checked */
   SOL*             sol                 /**< primal CIP solution */
   )
{
   VAR** vars;
   Real solval;
   Real sum;
   int nvars;
   int v;
   
   vars = consdata->vars;
   nvars = consdata->nvars;

   /* check the watched solution variable, if it already makes the constraint feasible */
   if( consdata->watchedsolvar != -1 )
   {
      assert(consdata->watchedsolvar < nvars);
      solval = SCIPgetSolVal(scip, sol, vars[consdata->watchedsolvar]);
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
      if( SCIPisFeasEQ(scip, solval, 1.0) )
         return TRUE;
   }

   /* calculate the constraint's activity */
   sum = 0.0;
   solval = 0.0;
   assert(SCIPfeastol(scip) < 0.1); /* to make the comparison against 1.1 working */
   for( v = 0; v < nvars && sum < 1.0; ++v )
   {
      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
      solval = SCIPgetSolVal(scip, sol, vars[v]);
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
      sum += solval;
   }

   /* if the last compared variable has solution value of 1.0, use it as new watched solution variable */
   if( SCIPisGE(scip, solval, 1.0) )
   {
      assert(0 < v && v <= nvars);
      consdata->watchedsolvar = v-1;
   }

   return SCIPisFeasGE(scip, sum, 1.0);
}

/** creates an LP row in a logic or constraint data object */
static
RETCODE createRow(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< logic or constraint */
   )
{
   CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), 1.0, SCIPinfinity(scip),
                  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   
   for( v = 0; v < consdata->nvars; ++v )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->row, consdata->vars[v], 1.0) );
   }

   return SCIP_OKAY;
}

/** adds logic or constraint as cut to the LP */
static
RETCODE addCut(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint */
   Real             violation           /**< absolute violation of the constraint */
   )
{
   CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   if( consdata->row == NULL )
   {
      /* convert logic or constraint data into LP row */
      CHECK_OKAY( createRow(scip, cons) );
   }
   assert(consdata->row != NULL);
   assert(!SCIProwIsInLP(consdata->row));
            
   /* insert LP row as cut */
   CHECK_OKAY( SCIPaddCut(scip, consdata->row, violation/(SCIProwGetNNonz(consdata->row)+1)) );

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
RETCODE separateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint to be separated */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            separated,          /**< pointer to store TRUE, if a cut was found */
   Bool*            reduceddom          /**< pointer to store TRUE, if a domain reduction was found */
   )
{
   CONSDATA* consdata;
   Bool addcut;
   Bool mustcheck;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(separated != NULL);
   assert(reduceddom != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);

   /* skip constraints already in the LP */
   if( consdata->row != NULL && SCIProwIsInLP(consdata->row) )
      return SCIP_OKAY;

   /* update and check the watched variables */
   CHECK_OKAY( processWatchedVars(scip, cons, eventhdlr, cutoff, reduceddom, &addcut, &mustcheck) );

   if( mustcheck )
   {
      assert(!addcut);

      /* variable's fixings didn't give us any information -> we have to check the constraint */
      if( consdata->row != NULL )
      {
         Real feasibility;

         assert(!SCIProwIsInLP(consdata->row));
         feasibility = SCIPgetRowLPFeasibility(scip, consdata->row);
         addcut = !SCIPisFeasible(scip, feasibility);
      }
      else
         addcut = !checkCons(scip, consdata, NULL);
   }

   if( addcut )
   {
      /* insert LP row as cut */
      CHECK_OKAY( addCut(scip, cons, 1.0) );
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      *separated = TRUE;
   }
   else
   {
      /* constraint was feasible -> increase age */
      CHECK_OKAY( SCIPincConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** enforces the pseudo solution on the given constraint */
static
RETCODE enforcePseudo(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< logic or constraint to be separated */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            infeasible,         /**< pointer to store TRUE, if the constraint was infeasible */
   Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   Bool*            solvelp             /**< pointer to store TRUE, if the LP has to be solved */
   )
{
   Bool addcut;
   Bool mustcheck;

   assert(!SCIPhasActnodeLP(scip));
   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(reduceddom != NULL);
   assert(solvelp != NULL);

   /* update and check the watched variables */
   CHECK_OKAY( processWatchedVars(scip, cons, eventhdlr, cutoff, reduceddom, &addcut, &mustcheck) );

   if( mustcheck )
   {
      CONSDATA* consdata;

      assert(!addcut);

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( checkCons(scip, consdata, NULL) )
      {
         /* constraint was feasible -> increase age */
         CHECK_OKAY( SCIPincConsAge(scip, cons) );
      }
      else
      {
         /* constraint was infeasible -> reset age */
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *infeasible = TRUE;
      }
   }
   else if( addcut )
   {
      /* a cut must be added to the LP -> we have to solve the LP immediately */
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      *solvelp = TRUE;
   }
   else
   {
      /* constraint was feasible -> increase age */
      CHECK_OKAY( SCIPincConsAge(scip, cons) );
   }
   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
DECL_CONSFREE(consFreeLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   CHECK_OKAY( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called when problem solving starts) */
#define consInitLogicor NULL


/** deinitialization method of constraint handler (called when problem solving exits) */
#define consExitLogicor NULL


/** solving start notification method of constraint handler (called when presolving was finished) */
#define consSolstartLogicor NULL


/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* free LP row and logic or constraint */
   CHECK_OKAY( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransLogicor)
{  /*lint --e{715}*/
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   /*debugMessage("Trans method of logic or constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPstage(scip) == SCIP_STAGE_INITSOLVE);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */

   /* create constraint data for target constraint */
   CHECK_OKAY( consdataCreateTransformed(scip, &targetdata, sourcedata->nvars, sourcedata->vars) );

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                  SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                  SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                  SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpLogicor)
{  /*lint --e{715}*/
   int c;

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsInitial(conss[c]) )
      {
         CHECK_OKAY( addCut(scip, conss[c], 0.0) );
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler */
static
DECL_CONSSEPA(consSepaLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool cutoff;
   Bool separated;
   Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("separating %d/%d logic or constraints\n", nusefulconss, nconss);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* step 1: check all useful logic or constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
   }

   /* step 2: combine logic or constraints to get more cuts */
   todoMessage("further cuts of logic or constraints");

   /* step 3: if no cuts were found and we are in the root node, separate remaining constraints */
   if( SCIPgetActDepth(scip) == 0 )
   {
      for( c = nusefulconss; c < nconss && !cutoff && !separated && !reduceddom; ++c )
      {
         CHECK_OKAY( separateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
      }
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( separated )
      *result = SCIP_SEPARATED;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}

#ifdef BRANCHLP
/** if fractional variables exist, chooses a set S of them and branches on (i) x(S) >= 1, and (ii) x(S) >= 0 */
static
RETCODE branchLP(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< logic or constraint handler */
   RESULT*          result              /**< pointer to store the result SCIP_BRANCHED, if branching was applied */
   )
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* posvaruses;
   INTARRAY* negvaruses;
   VAR** lpcands;
   VAR** branchcands;
   VAR* var;
   Real* usescores;
   Real maxvarusefac;
   Real minvarusefac;
   Real actusescore;
   Real branchweight;
   Real solval;
   int varindex;
   int nlpcands;
   int nbranchcands;
   int nselcands;
   int posuse;
   int neguse;
   int i;
   int j;

   todoMessage("use a better logicor branching on LP solution");

   assert(conshdlr != NULL);
   assert(result != NULL);

   /* get fractional variables */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, NULL, NULL, &nlpcands) );
   if( nlpcands == 0 )
      return SCIP_OKAY;

   /* get temporary memory */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &branchcands, nlpcands) );
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &usescores, nlpcands) );
   
   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   posvaruses = conshdlrdata->posvaruses;
   negvaruses = conshdlrdata->negvaruses;
   maxvarusefac = conshdlrdata->maxvarusefac;
   minvarusefac = conshdlrdata->minvarusefac;
   assert(posvaruses != NULL);
   assert(negvaruses != NULL);

   /* sort fractional variables by weighted number of uses in enabled logic or constraints */
   nbranchcands = 0;
   for( i = 0; i < nlpcands; ++i )
   {
      var = lpcands[i];
      varindex = SCIPvarGetIndex(var);
      posuse = SCIPgetIntarrayVal(scip, posvaruses, varindex);
      neguse = SCIPgetIntarrayVal(scip, negvaruses, varindex);
      if( posuse + neguse > 0 )
      {
         actusescore = maxvarusefac * MAX(posuse, neguse) + minvarusefac * MIN(posuse, neguse);
         for( j = nbranchcands; j > 0 && actusescore > usescores[j-1]; --j )
         {
            branchcands[j] = branchcands[j-1];
            usescores[j] = usescores[j-1];
         }
         assert(0 <= j && j <= nbranchcands);

         /* choose between normal and negated variable in branching, such that setting selected variable
          * to zero is hopefully leading to a feasible solution
          */
         if( posuse > neguse )
         {
            CHECK_OKAY( SCIPgetNegatedVar(scip, var, &branchcands[j]) );
         }
         else
            branchcands[j] = var;
         usescores[j] = actusescore;
         nbranchcands++;
      }
   }
   assert(nbranchcands <= nlpcands);

   /* if none of the fractional variables is member of a logic or constraint,
    * we are not responsible for doing the branching
    */
   if( nbranchcands > 0 )
   {
      /* select the first variables from the sorted candidate list, until MAXBRANCHWEIGHT is reached;
       * then choose one less
       */
      branchweight = 0.0;
      solval = 0.0;
      for( nselcands = 0; nselcands < nbranchcands && branchweight <= MAXBRANCHWEIGHT; ++nselcands )
      {
         solval = SCIPgetVarSol(scip, branchcands[nselcands]);
         assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
         branchweight += solval;
      }
      assert(nselcands > 0);
      nselcands--;
      branchweight -= solval;

      /* check, if we accumulated at least MIN and at most MAXBRANCHWEIGHT weight */
      if( MINBRANCHWEIGHT <= branchweight && branchweight <= MAXBRANCHWEIGHT )
      {
         NODE* node;

         /* perform the binary set branching on the selected variables */
         assert(nselcands <= nlpcands);
         
         /* create left child, fix x_i = 0 for all i \in S */
         CHECK_OKAY( SCIPcreateChild(scip, &node) );
         for( i = 0; i < nselcands; ++i )
         {
            CHECK_OKAY( SCIPchgVarUbNode(scip, node, branchcands[i], 0.0) );
         }

         /* create right child: add constraint x(S) >= 1 */
         CHECK_OKAY( SCIPcreateChild(scip, &node) );
         if( nselcands == 1 )
         {
            /* only one candidate selected: fix it to 1.0 */
            debugMessage("fixing variable <%s> to 1.0 in right child node\n", SCIPvarGetName(branchcands[0]));
            CHECK_OKAY( SCIPchgVarLbNode(scip, node, branchcands[0], 1.0) );
         }
         else
         {
            CONS* newcons;
            char name[MAXSTRLEN];
         
            /* add logic or constraint x(S) >= 1 */
            sprintf(name, "LOB%lld", SCIPgetNodenum(scip));

            CHECK_OKAY( SCIPcreateConsLogicor(scip, &newcons, name, nselcands, branchcands,
                           FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE) );
            CHECK_OKAY( SCIPaddConsNode(scip, node, newcons) );
            CHECK_OKAY( SCIPreleaseCons(scip, &newcons) );
         }
      
         *result = SCIP_BRANCHED;
         
#ifdef DEBUG
         printf("logic or LP branching: nselcands=%d/%d, weight(S)=%g, S={", nselcands, nlpcands, branchweight);
         for( i = 0; i < nselcands; ++i )
            printf(" %s[%g, %g]", SCIPvarGetName(branchcands[i]), usescores[i], SCIPgetSolVal(scip, NULL, branchcands[i]));
         printf(" }\n");
#endif
      }
   }

   /* free temporary memory */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &usescores) );
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &branchcands) );

   return SCIP_OKAY;
}
#endif

/** if unfixed variables exist, chooses a set S of them and creates |S|+1 child nodes:
 *   - for each variable i from S, create child node with x_0 = ... = x_i-1 = 0, x_i = 1
 *   - create an additional child node x_0 = ... = x_n-1 = 0
 */
static
RETCODE branchPseudo(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< logic or constraint handler */
   RESULT*          result              /**< pointer to store the result SCIP_BRANCHED, if branching was applied */
   )
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* posvaruses;
   INTARRAY* negvaruses;
   VAR** pseudocands;
   VAR** branchcands;
   VAR* var;
   NODE* node;
   Real* usescores;
   Real maxvarusefac;
   Real minvarusefac;
   Real actusescore;
   int varindex;
   int npseudocands;
   int maxnbranchcands;
   int nbranchcands;
   int posuse;
   int neguse;
   int i;
   int j;

   todoMessage("use a better logic or branching on pseudo solution");

   assert(conshdlr != NULL);
   assert(result != NULL);

   /* get fractional variables */
   CHECK_OKAY( SCIPgetPseudoBranchCands(scip, &pseudocands, &npseudocands) );
   if( npseudocands == 0 )
      return SCIP_OKAY;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   posvaruses = conshdlrdata->posvaruses;
   negvaruses = conshdlrdata->negvaruses;
   maxvarusefac = conshdlrdata->maxvarusefac;
   minvarusefac = conshdlrdata->minvarusefac;
   maxnbranchcands = conshdlrdata->npseudobranches-1;
   assert(posvaruses != NULL);
   assert(negvaruses != NULL);
   assert(maxnbranchcands >= 1);

   /* get temporary memory */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &branchcands, maxnbranchcands) );
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &usescores, maxnbranchcands) );
   
   /* sort unfixed variables by number of uses in enabled logic or constraints */
   nbranchcands = 0;
   for( i = 0; i < npseudocands; ++i )
   {
      var = pseudocands[i];
      varindex = SCIPvarGetIndex(var);
      posuse = SCIPgetIntarrayVal(scip, posvaruses, varindex);
      neguse = SCIPgetIntarrayVal(scip, negvaruses, varindex);
      if( posuse + neguse > 0 )
      {
         actusescore = maxvarusefac * MAX(posuse, neguse) + minvarusefac * MIN(posuse, neguse);
         if( nbranchcands < maxnbranchcands || actusescore > usescores[nbranchcands-1] )
         {
            for( j = MIN(nbranchcands, maxnbranchcands-1); j > 0 && actusescore > usescores[j-1]; --j )
            {
               branchcands[j] = branchcands[j-1];
               usescores[j] = usescores[j-1];
            }
            assert(0 <= j && j <= nbranchcands && j < maxnbranchcands);

            /* choose between normal and negated variable in branching, such that setting selected variable
             * to zero is hopefully leading to a feasible solution
             */
            if( posuse > neguse )
            {
               CHECK_OKAY( SCIPgetNegatedVar(scip, var, &branchcands[j]) );
            }
            else
               branchcands[j] = var;
            usescores[j] = actusescore;
            if( nbranchcands < maxnbranchcands )
               nbranchcands++;
         }
      }
   }
   assert(nbranchcands <= maxnbranchcands);

   /* if none of the unfixed variables is member of a logic or constraint,
    * we are not responsible for doing the branching
    */
   if( nbranchcands > 0 )
   {
      /* branch on the first part of the sorted candidates:
       * - for each of these variables i, create a child node x_0 = ... = x_i-1 = 0, x_i = 1
       * - create an additional child node x_0 = ... = x_n-1 = 0
       */

      /* create child with x_0 = ... = x_n-1 = 0 */
      CHECK_OKAY( SCIPcreateChild(scip, &node) );
      for( i = 0; i < nbranchcands; ++i )
      {
         CHECK_OKAY( SCIPchgVarUbNode(scip, node, branchcands[i], 0.0) );
      }

      /* create children with x_0 = ... = x_i-1 = 0, x_i = 1, i = n-1,...,0 */
      for( i = nbranchcands-1; i >= 0; --i )
      {            
         CHECK_OKAY( SCIPcreateChild(scip, &node) );
         for( j = 0; j < i; ++j )
         {
            CHECK_OKAY( SCIPchgVarUbNode(scip, node, branchcands[j], 0.0) );
         }
         CHECK_OKAY( SCIPchgVarLbNode(scip, node, branchcands[i], 1.0) );
      }

      *result = SCIP_BRANCHED;

#ifdef DEBUG
      printf("logic or pseudo branching: nbranchcands=%d/%d, S={", nbranchcands, maxnbranchcands);
      for( i = 0; i < nbranchcands; ++i )
         printf(" %s [%g]", SCIPvarGetName(branchcands[i]), usescores[i]);
      printf(" }\n");
#endif
   }

   /* free temporary memory */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &usescores) );
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &branchcands) );

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool cutoff;
   Bool separated;
   Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("LP enforcing %d logic or constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* step 1: check all useful logic or constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
   }

   /* step 2: check all obsolete logic or constraints for feasibility */
   for( c = nusefulconss; c < nconss && !cutoff && !separated && !reduceddom; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
   }

#ifdef BRANCHLP
   /* step 3: if solution is not integral, choose a variable set to branch on */
   if( !cutoff && !separated && !reduceddom )
   {
      CHECK_OKAY( branchLP(scip, conshdlr, result) );
      if( *result != SCIP_FEASIBLE )
         return SCIP_OKAY;
   }
#endif

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( separated )
      *result = SCIP_SEPARATED;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool cutoff;
   Bool infeasible;
   Bool reduceddom;
   Bool solvelp;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   /* if the solution is infeasible anyway due to objective value, skip the constraint processing and branch directly */
   if( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      CHECK_OKAY( branchPseudo(scip, conshdlr, result) );
      return SCIP_OKAY;
   }

   debugMessage("pseudo enforcing %d logic or constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   infeasible = FALSE;
   reduceddom = FALSE;
   solvelp = FALSE;

   /* check all logic or constraints for feasibility */
   for( c = 0; c < nconss && !cutoff && !reduceddom && !solvelp; ++c )
   {
      CHECK_OKAY( enforcePseudo(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &infeasible, &reduceddom, &solvelp) );
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else if( solvelp )
      *result = SCIP_SOLVELP;
   else if( infeasible )
   {
      *result = SCIP_INFEASIBLE;
      
      /* at least one constraint is violated by pseudo solution and we didn't find a better way to resolve this:
       * -> branch on pseudo solution
       */
      CHECK_OKAY( branchPseudo(scip, conshdlr, result) );
   }
   
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckLogicor)
{  /*lint --e{715}*/
   CONS* cons;
   CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* check all logic or constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
      {
         if( checkCons(scip, consdata, sol) )
         {
            /* constraint was feasible -> increase age */
            CHECK_OKAY( SCIPincConsAge(scip, cons) );
         }
         else
         {
            /* constraint is violated */
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            *result = SCIP_INFEASIBLE;
            return SCIP_OKAY;
         }
      }
   }
   
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
DECL_CONSPROP(consPropLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool cutoff;
   Bool reduceddom;
   Bool addcut;
   Bool mustcheck;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   reduceddom = FALSE;

   /* step 1: propagate all useful logic or constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      CHECK_OKAY( processWatchedVars(scip, conss[c], conshdlrdata->eventhdlr,
                     &cutoff, &reduceddom, &addcut, &mustcheck) );
   }

   /* step 2: if no reduction was found, propagate all obsolete logic or constraints */
   if( !cutoff && !reduceddom )
   {
      for( c = nusefulconss; c < nconss && !cutoff; ++c )
      {
         CHECK_OKAY( processWatchedVars(scip, conss[c], conshdlrdata->eventhdlr,
                        &cutoff, &reduceddom, &addcut, &mustcheck) );
      }
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
DECL_CONSPRESOL(consPresolLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONS* cons;
   CONSDATA* consdata;
   Bool infeasible;
   Bool redundant;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      debugMessage("presolving logic or constraint <%s>\n", SCIPconsGetName(cons));

      /* remove all variables that are fixed to zero, check redundancy due to fixed-to-one variable */
      CHECK_OKAY( applyFixings(scip, cons, conshdlrdata->eventhdlr, &redundant) );

      if( redundant )
      {
         debugMessage("logic or constraint <%s> is redundant\n", SCIPconsGetName(cons));
         CHECK_OKAY( SCIPdelCons(scip, cons) );
         (*ndelconss)++;
         *result = SCIP_SUCCESS;
         continue;
      }
      else if( !SCIPconsIsModifiable(cons) )
      {
         /* if unmodifiable constraint has no variables, it is infeasible,
          * if unmodifiable constraint has only one variable, this one can be fixed and the constraint deleted
          */
         if( consdata->nvars == 0 )
         {
            debugMessage("logic or constraint <%s> is infeasible\n", SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if( consdata->nvars == 1 )
         {
            debugMessage("logic or constraint <%s> has only one variable not fixed to 0.0\n",
               SCIPconsGetName(cons));
            
            assert(consdata->vars != NULL);
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
            
            CHECK_OKAY( SCIPfixVar(scip, consdata->vars[0], 1.0, &infeasible) );
            if( infeasible )
            {
               debugMessage(" -> infeasible fixing\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            CHECK_OKAY( SCIPdelCons(scip, cons) );
            (*nfixedvars)++;
            (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }
      }
   }
   
   return SCIP_OKAY;
}


/** conflict variable resolving method of constraint handler */
static
DECL_CONSRESCVAR(consRescvarLogicor)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   Bool infervarfound;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(infervar != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("conflict resolving method of logic or constraint handler\n");

   /* the only deductions are variables infered to 1.0 on logic or constraints where all other variables
    * are assigned to zero
    */
   assert(SCIPvarGetLbLocal(infervar) > 0.5); /* the inference variable must be assigned to one */

   infervarfound = FALSE;
   for( v = 0; v < consdata->nvars; ++v )
   {
      if( consdata->vars[v] != infervar )
      {
         assert(SCIPvarGetUbLocal(consdata->vars[v]) < 0.5); /* the reason variable must be assigned to zero */
         CHECK_OKAY( SCIPaddConflictVar(scip, consdata->vars[v]) );
      }
      else
      {
         assert(!infervarfound);
         infervarfound = TRUE;
      }
   }
   assert(infervarfound);

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockLogicor)
{  /*lint --e{715}*/
   consdataLockAllRoundings(SCIPconsGetData(cons), nlockspos, nlocksneg);

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockLogicor)
{  /*lint --e{715}*/
   consdataUnlockAllRoundings(SCIPconsGetData(cons), nunlockspos, nunlocksneg);

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
DECL_CONSENABLE(consActiveLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("activation information method of logic or constraint handler\n");

   /* increase the number of uses for each variable in the constraint */
   for( v = 0; v < consdata->nvars; ++v )
   {
      CHECK_OKAY( conshdlrdataIncVaruses(scip, conshdlrdata, consdata->vars[v]) );
   }

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
DECL_CONSDISABLE(consDeactiveLogicor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("deactivation information method of logic or constraint handler\n");

   /* decrease the number of uses for each variable in the constraint */
   for( v = 0; v < consdata->nvars; ++v )
   {
      CHECK_OKAY( conshdlrdataDecVaruses(scip, conshdlrdata, consdata->vars[v]) );
   }

   return SCIP_OKAY;
}

/** constraint enabling notification method of constraint handler */
#define consEnableLogicor NULL


/** constraint disabling notification method of constraint handler */
#define consDisableLogicor NULL




/*
 * upgrading of linear constraints
 */

/** creates and captures a normalized (with all coefficients +1) logic or constraint */
static
RETCODE createNormalizedLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients (+1.0 or -1.0) */
   int              mult,               /**< multiplier on the coefficients(+1 or -1) */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   VAR** transvars;
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);
   assert(mult == +1 || mult == -1);

   /* get temporary memory */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &transvars, nvars) );

   /* negate positive or negative variables */
   for( v = 0; v < nvars; ++v )
   {
      if( mult * vals[v] > 0.0 )
         transvars[v] = vars[v];
      else
      {
         CHECK_OKAY( SCIPgetNegatedVar(scip, vars[v], &transvars[v]) );
      }
      assert(transvars[v] != NULL);
   }

   /* create the constraint */
   CHECK_OKAY( SCIPcreateConsLogicor(scip, cons, name, nvars, transvars,
                  initial, separate, enforce, check, propagate, local, modifiable, removeable) );

   /* release temporary memory */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &transvars) );

   return SCIP_OKAY;
}

static
DECL_LINCONSUPGD(linconsUpgdLogicor)
{  /*lint --e{715}*/
   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to logic or constraint
    * - logic or constraints consist only of binary variables with a
    *   coefficient of +1.0 or -1.0 (variables with -1.0 coefficients can be negated):
    *        lhs     <= x1 + ... + xp - y1 - ... - yn <= rhs
    * - negating all variables y = (1-Y) with negative coefficients gives:
    *        lhs + n <= x1 + ... + xp + Y1 + ... + Yn <= rhs + n
    * - negating all variables x = (1-X) with positive coefficients and multiplying with -1 gives:
    *        p - rhs <= X1 + ... + Xp + y1 + ... + yn <= p - lhs
    * - logic or constraints have left hand side of +1.0, and right hand side of +infinity: x(S) >= 1.0
    *    -> without negations:  (lhs == 1 - n  and  rhs == +inf)  or  (lhs == -inf  and  rhs = p - 1)
    */
   if( nposbin + nnegbin == nvars && ncoeffspone + ncoeffsnone == nvars
      && ((SCIPisEQ(scip, lhs, 1.0 - ncoeffsnone) && SCIPisInfinity(scip, rhs))
         || (SCIPisInfinity(scip, -lhs) && SCIPisEQ(scip, rhs, ncoeffspone - 1.0))) )
   {
      int mult;

      debugMessage("upgrading constraint <%s> to logic or constraint\n", SCIPconsGetName(cons));
      
      /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
      mult = SCIPisInfinity(scip, rhs) ? +1 : -1;
      
      /* create the logic or constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      CHECK_OKAY( createNormalizedLogicor(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, mult,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecLogicor)
{  /*lint --e{715}*/
   CONSDATA* consdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   debugMessage("Exec method of bound tighten event handler for logic or constraints\n");

   consdata = (CONSDATA*)eventdata;
   assert(consdata != NULL);

   consdata->changed = TRUE;

   return SCIP_OKAY;
}




/*
 * Callback methods of conflict handler
 */

static
DECL_CONFLICTEXEC(conflictExecLogicor)
{  /*lint --e{715}*/
   CONS* cons;
   char consname[MAXSTRLEN];
   
   assert(conflicthdlr != NULL);
   assert(strcmp(SCIPconflicthdlrGetName(conflicthdlr), CONFLICTHDLR_NAME) == 0);
   assert(conflictvars != NULL || nconflictvars == 0);
   assert(result != NULL);

   /* don't already resolved conflicts */
   if( resolved )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* create a constraint out of the conflict set */
   sprintf(consname, "cf%d", SCIPgetNConss(scip));
   CHECK_OKAY( SCIPcreateConsLogicor(scip, &cons, consname, nconflictvars, conflictvars, 
                  FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE) );
   CHECK_OKAY( SCIPaddCons(scip, cons) );
   CHECK_OKAY( SCIPreleaseCons(scip, &cons) );

   *result = SCIP_CONSADDED;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for logic or constraints and includes it in SCIP */
RETCODE SCIPincludeConsHdlrLogicor(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound tighten events on watched variables */
   CHECK_OKAY( SCIPincludeEventHdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
                  NULL, NULL, NULL,
                  NULL, eventExecLogicor,
                  NULL) );

   /* create conflict handler for logic or constraints */
   CHECK_OKAY( SCIPincludeConflictHdlr(scip, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
                  NULL, NULL, NULL, conflictExecLogicor,
                  NULL) );

   /* create constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ,
                  CONSHDLR_NEEDSCONS,
                  consFreeLogicor, consInitLogicor, consExitLogicor, consSolstartLogicor,
                  consDeleteLogicor, consTransLogicor, 
                  consInitlpLogicor, consSepaLogicor, 
                  consEnfolpLogicor, consEnfopsLogicor, consCheckLogicor, 
                  consPropLogicor, consPresolLogicor, consRescvarLogicor,
                  consLockLogicor, consUnlockLogicor,
                  consActiveLogicor, consDeactiveLogicor,
                  consEnableLogicor, consDisableLogicor,
                  conshdlrdata) );

   /* include the linear constraint to logicor constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdLogicor, LINCONSUPGD_PRIORITY) );

   /* logic or constraint handler parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "conshdlr/logicor/npseudobranches", 
                  "number of children created in pseudo branching",
                  &conshdlrdata->npseudobranches, DEFAULT_NPSEUDOBRANCHES, 2, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "conshdlr/logicor/maxvarusefac", 
                  "branching factor to weigh maximum of positive and negative variable uses",
                  &conshdlrdata->maxvarusefac, DEFAULT_MAXVARUSEFAC, REAL_MIN, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "conshdlr/logicor/minvarusefac", 
                  "branching factor to weigh minimum of positive and negative variable uses",
                  &conshdlrdata->minvarusefac, DEFAULT_MINVARUSEFAC, REAL_MIN, REAL_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}


/** creates and captures a logic or constraint */
RETCODE SCIPcreateConsLogicor(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   assert(scip != NULL);

   /* find the logicor constraint handler */
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("logic or constraint handler not found");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   if( SCIPstage(scip) == SCIP_STAGE_PROBLEM )
   {
      /* create constraint in original problem */
      CHECK_OKAY( consdataCreate(scip, &consdata, nvars, vars) );
   }
   else
   {
      /* create constraint in transformed problem */
      CHECK_OKAY( consdataCreateTransformed(scip, &consdata, nvars, vars) );
   }

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
                  local, modifiable, removeable) );

   return SCIP_OKAY;
}
