/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (c) 2002-2004 Tobias Achterberg                              */
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
#pragma ident "@(#) $Id: cons_and.c,v 1.30 2004/08/13 09:16:40 bzfpfend Exp $"

/**@file   cons_and.c
 * @brief  constraint handler for and constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cons_and.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "and"
#define CONSHDLR_DESC          "constraint handler for and constraints: r = and(x1, ..., xn)"
#define CONSHDLR_SEPAPRIORITY   +850000
#define CONSHDLR_ENFOPRIORITY   -850000
#define CONSHDLR_CHECKPRIORITY  -850000
#define CONSHDLR_SEPAFREQ             1
#define CONSHDLR_PROPFREQ             1
#define CONSHDLR_EAGERFREQ          100
#define CONSHDLR_MAXPREROUNDS        -1
#define CONSHDLR_NEEDSCONS         TRUE

#define EVENTHDLR_NAME         "and"
#define EVENTHDLR_DESC         "bound tighten event handler for and constraints"




/*
 * Data structures
 */

/** constraint data for and constraints */
struct ConsData
{
   VAR**            vars;               /**< variables in the and operation */
   VAR*             resvar;             /**< resultant variable */
   ROW**            rows;               /**< rows for linear relaxation of and constraint */
   int              nvars;              /**< number of variables in and operation */
   int              varssize;           /**< size of vars array */
   int              rowssize;           /**< size of rows array */
   int              watchedvar1;        /**< position of first watched operator variable */
   int              watchedvar2;        /**< position of second watched operator variable */
   unsigned int     propagated:1;       /**< is constraint already preprocessed/propagated? */
   unsigned int     nofixedzero:1;      /**< is none of the opereator variables fixed to FALSE? */
};

/** constraint handler data */
struct ConshdlrData
{
   EVENTHDLR*       eventhdlr;          /**< event handler for bound tighten events on watched variables */
};




/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_1,                          /**< v_i = FALSE                                  =>  r   = FALSE          */
   PROPRULE_2,                          /**< r   = TRUE                                   =>  v_i = TRUE for all i */
   PROPRULE_3,                          /**< v_i = TRUE for all i                         =>  r   = TRUE           */
   PROPRULE_4,                          /**< r   = FALSE, v_i = TRUE for all i except j   =>  v_j = FALSE          */
   PROPRULE_INVALID                     /**< propagation was applied without a specific propagation rule */
};
typedef enum Proprule PROPRULE;




/*
 * Local methods
 */

/** creates constaint handler data */
static
RETCODE conshdlrdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, conshdlrdata) );

   /* get event handler for catching bound tighten events on watched variables */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      errorMessage("event handler for and constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   
   return SCIP_OKAY;
}

/** frees constraint handler data */
static
RETCODE conshdlrdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** gets number of LP rows needed for the LP relaxation of the constraint */
static
int consdataGetNRows(
   CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   return consdata->nvars + 1;
}

/** catches events for tightening lower bound of variable at given position */
static
RETCODE consdataCatchLbEvent(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< and constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   /* catch lower bound tighten events on variable */
   CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_LBTIGHTENED, eventhdlr,
         (EVENTDATA*)consdata) );
   
   return SCIP_OKAY;
}

/** catches events for tightening upper bound of variable at given position */
static
RETCODE consdataCatchUbEvent(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< and constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   /* catch upper bound tighten events on variable */
   CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_UBTIGHTENED, eventhdlr,
         (EVENTDATA*)consdata) );
   
   return SCIP_OKAY;
}

/** drops events for tightening lower bound of variable at given position */
static
RETCODE consdataDropLbEvent(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< and constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   /* drop lower bound events on variable */
   CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_LBTIGHTENED, eventhdlr,
         (EVENTDATA*)consdata) );

   return SCIP_OKAY;
}

/** drops events for tightening upper bound of variable at given position */
static
RETCODE consdataDropUbEvent(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< and constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   /* drop upper bound events on variable */
   CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_UBTIGHTENED, eventhdlr,
         (EVENTDATA*)consdata) );

   return SCIP_OKAY;
}

/** catches needed bound tightening events on all variables of constraint, except the special ones for watched variables */
static
RETCODE consdataCatchEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< and constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* catch tightening events for both bounds on resultant variable */
   CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->resvar, SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
         (EVENTDATA*)consdata) );

   /* catch tightening events for upper bound on operator variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      CHECK_OKAY( consdataCatchUbEvent(scip, consdata, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops needed bound tightening events on all variables of constraint, except the special ones for watched variables */
static
RETCODE consdataDropEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< and constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* drop tightening events for both bounds on resultant variable */
   CHECK_OKAY( SCIPdropVarEvent(scip, consdata->resvar, SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
         (EVENTDATA*)consdata) );

   /* drop tightening events for upper bound on operator variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      CHECK_OKAY( consdataDropUbEvent(scip, consdata, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** locks rounding for variable in transformed and constraint */
static
void lockRounding(
   VAR*             var,                /**< problem variable */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   /* forbid rounding of variable */
   SCIPvarLock(var, nlockspos + nlocksneg, nlockspos + nlocksneg);
}

/** unlocks rounding for variable in transformed and constraint */
static
void unlockRounding(
   VAR*             var,                /**< problem variable */
   int              nunlockspos,        /**< decrease in number of rounding locks for constraint */
   int              nunlocksneg         /**< decrease in number of rounding locks for constraint's negation */
   )
{
   /* allow rounding of variable */
   SCIPvarUnlock(var, nunlockspos + nunlocksneg, nunlockspos + nunlocksneg);
}

/** locks rounding for all variables in transformed and constraint */
static
void consdataLockAllRoundings(
   CONSDATA*        consdata,           /**< and constraint data */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   int i;

   assert(consdata != NULL);

   /* lock resultant variable */
   lockRounding(consdata->resvar, nlockspos, nlocksneg);

   /* lock all operand variables */
   for( i = 0; i < consdata->nvars; ++i )
      lockRounding(consdata->vars[i], nlockspos, nlocksneg);
}

/** unlocks rounding for all variables in transformed and constraint */
static
void consdataUnlockAllRoundings(
   CONSDATA*        consdata,           /**< and constraint data */
   int              nunlockspos,        /**< decrease in number of rounding locks for constraint */
   int              nunlocksneg         /**< decrease in number of rounding locks for constraint's negation */
   )
{
   int i;

   assert(consdata != NULL);

   /* unlock resultant variable */
   unlockRounding(consdata->resvar, nunlockspos, nunlocksneg);

   /* unlock all operand variables */
   for( i = 0; i < consdata->nvars; ++i )
      unlockRounding(consdata->vars[i], nunlockspos, nunlocksneg);
}

/** stores the given variable numbers as watched variables, and updates the event processing */
static
RETCODE consdataSwitchWatchedvars(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< and constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              watchedvar1,        /**< new first watched variable */
   int              watchedvar2         /**< new second watched variable */
   )
{
   assert(consdata != NULL);
   assert(watchedvar1 == -1 || watchedvar1 != watchedvar2);
   assert(watchedvar1 != -1 || watchedvar2 == -1);
   assert(watchedvar1 == -1 || (0 <= watchedvar1 && watchedvar1 < consdata->nvars));
   assert(watchedvar2 == -1 || (0 <= watchedvar2 && watchedvar2 < consdata->nvars));

   /* drop lower bound tighten events on old watched variables */
   if( consdata->watchedvar1 != -1 && consdata->watchedvar1 != watchedvar1 && consdata->watchedvar1 != watchedvar2 )
   {
      CHECK_OKAY( consdataDropLbEvent(scip, consdata, eventhdlr, consdata->watchedvar1) );
   }
   if( consdata->watchedvar2 != -1 && consdata->watchedvar2 != watchedvar1 && consdata->watchedvar2 != watchedvar2 )
   {
      CHECK_OKAY( consdataDropLbEvent(scip, consdata, eventhdlr, consdata->watchedvar2) );
   }

   /* catch lower bound tighten events on new watched variables */
   if( watchedvar1 != -1 && watchedvar1 != consdata->watchedvar1 && watchedvar1 != consdata->watchedvar2 )
   {
      CHECK_OKAY( consdataCatchLbEvent(scip, consdata, eventhdlr, watchedvar1) );
   }
   if( watchedvar2 != -1 && watchedvar2 != consdata->watchedvar1 && watchedvar2 != consdata->watchedvar2 )
   {
      CHECK_OKAY( consdataCatchLbEvent(scip, consdata, eventhdlr, watchedvar2) );
   }

   /* set the new watched variables */
   consdata->watchedvar1 = watchedvar1;
   consdata->watchedvar2 = watchedvar2;

   return SCIP_OKAY;
}

/** creates constraint data for and constraint */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              nvars,              /**< number of variables in the and operation */
   VAR**            vars,               /**< variables in and operation */
   VAR*             resvar              /**< resultant variable */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(resvar != NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );
   CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
   (*consdata)->resvar = resvar;
   (*consdata)->rows = NULL;
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->rowssize = 0;
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->propagated = FALSE;
   (*consdata)->nofixedzero = FALSE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      CHECK_OKAY( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
      CHECK_OKAY( SCIPgetTransformedVar(scip, (*consdata)->resvar, &(*consdata)->resvar) );

      /* catch needed bound tightening events on variables */
      CHECK_OKAY( consdataCatchEvents(scip, *consdata, eventhdlr) );
   }

   return SCIP_OKAY;
}

/** releases LP rows of constraint data and frees rows array */
static
RETCODE consdataFreeRows(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< constraint data */
   )
{
   int r;

   assert(consdata != NULL);

   if( consdata->rows != NULL )
   {
      for( r = 0; r < consdataGetNRows(consdata); ++r )
      {
         CHECK_OKAY( SCIPreleaseRow(scip, &consdata->rows[r]) );
      }
      SCIPfreeBlockMemoryArray(scip, &consdata->rows, consdata->rowssize);
   }

   return SCIP_OKAY;
}

/** frees constraint data for and constraint */
static
RETCODE consdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to the constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   if( SCIPisTransformed(scip) )
   {
      /* drop bound tighten events for watched variables */
      CHECK_OKAY( consdataSwitchWatchedvars(scip, *consdata, eventhdlr, -1, -1) );

      /* drop all other bound tightening events on variables */
      CHECK_OKAY( consdataDropEvents(scip, *consdata, eventhdlr) );
   }
   else
   {
      assert((*consdata)->watchedvar1 == -1);
      assert((*consdata)->watchedvar2 == -1);
   }

   /* release and free the rows */
   CHECK_OKAY( consdataFreeRows(scip, *consdata) );

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);
 
   return SCIP_OKAY;
}

/** prints and constraint to file stream */
static
void consdataPrint(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< and constraint data */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);

   if( file == NULL )
      file = stdout;

   /* print coefficients */
   fprintf(file, "<%s> == and(", SCIPvarGetName(consdata->resvar));
   for( v = 0; v < consdata->nvars; ++v )
   {
      if( v > 0 )
         fprintf(file, ", ");
      fprintf(file, "<%s>", SCIPvarGetName(consdata->vars[v]));
   }
   fprintf(file, ")\n");
}

/** deletes coefficient at given position from and constraint data */
static
RETCODE delCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< and constraint */
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
      if( consdata->watchedvar1 == pos )
      {
         CHECK_OKAY( consdataSwitchWatchedvars(scip, consdata, eventhdlr, consdata->watchedvar2, -1) );
      }
      if( consdata->watchedvar2 == pos )
      {
         CHECK_OKAY( consdataSwitchWatchedvars(scip, consdata, eventhdlr, consdata->watchedvar1, -1) );
      }
   }
   assert(pos != consdata->watchedvar1);
   assert(pos != consdata->watchedvar2);

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->nvars--;

   /* if the last variable (that moved) was watched, update the watched position */
   if( consdata->watchedvar1 == consdata->nvars )
      consdata->watchedvar1 = pos;
   if( consdata->watchedvar2 == consdata->nvars )
      consdata->watchedvar2 = pos;

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}

/** deletes all one-fixed variables */
static
RETCODE applyFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< and constraint */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   CONSDATA* consdata;
   VAR* var;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);

   v = 0;
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

      if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         CHECK_OKAY( delCoefPos(scip, cons, eventhdlr, v) );
      }
      else
         ++v;
   }

   debugMessage("after fixings: ");
   debug(consdataPrint(scip, consdata, NULL));

   return SCIP_OKAY;
}

/** creates LP rows corresponding to and constraint:
 *   - for each operator variable vi:  resvar - vi            <= 0
 *   - one additional row:             resvar - v1 - ... - vn >= n-1
 */
static 
RETCODE createRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to check */
   )
{
   CONSDATA* consdata;
   char rowname[MAXSTRLEN];
   int nvars;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows == NULL);

   nvars = consdata->nvars;

   /* get memory for rows */
   consdata->rowssize = consdataGetNRows(consdata);
   CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &consdata->rows, consdata->rowssize) );
   assert(consdata->rowssize == nvars+1);

   /* create operator rows */
   for( i = 0; i < nvars; ++i )
   {
      sprintf(rowname, "%s_%d", SCIPconsGetName(cons), i);
      CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->rows[i], rowname, -SCIPinfinity(scip), 0.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->rows[i], consdata->resvar, 1.0) );
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->rows[i], consdata->vars[i], -1.0) );
   }

   /* create additional row */
   sprintf(rowname, "%s_add", SCIPconsGetName(cons));
   CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->rows[nvars], rowname, -consdata->nvars + 1.0, SCIPinfinity(scip),
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   CHECK_OKAY( SCIPaddVarToRow(scip, consdata->rows[nvars], consdata->resvar, 1.0) );
   CHECK_OKAY( SCIPaddVarsToRowSameCoef(scip, consdata->rows[nvars], nvars, consdata->vars, -1.0) );

   return SCIP_OKAY;
}  

/** adds linear relaxation of and constraint to the LP */
static 
RETCODE addRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to check */
   )
{
   CONSDATA* consdata;
   int r;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rows == NULL )
   {
      CHECK_OKAY( createRelaxation(scip, cons) );
   }

   for( r = 0; r < consdataGetNRows(consdata); ++r )
   {
      CHECK_OKAY( SCIPaddCut(scip, consdata->rows[r], 1.0/(SCIProwGetNNonz(consdata->rows[r])+1)) );
   }

   return SCIP_OKAY;
}

/** checks and constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
RETCODE checkCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint to check */
   SOL*             sol,                /**< solution to check, NULL for current solution */
   Bool             checklprows,        /**< should LP rows be checked? */
   Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   CONSDATA* consdata;
   Bool mustcheck;
   int r;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   /* check, if we can skip this feasibility check, because all rows are in the LP and doesn't have to be checked */
   mustcheck = checklprows;
   mustcheck = mustcheck || (consdata->rows == NULL);
   if( !mustcheck )
   {
      assert(consdata->rows != NULL);
      for( r = 0; r < consdataGetNRows(consdata); ++r )
      {
         mustcheck = !SCIProwIsInLP(consdata->rows[r]);
         if( mustcheck )
            break;
      }         
   }

   /* check feasibility of constraint if necessary */
   if( mustcheck )
   {
      Real solval;
      int i;

      /* increase age of constraint; age is reset to zero, if a violation was found */
      CHECK_OKAY( SCIPincConsAge(scip, cons) );

      /* check, if all operator variables are TRUE */
      for( i = 0; i < consdata->nvars; ++i )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[i]);
         assert(SCIPisIntegral(scip, solval));
         if( solval < 0.5 )
            break;
      }

      /* if all operator variables are TRUE, the resultant has to be TRUE, otherwise, the resultant has to be FALSE */
      solval = SCIPgetSolVal(scip, sol, consdata->resvar);
      assert(SCIPisIntegral(scip, solval));

      if( (i == consdata->nvars) != (solval > 0.5) )
      {
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *violated = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** separates current LP solution */
static
RETCODE separateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint to check */
   Bool*            separated           /**< pointer to store whether a cut was found */
   )
{
   CONSDATA* consdata;
   Real feasibility;
   int r;

   assert(separated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *separated = FALSE;

   /* create all necessary rows for the linear relaxation */
   if( consdata->rows == NULL )
   {
      CHECK_OKAY( createRelaxation(scip, cons) );
   }
   assert(consdata->rows != NULL);

   /* test all rows for feasibility and add infeasible rows */
   for( r = 0; r < consdataGetNRows(consdata); ++r )
   {
      if( !SCIProwIsInLP(consdata->rows[r]) )
      {
         feasibility = SCIPgetRowLPFeasibility(scip, consdata->rows[r]);
         if( !SCIPisFeasible(scip, feasibility) )
         {
            CHECK_OKAY( SCIPaddCut(scip, consdata->rows[r], -feasibility/(SCIProwGetNNonz(consdata->rows[r])+1)) );
            *separated = TRUE;
         }
      }            
   }

   return SCIP_OKAY;
}

/** analyzes conflicting TRUE assignment to resultant of given constraint, and adds conflict clause to problem */
static
RETCODE analyzeConflictOne(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< and constraint that detected the conflict */
   int              falsepos            /**< position of operand that is fixed to FALSE */
   )
{
   CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarGetLbLocal(consdata->resvar) > 0.5);
   assert(0 <= falsepos && falsepos < consdata->nvars);
   assert(SCIPvarGetUbLocal(consdata->vars[falsepos]) < 0.5);

   /* initialize conflict analysis, and add resultant and single operand variable to conflict candidate queue */
   CHECK_OKAY( SCIPinitConflictAnalysis(scip) );
   CHECK_OKAY( SCIPaddConflictVar(scip, consdata->resvar) );
   CHECK_OKAY( SCIPaddConflictVar(scip, consdata->vars[falsepos]) );

   /* analyze the conflict */
   CHECK_OKAY( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** analyzes conflicting FALSE assignment to resultant of given constraint, and adds conflict clause to problem */
static
RETCODE analyzeConflictZero(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< or constraint that detected the conflict */
   )
{
   CONSDATA* consdata;
   int v;

   assert(!SCIPconsIsModifiable(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarGetUbLocal(consdata->resvar) < 0.5);

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   CHECK_OKAY( SCIPinitConflictAnalysis(scip) );
   CHECK_OKAY( SCIPaddConflictVar(scip, consdata->resvar) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(SCIPvarGetLbLocal(consdata->vars[v]) > 0.5);
      CHECK_OKAY( SCIPaddConflictVar(scip, consdata->vars[v]) );
   }

   /* analyze the conflict */
   CHECK_OKAY( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** propagates constraint with the following rules:
 *   (1) v_i = FALSE                                  =>  r   = FALSE
 *   (2) r   = TRUE                                   =>  v_i = TRUE for all i
 *   (3) v_i = TRUE for all i                         =>  r   = TRUE
 *   (4) r   = FALSE, v_i = TRUE for all i except j   =>  v_j = FALSE
 */
static
RETCODE propagateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< and constraint to be processed */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*             nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   CONSDATA* consdata;
   VAR* resvar;
   VAR** vars;
   int nvars;
   int watchedvar1;
   int watchedvar2;
   int i;
   Bool infeasible;
   Bool tightened;

   assert(cutoff != NULL);
   assert(nfixedvars != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   resvar = consdata->resvar;
   vars = consdata->vars;
   nvars = consdata->nvars;

   /* don't process the constraint, if none of the operator variables was fixed to FALSE, and if the watched variables
    * and the resultant weren't fixed to any value since last propagation call
    */
   if( consdata->propagated )
   {
      assert(consdata->nofixedzero);
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(resvar), 0.0));
      return SCIP_OKAY;
   }

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   CHECK_OKAY( SCIPincConsAge(scip, cons) );

   /* if one of the operator variables was fixed to FALSE, the resultant can be fixed to FALSE (rule (1)) */
   if( !consdata->nofixedzero )
   {
      for( i = 0; i < nvars && SCIPvarGetUbLocal(vars[i]) > 0.5; ++i ) /* search fixed operator */
      {}
      if( i < nvars )
      {
         debugMessage("constraint <%s>: operator var <%s> fixed to 0.0 -> fix resultant <%s> to 0.0\n",
            SCIPconsGetName(cons), SCIPvarGetName(vars[i]), SCIPvarGetName(resvar));
         CHECK_OKAY( SCIPinferBinVar(scip, resvar, FALSE, cons, PROPRULE_1, &infeasible, &tightened) );
         if( infeasible )
         {
            /* use conflict analysis to get a conflict clause out of the conflicting assignment */
            CHECK_OKAY( analyzeConflictOne(scip, cons, i) );
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else
         {
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
            if( tightened )
            {
               CHECK_OKAY( SCIPresetConsAge(scip, cons) );
               (*nfixedvars)++;
            }
         }

         return SCIP_OKAY;
      }
      else
         consdata->nofixedzero = TRUE;
   }
   assert(consdata->nofixedzero);

   /* if resultant is fixed to TRUE, all operator variables can be fixed to TRUE (rule (2)) */
   if( SCIPvarGetLbLocal(resvar) > 0.5 )
   {
      for( i = 0; i < nvars && !(*cutoff); ++i )
      {
         debugMessage("constraint <%s>: resultant var <%s> fixed to 1.0 -> fix operator var <%s> to 1.0\n",
            SCIPconsGetName(cons), SCIPvarGetName(resvar), SCIPvarGetName(vars[i]));
         CHECK_OKAY( SCIPinferBinVar(scip, vars[i], TRUE, cons, PROPRULE_2, &infeasible, &tightened) );
         if( infeasible )
         {
            /* use conflict analysis to get a conflict clause out of the conflicting assignment */
            CHECK_OKAY( analyzeConflictOne(scip, cons, i) );
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( tightened )
         {
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            (*nfixedvars)++;
         }
      }

      if( !(*cutoff) )
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      }

      return SCIP_OKAY;
   }

   /* rules (3) and (4) can only be applied, if we know all operator variables */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* rules (3) and (4) can not be applied, if we have at least two unfixed variables left;
    * that means, we only have to watch (i.e. capture events) of two variables, and switch to other variables
    * if these ones get fixed
    */
   watchedvar1 = consdata->watchedvar1;
   watchedvar2 = consdata->watchedvar2;

   /* check, if watched variables are still unfixed */
   if( watchedvar1 != -1 )
   {
      assert(SCIPvarGetUbLocal(vars[watchedvar1]) > 0.5); /* otherwise, rule (1) could be applied */
      if( SCIPvarGetLbLocal(vars[watchedvar1]) > 0.5 )
         watchedvar1 = -1;
   }
   if( watchedvar2 != -1 )
   {
      assert(SCIPvarGetUbLocal(vars[watchedvar2]) > 0.5); /* otherwise, rule (1) could be applied */
      if( SCIPvarGetLbLocal(vars[watchedvar2]) > 0.5 )
         watchedvar2 = -1;
   }
   
   /* if only one watched variable is still unfixed, make it the first one */
   if( watchedvar1 == -1 )
   {
      watchedvar1 = watchedvar2;
      watchedvar2 = -1;
   }
   assert(watchedvar1 != -1 || watchedvar2 == -1);

   /* if the watched variables are invalid (fixed), find new ones if existing */
   if( watchedvar2 == -1 )
   {
      for( i = 0; i < nvars; ++i )
      {
         assert(SCIPvarGetUbLocal(vars[i]) > 0.5); /* otherwise, rule (1) could be applied */
         if( SCIPvarGetLbLocal(vars[i]) < 0.5 )
         {
            if( watchedvar1 == -1 )
            {
               assert(watchedvar2 == -1);
               watchedvar1 = i;
            }
            else if( watchedvar1 != i )
            {
               watchedvar2 = i;
               break;
            }
         }
      }
   }
   assert(watchedvar1 != -1 || watchedvar2 == -1);

   /* if all variables are fixed to TRUE, the resultant can also be fixed to TRUE (rule (3)) */
   if( watchedvar1 == -1 )
   {
      assert(watchedvar2 == -1);
      
      debugMessage("constraint <%s>: all operator vars fixed to 1.0 -> fix resultant <%s> to 1.0\n",
         SCIPconsGetName(cons), SCIPvarGetName(resvar));
      CHECK_OKAY( SCIPinferBinVar(scip, resvar, TRUE, cons, PROPRULE_3, &infeasible, &tightened) );
      if( infeasible )
      {
         /* use conflict analysis to get a conflict clause out of the conflicting assignment */
         CHECK_OKAY( analyzeConflictZero(scip, cons) );
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;
      }
      else
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
         if( tightened )
         {
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            (*nfixedvars)++;
         }
      }

      return SCIP_OKAY;
   }

   /* if resultant is fixed to FALSE, and only one operator variable is not fixed to TRUE, this operator variable
    * can be fixed to FALSE (rule (4))
    */
   if( SCIPvarGetUbLocal(resvar) < 0.5 && watchedvar2 == -1 )
   {
      assert(watchedvar1 != -1);
      
      debugMessage("constraint <%s>: resultant <%s> fixed to 0.0, only one unfixed operand -> fix operand <%s> to 0.0\n",
         SCIPconsGetName(cons), SCIPvarGetName(resvar), SCIPvarGetName(vars[watchedvar1]));
      CHECK_OKAY( SCIPinferBinVar(scip, vars[watchedvar1], FALSE, cons, PROPRULE_4, &infeasible, &tightened) );
      if( infeasible )
      {
         /* use conflict analysis to get a conflict clause out of the conflicting assignment */
         CHECK_OKAY( analyzeConflictZero(scip, cons) );
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;
      }
      else
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
         if( tightened )
         {
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            (*nfixedvars)++;
         }
      }

      return SCIP_OKAY;
   }

   /* switch to the new watched variables */
   CHECK_OKAY( consdataSwitchWatchedvars(scip, consdata, eventhdlr, watchedvar1, watchedvar2) );

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** resolves a conflict on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) v_i = FALSE                                  =>  r   = FALSE
 *   (2) r   = TRUE                                   =>  v_i = TRUE for all i
 *   (3) v_i = TRUE for all i                         =>  r   = TRUE
 *   (4) r   = FALSE, v_i = TRUE for all i except j   =>  v_j = FALSE
 */
static
RETCODE resolveConflict(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< and constraint to be processed */
   VAR*             infervar,           /**< variable that was deduced */
   PROPRULE         proprule,           /**< propagation rule that deduced the value */
   RESULT*          result              /**< pointer to store the result of the conflict variable resolving call */
   )
{
   CONSDATA* consdata;
   VAR** vars;
   int nvars;
   int i;

   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   vars = consdata->vars;
   nvars = consdata->nvars;

   switch( proprule )
   {
   case PROPRULE_1:
      /* the resultant was infered to FALSE, because one operand variable was FALSE */
      assert(SCIPvarGetUbLocal(infervar) < 0.5);
      assert(infervar == consdata->resvar);
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPvarWasFixedEarlier(vars[i], consdata->resvar) && SCIPvarGetUbLocal(vars[i]) < 0.5 )
         {
            CHECK_OKAY( SCIPaddConflictVar(scip, vars[i]) );
            break;
         }
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_2:
      /* the operand variable was infered to TRUE, because the resultant was TRUE */
      assert(SCIPvarGetLbLocal(infervar) > 0.5);
      assert(SCIPvarWasFixedEarlier(consdata->resvar, infervar));
      assert(SCIPvarGetLbLocal(consdata->resvar) > 0.5);
      CHECK_OKAY( SCIPaddConflictVar(scip, consdata->resvar) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_3:
      /* the resultant was infered to TRUE, because all operand variables were TRUE */
      assert(SCIPvarGetLbLocal(infervar) > 0.5);
      assert(infervar == consdata->resvar);
      for( i = 0; i < nvars; ++i )
      {
         assert(SCIPvarWasFixedEarlier(vars[i], consdata->resvar));
         assert(SCIPvarGetLbLocal(vars[i]) > 0.5);
         CHECK_OKAY( SCIPaddConflictVar(scip, vars[i]) );
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_4:
      /* the operand variable was infered to FALSE, because the resultant was FALSE and all other operands were TRUE */
      assert(SCIPvarGetUbLocal(infervar) < 0.5);
      assert(SCIPvarWasFixedEarlier(consdata->resvar, infervar));
      assert(SCIPvarGetUbLocal(consdata->resvar) < 0.5);
      CHECK_OKAY( SCIPaddConflictVar(scip, consdata->resvar) );
      for( i = 0; i < nvars; ++i )
      {
         if( vars[i] != infervar )
         {
            assert(SCIPvarWasFixedEarlier(vars[i], infervar));
            assert(SCIPvarGetLbLocal(vars[i]) > 0.5);
            CHECK_OKAY( SCIPaddConflictVar(scip, vars[i]) );
         }
      }
      *result = SCIP_SUCCESS;
      break;

   default:
      errorMessage("invalid inference information %d in and constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}





/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
DECL_CONSFREE(consFreeAnd)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   CHECK_OKAY( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitAnd NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitAnd NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreAnd NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreAnd NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolAnd NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
DECL_CONSEXITSOL(consExitsolAnd)
{
   CONSDATA* consdata;
   int c;

   /* release and free the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      CHECK_OKAY( consdataFreeRows(scip, consdata) );
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteAnd)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   CHECK_OKAY( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransAnd)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create target constraint data */
   CHECK_OKAY( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr,
         sourcedata->nvars, sourcedata->vars, sourcedata->resvar) );

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpAnd)
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
DECL_CONSSEPA(consSepaAnd)
{  /*lint --e{715}*/
   Bool separated;
   int c;

   *result = SCIP_DIDNOTFIND;

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   } 

   /* combine constraints to get more cuts */
   /**@todo combine constraints to get further cuts */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpAnd)
{  /*lint --e{715}*/
   Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      CHECK_OKAY( checkCons(scip, conss[i], NULL, FALSE, &violated) );
      if( violated )
      {
         Bool separated;

         CHECK_OKAY( separateCons(scip, conss[i], &separated) );
         assert(separated); /* because the solution is integral, the separation always finds a cut */
         *result = SCIP_SEPARATED;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsAnd)
{  /*lint --e{715}*/
   Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      CHECK_OKAY( checkCons(scip, conss[i], NULL, TRUE, &violated) );
      if( violated )
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
DECL_CONSCHECK(consCheckAnd)
{  /*lint --e{715}*/
   Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      CHECK_OKAY( checkCons(scip, conss[i], sol, TRUE, &violated) );
      if( violated )
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
DECL_CONSPROP(consPropAnd)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool cutoff;
   int nfixedvars;
   int c;
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   nfixedvars = 0;

   /* propagate all useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      CHECK_OKAY( propagateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &nfixedvars) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
DECL_CONSPRESOL(consPresolAnd)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONS* cons;
   CONSDATA* consdata;
   Bool cutoff;
   Bool redundant;
   Bool aggregated;
   int oldnfixedvars;
   int oldnaggrvars;
   int c;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;
   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   cutoff = FALSE;
   for( c = 0; c < nconss && !cutoff; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->propagated = FALSE;

      /* propagate constraint */
      CHECK_OKAY( propagateCons(scip, cons, conshdlrdata->eventhdlr, &cutoff, nfixedvars) );

      /* remove all variables that are fixed to one */
      CHECK_OKAY( applyFixings(scip, cons, conshdlrdata->eventhdlr) );

      if( !cutoff && !SCIPconsIsDeleted(cons) && !SCIPconsIsModifiable(cons) )
      {
         assert(consdata->nvars >= 1); /* otherwise, propagateCons() has deleted the constraint */

         /* if only one variable is left, the resultant has to be equal to this single variable */
         if( consdata->nvars == 1 )
         {
            debugMessage("and constraint <%s> has only one variable not fixed to 1.0\n", SCIPconsGetName(cons));
            
            assert(consdata->vars != NULL);
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
            
            /* aggregate variables: resultant - operand == 0 */
            CHECK_OKAY( SCIPaggregateVars(scip, consdata->resvar, consdata->vars[0], 1.0, -1.0, 0.0,
                  &cutoff, &redundant, &aggregated) );
            assert(redundant);
            if( aggregated )
               (*naggrvars)++;

            /* delete constraint */
            CHECK_OKAY( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
         }
      }
   }

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** conflict variable resolving method of constraint handler */
static
DECL_CONSRESCVAR(consRescvarAnd)
{  /*lint --e{715}*/
   CHECK_OKAY( resolveConflict(scip, cons, infervar, (PROPRULE)inferinfo, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockAnd)
{  /*lint --e{715}*/
   consdataLockAllRoundings(SCIPconsGetData(cons), nlockspos, nlocksneg);

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockAnd)
{  /*lint --e{715}*/
   consdataUnlockAllRoundings(SCIPconsGetData(cons), nunlockspos, nunlocksneg);

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveAnd NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveAnd NULL


/** constraint enabling notification method of constraint handler */
#define consEnableAnd NULL


/** constraint disabling notification method of constraint handler */
#define consDisableAnd NULL


/** constraint display method of constraint handler */
static
DECL_CONSPRINT(consPrintAnd)
{
   consdataPrint(scip, SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecAnd)
{  /*lint --e{715}*/
   CONSDATA* consdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(event != NULL);

   consdata = (CONSDATA*)eventdata;
   assert(consdata != NULL);

   /* check, if the variable was fixed to zero */
   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_UBTIGHTENED )
      consdata->nofixedzero = FALSE;

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for and constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrAnd(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound tighten events */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL,
         NULL, eventExecAnd,
         NULL) );

   /* create constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, CONSHDLR_NEEDSCONS,
         consFreeAnd, consInitAnd, consExitAnd, 
         consInitpreAnd, consExitpreAnd, consInitsolAnd, consExitsolAnd,
         consDeleteAnd, consTransAnd, consInitlpAnd,
         consSepaAnd, consEnfolpAnd, consEnfopsAnd, consCheckAnd, 
         consPropAnd, consPresolAnd, consRescvarAnd,
         consLockAnd, consUnlockAnd,
         consActiveAnd, consDeactiveAnd, 
         consEnableAnd, consDisableAnd,
         consPrintAnd,
         conshdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a and constraint */
RETCODE SCIPcreateConsAnd(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   VAR*             resvar,             /**< resultant variable of the operation */
   int              nvars,              /**< number of operator variables in the constraint */
   VAR**            vars,               /**< array with operator variables of constraint */
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
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;

   /* find the and constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("and constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, resvar) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, removeable) );

   return SCIP_OKAY;
}
