/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (c) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_and.c,v 1.62 2005/09/22 14:43:47 bzfpfend Exp $"

/**@file   cons_and.c
 * @brief  constraint handler for and constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_and.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "and"
#define CONSHDLR_DESC          "constraint handler for and constraints: r = and(x1, ..., xn)"
#define CONSHDLR_SEPAPRIORITY   +850000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -850000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -850000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define EVENTHDLR_NAME         "and"
#define EVENTHDLR_DESC         "bound change event handler for and constraints"




/*
 * Data structures
 */

/** constraint data for and constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables in the and operation */
   SCIP_VAR*             resvar;             /**< resultant variable */
   SCIP_ROW**            rows;               /**< rows for linear relaxation of and constraint */
   int                   nvars;              /**< number of variables in and operation */
   int                   varssize;           /**< size of vars array */
   int                   rowssize;           /**< size of rows array */
   int                   watchedvar1;        /**< position of first watched operator variable */
   int                   watchedvar2;        /**< position of second watched operator variable */
   int                   filterpos1;         /**< event filter position of first watched operator variable */
   int                   filterpos2;         /**< event filter position of second watched operator variable */
   unsigned int          propagated:1;       /**< is constraint already preprocessed/propagated? */
   unsigned int          nofixedzero:1;      /**< is none of the opereator variables fixed to FALSE? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events on watched variables */
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

#if 0
/** installs rounding locks for the given variable in the given and constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< and constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}
#endif

/** removes rounding locks for the given variable in the given and constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< and constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** creates constaint handler data */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   /* get event handler for catching bound change events on variables */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for and constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   
   return SCIP_OKAY;
}

/** frees constraint handler data */
static
SCIP_RETCODE conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
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
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   return consdata->nvars + 1;
}

/** catches events for the watched variable at given position */
static
SCIP_RETCODE consdataCatchWatchedEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< and constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos,                /**< array position of variable to catch bound change events for */
   int*                  filterpos           /**< pointer to store position of event filter entry */
   )
{
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(filterpos != NULL);

   /* catch tightening events for lower bound and relaxed events for upper bounds on watched variable */
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
         eventhdlr, (SCIP_EVENTDATA*)consdata, filterpos) );

   return SCIP_OKAY;
}


/** drops events for the watched variable at given position */
static
SCIP_RETCODE consdataDropWatchedEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< and constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos,                /**< array position of watched variable to drop bound change events for */
   int                   filterpos           /**< position of event filter entry */
   )
{
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(filterpos >= 0);

   /* drop tightening events for lower bound and relaxed events for upper bounds on watched variable */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
         eventhdlr, (SCIP_EVENTDATA*)consdata, filterpos) );

   return SCIP_OKAY;
}

/** catches needed events on all variables of constraint, except the special ones for watched variables */
static
SCIP_RETCODE consdataCatchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< and constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* catch bound change events for both bounds on resultant variable */
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->resvar, SCIP_EVENTTYPE_BOUNDCHANGED, 
         eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

   /* catch tightening events for upper bound and relaxed events for lower bounds on operator variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[i], SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );
   }

   return SCIP_OKAY;
}

/** drops events on all variables of constraint, except the special ones for watched variables */
static
SCIP_RETCODE consdataDropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< and constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* drop bound change events for both bounds on resultant variable */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->resvar, SCIP_EVENTTYPE_BOUNDCHANGED,
         eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

   /* drop tightening events for upper bound and relaxed events for lower bounds on operator variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[i], SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );
   }

   return SCIP_OKAY;
}

/** stores the given variable numbers as watched variables, and updates the event processing */
static
SCIP_RETCODE consdataSwitchWatchedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< and constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   watchedvar1,        /**< new first watched variable */
   int                   watchedvar2         /**< new second watched variable */
   )
{
   assert(consdata != NULL);
   assert(watchedvar1 == -1 || watchedvar1 != watchedvar2);
   assert(watchedvar1 != -1 || watchedvar2 == -1);
   assert(watchedvar1 == -1 || (0 <= watchedvar1 && watchedvar1 < consdata->nvars));
   assert(watchedvar2 == -1 || (0 <= watchedvar2 && watchedvar2 < consdata->nvars));

   /* if one watched variable is equal to the old other watched variable, just switch positions */
   if( watchedvar1 == consdata->watchedvar2 || watchedvar2 == consdata->watchedvar1 )
   {
      int tmp;
      
      tmp = consdata->watchedvar1;
      consdata->watchedvar1 = consdata->watchedvar2;
      consdata->watchedvar2 = tmp;
      tmp = consdata->filterpos1;
      consdata->filterpos1 = consdata->filterpos2;
      consdata->filterpos2 = tmp;
   }
   assert(watchedvar1 == -1 || watchedvar1 != consdata->watchedvar2);
   assert(watchedvar2 == -1 || watchedvar2 != consdata->watchedvar1);

   /* drop events on old watched variables */
   if( consdata->watchedvar1 != -1 && consdata->watchedvar1 != watchedvar1 )
   {
      assert(consdata->filterpos1 != -1);
      SCIP_CALL( consdataDropWatchedEvents(scip, consdata, eventhdlr, consdata->watchedvar1, consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 && consdata->watchedvar2 != watchedvar2 )
   {
      assert(consdata->filterpos2 != -1);
      SCIP_CALL( consdataDropWatchedEvents(scip, consdata, eventhdlr, consdata->watchedvar2, consdata->filterpos2) );
   }

   /* catch events on new watched variables */
   if( watchedvar1 != -1 && watchedvar1 != consdata->watchedvar1 )
   {
      SCIP_CALL( consdataCatchWatchedEvents(scip, consdata, eventhdlr, watchedvar1, &consdata->filterpos1) );
   }
   if( watchedvar2 != -1 && watchedvar2 != consdata->watchedvar2 )
   {
      SCIP_CALL( consdataCatchWatchedEvents(scip, consdata, eventhdlr, watchedvar2, &consdata->filterpos2) );
   }

   /* set the new watched variables */
   consdata->watchedvar1 = watchedvar1;
   consdata->watchedvar2 = watchedvar2;

   return SCIP_OKAY;
}

/** creates constraint data for and constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   nvars,              /**< number of variables in the and operation */
   SCIP_VAR**            vars,               /**< variables in and operation */
   SCIP_VAR*             resvar              /**< resultant variable */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(resvar != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
   (*consdata)->resvar = resvar;
   (*consdata)->rows = NULL;
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->rowssize = 0;
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->filterpos1 = -1;
   (*consdata)->filterpos2 = -1;
   (*consdata)->propagated = FALSE;
   (*consdata)->nofixedzero = FALSE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->resvar, &(*consdata)->resvar) );

      /* catch needed events on variables */
      SCIP_CALL( consdataCatchEvents(scip, *consdata, eventhdlr) );
   }

   return SCIP_OKAY;
}

/** releases LP rows of constraint data and frees rows array */
static
SCIP_RETCODE consdataFreeRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int r;

   assert(consdata != NULL);

   if( consdata->rows != NULL )
   {
      for( r = 0; r < consdataGetNRows(consdata); ++r )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rows[r]) );
      }
      SCIPfreeBlockMemoryArray(scip, &consdata->rows, consdata->rowssize);
   }

   return SCIP_OKAY;
}

/** frees constraint data for and constraint */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to the constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   if( SCIPisTransformed(scip) )
   {
      /* drop events for watched variables */
      SCIP_CALL( consdataSwitchWatchedvars(scip, *consdata, eventhdlr, -1, -1) );

      /* drop all other events on variables */
      SCIP_CALL( consdataDropEvents(scip, *consdata, eventhdlr) );
   }
   else
   {
      assert((*consdata)->watchedvar1 == -1);
      assert((*consdata)->watchedvar2 == -1);
   }

   /* release and free the rows */
   SCIP_CALL( consdataFreeRows(scip, *consdata) );

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);
 
   return SCIP_OKAY;
}

/** prints and constraint to file stream */
static
void consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< and constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);

   /* print coefficients */
   SCIPinfoMessage(scip, file, "<%s> == and(", SCIPvarGetName(consdata->resvar));
   for( v = 0; v < consdata->nvars; ++v )
   {
      if( v > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "<%s>", SCIPvarGetName(consdata->vars[v]));
   }
   SCIPinfoMessage(scip, file, ")\n");
}

/** deletes coefficient at given position from and constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< and constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(consdata->vars[pos]));

   /* remove the rounding locks of the variable */
   SCIP_CALL( unlockRounding(scip, cons, consdata->vars[pos]) );

   if( SCIPconsIsTransformed(cons) )
   {
      /* if the position is watched, stop watching the position */
      if( consdata->watchedvar1 == pos )
      {
         SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, consdata->watchedvar2, -1) );
      }
      if( consdata->watchedvar2 == pos )
      {
         SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, consdata->watchedvar1, -1) );
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
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< and constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
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
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
      }
      else
         ++v;
   }

   SCIPdebugMessage("after fixings: ");
   SCIPdebug(consdataPrint(scip, consdata, NULL));

   return SCIP_OKAY;
}

/** creates LP rows corresponding to and constraint:
 *   - for each operator variable vi:  resvar - vi            <= 0
 *   - one additional row:             resvar - v1 - ... - vn >= n-1
 */
static 
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSDATA* consdata;
   char rowname[SCIP_MAXSTRLEN];
   int nvars;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows == NULL);

   nvars = consdata->nvars;

   /* get memory for rows */
   consdata->rowssize = consdataGetNRows(consdata);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->rows, consdata->rowssize) );
   assert(consdata->rowssize == nvars+1);

   /* create operator rows */
   for( i = 0; i < nvars; ++i )
   {
      sprintf(rowname, "%s_%d", SCIPconsGetName(cons), i);
      SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[i], rowname, -SCIPinfinity(scip), 0.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[i], consdata->resvar, 1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[i], consdata->vars[i], -1.0) );
   }

   /* create additional row */
   sprintf(rowname, "%s_add", SCIPconsGetName(cons));
   SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[nvars], rowname, -consdata->nvars + 1.0, SCIPinfinity(scip),
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[nvars], consdata->resvar, 1.0) );
   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[nvars], nvars, consdata->vars, -1.0) );

   return SCIP_OKAY;
}  

/** adds linear relaxation of and constraint to the LP */
static 
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSDATA* consdata;
   int r;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }

   for( r = 0; r < consdataGetNRows(consdata); ++r )
   {
      SCIP_CALL( SCIPaddCut(scip, consdata->rows[r], FALSE) );
   }

   return SCIP_OKAY;
}

/** checks and constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool             checklprows,        /**< should LP rows be checked? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool mustcheck;
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
      SCIP_Real solval;
      int i;

      /* increase age of constraint; age is reset to zero, if a violation was found */
      SCIP_CALL( SCIPincConsAge(scip, cons) );

      /* check, if all operator variables are TRUE */
      for( i = 0; i < consdata->nvars; ++i )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[i]);
         assert(SCIPisFeasIntegral(scip, solval));
         if( solval < 0.5 )
            break;
      }

      /* if all operator variables are TRUE, the resultant has to be TRUE, otherwise, the resultant has to be FALSE */
      solval = SCIPgetSolVal(scip, sol, consdata->resvar);
      assert(SCIPisFeasIntegral(scip, solval));

      if( (i == consdata->nvars) != (solval > 0.5) )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *violated = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** separates current LP solution */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_Bool*            separated           /**< pointer to store whether a cut was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real feasibility;
   int r;

   assert(separated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *separated = FALSE;

   /* create all necessary rows for the linear relaxation */
   if( consdata->rows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->rows != NULL);

   /* test all rows for feasibility and add infeasible rows */
   for( r = 0; r < consdataGetNRows(consdata); ++r )
   {
      if( !SCIProwIsInLP(consdata->rows[r]) )
      {
         feasibility = SCIPgetRowLPFeasibility(scip, consdata->rows[r]);
         if( SCIPisFeasNegative(scip, feasibility) )
         {
            SCIP_CALL( SCIPaddCut(scip, consdata->rows[r], FALSE) );
            *separated = TRUE;
         }
      }            
   }

   return SCIP_OKAY;
}

/** analyzes conflicting TRUE assignment to resultant of given constraint, and adds conflict clause to problem */
static
SCIP_RETCODE analyzeConflictOne(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< and constraint that detected the conflict */
   int                   falsepos            /**< position of operand that is fixed to FALSE */
   )
{
   SCIP_CONSDATA* consdata;

   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarGetLbLocal(consdata->resvar) > 0.5);
   assert(0 <= falsepos && falsepos < consdata->nvars);
   assert(SCIPvarGetUbLocal(consdata->vars[falsepos]) < 0.5);

   /* initialize conflict analysis, and add resultant and single operand variable to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );
   SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
   SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[falsepos]) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** analyzes conflicting FALSE assignment to resultant of given constraint, and adds conflict clause to problem */
static
SCIP_RETCODE analyzeConflictZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< or constraint that detected the conflict */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   assert(!SCIPconsIsModifiable(cons));

   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPvarGetUbLocal(consdata->resvar) < 0.5);

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );
   SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(SCIPvarGetLbLocal(consdata->vars[v]) > 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** propagates constraint with the following rules:
 *   (1) v_i = FALSE                                  =>  r   = FALSE
 *   (2) r   = TRUE                                   =>  v_i = TRUE for all i
 *   (3) v_i = TRUE for all i                         =>  r   = TRUE
 *   (4) r   = FALSE, v_i = TRUE for all i except j   =>  v_j = FALSE
 */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< and constraint to be processed */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* resvar;
   SCIP_VAR** vars;
   int nvars;
   int watchedvar1;
   int watchedvar2;
   int i;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

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
   SCIP_CALL( SCIPincConsAge(scip, cons) );

   /* if one of the operator variables was fixed to FALSE, the resultant can be fixed to FALSE (rule (1)) */
   if( !consdata->nofixedzero )
   {
      for( i = 0; i < nvars && SCIPvarGetUbLocal(vars[i]) > 0.5; ++i ) /* search fixed operator */
      {}
      if( i < nvars )
      {
         SCIPdebugMessage("constraint <%s>: operator var <%s> fixed to 0.0 -> fix resultant <%s> to 0.0\n",
            SCIPconsGetName(cons), SCIPvarGetName(vars[i]), SCIPvarGetName(resvar));
         SCIP_CALL( SCIPinferBinvarCons(scip, resvar, FALSE, cons, (int)PROPRULE_1, &infeasible, &tightened) );
         if( infeasible )
         {
            /* use conflict analysis to get a conflict clause out of the conflicting assignment */
            SCIP_CALL( analyzeConflictOne(scip, cons, i) );
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else
         {
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
            if( tightened )
            {
               SCIP_CALL( SCIPresetConsAge(scip, cons) );
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
         SCIPdebugMessage("constraint <%s>: resultant var <%s> fixed to 1.0 -> fix operator var <%s> to 1.0\n",
            SCIPconsGetName(cons), SCIPvarGetName(resvar), SCIPvarGetName(vars[i]));
         SCIP_CALL( SCIPinferBinvarCons(scip, vars[i], TRUE, cons, (int)PROPRULE_2, &infeasible, &tightened) );
         if( infeasible )
         {
            /* use conflict analysis to get a conflict clause out of the conflicting assignment */
            SCIP_CALL( analyzeConflictOne(scip, cons, i) );
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( tightened )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            (*nfixedvars)++;
         }
      }

      if( !(*cutoff) )
      {
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
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
      
      SCIPdebugMessage("constraint <%s>: all operator vars fixed to 1.0 -> fix resultant <%s> to 1.0\n",
         SCIPconsGetName(cons), SCIPvarGetName(resvar));
      SCIP_CALL( SCIPinferBinvarCons(scip, resvar, TRUE, cons, (int)PROPRULE_3, &infeasible, &tightened) );
      if( infeasible )
      {
         /* use conflict analysis to get a conflict clause out of the conflicting assignment */
         SCIP_CALL( analyzeConflictZero(scip, cons) );
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         if( tightened )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
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
      
      SCIPdebugMessage("constraint <%s>: resultant <%s> fixed to 0.0, only one unfixed operand -> fix operand <%s> to 0.0\n",
         SCIPconsGetName(cons), SCIPvarGetName(resvar), SCIPvarGetName(vars[watchedvar1]));
      SCIP_CALL( SCIPinferBinvarCons(scip, vars[watchedvar1], FALSE, cons, (int)PROPRULE_4, &infeasible, &tightened) );
      if( infeasible )
      {
         /* use conflict analysis to get a conflict clause out of the conflicting assignment */
         SCIP_CALL( analyzeConflictZero(scip, cons) );
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         if( tightened )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            (*nfixedvars)++;
         }
      }

      return SCIP_OKAY;
   }

   /* switch to the new watched variables */
   SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, watchedvar1, watchedvar2) );

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
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   PROPRULE         proprule,           /**< propagation rule that deduced the value */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
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
      assert(SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < 0.5);
      assert(infervar == consdata->resvar);
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPvarGetUbAtIndex(vars[i], bdchgidx, FALSE) < 0.5 )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
            break;
         }
      }
      assert(i < nvars);
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_2:
      /* the operand variable was infered to TRUE, because the resultant was TRUE */
      assert(SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) > 0.5);
      assert(SCIPvarGetLbAtIndex(consdata->resvar, bdchgidx, FALSE) > 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_3:
      /* the resultant was infered to TRUE, because all operand variables were TRUE */
      assert(SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) > 0.5);
      assert(infervar == consdata->resvar);
      for( i = 0; i < nvars; ++i )
      {
         assert(SCIPvarGetLbAtIndex(vars[i], bdchgidx, FALSE) > 0.5);
         SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_4:
      /* the operand variable was infered to FALSE, because the resultant was FALSE and all other operands were TRUE */
      assert(SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < 0.5);
      assert(SCIPvarGetUbAtIndex(consdata->resvar, bdchgidx, FALSE) < 0.5);
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->resvar) );
      for( i = 0; i < nvars; ++i )
      {
         if( vars[i] != infervar )
         {
            assert(SCIPvarGetLbAtIndex(vars[i], bdchgidx, FALSE) > 0.5);
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
         }
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in and constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}





/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

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
SCIP_DECL_CONSEXITSOL(consExitsolAnd)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   /* release and free the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      SCIP_CALL( consdataFreeRows(scip, consdata) );
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create target constraint data */
   SCIP_CALL( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr,
         sourcedata->nvars, sourcedata->vars, sourcedata->resvar) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpAnd)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( SCIPconsIsInitial(conss[i]) )
      {
         SCIP_CALL( addRelaxation(scip, conss[i]) );
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpAnd)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int c;

   *result = SCIP_DIDNOTFIND;

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   } 

   /* combine constraints to get more cuts */
   /**@todo combine constraints to get further cuts */

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
#define consSepasolAnd NULL   /*??????????????????*/


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpAnd)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, FALSE, &violated) );
      if( violated )
      {
         SCIP_Bool separated;

         SCIP_CALL( separateCons(scip, conss[i], &separated) );
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
SCIP_DECL_CONSENFOPS(consEnfopsAnd)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, TRUE, &violated) );
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
SCIP_DECL_CONSCHECK(consCheckAnd)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, checklprows, &violated) );
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
SCIP_DECL_CONSPROP(consPropAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   int nfixedvars;
   int c;
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   nfixedvars = 0;

   /* propagate all useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( propagateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &nfixedvars) );
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
SCIP_DECL_CONSPRESOL(consPresolAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;
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
      SCIP_CALL( propagateCons(scip, cons, conshdlrdata->eventhdlr, &cutoff, nfixedvars) );

      /* remove all variables that are fixed to one */
      SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr) );

      if( !cutoff && !SCIPconsIsDeleted(cons) && !SCIPconsIsModifiable(cons) )
      {
         assert(consdata->nvars >= 1); /* otherwise, propagateCons() has deleted the constraint */

         /* if only one variable is left, the resultant has to be equal to this single variable */
         if( consdata->nvars == 1 )
         {
            SCIPdebugMessage("and constraint <%s> has only one variable not fixed to 1.0\n", SCIPconsGetName(cons));
            
            assert(consdata->vars != NULL);
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
            
            /* aggregate variables: resultant - operand == 0 */
            SCIP_CALL( SCIPaggregateVars(scip, consdata->resvar, consdata->vars[0], 1.0, -1.0, 0.0,
                  &cutoff, &redundant, &aggregated) );
            assert(redundant);
            if( aggregated )
               (*naggrvars)++;

            /* delete constraint */
            SCIP_CALL( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
         }
      }
   }

   /**@todo preprocess pairs of and constraints */

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropAnd)
{  /*lint --e{715}*/
   SCIP_CALL( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockAnd)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* resultant variable */
   SCIP_CALL( SCIPaddVarLocks(scip, consdata->resvar, nlockspos + nlocksneg, nlockspos + nlocksneg) );

   /* operand variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

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
SCIP_DECL_CONSPRINT(consPrintAnd)
{  /*lint --e{715}*/
   consdataPrint(scip, SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecAnd)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(event != NULL);

   consdata = (SCIP_CONSDATA*)eventdata;
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
SCIP_RETCODE SCIPincludeConshdlrAnd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create event handler for events on variables */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecAnd,
         NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeAnd, consInitAnd, consExitAnd, 
         consInitpreAnd, consExitpreAnd, consInitsolAnd, consExitsolAnd,
         consDeleteAnd, consTransAnd, consInitlpAnd,
         consSepalpAnd, consSepasolAnd, consEnfolpAnd, consEnfopsAnd, consCheckAnd, 
         consPropAnd, consPresolAnd, consRespropAnd, consLockAnd,
         consActiveAnd, consDeactiveAnd, 
         consEnableAnd, consDisableAnd,
         consPrintAnd,
         conshdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a and constraint */
SCIP_RETCODE SCIPcreateConsAnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             resvar,             /**< resultant variable of the operation */
   int                   nvars,              /**< number of operator variables in the constraint */
   SCIP_VAR**            vars,               /**< array with operator variables of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removeable          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   /* find the and constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("and constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, resvar) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removeable) );

   return SCIP_OKAY;
}

