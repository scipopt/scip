/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (c) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_xor.c,v 1.45 2006/04/10 16:15:24 bzfpfend Exp $"

/**@file   cons_xor.c
 * @brief  constraint handler for xor constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_xor.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "xor"
#define CONSHDLR_DESC          "constraint handler for xor constraints: r = xor(x1, ..., xn)"
#define CONSHDLR_SEPAPRIORITY   +850200 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -850200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -850200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define EVENTHDLR_NAME         "xor"
#define EVENTHDLR_DESC         "event handler for xor constraints"

#define NROWS 4



/*
 * Data structures
 */

/** constraint data for xor constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables (including resultant) in the xor operation */
   SCIP_VAR*             intvar;             /**< internal variable for LP relaxation */
   SCIP_ROW*             rows[NROWS];        /**< rows for linear relaxation of xor constraint */
   int                   nvars;              /**< number of variables (including resultant) in xor operation */
   int                   varssize;           /**< size of vars array */
   int                   watchedvar1;        /**< position of first watched operator variable */
   int                   watchedvar2;        /**< position of second watched operator variable */
   int                   filterpos1;         /**< event filter position of first watched operator variable */
   int                   filterpos2;         /**< event filter position of second watched operator variable */
   unsigned int          propagated:1;       /**< is constraint already preprocessed/propagated? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for events on watched variables */
};




/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_1,                          /**< all except one variable fixed  =>  fix remaining variable */
   PROPRULE_INVALID                     /**< propagation was applied without a specific propagation rule */
};
typedef enum Proprule PROPRULE;




/*
 * Local methods
 */

/** installs rounding locks for the given variable in the given xor constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given xor constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
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

   /* get event handler for catching events on watched variables */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for xor constraints not found\n");
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

/** stores the given variable numbers as watched variables, and updates the event processing */
static
SCIP_RETCODE consdataSwitchWatchedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< xor constraint data */
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
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar1], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 && consdata->watchedvar2 != watchedvar2 )
   {
      assert(consdata->filterpos2 != -1);
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar2], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, consdata->filterpos2) );
   }

   /* catch events on new watched variables */
   if( watchedvar1 != -1 && watchedvar1 != consdata->watchedvar1 )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[watchedvar1], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, &consdata->filterpos1) );
   }
   if( watchedvar2 != -1 && watchedvar2 != consdata->watchedvar2 )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[watchedvar2], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, &consdata->filterpos2) );
   }

   /* set the new watched variables */
   consdata->watchedvar1 = watchedvar1;
   consdata->watchedvar2 = watchedvar2;

   return SCIP_OKAY;
}

/** creates constraint data for xor constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
   int                   nvars,              /**< number of variables in the xor operation */
   SCIP_VAR**            vars,               /**< variables in xor operation */
   SCIP_VAR*             resvar              /**< resultant variable */
   )
{
   int r;

   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(resvar != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   /* store the resultant variable as first variable in the array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vars, nvars+1) );
   (*consdata)->vars[0] = resvar;
   BMScopyMemoryArray(&(*consdata)->vars[1], vars, nvars);

   (*consdata)->intvar = NULL;
   for( r = 0; r < NROWS; ++r )
      (*consdata)->rows[r] = NULL;
   (*consdata)->nvars = nvars+1;
   (*consdata)->varssize = nvars+1;
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->filterpos1 = -1;
   (*consdata)->filterpos2 = -1;
   (*consdata)->propagated = FALSE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
   }

   return SCIP_OKAY;
}

/** releases LP row of constraint data */
static
SCIP_RETCODE consdataFreeRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int r;

   assert(consdata != NULL);

   for( r = 0; r < NROWS; ++r )
   {
      if( consdata->rows[r] != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rows[r]) );
      }
   }

   return SCIP_OKAY;
}

/** frees constraint data for xor constraint */
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
   }
   else
   {
      assert((*consdata)->watchedvar1 == -1);
      assert((*consdata)->watchedvar2 == -1);
   }

   /* release LP row */
   SCIP_CALL( consdataFreeRows(scip, *consdata) );

   /* release internal variable */
   if( (*consdata)->intvar != NULL )
   {
      /* if internal variable is not defined anymore by any checked constraint, delete it from the problem;
       * otherwise, it could be fixed to a wrong value in dual presolving
       */
      if( SCIPvarMayRoundDown((*consdata)->intvar) && SCIPvarMayRoundUp((*consdata)->intvar) )
      {
         SCIP_CALL( SCIPdelVar(scip, (*consdata)->intvar) );
      }
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->intvar) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);
 
   return SCIP_OKAY;
}

/** prints xor constraint to file stream */
static
void consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< xor constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);
   assert(consdata->nvars >= 1);

   /* print coefficients */
   SCIPinfoMessage(scip, file, "<%s> == xor(", SCIPvarGetName(consdata->vars[0]));
   for( v = 1; v < consdata->nvars; ++v )
   {
      if( v > 1 )
         SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "<%s>", SCIPvarGetName(consdata->vars[v]));
   }
   SCIPinfoMessage(scip, file, ")\n");
}

/** deletes coefficient at given position from xor constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
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

   /* remove the rounding locks of the deleted variable */
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

/** deletes all zero-fixed variables and all pairs of one-fixed variables */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int v;
   int lastfixedone;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);

   v = 0;
   lastfixedone = -1;
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

      if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
      }
      else if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         if( lastfixedone >= 0 )
         {
            /* delete variable at position v > lastfixedone first, because this doesn't affect previous array positions */
            assert(v > lastfixedone);
            SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
            SCIP_CALL( delCoefPos(scip, cons, eventhdlr, lastfixedone) );
            lastfixedone = -1;
         }
         else
         {
            lastfixedone = v;
            ++v;
         }
      }
      else
         ++v;
   }

   SCIPdebugMessage("after fixings: ");
   SCIPdebug(consdataPrint(scip, consdata, NULL));

   return SCIP_OKAY;
}

/** creates LP row corresponding to xor constraint: 
 *    resvar + v1 + ... + vn - 2q == 0
 *  with internal integer variable q;
 *  in the special case of 2 operand variables r == x xor y, the following linear system is created:
 *    + x - y - r <= 0
 *    - x + y - r <= 0
 *    - x - y + r <= 0
 *    + x + y + r <= 2
 */
static 
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSDATA* consdata;
   char varname[SCIP_MAXSTRLEN];

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows[0] == NULL);

   if( SCIPconsIsModifiable(cons) || consdata->nvars != 3 )  /* resultant is member of vars array */
   {
      /* create internal variable, if not yet existing */
      if( consdata->intvar == NULL )
      {
         int ub;

         sprintf(varname, "%s_int", SCIPconsGetName(cons));
         ub = consdata->nvars/2;
         SCIP_CALL( SCIPcreateVar(scip, &consdata->intvar, varname, 0.0, (SCIP_Real)ub, 0.0,
               consdata->nvars >= 4 ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_BINARY,
               SCIPconsIsInitial(cons), SCIPconsIsRemoveable(cons), NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, consdata->intvar) );

         /* install the rounding locks for the internal variable */
         SCIP_CALL( lockRounding(scip, cons, consdata->intvar) );
      }

      /* create LP row (resultant variable is also stored in vars array) */
      SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[0], SCIPconsGetName(cons), 0.0, 0.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[0], consdata->intvar, -2.0) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[0], consdata->nvars, consdata->vars, 1.0) );
   }
   else
   {
      char rowname[SCIP_MAXSTRLEN];
      int r;

      /* create the <= 0 rows with one positive sign */
      for( r = 0; r < 3; ++r )
      {
         int v;

         sprintf(rowname, "%s_%d", SCIPconsGetName(cons), r);
         SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[r], rowname, -SCIPinfinity(scip), 0.0,
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
         for( v = 0; v < 3; ++v )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[r], consdata->vars[v], v == r ? +1.0 : -1.0) );
         }
      }

      /* create the <= 2 row with all positive signs */
      sprintf(rowname, "%s_3", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[3], rowname, -SCIPinfinity(scip), 2.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[3], consdata->nvars, consdata->vars, 1.0) );
   }

   return SCIP_OKAY;
}  

/** adds linear relaxation of or constraint to the LP */
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

   if( consdata->rows[0] == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->rows[0] != NULL);
   for( r = 0; r < NROWS; ++r )
   {
      if( consdata->rows[r] != NULL )
      {
         SCIP_CALL( SCIPaddCut(scip, NULL, consdata->rows[r], FALSE) );
      }
   }

   return SCIP_OKAY;
}

/** returns whether all rows of the LP relaxation are in the current LP */
static
SCIP_Bool allRowsInLP(
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->rows[0] == NULL )      /* LP relaxation does not exist */
      return FALSE;
   else
   {
      int r;
      for( r = 0; r < NROWS; ++r )
      {
         if( consdata->rows[r] != NULL && !SCIProwIsInLP(consdata->rows[r]) )
            return FALSE;
      }
      return TRUE;
   }
}

/** checks xor constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
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

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   /* check feasibility of constraint if necessary */
   if( checklprows || !allRowsInLP(consdata) )
   {
      SCIP_Real solval;
      int i;
      SCIP_Bool odd;

      /* increase age of constraint; age is reset to zero, if a violation was found */
      SCIP_CALL( SCIPincConsAge(scip, cons) );
      
      /* check, if all variables (including resultant) sum up to an even value */
      odd = FALSE;
      for( i = 0; i < consdata->nvars; ++i )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[i]);
         assert(SCIPisFeasIntegral(scip, solval));
         odd = (odd != (solval > 0.5));
      }
      if( odd )
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
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
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

   /* create row for the linear relaxation */
   if( consdata->rows[0] == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->rows[0] != NULL);

   /* test rows for feasibility and add it, if it is infeasible */
   for( r = 0; r < NROWS; ++r )
   {
      if( sol != NULL || (consdata->rows[r] != NULL && !SCIProwIsInLP(consdata->rows[r])) )
      {
         feasibility = SCIPgetRowSolFeasibility(scip, consdata->rows[r], sol);
         if( SCIPisFeasNegative(scip, feasibility) )
         {
            SCIP_CALL( SCIPaddCut(scip, sol, consdata->rows[r], FALSE) );
            *separated = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< xor constraint that detected the conflict */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** propagates constraint with the following rule:
 *   (1) all except one variable fixed  =>  fix remaining variable
 */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint to be processed */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int watchedvar1;
   int watchedvar2;
   int i;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Bool odd;

   assert(cutoff != NULL);
   assert(nfixedvars != NULL);

   /* propagation can only be applied, if we know all operator variables */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   nvars = consdata->nvars;

   /* don't process the constraint, if the watched variables weren't fixed to any value since last propagation call */
   if( consdata->propagated )
      return SCIP_OKAY;

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   if( !SCIPinRepropagation(scip) )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   /* propagation can not be applied, if we have at least two unfixed variables left;
    * that means, we only have to watch (i.e. capture events) of two variables, and switch to other variables
    * if these ones get fixed
    */
   watchedvar1 = consdata->watchedvar1;
   watchedvar2 = consdata->watchedvar2;

   /* check, if watched variables are still unfixed */
   if( watchedvar1 != -1 )
   {
      if( SCIPvarGetLbLocal(vars[watchedvar1]) > 0.5 || SCIPvarGetUbLocal(vars[watchedvar1]) < 0.5 )
         watchedvar1 = -1;
   }
   if( watchedvar2 != -1 )
   {
      if( SCIPvarGetLbLocal(vars[watchedvar2]) > 0.5 || SCIPvarGetUbLocal(vars[watchedvar2]) < 0.5 )
         watchedvar2 = -1;
   }
   
   /* if only one watched variable is still unfixed, make it the first one */
   if( watchedvar1 == -1 )
   {
      watchedvar1 = watchedvar2;
      watchedvar2 = -1;
   }
   assert(watchedvar1 != -1 || watchedvar2 == -1);

   /* if the watched variables are invalid (fixed), find new ones if existing; count the parity */
   odd = FALSE;
   if( watchedvar2 == -1 )
   {
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPvarGetLbLocal(vars[i]) > 0.5 )
            odd = !odd;
         else if( SCIPvarGetUbLocal(vars[i]) > 0.5 )
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

   /* if all variables are fixed, we can decide the feasibility of the constraint */
   if( watchedvar1 == -1 )
   {
      assert(watchedvar2 == -1);
      
      SCIPdebugMessage("constraint <%s>: all vars fixed -> constraint is %s\n",
         SCIPconsGetName(cons), odd ? "infeasible" : "feasible");
      if( odd )
      {
         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflict(scip, cons) );

         *cutoff = TRUE;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );

      return SCIP_OKAY;
   }

   /* if only one variable is not fixed, this variable can be deduced */
   if( watchedvar2 == -1 )
   {
      assert(watchedvar1 != -1);
      
      SCIPdebugMessage("constraint <%s>: only one unfixed variable -> fix <%s> to %d\n",
         SCIPconsGetName(cons), SCIPvarGetName(vars[watchedvar1]), odd);
      SCIP_CALL( SCIPinferBinvarCons(scip, vars[watchedvar1], odd, cons, (int)PROPRULE_1, &infeasible, &tightened) );
      assert(!infeasible);
      assert(tightened);
      (*nfixedvars)++;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );

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
 *   (1) all except one variable fixed  =>  fix remaining variable
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
      /* the variable was infered, because all other variables were fixed */
      for( i = 0; i < nvars; ++i )
      {
         if( vars[i] != infervar )
         {
            assert(SCIPisEQ(scip, SCIPvarGetLbAtIndex(vars[i], bdchgidx, FALSE), 
                  SCIPvarGetUbAtIndex(vars[i], bdchgidx, FALSE)));
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
         }
         else
            assert(SCIPisEQ(scip, SCIPvarGetLbAtIndex(vars[i], bdchgidx, TRUE), 
                  SCIPvarGetUbAtIndex(vars[i], bdchgidx, TRUE)));
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in xor constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}





/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeXor)
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
#define consInitXor NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitXor NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreXor NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreXor NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolXor NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolXor)
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
SCIP_DECL_CONSDELETE(consDeleteXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->nvars >= 1);
   assert(sourcedata->vars != NULL);

   /* create target constraint data */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->nvars-1, &sourcedata->vars[1], sourcedata->vars[0]) );

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
SCIP_DECL_CONSINITLP(consInitlpXor)
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
SCIP_DECL_CONSSEPALP(consSepalpXor)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int c;

   *result = SCIP_DIDNOTFIND;

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   } 

   /* combine constraints to get more cuts */
   /**@todo combine constraints to get further cuts */

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolXor)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int c;

   *result = SCIP_DIDNOTFIND;

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   } 

   /* combine constraints to get more cuts */
   /**@todo combine constraints to get further cuts */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpXor)
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

         SCIP_CALL( separateCons(scip, conss[i], NULL, &separated) );
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
SCIP_DECL_CONSENFOPS(consEnfopsXor)
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
SCIP_DECL_CONSCHECK(consCheckXor)
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
SCIP_DECL_CONSPROP(consPropXor)
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
SCIP_DECL_CONSPRESOL(consPresolXor)
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

      /* remove all variables that are fixed to zero and all pairs of variables fixed to one */
      SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr) );

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons, conshdlrdata->eventhdlr, &cutoff, nfixedvars) );

      if( !cutoff && !SCIPconsIsDeleted(cons) && !SCIPconsIsModifiable(cons) )
      {
         assert(consdata->nvars >= 2); /* otherwise, propagateCons() has deleted the constraint */

         /* if only two variables are left, both have to be equal */
         if( consdata->nvars == 2 )
         {
            SCIPdebugMessage("xor constraint <%s> has only two unfixed variables, remaining vars are even\n",
               SCIPconsGetName(cons));
            
            assert(consdata->vars != NULL);
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[1]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[1]), 1.0));
            
            /* aggregate variables: vars[0] - vars[1] == 0 */
            SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], 1.0, -1.0, 0.0,
                  &cutoff, &redundant, &aggregated) );
            assert(redundant);
            if( aggregated )
               (*naggrvars)++;

            /* delete constraint */
            SCIP_CALL( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
         }

         /* if only three variables are left, but one of them is fixed to one, the others have to be opposite */
         if( consdata->nvars == 3 )
         {
            assert(consdata->vars != NULL);

            if( SCIPvarGetLbGlobal(consdata->vars[0]) > 0.5 )
            {
               SCIPdebugMessage("xor constraint <%s> has only two unfixed variables, remaining vars are odd\n",
                  SCIPconsGetName(cons));
            
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[1]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[1]), 1.0));
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[2]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[2]), 1.0));
               
               /* aggregate variables: vars[1] + vars[2] == 1 */
               SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[1], consdata->vars[2], 1.0, 1.0, 1.0,
                     &cutoff, &redundant, &aggregated) );
               assert(redundant);
               if( aggregated )
                  (*naggrvars)++;

               /* delete constraint */
               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
            }
            else if( SCIPvarGetLbGlobal(consdata->vars[1]) > 0.5 )
            {
               SCIPdebugMessage("xor constraint <%s> has only two unfixed variables, remaining vars are odd\n",
                  SCIPconsGetName(cons));
            
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[2]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[2]), 1.0));
               
               /* aggregate variables: vars[0] + vars[2] == 1 */
               SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[2], 1.0, 1.0, 1.0,
                     &cutoff, &redundant, &aggregated) );
               assert(redundant);
               if( aggregated )
                  (*naggrvars)++;

               /* delete constraint */
               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
            }
            else if( SCIPvarGetLbGlobal(consdata->vars[2]) > 0.5 )
            {
               SCIPdebugMessage("xor constraint <%s> has only two unfixed variables, remaining vars are odd\n",
                  SCIPconsGetName(cons));
            
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[1]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[1]), 1.0));
               
               /* aggregate variables: vars[0] + vars[1] == 1 */
               SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], 1.0, 1.0, 1.0,
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
   }

   /**@todo preprocess pairs of xor constraints */

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropXor)
{  /*lint --e{715}*/
   SCIP_CALL( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* external variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   /* internal variable */
   if( consdata->intvar != NULL )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->intvar, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveXor NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveXor NULL


/** constraint enabling notification method of constraint handler */
#define consEnableXor NULL


/** constraint disabling notification method of constraint handler */
#define consDisableXor NULL


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintXor)
{  /*lint --e{715}*/
   consdataPrint(scip, SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(event != NULL);

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for xor constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrXor(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create event handler for events on variables */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecXor,
         NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeXor, consInitXor, consExitXor, 
         consInitpreXor, consExitpreXor, consInitsolXor, consExitsolXor,
         consDeleteXor, consTransXor, consInitlpXor,
         consSepalpXor, consSepasolXor, consEnfolpXor, consEnfopsXor, consCheckXor, 
         consPropXor, consPresolXor, consRespropXor, consLockXor,
         consActiveXor, consDeactiveXor, 
         consEnableXor, consDisableXor,
         consPrintXor,
         conshdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a xor constraint */
SCIP_RETCODE SCIPcreateConsXor(
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

   /* find the xor constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("xor constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, nvars, vars, resvar) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removeable) );

   return SCIP_OKAY;
}

/** gets number of variables in xor constraint */
int SCIPgetNVarsXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an xor constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets array of variables in xor constraint */
SCIP_VAR** SCIPgetVarsXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an xor constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

