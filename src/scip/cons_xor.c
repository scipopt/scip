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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_xor.c,v 1.28 2005/03/21 11:37:30 bzfpfend Exp $"

/**@file   cons_xor.c
 * @brief  constraint handler for xor constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_xor.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "xor"
#define CONSHDLR_DESC          "constraint handler for xor constraints: r = xor(x1, ..., xn)"
#define CONSHDLR_SEPAPRIORITY   +850200 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -850200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -850200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             5 /**< frequency for separating cuts; zero means to separate only in the root node */
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




/*
 * Data structures
 */

/** constraint data for xor constraints */
struct ConsData
{
   VAR**            vars;               /**< variables (including resultant) in the xor operation */
   VAR*             intvar;             /**< internal variable for LP relaxation */
   ROW*             row;                /**< row for linear relaxation of xor constraint */
   int              nvars;              /**< number of variables (including resultant) in xor operation */
   int              varssize;           /**< size of vars array */
   int              watchedvar1;        /**< position of first watched operator variable */
   int              watchedvar2;        /**< position of second watched operator variable */
   int              filterpos1;         /**< event filter position of first watched operator variable */
   int              filterpos2;         /**< event filter position of second watched operator variable */
   unsigned int     propagated:1;       /**< is constraint already preprocessed/propagated? */
};

/** constraint handler data */
struct ConshdlrData
{
   EVENTHDLR*       eventhdlr;          /**< event handler for events on watched variables */
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
RETCODE lockRounding(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< xor constraint */
   VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   CHECK_OKAY( SCIPlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given xor constraint */
static
RETCODE unlockRounding(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< xor constraint */
   VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   CHECK_OKAY( SCIPunlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** creates constaint handler data */
static
RETCODE conshdlrdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, conshdlrdata) );

   /* get event handler for catching events on watched variables */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      errorMessage("event handler for xor constraints not found\n");
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

/** stores the given variable numbers as watched variables, and updates the event processing */
static
RETCODE consdataSwitchWatchedvars(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< xor constraint data */
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
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar1], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (EVENTDATA*)consdata, consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 && consdata->watchedvar2 != watchedvar2 )
   {
      assert(consdata->filterpos2 != -1);
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar2], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (EVENTDATA*)consdata, consdata->filterpos2) );
   }

   /* catch events on new watched variables */
   if( watchedvar1 != -1 && watchedvar1 != consdata->watchedvar1 )
   {
      CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[watchedvar1], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (EVENTDATA*)consdata, &consdata->filterpos1) );
   }
   if( watchedvar2 != -1 && watchedvar2 != consdata->watchedvar2 )
   {
      CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[watchedvar2], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (EVENTDATA*)consdata, &consdata->filterpos2) );
   }

   /* set the new watched variables */
   consdata->watchedvar1 = watchedvar1;
   consdata->watchedvar2 = watchedvar2;

   return SCIP_OKAY;
}

/** creates constraint data for xor constraint */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              nvars,              /**< number of variables in the xor operation */
   VAR**            vars,               /**< variables in xor operation */
   VAR*             resvar              /**< resultant variable */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(resvar != NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );

   /* store the resultant variable as first variable in the array */
   CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &(*consdata)->vars, nvars+1) );
   (*consdata)->vars[0] = resvar;
   copyMemoryArray(&(*consdata)->vars[1], vars, nvars);

   (*consdata)->intvar = NULL;
   (*consdata)->row = NULL;
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
      CHECK_OKAY( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
   }

   return SCIP_OKAY;
}

/** releases LP row of constraint data */
static
RETCODE consdataFreeRows(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->row != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &consdata->row) );
   }

   return SCIP_OKAY;
}

/** frees constraint data for xor constraint */
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
      /* drop events for watched variables */
      CHECK_OKAY( consdataSwitchWatchedvars(scip, *consdata, eventhdlr, -1, -1) );
   }
   else
   {
      assert((*consdata)->watchedvar1 == -1);
      assert((*consdata)->watchedvar2 == -1);
   }

   /* release LP row */
   CHECK_OKAY( consdataFreeRows(scip, *consdata) );

   /* release internal variable */
   if( (*consdata)->intvar != NULL )
   {
      CHECK_OKAY( SCIPreleaseVar(scip, &(*consdata)->intvar) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);
 
   return SCIP_OKAY;
}

/** prints xor constraint to file stream */
static
void consdataPrint(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< xor constraint data */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);
   assert(consdata->nvars >= 1);

   if( file == NULL )
      file = stdout;

   /* print coefficients */
   fprintf(file, "<%s> == xor(", SCIPvarGetName(consdata->vars[0]));
   for( v = 1; v < consdata->nvars; ++v )
   {
      if( v > 1 )
         fprintf(file, ", ");
      fprintf(file, "<%s>", SCIPvarGetName(consdata->vars[v]));
   }
   fprintf(file, ")\n");
}

/** deletes coefficient at given position from xor constraint data */
static
RETCODE delCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< xor constraint */
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

   /* remove the rounding locks of the deleted variable */
   CHECK_OKAY( unlockRounding(scip, cons, consdata->vars[pos]) );

   if( SCIPconsIsTransformed(cons) )
   {
      /* if the position is watched, stop watching the position */
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

/** deletes all zero-fixed variables and all pairs of one-fixed variables */
static
RETCODE applyFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< xor constraint */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   CONSDATA* consdata;
   VAR* var;
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
         CHECK_OKAY( delCoefPos(scip, cons, eventhdlr, v) );
      }
      else if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         if( lastfixedone >= 0 )
         {
            /* delete variable at position v > lastfixedone first, because this doesn't affect previous array positions */
            assert(v > lastfixedone);
            CHECK_OKAY( delCoefPos(scip, cons, eventhdlr, v) );
            CHECK_OKAY( delCoefPos(scip, cons, eventhdlr, lastfixedone) );
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

   debugMessage("after fixings: ");
   debug(consdataPrint(scip, consdata, NULL));

   return SCIP_OKAY;
}

/** creates LP row corresponding to xor constraint: 
 *    resvar + v1 + ... + vn - 2q == 0
 *  with internal integer variable q
 */
static 
RETCODE createRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to check */
   )
{
   CONSDATA* consdata;
   char varname[MAXSTRLEN];

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   /* create internal variable, if not yet existing */
   if( consdata->intvar == NULL )
   {
      sprintf(varname, "%s_int", SCIPconsGetName(cons));
      CHECK_OKAY( SCIPcreateVar(scip, &consdata->intvar, varname, 0.0, SCIPfloor(scip, consdata->nvars/2.0), 0.0,
            SCIP_VARTYPE_INTEGER, SCIPconsIsInitial(cons), SCIPconsIsRemoveable(cons), NULL, NULL, NULL, NULL) );
      CHECK_OKAY( SCIPaddVar(scip, consdata->intvar) );

      /* install the rounding locks for the internal variable */
      CHECK_OKAY( lockRounding(scip, cons, consdata->intvar) );
   }

   /* create LP row (resultant variable is also stored in vars array) */
   /**@todo change LP relaxation! */
   CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), 0.0, 0.0,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   CHECK_OKAY( SCIPaddVarToRow(scip, consdata->row, consdata->intvar, -2.0) );
   CHECK_OKAY( SCIPaddVarsToRowSameCoef(scip, consdata->row, consdata->nvars, consdata->vars, 1.0) );

   return SCIP_OKAY;
}  

/** adds linear relaxation of or constraint to the LP */
static 
RETCODE addRelaxation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to check */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      CHECK_OKAY( createRelaxation(scip, cons) );
   }
   CHECK_OKAY( SCIPaddCut(scip, consdata->row, FALSE) );

   return SCIP_OKAY;
}

/** checks xor constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
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

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   /* check feasibility of constraint if necessary */
   if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
   {
      Real solval;
      int i;
      Bool odd;

      /* increase age of constraint; age is reset to zero, if a violation was found */
      CHECK_OKAY( SCIPincConsAge(scip, cons) );
      
      /* check, if all variables (including resultant) sum up to an even value */
      odd = FALSE;
      for( i = 0; i < consdata->nvars; ++i )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[i]);
         assert(SCIPisFeasIntegral(scip, solval));
         odd = odd ^ (solval > 0.5);
      }
      if( odd )
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

   assert(separated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *separated = FALSE;

   /* create row for the linear relaxation */
   if( consdata->row == NULL )
   {
      CHECK_OKAY( createRelaxation(scip, cons) );
   }
   assert(consdata->row != NULL);

   /* test row for feasibility and add it, if it is infeasible */
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

/** analyzes conflicting assignment on given constraint, and adds conflict clause to problem */
static
RETCODE analyzeConflict(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< xor constraint that detected the conflict */
   )
{
   CONSDATA* consdata;
   int v;

   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   CHECK_OKAY( SCIPinitConflictAnalysis(scip) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      CHECK_OKAY( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
   }

   /* analyze the conflict */
   CHECK_OKAY( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** propagates constraint with the following rule:
 *   (1) all except one variable fixed  =>  fix remaining variable
 */
static
RETCODE propagateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< xor constraint to be processed */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*             nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   CONSDATA* consdata;
   VAR** vars;
   int nvars;
   int watchedvar1;
   int watchedvar2;
   int i;
   Bool infeasible;
   Bool tightened;
   Bool odd;

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
   CHECK_OKAY( SCIPincConsAge(scip, cons) );

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
      
      debugMessage("constraint <%s>: all vars fixed -> constraint is %s\n",
         SCIPconsGetName(cons), odd ? "infeasible" : "feasible");
      if( odd )
      {
         /* use conflict analysis to get a conflict clause out of the conflicting assignment */
         CHECK_OKAY( analyzeConflict(scip, cons) );

         *cutoff = TRUE;
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      }
      CHECK_OKAY( SCIPdelConsLocal(scip, cons) );

      return SCIP_OKAY;
   }

   /* if only one variable is not fixed, this variable can be deduced */
   if( watchedvar2 == -1 )
   {
      assert(watchedvar1 != -1);
      
      debugMessage("constraint <%s>: only one unfixed variable -> fix <%s> to %d\n",
         SCIPconsGetName(cons), SCIPvarGetName(vars[watchedvar1]), odd);
      CHECK_OKAY( SCIPinferBinvarCons(scip, vars[watchedvar1], odd, cons, PROPRULE_1, &infeasible, &tightened) );
      assert(!infeasible);
      assert(tightened);
      (*nfixedvars)++;
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      CHECK_OKAY( SCIPdelConsLocal(scip, cons) );

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
 *   (1) all except one variable fixed  =>  fix remaining variable
 */
static
RETCODE resolvePropagation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint that inferred the bound change */
   VAR*             infervar,           /**< variable that was deduced */
   PROPRULE         proprule,           /**< propagation rule that deduced the value */
   BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
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
      /* the variable was infered, because all other variables were fixed */
      for( i = 0; i < nvars; ++i )
      {
         if( vars[i] != infervar )
         {
            assert(SCIPisEQ(scip, SCIPvarGetLbAtIndex(vars[i], bdchgidx, FALSE), 
                  SCIPvarGetUbAtIndex(vars[i], bdchgidx, FALSE)));
            CHECK_OKAY( SCIPaddConflictBinvar(scip, vars[i]) );
         }
         else
            assert(SCIPisEQ(scip, SCIPvarGetLbAtIndex(vars[i], bdchgidx, TRUE), 
                  SCIPvarGetUbAtIndex(vars[i], bdchgidx, TRUE)));
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_INVALID:
   default:
      errorMessage("invalid inference information %d in xor constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}





/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
DECL_CONSFREE(consFreeXor)
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
DECL_CONSEXITSOL(consExitsolXor)
{  /*lint --e{715}*/
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
DECL_CONSDELETE(consDeleteXor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   CHECK_OKAY( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransXor)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->nvars >= 1);
   assert(sourcedata->vars != NULL);

   /* create target constraint data */
   CHECK_OKAY( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr,
         sourcedata->nvars-1, &sourcedata->vars[1], sourcedata->vars[0]) );

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpXor)
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
DECL_CONSSEPA(consSepaXor)
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
DECL_CONSENFOLP(consEnfolpXor)
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
DECL_CONSENFOPS(consEnfopsXor)
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
DECL_CONSCHECK(consCheckXor)
{  /*lint --e{715}*/
   Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      CHECK_OKAY( checkCons(scip, conss[i], sol, checklprows, &violated) );
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
DECL_CONSPROP(consPropXor)
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
DECL_CONSPRESOL(consPresolXor)
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

      /* remove all variables that are fixed to zero and all pairs of variables fixed to one */
      CHECK_OKAY( applyFixings(scip, cons, conshdlrdata->eventhdlr) );

      /* propagate constraint */
      CHECK_OKAY( propagateCons(scip, cons, conshdlrdata->eventhdlr, &cutoff, nfixedvars) );

      if( !cutoff && !SCIPconsIsDeleted(cons) && !SCIPconsIsModifiable(cons) )
      {
         assert(consdata->nvars >= 2); /* otherwise, propagateCons() has deleted the constraint */

         /* if only two variables are left, both have to be equal */
         if( consdata->nvars == 2 )
         {
            debugMessage("xor constraint <%s> has only two unfixed variables, remaining vars are even\n",
               SCIPconsGetName(cons));
            
            assert(consdata->vars != NULL);
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[1]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[1]), 1.0));
            
            /* aggregate variables: vars[0] - vars[1] == 0 */
            CHECK_OKAY( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], 1.0, -1.0, 0.0,
                  &cutoff, &redundant, &aggregated) );
            assert(redundant);
            if( aggregated )
               (*naggrvars)++;

            /* delete constraint */
            CHECK_OKAY( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
         }

         /* if only three variables are left, but one of them is fixed to one, the others have to be opposite */
         if( consdata->nvars == 3 )
         {
            assert(consdata->vars != NULL);

            if( SCIPvarGetLbGlobal(consdata->vars[0]) > 0.5 )
            {
               debugMessage("xor constraint <%s> has only two unfixed variables, remaining vars are odd\n",
                  SCIPconsGetName(cons));
            
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[1]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[1]), 1.0));
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[2]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[2]), 1.0));
               
               /* aggregate variables: vars[1] + vars[2] == 1 */
               CHECK_OKAY( SCIPaggregateVars(scip, consdata->vars[1], consdata->vars[2], 1.0, 1.0, 1.0,
                     &cutoff, &redundant, &aggregated) );
               assert(redundant);
               if( aggregated )
                  (*naggrvars)++;

               /* delete constraint */
               CHECK_OKAY( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
            }
            else if( SCIPvarGetLbGlobal(consdata->vars[1]) > 0.5 )
            {
               debugMessage("xor constraint <%s> has only two unfixed variables, remaining vars are odd\n",
                  SCIPconsGetName(cons));
            
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[2]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[2]), 1.0));
               
               /* aggregate variables: vars[0] + vars[2] == 1 */
               CHECK_OKAY( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[2], 1.0, 1.0, 1.0,
                     &cutoff, &redundant, &aggregated) );
               assert(redundant);
               if( aggregated )
                  (*naggrvars)++;

               /* delete constraint */
               CHECK_OKAY( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
            }
            else if( SCIPvarGetLbGlobal(consdata->vars[2]) > 0.5 )
            {
               debugMessage("xor constraint <%s> has only two unfixed variables, remaining vars are odd\n",
                  SCIPconsGetName(cons));
            
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
               assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[1]), 0.0));
               assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[1]), 1.0));
               
               /* aggregate variables: vars[0] + vars[1] == 1 */
               CHECK_OKAY( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], 1.0, 1.0, 1.0,
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
DECL_CONSRESPROP(consRespropXor)
{  /*lint --e{715}*/
   CHECK_OKAY( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockXor)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* external variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      CHECK_OKAY( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   /* internal variable */
   if( consdata->intvar != NULL )
   {
      CHECK_OKAY( SCIPaddVarLocks(scip, consdata->intvar, nlockspos + nlocksneg, nlockspos + nlocksneg) );
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
DECL_CONSPRINT(consPrintXor)
{
   consdataPrint(scip, SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecXor)
{  /*lint --e{715}*/
   CONSDATA* consdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(event != NULL);

   consdata = (CONSDATA*)eventdata;
   assert(consdata != NULL);

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for xor constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrXor(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for events on variables */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecXor,
         NULL) );

   /* create constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeXor, consInitXor, consExitXor, 
         consInitpreXor, consExitpreXor, consInitsolXor, consExitsolXor,
         consDeleteXor, consTransXor, consInitlpXor,
         consSepaXor, consEnfolpXor, consEnfopsXor, consCheckXor, 
         consPropXor, consPresolXor, consRespropXor, consLockXor,
         consActiveXor, consDeactiveXor, 
         consEnableXor, consDisableXor,
         consPrintXor,
         conshdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a xor constraint */
RETCODE SCIPcreateConsXor(
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
   Bool             dynamic,            /**< is constraint subject to aging? */
   Bool             removeable          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;

   /* find the xor constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("xor constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, resvar) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removeable) );

   return SCIP_OKAY;
}
