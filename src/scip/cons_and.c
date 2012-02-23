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

/**@file   cons_and.c
 * @brief  Constraint handler for "and" constraints,  \f$r = x_1 \wedge x_2 \wedge \dots  \wedge x_n\f$
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Michael Winkler
 *
 * This constraint handler deals with "and" constraint. These are constraint of the form:
 *
 * \f[
 *    r = x_1 \wedge x_2 \wedge \dots  \wedge x_n
 * \f]
 *
 * where \f$x_i\f$ is a binary variable for all \f$i\f$. Hence, \f$r\f$ is also of binary type. The variable \f$r\f$ is
 * called resultant and the \f$x\f$'s operators.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_and.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_nonlinear.h"
#include "scip/pub_misc.h"
#include "scip/debug.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "and"
#define CONSHDLR_DESC          "constraint handler for and constraints: r = and(x1, ..., xn)"
#define CONSHDLR_SEPAPRIORITY   +850100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -850100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -850100 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "and"
#define EVENTHDLR_DESC         "bound change event handler for and constraints"

#define DEFAULT_PRESOLPAIRWISE     TRUE /**< should pairwise constraint comparison be performed in presolving? */
#define DEFAULT_LINEARIZE         FALSE /**< should constraint get linearize and removed? */
#define DEFAULT_ENFORCECUTS        TRUE /**< should cuts be separated during LP enforcing? */
#define DEFAULT_AGGRLINEARIZATION FALSE /**< should an aggregated linearization be used? */

#define HASHSIZE_ANDCONS         131101 /**< minimal size of hash table in and constraint tables */
#define DEFAULT_PRESOLUSEHASHING   TRUE /**< should hash table be used for detecting redundant constraints in advance */
#define NMINCOMPARISONS          200000 /**< number for minimal pairwise presolving comparisons */
#define MINGAINPERNMINCOMPARISONS 1e-06 /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise comparison round */
#define EXPRGRAPHREFORM_PRIORITY 100000 /**< priority of expression graph node reformulation method */

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
   int                   nrows;              /**< number of rows for linear relaxation of and constraint */
   int                   watchedvar1;        /**< position of first watched operator variable */
   int                   watchedvar2;        /**< position of second watched operator variable */
   int                   filterpos1;         /**< event filter position of first watched operator variable */
   int                   filterpos2;         /**< event filter position of second watched operator variable */
   unsigned int          propagated:1;       /**< is constraint already preprocessed/propagated? */
   unsigned int          nofixedzero:1;      /**< is none of the operator variables fixed to FALSE? */
   unsigned int          impladded:1;        /**< were the implications of the constraint already added? */
   unsigned int          opimpladded:1;      /**< was the implication for 2 operands with fixed resultant added? */
   unsigned int          sorted:1;           /**< are the constraint's variables sorted? */
   unsigned int          changed:1;          /**< was constraint changed since last pair preprocessing round? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events on watched variables */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Bool             presolusehashing;   /**< should hash table be used for detecting redundant constraints in advance */
   SCIP_Bool             linearize;          /**< should constraint get linearize and removed? */
   SCIP_Bool             enforcecuts;        /**< should cuts be separated during LP enforcing? */
   SCIP_Bool             aggrlinearization;  /**< should an aggregated linearization be used?  */
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

/** creates constraint handler data */
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

/** ensures, that the vars array can store at least num entries */
static
SCIP_RETCODE consdataEnsureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);
   
   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

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
   int v;

   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(resvar != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
   (*consdata)->resvar = resvar;
   (*consdata)->rows = NULL;
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->nrows = 0;
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->filterpos1 = -1;
   (*consdata)->filterpos2 = -1;
   (*consdata)->propagated = FALSE;
   (*consdata)->nofixedzero = FALSE;
   (*consdata)->impladded = FALSE;
   (*consdata)->opimpladded = FALSE;
   (*consdata)->sorted = FALSE;
   (*consdata)->changed = TRUE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->resvar, &(*consdata)->resvar) );

      /* catch needed events on variables */
      SCIP_CALL( consdataCatchEvents(scip, *consdata, eventhdlr) );
   }

   /* capture vars */
   SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->resvar) );
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[v]) );
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
      for( r = 0; r < consdata->nrows; ++r )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rows[r]) );
      }
      SCIPfreeBlockMemoryArray(scip, &consdata->rows, consdata->nrows);
      
      consdata->nrows = 0;
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
   int v;

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

   /* release vars */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[v])) );
   }
   SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->resvar)) );


   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);
 
   return SCIP_OKAY;
}

/** prints and constraint to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< and constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(consdata != NULL);

   /* print resultant */
   SCIP_CALL( SCIPwriteVarName(scip, file, consdata->resvar, TRUE) );

   /* start the variable list */
   SCIPinfoMessage(scip, file, " == and(");

   /* print variable list */
   SCIP_CALL( SCIPwriteVarsList(scip, file, consdata->vars, consdata->nvars, TRUE, ',') );

   /* close the variable list */
   SCIPinfoMessage(scip, file, ")");

   return SCIP_OKAY;
}

/** adds coefficient to and constraint */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows == NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   SCIP_CALL( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1) );
   consdata->vars[consdata->nvars] = var;
   consdata->nvars++;
   consdata->sorted = (consdata->nvars == 1);
   consdata->changed = TRUE;

   /* capture variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   /* if we are in transformed problem, catch the variable's events */
   if( transformed )
   {
      /* catch bound change events of variable */
      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );
   }

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockRounding(scip, cons, var) );

   /**@todo update LP rows */
   if( consdata->rows != NULL )
   {
      SCIPerrorMessage("cannot add coefficients to and constraint after LP relaxation was created\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
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
      /* drop bound change events of variable */
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );
   }

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

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &(consdata->vars[pos])) );

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->nvars--;

   /* if the last variable (that moved) was watched, update the watched position */
   if( consdata->watchedvar1 == consdata->nvars )
      consdata->watchedvar1 = pos;
   if( consdata->watchedvar2 == consdata->nvars )
      consdata->watchedvar2 = pos;

   consdata->propagated = FALSE;
   consdata->sorted = FALSE;
   consdata->changed = TRUE;

   return SCIP_OKAY;
}

/** index comparison method of and constraints: compares two indices of the variable set in the and constraint */
static
SCIP_DECL_SORTINDCOMP(consdataCompVar)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nvars);
   assert(0 <= ind2 && ind2 < consdata->nvars);

   return SCIPvarCompare(consdata->vars[ind1], consdata->vars[ind2]);
}

/** sorts and constraint's variables by non-decreasing variable index */
static
SCIP_RETCODE consdataSort(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->nvars == 0 )
      consdata->sorted = TRUE;
   else if( !consdata->sorted )
   {
      SCIP_VAR* varv;
      int* perm;
      int v;
      int i;
      int nexti;

      /* get temporary memory to store the sorted permutation */
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, consdata->nvars) );

      /* call bubble sort */
      SCIPsort(perm, consdataCompVar, (void*)consdata, consdata->nvars);

      /* permute the variables in the constraint according to the resulting permutation */
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( perm[v] != v )
         {
            SCIP_Bool iswatchedvar1;
            SCIP_Bool iswatchedvar2;

            varv = consdata->vars[v];
            iswatchedvar1 = (consdata->watchedvar1 == v);
            iswatchedvar2 = (consdata->watchedvar2 == v);
            i = v;
            do
            {
               assert(0 <= perm[i] && perm[i] < consdata->nvars);
               assert(perm[i] != i);
               consdata->vars[i] = consdata->vars[perm[i]];
               if( consdata->watchedvar1 == perm[i] )
                  consdata->watchedvar1 = i;
               if( consdata->watchedvar2 == perm[i] )
                  consdata->watchedvar2 = i;
               nexti = perm[i];
               perm[i] = i;
               i = nexti;
            }
            while( perm[i] != v );
            consdata->vars[i] = varv;
            if( iswatchedvar1 )
               consdata->watchedvar1 = i;
            if( iswatchedvar2 )
               consdata->watchedvar2 = i;
            perm[i] = i;
         }
      }
      consdata->sorted = TRUE;

#ifdef SCIP_DEBUG
      /* check sorting */
      for( v = 0; v < consdata->nvars; ++v )
      {
         assert(v == consdata->nvars-1 || SCIPvarCompare(consdata->vars[v], consdata->vars[v+1]) <= 0);
         assert(perm[v] == v);
      }
#endif

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &perm);
   }
   assert(consdata->sorted);

   return SCIP_OKAY;
}

/** deletes all one-fixed variables and removes multiple entries */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< and constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the problem is infeasible */
   int*                  nfixedvars,         /**< pointer to add up the number of found domain reductions */
   int*                  nchgcoefs           /**< pointer to add up the number of changed coefficients */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Bool* contained;
   int nprobvars;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(nchgcoefs != NULL);

   v = 0;
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      assert(SCIPvarIsBinary(var));

      if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         (*nchgcoefs)++;
      }
      else
      {
         SCIP_VAR* repvar;
         SCIP_Bool negated;
         
         /* get binary representative of variable */
         SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &repvar, &negated) );

         /* check, if the variable should be replaced with the representative */
         if( repvar != var )
         {
            /* delete old (aggregated) variable */
            SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );

            /* add representative instead */
            SCIP_CALL( addCoef(scip, cons, eventhdlr, repvar) );
         }
         else
            ++v;
      }
   }

   /* sort the variables in the constraint */
   SCIP_CALL( consdataSort(scip, consdata) );

   /* search for multiple variables; scan from back to front because deletion doesn't affect the order of the front
    * variables
    */
   nprobvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &contained, nprobvars) );
   BMSclearMemoryArray(contained, nprobvars);
   var = NULL;
   for( v = consdata->nvars-1; v >= 0; --v )
   {
      if( consdata->vars[v] == var )
      {
         /* delete the multiple variable */
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         (*nchgcoefs)++;
      }
      else
      {
         SCIP_VAR* probvar;
         int probidx;

         /* we found a new variable */
         var = consdata->vars[v];
         probvar = SCIPvarGetProbvar(var);
         probidx = SCIPvarGetProbindex(probvar);
         assert(0 <= probidx && probidx < nprobvars);
         if( contained[probidx] )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            SCIPdebugMessage("and constraint <%s>: variable <%s> and its negation are present -> fix <%s> = 0\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), SCIPvarGetName(consdata->resvar));

            /* negation of the variable is already present in the constraint: fix resultant to zero */
#ifndef NDEBUG
            {
               int i;
               for( i = consdata->nvars - 1; i > v && var != SCIPvarGetNegatedVar(consdata->vars[i]); --i )
               {}
               assert(i > v);
            }
#endif
            SCIP_CALL( SCIPfixVar(scip, consdata->resvar, 0.0, &infeasible, &fixed) );
            *cutoff = *cutoff && infeasible;
            if( fixed )
               (*nfixedvars)++;

            SCIP_CALL( SCIPdelCons(scip, cons) );
            break;
         }
         contained[probidx] = TRUE;
      }
   }
   SCIPfreeBufferArray(scip, &contained);

   SCIPdebugMessage("after fixings: ");
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL)) );
   SCIPdebugPrintf("\n");

   return SCIP_OKAY;
}

/** creates a linearization of the and constraint */
static 
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons               /**< constraint to check */
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
   consdata->nrows = consdataGetNRows(consdata);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->rows, consdata->nrows) );
   
   /* creates LP rows corresponding to and constraint:
    *   - one additional row:             resvar - v1 - ... - vn >= 1-n
    *   - for each operator variable vi:  resvar - vi            <= 0
    */
   
   /* create additional row */
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_add", SCIPconsGetName(cons));
   SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[0], rowname, 
         -consdata->nvars + 1.0, SCIPinfinity(scip),
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[0], consdata->resvar, 1.0) );
   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[0], nvars, consdata->vars, -1.0) );

   /* create operator rows */
   for( i = 0; i < nvars; ++i )
   {
      (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_%d", SCIPconsGetName(cons), i);
      SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[i+1], rowname, -SCIPinfinity(scip), 0.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[i+1], consdata->resvar, 1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[i+1], consdata->vars[i], -1.0) );
   }
   
   return SCIP_OKAY;
}  

/** adds linear relaxation of and constraint to the LP */
static 
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_ROW* aggrrow;
   SCIP_CONSDATA* consdata;

   char rowname[SCIP_MAXSTRLEN];
   
   /* in the root LP we only add the weaker relaxation which consists of two rows:
    *   - one additional row:             resvar - v1 - ... - vn >= 1-n
    *   - aggregated row:               n*resvar - v1 - ... - vn <= 0.0
    *
    * during separation we separate the stronger relaxation which consists of n+1 row:
    *   - one additional row:             resvar - v1 - ... - vn >= 1-n
    *   - for each operator variable vi:  resvar - vi            <= 0
    */
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   if( consdata->rows == NULL )
   {
      /* create the n+1 row relaxation */
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   
   /* create/add/releas the row aggregated row */
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_operators", SCIPconsGetName(cons));
   SCIP_CALL( SCIPcreateEmptyRow(scip, &aggrrow, rowname, -SCIPinfinity(scip), 0.0,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
   SCIP_CALL( SCIPaddVarToRow(scip, aggrrow, consdata->resvar, (SCIP_Real) consdata->nvars) );
   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, aggrrow, consdata->nvars, consdata->vars, -1.0) );
   SCIP_CALL( SCIPaddCut(scip, NULL, aggrrow, FALSE) );
   SCIP_CALL( SCIPreleaseRow(scip, &aggrrow) );

   /* add additional row */
   if( !SCIProwIsInLP(consdata->rows[0]) )
   {
      SCIP_CALL( SCIPaddCut(scip, NULL, consdata->rows[0], FALSE) );
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
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
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
      
      for( r = 0; r < consdata->nrows; ++r )
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

         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );

            SCIPinfoMessage(scip, NULL, "violation:");
            if( i == consdata->nvars )
            {
               SCIPinfoMessage(scip, NULL, " all operands are TRUE and resultant <%s> = FALSE\n",                   
                  SCIPvarGetName(consdata->resvar)); 
            }
            else
            {
               SCIPinfoMessage(scip, NULL, " operand <%s> = FALSE and resultant <%s> = TRUE\n",
                  SCIPvarGetName(consdata->vars[i-1]), SCIPvarGetName(consdata->resvar)); 
            }
         }
      }
   }
   
   return SCIP_OKAY;
}

/** separates given primal solution */
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

   /* create all necessary rows for the linear relaxation */
   if( consdata->rows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->rows != NULL);

   /* test all rows for feasibility and add infeasible rows */
   for( r = 0; r < consdata->nrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->rows[r]) )
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

/** analyzes conflicting TRUE assignment to resultant of given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflictOne(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< and constraint that detected the conflict */
   int                   falsepos            /**< position of operand that is fixed to FALSE */
   )
{
   SCIP_CONSDATA* consdata;

   /* conflict analysis can only be applied in solving stage and if it turned on */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
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

/** analyzes conflicting FALSE assignment to resultant of given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflictZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< or constraint that detected the conflict */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   assert(!SCIPconsIsModifiable(cons));

   /* conflict analysis can only be applied in solving stage and if it is applicable */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
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
 *
 *  additional if the resultant is fixed to zero during presolving or in the root node (globally), then the "and"
 *  constraint is collapsed to a linear (logicor) constraint of the form 
 *  -> sum_{i=0}^{n-1} ~v_i >= 1
 */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< and constraint to be processed */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars,         /**< pointer to add up the number of found domain reductions */
   int*                  nupgdconss          /**< pointer to add up the number of upgraded constraints */
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
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(resvar), 0.0));
      return SCIP_OKAY;
   }

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   if( !SCIPinRepropagation(scip) )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

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
            /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
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
            /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
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

   /* rules (3) and (4) cannot be applied, if we have at least two unfixed variables left;
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
         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
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
   if( SCIPvarGetUbLocal(resvar) < 0.5 )
   {
      if( watchedvar2 == -1 )
      {
         assert(watchedvar1 != -1);
         
         SCIPdebugMessage("constraint <%s>: resultant <%s> fixed to 0.0, only one unfixed operand -> fix operand <%s> to 0.0\n",
            SCIPconsGetName(cons), SCIPvarGetName(resvar), SCIPvarGetName(vars[watchedvar1]));
         SCIP_CALL( SCIPinferBinvarCons(scip, vars[watchedvar1], FALSE, cons, (int)PROPRULE_4, &infeasible, &tightened) );
         if( infeasible )
         {
            /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
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
      else if( SCIPgetDepth(scip) <= 0 )
      {
         /* since the resultant variable is globally fixed to zero the and constraint collapses to linear constraint of
          * the form 
          * -> \sum_{i=0}^{n-1} v_i <= n-1
          *
          * this can be transformed into a logicor constraint of the form
          * -> \sum_{i=0}^{n-1} ~v_i >= 1
          *
          * create, add, and release the logicor constraint and remove the and constraint globally 
          */
         
         SCIP_VAR** consvars;
         SCIP_CONS* lincons;
         
         assert(SCIPvarGetUbGlobal(resvar) < 0.5);
         
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );
         
         /* collect negated variables */
         for( i = 0; i < nvars; ++i )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, vars[i], &consvars[i]) );
         }

         /* create, add, and release the logicor constraint */
         SCIP_CALL( SCIPcreateConsLogicor(scip, &lincons, SCIPconsGetName(cons), nvars, consvars,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), 
               SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
         SCIP_CALL( SCIPaddCons(scip, lincons) );
         SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

         /* remove the "and" constraint globally */
         SCIP_CALL( SCIPdelCons(scip, cons) );

         (*nupgdconss)++;

         SCIPfreeBufferArray(scip, &consvars);
      }
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
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
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

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyAndcons)
{  /*lint --e{715}*/
   /* the key is the element itself */ 
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqAndcons)
{
   SCIP* scip;
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
   SCIP_Bool coefsequal;
   int i;

   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);
   scip = (SCIP*)userptr; 
   assert(scip != NULL);
   
   /* checks trivial case */
   if( consdata1->nvars != consdata2->nvars )
      return FALSE;

   /* sorts the constraints */
   SCIP_CALL_ABORT( consdataSort(scip, consdata1) );
   SCIP_CALL_ABORT( consdataSort(scip, consdata2) );

   coefsequal = TRUE;

   for( i = 0; i < consdata1->nvars ; ++i )
   {
      /* tests if variables are equal */
      if( consdata1->vars[i] != consdata2->vars[i] )
      {
         assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 1 || 
            SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == -1);
         coefsequal = FALSE;
         break;
      }
      assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 0); 
   } 
   
   return coefsequal;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValAndcons)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   unsigned int hashval;
   int minidx;
   int mididx;
   int maxidx;
   
   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);
   assert(consdata->sorted);
   assert(consdata->nvars > 0);

   minidx = SCIPvarGetIndex(consdata->vars[0]);
   mididx = SCIPvarGetIndex(consdata->vars[consdata->nvars / 2]);
   maxidx = SCIPvarGetIndex(consdata->vars[consdata->nvars - 1]);
   assert(minidx >= 0 && minidx <= maxidx);

   hashval = (consdata->nvars << 29) + (minidx << 22) + (mididx << 11) + maxidx; /*lint !e701*/

   return hashval;
}

/** updates the flags of the first constraint according to the ones of the second constraint */
static
SCIP_RETCODE updateFlags(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< constraint that should stay */
   SCIP_CONS*            cons1               /**< constraint that should be deleted */
   )
{
   if( SCIPconsIsInitial(cons1) )
   {
      SCIP_CALL( SCIPsetConsInitial(scip, cons0, TRUE) );
   }
   if( SCIPconsIsSeparated(cons1) )
   {
      SCIP_CALL( SCIPsetConsSeparated(scip, cons0, TRUE) );
   }
   if( SCIPconsIsEnforced(cons1) )
   {
      SCIP_CALL( SCIPsetConsEnforced(scip, cons0, TRUE) );
   }
   if( SCIPconsIsChecked(cons1) )
   {
      SCIP_CALL( SCIPsetConsChecked(scip, cons0, TRUE) );
   }
   if( SCIPconsIsPropagated(cons1) )
   {
      SCIP_CALL( SCIPsetConsPropagated(scip, cons0, TRUE) );
   }
   if( !SCIPconsIsDynamic(cons1) )
   {
      SCIP_CALL( SCIPsetConsDynamic(scip, cons0, FALSE) );
   }
   if( !SCIPconsIsRemovable(cons1) )
   {
      SCIP_CALL( SCIPsetConsRemovable(scip, cons0, FALSE) );
   }
   if( SCIPconsIsStickingAtNode(cons1) )
   {
      SCIP_CALL( SCIPsetConsStickingAtNode(scip, cons0, TRUE) );
   }

   return SCIP_OKAY;
}


/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint 
 *  accordingly; in contrast to removeRedundantConstraints(), it uses a hash table 
 */
static
SCIP_RETCODE detectRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   int*                  firstchange,        /**< pointer to store first changed constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
)
{
   SCIP_HASHTABLE* hashtable;
   int hashtablesize;
   int c;

   assert(conss != NULL);
   assert(ndelconss != NULL);

   /* create a hash table for the constraint set */
   hashtablesize = SCIPcalcHashtableSize(10*nconss);
   hashtablesize = MAX(hashtablesize, HASHSIZE_ANDCONS);
   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashGetKeyAndcons, hashKeyEqAndcons, hashKeyValAndcons, (void*) scip) );

   *cutoff = FALSE;

   /* check all constraints in the given set for redundancy */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons0;
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata0;

      cons0 = conss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;

      consdata0 = SCIPconsGetData(cons0);
      /* sort the constraint */
      SCIP_CALL( consdataSort(scip, consdata0) );

      /* get constraint from current hash table with same variables as cons0 */
      cons1 = (SCIP_CONS*)(SCIPhashtableRetrieve(hashtable, (void*)cons0));
 
      if( cons1 != NULL )
      {
         SCIP_CONSDATA* consdata1;
         SCIP_Bool redundant;


         assert(SCIPconsIsActive(cons1));
         assert(!SCIPconsIsModifiable(cons1));
      
         consdata1 = SCIPconsGetData(cons1);
         
         assert(consdata0 != NULL && consdata1 != NULL);
         assert(consdata0->nvars >= 1 && consdata0->nvars == consdata1->nvars);
         
         assert(consdata0->sorted && consdata1->sorted);
         assert(consdata0->vars[0] == consdata1->vars[0]);

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( updateFlags(scip, cons1, cons0) ); 
         redundant = FALSE;

         if( consdata0->resvar != consdata1->resvar )
         {
            SCIP_Bool aggregated;
            
            assert(SCIPvarCompare(consdata0->resvar, consdata1->resvar) != 0); 
         
            /* aggregate resultants */
            SCIP_CALL( SCIPaggregateVars(scip, consdata0->resvar, consdata1->resvar, 1.0, -1.0, 0.0,
                  cutoff, &redundant, &aggregated) );
            assert(redundant || SCIPdoNotAggr(scip));

            if( aggregated )
               ++(*naggrvars);
            if( *cutoff )
               goto TERMINATE;
         }
         else 
            redundant = TRUE;

         /* delete consdel */
         if( redundant )
         {
            SCIP_CALL( SCIPdelCons(scip, cons0) );
            (*ndelconss)++;
         }

         /* update the first changed constraint to begin the next aggregation round with */
         if( consdata0->changed && SCIPconsGetPos(cons1) < *firstchange )
            *firstchange = SCIPconsGetPos(cons1);

         assert(SCIPconsIsActive(cons1));
      }
      else
      {
         /* no such constraint in current hash table: insert cons0 into hash table */  
         SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) cons0) );
      }
   }
 TERMINATE:
   /* free hash table */
   SCIPhashtableFree(&hashtable);

   return SCIP_OKAY;
}


/** compares constraint with all prior constraints for possible redundancy or aggregation,
 *  and removes or changes constraint accordingly
 */
static
SCIP_RETCODE preprocessConstraintPairs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   firstchange,        /**< first constraint that changed since last pair preprocessing round */
   int                   chkind,             /**< index of constraint to check against all prior indices upto startind */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  nbdchgs,            /**< pointer to count the number of performed bound changes, or NULL */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   SCIP_Bool cons0changed;
   int c;

   assert(conss != NULL);
   assert(firstchange <= chkind);
   assert(cutoff != NULL);
   assert(naggrvars != NULL);
   assert(nbdchgs != NULL);
   assert(ndelconss != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);

   /* sort the constraint */
   SCIP_CALL( consdataSort(scip, consdata0) );

   /* check constraint against all prior constraints */
   cons0changed = consdata0->changed;
   
   if( SCIPconsIsActive(cons0) )
   {
      for( c = (cons0changed ? 0 : firstchange); c < chkind && !(*cutoff) && !SCIPisStopped(scip); ++c )
      {
         SCIP_CONS* cons1;
         SCIP_CONSDATA* consdata1;
         SCIP_Bool cons0superset;
         SCIP_Bool cons1superset;
         int v0;
         int v1;

         cons1 = conss[c];

         /* ignore inactive and modifiable constraints */
         if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
            continue;

         consdata1 = SCIPconsGetData(cons1);
         assert(consdata1 != NULL);

#if 0
         SCIPdebugMessage("preprocess and constraint pair <%s>[chg:%d] and <%s>[chg:%d]\n",
            SCIPconsGetName(cons0), cons0changed, SCIPconsGetName(cons1), consdata1->changed);
#endif

         /* if both constraints were not changed since last round, we can ignore the pair */
         if( !cons0changed && !consdata1->changed )
            continue;

         assert(consdata1->nvars >= 1);

         /* sort the constraint */
         SCIP_CALL( consdataSort(scip, consdata1) );

         /* check consdata0 against consdata1:
          * - if they consist of the same operands, the resultants can be aggregated
          * - if one operand list is a subset of the other, add implication r0 = 1 -> r1 = 1, or r1 = 1 -> r0 = 1
          */
         v0 = 0;
         v1 = 0;
         cons0superset = TRUE;
         cons1superset = TRUE;
         while( (v0 < consdata0->nvars || v1 < consdata1->nvars) && (cons0superset || cons1superset) )
         {
            int varcmp;

            /* test, if variable appears in only one or in both constraints */
            if( v0 < consdata0->nvars && v1 < consdata1->nvars )
               varcmp = SCIPvarCompare(consdata0->vars[v0], consdata1->vars[v1]);
            else if( v0 < consdata0->nvars )
               varcmp = -1;
            else
               varcmp = +1;

            switch( varcmp )
            {
            case -1:
               /* variable doesn't appear in consdata1 */
               cons1superset = FALSE;
               v0++;
               break;

            case +1:
               /* variable doesn't appear in consdata0 */
               cons0superset = FALSE;
               v1++;
               break;

            case 0:
               /* variable appears in both constraints */
               v0++;
               v1++;
               break;

            default:
               SCIPerrorMessage("invalid comparison result\n");
               SCIPABORT();
            }
         }

         /* check for equivalence and domination */
         if( cons0superset && cons1superset )
         {
            SCIP_Bool infeasible;
            SCIP_Bool redundant;
            SCIP_Bool aggregated;

            /* constraints are equivalent */
            SCIPdebugMessage("equivalent and constraints <%s> and <%s>: aggregate resultants <%s> == <%s>\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), SCIPvarGetName(consdata0->resvar),
               SCIPvarGetName(consdata1->resvar));
         
            /* aggregate resultants */
            SCIP_CALL( SCIPaggregateVars(scip, consdata0->resvar, consdata1->resvar, 1.0, -1.0, 0.0,
                  &infeasible, &redundant, &aggregated) );
            assert(redundant || SCIPdoNotAggr(scip));

            if( aggregated )
            {
               assert(redundant);
               (*naggrvars)++;
            }
            
            if( redundant )
            {
               /* delete constraint */
               SCIP_CALL( SCIPdelCons(scip, cons1) );
               (*ndelconss)++;
            }

            *cutoff = *cutoff || infeasible;
         }
         else if( cons0superset )
         {
            SCIP_Bool infeasible;
            int nboundchgs;

            /* the conjunction of cons0 is a superset of the conjunction of cons1 */
            SCIPdebugMessage("and constraint <%s> is superset of <%s>: add implication <%s> = 1 -> <%s> = 1\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), SCIPvarGetName(consdata0->resvar),
               SCIPvarGetName(consdata1->resvar));

            /* add implication */
            SCIP_CALL( SCIPaddVarImplication(scip, consdata0->resvar, TRUE, consdata1->resvar, SCIP_BOUNDTYPE_LOWER, 1.0,
                  &infeasible, &nboundchgs) );
            *cutoff = *cutoff || infeasible;
            (*nbdchgs) += nboundchgs;
         }
         else if( cons1superset )
         {
            SCIP_Bool infeasible;
            int nboundchgs;

            /* the conjunction of cons1 is a superset of the conjunction of cons0 */
            SCIPdebugMessage("and constraint <%s> is superset of <%s>: add implication <%s> = 1 -> <%s> = 1\n",
               SCIPconsGetName(cons1), SCIPconsGetName(cons0), SCIPvarGetName(consdata1->resvar),
               SCIPvarGetName(consdata0->resvar));

            /* add implication */
            SCIP_CALL( SCIPaddVarImplication(scip, consdata1->resvar, TRUE, consdata0->resvar, SCIP_BOUNDTYPE_LOWER, 1.0,
                  &infeasible, &nboundchgs) );
            *cutoff = *cutoff || infeasible;
            (*nbdchgs) += nboundchgs;
         }
      }
   }
   consdata0->changed = FALSE;

   return SCIP_OKAY;
}

/** tries to reformulate an expression graph node that is a product of binary variables via introducing an and constraint */
static
SCIP_DECL_EXPRGRAPHNODEREFORM(exprgraphnodeReformAnd)
{
   SCIP_EXPRGRAPHNODE* child;
   char name[SCIP_MAXSTRLEN];
   int nchildren;
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   int c;

   assert(scip != NULL);
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(naddcons != NULL);
   assert(reformnode != NULL);

   *reformnode = NULL;

   /* allow only products given as EXPR_PRODUCT or EXPR_POLYNOMIAL with only 1 monomial */
   if( SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_PRODUCT &&
       (SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_POLYNOMIAL || SCIPexprgraphGetNodePolynomialNMonomials(node) > 1)
     )
      return SCIP_OKAY;

   nchildren = SCIPexprgraphGetNodeNChildren(node);

   /* for a polynomial with only one monomial, all children should appear as factors in the monomial
    * since we assume that the factors have been merged, this means that the number of factors in the monomial should equal the number of children of the node
    */
   assert(SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_POLYNOMIAL || SCIPexprGetMonomialNFactors(SCIPexprgraphGetNodePolynomialMonomials(node)[0]) == nchildren);

   /* check only products with at least 3 variables (2 variables are taken of by cons_quadratic) */
   if( nchildren <= 2 )
      return SCIP_OKAY;

   /* check if all factors correspond to binary variables, and if so, setup vars array */
   for( c = 0; c < nchildren; ++c )
   {
      child = SCIPexprgraphGetNodeChildren(node)[c];

      if( SCIPexprgraphGetNodeOperator(child) != SCIP_EXPR_VARIDX )
         return SCIP_OKAY;

      var = (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, child);
      if( !SCIPvarIsBinary(var) )
         return SCIP_OKAY;
   }

   /* node corresponds to product of binary variables (maybe with coefficient and constant, if polynomial) */
   SCIPdebugMessage("reformulate node %p via and constraint\n", (void*)node);

   /* collect variables in product */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nchildren) );
   for( c = 0; c < nchildren; ++c )
   {
      child = SCIPexprgraphGetNodeChildren(node)[c];
      vars[c] = (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, child);
   }

   /* create variable for resultant
    * cons_and wants to add implications for resultant, which is only possible for binary variables currently
    * so choose binary as vartype, even though implicit integer had been sufficient
    */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%dand", *naddcons);
   SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
      TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, var) );

#ifdef SCIP_DEBUG_SOLUTION
   {
      SCIP_Bool debugval;
      SCIP_Real varval;

      debugval = TRUE;
      for( c = 0; c < nchildren; ++c )
      {
         SCIP_CALL( SCIPdebugGetSolVal(scip, vars[c], &varval) );
         debugval = debugval && (varval > 0.5);
      }
      SCIP_CALL( SCIPdebugAddSolVal(scip, var, debugval ? 1.0 : 0.0) );
   }
#endif

   /* create and constraint */
   SCIP_CALL( SCIPcreateConsAnd(scip, &cons, name, var, nchildren, vars,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   ++*naddcons;

   SCIPfreeBufferArray(scip, &vars);

   /* add var to exprgraph */
   SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&var, reformnode) );
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   /* if we have coefficient and constant, then replace reformnode by linear expression in reformnode */
   if( SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_POLYNOMIAL )
   {
      SCIP_Real coef;
      SCIP_Real constant;

      coef = SCIPexprGetMonomialCoef(SCIPexprgraphGetNodePolynomialMonomials(node)[0]);
      constant = SCIPexprgraphGetNodePolynomialConstant(node);

      if( coef != 1.0 || constant != 0.0 )
      {
         SCIP_EXPRGRAPHNODE* linnode;
         SCIP_CALL( SCIPexprgraphCreateNodeLinear(SCIPblkmem(scip), &linnode, 1, &coef, constant) );
         SCIP_CALL( SCIPexprgraphAddNode(exprgraph, linnode, -1, 1, reformnode) );
         *reformnode = linnode;
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyAnd)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrAnd(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}

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
static
SCIP_DECL_CONSINITPRE(consInitpreAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( nconss == 0 || conss != NULL );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   if( conshdlrdata->linearize )
   {
      /* linearize all "and" constraints  and remove the "and" constraints */
      SCIP_CONS* newcons;
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      char consname[SCIP_MAXSTRLEN];
      
      SCIP_VAR** vars;
      SCIP_Real* vals;
      
      int nvars;
      int c, v;
      
      /* alloc buffer array */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );
      
      for( c = 0; c < nconss; ++c )
      {
         cons = conss[c];
         assert( cons != NULL );
         
         consdata = SCIPconsGetData(cons);
         assert( consdata != NULL );
         assert( consdata->resvar != NULL );
         
         nvars = consdata->nvars;
         
         if( !conshdlrdata->aggrlinearization )
         {
            vars[0] = consdata->resvar;
            vals[0] = 1.0;
            vals[1] = -1.0;
         
            /* create operator linear constraints */
            for( v = 0; v < nvars; ++v )
            {
               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_%d", SCIPconsGetName(cons), v);
               vars[1] = consdata->vars[v];
            
               SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, consname, 2, vars, vals, -SCIPinfinity(scip), 0.0,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), 
                     SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                     SCIPconsIsStickingAtNode(cons)) );
               
               /* add constraint */
               SCIP_CALL( SCIPaddCons(scip, newcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            }
         }
         
         /* realloc  buffer array */
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, nvars + 1) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &vals, nvars + 1) );

         for( v = 0; v < nvars; ++v )
         {
            vars[v] = consdata->vars[v];
            vals[v] = -1.0;
         }
         
         vars[nvars] = consdata->resvar;

         if( conshdlrdata->aggrlinearization )
         {
            /* create additional linear constraint */
            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_operators", SCIPconsGetName(cons));

            vals[nvars] = (SCIP_Real) nvars;

            SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, consname, nvars + 1, vars, vals, -SCIPinfinity(scip), 0.0,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), 
                  SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                  SCIPconsIsStickingAtNode(cons)) );

            /* add constraint */
            SCIP_CALL( SCIPaddCons(scip, newcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         }

         /* create additional linear constraint */
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_add", SCIPconsGetName(cons));

         vals[nvars] = 1.0;

         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, consname, nvars + 1, vars, vals, -nvars + 1.0, SCIPinfinity(scip),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), 
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
               SCIPconsIsStickingAtNode(cons)) );
         
         /* add constraint */
         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         
         /* delete constraint */
         SCIP_CALL( SCIPdelCons(scip, cons) );
      }
      
      /* free buffer array */
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &vals);
   }
   
   return SCIP_OKAY;
}



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
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpAnd)
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
SCIP_DECL_CONSSEPALP(consSepalpAnd)
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
SCIP_DECL_CONSSEPASOL(consSepasolAnd)
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
SCIP_DECL_CONSENFOLP(consEnfolpAnd)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool separated;
   SCIP_Bool violated;
   int i;

   separated = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, FALSE, FALSE, &violated) );
      if( violated )
      {
         if( conshdlrdata->enforcecuts )
         {
            
            SCIP_CALL( separateCons(scip, conss[i], NULL, &separated) );
            assert(separated); /* because the solution is integral, the separation always finds a cut */
         }
         else
         {
            *result = SCIP_INFEASIBLE;
            return SCIP_OKAY;
         }
      }
   } 

   if( separated )
      *result = SCIP_SEPARATED;
   else
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
      SCIP_CALL( checkCons(scip, conss[i], NULL, TRUE, FALSE, &violated) );
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
      SCIP_CALL( checkCons(scip, conss[i], sol, checklprows, printreason, &violated) );
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
   int nupgdconss;
   int c;
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   nfixedvars = 0;
   nupgdconss = 0;

   /* propagate all useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( propagateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &nfixedvars, &nupgdconss) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 || nupgdconss > 0 )
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
   SCIP_Bool delay;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnupgdconss;
   int firstchange;
   int c;

   assert(result != NULL);

   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   cutoff = FALSE;
   delay = FALSE;
   firstchange = INT_MAX;
   for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->propagated = FALSE;

      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == INT_MAX && consdata->changed )
         firstchange = c;

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons, conshdlrdata->eventhdlr, &cutoff, nfixedvars, nupgdconss) );

      /* remove all variables that are fixed to one; merge multiple entries of the same variable;
       * fix resultant to zero if a pair of negated variables is contained in the operand variables
       */
      if( !cutoff && !SCIPconsIsDeleted(cons) )
      {
         SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr, &cutoff, nfixedvars, nchgcoefs) );
      }

      if( !cutoff && !SCIPconsIsDeleted(cons) && !SCIPconsIsModifiable(cons) )
      {
         assert(consdata->nvars >= 1); /* otherwise, propagateCons() has deleted the constraint */

         /* if only one variable is left, the resultant has to be equal to this single variable */
         if( consdata->nvars == 1 )
         {
            SCIP_Bool redundant;
            SCIP_Bool aggregated;

            SCIPdebugMessage("and constraint <%s> has only one variable not fixed to 1.0\n", SCIPconsGetName(cons));
            
            assert(consdata->vars != NULL);
            assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
            
            /* aggregate variables: resultant - operand == 0 */
            SCIP_CALL( SCIPaggregateVars(scip, consdata->resvar, consdata->vars[0], 1.0, -1.0, 0.0,
                  &cutoff, &redundant, &aggregated) );
            assert(redundant || SCIPdoNotAggr(scip));

            if( aggregated )
            {
               assert(redundant);
               (*naggrvars)++;
            }
            
            if( redundant )
            {
               /* delete constraint */
               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
            }
         }
         else if( !consdata->impladded )
         {
            int i;

            /* add implications: resultant == 1 -> all operands == 1 */
            for( i = 0; i < consdata->nvars && !cutoff; ++i )
            {
               int nimplbdchgs;

               SCIP_CALL( SCIPaddVarImplication(scip, consdata->resvar, TRUE, consdata->vars[i],
                     SCIP_BOUNDTYPE_LOWER, 1.0, &cutoff, &nimplbdchgs) );
               (*nchgbds) += nimplbdchgs;
            }
            consdata->impladded = TRUE;
         }

         /* if in r = x and y, the resultant is fixed to zero, add implication x = 1 -> y = 0 */
         if( !cutoff && SCIPconsIsActive(cons) && consdata->nvars == 2 && !consdata->opimpladded
            && SCIPvarGetUbGlobal(consdata->resvar) < 0.5 )
         {
            int nimplbdchgs;
            
            SCIP_CALL( SCIPaddVarImplication(scip, consdata->vars[0], TRUE, consdata->vars[1],
                  SCIP_BOUNDTYPE_UPPER, 0.0, &cutoff, &nimplbdchgs) );
            (*nchgbds) += nimplbdchgs;
            consdata->opimpladded = TRUE;
         }
      }
   }

   /* process pairs of constraints: check them for equal operands in order to aggregate resultants;
    * only apply this expensive procedure, if the single constraint preprocessing did not find any reductions
    * (otherwise, we delay the presolving to be called again next time)
    */
   if( !cutoff && conshdlrdata->presolusehashing )
   {
      if( *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars )
      {
         if( firstchange < nconss ) 
         {
            /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
            SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &firstchange, &cutoff, naggrvars, ndelconss) );
            oldnaggrvars = *naggrvars;
         }
      }
      else
         delay = TRUE;
   }

   if( !cutoff && conshdlrdata->presolpairwise )
   {
      if( *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars )
      {
         SCIP_Longint npaircomparisons;
         npaircomparisons = 0;
         oldndelconss = *ndelconss;
         
         for( c = firstchange; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
         {
            if( SCIPconsIsActive(conss[c]) && !SCIPconsIsModifiable(conss[c]) )
            {
               npaircomparisons += ((SCIPconsGetData(conss[c])->changed) ? (SCIP_Longint) c : ((SCIP_Longint) c - (SCIP_Longint) firstchange));
               
               SCIP_CALL( preprocessConstraintPairs(scip, conss, firstchange, c,
                                                    &cutoff, naggrvars, nchgbds, ndelconss) );
               
               if( npaircomparisons > NMINCOMPARISONS )
               {
                  if( ((*ndelconss - oldndelconss) + (*naggrvars - oldnaggrvars) + (*nchgbds - oldnchgbds)/2.0) / ((SCIP_Real) npaircomparisons) < MINGAINPERNMINCOMPARISONS )
                     break;
                  oldndelconss = *ndelconss;
                  oldnaggrvars = *naggrvars;
                  oldnchgbds = *nchgbds;

                  npaircomparisons = 0;
               }
            }
         }
      }
      else
         delay = TRUE;
   }

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( delay )
      *result = SCIP_DELAYED;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nchgbds > oldnchgbds
            || *ndelconss > oldndelconss || *nupgdconss > oldnupgdconss )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

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


/** variable deletion method of constraint handler */
#define consDelvarsAnd NULL


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintAnd)
{  /*lint --e{715}*/
   
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file) );
      
   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyAnd)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   SCIP_VAR** vars;
   SCIP_VAR* sourceresvar;
   SCIP_VAR* resvar;
   const char* consname;
   int nvars;
   int v;

   assert(valid != NULL);
   (*valid) = TRUE;
   
   sourceresvar = SCIPgetResultantAnd(sourcescip, sourcecons);

   /* map resultant to active variable of the target SCIP  */
   SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourceresvar, &resvar, varmap, consmap, global, valid) );
   assert(!(*valid) || resvar != NULL);

   /* we do not copy, if a variable is missing */
   if( !(*valid) )
      return SCIP_OKAY;

   /* map operand variables to active variables of the target SCIP  */
   sourcevars = SCIPgetVarsAnd(sourcescip, sourcecons);
   nvars = SCIPgetNVarsAnd(sourcescip, sourcecons);

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &vars[v], varmap, consmap, global, valid) );
      assert(!(*valid) || vars[v] != NULL);

      /* we do not copy, if a variable is missing */
      if( !(*valid) )
         goto TERMINATE;
   }
   
   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);
 
   /* creates and captures a and constraint */
   SCIP_CALL( SCIPcreateConsAnd(scip, cons, consname, resvar, nvars, vars, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

 TERMINATE:   
   /* free buffer array */
   SCIPfreeBufferArray(scip, &vars);
   
   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseAnd)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_VAR* resvar;
   char* endptr;
   int requiredsize;
   int varssize;
   int nvars;
   
   SCIPdebugMessage("parse <%s> as and constraint\n", str);

   /* parse variable name */ 
   SCIP_CALL( SCIPparseVarName(scip, str, &resvar, &endptr) );
   str = endptr;

   if( resvar == NULL )
   {
      SCIPdebugMessage("resultant variable does not exist \n");
      *success = FALSE;
   }
   else
   {
      char* strcopy;
      char* token;
      char* saveptr;

      /* copy string for truncating it */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &strcopy, str, (int)(strlen(str)+1)));

      /* cutoff "== and(" form the constraint string */
      (void) SCIPstrtok(strcopy, "(", &saveptr );

      /* cutoff ")" form the constraint string */
      token = SCIPstrtok(NULL, ")", &saveptr );

      varssize = 100;
      nvars = 0;

      /* allocate buffer array for variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );

      /* parse string */
      SCIP_CALL( SCIPparseVarsList(scip, token, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );
      token = endptr;

      if( *success )
      {
         /* check if the size of the variable array was great enough */
         if( varssize < requiredsize )
         {
            /* reallocate memory */
            varssize = requiredsize;
            SCIP_CALL( SCIPreallocBufferArray(scip, &vars, varssize) );
            
            /* parse string again with the correct size of the variable array */
            SCIP_CALL( SCIPparseVarsList(scip, token, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );
         }
         
         assert(*success);
         assert(varssize >= requiredsize);
         
         /* create and constraint */
         SCIP_CALL( SCIPcreateConsAnd(scip, cons, name, resvar, nvars, vars, 
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
      }

      /* free variable buffer */
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &strcopy);
   }
   
   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsAnd)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars + 1 )
      (*success) = FALSE;
   else
   {
      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      vars[consdata->nvars] = consdata->resvar;
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variable (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsAnd)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars + 1;
   (*success) = TRUE;

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
   SCIP*                 scip               /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create event handler for events on variables */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecAnd,
         NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyAnd,
         consFreeAnd, consInitAnd, consExitAnd,
         consInitpreAnd, consExitpreAnd, consInitsolAnd, consExitsolAnd,
         consDeleteAnd, consTransAnd, consInitlpAnd,
         consSepalpAnd, consSepasolAnd, consEnfolpAnd, consEnfopsAnd, consCheckAnd,
         consPropAnd, consPresolAnd, consRespropAnd, consLockAnd,
         consActiveAnd, consDeactiveAnd,
         consEnableAnd, consDisableAnd,
         consDelvarsAnd, consPrintAnd, consCopyAnd, consParseAnd,
         consGetVarsAnd, consGetNVarsAnd, conshdlrdata) );

   /* add and constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/and/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance",
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/linearize",
         "should the \"and\" constraint get linearized and removed (in presolving)?",
         &conshdlrdata->linearize, TRUE, DEFAULT_LINEARIZE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/enforcecuts",
         "should cuts be separated during LP enforcing?",
         &conshdlrdata->enforcecuts, TRUE, DEFAULT_ENFORCECUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/aggrlinearization",
         "should an aggregated linearization be used?",
         &conshdlrdata->aggrlinearization, TRUE, DEFAULT_AGGRLINEARIZATION, NULL, NULL) );

   if( SCIPfindConshdlr(scip, "nonlinear") != NULL )
   {
      /* include the and-constraint upgrade in the nonlinear constraint handler */
      SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, NULL, exprgraphnodeReformAnd, EXPRGRAPHREFORM_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   return SCIP_OKAY;
}

/** creates and captures a and constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsAnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             resvar,             /**< resultant variable of the operation */
   int                   nvars,              /**< number of operator variables in the constraint */
   SCIP_VAR**            vars,               /**< array with operator variables of constraint */
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
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** gets number of variables in and constraint */
int SCIPgetNVarsAnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an and constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets array of variables in and constraint */
SCIP_VAR** SCIPgetVarsAnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an and constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}


/** gets the resultant variable in and constraint */
SCIP_VAR* SCIPgetResultantAnd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;
   
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an and constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   return consdata->resvar;
}
