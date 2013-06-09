/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_logicor.c
 * @brief  Constraint handler for logic or constraints \f$1^T x \ge 1\f$
 *         (equivalent to set covering, but algorithms are suited for depth first search).
 * @author Tobias Achterberg
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_logicor.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/pub_misc.h"


#define CONSHDLR_NAME          "logicor"
#define CONSHDLR_DESC          "logic or constraints"
#define CONSHDLR_SEPAPRIORITY    +10000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -2000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -2000000 /**< priority of the constraint handler for checking feasibility */
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

#define LINCONSUPGD_PRIORITY    +800000 /**< priority of the constraint handler for upgrading of linear constraints */

#define EVENTHDLR_NAME         "logicor"
#define EVENTHDLR_DESC         "event handler for logic or constraints"

#define CONFLICTHDLR_NAME      "logicor"
#define CONFLICTHDLR_DESC      "conflict handler creating logic or constraints"
#define CONFLICTHDLR_PRIORITY  LINCONSUPGD_PRIORITY

#define DEFAULT_PRESOLPAIRWISE     TRUE /**< should pairwise constraint comparison be performed in presolving? */

#define HASHSIZE_LOGICORCONS     131101 /**< minimal size of hash table in logicor constraint tables */
#define DEFAULT_PRESOLUSEHASHING   TRUE /**< should hash table be used for detecting redundant constraints in advance */
#define NMINCOMPARISONS          200000 /**< number for minimal pairwise presolving comparisons */
#define MINGAINPERNMINCOMPARISONS 1e-06 /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise comparison round */
#define DEFAULT_DUALPRESOLVING     TRUE /**< should dual presolving steps be performed? */
#define DEFAULT_NEGATEDCLIQUE      TRUE /**< should negated clique information be used in presolving */

/* @todo make this a parameter setting */
#if 1 /* @todo test which AGEINCREASE formula is better! */
#define AGEINCREASE(n) (1.0 + 0.2*n)
#else
#define AGEINCREASE(n) (0.1*n)
#endif


/* @todo maybe use event SCIP_EVENTTYPE_VARUNLOCKED to decide for another dual-presolving run on a constraint */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for events on watched variables */
   SCIP_CONSHDLR*        conshdlrlinear;     /**< pointer to linear constraint handler or NULL if not included */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Bool             presolusehashing;   /**< should hash table be used for detecting redundant constraints in advance */
   SCIP_Bool             dualpresolving;     /**< should dual presolving steps be performed? */
   SCIP_Bool             usenegatedclique;   /**< should negated clique information be used in presolving */
   int                   nlastcliques;       /**< number of cliques after last negated clique presolving round */
   int                   nlastimpls;         /**< number of implications after last negated clique presolving round */
};

/** logic or constraint data */
struct SCIP_ConsData
{
   SCIP_ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   SCIP_VAR**            vars;               /**< variables of the constraint */
   int                   varssize;           /**< size of vars array */
   int                   nvars;              /**< number of variables in the constraint */
   int                   watchedvar1;        /**< position of the first watched variable */
   int                   watchedvar2;        /**< position of the second watched variable */
   int                   filterpos1;         /**< event filter position of first watched variable */
   int                   filterpos2;         /**< event filter position of second watched variable */
   unsigned int          impladded:1;        /**< was the 2-variable logic or constraint already added as implication? */
   unsigned int          sorted:1;           /**< are the constraint's variables sorted? */
   unsigned int          changed:1;          /**< was constraint changed since last redundancy round in preprocessing? */
   unsigned int          merged:1;           /**< are the constraint's equal/negated variables already merged? */
};


/*
 * Local methods
 */

/** installs rounding locks for the given variable in the given logic or constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding down may violate the constraint */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, FALSE) );

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given logic or constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding down may violate the constraint */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, FALSE) );

   return SCIP_OKAY;
}

/** creates constraint handler data for logic or constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata,       /**< pointer to store the constraint handler data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   (*conshdlrdata)->nlastcliques = 0;
   (*conshdlrdata)->nlastimpls = 0;

   /* set event handler for catching events on watched variables */
   (*conshdlrdata)->eventhdlr = eventhdlr;

   return SCIP_OKAY;
}

/** frees constraint handler data for logic or constraint handler */
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

/** ensures, that the vars array can store at least num entries */
static
SCIP_RETCODE consdataEnsureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< logicor constraint data */
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

/** creates a logic or constraint data object */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the logic or constraint data */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars                /**< variables of the constraint */
   )
{
   int v;

   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->row = NULL;
   if( nvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
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
   (*consdata)->filterpos1 = -1;
   (*consdata)->filterpos2 = -1;
   (*consdata)->impladded = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->sorted = (nvars <= 1);
   (*consdata)->merged = (nvars <= 1);

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
   }

   /* capture variables */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[v]) );
   }

   return SCIP_OKAY;
}   

/** frees a logic or constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to the logic or constraint */
   )
{
   int v;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* release variables */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[v])) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints logic or constraint to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< logic or constraint data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             endline             /**< should an endline be set? */
   )
{
   assert(consdata != NULL);

   /* print constraint type */
   SCIPinfoMessage(scip, file, "logicor(");

   /* print variable list */
   SCIP_CALL( SCIPwriteVarsList(scip, file, consdata->vars, consdata->nvars, TRUE, ',') );
   
   /* close bracket */
   SCIPinfoMessage(scip, file, ")");
   
   if( endline )
      SCIPinfoMessage(scip, file, "\n");

   return SCIP_OKAY;
}

/** stores the given variable numbers as watched variables, and updates the event processing */
static
SCIP_RETCODE switchWatchedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   watchedvar1,        /**< new first watched variable */
   int                   watchedvar2         /**< new second watched variable */
   )
{
   SCIP_CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
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
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar1],
            SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED, eventhdlr, (SCIP_EVENTDATA*)cons,
            consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 && consdata->watchedvar2 != watchedvar2 )
   {
      assert(consdata->filterpos2 != -1);
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar2],
            SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED, eventhdlr, (SCIP_EVENTDATA*)cons, 
            consdata->filterpos2) );
   }

   /* catch events on new watched variables */
   if( watchedvar1 != -1 && watchedvar1 != consdata->watchedvar1 )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[watchedvar1],
            SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED, eventhdlr, (SCIP_EVENTDATA*)cons,
            &consdata->filterpos1) );
   }
   if( watchedvar2 != -1 && watchedvar2 != consdata->watchedvar2 )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[watchedvar2],
            SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED, eventhdlr, (SCIP_EVENTDATA*)cons,
            &consdata->filterpos2) );
   }

   /* set the new watched variables */
   consdata->watchedvar1 = watchedvar1;
   consdata->watchedvar2 = watchedvar2;
   
   return SCIP_OKAY;
}

/** adds coefficient in logicor constraint */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logicor constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   SCIP_CALL( consdataEnsureVarsSize(scip, consdata, consdata->nvars + 1) );
   consdata->vars[consdata->nvars] = var;
   SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[consdata->nvars]) );
   consdata->nvars++;

   /* we only catch this event in presolving stage */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_CONSHDLR* conshdlr;

      conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
      assert(conshdlr != NULL);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
            (SCIP_EVENTDATA*)cons, NULL) );
   }

   consdata->sorted = (consdata->nvars == 1);
   consdata->changed = TRUE;

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockRounding(scip, cons, var) );

   /* add the new coefficient to the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, 1.0) );
   }

   consdata->merged = FALSE;

   return SCIP_OKAY;
}

/** deletes coefficient at given position from logic or constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint */
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

   /* remove the rounding locks of variable */
   SCIP_CALL( unlockRounding(scip, cons, consdata->vars[pos]) );

   /* we only catch this event in presolving stage, so we need to only drop it there */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_VARFIXED, eventhdlr,
            (SCIP_EVENTDATA*)cons, -1) );
   }

   if( SCIPconsIsTransformed(cons) )
   {
      /* if the position is watched, stop watching the position */
      if( consdata->watchedvar1 == pos )
      {
         SCIP_CALL( switchWatchedvars(scip, cons, eventhdlr, consdata->watchedvar2, -1) );
      }
      if( consdata->watchedvar2 == pos )
      {
         SCIP_CALL( switchWatchedvars(scip, cons, eventhdlr, consdata->watchedvar1, -1) );
      }
   }
   assert(pos != consdata->watchedvar1);
   assert(pos != consdata->watchedvar2);

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &consdata->vars[pos]) );

   /* move the last variable to the free slot */
   if( pos != consdata->nvars - 1 )
   {
      consdata->vars[pos] = consdata->vars[consdata->nvars-1];
      consdata->sorted = FALSE;
   }
   consdata->nvars--;

   /* if the last variable (that moved) was watched, update the watched position */
   if( consdata->watchedvar1 == consdata->nvars )
      consdata->watchedvar1 = pos;
   if( consdata->watchedvar2 == consdata->nvars )
      consdata->watchedvar2 = pos;

   consdata->changed = TRUE;

   SCIP_CALL( SCIPenableConsPropagation(scip, cons) );

   return SCIP_OKAY;
}
/** in case a part (more than one variable) in the logic or constraint is independent of every else, we can perform dual
 *  reductions;
 *  - fix the variable with the smallest object coefficient to one if the constraint is not modifiable and all
 *    variable are independant
 *  - fix all independant variables with negative object coefficient to one
 *  - fix all remaining independant variables to zero 
 *
 * Note: the following dual reduction for logic or constraints is already performed by the presolver "dualfix"
 *       - if a variable in a set covering constraint is only locked by that constraint and has negative or zero
 *         objective coefficient than it can be fixed to one
 */
static
SCIP_RETCODE dualPresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< setppc constraint */
   int*                  nfixedvars,         /**< pointer to count number of fixings */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints  */
   SCIP_RESULT*          result              /**< pointer to store the result SCIP_SUCCESS, if presolving was performed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_VAR* activevar;
   SCIP_Real bestobjval;
   SCIP_Real objval;
   SCIP_Real fixval;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;
   SCIP_Bool negated;
   int nfixables;
   int nvars;
   int idx;
   int v;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);
   assert(result != NULL);

   /* constraints for which the check flag is set to FALSE, did not contribute to the lock numbers; therefore, we cannot
    * use the locks to decide for a dual reduction using this constraint; for example after a restart the cuts which are
    * added to the problems have the check flag set to FALSE 
    */
   if( !SCIPconsIsChecked(cons) )
      return SCIP_OKAY;
   
   assert(SCIPconsIsActive(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   
   /* we don't want to consider small constraints (note that the constraints can be modifiable, so we can't delete this
    * constraint) 
    */
   if( nvars < 2 )
      return SCIP_OKAY;

   vars = consdata->vars;
   idx = -1;
   bestobjval = SCIP_INVALID;
   
   nfixables = 0;

   /* check if we can apply the dual reduction; therefore count the number of variables where the logic or has the only
    * locks on 
    */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);

      /* the variable should not be (globally) fixed */
      assert(SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5);
      
      /* in case an other constraints has also locks on that variable we cannot perform a dual reduction on these
       * variables
       */
      if( SCIPvarGetNLocksDown(var) > 1 || SCIPvarGetNLocksUp(var) > 0 ) 
         continue;
      
      ++nfixables;
      negated = FALSE;
      
      /* get the active variable */
      SCIP_CALL( SCIPvarGetProbvarBinary(&var, &negated) );
      assert(SCIPvarIsActive(var));

      if( negated )
         objval = -SCIPvarGetObj(var);
      else
         objval = SCIPvarGetObj(var);
      
      /* check if the current variable has a smaller objective coefficient */
      if( SCIPisLT(scip, objval, bestobjval) )
      {
         idx = v;
         bestobjval = objval;
      }
   }

   if( nfixables < 2 )
      return SCIP_OKAY;

   assert(idx >= 0 && idx < nvars);
   assert(bestobjval < SCIPinfinity(scip));
   
   *result = SCIP_SUCCESS;

   /* fix all redundant variables to their best bound */
   
   /* first part of all variables */
   for( v = 0; v < idx; ++v ) 
   {
      var = vars[v];
      assert(var != NULL);
         
      /* in case an other constraints has also locks on that variable we cannot perform a dual reduction on these
       * variables
       */
      if( SCIPvarGetNLocksDown(var) > 1 || SCIPvarGetNLocksUp(var) > 0 ) 
         continue;

      activevar = var;
      negated = FALSE;

      /* get the active variable */
      SCIP_CALL( SCIPvarGetProbvarBinary(&activevar, &negated) );
      assert(SCIPvarIsActive(activevar));
         
      if( negated )
         objval = -SCIPvarGetObj(activevar);
      else
         objval = SCIPvarGetObj(activevar);

      if( objval > 0.0 )
         fixval = 0.0;
      else
         fixval = 1.0;
            
      SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) ); 
      assert(!infeasible); 
      assert(fixed);
   
      SCIPdebugMessage(" -> fixed <%s> == %g\n", SCIPvarGetName(var), fixval);
      ++(*nfixedvars);
   }

   /* second part of all variables */
   for( v = idx + 1; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);

      /* in case an other constraints has also locks on that variable we cannot perform a dual reduction on these
       * variables
       */
      if( SCIPvarGetNLocksDown(var) > 1 || SCIPvarGetNLocksUp(var) > 0 )
         continue;
      
      activevar = var;
      negated = FALSE;

      /* get the active variable */
      SCIP_CALL( SCIPvarGetProbvarBinary(&activevar, &negated) );
      assert(SCIPvarIsActive(activevar));
      
      if( negated )
         objval = -SCIPvarGetObj(activevar);
      else
         objval = SCIPvarGetObj(activevar);
      
      if( objval > 0.0 )
         fixval = 0.0;
      else
         fixval = 1.0;

      SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
      assert(!infeasible);
      assert(fixed);
      
      SCIPdebugMessage(" -> fixed <%s> == %g\n", SCIPvarGetName(var), fixval);
      ++(*nfixedvars);
   }
   
   /* if all variable have our appreciated number of locks and the constraint is not modifiable, or if the bestobjval is
    * less than or equal to zero, we can fix the variable with the smallest objective coefficient to one and the
    * constraint gets redundant
    */
   if( (nfixables == nvars && !SCIPconsIsModifiable(cons)) || bestobjval <= 0.0 )
   {
      SCIP_CALL( SCIPfixVar(scip, vars[idx], 1.0, &infeasible, &fixed) );
      assert(!infeasible);
      assert(fixed);
      
      SCIPdebugMessage(" -> fixed <%s> == 1.0\n", SCIPvarGetName(vars[idx]));
      ++(*nfixedvars);

      /* remove constraint since i*/
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
   }

   return SCIP_OKAY;
}

/** deletes all zero-fixed variables, checks for variables fixed to one, replace all variables which are not active or
 *  not a negation of an active variable by there active or negation of an active counterpart
 */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            redundant           /**< returns whether a variable fixed to one exists in the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int v;
   SCIP_VAR** vars;
   SCIP_Bool* negarray;
   int nvars;

   assert(eventhdlr != NULL);
   assert(redundant != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);

   *redundant = FALSE;
   v = 0;

   /* remove zeros and mark constraint redundant when found one variable fixed to one */
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      assert(SCIPvarIsBinary(var));

      if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         *redundant = TRUE;
         return SCIP_OKAY;
      }
      else if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
      }
      else
         ++v;
   }

   if( consdata->nvars == 0 )
      return SCIP_OKAY;

   nvars = consdata->nvars;

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &negarray, nvars) );

   /* get active or negation of active variables */
   SCIP_CALL( SCIPgetBinvarRepresentatives(scip, nvars, consdata->vars, vars, negarray) );

   /* renew all variables, important that we do a backwards loop because deletion only affect rear items */
   for( v = nvars - 1; v >= 0; --v )
   {
      if( vars[v] != consdata->vars[v] )
      {
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         SCIP_CALL( addCoef(scip, cons, vars[v]) );
      }
   }
   assert(consdata->nvars == nvars);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &negarray);
   SCIPfreeBufferArray(scip, &vars);

   SCIPdebugMessage("after fixings: ");
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL, TRUE)) );

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< logic or constraint that detected the conflict */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   /* conflict analysis can only be applied in solving stage and if it is applicable */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
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

/** disables or deletes the given constraint, depending on the current depth */
static
SCIP_RETCODE disableCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< bound disjunction constraint to be disabled */
   )
{
   if( SCIPgetDepth(scip) == 0 )
   {
      SCIP_CALL( SCIPdelCons(scip, cons) );
   }
   else
   {
      SCIPdebugMessage("disabling constraint cons <%s> at depth %d\n", SCIPconsGetName(cons), SCIPgetDepth(scip));
      SCIP_CALL( SCIPdisableCons(scip, cons) );
   }

   return SCIP_OKAY;
}

/** find pairs of negated variables in constraint: constraint is redundant */
/** find sets of equal variables in constraint: multiple entries of variable can be replaced by single entry */
static
SCIP_RETCODE mergeMultiples(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   unsigned char**       entries,            /**< array to store whether two positions in constraints represent the same variable */
   int*                  nentries,           /**< pointer for array size, if array will be to small it's corrected */
   SCIP_Bool*            redundant,          /**< returns whether a variable fixed to one exists in the constraint */
   int*                  nchgcoefs           /**< pointer to count number of changed/deleted coefficients */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   SCIP_Bool* negarray;
   SCIP_VAR* var;
   int v;
   int pos;
   int nintvars;
#ifndef NDEBUG
   int nbinvars;
   int nimplvars;
#endif

   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);
   assert(*entries != NULL);
   assert(nentries != NULL);
   assert(redundant != NULL);
   assert(nchgcoefs != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   *redundant = FALSE;

   if( consdata->merged )
      return SCIP_OKAY;

   if( consdata->nvars <= 1 )
   {
      consdata->merged = TRUE;
      return SCIP_OKAY;
   }

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &negarray, nvars) );

   assert(consdata->vars != NULL && nvars > 0);

   nintvars = SCIPgetNIntVars(scip);
#ifndef NDEBUG
   nbinvars = SCIPgetNBinVars(scip);
   nimplvars = SCIPgetNImplVars(scip);
#endif

   assert(*nentries >= nbinvars + nimplvars);

   /* get active or negation of active variables */
   SCIP_CALL( SCIPgetBinvarRepresentatives(scip, nvars, consdata->vars, vars, negarray) );

   /* initialize entries array */
   for( v = nvars - 1; v >= 0; --v )
   {
      assert(negarray[v] ? SCIPvarIsNegated(vars[v]) : TRUE);
      var = negarray[v] ? SCIPvarGetNegationVar(vars[v]) : vars[v];

      pos = SCIPvarGetProbindex(var);

      assert(SCIPvarIsActive(var));
      assert((pos < nbinvars && SCIPvarGetType(var) == SCIP_VARTYPE_BINARY)
	 || (pos >= (nbinvars + nintvars) && pos < (nbinvars + nintvars + nimplvars)
	    && SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT && SCIPvarIsBinary(var)));

      /* subtract number of integer variables because we only allocated memory for all binary and implicit variables */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
	 pos -= nintvars;

      /* var is not active yet */
      (*entries)[pos] = 0;
   }

   /* check all vars for multiple entries, do necessary backwards loop because deletion only affect rear items */
   for( v = nvars - 1; v >= 0; --v )
   {
      var = negarray[v] ? SCIPvarGetNegationVar(vars[v]) : vars[v];

      pos = SCIPvarGetProbindex(var);

      /* subtract number of integer variables because we only allocated memory for all binary and implicit variables */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
	 pos -= nintvars;

      /* if var occurs first time in constraint init entries array */
      if( (*entries)[pos] == 0 )
         (*entries)[pos] = negarray[v] ? 2 : 1;
      /* if var occurs second time in constraint, first time it was not negated */
      else if( (*entries)[pos] == 1 )
      {
         if( negarray[v] )
         {
            *redundant = TRUE;
            goto TERMINATE;
         }
         else
         {
            SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
	    ++(*nchgcoefs);
         }
      }
      /* if var occurs second time in constraint, first time it was negated */
      else
      {
         if( !negarray[v] )
         {
            *redundant = TRUE;
            goto TERMINATE;
         }
         else
         {
            SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
	    ++(*nchgcoefs);
         }
      }
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &negarray);
   SCIPfreeBufferArray(scip, &vars);

   consdata->merged = TRUE;

   return SCIP_OKAY;
}

/** checks constraint for violation only looking at the watched variables, applies fixings if possible */
static
SCIP_RETCODE processWatchedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint to be processed */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   SCIP_Bool*            addcut,             /**< pointer to store whether this constraint must be added as a cut */
   SCIP_Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Longint nbranchings1;
   SCIP_Longint nbranchings2;
   int nvars;
   int watchedvar1;
   int watchedvar2;

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

   SCIPdebugMessage("processing watched variables of constraint <%s>\n", SCIPconsGetName(cons));

   vars = consdata->vars;
   nvars = consdata->nvars;
   assert(nvars == 0 || vars != NULL);

   /* check watched variables if they are fixed to one */
   if( consdata->watchedvar1 >= 0 && SCIPvarGetLbLocal(vars[consdata->watchedvar1]) > 0.5 )
   {
      /* the variable is fixed to one, making the constraint redundant -> disable the constraint */
      SCIPdebugMessage(" -> disabling constraint <%s> (watchedvar1 fixed to 1.0)\n", SCIPconsGetName(cons));
      SCIP_CALL( disableCons(scip, cons) );
      return SCIP_OKAY;
   }
   if( consdata->watchedvar2 >= 0 && SCIPvarGetLbLocal(vars[consdata->watchedvar2]) > 0.5 )
   {
      /* the variable is fixed to one, making the constraint redundant -> disable the constraint */
      SCIPdebugMessage(" -> disabling constraint <%s> (watchedvar2 fixed to 1.0)\n", SCIPconsGetName(cons));
      SCIP_CALL( disableCons(scip, cons) );
      return SCIP_OKAY;
   }

   /* check if watched variables are still unfixed */
   watchedvar1 = -1;
   watchedvar2 = -1;
   nbranchings1 = SCIP_LONGINT_MAX;
   nbranchings2 = SCIP_LONGINT_MAX;
   if( consdata->watchedvar1 >= 0 && SCIPvarGetUbLocal(vars[consdata->watchedvar1]) > 0.5 )
   {
      watchedvar1 = consdata->watchedvar1;
      nbranchings1 = -1; /* prefer keeping the watched variable */
   }
   if( consdata->watchedvar2 >= 0 && SCIPvarGetUbLocal(vars[consdata->watchedvar2]) > 0.5 )
   {
      if( watchedvar1 == -1 )
      {
         watchedvar1 = consdata->watchedvar2;
         nbranchings1 = -1; /* prefer keeping the watched variable */
      }
      else
      {
         watchedvar2 = consdata->watchedvar2;
         nbranchings2 = -1; /* prefer keeping the watched variable */
      }
   }
   assert(watchedvar1 >= 0 || watchedvar2 == -1);
   assert(nbranchings1 <= nbranchings2);

   /* search for new watched variables */
   if( watchedvar2 == -1 )
   {
      int v;

      for( v = 0; v < nvars; ++v )
      {
         SCIP_Longint nbranchings;

         /* don't process the watched variables again */
         if( v == consdata->watchedvar1 || v == consdata->watchedvar2 )
            continue;

         /* check, if the variable is fixed */
         if( SCIPvarGetUbLocal(vars[v]) < 0.5 )
            continue;
         
         /* check, if the literal is satisfied */
         if( SCIPvarGetLbLocal(vars[v]) > 0.5 )
         {
            assert(v != consdata->watchedvar1);
            assert(v != consdata->watchedvar2);
            
            /* the variable is fixed to one, making the constraint redundant;
             * make sure, the feasible variable is watched and disable the constraint
             */
            SCIPdebugMessage(" -> disabling constraint <%s> (variable <%s> fixed to 1.0)\n", 
               SCIPconsGetName(cons), SCIPvarGetName(vars[v]));
            if( consdata->watchedvar1 != -1 )
            {
               SCIP_CALL( switchWatchedvars(scip, cons, eventhdlr, consdata->watchedvar1, v) );
            }
            else
            {
               SCIP_CALL( switchWatchedvars(scip, cons, eventhdlr, v, consdata->watchedvar2) );
            }
            SCIP_CALL( disableCons(scip, cons) );
            return SCIP_OKAY;
         }
         
         /* the variable is unfixed and can be used as watched variable */
         nbranchings = SCIPvarGetNBranchingsCurrentRun(vars[v], SCIP_BRANCHDIR_DOWNWARDS);
         assert(nbranchings >= 0);
         if( nbranchings < nbranchings2 )
         {
            if( nbranchings < nbranchings1 )
            {
               watchedvar2 = watchedvar1;
               nbranchings2 = nbranchings1;
               watchedvar1 = v;
               nbranchings1 = nbranchings;
            }
            else
            {
               watchedvar2 = v;
               nbranchings2 = nbranchings;
            }
         }
      }
   }
   assert(nbranchings1 <= nbranchings2);
   assert(watchedvar1 >= 0 || watchedvar2 == -1);

   if( watchedvar1 == -1 )
   {
      /* there is no unfixed variable left -> the constraint is infeasible
       *  - a modifiable constraint must be added as a cut and further pricing must be performed in the LP solving loop
       *  - an unmodifiable constraint is infeasible and the node can be cut off
       */
      assert(watchedvar2 == -1);

      SCIPdebugMessage(" -> constraint <%s> is infeasible\n", SCIPconsGetName(cons));

      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      if( SCIPconsIsModifiable(cons) )
         *addcut = TRUE;
      else
      {
         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflict(scip, cons) );

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
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(vars[watchedvar1]), 0.0));
      assert(SCIPisFeasEQ(scip, SCIPvarGetUbLocal(vars[watchedvar1]), 1.0));
      if( SCIPconsIsModifiable(cons) )
         *mustcheck = TRUE;
      else
      {
         SCIP_Bool infbdchg;

         /* fixed remaining variable to one and disable constraint; make sure, the fixed-to-one variable is watched */
         SCIPdebugMessage(" -> single-literal constraint <%s> (fix <%s> to 1.0) at depth %d\n", 
            SCIPconsGetName(cons), SCIPvarGetName(vars[watchedvar1]), SCIPgetDepth(scip));
         SCIP_CALL( SCIPinferBinvarCons(scip, vars[watchedvar1], TRUE, cons, 0, &infbdchg, NULL) );
         assert(!infbdchg);
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         if( watchedvar1 != consdata->watchedvar1 ) /* keep one of the watched variables */
         {
            SCIP_CALL( switchWatchedvars(scip, cons, eventhdlr, watchedvar1, consdata->watchedvar1) );
         }
         SCIP_CALL( disableCons(scip, cons) );
         *reduceddom = TRUE;
      }
   }
   else
   {
      SCIPdebugMessage(" -> new watched variables <%s> and <%s> of constraint <%s> are still unfixed\n",
         SCIPvarGetName(vars[watchedvar1]), SCIPvarGetName(vars[watchedvar2]), SCIPconsGetName(cons));

      /* switch to the new watched variables */
      SCIP_CALL( switchWatchedvars(scip, cons, eventhdlr, watchedvar1, watchedvar2) );

      /* there are at least two unfixed variables -> the constraint must be checked manually */
      *mustcheck = TRUE;

      /* disable propagation of constraint until a watched variable gets fixed */
      SCIP_CALL( SCIPdisableConsPropagation(scip, cons) );

      /* increase aging counter */
      SCIP_CALL( SCIPaddConsAge(scip, cons, AGEINCREASE(consdata->nvars)) );
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, returns TRUE iff constraint is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint to be checked */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            violated            /**< pointer to store whether the given solution violates the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real solval;
   SCIP_Real sum;
   int nvars;
   int v;

   assert(violated != NULL);

   *violated = FALSE;
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   nvars = consdata->nvars;
   
   /* calculate the constraint's activity */
   sum = 0.0;
   solval = 0.0;
   for( v = 0; v < nvars && sum < 1.0; ++v )
   {
      assert(SCIPvarIsBinary(vars[v]));
      solval = SCIPgetSolVal(scip, sol, vars[v]);
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
      sum += solval;
   }

   *violated = SCIPisFeasLT(scip, sum, 1.0);

   return SCIP_OKAY;
}

/** creates an LP row in a logic or constraint data object */
static
SCIP_RETCODE createRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< logic or constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->row, SCIPconsGetHdlr(cons), SCIPconsGetName(cons), 1.0, SCIPinfinity(scip),
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->row, consdata->nvars, consdata->vars, 1.0) );

   return SCIP_OKAY;
}

/** adds logic or constraint as cut to the LP */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint */
   SCIP_SOL*             sol                 /**< primal CIP solution, NULL for current LP solution */
   )
{
   SCIP_CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   if( consdata->row == NULL )
   {
      /* convert logic or constraint data into LP row */
      SCIP_CALL( createRow(scip, cons) );
   }
   assert(consdata->row != NULL);

   /* insert LP row as cut */
   if( !SCIProwIsInLP(consdata->row) )
   {
      SCIPdebugMessage("adding constraint <%s> as cut to the LP\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPaddCut(scip, sol, consdata->row, FALSE) );
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint to be separated */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            separated,          /**< pointer to store TRUE, if a cut was found */
   SCIP_Bool*            reduceddom          /**< pointer to store TRUE, if a domain reduction was found */
   )
{
   SCIP_Bool addcut;
   SCIP_Bool mustcheck;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(separated != NULL);
   assert(reduceddom != NULL);

   SCIPdebugMessage("separating constraint <%s>\n", SCIPconsGetName(cons));

   /* update and check the watched variables, if they were changed since last processing */
   if( sol == NULL && SCIPconsIsPropagationEnabled(cons) )
   {
      SCIP_CALL( processWatchedVars(scip, cons, eventhdlr, cutoff, reduceddom, &addcut, &mustcheck) );
   }
   else
   {
      addcut = FALSE;
      mustcheck = TRUE;
   }

   if( mustcheck )
   {
      SCIP_CONSDATA* consdata;

      assert(!addcut);
      
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* variable's fixings didn't give us any information -> we have to check the constraint */
      if( sol == NULL && consdata->row != NULL )
      {
         /* skip constraints already in the LP */
         if( SCIProwIsInLP(consdata->row) )
            return SCIP_OKAY;
         else
         {
            SCIP_Real feasibility;
            
            assert(!SCIProwIsInLP(consdata->row));
            feasibility = SCIPgetRowLPFeasibility(scip, consdata->row);
            addcut = SCIPisFeasNegative(scip, feasibility);
         }
      }
      else
      {
         SCIP_CALL( checkCons(scip, cons, sol, &addcut) );
      }
   }

   if( addcut )
   {
      /* insert LP row as cut */
      SCIP_CALL( addCut(scip, cons, sol) );
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *separated = TRUE;
   }

   return SCIP_OKAY;
}

/** enforces the pseudo solution on the given constraint */
static
SCIP_RETCODE enforcePseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint to be separated */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the constraint was infeasible */
   SCIP_Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   SCIP_Bool*            solvelp             /**< pointer to store TRUE, if the LP has to be solved */
   )
{
   SCIP_Bool addcut;
   SCIP_Bool mustcheck;

   assert(!SCIPhasCurrentNodeLP(scip));
   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(reduceddom != NULL);
   assert(solvelp != NULL);

   /* update and check the watched variables, if they were changed since last processing */
   if( SCIPconsIsPropagationEnabled(cons) )
   {
      SCIP_CALL( processWatchedVars(scip, cons, eventhdlr, cutoff, reduceddom, &addcut, &mustcheck) );
   }
   else
   {
      addcut = FALSE;
      mustcheck = TRUE;
   }

   if( mustcheck )
   {
      SCIP_Bool violated;

      assert(!addcut);

      SCIP_CALL( checkCons(scip, cons, NULL, &violated) );
      if( violated )
      {
         /* constraint was infeasible -> reset age */
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *infeasible = TRUE;
      }
   }
   else if( addcut )
   {
      /* a cut must be added to the LP -> we have to solve the LP immediately */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *solvelp = TRUE;
   }

   return SCIP_OKAY;
}

/** sorts logicor constraint's variables by non-decreasing variable index */
static
void consdataSort(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->sorted )
   {
      if( consdata->nvars <= 1 )
	 consdata->sorted = TRUE;
      else
      {
	 SCIP_VAR* var1 = NULL;
	 SCIP_VAR* var2 = NULL;

	 /* remember watch variables */
	 if( consdata->watchedvar1 != -1 )
	 {
	    var1 = consdata->vars[consdata->watchedvar1];
	    assert(var1 != NULL);
	    consdata->watchedvar1 = -1;
	    if( consdata->watchedvar2 != -1 )
	    {
	       var2 = consdata->vars[consdata->watchedvar2];
	       assert(var2 != NULL);
	       consdata->watchedvar2 = -1;
	    }
	 }
	 assert(consdata->watchedvar1 == -1);
	 assert(consdata->watchedvar2 == -1);
	 assert(var1 != NULL || var2 == NULL);

	 /* sort variables after index */
	 SCIPsortPtr((void**)consdata->vars, SCIPvarComp, consdata->nvars);
	 consdata->sorted = TRUE;

	 /* correct watched variables */
	 if( var1 != NULL )
	 {
	    int pos;
#ifndef NDEBUG
	    SCIP_Bool found;

	    found = SCIPsortedvecFindPtr((void**)consdata->vars, SCIPvarComp, (void*)var1, consdata->nvars, &pos);
	    assert(found);
#else
	    SCIPsortedvecFindPtr((void**)consdata->vars, SCIPvarComp, (void*)var1, consdata->nvars, &pos);
#endif
	    assert(pos >= 0 && pos < consdata->nvars);
	    consdata->watchedvar1 = pos;

	    if( var2 != NULL )
	    {
#ifndef NDEBUG
	       found = SCIPsortedvecFindPtr((void**)consdata->vars, SCIPvarComp, (void*)var2, consdata->nvars, &pos);
	       assert(found);
#else
	       SCIPsortedvecFindPtr((void**)consdata->vars, SCIPvarComp, (void*)var2, consdata->nvars, &pos);
#endif
	       assert(pos >= 0 && pos < consdata->nvars);
	       consdata->watchedvar2 = pos;
	    }
	 }
      }
   }

#ifdef SCIP_DEBUG
   /* check sorting */
   {
      int v;

      for( v = 0; v < consdata->nvars; ++v )
      {
         assert(v == consdata->nvars-1 || SCIPvarCompare(consdata->vars[v], consdata->vars[v+1]) <= 0);
      }
   }
#endif
}

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyLogicorcons)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqLogicorcons)
{
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
   SCIP_Bool coefsequal;
   int i;
#ifndef NDEBUG
   SCIP* scip;

   scip = (SCIP*)userptr;
   assert(scip != NULL);
#endif

   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);

   /* checks trivial case */
   if( consdata1->nvars != consdata2->nvars )
      return FALSE;

   /* sorts the constraints */
   consdataSort(consdata1);
   consdataSort(consdata2);
   assert(consdata1->sorted);
   assert(consdata2->sorted);

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
SCIP_DECL_HASHKEYVAL(hashKeyValLogicorcons)
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
   hashtablesize = MAX(hashtablesize, HASHSIZE_LOGICORCONS);
   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashGetKeyLogicorcons, hashKeyEqLogicorcons, hashKeyValLogicorcons, (void*) scip) );

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
      consdataSort(consdata0);
      assert(consdata0->sorted);

      /* get constraint from current hash table with same variables as cons0 */
      cons1 = (SCIP_CONS*)(SCIPhashtableRetrieve(hashtable, (void*)cons0));
 
      if( cons1 != NULL )
      {
#ifndef NDEBUG
         SCIP_CONSDATA* consdata1;
#endif

         assert(SCIPconsIsActive(cons1));
         assert(!SCIPconsIsModifiable(cons1));
      
#ifndef NDEBUG
         consdata1 = SCIPconsGetData(cons1);
#endif
         assert(consdata0 != NULL && consdata1 != NULL);
         assert(consdata0->nvars >= 1 && consdata0->nvars == consdata1->nvars);
         
         assert(consdata0->sorted && consdata1->sorted);
         assert(consdata0->vars[0] == consdata1->vars[0]);

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( updateFlags(scip, cons1, cons0) ); 

         /* delete consdel */
         SCIP_CALL( SCIPdelCons(scip, cons0) );
         (*ndelconss)++;

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

   /* free hash table */
   SCIPhashtableFree(&hashtable);

   return SCIP_OKAY;
}

/** removes the redundant second constraint and updates the flags of the first one */
static
SCIP_RETCODE removeRedundantCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< constraint that should stay */
   SCIP_CONS*            cons1,              /**< constraint that should be deleted */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   assert(ndelconss != NULL);

   SCIPdebugMessage(" -> removing logicor constraint <%s> which is redundant to <%s>\n",
      SCIPconsGetName(cons1), SCIPconsGetName(cons0));
   SCIPdebugPrintCons(scip, cons0, NULL);
   SCIPdebugPrintCons(scip, cons1, NULL);

   /* update flags of cons0 */
   SCIP_CALL( updateFlags(scip, cons0, cons1) );

   /* delete cons1 */
   SCIP_CALL( SCIPdelCons(scip, cons1) );
   (*ndelconss)++;

   return SCIP_OKAY;
}


/** deletes redundant constraints */
static
SCIP_RETCODE removeRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int*                  firstchange,        /**< first constraint that changed since last pair preprocessing round */
   int                   chkind,             /**< index of constraint to check against all prior indices upto startind */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   SCIP_Bool cons0changed;
   int c;

   assert(ndelconss != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);

   /* sort the constraint */
   consdataSort(consdata0);
   assert(consdata0->sorted);

   /* check constraint against all prior constraints */
   cons0changed = consdata0->changed;
   consdata0->changed = FALSE;
   for( c = (cons0changed ? 0 : *firstchange); c < chkind && SCIPconsIsActive(cons0); ++c )
   {
      SCIP_CONS* consstay;
      SCIP_CONS* consdel;
      SCIP_CONSDATA* consdatastay;
      SCIP_CONSDATA* consdatadel;
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata1;
      SCIP_Bool consdelisredundant;
      int v0;
      int v1;

      cons1 = conss[c];
      assert(SCIPconsIsActive(cons0));

      /* ignore inactive and modifiable constraints */
      if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
         continue;

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

      /* sort the constraint */
      consdataSort(consdata1);
      assert(consdata1->sorted);

      if( consdata0->nvars <= consdata1->nvars )
      {
         consstay = cons0;
         consdel = cons1;
         consdatastay = consdata0;
         consdatadel = consdata1;
      }
      else
      {
         consstay = cons1;
         consdel = cons0;
         consdatastay = consdata1;
         consdatadel = consdata0;
      }

      v0 = 0;
      v1 = 0;
      consdelisredundant = TRUE;

      while( v0 < consdatastay->nvars && v1 < consdatadel->nvars )
      {
         int index0;
         int index1;

         index0 = SCIPvarGetIndex(consdatastay->vars[v0]);
         index1 = SCIPvarGetIndex(consdatadel->vars[v1]);
         if( index1 < index0 )
         {
            for( ++v1; v1 < consdatadel->nvars; ++v1 )
            {
               index1 = SCIPvarGetIndex(consdatadel->vars[v1]);
               if( index1 >= index0 )
                  break;
            }
         }
         if( index0 == index1 )
         {
            v0++;
            v1++;
         }
         else
         {
            consdelisredundant = FALSE;
            break;
         }
      }

      if (v0 < consdatastay->nvars)
            consdelisredundant = FALSE;

      if( consdelisredundant )
      {
         /* delete consdel */
         SCIPdebugMessage("logicor constraint <%s> is contained in <%s>\n", SCIPconsGetName(consdel), SCIPconsGetName(consstay));
         SCIPdebugPrintCons(scip, consstay, NULL);
         SCIPdebugPrintCons(scip, consdel, NULL);
         SCIP_CALL( removeRedundantCons(scip, consstay, consdel, ndelconss) );

         /* update the first changed constraint to begin the next aggregation round with */
         if( consdatastay->changed && SCIPconsGetPos(consstay) < *firstchange )
            *firstchange = SCIPconsGetPos(consstay);
      }
   }

   return SCIP_OKAY;
}

#define MAXCOMPARISONS 1000000

/* try to find a negated clique in a constraint which makes this constraint redundant but we need to keep the negated
 * clique information alive, so we create a corresponding set-packing constraint
 */
static
SCIP_RETCODE removeConstraintsDueToNegCliques(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< logicor constraint handler */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_CONS**           conss,              /**< all constraints */
   int                   nconss,             /**< number of constraints */
   int*                  nupgdconss,         /**< pointer to count number of upgraded constraints */
   int*                  nchgcoefs           /**< pointer to count number of changed/deleted coefficients */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** repvars;
   SCIP_Bool* negated;
   SCIP_VAR* var1;
   int c;
   int size;
   int maxcomppercons;
   int comppercons;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(eventhdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nupgdconss != NULL);
   assert(nchgcoefs != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( nconss == 0 )
      return SCIP_OKAY;

   if( SCIPgetNCliques(scip) == conshdlrdata->nlastcliques && SCIPgetNImplications(scip) == conshdlrdata->nlastimpls )
      return SCIP_OKAY;

   /* estimate the maximal number ob variables in a logicor constraint */
   size = SCIPgetNBinVars(scip) + SCIPgetNImplVars(scip);
   assert(size > 0);

   /* temporary memory for active/negation of active variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &repvars, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &negated, size) );

   /* iterate over all constraints and try to find negated cliques in logicors */
   for( c = nconss - 1; c >= 0 && !SCIPisStopped(scip); --c )
   {
      int v;

      assert(conss != NULL); /* for flexelint */

      cons = conss[c];
      assert(cons != NULL);

      if( !SCIPconsIsActive(cons) )
         continue;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      assert(consdata->nvars > 1);
      assert(consdata->nvars <= size);

      if( SCIPconsIsModifiable(cons) && consdata->nvars == 2 )
         continue;

      maxcomppercons = MAXCOMPARISONS / nconss;
      comppercons = 0;

      /* get binary representations of constraint variables */
      SCIP_CALL( SCIPgetBinvarRepresentatives(scip, consdata->nvars, consdata->vars, repvars, negated) );
      for( v = consdata->nvars - 1; v > 0; --v )
      {
         SCIP_Bool breakloop;
         SCIP_Bool neg1;
         int w;

	 var1 = repvars[v];
	 neg1 = negated[v];

         /* if there is no negated variable, there can't be a negated clique */
         if( SCIPvarGetNegatedVar(var1) == NULL )
            continue;

         /* get active counterpart to check for common cliques */
         if( SCIPvarGetStatus(var1) == SCIP_VARSTATUS_NEGATED )
         {
            var1 = SCIPvarGetNegatedVar(var1);
            neg1 = TRUE;
         }
         else
            neg1 = FALSE;

         if( !SCIPvarIsActive(var1) )
            continue;

	 /* no cliques available */
	 if( SCIPvarGetNCliques(var1, neg1) == 0 && SCIPvarGetNImpls(var1, neg1) == 0 )
	    continue;

	 comppercons += (v - 1);

         breakloop = FALSE;

         for( w = v - 1; w >= 0; --w )
         {
	    SCIP_VAR* var2;
            SCIP_Bool neg2;

	    var2 = repvars[w];
	    neg2 = negated[w];

            /* if there is no negated variable, there can't be a negated clique */
            if( SCIPvarGetNegatedVar(var2) == NULL )
               continue;

            if( SCIPvarGetStatus(var2) == SCIP_VARSTATUS_NEGATED )
            {
               var2 = SCIPvarGetNegatedVar(var2);
               neg2 = TRUE;
            }
            else
               neg2 = FALSE;

            if( !SCIPvarIsActive(var2) )
               continue;

	    /* no cliques available */
	    if( SCIPvarGetNCliques(var2, neg2) == 0 && SCIPvarGetNImpls(var2, neg2) == 0 )
	       continue;

	    /* check if both active variable are the same */
	    if( var1 == var2 )
	    {
	       if( neg1 != neg2 )
	       {
		  SCIPdebugMessage("logicor constraint <%s> is redundant, because variable <%s> and its negation <%s> exist\n", SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));

		  SCIP_CALL( SCIPdelCons(scip, cons) );

		  breakloop = TRUE;
	       }
	       else
	       {
#ifndef NDEBUG
		  SCIP_VAR* lastvar = consdata->vars[consdata->nvars - 1];
#endif
		  SCIPdebugMessage("in logicor constraint <%s>, active variable of <%s> and active variable of <%s> are the same, removing the first\n", SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[v]), SCIPvarGetName(consdata->vars[w]));

		  SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );

		  if( v < consdata->nvars )
		  {
		     /* delCoefPos replaces the variable on position v with the last one, so w also need to correct the
		      * negated array the same way, and because of deletion the number of variables is already decreased
		      */
		     assert(consdata->vars[v] == lastvar);
		     negated[v] = negated[consdata->nvars];
		  }
		  ++(*nchgcoefs);
	       }
	       break;
	    }

            if( SCIPvarsHaveCommonClique(var1, neg1, var2, neg2, TRUE) )
            {
               SCIP_CONS* newcons;
               SCIP_VAR* vars[2];

               /* this negated clique information could be created out of this logicor constraint even if there are more
                * than two variables left (, for example by probing), we need to keep this information by creating a
                * setppc constraint instead 
                */

               /* get correct variables */
               if( !neg1 )
                  vars[0] = SCIPvarGetNegatedVar(var1);
               else
                  vars[0] = var1;

               if( !neg2 )
                  vars[1] = SCIPvarGetNegatedVar(var2);
               else
                  vars[1] = var2;

               SCIP_CALL( SCIPcreateConsSetpack(scip, &newcons, SCIPconsGetName(cons), 2, vars, 
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), 
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

               SCIP_CALL( SCIPaddCons(scip, newcons) );
               SCIPdebugPrintCons(scip, newcons, NULL);

               SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

               SCIPdebugMessage("logicor constraint <%s> is redundant due to negated clique information and will be replaced by a setppc constraint \n", SCIPconsGetName(cons));
               SCIPdebugMessage("variable <%s> and variable <%s> are in a negated clique\n", SCIPvarGetName(consdata->vars[v]), SCIPvarGetName(consdata->vars[w]));

               SCIP_CALL( SCIPdelCons(scip, cons) );
               ++(*nupgdconss);

               breakloop = TRUE;
               break;
            }
         }
         if( breakloop )
            break;

	 /* do not do to many comparisons */
	 if( comppercons > maxcomppercons )
	    break;
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &negated);
   SCIPfreeBufferArray(scip, &repvars);

   return SCIP_OKAY;
}

/*
 * upgrading of linear constraints
 */

/** creates and captures a normalized (with all coefficients +1) logic or constraint */
static
SCIP_RETCODE createNormalizedLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients (+1.0 or -1.0) */
   int                   mult,               /**< multiplier on the coefficients(+1 or -1) */
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
   SCIP_VAR** transvars;
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);
   assert(mult == +1 || mult == -1);

   /* get temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &transvars, nvars) );

   /* negate positive or negative variables */
   for( v = 0; v < nvars; ++v )
   {
      if( mult * vals[v] > 0.0 )
         transvars[v] = vars[v];
      else
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, vars[v], &transvars[v]) );
      }
      assert(transvars[v] != NULL);
   }

   /* create the constraint */
   SCIP_CALL( SCIPcreateConsLogicor(scip, cons, name, nvars, transvars,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &transvars);

   return SCIP_OKAY;
}

static
SCIP_DECL_LINCONSUPGD(linconsUpgdLogicor)
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
   if( nvars > 2 && nposbin + nnegbin == nvars && ncoeffspone + ncoeffsnone == nvars
      && ((SCIPisEQ(scip, lhs, 1.0 - ncoeffsnone) && SCIPisInfinity(scip, rhs))
         || (SCIPisInfinity(scip, -lhs) && SCIPisEQ(scip, rhs, ncoeffspone - 1.0))) )
   {
      int mult;

      SCIPdebugMessage("upgrading constraint <%s> to logic or constraint\n", SCIPconsGetName(cons));
      
      /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
      mult = SCIPisInfinity(scip, rhs) ? +1 : -1;
      
      /* create the logic or constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( createNormalizedLogicor(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, mult,
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
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyLogicor)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrLogicor(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free constraint handler data */
   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;
   int v;

   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* catch all variable event for deleted variables, which is only used in presolving */
   for( c = nconss - 1; c >= 0; --c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      for( v = consdata->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[v], SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
               (SCIP_EVENTDATA*)conss[c], NULL) );
      }
   }

   return SCIP_OKAY;
}
/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;
   int v;

   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* drop all variable event for deleted variables, which was only used in presolving */
   for( c = nconss - 1; c >= 0; --c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      for( v = consdata->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[v], SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
               (SCIP_EVENTDATA*)conss[c], -1) );
      }
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolLogicor)
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
SCIP_DECL_CONSDELETE(consDeleteLogicor)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int v;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      for( v = (*consdata)->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->vars[v], SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
               (SCIP_EVENTDATA*)cons, -1) );
      }
   }

   /* free LP row and logic or constraint */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransLogicor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   /*debugMessage("Trans method of logic or constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */

   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->nvars, sourcedata->vars) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int v;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      for( v = targetdata->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, targetdata->vars[v], SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
               (SCIP_EVENTDATA*)(*targetcons), NULL) );
      }
   }

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpLogicor)
{  /*lint --e{715}*/
   int c;

   for( c = 0; c < nconss; ++c )
   {
      assert(SCIPconsIsInitial(conss[c]));
      SCIP_CALL( addCut(scip, conss[c], NULL) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   SCIP_Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("separating %d/%d logic or constraints\n", nusefulconss, nconss);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* check all useful logic or constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
   }

   /* combine logic or constraints to get more cuts */
   /**@todo further cuts of logic or constraints */

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   SCIP_Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("separating %d/%d logic or constraints\n", nusefulconss, nconss);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* check all useful logic or constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
   }

   /* combine logic or constraints to get more cuts */
   /**@todo further cuts of logic or constraints */

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool separated;
   SCIP_Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("LP enforcing %d logic or constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* check all useful logic or constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
   }

   /* check all obsolete logic or constraints for feasibility */
   for( c = nusefulconss; c < nconss && !cutoff && !separated && !reduceddom; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, conshdlrdata->eventhdlr, &cutoff, &separated, &reduceddom) );
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


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_Bool reduceddom;
   SCIP_Bool solvelp;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("pseudo enforcing %d logic or constraints\n", nconss);

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
      SCIP_CALL( enforcePseudo(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &infeasible, &reduceddom, &solvelp) );
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else if( solvelp )
      *result = SCIP_SOLVELP;
   else if( infeasible )
      *result = SCIP_INFEASIBLE;
   
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckLogicor)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
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
         SCIP_Bool violated;

         SCIP_CALL( checkCons(scip, cons, sol, &violated) );
         if( violated )
         {
            /* constraint is violated */
            *result = SCIP_INFEASIBLE;

            if( printreason )
            {
#ifndef NDEBUG
               int v;
               for( v = 0; v < consdata->nvars; ++v )
               {
                  assert( consdata->vars[v] != NULL);
                  assert( SCIPvarIsBinary(consdata->vars[v]) );
                  assert( SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[v])) );
               }
#endif
               SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
               SCIPinfoMessage(scip, NULL, ";\nviolation: all variables are set to zero\n");
            }

            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool reduceddom;
   SCIP_Bool addcut;
   SCIP_Bool mustcheck;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   reduceddom = FALSE;

   /* propagate all useful logic or constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIPdebugMessage(" propagate constraint %s\n", SCIPconsGetName(conss[c]));
      SCIP_CALL( processWatchedVars(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &reduceddom, &addcut, &mustcheck) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   unsigned char* entries;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   SCIP_Bool fixed;
   int c;
   int firstchange;
   int nentries;
   int oldnfixedvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnupgdconss;
   int oldnchgcoefs;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   nentries = SCIPgetNBinVars(scip) + SCIPgetNImplVars(scip);

   oldnfixedvars = *nfixedvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;
   oldnchgcoefs = *nchgcoefs;

   firstchange = INT_MAX;

   SCIP_CALL( SCIPallocBufferArray(scip, &entries, nentries) );

   /* process constraints */
   for( c = 0; c < nconss && *result != SCIP_CUTOFF && !SCIPisStopped(scip); ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      SCIPdebugMessage("presolving logic or constraint <%s>\n", SCIPconsGetName(cons));

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
      {
         SCIP_CALL( SCIPenableConsPropagation(scip, cons) );
      }

      /* remove all variables that are fixed to zero, check redundancy due to fixed-to-one variable */
      SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr, &redundant) );

      /* find pairs of negated variables in constraint: constraint is redundant */
      /* find sets of equal variables in constraint: multiple entries of variable can be replaced by single entry */
      if( !redundant )
      {
         SCIP_CALL( mergeMultiples(scip, cons, conshdlrdata->eventhdlr, &entries, &nentries, &redundant, nchgcoefs) );
      }

      if( redundant )
      {
         SCIPdebugMessage("logic or constraint <%s> is redundant\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );
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
            SCIPdebugMessage("logic or constraint <%s> is infeasible\n", SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            goto TERMINATE;
         }
         else if( consdata->nvars == 1 )
         {
            SCIPdebugMessage("logic or constraint <%s> has only one variable not fixed to 0.0\n",
               SCIPconsGetName(cons));

            assert(consdata->vars != NULL);
            assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));

            if( SCIPvarGetStatus(consdata->vars[0]) != SCIP_VARSTATUS_MULTAGGR )
            {
               SCIPdebugMessage(" -> fix variable and delete constraint\n");

               SCIP_CALL( SCIPfixVar(scip, consdata->vars[0], 1.0, &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMessage(" -> infeasible fixing\n");
                  *result = SCIP_CUTOFF;
                  goto TERMINATE;
               }
               assert(fixed);
               (*nfixedvars)++;

               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
            else if( conshdlrdata->conshdlrlinear != NULL )
            {
               SCIP_Real coef;
               SCIP_CONS* conslinear;
               char consname[SCIP_MAXSTRLEN];

               SCIPdebugMessage(" -> variable is multi-aggregated, upgrade to linear constraint <%s> == 1 \n",
                  SCIPvarGetName(consdata->vars[0]));

               coef = 1.0;
               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "fixmaggr_%s_%s", SCIPconsGetName(cons),SCIPvarGetName(consdata->vars[0]) );
               SCIP_CALL( SCIPcreateConsLinear(scip, &conslinear, consname, 1, consdata->vars, &coef, 1.0, 1.0,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), 
                     SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                     SCIPconsIsStickingAtNode(cons)) );

               /* add constraint */
               SCIP_CALL( SCIPaddCons(scip, conslinear) );
               SCIP_CALL( SCIPreleaseCons(scip, &conslinear) );
               SCIP_CALL( SCIPdelCons(scip, cons) );       

               (*nupgdconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
         }
         else if( consdata->nvars == 2 && !consdata->impladded )
         {
            SCIP_Bool implinfeasible;
            int nimplbdchgs;
            SCIP_Bool values[2];

            values[0] = FALSE;
            values[1] = FALSE;
            /* a two-variable logicor constraint x + y >= 1 yields the implication x == 0 -> y == 1, and is represented
             * by the clique inequality ~x + ~y <= 1 
             */
            SCIP_CALL( SCIPaddClique(scip, consdata->vars, values, consdata->nvars, &implinfeasible, &nimplbdchgs) );
            *nchgbds += nimplbdchgs;
            if( implinfeasible )
            {
               *result = SCIP_CUTOFF;
               goto TERMINATE;
            }

            /* adding the above implication could lead to fixings, which render the constraint redundant */
            if ( nimplbdchgs > 0 )
            {
               /* remove all variables that are fixed to zero, check redundancy due to fixed-to-one variable */
               SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr, &redundant) );

               if( redundant )
               {
                  SCIPdebugMessage("logic or constraint <%s> is redundant\n", SCIPconsGetName(cons));
                  SCIP_CALL( SCIPdelCons(scip, cons) );
                  (*ndelconss)++;
                  *result = SCIP_SUCCESS;
                  continue;
               }
            }
            consdata->impladded = TRUE;
         }
      }

      /* perform dual reductions */
      if( conshdlrdata->dualpresolving )
      {
         SCIP_CALL( dualPresolving(scip, cons, nfixedvars, ndelconss, result) );

         /* if dual reduction deleted the constraint we take the next */
         if( !SCIPconsIsActive(cons) )
            continue;
      }

      /* remember the first changed constraint to begin the next redundancy round with */
      if( firstchange == INT_MAX && consdata->changed )
         firstchange = c;
   }

   assert(*result != SCIP_CUTOFF);

   /* fast preprocessing of pairs of logic or constraints, used for equal constraints */
   if( firstchange < nconss && conshdlrdata->presolusehashing )
   {
      /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
      SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &firstchange, ndelconss) );
   }

   /* preprocess pairs of logic or constraints and apply negated clique presolving */
   if( oldnfixedvars == *nfixedvars && oldnchgbds == *nchgbds && oldndelconss == *ndelconss && oldnupgdconss == *nupgdconss && oldnchgcoefs == *nchgcoefs )
   {
      /* check constraints for redundancy */
      if( conshdlrdata->presolpairwise )
      {
        SCIP_Longint npaircomparisons;
        npaircomparisons = 0;
        oldndelconss = *ndelconss;

        for( c = firstchange; c < nconss && !SCIPisStopped(scip); ++c )
         {
            if( SCIPconsIsActive(conss[c]) && !SCIPconsIsModifiable(conss[c]) )
            {
               npaircomparisons += (SCIPconsGetData(conss[c])->changed) ? (SCIP_Longint) c : ((SCIP_Longint) c - (SCIP_Longint) firstchange);

               SCIP_CALL( removeRedundantConstraints(scip, conss, &firstchange, c, ndelconss) );

               if( npaircomparisons > NMINCOMPARISONS )
               {
                  if( (*ndelconss - oldndelconss) / ((SCIP_Real)npaircomparisons) < MINGAINPERNMINCOMPARISONS )
                     break;
                  oldndelconss = *ndelconss;
                  npaircomparisons = 0;
               }
            }
         }
      }

      /* check for redundant constraints due to negated clique information */
      if( conshdlrdata->usenegatedclique )
      {
         SCIP_CALL( removeConstraintsDueToNegCliques(scip, conshdlr, conshdlrdata->eventhdlr, conss, nconss, nupgdconss, nchgcoefs) );

         conshdlrdata->nlastcliques = SCIPgetNCliques(scip);
         conshdlrdata->nlastimpls = SCIPgetNImplications(scip);
      }
   }

 TERMINATE:

   SCIPfreeBufferArray(scip, &entries);

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropLogicor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
#ifndef NDEBUG
   SCIP_Bool infervarfound;
#endif
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(infervar != NULL);
   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("conflict resolving method of logic or constraint handler\n");

   /* the only deductions are variables infered to 1.0 on logic or constraints where all other variables
    * are assigned to zero
    */
   assert(SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) > 0.5); /* the inference variable must be assigned to one */

#ifndef NDEBUG
   infervarfound = FALSE;
#endif
   for( v = 0; v < consdata->nvars; ++v )
   {
      if( consdata->vars[v] != infervar )
      {
         /* the reason variable must have been assigned to zero */
         assert(SCIPvarGetUbAtIndex(consdata->vars[v], bdchgidx, FALSE) < 0.5);
         SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
      }
#ifndef NDEBUG
      else
      {
         assert(!infervarfound);
         infervarfound = TRUE;
      }
#endif
   }
   assert(infervarfound);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockLogicor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* lock every single coefficient */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos, nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->watchedvar1 == -1 || consdata->watchedvar1 != consdata->watchedvar2);

   SCIPdebugMessage("activating information for logic or constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL, TRUE)) );

   /* catch events on watched variables */
   if( consdata->watchedvar1 != -1 )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[consdata->watchedvar1],
            SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)cons,
            &consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[consdata->watchedvar2],
            SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)cons,
            &consdata->filterpos2) );
   }

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveLogicor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->watchedvar1 == -1 || consdata->watchedvar1 != consdata->watchedvar2);

   SCIPdebugMessage("deactivating information for logic or constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL, TRUE)) );

   /* drop events on watched variables */
   if( consdata->watchedvar1 != -1 )
   {
      assert(consdata->filterpos1 != -1);
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar1],
            SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)cons,
            consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 )
   {
      assert(consdata->filterpos2 != -1);
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar2],
            SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)cons,
            consdata->filterpos2) );
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintLogicor)
{  /*lint --e{715}*/

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file, FALSE) );
    
   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyLogicor)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   const char* consname;
   int nvars;

   /* get variables and coefficients of the source constraint */
   sourcevars = SCIPgetVarsLogicor(sourcescip, sourcecons);
   nvars = SCIPgetNVarsLogicor(sourcescip, sourcecons);
   
   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   /* copy the logic using the linear constraint copy method */
   SCIP_CALL( SCIPcopyConsLinear(scip, cons, sourcescip, consname, nvars, sourcevars, NULL,
         1.0, SCIPinfinity(scip), varmap, consmap,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   assert(cons != NULL);
   
   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseLogicor)
{  /*lint --e{715}*/
   SCIP_VAR** vars;

   char* strcopy;
   char* token;
   char* saveptr;
   char* endptr;
   int requiredsize;
   int varssize;
   int nvars;
   
   SCIPdebugMessage("parse <%s> as logicor constraint\n", str);

   /* copy string for truncating it */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &strcopy, str, (int)(strlen(str)+1)));

   /* cutoff "logicor" form the constraint string */
   (void) SCIPstrtok(strcopy, "(", &saveptr ); 

   /* cutoff ")" form the constraint string */
   token = SCIPstrtok(NULL, ")", &saveptr ); 
   
   varssize = 100;
   nvars = 0;

   /* allocate buffer array for variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );

   /* parse string */
   SCIP_CALL( SCIPparseVarsList(scip, token, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );
   
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

      /* create logicor constraint */
      SCIP_CALL( SCIPcreateConsLogicor(scip, cons, name, nvars, vars,  
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   /* free buffers */
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &strcopy);
   
   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsLogicor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsLogicor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecLogicor)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   SCIPdebugMessage("exec method of event handler for logic or constraints\n");

   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_LBRELAXED )
   {
      SCIPdebugMessage("enabling constraint cons <%s> at depth %d\n", SCIPconsGetName((SCIP_CONS*)eventdata), SCIPgetDepth(scip));

      SCIP_CALL( SCIPenableCons(scip, (SCIP_CONS*)eventdata) );
      SCIP_CALL( SCIPenableConsPropagation(scip, (SCIP_CONS*)eventdata) );
   }
   else if( SCIPeventGetType(event) == SCIP_EVENTTYPE_UBTIGHTENED )
   {
      SCIP_CALL( SCIPenableConsPropagation(scip, (SCIP_CONS*)eventdata) );
   }

   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_VARFIXED )
   {
      SCIP_VAR* var = SCIPeventGetVar(event);

      /* we only catch this event in presolving stage */
      assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
      assert(var != NULL);

      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_FIXED && SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5 )
      {
	 SCIP_CONS* cons = (SCIP_CONS*)eventdata;

         if( SCIPconsIsActive(cons) )
	 {
	    SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
	    assert(consdata != NULL);
	    consdata->merged = FALSE;
	 }
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of conflict handler
 */

static
SCIP_DECL_CONFLICTEXEC(conflictExecLogicor)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int i;

   assert(conflicthdlr != NULL);
   assert(strcmp(SCIPconflicthdlrGetName(conflicthdlr), CONFLICTHDLR_NAME) == 0);
   assert(bdchginfos != NULL || nbdchginfos == 0);
   assert(result != NULL);

   /* don't process already resolved conflicts */
   if( resolved )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* create array of variables in conflict constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nbdchginfos) );
   for( i = 0; i < nbdchginfos; ++i )
   {
      assert(bdchginfos != NULL); /* for flexelint */
      assert(bdchginfos[i] != NULL);

      vars[i] = SCIPbdchginfoGetVar(bdchginfos[i]);

      /* we can only treat binary variables */
      if( !SCIPvarIsBinary(vars[i]) )
         break;

      /* if the variable is fixed to one in the conflict set, we have to use its negation */
      if( SCIPbdchginfoGetNewbound(bdchginfos[i]) > 0.5 )
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, vars[i], &vars[i]) );
      }
   }

   if( i == nbdchginfos )
   {
      SCIP_CONS* cons;
      char consname[SCIP_MAXSTRLEN];
      
      /* create a constraint out of the conflict set */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cf%d_%"SCIP_LONGINT_FORMAT, SCIPgetNRuns(scip), SCIPgetNConflictConssApplied(scip));
      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, consname, nbdchginfos, vars, 
            FALSE, separate, FALSE, FALSE, TRUE, local, FALSE, dynamic, removable, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, validnode) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      
      *result = SCIP_CONSADDED;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for logic or constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrLogicor(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONFLICTHDLR* conflicthdlr;
   SCIP_EVENTHDLR* eventhdlr;

   /* create event handler for events on watched variables */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecLogicor, NULL) );

   /* create conflict handler for logic or constraints */
   SCIP_CALL( SCIPincludeConflicthdlrBasic(scip, &conflicthdlr, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
         conflictExecLogicor, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpLogicor, consEnfopsLogicor, consCheckLogicor, consLockLogicor,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveLogicor) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyLogicor, consCopyLogicor) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveLogicor) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteLogicor) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreLogicor) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolLogicor) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeLogicor) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsLogicor) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsLogicor) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreLogicor) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpLogicor) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseLogicor) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolLogicor,CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintLogicor) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropLogicor, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropLogicor) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpLogicor, consSepasolLogicor, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransLogicor) );

   conshdlrdata->conshdlrlinear = SCIPfindConshdlr(scip,"linear");

   if( conshdlrdata->conshdlrlinear != NULL )
   {
      /* include the linear constraint to logicor constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdLogicor, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }

   /* logic or constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/logicor/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/logicor/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance",
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/logicor/dualpresolving",
         "should dual presolving steps be performed?",
         &conshdlrdata->dualpresolving, TRUE, DEFAULT_DUALPRESOLVING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/logicor/negatedclique",
         "should negated clique information be used in presolving",
         &conshdlrdata->usenegatedclique, TRUE, DEFAULT_NEGATEDCLIQUE, NULL, NULL) );

   return SCIP_OKAY;
}


/** creates and captures a logic or constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
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

   assert(scip != NULL);

   /* find the logicor constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("logic or constraint handler not found\n");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   SCIP_CALL( consdataCreate(scip, &consdata, nvars, vars) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   if( SCIPisTransformed(scip) && SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int v;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      for( v = consdata->nvars - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[v], SCIP_EVENTTYPE_VARFIXED, conshdlrdata->eventhdlr,
               (SCIP_EVENTDATA*)(*cons), NULL) );
      }
   }

   return SCIP_OKAY;
}

/** creates and captures a logicor constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsLogicor(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsLogicor() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars                /**< array with variables of constraint entries */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsLogicor(scip, cons, name, nvars, vars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** adds coefficient in logic or constraint */
SCIP_RETCODE SCIPaddCoefLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logicor constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   assert(var != NULL);

   /*debugMessage("adding variable <%s> to logicor constraint <%s>\n",
     SCIPvarGetName(var), SCIPconsGetName(cons));*/

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a logic or constraint\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( addCoef(scip, cons, var) );

   return SCIP_OKAY;
}

/** gets number of variables in logic or constraint */
int SCIPgetNVarsLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a logic or constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets array of variables in logic or constraint */
SCIP_VAR** SCIPgetVarsLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a logic or constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets the dual solution of the logic or constraint in the current LP */
SCIP_Real SCIPgetDualsolLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a logic or constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets the dual Farkas value of the logic or constraint in the current infeasible LP */
SCIP_Real SCIPgetDualfarkasLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a logic or constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}

/** returns the linear relaxation of the given logic or constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a logic or constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->row;
}

