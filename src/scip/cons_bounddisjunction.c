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

/**@file   cons_bounddisjunction.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for bound disjunction constraints
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "scip/cons_bounddisjunction.h"
#include "scip/cons_linear.h"
#include "scip/pub_misc.h"


#define CONSHDLR_NAME          "bounddisjunction"
#define CONSHDLR_DESC          "bound disjunction constraints"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -3000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -3000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "bounddisjunction"
#define EVENTHDLR_DESC         "event handler for bound disjunction constraints"

#define CONFLICTHDLR_NAME      "bounddisjunction"
#define CONFLICTHDLR_DESC      "conflict handler creating bound disjunction constraints"
#define CONFLICTHDLR_PRIORITY  -3000000


/* @todo make this a parameter setting */
#if 1 /* @todo test which AGEINCREASE formula is better! */
#define AGEINCREASE(n) (1.0 + 0.2*n)
#else
#define AGEINCREASE(n) (0.1*n)
#endif


/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for events on watched variables */
};

/** bound disjunction constraint data */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables of the literals in the constraint */
   SCIP_BOUNDTYPE*       boundtypes;         /**< types of bounds of the literals (lower or upper bounds) */
   SCIP_Real*            bounds;             /**< bounds of the literals */
   int                   varssize;           /**< size of vars, boundtypes, and bounds arrays */
   int                   nvars;              /**< number of variables in the constraint */
   int                   watchedvar1;        /**< position of the first watched variable */
   int                   watchedvar2;        /**< position of the second watched variable */
   int                   filterpos1;         /**< event filter position of first watched variable */
   int                   filterpos2;         /**< event filter position of second watched variable */
};




/*
 * Local methods
 */

/** adds rounding locks for the given variable in the given bound disjunction constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint */
   SCIP_CONSDATA*        consdata,           /**< bound disjunction constraint data */
   int                   pos                 /**< position of the variable in the constraint */
   )
{
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   if( consdata->boundtypes[pos] == SCIP_BOUNDTYPE_LOWER )
   {
      /* rounding down may violate the constraint */
      SCIP_CALL( SCIPlockVarCons(scip, consdata->vars[pos], cons, TRUE, FALSE) );
   }
   else
   {
      /* rounding up may violate the constraint */
      SCIP_CALL( SCIPlockVarCons(scip, consdata->vars[pos], cons, FALSE, TRUE) );
   }

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given bound disjunction constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint */
   SCIP_CONSDATA*        consdata,           /**< bound disjunction constraint data */
   int                   pos                 /**< position of the variable in the constraint */
   )
{
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   if( consdata->boundtypes[pos] == SCIP_BOUNDTYPE_LOWER )
   {
      /* rounding down may violate the constraint */
      SCIP_CALL( SCIPunlockVarCons(scip, consdata->vars[pos], cons, TRUE, FALSE) );
   }
   else
   {
      /* rounding up may violate the constraint */
      SCIP_CALL( SCIPunlockVarCons(scip, consdata->vars[pos], cons, FALSE, TRUE) );
   }

   return SCIP_OKAY;
}

/** catches the events on a single variable of the bound disjunction constraint */
static
SCIP_RETCODE catchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint */
   SCIP_CONSDATA*        consdata,           /**< bound disjunction constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos,                /**< position of the variable in the constraint */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   if( consdata->boundtypes[pos] == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)cons, filterpos) );
   }
   else
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)cons, filterpos) );
   }

   return SCIP_OKAY;
}

/** drops the events on a single variable of the bound disjunction constraint */
static
SCIP_RETCODE dropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint */
   SCIP_CONSDATA*        consdata,           /**< bound disjunction constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos,                /**< position of the variable in the constraint */
   int                   filterpos           /**< position of event filter entry returned by SCIPcatchVarEvent(), or -1 */
   )
{
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   if( consdata->boundtypes[pos] == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)cons, filterpos) );
   }
   else
   {
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED,
            eventhdlr, (SCIP_EVENTDATA*)cons, filterpos) );
   }

   return SCIP_OKAY;
}

/** creates constaint handler data for bound disjunction constraint handler */
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
      SCIPerrorMessage("event handler for bound disjunction constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   
   return SCIP_OKAY;
}

/** frees constraint handler data for bound disjunction constraint handler */
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

/** creates a bound disjunction constraint data object */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the bound disjunction constraint data */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< variables of the literals in the constraint */
   SCIP_BOUNDTYPE*       boundtypes,         /**< types of bounds of the literals (lower or upper bounds) */
   SCIP_Real*            bounds              /**< bounds of the literals */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || boundtypes != NULL);
   assert(nvars == 0 || bounds != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   if( nvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->boundtypes, boundtypes, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->bounds, bounds, nvars) );
      (*consdata)->varssize = nvars;
      (*consdata)->nvars = nvars;
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->boundtypes = NULL;
      (*consdata)->bounds = NULL;
      (*consdata)->varssize = 0;
      (*consdata)->nvars = 0;
   }
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->filterpos1 = -1;
   (*consdata)->filterpos2 = -1;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
   }

   return SCIP_OKAY;
}   

/** frees a bound disjunction constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to the bound disjunction constraint */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->boundtypes, (*consdata)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->bounds, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints bound disjunction constraint to file stream */
static
void consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< bound disjunction constraint data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             endline             /**< should an endline be set? */
   )
{
   int v;

   assert(consdata != NULL);

   /* print coefficients */
   SCIPinfoMessage(scip, file, "bounddisjunction(");
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      if( v > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "<%s> %s %.15g", SCIPvarGetName(consdata->vars[v]),
         consdata->boundtypes[v] == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", consdata->bounds[v]);
   }
   SCIPinfoMessage(scip, file, ")");

   if( endline )
      SCIPinfoMessage(scip, file, "\n");
}

/** stores the given variable numbers as watched variables, and updates the event processing */
static
SCIP_RETCODE switchWatchedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint */
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
      SCIP_CALL( dropEvents(scip, cons, consdata, eventhdlr, consdata->watchedvar1, consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 && consdata->watchedvar2 != watchedvar2 )
   {
      assert(consdata->filterpos2 != -1);
      SCIP_CALL( dropEvents(scip, cons, consdata, eventhdlr, consdata->watchedvar2, consdata->filterpos2) );
   }

   /* catch events on new watched variables */
   if( watchedvar1 != -1 && watchedvar1 != consdata->watchedvar1 )
   {
      SCIP_CALL( catchEvents(scip, cons, consdata, eventhdlr, watchedvar1, &consdata->filterpos1) );
   }
   if( watchedvar2 != -1 && watchedvar2 != consdata->watchedvar2 )
   {
      SCIP_CALL( catchEvents(scip, cons, consdata, eventhdlr, watchedvar2, &consdata->filterpos2) );
   }

   /* set the new watched variables */
   consdata->watchedvar1 = watchedvar1;
   consdata->watchedvar2 = watchedvar2;
   
   return SCIP_OKAY;
}

/** deletes coefficient at given position from bound disjunction constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint */
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
   SCIP_CALL( unlockRounding(scip, cons, consdata, pos) );

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

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->boundtypes[pos] = consdata->boundtypes[consdata->nvars-1];
   consdata->bounds[pos] = consdata->bounds[consdata->nvars-1];
   consdata->nvars--;

   /* if the last variable (that moved) was watched, update the watched position */
   if( consdata->watchedvar1 == consdata->nvars )
      consdata->watchedvar1 = pos;
   if( consdata->watchedvar2 == consdata->nvars )
      consdata->watchedvar2 = pos;

   SCIP_CALL( SCIPenableConsPropagation(scip, cons) );

   return SCIP_OKAY;
}

/** adds literal to bound disjunction constraint data */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_VAR*             var,                /**< variable in literal */
   SCIP_BOUNDTYPE        boundtype,          /**< boundtype of literal */
   SCIP_Real             bound               /**< bound of literal */
   )
{
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(var != NULL);
   assert(!SCIPisInfinity(scip, REALABS(bound)));
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* ensure enough memory in consdata arrays */
   if( consdata->varssize == consdata->nvars )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, consdata->nvars + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars,       consdata->varssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->boundtypes, consdata->varssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->bounds,     consdata->varssize, newsize) );
      consdata->varssize = newsize;
   }
   assert(consdata->varssize > consdata->nvars);

   /* add the variable to the end of the array */
   consdata->vars[consdata->nvars] = var;
   consdata->boundtypes[consdata->nvars] = boundtype;
   consdata->bounds[consdata->nvars] = bound;
   consdata->nvars++;

   if( SCIPconsIsTransformed(cons) )
   {
      /* add rounding lock of variable */
      SCIP_CALL( lockRounding(scip, cons, consdata, consdata->nvars-1) );

      /* if less than 2 variables are watched, add the new one to the watched variables */
      if( consdata->watchedvar1 == -1 )
      {
         assert(consdata->watchedvar2 == -1);
         SCIP_CALL( switchWatchedvars(scip, cons, eventhdlr, consdata->nvars-1, -1) );
      }
      else if( consdata->watchedvar2 == -1 )
      {
         SCIP_CALL( switchWatchedvars(scip, cons, eventhdlr, consdata->watchedvar1, consdata->nvars-1) );
      }
   }

   SCIP_CALL( SCIPenableConsPropagation(scip, cons) );

   return SCIP_OKAY;
}

/** deletes all variables with global bounds violating the literal, checks for global bounds satisfying the literal */
static
SCIP_RETCODE applyGlobalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            redundant           /**< returns whether a variable fixed to one exists in the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int v;
   SCIP_Real bnd;

   assert(eventhdlr != NULL);
   assert(redundant != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->vars != NULL);

   *redundant = FALSE;
   v = 0;
   while( v < consdata->nvars )
   {
      SCIP_VAR* var;

      var = consdata->vars[v];

      if( consdata->boundtypes[v] == SCIP_BOUNDTYPE_LOWER )
      {
         bnd = SCIPcomputeVarLbGlobal(scip, var);
         if( SCIPisFeasGE(scip, bnd, consdata->bounds[v]) )
         {
            *redundant = TRUE;
            return SCIP_OKAY;
         }
         else
         {
            bnd = SCIPcomputeVarUbGlobal(scip, var);
            if( SCIPisFeasLT(scip, bnd, consdata->bounds[v]) )
            {
               SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
            }
            else
               ++v;
         }
      }
      else
      {
         assert(consdata->boundtypes[v] == SCIP_BOUNDTYPE_UPPER);
         bnd = SCIPcomputeVarUbGlobal(scip, var);
         if( SCIPisFeasLE(scip, bnd, consdata->bounds[v]) )
         {
            *redundant = TRUE;
            return SCIP_OKAY;
         }
         else
         {
            bnd = SCIPcomputeVarLbGlobal(scip, var);
            if( SCIPisFeasGT(scip, bnd, consdata->bounds[v]) )
            {
               SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
            }
            else
               ++v;
         }
      }
   }

   SCIPdebugMessage("after global bounds: ");
   SCIPdebug(consdataPrint(scip, consdata, NULL, TRUE));

   return SCIP_OKAY;
}

/** replace variables by their representative active (or multi-aggregated) variables */
static
SCIP_RETCODE removeFixedVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_Bool*            redundant           /**< flag to indicate whether constraint has been bound redundant */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_BOUNDTYPE boundtype;
   SCIP_Real bound;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   v = 0;
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      assert(var != NULL);

      if( SCIPvarIsActive(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      {
         ++v;
         continue;
      }

      /* get active/fixed/multiaggr equivalent of v'th literal */
      bound = consdata->bounds[v];
      boundtype = consdata->boundtypes[v];
      SCIP_CALL( SCIPvarGetProbvarBound(&var, &bound, &boundtype) );

      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
      {
         /* if literal is satisfied, then constraint is redundant and we can stop */
         if( (boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisFeasLE(scip, bound, SCIPvarGetLbGlobal(var))) ||
             (boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisFeasGE(scip, bound, SCIPvarGetUbGlobal(var))) )
         {
            *redundant = TRUE;
            break;
         }
         /* if literal is not satisfied, then it just be removed */
      }
      else
      {
         /* add new literal */
         SCIP_CALL( addCoef(scip, cons, eventhdlr, var, boundtype, bound) );
      }

      /* remove old literal */
      SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
   }

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< bound disjunction constraint that detected the conflict */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* initialize conflict analysis, and add all bounds of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      /* the opposite bound is in conflict with this literal */
      SCIP_CALL( SCIPaddConflictBd(scip, consdata->vars[v], SCIPboundtypeOpposite(consdata->boundtypes[v]), NULL) );
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** returns whether literal at the given position is satisfied in the local bounds */
static
SCIP_Bool isLiteralSatisfied(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< bound disjunction constraint data */
   int                   pos                 /**< position of the literal */
   )
{
   SCIP_Real bnd;

   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   if( consdata->boundtypes[pos] == SCIP_BOUNDTYPE_LOWER )
   {
      bnd = SCIPcomputeVarLbLocal(scip, consdata->vars[pos]);
      return SCIPisFeasGE(scip, bnd, consdata->bounds[pos]);
   }
   else
   {
      bnd = SCIPcomputeVarUbLocal(scip, consdata->vars[pos]);
      return SCIPisFeasLE(scip, bnd, consdata->bounds[pos]);
   }
}

/** returns whether literal at the given position is violated in the local bounds */
static
SCIP_Bool isLiteralViolated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< bound disjunction constraint data */
   int                   pos                 /**< position of the literal */
   )
{
   SCIP_Real bnd;

   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   if( consdata->boundtypes[pos] == SCIP_BOUNDTYPE_LOWER )
   {
      bnd = SCIPcomputeVarUbLocal(scip, consdata->vars[pos]);
      return SCIPisFeasLT(scip, bnd, consdata->bounds[pos]);
   }
   else
   {
      bnd = SCIPcomputeVarLbLocal(scip, consdata->vars[pos]);
      return SCIPisFeasGT(scip, bnd, consdata->bounds[pos]);
   }
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
      SCIP_CALL( SCIPdisableCons(scip, cons) );
   }

   return SCIP_OKAY;
}

/** checks constraint for violation only looking at the watched variables, applies bound changes if possible */
static
SCIP_RETCODE processWatchedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint to be processed */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the constraint is infeasible in current bounds */
   SCIP_Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   SCIP_Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;
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
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->watchedvar1 == -1 || consdata->watchedvar1 != consdata->watchedvar2);

   *mustcheck = FALSE;

   SCIPdebugMessage("processing watched variables of constraint <%s>\n", SCIPconsGetName(cons));

   nvars = consdata->nvars;
   vars = consdata->vars;
   boundtypes = consdata->boundtypes;
   bounds = consdata->bounds;
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || boundtypes != NULL);
   assert(nvars == 0 || bounds != NULL);

   /* check watched variables if they are satisfying the literal */
   if( consdata->watchedvar1 >= 0 && isLiteralSatisfied(scip, consdata, consdata->watchedvar1) )
   {
      /* the literal is satisfied, making the constraint redundant */
      SCIPdebugMessage(" -> disabling constraint <%s> (watchedvar1 satisfied)\n", SCIPconsGetName(cons));
      SCIP_CALL( disableCons(scip, cons) );
      return SCIP_OKAY;
   }
   if( consdata->watchedvar2 >= 0 && isLiteralSatisfied(scip, consdata, consdata->watchedvar2) )
   {
      /* the literal is satisfied, making the constraint redundant */
      SCIPdebugMessage(" -> disabling constraint <%s> (watchedvar2 satisfied)\n", SCIPconsGetName(cons));
      SCIP_CALL( disableCons(scip, cons) );
      return SCIP_OKAY;
   }

   /* check if watched variables are still undecided */
   watchedvar1 = -1;
   watchedvar2 = -1;
   nbranchings1 = SCIP_LONGINT_MAX;
   nbranchings2 = SCIP_LONGINT_MAX;
   if( consdata->watchedvar1 >= 0 && !isLiteralViolated(scip, consdata, consdata->watchedvar1) )
   {
      watchedvar1 = consdata->watchedvar1;
      nbranchings1 = -1; /* prefer keeping the watched variable */
   }
   if( consdata->watchedvar2 >= 0 && !isLiteralViolated(scip, consdata, consdata->watchedvar2) )
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

         /* check, if the literal is violated */
         if( isLiteralViolated(scip, consdata, v) )
            continue;

         /* check, if the literal is satisfied */
         if( isLiteralSatisfied(scip, consdata, v) )
         {
            assert(v != consdata->watchedvar1);
            assert(v != consdata->watchedvar2);
            
            /* the literal is satisfied, making the constraint redundant;
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
      
         /* the literal is still undecided and can be used as watched variable */
         nbranchings = SCIPvarGetNBranchingsCurrentRun(vars[v], 
            boundtypes[v] == SCIP_BOUNDTYPE_LOWER ? SCIP_BRANCHDIR_DOWNWARDS : SCIP_BRANCHDIR_UPWARDS);
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
      /* there is no undecided literal left -> the constraint is infeasible
       *  - a modifiable constraint is infeasible
       *  - an unmodifiable constraint is infeasible and the node can be cut off
       */
      assert(watchedvar2 == -1);

      SCIPdebugMessage(" -> constraint <%s> is infeasible\n", SCIPconsGetName(cons));
      *infeasible = TRUE;

      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      if( !SCIPconsIsModifiable(cons) )
      {
         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflict(scip, cons) );

         /* mark the node to be cut off */
         *cutoff = TRUE;
      }
   }
   else if( watchedvar2 == -1 )
   {
      /* there is only one undecided literal:
       * - a modifiable constraint must be checked manually
       * - we cannot change bounds of multi-aggregated variables and have to check manually
       * - an unmodifiable constraint is feasible and can be disabled after the remaining literal is satisfied
       */
      assert(0 <= watchedvar1 && watchedvar1 < nvars);
      assert(!isLiteralViolated(scip, consdata, watchedvar1));
      assert(!isLiteralSatisfied(scip, consdata, watchedvar1));
      if( SCIPconsIsModifiable(cons)
         || SCIPvarGetStatus(SCIPvarGetProbvar(vars[watchedvar1])) == SCIP_VARSTATUS_MULTAGGR )
         *mustcheck = TRUE;
      else
      {
         SCIP_Bool infbdchg;

         /* satisfy remaining literal and disable constraint; make sure, the fixed-to-one variable is watched */
         SCIPdebugMessage(" -> single-literal constraint <%s> (change bound <%s> %s %g) at depth %d\n", 
            SCIPconsGetName(cons), SCIPvarGetName(vars[watchedvar1]), 
            boundtypes[watchedvar1] == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", bounds[watchedvar1], SCIPgetDepth(scip));
         if( boundtypes[watchedvar1] == SCIP_BOUNDTYPE_LOWER )
         {
            SCIP_CALL( SCIPinferVarLbCons(scip, vars[watchedvar1], bounds[watchedvar1], cons, watchedvar1, FALSE,
                  &infbdchg, NULL) );
         }
         else
         {
            SCIP_CALL( SCIPinferVarUbCons(scip, vars[watchedvar1], bounds[watchedvar1], cons, watchedvar1, FALSE,
                  &infbdchg, NULL) );
         }
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
      SCIPdebugMessage(" -> new watched variables <%s> and <%s> of constraint <%s> are still undecided\n",
         SCIPvarGetName(vars[watchedvar1]), SCIPvarGetName(vars[watchedvar2]), SCIPconsGetName(cons));

      /* switch to the new watched variables */
      SCIP_CALL( switchWatchedvars(scip, cons, eventhdlr, watchedvar1, watchedvar2) );

      /* there are at least two undecided variables -> the constraint must be checked manually */
      *mustcheck = TRUE;

      /* disable propagation of constraint until the corresponding bound of a watched variable changed */
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
   SCIP_CONS*            cons,               /**< bound disjunction constraint to be checked */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            violated            /**< pointer to store whether the given solution violates the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;
   SCIP_Real solval;
   int nvars;
   int v;

   assert(violated != NULL);

   *violated = FALSE;
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   vars = consdata->vars;
   boundtypes = consdata->boundtypes;
   bounds = consdata->bounds;
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || boundtypes != NULL);
   assert(nvars == 0 || bounds != NULL);

   /* check the given solution */
   *violated = TRUE;
   for( v = 0; v < nvars; ++v )
   {
      solval = SCIPgetSolVal(scip, sol, vars[v]);
      if( (boundtypes[v] == SCIP_BOUNDTYPE_LOWER && SCIPisFeasGE(scip, solval, bounds[v]))
         || (boundtypes[v] == SCIP_BOUNDTYPE_UPPER && SCIPisFeasLE(scip, solval, bounds[v])) )
      {
         *violated = FALSE;
         break;
      }
   }
   return SCIP_OKAY;
}

/* registers variables of a constraint as branching candidates 
 * indicates whether an n-ary branch is necessary to enforce this constraint, 
 * because all active literals are w.r.t. continuous variables which bound (in the literal) is at the variable's bound 
 */
static
SCIP_RETCODE registerBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint which variables should be registered for branching */
   SCIP_Bool*            neednarybranch      /**< pointer to store TRUE, if n-ary branching is necessary to enforce this constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;
   SCIP_Real violation;
   SCIP_Real varlb;
   SCIP_Real varub;
   int nvars;
   int v;
   
   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(neednarybranch != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   nvars = consdata->nvars;
   vars = consdata->vars;
   boundtypes = consdata->boundtypes;
   bounds = consdata->bounds;
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || boundtypes != NULL);
   assert(nvars == 0 || bounds != NULL);
   
   *neednarybranch = TRUE;
   
   for( v = 0; v < nvars; ++v )
   {
      /* constraint should be violated, so all bounds in the constraint have to be violated */
      assert( !(boundtypes[v] == SCIP_BOUNDTYPE_LOWER && SCIPisFeasGE(scip, SCIPgetSolVal(scip, NULL, vars[v]), bounds[v])) &&
         !(boundtypes[v] == SCIP_BOUNDTYPE_UPPER && SCIPisFeasLE(scip, SCIPgetSolVal(scip, NULL, vars[v]), bounds[v])) );

      varlb = SCIPcomputeVarLbLocal(scip, vars[v]);
      varub = SCIPcomputeVarUbLocal(scip, vars[v]);
      /* if literal is x >= varlb, but upper bound on x is < varlb, then this literal can never be satisfied,
       * thus there is no use for branching */
      if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER && SCIPisFeasLT(scip, varub, bounds[v]) )
         continue;
      /* if literal is x <= varub, but lower bound on x is > varub, then this literal can never be satisfied,
       * thus there is no use for branching */
      if( boundtypes[v] == SCIP_BOUNDTYPE_UPPER && SCIPisFeasGT(scip, varlb, bounds[v]) )
         continue;

      violation = SCIPgetSolVal(scip, NULL, vars[v]) - bounds[v];

      /* if variable is continuous, then we cannot branch on one of the variable bounds */
      if( SCIPvarGetType(vars[v]) != SCIP_VARTYPE_CONTINUOUS ||
         ((SCIPisInfinity(scip, -varlb) || !SCIPisFeasEQ(scip, bounds[v], varlb)) &&
          (SCIPisInfinity(scip,  varub) || !SCIPisFeasEQ(scip, bounds[v], varub))) )
      {
         SCIP_CALL( SCIPaddExternBranchCand(scip, vars[v], REALABS(violation), bounds[v]) );
         *neednarybranch = FALSE;
      }
   }
   
   return SCIP_OKAY;
}

/** enforces the pseudo or LP solution on the given constraint */
static
SCIP_RETCODE enforceCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< bound disjunction constraint to be separated */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the constraint was infeasible */
   SCIP_Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   SCIP_Bool*            registeredbrcand    /**< pointer to store TRUE, if branching variable candidates were registered or was already true */
   )
{
   SCIP_Bool mustcheck;
   SCIP_Bool neednarybranch;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(reduceddom != NULL);
   assert(registeredbrcand != NULL);

   /* update and check the watched variables, if they were changed since last processing */
   if( SCIPconsIsPropagationEnabled(cons) )
   {
      SCIP_CALL( processWatchedVars(scip, cons, eventhdlr, cutoff, infeasible, reduceddom, &mustcheck) );
   }
   else
      mustcheck = TRUE;

   if( mustcheck )
   {
      SCIP_Bool violated;

      SCIP_CALL( checkCons(scip, cons, NULL, &violated) );
      if( violated )
      {
         /* constraint was infeasible -> reset age */
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *infeasible = TRUE;
         
         /* register branching candidates */
         SCIP_CALL( registerBranchingCandidates(scip, cons, &neednarybranch) );
         
         if( !neednarybranch )
            *registeredbrcand = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** enforces a constraint by creating an n-ary branch consisting of a set of child nodes, each enforcing one literal
 */
static
SCIP_RETCODE createNAryBranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< bound disjunction constraint to branch on */
)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;
   SCIP_Real varlb;
   SCIP_Real varub;
   int nvars;
   int v;
   
   SCIP_Real  priority;
   SCIP_Real  estimate;
   SCIP_NODE* node;
   
   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   nvars = consdata->nvars;
   vars = consdata->vars;
   boundtypes = consdata->boundtypes;
   bounds = consdata->bounds;
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || boundtypes != NULL);
   assert(nvars == 0 || bounds != NULL);
   
   for( v = 0; v < nvars; ++v )
   {
      /* constraint should be violated, so all bounds in the constraint have to be violated */
      assert( !(boundtypes[v] == SCIP_BOUNDTYPE_LOWER && SCIPisFeasGE(scip, SCIPgetSolVal(scip, NULL, vars[v]), bounds[v])) &&
         !(boundtypes[v] == SCIP_BOUNDTYPE_UPPER && SCIPisFeasLE(scip, SCIPgetSolVal(scip, NULL, vars[v]), bounds[v])) );

      varlb = SCIPcomputeVarLbLocal(scip, vars[v]);
      varub = SCIPcomputeVarUbLocal(scip, vars[v]);
      /* if literal is x >= varlb, but upper bound on x is < varlb, then this literal can never be satisfied,
       * thus there is no use in creating an extra child for it */
      if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER && SCIPisFeasLT(scip, varub, bounds[v]) )
         continue;
      /* if literal is x <= varub, but lower bound on x is > varub, then this literal can never be satisfied,
       * thus there is no use in creating an extra child for it */
      if( boundtypes[v] == SCIP_BOUNDTYPE_UPPER && SCIPisFeasGT(scip, varlb, bounds[v]) )
         continue;

      /* create a child that enforces the current literal */      
      priority = SCIPcalcNodeselPriority(scip, vars[v], boundtypes[v] == SCIP_BOUNDTYPE_LOWER ? 
         SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS, bounds[v]);
      estimate = SCIPcalcChildEstimate  (scip, vars[v], bounds[v]);

      SCIPdebugMessage(" -> creating child to enforce: <%s> %c= %g (priority: %g, estimate: %g)\n",
         SCIPvarGetName(vars[v]), boundtypes[v] == SCIP_BOUNDTYPE_LOWER ? '>' : '<', bounds[v], priority, estimate);
      
      SCIP_CALL( SCIPcreateChild(scip, &node, priority, estimate) );

      /* enforce current literal */
      if( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_MULTAGGR )
      {
         SCIP_CONS* brcons;
         SCIP_Real  one;

         one = 1.0;
         
         if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER )
         {
            SCIP_CALL( SCIPcreateConsLinear(scip, &brcons, "bounddisjbranch", 1, &vars[v], &one, bounds[v], SCIPinfinity(scip),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
               SCIPconsIsStickingAtNode(cons)) );
         }
         else
         {
            SCIP_CALL( SCIPcreateConsLinear(scip, &brcons, "bounddisjbranch", 1, &vars[v], &one, -SCIPinfinity(scip), bounds[v],
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
               SCIPconsIsStickingAtNode(cons)) );
         }
         SCIP_CALL( SCIPaddConsNode(scip, node, brcons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &brcons) );
      }
      else
      {
         assert(SCIPvarIsActive(vars[v]));
         if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER )
         {
            SCIP_CALL( SCIPchgVarLbNode(scip, node, vars[v], bounds[v]) );
         }
         else
         {
            SCIP_CALL( SCIPchgVarUbNode(scip, node, vars[v], bounds[v]) );
         }
      }
      
      /* delete bound disjunction constraint from child node */
      SCIP_CALL( SCIPdelConsNode(scip, node, cons) );      
   }
   
   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyBounddisjunction)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBounddisjunction(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeBounddisjunction)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitBounddisjunction NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitBounddisjunction NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreBounddisjunction NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreBounddisjunction)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Bool redundant;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* fast processing of constraints, apply global bounds and remove fixed variables */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      SCIPdebugMessage("exit-presolving bound disjunction constraint <%s>\n", SCIPconsGetName(cons));

      /* remove all literals that are violated in global bounds, check redundancy due to global bounds */
      SCIP_CALL( applyGlobalBounds(scip, cons, conshdlrdata->eventhdlr, &redundant) );

      if( !redundant )
      {
         /* replace variables by their representative active (or multi-aggregated) variables */
         SCIP_CALL( removeFixedVariables(scip, cons, conshdlrdata->eventhdlr, &redundant) );
      }

      if( redundant )
      {
         SCIPdebugMessage("bound disjunction constraint <%s> is redundant\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );
         continue;
      }
      else if( !SCIPconsIsModifiable(cons) && consdata->nvars == 0 )
      {
         /* if unmodifiable constraint has no variables, it is infeasible */
         SCIPdebugMessage("bound disjunction constraint <%s> is infeasible\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolBounddisjunction NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#define consExitsolBounddisjunction NULL


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteBounddisjunction)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* free LP row and bound disjunction constraint */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransBounddisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   /*debugMessage("Trans method of bound disjunction constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->nvars, sourcedata->vars,
         sourcedata->boundtypes, sourcedata->bounds) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
#define consInitlpBounddisjunction NULL


/** separation method of constraint handler for LP solutions */
#define consSepalpBounddisjunction NULL


/** separation method of constraint handler for arbitrary primal solutions */
#define consSepasolBounddisjunction NULL


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBounddisjunction)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_Bool reduceddom;
   SCIP_Bool registeredbrcand;
   SCIP_Bool infeasiblecons;
   int c;
   int nnarybranchconsvars;
   SCIP_CONS* narybranchcons; /* constraint that is a candidate for an n-ary branch */

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("LP enforcing %d bound disjunction constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   infeasible = FALSE;
   reduceddom = FALSE;
   registeredbrcand = FALSE;
   narybranchcons = NULL;
   nnarybranchconsvars = INT_MAX;

   /* check all bound disjunction constraints for feasibility */
   for( c = 0; c < nconss && !cutoff && !reduceddom; ++c )
   {
      infeasiblecons = FALSE;
      SCIP_CALL( enforceCurrentSol(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &infeasiblecons, &reduceddom, &registeredbrcand) );
      infeasible |= infeasiblecons;
      if( infeasiblecons && !registeredbrcand )
      {
         /* if cons. c has less literals than the previous candidate for an n-ary branch, then keep cons. c as candidate for n-ary branch */
         if( narybranchcons == NULL || SCIPconsGetData(conss[c])->nvars < nnarybranchconsvars )
         {
            narybranchcons = conss[c];
            nnarybranchconsvars = SCIPconsGetData(narybranchcons)->nvars;
            assert(nnarybranchconsvars > 0);
         }
      }
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else if( infeasible )
   {
      if( registeredbrcand )
      {
         *result = SCIP_INFEASIBLE;
      }
      else
      {
         SCIP_CALL( createNAryBranch(scip, narybranchcons) );
         *result = SCIP_BRANCHED;
      }
   }
   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBounddisjunction)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_Bool reduceddom;
   SCIP_Bool registeredbrcand;
   int c;
   SCIP_CONS* narybranchcons; /* constraint that is a candidate for an n-ary branch */

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("pseudo enforcing %d bound disjunction constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   infeasible = FALSE;
   reduceddom = FALSE;
   registeredbrcand = FALSE;
   narybranchcons = NULL;

   /* check all bound disjunction constraints for feasibility */
   for( c = 0; c < nconss && !cutoff && !reduceddom; ++c )
   {
      SCIP_CALL( enforceCurrentSol(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &infeasible, &reduceddom, &registeredbrcand) );
      if( infeasible && !registeredbrcand )
      {
         /* if cons. c has less literals than the previous candidate for an n-ary branch, then keep cons. c as candidate for n-ary branch */
         if( !narybranchcons || SCIPconsGetData(conss[c])->nvars < SCIPconsGetData(narybranchcons)->nvars )
            narybranchcons = conss[c];
      }
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reduceddom )
      *result = SCIP_REDUCEDDOM;
   else if( infeasible )
   {
      if( registeredbrcand )
      {
         *result = SCIP_INFEASIBLE;
      }
      else
      {
         SCIP_CALL( createNAryBranch(scip, narybranchcons) );
         *result = SCIP_BRANCHED;
      }
   }
   
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckBounddisjunction)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* check all bound disjunction constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_Bool violated;

      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      SCIP_CALL( checkCons(scip, cons, sol, &violated) );
      if( violated )
      {
         if( printreason )
         {
            int v;

            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "violation: ");
            for( v = 0; v < consdata->nvars; ++v )
            {
               assert(consdata->vars[v] != NULL);
               if( v > 0 )
                  SCIPinfoMessage(scip, NULL, ", ");
               SCIPinfoMessage(scip, NULL, "<%s> = %.15g", 
                  SCIPvarGetName(consdata->vars[v]), SCIPgetSolVal(scip, sol, consdata->vars[v]));
            }
            SCIPinfoMessage(scip, NULL, ")\n");
         }

         /* constraint is violated */
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }
   
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropBounddisjunction)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool infeasible;
   SCIP_Bool reduceddom;
   SCIP_Bool mustcheck;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   infeasible = FALSE;
   reduceddom = FALSE;

   /* propagate all useful bound disjunction constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( processWatchedVars(scip, conss[c], conshdlrdata->eventhdlr,
            &cutoff, &infeasible, &reduceddom, &mustcheck) );
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
SCIP_DECL_CONSPRESOL(consPresolBounddisjunction)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   SCIP_Bool tightened;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   for( c = 0; c < nconss && *result != SCIP_CUTOFF && !SCIPisStopped(scip); ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      SCIPdebugMessage("presolving bound disjunction constraint <%s>\n", SCIPconsGetName(cons));

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
      {
         SCIP_CALL( SCIPenableConsPropagation(scip, cons) );
      }

      /* remove all literals that are violated in global bounds, check redundancy due to global bounds */
      SCIP_CALL( applyGlobalBounds(scip, cons, conshdlrdata->eventhdlr, &redundant) );

      if( !redundant )
      {
         /* replace variables by their representative active (or multi-aggregated) variables */
         SCIP_CALL( removeFixedVariables(scip, cons, conshdlrdata->eventhdlr, &redundant) );
      }

      /**@todo find pairs of negated variables in constraint: constraint is redundant */
      /**@todo find sets of equal variables in constraint: multiple entries of variable can be replaced by single entry */

      if( redundant )
      {
         SCIPdebugMessage("bound disjunction constraint <%s> is redundant\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss)++;
         *result = SCIP_SUCCESS;
         continue;
      }
      else if( !SCIPconsIsModifiable(cons) )
      {
         /* if unmodifiable constraint has no variables, it is infeasible,
          * if unmodifiable constraint has only one variable, the literal can be satisfied and the constraint deleted
          */
         if( consdata->nvars == 0 )
         {
            SCIPdebugMessage("bound disjunction constraint <%s> is infeasible\n", SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if( consdata->nvars == 1 )
         {
            SCIPdebugMessage("bound disjunction constraint <%s> has only one undecided literal\n",
               SCIPconsGetName(cons));
            
            assert(consdata->vars != NULL);
            assert(!isLiteralSatisfied(scip, consdata, 0));
            assert(!isLiteralViolated(scip, consdata, 0));

            if( SCIPvarIsActive(consdata->vars[0]) )
            {
               if( consdata->boundtypes[0] == SCIP_BOUNDTYPE_LOWER )
               {
                  SCIP_CALL( SCIPtightenVarLb(scip, consdata->vars[0], consdata->bounds[0], TRUE, &infeasible, &tightened) );
               }
               else
               {
                  SCIP_CALL( SCIPtightenVarUb(scip, consdata->vars[0], consdata->bounds[0], TRUE, &infeasible, &tightened) );
               }
               if( infeasible )
               {
                  SCIPdebugMessage(" -> infeasible fixing\n");
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               assert(tightened);
               (*nchgbds)++;
            }
            else
            {
               /* upgrade to a linear constraint, if vars[0] is multiaggregated */
               SCIP_CONS* lincons;
               SCIP_Real one;

               assert(SCIPvarGetStatus(consdata->vars[0]) == SCIP_VARSTATUS_MULTAGGR);

               one = 1.0;
               if( consdata->boundtypes[0] == SCIP_BOUNDTYPE_LOWER )
               {
                  SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons),
                     1, &consdata->vars[0], &one, consdata->bounds[0], SCIPinfinity(scip),
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
                     SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                     SCIPconsIsStickingAtNode(cons)) );
               }
               else
               {
                  SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons),
                     1, &consdata->vars[0], &one, -SCIPinfinity(scip), consdata->bounds[0],
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
                     SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
                     SCIPconsIsStickingAtNode(cons)) );
               }
               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
               (*nupgdconss)++;
            }

            SCIP_CALL( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }
      }
   }
   
   /**@todo preprocess pairs of bound disjunction constraints */

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropBounddisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_BOUNDTYPE* boundtypes;
#ifndef NDEBUG
   SCIP_Real* bounds;
#endif
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(infervar != NULL);
   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(0 <= inferinfo && inferinfo < consdata->nvars);
   assert(consdata->vars[inferinfo] == infervar);

   vars = consdata->vars;
   boundtypes = consdata->boundtypes;
#ifndef NDEBUG
   bounds = consdata->bounds;
   assert(bounds != NULL);
#endif
   assert(boundtypes != NULL);

   SCIPdebugMessage("conflict resolving method of bound disjunction constraint handler\n");

   /* the only deductions are bounds tightened to a literal's bound on bound disjunction constraints where all other
    * literals are violated
    */
   assert((boundtypes[inferinfo] == SCIP_BOUNDTYPE_LOWER
         && SCIPisFeasGE(scip, SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE), bounds[inferinfo]))
      || (boundtypes[inferinfo] == SCIP_BOUNDTYPE_UPPER
         && SCIPisFeasLE(scip, SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE), bounds[inferinfo])));

   for( v = 0; v < consdata->nvars; ++v )
   {
      if( v != inferinfo )
      {
         assert(consdata->vars[v] != infervar || consdata->boundtypes[v] != consdata->boundtypes[inferinfo]);

         /* the reason literal must have been violated
          * we do not check for multiaggregated variables, since SCIPvarGetXbAtIndex is not implemented for multiaggr. variables */
         assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_MULTAGGR
            || (boundtypes[v] == SCIP_BOUNDTYPE_LOWER
               && SCIPisFeasLT(scip, SCIPvarGetUbAtIndex(vars[v], bdchgidx, TRUE), bounds[v]))
            || (boundtypes[v] == SCIP_BOUNDTYPE_UPPER
               && SCIPisFeasGT(scip, SCIPvarGetLbAtIndex(vars[v], bdchgidx, TRUE), bounds[v])));
         SCIP_CALL( SCIPaddConflictBd(scip, vars[v], SCIPboundtypeOpposite(boundtypes[v]), bdchgidx) );
      }
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBounddisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* lock every single coefficient */
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( consdata->boundtypes[i] == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos, nlocksneg) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveBounddisjunction)
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

   SCIPdebugMessage("activating information for bound disjunction constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebug(consdataPrint(scip, consdata, NULL, TRUE));

   /* catch events on watched variables */
   if( consdata->watchedvar1 != -1 )
   {
      SCIP_CALL( catchEvents(scip, cons, consdata, conshdlrdata->eventhdlr, consdata->watchedvar1,
            &consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 )
   {
      SCIP_CALL( catchEvents(scip, cons, consdata, conshdlrdata->eventhdlr, consdata->watchedvar2,
            &consdata->filterpos2) );
   }

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveBounddisjunction)
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

   SCIPdebugMessage("deactivating information for bound disjunction constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebug(consdataPrint(scip, consdata, NULL, TRUE));

   /* drop events on watched variables */
   if( consdata->watchedvar1 != -1 )
   {
      assert(consdata->filterpos1 != -1);
      SCIP_CALL( dropEvents(scip, cons, consdata, conshdlrdata->eventhdlr, consdata->watchedvar1, consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 )
   {
      assert(consdata->filterpos2 != -1);
      SCIP_CALL( dropEvents(scip, cons, consdata, conshdlrdata->eventhdlr, consdata->watchedvar2, consdata->filterpos2) );
   }

   return SCIP_OKAY;
}


/** constraint enabling notification method of constraint handler */
#define consEnableBounddisjunction NULL


/** constraint disabling notification method of constraint handler */
#define consDisableBounddisjunction NULL


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintBounddisjunction)
{  /*lint --e{715}*/

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   consdataPrint(scip, SCIPconsGetData(cons), file, FALSE);
  
   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyBounddisjunction)
{
   SCIP_VAR** sourcevars;
   SCIP_VAR** targetvars;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;

   int nvars;
   int v;
   
   assert(valid != NULL);
   
   *valid = TRUE;

   /* get source data */
   sourcevars = SCIPgetVarsBounddisjunction(sourcescip, sourcecons);
   nvars = SCIPgetNVarsBounddisjunction(sourcescip, sourcecons);
   boundtypes = SCIPgetBoundtypesBounddisjunction(sourcescip, sourcecons);
   bounds = SCIPgetBoundsBounddisjunction(sourcescip, sourcecons);

   SCIP_CALL( SCIPallocBufferArray(scip, &targetvars, nvars) );
   
   /* map source variables to active variables of the target SCIP */
   for( v = 0; v < nvars && *valid; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &targetvars[v], varmap, consmap, global, valid) );
      assert(!(*valid) || targetvars[v] != NULL);
   }

   /* only create the target constraint, if all variables could be copied */
   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, cons, name ? name : SCIPconsGetName(sourcecons), nvars, targetvars, boundtypes,
         bounds, initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );      
   }

   SCIPfreeBufferArray(scip, &targetvars);

   return SCIP_OKAY;   
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseBounddisjunction)
{  /*lint --e{715}*/
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;
   SCIP_VAR** vars;
   const char* s;
   char* t;
   int varssize;
   int nvars;
   int pos;

   assert( success != NULL );
   *success = TRUE;

   SCIPdebugMessage("parse <%s> as bounddisjunction constraint\n", str);

   /* skip white space */
   s = str;
   while ( *s != '\0' && isspace(*s) )
      ++s;

   /* check for string "bounddisjunction" */
   if ( strncmp(s, "bounddisjunction(", 16) != 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "error during parsing: expected \"bounddisjunction(\" in <%s>.\n", s);
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* skip "bounddisjunction(" */
   s += 17;

   varssize = 100;
   nvars = 0;

   /* allocate buffer array for variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, varssize) );

   /* parse string until ")" */
   while ( *s != '\0' && *s != ')' )
   {
      SCIP_VAR* var;

      /* parse variable name */ 
      pos = 0;
      SCIP_CALL( SCIPparseVarName(scip, s, pos, &var, &pos) );

      if ( var == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "variable with name <%s> does not exist\n", SCIPvarGetName(var));
         *success = FALSE;
         goto TERMINATE;
      }

      /* skip white space */
      s = &s[pos];
      while ( *s != '\0' && isspace(*s) && *s != '>' && *s != '<' )
         ++s;

      /* parse bound type */
      switch ( *s )
      {
      case '<':
         boundtypes[nvars] = SCIP_BOUNDTYPE_UPPER;
         break;
      case '>':
         boundtypes[nvars] = SCIP_BOUNDTYPE_LOWER;
         break;
      default:
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "variable with name <%s> does not exist\n", SCIPvarGetName(var));
         *success = FALSE;
         goto TERMINATE;
      }

      ++s;
      if ( *s != '=' )
      {
         SCIPdebugMessage("expected '=': %s\n", s);
         *success = FALSE;
         goto TERMINATE;
      }

      /* skip '=' */
      ++s;

      /* skip white space */
      while ( *s != '\0' && isspace(*s) )
         ++s;

      /* parse bound value */
      bounds[nvars] = strtod(s, &t);
      if ( t == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error during parsing of the weight: %s\n", s);
         *success = FALSE;
         goto TERMINATE;
      }

      /* skip white space */
      s = t;
      while ( (*s != '\0' && isspace(*s)) || *s == ',' )
         ++s;

      /* set variable */
      vars[nvars++] = var;
      
      /* check if the size of the variable array was big enough */
      if ( nvars > varssize )
      {
         /* reallocate memory */
         varssize *= 2;
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, varssize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &boundtypes, varssize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &bounds, varssize) );
      }
   }
   /* ignore if the string ended without ")" */

   /* add bounddisjunction */
   if ( *success && nvars > 0 )
   {
      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, cons, name, nvars, vars, boundtypes, bounds, 
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

 TERMINATE:
   /* free variable buffer */
   SCIPfreeBufferArray(scip, &bounds);
   SCIPfreeBufferArray(scip, &boundtypes);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecBounddisjunction)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   /*SCIPdebugMessage("exec method of event handler for bound disjunction constraints\n");*/

   if( (SCIPeventGetType(event) & SCIP_EVENTTYPE_BOUNDRELAXED) != 0 )
   {
      SCIP_CALL( SCIPenableCons(scip, (SCIP_CONS*)eventdata) );
   }
   else
      assert((SCIPeventGetType(event) & SCIP_EVENTTYPE_BOUNDTIGHTENED) != 0);

   SCIP_CALL( SCIPenableConsPropagation(scip, (SCIP_CONS*)eventdata) );

   return SCIP_OKAY;
}




/*
 * Callback methods of conflict handler
 */

static
SCIP_DECL_CONFLICTEXEC(conflictExecBounddisjunction)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;
   SCIP_CONS* cons;
   char consname[SCIP_MAXSTRLEN];
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

   /* create array of variables, boundtypes, and bound values in conflict constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nbdchginfos) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, nbdchginfos) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, nbdchginfos) );
   for( i = 0; i < nbdchginfos; ++i )
   {
      assert(bdchginfos != NULL);

      vars[i] = SCIPbdchginfoGetVar(bdchginfos[i]);
      boundtypes[i] = SCIPboundtypeOpposite(SCIPbdchginfoGetBoundtype(bdchginfos[i]));
      bounds[i] = SCIPbdchginfoGetNewbound(bdchginfos[i]);
      /* for continuous variables, we can only use the relaxed version of the bounds negation: !(x <= u) -> x >= u */
      if( SCIPvarIsIntegral(vars[i]) )
      {
         assert(SCIPisIntegral(scip, bounds[i]));
         bounds[i] += (boundtypes[i] == SCIP_BOUNDTYPE_LOWER ? +1.0 : -1.0);
      }
      else if( (boundtypes[i] == SCIP_BOUNDTYPE_LOWER && SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(vars[i]), bounds[i]))
         || (boundtypes[i] == SCIP_BOUNDTYPE_UPPER && SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(vars[i]), bounds[i])) )
      {
         /* the literal is satisfied in global bounds (may happen due to weak "negation" of continuous variables)
          * -> discard the conflict constraint
          */
         break;
      }
   }
      
   /* create a constraint out of the conflict set */
   if( i == nbdchginfos )
   {
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cf%d_%"SCIP_LONGINT_FORMAT, SCIPgetNRuns(scip), SCIPgetNConflictConssApplied(scip));
      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, consname, nbdchginfos, vars, boundtypes, bounds,
            FALSE, TRUE, FALSE, FALSE, TRUE, local, FALSE, dynamic, removable, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, validnode) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      *result = SCIP_CONSADDED;
   }
   else
      *result = SCIP_DIDNOTFIND;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &bounds);
   SCIPfreeBufferArray(scip, &boundtypes);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for bound disjunction constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBounddisjunction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create event handler for events on watched variables */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, eventExecBounddisjunction, NULL) );

   /* create conflict handler for bound disjunction constraints */
   SCIP_CALL( SCIPincludeConflicthdlr(scip, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
         NULL, NULL, NULL, NULL, NULL, NULL, conflictExecBounddisjunction, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyBounddisjunction,
         consFreeBounddisjunction, consInitBounddisjunction, consExitBounddisjunction, 
         consInitpreBounddisjunction, consExitpreBounddisjunction, consInitsolBounddisjunction, consExitsolBounddisjunction,
         consDeleteBounddisjunction, consTransBounddisjunction, 
         consInitlpBounddisjunction, consSepalpBounddisjunction, consSepasolBounddisjunction, 
         consEnfolpBounddisjunction, consEnfopsBounddisjunction, consCheckBounddisjunction, 
         consPropBounddisjunction, consPresolBounddisjunction, consRespropBounddisjunction, consLockBounddisjunction,
         consActiveBounddisjunction, consDeactiveBounddisjunction,
         consEnableBounddisjunction, consDisableBounddisjunction,
         consPrintBounddisjunction, consCopyBounddisjunction, consParseBounddisjunction,
         conshdlrdata) );

   return SCIP_OKAY;
}


/** creates and captures a bound disjunction constraint */
SCIP_RETCODE SCIPcreateConsBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< variables of the literals in the constraint */
   SCIP_BOUNDTYPE*       boundtypes,         /**< types of bounds of the literals (lower or upper bounds) */
   SCIP_Real*            bounds,             /**< bounds of the literals */
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
                                              *   are seperated as constraints. */
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

   /* find the bounddisjunction constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("bound disjunction constraint handler not found\n");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   SCIP_CALL( consdataCreate(scip, &consdata, nvars, vars, boundtypes, bounds) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** gets number of variables in bound disjunction constraint */
int SCIPgetNVarsBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a bound disjunction constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets array of variables in bound disjunction constraint */
SCIP_VAR** SCIPgetVarsBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a bound disjunction constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets array of bound types in bound disjunction constraint */
SCIP_BOUNDTYPE* SCIPgetBoundtypesBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a bound disjunction constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->boundtypes;
}

/** gets array of bounds in bound disjunction constraint */
SCIP_Real* SCIPgetBoundsBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a bound disjunction constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->bounds;
}

