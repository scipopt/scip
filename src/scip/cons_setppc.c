/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
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
#pragma ident "@(#) $Id: cons_setppc.c,v 1.51 2004/07/07 09:52:41 bzfwolte Exp $"

/**@file   cons_setppc.c
 * @brief  constraint handler for the set partitioning / packing / covering constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_setppc.h"
#include "cons_linear.h"


#define CONSHDLR_NAME          "setppc"
#define CONSHDLR_DESC          "set partitioning / packing / covering constraints"
#define CONSHDLR_SEPAPRIORITY   +700000
#define CONSHDLR_ENFOPRIORITY   +700000
#define CONSHDLR_CHECKPRIORITY  -700000
#define CONSHDLR_SEPAFREQ             5
#define CONSHDLR_PROPFREQ             1
#define CONSHDLR_EAGERFREQ          100
#define CONSHDLR_MAXPREROUNDS        -1
#define CONSHDLR_NEEDSCONS         TRUE

#define LINCONSUPGD_PRIORITY    +700000

#define EVENTHDLR_NAME         "setppc"
#define EVENTHDLR_DESC         "bound change event handler for set partitioning / packing / covering constraints"

#define CONFLICTHDLR_NAME      "setppc"
#define CONFLICTHDLR_DESC      "conflict handler creating set covering constraints"
#define CONFLICTHDLR_PRIORITY  LINCONSUPGD_PRIORITY

#ifdef BRANCHLP
#define MINBRANCHWEIGHT             0.3  /**< minimum weight of both sets in binary set branching */
#define MAXBRANCHWEIGHT             0.7  /**< maximum weight of both sets in binary set branching */
#endif
#define DEFAULT_NPSEUDOBRANCHES       2  /**< number of children created in pseudo branching (0: disable branching) */


/** type of setppc constraint: set partitioning, set packing, or set covering */
enum SetppcType
{
   SCIP_SETPPCTYPE_PARTITIONING = 0,     /**< constraint is a set partitioning constraint: sum(x) == 1 */
   SCIP_SETPPCTYPE_PACKING      = 1,     /**< constraint is a set packing constraint:      sum(x) <= 1 */
   SCIP_SETPPCTYPE_COVERING     = 2      /**< constraint is a set covering constraint:     sum(x) >= 1 */
};
typedef enum SetppcType SETPPCTYPE;

/** constraint handler data */
struct ConshdlrData
{
   EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   INTARRAY*        varuses;            /**< number of times a var is used in the active set ppc constraints */
   int              npseudobranches;    /**< number of children created in pseudo branching (0 to disable branching) */
};

/** constraint data for set partitioning / packing / covering constraints */
struct ConsData
{
   ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   VAR**            vars;               /**< variables of the constraint */
   int              varssize;           /**< size of vars array */
   int              nvars;              /**< number of variables in the constraint */
   int              nfixedzeros;        /**< current number of variables fixed to zero in the constraint */
   int              nfixedones;         /**< current number of variables fixed to one in the constraint */
   unsigned int     setppctype:2;       /**< type of constraint: set partitioning, packing or covering */
   unsigned int     propagated:1;       /**< is constraint already preprocessed/propagated? */
};




/*
 * Local methods
 */

/** creates constaint handler data for set partitioning / packing / covering constraint handler */
static
RETCODE conshdlrdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, conshdlrdata) );
   CHECK_OKAY( SCIPcreateIntarray(scip, &(*conshdlrdata)->varuses) );
   (*conshdlrdata)->npseudobranches = DEFAULT_NPSEUDOBRANCHES;

   /* get event handler for bound change events */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      errorMessage("event handler for set partitioning / packing / covering constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** frees constraint handler data for set partitioning / packing / covering constraint handler */
static
RETCODE conshdlrdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   CHECK_OKAY( SCIPfreeIntarray(scip, &(*conshdlrdata)->varuses) );
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
   INTARRAY* varuses;

   assert(conshdlrdata != NULL);
   assert(var != NULL);

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   /* if the variable is the negation of a problem variable, count the varuses in the problem variable */
   if( SCIPvarIsNegated(var) )
   {
      VAR* negvar;

      CHECK_OKAY( SCIPgetNegatedVar(scip, var, &negvar) );
      var = negvar;
   }

   /* increase varuses counter */
   CHECK_OKAY( SCIPincIntarrayVal(scip, varuses, SCIPvarGetIndex(var), +1) );

   /*debugMessage("varuses of <%s>: %d\n", SCIPvarGetName(var), SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var)));*/

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
   INTARRAY* varuses;

   assert(conshdlrdata != NULL);
   assert(var != NULL);

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   /* if the variable is the negation of a problem variable, count the varuses in the problem variable */
   if( SCIPvarIsNegated(var) )
   {
      VAR* negvar;

      CHECK_OKAY( SCIPgetNegatedVar(scip, var, &negvar) );
      var = negvar;
   }

   /* decrease varuses counter */
   CHECK_OKAY( SCIPincIntarrayVal(scip, varuses, SCIPvarGetIndex(var), -1) );
   /*debugMessage("varuses of <%s>: %d\n", SCIPvarGetName(var), SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var)));*/
   assert(SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var)) >= 0);

   return SCIP_OKAY;
}

/** ensures, that the vars array can store at least num entries */
static
RETCODE consdataEnsureVarsSize(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);
   
   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      CHECK_OKAY( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}

/** catches events for variable at given position */
static
RETCODE consdataCatchEvent(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< set partitioning / packing / covering constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   VAR* var;

   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars != NULL);

   var = consdata->vars[pos];
   assert(var != NULL);

   /* catch bound change events on variable */
   CHECK_OKAY( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (EVENTDATA*)consdata) );
   
   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      consdata->nfixedzeros++;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      consdata->nfixedones++;

   return SCIP_OKAY;
}

/** drops events for variable at given position */
static
RETCODE consdataDropEvent(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< set partitioning / packing / covering constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   VAR* var;

   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars != NULL);

   var = consdata->vars[pos];
   assert(var != NULL);
   
   /* drop events on variable */
   CHECK_OKAY( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (EVENTDATA*)consdata) );

   /* update the fixed variables counters for this variable */
   if( SCIPisEQ(scip, SCIPvarGetUbLocal(var), 0.0) )
      consdata->nfixedzeros--;
   else if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), 1.0) )
      consdata->nfixedones--;

   return SCIP_OKAY;
}

/** catches bound change events for all variables in transformed setppc constraint */
static
RETCODE consdataCatchAllEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< setppc constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* catch event for every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      CHECK_OKAY( consdataCatchEvent(scip, consdata, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for all variables in transformed setppc constraint */
static
RETCODE consdataDropAllEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< setppc constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);

   /* drop event of every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      CHECK_OKAY( consdataDropEvent(scip, consdata, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** locks the rounding locks associated to the given variable in the setppc constraint */
static
void consdataLockRounding(
   CONSDATA*        consdata,           /**< linear constraint data */
   VAR*             var,                /**< variable of constraint entry */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   assert(consdata != NULL);

   /* forbid rounding of variable */
   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      SCIPvarLock(var, nlockspos + nlocksneg, nlockspos + nlocksneg);
      break;
   case SCIP_SETPPCTYPE_PACKING:
      SCIPvarLock(var, nlocksneg, nlockspos);
      break;
   case SCIP_SETPPCTYPE_COVERING:
      SCIPvarLock(var, nlockspos, nlocksneg);
      break;
   default:
      errorMessage("unknown setppc type\n");
      abort();
   }
}

/** unlocks the rounding locks associated to the given variable in the setppc constraint */
static
void consdataUnlockRounding(
   CONSDATA*        consdata,           /**< linear constraint data */
   VAR*             var,                /**< variable of constraint entry */
   int              nunlockspos,        /**< decrease in number of rounding locks for constraint */
   int              nunlocksneg         /**< decrease in number of rounding locks for constraint's negation */
   )
{
   assert(consdata != NULL);

   /* allow rounding of variable */
   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      SCIPvarUnlock(var, nunlockspos + nunlocksneg, nunlockspos + nunlocksneg);
      break;
   case SCIP_SETPPCTYPE_PACKING:
      SCIPvarUnlock(var, nunlocksneg, nunlockspos);
      break;
   case SCIP_SETPPCTYPE_COVERING:
      SCIPvarUnlock(var, nunlockspos, nunlocksneg);
      break;
   default:
      errorMessage("unknown setppc type\n");
      abort();
   }
}

/** locks the rounding locks of all variables in the setppc constraint */
static
void consdataLockAllRoundings(
   CONSDATA*        consdata,           /**< linear constraint data */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   int i;

   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; ++i )
      consdataLockRounding(consdata, consdata->vars[i], nlockspos, nlocksneg);
}

/** unlocks the rounding locks of all variables in the setppc constraint */
static
void consdataUnlockAllRoundings(
   CONSDATA*        consdata,           /**< linear constraint data */
   int              nunlockspos,        /**< decrease in number of rounding locks for constraint */
   int              nunlocksneg         /**< decrease in number of rounding locks for constraint's negation */
   )
{
   int i;

   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; ++i )
      consdataUnlockRounding(consdata, consdata->vars[i], nunlockspos, nunlocksneg);
}

/** creates a set partitioning / packing / covering constraint data object */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the set partitioning / packing / covering constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< variables of the constraint */
   SETPPCTYPE       setppctype          /**< type of constraint: set partitioning, packing, or covering constraint */
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
   (*consdata)->nfixedzeros = 0;
   (*consdata)->nfixedones = 0;
   (*consdata)->setppctype = setppctype; /*lint !e641*/
   (*consdata)->propagated = FALSE;

   return SCIP_OKAY;
}   

/** creates a transformed set partitioning / packing / covering constraint data object */
static
RETCODE consdataCreateTransformed(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the set partitioning / packing / covering constraint */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< variables of the constraint */
   SETPPCTYPE       setppctype          /**< type of constraint: set partitioning, packing, or covering constraint */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, consdata, nvars, vars, setppctype) );

   /* transform the variables */
   CHECK_OKAY( SCIPgetTransformedVars(scip, nvars, (*consdata)->vars, (*consdata)->vars) );

   /* catch bound change events of variables */
   CHECK_OKAY( consdataCatchAllEvents(scip, *consdata, eventhdlr) );

   return SCIP_OKAY;
}

/** frees a set partitioning / packing / covering constraint data */
static
RETCODE consdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the set partitioning / packing / covering constraint */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* if constraint belongs to transformed problem space, drop bound change events on variables */
   if( (*consdata)->nvars > 0 && SCIPvarIsTransformed((*consdata)->vars[0]) )
   {
      CHECK_OKAY( consdataDropAllEvents(scip, *consdata, eventhdlr) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints set partitioning / packing / covering constraint to file stream */
static
void consdataPrint(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< set partitioning / packing / covering constraint data */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);

   if( file == NULL )
      file = stdout;

   /* print coefficients */
   if( consdata->nvars == 0 )
      fprintf(file, "0 ");
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      fprintf(file, "+<%s> ", SCIPvarGetName(consdata->vars[v]));
   }

   /* print right hand side */
   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      fprintf(file, "== 1\n");
      break;
   case SCIP_SETPPCTYPE_PACKING:
      fprintf(file, ">= 1\n");
      break;
   case SCIP_SETPPCTYPE_COVERING:
      fprintf(file, "<= 1\n");
      break;
   default:
      errorMessage("unknown setppc type\n");
      abort();
   }
}

/** adds coefficient in setppc constraint */
static
RETCODE addCoef(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   VAR*             var                 /**< variable to add to the constraint */
   )
{
   CONSDATA* consdata;
   Bool transformed;

   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      CHECK_OKAY( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   CHECK_OKAY( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1) );
   consdata->vars[consdata->nvars] = var;
   consdata->nvars++;

   /* if we are in transformed problem, catch the variable's events */
   if( transformed )
   {
      CONSHDLR* conshdlr;
      CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* catch bound change events of variable */
      CHECK_OKAY( consdataCatchEvent(scip, consdata, conshdlrdata->eventhdlr, consdata->nvars-1) );
   }

   /* if necessary, update the rounding locks of variable */
   if( SCIPconsIsLocked(cons) )
   {
      assert(transformed);
      consdataLockRounding(consdata, var, (int)SCIPconsIsLockedPos(cons), (int)SCIPconsIsLockedNeg(cons));
   }

   consdata->propagated = FALSE;

   /* add the new coefficient to the LP row */
   if( consdata->row != NULL )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->row, var, 1.0) );
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from setppc constraint data */
static
RETCODE delCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< set partitioning / packing / covering constraint */
   int              pos                 /**< position of coefficient to delete */
   )
{
   CONSDATA* consdata;
   VAR* var;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   var = consdata->vars[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* if necessary, update the rounding locks of variable */
   if( SCIPconsIsActive(cons) && SCIPconsIsGlobal(cons) )
   {
      assert(SCIPconsIsTransformed(cons));
      consdataUnlockRounding(consdata, var, (int)SCIPconsIsLockedPos(cons), (int)SCIPconsIsLockedNeg(cons));
   }

   if( SCIPconsIsTransformed(cons) )
   {
      CONSHDLR* conshdlr;
      CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* drop bound change events of variable */
      CHECK_OKAY( consdataDropEvent(scip, consdata, conshdlrdata->eventhdlr, pos) );
   }

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->nvars--;

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}

/** deletes all zero-fixed variables */
static
RETCODE applyFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< set partitioning / packing / covering constraint */
   )
{
   CONSDATA* consdata;
   VAR* var;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nfixedzeros >= 1 )
   {
      assert(consdata->vars != NULL);

      v = 0;
      while( v < consdata->nvars )
      {
         var = consdata->vars[v];
         if( SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
         {
            CHECK_OKAY( delCoefPos(scip, cons, v) );
         }
         else
            ++v;
      }
   }

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint where all of the variables where assigned to zero,
 *  and adds conflict clause to problem
 */
static
RETCODE analyzeConflictZero(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< set partitioning / packing / covering constraint data */
   )
{
   int v;

   assert(consdata != NULL);
   assert(consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING
      || consdata->setppctype == SCIP_SETPPCTYPE_COVERING); /*lint !e641*/

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

/** analyzes conflicting assignment on given constraint where two of the variables where assigned to one,
 *  and adds conflict clause to problem
 */
static
RETCODE analyzeConflictOne(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< set partitioning / packing / covering constraint data */
   )
{
   int v;
   int n;

   assert(consdata != NULL);
   assert(consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING
      || consdata->setppctype == SCIP_SETPPCTYPE_PACKING); /*lint !e641*/

   /* initialize conflict analysis, and add the two variables assigned to one to conflict candidate queue */
   CHECK_OKAY( SCIPinitConflictAnalysis(scip) );
   n = 0;
   for( v = 0; v < consdata->nvars && n < 2; ++v )
   {
      if( SCIPvarGetLbLocal(consdata->vars[v]) > 0.5 )
      {
         CHECK_OKAY( SCIPaddConflictVar(scip, consdata->vars[v]) );
         n++;
      }
   }
   assert(n == 2);

   /* analyze the conflict */
   CHECK_OKAY( SCIPanalyzeConflict(scip, NULL) );

   return SCIP_OKAY;
}

/** checks constraint for violation only looking at the fixed variables, applies further fixings if possible */
static
RETCODE processFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< set partitioning / packing / covering constraint to be processed */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   Bool*            addcut,             /**< pointer to store whether this constraint must be added as a cut */
   Bool*            mustcheck           /**< pointer to store whether this constraint must be checked for feasibility */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(reduceddom != NULL);
   assert(addcut != NULL);
   assert(mustcheck != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(0 <= consdata->nfixedzeros && consdata->nfixedzeros <= consdata->nvars);
   assert(0 <= consdata->nfixedones && consdata->nfixedones <= consdata->nvars);

   *addcut = FALSE;
   *mustcheck = FALSE;

   if( consdata->nfixedones >= 2 )
   {
      /* at least two variables are fixed to 1:
       * - a set covering constraint is feasible anyway and can be disabled
       * - a set partitioning or packing constraint is infeasible
       */
      if( consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      }
      else
      {
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         if( SCIPconsIsGlobal(cons) )
         {
            /* use conflict analysis to get a conflict clause out of the conflicting assignment */
            CHECK_OKAY( analyzeConflictOne(scip, consdata) );
         }
         *cutoff = TRUE;
      }
   }
   else if( consdata->nfixedones == 1 )
   {
      /* exactly one variable is fixed to 1:
       * - a set covering constraint is feasible anyway and can be disabled
       * - all other variables in a set partitioning or packing constraint must be zero
       */
      if( consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      }
      else
      {
         if( consdata->nfixedzeros < consdata->nvars - 1 )
         {
            VAR** vars;
            VAR* var;
            Bool fixedonefound;
            Bool fixed;
            Bool infeasible;
            Bool tightened;
            int nvars;
            int v;

            /* unfixed variables exist: fix them to zero */
            vars = consdata->vars;
            nvars = consdata->nvars;
            fixedonefound = FALSE;
            fixed = FALSE;
            for( v = 0; v < nvars; ++v )
            {
               var = vars[v];
               assert(!fixedonefound || SCIPisZero(scip, SCIPvarGetLbLocal(var)));
               assert(SCIPisZero(scip, SCIPvarGetUbLocal(var)) || SCIPisEQ(scip, SCIPvarGetUbLocal(var), 1.0));
               if( SCIPvarGetLbLocal(var) < 0.5 )
               {
                  CHECK_OKAY( SCIPinferBinVar(scip, var, FALSE, cons, 0, &infeasible, &tightened) );
                  assert(!infeasible);
                  fixed = fixed || tightened;
               }
               else
                  fixedonefound = TRUE;
            }
            /* the fixed to one variable must have been found, and at least one variable must have been fixed */
            assert(fixedonefound && fixed);

            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            *reduceddom = TRUE;
         }

         /* now all other variables are fixed to zero:
          * the constraint is feasible, and if it's not modifiable, it is redundant
          */
         if( !SCIPconsIsModifiable(cons) )
         {
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
         }
      }
   }
   else if( consdata->nfixedzeros == consdata->nvars )
   {
      /* all variables are fixed to zero:
       * - a set packing constraint is feasible anyway, and if it's unmodifiable, it can be disabled
       * - a set partitioning or covering constraint is infeasible, and if it's unmodifiable, the node
       *   can be cut off -- otherwise, the constraint must be added as a cut and further pricing must
       *   be performed
       */
      assert(consdata->nfixedones == 0);
      
      if( consdata->setppctype == SCIP_SETPPCTYPE_PACKING ) /*lint !e641*/
      {
         if( !SCIPconsIsModifiable(cons) )
         {
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
         }
      }
      else
      {
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         if( SCIPconsIsModifiable(cons) )
            *addcut = TRUE;
         else
         {
            if( SCIPconsIsGlobal(cons) )
            {
               /* use conflict analysis to get a conflict clause out of the conflicting assignment */
               CHECK_OKAY( analyzeConflictZero(scip, consdata) );
            }
            *cutoff = TRUE;
         }
      }
   }
   else if( consdata->nfixedzeros == consdata->nvars - 1 )
   {
      /* all variables except one are fixed to zero:
       * - a set packing constraint is feasible anyway, and if it's unmodifiable, it can be disabled
       * - an unmodifiable set partitioning or covering constraint is feasible and can be disabled after the
       *   remaining variable is fixed to one
       * - a modifiable set partitioning or covering constraint must be checked manually
       */
      assert(consdata->nfixedones == 0);
      
      if( consdata->setppctype == SCIP_SETPPCTYPE_PACKING ) /*lint !e641*/
      {
         if( !SCIPconsIsModifiable(cons) )
         {
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
         }
      }
      else if( !SCIPconsIsModifiable(cons) )
      {
         VAR** vars;
         VAR* var;
         Bool infeasible;
         Bool tightened;
         int nvars;
         int v;
         
         /* search the single variable that can be fixed */
         vars = consdata->vars;
         nvars = consdata->nvars;
         for( v = 0; v < nvars; ++v )
         {
            var = vars[v];
            assert(SCIPisZero(scip, SCIPvarGetLbLocal(var)));
            assert(SCIPisZero(scip, SCIPvarGetUbLocal(var)) || SCIPisEQ(scip, SCIPvarGetUbLocal(var), 1.0));
            if( SCIPvarGetUbLocal(var) > 0.5 )
            {
               CHECK_OKAY( SCIPinferBinVar(scip, var, TRUE, cons, 0, &infeasible, &tightened) );
               assert(!infeasible);
               assert(tightened);
               break;
            }
         }
         assert(v < nvars);

         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
         *reduceddom = TRUE;
      }
      else
         *mustcheck = TRUE;
   }
   else
   {
      /* no variable is fixed to one, and at least two variables are not fixed to zero:
       * - the constraint must be checked manually
       */
      assert(consdata->nfixedones == 0);
      assert(consdata->nfixedzeros < consdata->nvars - 1);

      *mustcheck = TRUE;
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, returns TRUE iff constraint is feasible */
static
Bool checkCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< set partitioning / packing / covering constraint to be checked */
   SOL*             sol                 /**< primal CIP solution */
   )
{
   VAR** vars;
   Real solval;
   Real sum;
   int nvars;
   int v;
   
   /* calculate the constraint's activity */
   vars = consdata->vars;
   nvars = consdata->nvars;
   sum = 0.0;
   assert(SCIPfeastol(scip) < 0.1); /* to make the comparison against 1.1 working */
   for( v = 0; v < nvars && sum < 1.1; ++v )  /* if sum >= 1.1, the feasibility is clearly decided */
   {
      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY);
      solval = SCIPgetSolVal(scip, sol, vars[v]);
      assert(SCIPisFeasGE(scip, solval, 0.0) && SCIPisFeasLE(scip, solval, 1.0));
      sum += solval;
   }

   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      return SCIPisFeasEQ(scip, sum, 1.0);
   case SCIP_SETPPCTYPE_PACKING:
      return SCIPisFeasLE(scip, sum, 1.0);
   case SCIP_SETPPCTYPE_COVERING:
      return SCIPisFeasGE(scip, sum, 1.0);
   default:
      errorMessage("unknown setppc type\n");
      abort();
   }
}

/** creates an LP row in a set partitioning / packing / covering constraint data object */
static
RETCODE createRow(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< set partitioning / packing / covering constraint */
   )
{
   CONSDATA* consdata;
   Real lhs;
   Real rhs;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   switch( consdata->setppctype )
   {
   case SCIP_SETPPCTYPE_PARTITIONING:
      lhs = 1.0;
      rhs = 1.0;
      break;
   case SCIP_SETPPCTYPE_PACKING:
      lhs = -SCIPinfinity(scip);
      rhs = 1.0;
      break;
   case SCIP_SETPPCTYPE_COVERING:
      lhs = 1.0;
      rhs = SCIPinfinity(scip);
      break;
   default:
      errorMessage("unknown setppc type\n");
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), lhs, rhs,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );

   CHECK_OKAY( SCIPaddVarsToRowSameCoef(scip, consdata->row, consdata->nvars, consdata->vars, 1.0) );

   return SCIP_OKAY;
}

/** adds setppc constraint as cut to the LP */
static
RETCODE addCut(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< setppc constraint */
   Real             violation           /**< absolute violation of the constraint */
   )
{
   CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   if( consdata->row == NULL )
   {
      /* convert set partitioning constraint data into LP row */
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
   CONS*            cons,               /**< set partitioning / packing / covering constraint to be separated */
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
   assert(0 <= consdata->nfixedzeros && consdata->nfixedzeros <= consdata->nvars);
   assert(0 <= consdata->nfixedones && consdata->nfixedones <= consdata->nvars);

   /* skip constraints already in the LP */
   if( consdata->row != NULL && SCIProwIsInLP(consdata->row) )
      return SCIP_OKAY;

   /* check constraint for violation only looking at the fixed variables, apply further fixings if possible */
   CHECK_OKAY( processFixings(scip, cons, cutoff, reduceddom, &addcut, &mustcheck) );

   if( mustcheck )
   {
      assert(!addcut);

      /* variable's fixings didn't give us any information -> we have to check the constraint */
      if( consdata->row != NULL )
      {
         assert(!SCIProwIsInLP(consdata->row));
         addcut = !SCIPisFeasible(scip, SCIPgetRowLPFeasibility(scip, consdata->row));
      }
      else
         addcut = !checkCons(scip, consdata, NULL);

      if( !addcut )
      {
         /* constraint was feasible -> increase age */
         CHECK_OKAY( SCIPincConsAge(scip, cons) );
      }
   }

   if( addcut )
   {
      /* insert LP row as cut */
      CHECK_OKAY( addCut(scip, cons, 1.0) );
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      *separated = TRUE;
   }

   return SCIP_OKAY;
}

/** enforces the pseudo solution on the given constraint */
static
RETCODE enforcePseudo(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< set partitioning / packing / covering constraint to be separated */
   Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   Bool*            infeasible,         /**< pointer to store TRUE, if the constraint was infeasible */
   Bool*            reduceddom,         /**< pointer to store TRUE, if a domain reduction was found */
   Bool*            solvelp             /**< pointer to store TRUE, if the LP has to be solved */
   )
{
   Bool addcut;
   Bool mustcheck;

   assert(!SCIPhasActNodeLP(scip));
   assert(cons != NULL);
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(reduceddom != NULL);
   assert(solvelp != NULL);

   /* check constraint for violation only looking at the fixed variables, apply further fixings if possible */
   CHECK_OKAY( processFixings(scip, cons, cutoff, reduceddom, &addcut, &mustcheck) );

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

   if( addcut )
   {
      /* a cut must be added to the LP -> we have to solve the LP immediately */
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
      *solvelp = TRUE;
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
DECL_CONSFREE(consFreeSetppc)
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


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitSetppc NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitSetppc NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreSetppc NULL

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreSetppc NULL
/**@TODO enable implication detection */
#if 0
static
DECL_CONSEXITPRE(consExitpreSetppc)
{
   CONSDATA* consdata;
   int c;
   int i;
   int j;

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
#ifdef DEBUG
      debugMessage("\n setppc constraint <%s>:\n", SCIPconsGetName(conss[c]));
      for( i = 0; i < consdata->nvars; i++ )
         debugMessage("1<%s> + ", SCIPvarGetName(consdata->vars[i]));
      if( consdata->setppctype == 0 ) 
         debugMessage(" == 1\n");
      else if( consdata->setppctype == 1 )
         debugMessage(" <= 1\n");
      else
         debugMessage(" >= 1\n");
#endif 

      /* constraint is a set partitioning constraint: sum(x) == 1 or
       * constraint is a set packing constraint:      sum(x) <= 1 */
      if( consdata->setppctype == 0 || consdata->setppctype == 1 )
      {
         for( i = 0; i < consdata->nvars - 1; i++ )
         {
            for( j = i+1; j < consdata->nvars; j++ )
            {
               /* add lower bound implication (x_i >= 1  ==> x_j <= 0) to x_i */
               CHECK_OKAY( SCIPaddVarLbimpl(scip, consdata->vars[i], 1, consdata->vars[j], TRUE, 0) );
               
               /* add lower bound implication (x_j >= 1  ==> x_i <= 0) to x_j */
               CHECK_OKAY( SCIPaddVarLbimpl(scip, consdata->vars[j], 1, consdata->vars[i], TRUE, 0) );
            }
         }
      }
      /* constraint is a set covering constraint:     sum(x) >= 1 and has 2 variables x_j */
      else if( consdata->nvars == 2 )
      {
         assert(consdata->setppctype == 2);
         /* add upper bound implication (x_0 <= 0  ==> x_1 >= 1) to x_0 */
         CHECK_OKAY( SCIPaddVarUbimpl(scip, consdata->vars[0], 0, consdata->vars[1], FALSE, 1) );
         
         /* add upper bound implication (x_1 <= 0  ==> x_0 >= 1) to x_1 */
         CHECK_OKAY( SCIPaddVarUbimpl(scip, consdata->vars[1], 0, consdata->vars[0], FALSE, 1) );
      }
   }

   return SCIP_OKAY;
}
#endif

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolSetppc NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
DECL_CONSEXITSOL(consExitsolSetppc)
{
   CONSDATA* consdata;
   int c;

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         CHECK_OKAY( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteSetppc)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);
   
   /* free setppc constraint data */
   CHECK_OKAY( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransSetppc)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   /*debugMessage("Trans method of setppc constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPstage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);
   
   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */

   /* create constraint data for target constraint */
   CHECK_OKAY( consdataCreateTransformed(scip, &targetdata, conshdlrdata->eventhdlr,
         sourcedata->nvars, sourcedata->vars, (SETPPCTYPE)sourcedata->setppctype) );

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpSetppc)
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
DECL_CONSSEPA(consSepaSetppc)
{  /*lint --e{715}*/
   Bool cutoff;
   Bool separated;
   Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("separating %d/%d set partitioning / packing / covering constraints\n", nusefulconss, nconss);

   *result = SCIP_DIDNOTFIND;

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* check all useful set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], &cutoff, &separated, &reduceddom) );
   }

   /* combine set partitioning / packing / covering constraints to get more cuts */
   /**@todo further cuts of set partitioning / packing / covering constraints */

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
/** if fractional variables exist, chooses a set S of them and branches on (i) x(S) == 0, and (ii) x(S) >= 1 */
static
RETCODE branchLP(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< set partitioning / packing / covering constraint handler */
   RESULT*          result              /**< pointer to store the result SCIP_BRANCHED, if branching was applied */
   )
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* varuses;
   VAR** lpcands;
   VAR** sortcands;
   VAR* var;
   Real branchweight;
   Real solval;
   int* uses;
   int nlpcands;
   int nsortcands;
   int nselcands;
   int numuses;
   int i;
   int j;

   /**@todo use a better set partitioning / packing / covering branching on LP solution (use SOS branching) */

   assert(conshdlr != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   /* get fractional variables */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, NULL, NULL, &nlpcands) );
   if( nlpcands == 0 )
      return SCIP_OKAY;

   /* get temporary memory */
   CHECK_OKAY( SCIPallocBufferArray(scip, &sortcands, nlpcands) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &uses, nlpcands) );
   
   /* sort fractional variables by number of uses in enabled set partitioning / packing / covering constraints */
   nsortcands = 0;
   for( i = 0; i < nlpcands; ++i )
   {
      var = lpcands[i];
      numuses = SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var));
      if( numuses > 0 )
      {
         for( j = nsortcands; j > 0 && numuses > uses[j-1]; --j )
         {
            sortcands[j] = sortcands[j-1];
            uses[j] = uses[j-1];
         }
         assert(0 <= j && j <= nsortcands);
         sortcands[j] = var;
         uses[j] = numuses;
         nsortcands++;
      }
   }
   assert(nsortcands <= nlpcands);

   /* if none of the fractional variables is member of a set partitioning / packing / covering constraint,
    * we are not responsible for doing the branching
    */
   if( nsortcands > 0 )
   {
      /* select the first variables from the sorted candidate list, until MAXBRANCHWEIGHT is reached;
       * then choose one less
       */
      branchweight = 0.0;
      solval = 0.0;
      for( nselcands = 0; nselcands < nsortcands && branchweight <= MAXBRANCHWEIGHT; ++nselcands )
      {
         solval = SCIPgetVarSol(scip, sortcands[nselcands]);
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
            CHECK_OKAY( SCIPchgVarUbNode(scip, node, sortcands[i], 0.0) );
         }

         /* create right child: add constraint x(S) >= 1 */
         CHECK_OKAY( SCIPcreateChild(scip, &node) );
         if( nselcands == 1 )
         {
            /* only one candidate selected: fix it to 1.0 */
            debugMessage("fixing variable <%s> to 1.0 in right child node\n", SCIPvarGetName(sortcands[0]));
            CHECK_OKAY( SCIPchgVarLbNode(scip, node, sortcands[0], 1.0) );
         }
         else
         {
            CONS* newcons;
            char name[MAXSTRLEN];
         
            /* add set covering constraint x(S) >= 1 */
            sprintf(name, "BSB%lld", SCIPgetNTotalNodes(scip));

            CHECK_OKAY( SCIPcreateConsSetcover(scip, &newcons, name, nselcands, sortcands,
                  FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE) );
            CHECK_OKAY( SCIPaddConsNode(scip, node, newcons) );
            CHECK_OKAY( SCIPreleaseCons(scip, &newcons) );
         }
      
         *result = SCIP_BRANCHED;
         
#ifdef DEBUG
         debugMessage("binary set branching: nselcands=%d/%d, weight(S)=%g, A={", nselcands, nlpcands, branchweight);
         for( i = 0; i < nselcands; ++i )
            printf(" %s[%g]", SCIPvarGetName(sortcands[i]), SCIPgetSolVal(scip, NULL, sortcands[i]));
         printf(" }\n");
#endif
      }
   }

   /* free temporary memory */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &uses) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &sortcands) );

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
   CONSHDLR*        conshdlr,           /**< set partitioning / packing / covering constraint handler */
   RESULT*          result              /**< pointer to store the result SCIP_BRANCHED, if branching was applied */
   )
{
   CONSHDLRDATA* conshdlrdata;
   INTARRAY* varuses;
   VAR** pseudocands;
   VAR** branchcands;
   VAR* var;
   NODE* node;
   int* canduses;
   int npseudocands;
   int maxnbranchcands;
   int nbranchcands;
   int uses;
   int i;
   int j;

   /**@todo use a better set partitioning / packing / covering branching on pseudo solution (use SOS branching) */

   assert(conshdlr != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check, if pseudo branching is disabled */
   if( conshdlrdata->npseudobranches <= 1 )
      return SCIP_OKAY;

   /* get fractional variables */
   CHECK_OKAY( SCIPgetPseudoBranchCands(scip, &pseudocands, NULL, &npseudocands) );
   if( npseudocands == 0 )
      return SCIP_OKAY;

   varuses = conshdlrdata->varuses;
   assert(varuses != NULL);

   /* choose the maximal number of branching variables */
   maxnbranchcands = conshdlrdata->npseudobranches-1;
   assert(maxnbranchcands >= 1);

   /* get temporary memory */
   CHECK_OKAY( SCIPallocBufferArray(scip, &branchcands, maxnbranchcands) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &canduses, maxnbranchcands) );
   
   /* sort unfixed variables by number of uses in enabled set partitioning / packing / covering constraints */
   nbranchcands = 0;
   for( i = 0; i < npseudocands; ++i )
   {
      var = pseudocands[i];
      uses = SCIPgetIntarrayVal(scip, varuses, SCIPvarGetIndex(var));
      if( uses > 0 )
      {
         if( nbranchcands < maxnbranchcands || uses > canduses[nbranchcands-1] )
         {
            for( j = MIN(nbranchcands, maxnbranchcands-1); j > 0 && uses > canduses[j-1]; --j )
            {
               branchcands[j] = branchcands[j-1];
               canduses[j] = canduses[j-1];
            }
            assert(0 <= j && j <= nbranchcands && j < maxnbranchcands);
            branchcands[j] = var;
            canduses[j] = uses;
            if( nbranchcands < maxnbranchcands )
               nbranchcands++;
         }
      }
   }
   assert(nbranchcands <= maxnbranchcands);

   /* if none of the unfixed variables is member of a set partitioning / packing / covering constraint,
    * we are not responsible for doing the branching
    */
   if( nbranchcands > 0 )
   {
      /* branch on the first part of the sorted candidates:
       * - for each of these variables i, create a child node x_0 = ... = x_i-1 = 0, x_i = 1
       * - create an additional child node x_0 = ... = x_n-1 = 0
       */
      for( i = 0; i < nbranchcands; ++i )
      {            
         /* create child with x_0 = ... = x_i-1 = 0, x_i = 1 */
         CHECK_OKAY( SCIPcreateChild(scip, &node, (Real)nbranchcands) );
         for( j = 0; j < i; ++j )
         {
            CHECK_OKAY( SCIPchgVarUbNode(scip, node, branchcands[j], 0.0) );
         }
         CHECK_OKAY( SCIPchgVarLbNode(scip, node, branchcands[i], 1.0) );
      }
      /* create child with x_0 = ... = x_n = 0 */
      CHECK_OKAY( SCIPcreateChild(scip, &node, (Real)i) );
      for( i = 0; i < nbranchcands; ++i )
      {
         CHECK_OKAY( SCIPchgVarUbNode(scip, node, branchcands[i], 0.0) );
      }

      *result = SCIP_BRANCHED;

#ifdef DEBUG
      {
         int nchildren;
         CHECK_OKAY( SCIPgetChildren(scip, NULL, &nchildren) );
         debugMessage("branched on pseudo solution: %d children\n", nchildren);
      }
#endif
   }

   /* free temporary memory */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &canduses) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &branchcands) );

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpSetppc)
{  /*lint --e{715}*/
   Bool cutoff;
   Bool separated;
   Bool reduceddom;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   debugMessage("LP enforcing %d set partitioning / packing / covering constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   cutoff = FALSE;
   separated = FALSE;
   reduceddom = FALSE;

   /* check all useful set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nusefulconss && !cutoff && !reduceddom; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], &cutoff, &separated, &reduceddom) );
   }

   /* check all obsolete set partitioning / packing / covering constraints for feasibility */
   for( c = nusefulconss; c < nconss && !cutoff && !separated && !reduceddom; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], &cutoff, &separated, &reduceddom) );
   }

#ifdef BRANCHLP
   if( !cutoff && !separated && !reduceddom )
   {
      /* if solution is not integral, choose a variable set to branch on */
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
DECL_CONSENFOPS(consEnfopsSetppc)
{  /*lint --e{715}*/
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

   debugMessage("pseudo enforcing %d set partitioning / packing / covering constraints\n", nconss);

   *result = SCIP_FEASIBLE;

   cutoff = FALSE;
   infeasible = FALSE;
   reduceddom = FALSE;
   solvelp = FALSE;

   /* check all set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nconss && !cutoff && !reduceddom && !solvelp; ++c )
   {
      CHECK_OKAY( enforcePseudo(scip, conss[c], &cutoff, &infeasible, &reduceddom, &solvelp) );
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
DECL_CONSCHECK(consCheckSetppc)
{  /*lint --e{715}*/
   CONS* cons;
   CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* check all set partitioning / packing / covering constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
      {
         if( !checkCons(scip, consdata, sol) )
         {
            /* constraint is violated */
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
            *result = SCIP_INFEASIBLE;
            return SCIP_OKAY;
         }
         else
         {
            CHECK_OKAY( SCIPincConsAge(scip, cons) );
         }
      }
   }
   
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
DECL_CONSPROP(consPropSetppc)
{  /*lint --e{715}*/
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

   debugMessage("propagating %d/%d set partitioning / packing / covering constraints\n", nusefulconss, nconss);

   cutoff = FALSE;
   reduceddom = FALSE;

   /* propagate all useful set partitioning / packing / covering constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      CHECK_OKAY( processFixings(scip, conss[c], &cutoff, &reduceddom, &addcut, &mustcheck) );
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
DECL_CONSPRESOL(consPresolSetppc)
{  /*lint --e{715}*/
   CONS* cons;
   CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* process constraints */
   for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->propagated = FALSE;

      if( consdata->propagated )
         continue;

      debugMessage("presolving set partitioning / packing / covering constraint <%s>\n", SCIPconsGetName(cons));

      /* remove all variables that are fixed to zero */
      CHECK_OKAY( applyFixings(scip, cons) );

      if( consdata->nfixedones >= 2 )
      {
         /* at least two variables are fixed to 1:
          * - a set covering constraint is feasible anyway and can be deleted
          * - a set partitioning or packing constraint is infeasible
          */
         if( consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
         {
            debugMessage("set covering constraint <%s> is redundant\n", SCIPconsGetName(cons));
            CHECK_OKAY( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }
         else
         {
            debugMessage("set partitioning / packing constraint <%s> is infeasible\n", SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
      }
      else if( consdata->nfixedones == 1 )
      {
         /* exactly one variable is fixed to 1:
          * - a set covering constraint is feasible anyway and can be disabled
          * - all other variables in a set partitioning or packing constraint must be zero
          */
         if( consdata->setppctype == SCIP_SETPPCTYPE_COVERING ) /*lint !e641*/
         {
            debugMessage("set covering constraint <%s> is redundant\n", SCIPconsGetName(cons));
            CHECK_OKAY( SCIPdelCons(scip, cons) );
            (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }
         else
         {
            VAR* var;
            int v;
            
            debugMessage("set partitioning / packing constraint <%s> has a variable fixed to 1.0\n", SCIPconsGetName(cons));
            for( v = 0; v < consdata->nvars; ++v )
            {
               var = consdata->vars[v];
               if( SCIPisZero(scip, SCIPvarGetLbGlobal(var)) && !SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
               {
                  Bool infeasible;
                  Bool fixed;

                  CHECK_OKAY( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );
                  if( infeasible )
                  {
                     debugMessage("setppc constraint <%s>: infeasible fixing <%s> == 0\n",
                        SCIPconsGetName(cons), SCIPvarGetName(var));
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
                  assert(fixed);
                  (*nfixedvars)++;
                  *result = SCIP_SUCCESS;
               }
            }

            /* now all other variables are fixed to zero:
             * the constraint is feasible, and if it's not modifiable, it is redundant
             */
            if( !SCIPconsIsModifiable(cons) )
            {
               debugMessage("set partitioning / packing constraint <%s> is redundant\n", SCIPconsGetName(cons));
               CHECK_OKAY( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
         }
      }
      else if( !SCIPconsIsModifiable(cons) )
      {
         /* all other preprocessings can only be done on non-modifiable constraints */
         if( consdata->nfixedzeros == consdata->nvars )
         {
            /* all variables are fixed to zero:
             * - a set packing constraint is feasible anyway and can be deleted
             * - a set partitioning or covering constraint is infeasible, and so is the whole problem
             */
            assert(consdata->nfixedones == 0);
            
            if( consdata->setppctype == SCIP_SETPPCTYPE_PACKING ) /*lint !e641*/
            {
               debugMessage("set packing constraint <%s> is redundant\n", SCIPconsGetName(cons));
               CHECK_OKAY( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
            else
            {
               debugMessage("set partitioning / covering constraint <%s> is infeasible\n", SCIPconsGetName(cons));
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
         }
         else if( consdata->nfixedzeros == consdata->nvars - 1 )
         {
            /* all variables except one are fixed to zero:
             * - a set packing constraint is feasible anyway, and can be deleted
             * - a set partitioning or covering constraint is feasible and can be deleted after the
             *   remaining variable is fixed to one
             */
            assert(consdata->nfixedones == 0);
         
            if( consdata->setppctype == SCIP_SETPPCTYPE_PACKING ) /*lint !e641*/
            {
               debugMessage("set packing constraint <%s> is redundant\n", SCIPconsGetName(cons));
               CHECK_OKAY( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
            else
            {
               VAR* var;
               Bool infeasible;
               Bool fixed;
               Bool found;
               int v;
               
               debugMessage("set partitioning / covering constraint <%s> has only one variable not fixed to 0.0\n",
                  SCIPconsGetName(cons));
               
               /* search unfixed variable */
               found = FALSE;
               var = NULL;
               for( v = 0; v < consdata->nvars && !found; ++v )
               {
                  var = consdata->vars[v];
                  found = !SCIPisZero(scip, SCIPvarGetUbGlobal(var));
               }
               assert(found);

               CHECK_OKAY( SCIPfixVar(scip, var, 1.0, &infeasible, &fixed) );
               if( infeasible )
               {
                  debugMessage("setppc constraint <%s>: infeasible fixing <%s> == 1\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var));
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               assert(fixed);
               (*nfixedvars)++;

               CHECK_OKAY( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
               *result = SCIP_SUCCESS;
               continue;
            }
         }
         else if( consdata->nfixedzeros == consdata->nvars - 2
            && consdata->setppctype == SCIP_SETPPCTYPE_PARTITIONING ) /*lint !e641*/
         {
            VAR* var;
            VAR* var1;
            VAR* var2;
            Bool infeasible;
            Bool redundant;
            Bool aggregated;
            int v;

            /* aggregate variable and delete constraint, if set partitioning constraint consists only of two
             * non-fixed variables
             */
            
            /* search unfixed variable */
            var1 = NULL;
            var2 = NULL;
            for( v = 0; v < consdata->nvars && var2 == NULL; ++v )
            {
               var = consdata->vars[v];
               if( !SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
               {
                  if( var1 == NULL )
                     var1 = var;
                  else
                     var2 = var;
               }
            }
            assert(var1 != NULL && var2 != NULL);

            /* aggregate binary equality var1 + var2 == 1 */
            debugMessage("set partitioning constraint <%s>: aggregate <%s> + <%s> == 1\n",
               SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
            CHECK_OKAY( SCIPaggregateVars(scip, var1, var2, 1.0, 1.0, 1.0, &infeasible, &redundant, &aggregated) );
            if( infeasible )
            {
               debugMessage("set partitioning constraint <%s>: infeasible aggregation <%s> + <%s> == 1\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            if( aggregated )
               (*naggrvars)++;
            if( redundant )
            {
               CHECK_OKAY( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
            }
            *result = SCIP_SUCCESS;
            continue;
         }
      }

      consdata->propagated = TRUE;
   }
   
   return SCIP_OKAY;
}


/** conflict variable resolving method of constraint handler */
static
DECL_CONSRESCVAR(consRescvarSetppc)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(infervar != NULL);
   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("conflict resolving method of set partitioning / packing / covering constraint handler\n");

   if( SCIPvarGetLbLocal(infervar) > 0.5 )
   {
      Bool confvarfound;

      /* the inference variable was infered to 1.0:
       * the inference constraint has to be a set partitioning or covering constraint, and the reason for
       * the deduction is the assignment to zero of all other variables
       */
      confvarfound = FALSE;
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( consdata->vars[v] != infervar )
         {
            assert(SCIPvarGetUbLocal(consdata->vars[v]) < 0.5); /* the reason variable must be assigned to zero */
            CHECK_OKAY( SCIPaddConflictVar(scip, consdata->vars[v]) );
         }
         else
         {
            assert(!confvarfound);
            confvarfound = TRUE;
         }
      }
      assert(confvarfound);
   }
   else
   {
      Bool reasonvarfound;

      /* the inference variable was infered to 0.0:
       * the inference constraint has to be a set partitioning or packing constraint, and the reason for
       * the deduction is the assignment to 1.0 of a single different variable
       */
      assert(SCIPvarGetUbLocal(infervar) < 0.5);
      reasonvarfound = FALSE;
      for( v = 0; v < consdata->nvars && !reasonvarfound; ++v )
      {
         if( SCIPvarGetLbLocal(consdata->vars[v]) > 0.5 )
         {
            CHECK_OKAY( SCIPaddConflictVar(scip, consdata->vars[v]) );
            reasonvarfound = TRUE;
         }
      }
      assert(reasonvarfound);
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockSetppc)
{  /*lint --e{715}*/
   consdataLockAllRoundings(SCIPconsGetData(cons), nlockspos, nlocksneg);

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockSetppc)
{  /*lint --e{715}*/
   consdataUnlockAllRoundings(SCIPconsGetData(cons), nunlockspos, nunlocksneg);

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
DECL_CONSACTIVE(consActiveSetppc)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPconsIsTransformed(cons));

   debugMessage("activation information for set partitioning / packing / covering constraint <%s>\n",
      SCIPconsGetName(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* increase the number of uses for each variable in the constraint */
   for( v = 0; v < consdata->nvars; ++v )
   {
      CHECK_OKAY( conshdlrdataIncVaruses(scip, conshdlrdata, consdata->vars[v]) );
   }

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
DECL_CONSDEACTIVE(consDeactiveSetppc)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPconsIsTransformed(cons));

   debugMessage("deactivation information for set partitioning / packing / covering constraint <%s>\n",
      SCIPconsGetName(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* decrease the number of uses for each variable in the constraint */
   for( v = 0; v < consdata->nvars; ++v )
   {
      CHECK_OKAY( conshdlrdataDecVaruses(scip, conshdlrdata, consdata->vars[v]) );
   }

   return SCIP_OKAY;
}


/** constraint enabling notification method of constraint handler */
#define consEnableSetppc NULL


/** constraint disabling notification method of constraint handler */
#define consDisableSetppc NULL

/** constraint display method of constraint handler */
static
DECL_CONSPRINT(consPrintSetppc)
{
   consdataPrint(scip, SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}




/** creates and captures a set partitioning / packing / covering constraint */
static
RETCODE createConsSetppc(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   SETPPCTYPE       setppctype,         /**< type of constraint: set partitioning, packing, or covering constraint */
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

   /* find the set partitioning constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("set partitioning / packing / covering constraint handler not found\n");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   if( SCIPstage(scip) == SCIP_STAGE_PROBLEM )
   {
      /* create constraint in original problem */
      CHECK_OKAY( consdataCreate(scip, &consdata, nvars, vars, setppctype) );
   }
   else
   {
      CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* create constraint in transformed problem */
      CHECK_OKAY( consdataCreateTransformed(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, setppctype) );
   }

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, removeable) );

   return SCIP_OKAY;
}

/** creates and captures a normalized (with all coefficients +1) setppc constraint */
static
RETCODE createNormalizedSetppc(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients (+1.0 or -1.0) */
   int              mult,               /**< multiplier on the coefficients(+1 or -1) */
   SETPPCTYPE       setppctype,         /**< type of constraint: set partitioning, packing, or covering constraint */
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
   CHECK_OKAY( SCIPallocBufferArray(scip, &transvars, nvars) );

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
   CHECK_OKAY( createConsSetppc(scip, cons, name, nvars, transvars, setppctype,
         initial, separate, enforce, check, propagate, local, modifiable, removeable) );

   /* release temporary memory */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &transvars) );

   return SCIP_OKAY;
}

static
DECL_LINCONSUPGD(linconsUpgdSetppc)
{  /*lint --e{715}*/
   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to set partitioning, packing, or covering constraint
    * - all set partitioning / packing / covering constraints consist only of binary variables with a
    *   coefficient of +1.0 or -1.0 (variables with -1.0 coefficients can be negated):
    *        lhs     <= x1 + ... + xp - y1 - ... - yn <= rhs
    * - negating all variables y = (1-Y) with negative coefficients gives:
    *        lhs + n <= x1 + ... + xp + Y1 + ... + Yn <= rhs + n
    * - negating all variables x = (1-X) with positive coefficients and multiplying with -1 gives:
    *        p - rhs <= X1 + ... + Xp + y1 + ... + yn <= p - lhs
    * - a set partitioning constraint has left hand side of +1.0, and right hand side of +1.0 : x(S) == 1.0
    *    -> without negations:  lhs == rhs == 1 - n  or  lhs == rhs == p - 1
    * - a set packing constraint has left hand side of -infinity, and right hand side of +1.0 : x(S) <= 1.0
    *    -> without negations:  (lhs == -inf  and  rhs == 1 - n)  or  (lhs == p - 1  and  rhs = +inf)
    * - a set covering constraint has left hand side of +1.0, and right hand side of +infinity: x(S) >= 1.0
    *    -> without negations:  (lhs == 1 - n  and  rhs == +inf)  or  (lhs == -inf  and  rhs = p - 1)
    */
   if( nposbin + nnegbin == nvars && ncoeffspone + ncoeffsnone == nvars )
   {
      int mult;

      if( SCIPisEQ(scip, lhs, rhs) && (SCIPisEQ(scip, lhs, 1.0 - ncoeffsnone) || SCIPisEQ(scip, lhs, ncoeffspone - 1.0)) )
      {
         debugMessage("upgrading constraint <%s> to set partitioning constraint\n", SCIPconsGetName(cons));
         
         /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
         mult = SCIPisEQ(scip, lhs, 1.0 - ncoeffsnone) ? +1 : -1;

         /* create the set partitioning constraint (an automatically upgraded constraint is always unmodifiable) */
         assert(!SCIPconsIsModifiable(cons));
         CHECK_OKAY( createNormalizedSetppc(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, mult,
               SCIP_SETPPCTYPE_PARTITIONING,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
      }
      else if( (SCIPisInfinity(scip, -lhs) && SCIPisEQ(scip, rhs, 1.0 - ncoeffsnone))
         || (SCIPisEQ(scip, lhs, ncoeffspone - 1.0) && SCIPisInfinity(scip, rhs)) )
      {
         debugMessage("upgrading constraint <%s> to set packing constraint\n", SCIPconsGetName(cons));
         
         /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
         mult = SCIPisInfinity(scip, -lhs) ? +1 : -1;

         /* create the set packing constraint (an automatically upgraded constraint is always unmodifiable) */
         assert(!SCIPconsIsModifiable(cons));
         CHECK_OKAY( createNormalizedSetppc(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, mult,
               SCIP_SETPPCTYPE_PACKING,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
      }
      else if( (SCIPisEQ(scip, lhs, 1.0 - ncoeffsnone) && SCIPisInfinity(scip, rhs))
         || (SCIPisInfinity(scip, -lhs) && SCIPisEQ(scip, rhs, ncoeffspone - 1.0)) )
      {
         debugMessage("upgrading constraint <%s> to set covering constraint\n", SCIPconsGetName(cons));
         
         /* check, if we have to multiply with -1 (negate the positive vars) or with +1 (negate the negative vars) */
         mult = SCIPisInfinity(scip, rhs) ? +1 : -1;

         /* create the set covering constraint (an automatically upgraded constraint is always unmodifiable) */
         assert(!SCIPconsIsModifiable(cons));
         CHECK_OKAY( createNormalizedSetppc(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, mult,
               SCIP_SETPPCTYPE_COVERING,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
      }
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecSetppc)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   debugMessage("Exec method of bound change event handler for set partitioning / packing / covering constraints\n");

   consdata = (CONSDATA*)eventdata;
   assert(consdata != NULL);

   eventtype = SCIPeventGetType(event);
   switch( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      consdata->nfixedones++;
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      consdata->nfixedones--;
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      consdata->nfixedzeros++;
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      consdata->nfixedzeros--;
      break;
   default:
      errorMessage("invalid event type\n");
      abort();
   }
   assert(0 <= consdata->nfixedzeros && consdata->nfixedzeros <= consdata->nvars);
   assert(0 <= consdata->nfixedones && consdata->nfixedones <= consdata->nvars);

   consdata->propagated = FALSE;

   debugMessage(" -> constraint has %d zero-fixed and %d one-fixed of %d variables\n", 
      consdata->nfixedzeros, consdata->nfixedones, consdata->nvars);

   return SCIP_OKAY;
}




/*
 * Callback methods of conflict handler
 */

static
DECL_CONFLICTEXEC(conflictExecSetppc)
{  /*lint --e{715}*/
   CONS* cons;
   char consname[MAXSTRLEN];
   
   assert(conflicthdlr != NULL);
   assert(strcmp(SCIPconflicthdlrGetName(conflicthdlr), CONFLICTHDLR_NAME) == 0);
   assert(conflictvars != NULL || nconflictvars == 0);
   assert(result != NULL);

   /* don't process already resolved conflicts */
   if( resolved )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* create a constraint out of the conflict set */
   sprintf(consname, "cf%d", SCIPgetNGlobalConss(scip));
   CHECK_OKAY( SCIPcreateConsSetcover(scip, &cons, consname, nconflictvars, conflictvars, 
         FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE) );
   CHECK_OKAY( SCIPaddConsNode(scip, node, cons) );
   CHECK_OKAY( SCIPreleaseCons(scip, &cons) );

   *result = SCIP_CONSADDED;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for set partitioning / packing / covering constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrSetppc(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL,
         NULL, eventExecSetppc,
         NULL) );

   /* create conflict handler for setppc constraints */
   CHECK_OKAY( SCIPincludeConflicthdlr(scip, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
         NULL, NULL, NULL, conflictExecSetppc,
         NULL) );

   /* create constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, 
         CONSHDLR_MAXPREROUNDS, CONSHDLR_NEEDSCONS,
         consFreeSetppc, consInitSetppc, consExitSetppc, 
         consInitpreSetppc, consExitpreSetppc, consInitsolSetppc, consExitsolSetppc,
         consDeleteSetppc, consTransSetppc, 
         consInitlpSetppc, consSepaSetppc, consEnfolpSetppc, consEnfopsSetppc, consCheckSetppc, 
         consPropSetppc, consPresolSetppc, consRescvarSetppc,
         consLockSetppc, consUnlockSetppc,
         consActiveSetppc, consDeactiveSetppc, 
         consEnableSetppc, consDisableSetppc,
         consPrintSetppc,
         conshdlrdata) );

   /* include the linear constraint to set partitioning constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdSetppc, LINCONSUPGD_PRIORITY) );

   /* set partitioning constraint handler parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
         "constraints/setppc/npseudobranches", 
         "number of children created in pseudo branching (0: disable pseudo branching)",
         &conshdlrdata->npseudobranches, DEFAULT_NPSEUDOBRANCHES, 0, INT_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}


/** creates and captures a set partitioning constraint */
RETCODE SCIPcreateConsSetpart(
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
   return createConsSetppc(scip, cons, name, nvars, vars, SCIP_SETPPCTYPE_PARTITIONING,
      initial, separate, enforce, check, propagate, local, modifiable, removeable);
}

/** creates and captures a set packing constraint */
RETCODE SCIPcreateConsSetpack(
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
   return createConsSetppc(scip, cons, name, nvars, vars, SCIP_SETPPCTYPE_PACKING,
      initial, separate, enforce, check, propagate, local, modifiable, removeable);
}

/** creates and captures a set covering constraint */
RETCODE SCIPcreateConsSetcover(
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
   return createConsSetppc(scip, cons, name, nvars, vars, SCIP_SETPPCTYPE_COVERING,
      initial, separate, enforce, check, propagate, local, modifiable, removeable);
}

/** adds coefficient in set partitioning / packing / covering constraint */
RETCODE SCIPaddCoefSetppc(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   VAR*             var                 /**< variable to add to the constraint */
   )
{
   assert(var != NULL);

   /*debugMessage("adding variable <%s> to setppc constraint <%s>\n", 
     SCIPvarGetName(var), SCIPconsGetName(cons));*/

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      return SCIP_INVALIDDATA;
   }
   
   CHECK_OKAY( addCoef(scip, cons, var) );

   return SCIP_OKAY;
}

/** gets the dual solution of the set partitioning / packing / covering constraint in the current LP */
Real SCIPgetDualsolSetppc(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   )
{
   CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not a set partitioning / packing / covering constraint\n");
      abort();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

