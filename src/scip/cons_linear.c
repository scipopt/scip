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
#pragma ident "@(#) $Id: cons_linear.c,v 1.99 2004/06/01 16:40:14 bzfpfend Exp $"

/**@file   cons_linear.c
 * @brief  constraint handler for linear constraints
 * @author Tobias Achterberg
 *
 *  Linear constraints are separated with a high priority, because they are easy
 *  to separate. Instead of using the global cut pool, the same effect can be
 *  implemented by adding linear constraints to the root node, such that they are
 *  separated each time, the linear constraints are separated. A constraint
 *  handler, which generates linear constraints in this way should have a lower
 *  separation priority than the linear constraint handler, and it should have a
 *  separation frequency that is a multiple of the frequency of the linear
 *  constraint handler. In this way, it can be avoided to separate the same cut
 *  twice, because if a separation run of the handler is always preceded by a
 *  separation of the linear constraints, the priorily added constraints are
 *  always satisfied.
 *
 *  Linear constraints are enforced and checked with a very low priority. Checking
 *  of (many) linear constraints is much more involved than checking the solution
 *  values for integrality. Because we are separating the linear constraints quite
 *  often, it is only necessary to enforce them for integral solutions. A constraint
 *  handler which generates pool cuts in its enforcing method should have an
 *  enforcing priority smaller than that of the linear constraint handler to avoid
 *  regenerating constraints which already exist.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_linear.h"



#define CONSHDLR_NAME          "linear"
#define CONSHDLR_DESC          "linear constraints of the form  lhs <= a^T x <= rhs"
#define CONSHDLR_SEPAPRIORITY  +1000000
#define CONSHDLR_ENFOPRIORITY  -1000000
#define CONSHDLR_CHECKPRIORITY -1000000
#define CONSHDLR_SEPAFREQ             5
#define CONSHDLR_PROPFREQ             5
#define CONSHDLR_EAGERFREQ          100
#define CONSHDLR_NEEDSCONS         TRUE /**< the constraint handler should only be called, if linear constraints exist */

#define DEFAULT_TIGHTENBOUNDSFREQ     1 /**< multiplier on propagation frequency, how often the bounds are tightened */
#define DEFAULT_MAXAGGRNORMSCALE    5.0 /**< maximal allowed relative gain in maximum norm for constraint aggregation */

#define EVENTHDLR_NAME         "linear"
#define EVENTHDLR_DESC         "bound change event handler for linear constraints"

#define CONFLICTHDLR_NAME      "linear"
#define CONFLICTHDLR_DESC      "conflict handler creating linear constraints"
#define CONFLICTHDLR_PRIORITY  -1000000


/** constraint data for linear constraints */
struct ConsData
{
   ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   VAR**            vars;               /**< variables of constraint entries */
   Real*            vals;               /**< coefficients of constraint entries */
   EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   Real             lhs;                /**< left hand side of row (for ranged rows) */
   Real             rhs;                /**< right hand side of row */
   Real             maxabsval;          /**< maximum absolute value of all coefficients */
   Real             pseudoactivity;     /**< pseudo activity value in current pseudo solution */
   Real             minactivity;        /**< minimal value w.r.t. the variable's bounds for the constraint's activity,
                                         *   ignoring the coefficients contributing with infinite value */
   Real             maxactivity;        /**< maximal value w.r.t. the variable's bounds for the constraint's activity,
                                         *   ignoring the coefficients contributing with infinite value */
   int              minactivityinf;     /**< number of coefficients contributing with infinite value to minactivity */
   int              maxactivityinf;     /**< number of coefficients contributing with infinite value to maxactivity */
   int              varssize;           /**< size of the vars- and vals-arrays */
   int              nvars;              /**< number of nonzeros in constraint */
   unsigned int     validmaxabsval:1;   /**< is the maximum absolute value valid? */
   unsigned int     validactivities:1;  /**< are the pseudo activity and activity bounds valid? */
   unsigned int     propagated:1;       /**< is constraint already preprocessed/propagated? */
   unsigned int     changed:1;          /**< was constraint changed since last aggregation round in preprocessing? */
   unsigned int     normalized:1;       /**< is the constraint in normalized form? */
   unsigned int     upgraded:1;         /**< is the constraint upgraded and will it be removed after preprocessing? */
   unsigned int     sorted:1;           /**< are the constraint's variables sorted? */
};

/** event data for bound change event */
struct EventData
{
   CONSDATA*        consdata;           /**< linear constraint data to process the bound change for */
   int              varpos;             /**< position of variable in vars array */
};

/** constraint handler data */
struct ConshdlrData
{
   EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   LINCONSUPGRADE** linconsupgrades;    /**< linear constraint upgrade methods for specializing linear constraints */
   Real             maxaggrnormscale;   /**< maximal allowed relative gain in maximum norm for constraint aggregation */
   int              linconsupgradessize;/**< size of linconsupgrade array */
   int              nlinconsupgrades;   /**< number of linear constraint upgrade methods */
   int              tightenboundsfreq;  /**< multiplier on propagation frequency, how often the bounds are tightened */
};

/** linear constraint update method */
struct LinConsUpgrade
{
   DECL_LINCONSUPGD((*linconsupgd));    /**< method to call for upgrading linear constraint */
   int              priority;           /**< priority of upgrading method */
};




/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that linconsupgrades array can store at least num entries */
static
RETCODE conshdlrdataEnsureLinconsupgradesSize(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA*    conshdlrdata,       /**< linear constraint handler data */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->nlinconsupgrades <= conshdlrdata->linconsupgradessize);
   
   if( num > conshdlrdata->linconsupgradessize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      CHECK_OKAY( SCIPreallocMemoryArray(scip, &conshdlrdata->linconsupgrades, newsize) );
      conshdlrdata->linconsupgradessize = newsize;
   }
   assert(num <= conshdlrdata->linconsupgradessize);

   return SCIP_OKAY;
}

/** ensures, that vars and vals arrays can store at least num entries */
static
RETCODE consdataEnsureVarsSize(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   int              num,                /**< minimum number of entries to store */
   Bool             transformed         /**< is constraint from transformed problem? */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);
   
   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      CHECK_OKAY( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      CHECK_OKAY( SCIPreallocBlockMemoryArray(scip, &consdata->vals, consdata->varssize, newsize) );
      if( transformed )
      {
         CHECK_OKAY( SCIPreallocBlockMemoryArray(scip, &consdata->eventdatas, consdata->varssize, newsize) );
      }
      else
         assert(consdata->eventdatas == NULL);
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}




/*
 * local methods for managing linear constraint update methods
 */

/** creates a linear constraint upgrade data object */
static
RETCODE linconsupgradeCreate(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSUPGRADE** linconsupgrade,     /**< pointer to store the linear constraint upgrade */
   DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   int              priority            /**< priority of upgrading method */
   )
{
   assert(linconsupgrade != NULL);
   assert(linconsupgd != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, linconsupgrade) );
   (*linconsupgrade)->linconsupgd = linconsupgd;
   (*linconsupgrade)->priority = priority;

   return SCIP_OKAY;
}

/** frees a linear constraint upgrade data object */
static
void linconsupgradeFree(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSUPGRADE** linconsupgrade      /**< pointer to the linear constraint upgrade */
   )
{
   assert(linconsupgrade != NULL);
   assert(*linconsupgrade != NULL);

   SCIPfreeMemory(scip, linconsupgrade);
}

/** creates constaint handler data for linear constraint handler */
static
RETCODE conshdlrdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, conshdlrdata) );
   (*conshdlrdata)->linconsupgrades = NULL;
   (*conshdlrdata)->linconsupgradessize = 0;
   (*conshdlrdata)->nlinconsupgrades = 0;
   (*conshdlrdata)->tightenboundsfreq = DEFAULT_TIGHTENBOUNDSFREQ;

   /* get event handler for updating linear constraint activity bounds */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      errorMessage("event handler for linear constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** frees constraint handler data for linear constraint handler */
static
void conshdlrdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   int i;

   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   for( i = 0; i < (*conshdlrdata)->nlinconsupgrades; ++i )
   {
      linconsupgradeFree(scip, &(*conshdlrdata)->linconsupgrades[i]);
   }
   SCIPfreeMemoryArrayNull(scip, &(*conshdlrdata)->linconsupgrades);

   SCIPfreeMemory(scip, conshdlrdata);
}

/** adds a linear constraint update method to the constraint handler's data */
static
RETCODE conshdlrdataIncludeUpgrade(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   LINCONSUPGRADE*  linconsupgrade      /**< linear constraint upgrade method */
   )
{
   int i;

   assert(conshdlrdata != NULL);
   assert(linconsupgrade != NULL);

   CHECK_OKAY( conshdlrdataEnsureLinconsupgradesSize(scip, conshdlrdata, conshdlrdata->nlinconsupgrades+1) );

   for( i = conshdlrdata->nlinconsupgrades;
        i > 0 && conshdlrdata->linconsupgrades[i-1]->priority < linconsupgrade->priority; --i )
   {
      conshdlrdata->linconsupgrades[i] = conshdlrdata->linconsupgrades[i-1];
   }
   assert(0 <= i && i <= conshdlrdata->nlinconsupgrades);
   conshdlrdata->linconsupgrades[i] = linconsupgrade;
   conshdlrdata->nlinconsupgrades++;

   return SCIP_OKAY;
}




/*
 * local methods
 */

/** creates event data for variable at given position, and catches events */
static
RETCODE consdataCatchEvent(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars != NULL);
   assert(consdata->vars[pos] != NULL);
   assert(SCIPvarIsTransformed(consdata->vars[pos]));
   assert(consdata->eventdatas != NULL);
   assert(consdata->eventdatas[pos] == NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, &consdata->eventdatas[pos]) );
   consdata->eventdatas[pos]->consdata = consdata;
   consdata->eventdatas[pos]->varpos = pos;

   CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, 
                  consdata->eventdatas[pos]) );

   return SCIP_OKAY;
}

/** deletes event data for variable at given position, and drops events */
static
RETCODE consdataDropEvent(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars[pos] != NULL);
   assert(consdata->eventdatas[pos] != NULL);
   assert(consdata->eventdatas[pos]->consdata == consdata);
   assert(consdata->eventdatas[pos]->varpos == pos);
   
   CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
                  consdata->eventdatas[pos]) );

   SCIPfreeBlockMemory(scip, &consdata->eventdatas[pos]);

   return SCIP_OKAY;
}

/** catches bound change events for all variables in transformed linear constraint */
static
RETCODE consdataCatchAllEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
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

/** drops bound change events for all variables in transformed linear constraint */
static
RETCODE consdataDropAllEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
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

/** locks the rounding locks associated to the given coefficient in the linear constraint */
static
void consdataLockRounding(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   VAR*             var,                /**< variable of constraint entry */
   Real             val,                /**< coefficient of constraint entry */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));

   if( SCIPisPositive(scip, val) )
   {
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         SCIPvarLock(var, nlockspos, nlocksneg);
      if( !SCIPisInfinity(scip, consdata->rhs) )
         SCIPvarLock(var, nlocksneg, nlockspos);
   }
   else
   {
      if( !SCIPisInfinity(scip, consdata->rhs) )
         SCIPvarLock(var, nlockspos, nlocksneg);
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         SCIPvarLock(var, nlocksneg, nlockspos);
   }
}

/** unlocks the rounding locks associated to the given coefficient in the linear constraint */
static
void consdataUnlockRounding(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   VAR*             var,                /**< variable of constraint entry */
   Real             val,                /**< coefficient of constraint entry */
   int              nunlockspos,        /**< decrease in number of rounding locks for constraint */
   int              nunlocksneg         /**< decrease in number of rounding locks for constraint's negation */
   )
{
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));

   if( SCIPisPositive(scip, val) )
   {
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         SCIPvarUnlock(var, nunlockspos, nunlocksneg);
      if( !SCIPisInfinity(scip, consdata->rhs) )
         SCIPvarUnlock(var, nunlocksneg, nunlockspos);
   }
   else
   {
      if( !SCIPisInfinity(scip, consdata->rhs) )
         SCIPvarUnlock(var, nunlockspos, nunlocksneg);
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         SCIPvarUnlock(var, nunlocksneg, nunlockspos);
   }
}

/** locks the rounding locks of all coefficients in the linear constraint */
static
void consdataLockAllRoundings(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   Bool haslhs;
   Bool hasrhs;
   int i;

   assert(consdata != NULL);

   haslhs = !SCIPisInfinity(scip, -consdata->lhs);
   hasrhs = !SCIPisInfinity(scip, consdata->rhs);
   
   /* lock rounding of every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( SCIPisPositive(scip, consdata->vals[i]) )
      {
         if( haslhs )
            SCIPvarLock(consdata->vars[i], nlockspos, nlocksneg);
         if( hasrhs )
            SCIPvarLock(consdata->vars[i], nlocksneg, nlockspos);
      }
      else
      {
         if( haslhs )
            SCIPvarLock(consdata->vars[i], nlocksneg, nlockspos);
         if( hasrhs )
            SCIPvarLock(consdata->vars[i], nlockspos, nlocksneg);
      }
   }
}

/** unlocks the rounding locks of all coefficients in the linear constraint */
static
void consdataUnlockAllRoundings(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   int              nunlockspos,        /**< decrease in number of rounding locks for constraint */
   int              nunlocksneg         /**< decrease in number of rounding locks for constraint's negation */
   )
{
   Bool haslhs;
   Bool hasrhs;
   int i;

   assert(consdata != NULL);

   haslhs = !SCIPisInfinity(scip, -consdata->lhs);
   hasrhs = !SCIPisInfinity(scip, consdata->rhs);
   
   /* unlock rounding of every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( SCIPisPositive(scip, consdata->vals[i]) )
      {
         if( haslhs )
            SCIPvarUnlock(consdata->vars[i], nunlockspos, nunlocksneg);
         if( hasrhs )
            SCIPvarUnlock(consdata->vars[i], nunlocksneg, nunlockspos);
      }
      else
      {
         if( haslhs )
            SCIPvarUnlock(consdata->vars[i], nunlocksneg, nunlockspos);
         if( hasrhs )
            SCIPvarUnlock(consdata->vars[i], nunlockspos, nunlocksneg);
      }
   }
}

/** creates a linear constraint data of the original problem */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to linear constraint data */
   int              nvars,              /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs                 /**< right hand side of row */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   if( SCIPisGT(scip, lhs, rhs) )
   {
      errorMessage("left hand side of linear constraint greater than right hand side\n");
      errorMessage(" -> lhs=%f, rhs=%f", lhs, rhs);
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );

   if( nvars > 0 )
   {
      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vals, vals, nvars) );
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->vals = NULL;
   }
   (*consdata)->eventdatas = NULL;

   (*consdata)->row = NULL;
   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;
   (*consdata)->maxabsval = SCIP_INVALID;
   (*consdata)->pseudoactivity = SCIP_INVALID;
   (*consdata)->minactivity = SCIP_INVALID;
   (*consdata)->maxactivity = SCIP_INVALID;
   (*consdata)->minactivityinf = -1;
   (*consdata)->maxactivityinf = -1;
   (*consdata)->varssize = nvars;
   (*consdata)->nvars = nvars;
   (*consdata)->validmaxabsval = FALSE;
   (*consdata)->validactivities = FALSE;
   (*consdata)->propagated = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->normalized = FALSE;
   (*consdata)->upgraded = FALSE;
   (*consdata)->sorted = (nvars <= 1);
   
   return SCIP_OKAY;
}

/** creates a linear constraint data of the transformed problem */
static
RETCODE consdataCreateTransformed(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to linear constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              nvars,              /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs                 /**< right hand side of row */
   )
{
   assert(consdata != NULL);

   /* create linear constraint data */
   CHECK_OKAY( consdataCreate(scip, consdata, nvars, vars, vals, lhs, rhs) );

   /* allocate the additional needed eventdatas array */
   assert((*consdata)->eventdatas == NULL);
   CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &(*consdata)->eventdatas, (*consdata)->varssize) );

   /* get transformed variables */
   CHECK_OKAY( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

   /* initialize the eventdatas array */
   clearMemoryArray((*consdata)->eventdatas, (*consdata)->nvars);

   /* catch bound change events of variables */
   CHECK_OKAY( consdataCatchAllEvents(scip, *consdata, eventhdlr) );

   return SCIP_OKAY;
}

/** frees a linear constraint data */
static
RETCODE consdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to linear constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->varssize >= 0);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* free event datas */
   if( (*consdata)->eventdatas != NULL )
   {
      /* drop bound change events of variables */
      CHECK_OKAY( consdataDropAllEvents(scip, *consdata, eventhdlr) );

      /* free additional eventdatas array */
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->eventdatas, (*consdata)->varssize);
   }
   assert((*consdata)->eventdatas == NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vals, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints linear constraint to file stream */
static
void consdataPrint(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);

   if( file == NULL )
      file = stdout;

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      fprintf(file, "%+g <= ", consdata->lhs);

   /* print coefficients */
   if( consdata->nvars == 0 )
      fprintf(file, "0 ");
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      fprintf(file, "%+g%s ", consdata->vals[v], SCIPvarGetName(consdata->vars[v]));
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      fprintf(file, "= %+g\n", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      fprintf(file, "<= %+g\n", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      fprintf(file, ">= %+g\n", consdata->lhs);
   else
      fprintf(file, " [free]\n");
}

/** updates minimum and maximum activity for a change in lower bound */
static
void consdataUpdateChgLb(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   VAR*             var,                /**< variable that has been changed */
   Real             oldlb,              /**< old lower bound of variable */
   Real             newlb,              /**< new lower bound of variable */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      assert(consdata->pseudoactivity < SCIP_INVALID);
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(consdata->minactivityinf >= 0);
      assert(consdata->maxactivityinf >= 0);
      assert(!SCIPisInfinity(scip, oldlb));
      assert(!SCIPisInfinity(scip, newlb));

      if( SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_LOWER )
         consdata->pseudoactivity += val * (newlb - oldlb);

      if( val > 0.0 )
      {
         if( SCIPisInfinity(scip, -oldlb) )
         {
            assert(consdata->minactivityinf >= 1);
            consdata->minactivityinf--;
         }
         else
            consdata->minactivity -= val * oldlb;

         if( SCIPisInfinity(scip, -newlb) )
            consdata->minactivityinf++;
         else
            consdata->minactivity += val * newlb;
      }
      else
      {
         if( SCIPisInfinity(scip, -oldlb) )
         {
            assert(consdata->maxactivityinf >= 1);
            consdata->maxactivityinf--;
         }
         else
            consdata->maxactivity -= val * oldlb;

         if( SCIPisInfinity(scip, -newlb) )
            consdata->maxactivityinf++;
         else
            consdata->maxactivity += val * newlb;
      }
   }
}

/** updates minimum and maximum activity for a change in upper bound */
static
void consdataUpdateChgUb(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   VAR*             var,                /**< variable that has been changed */
   Real             oldub,              /**< old upper bound of variable */
   Real             newub,              /**< new upper bound of variable */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      assert(consdata->pseudoactivity < SCIP_INVALID);
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(!SCIPisInfinity(scip, -oldub));
      assert(!SCIPisInfinity(scip, -newub));

      if( SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_UPPER )
         consdata->pseudoactivity += val * (newub - oldub);

      if( val > 0.0 )
      {
         if( SCIPisInfinity(scip, oldub) )
         {
            assert(consdata->maxactivityinf >= 1);
            consdata->maxactivityinf--;
         }
         else
            consdata->maxactivity -= val * oldub;

         if( SCIPisInfinity(scip, newub) )
            consdata->maxactivityinf++;
         else
            consdata->maxactivity += val * newub;
      }
      else
      {
         if( SCIPisInfinity(scip, oldub) )
         {
            assert(consdata->minactivityinf >= 1);
            consdata->minactivityinf--;
         }
         else
            consdata->minactivity -= val * oldub;

         if( SCIPisInfinity(scip, newub) )
            consdata->minactivityinf++;
         else
            consdata->minactivity += val * newub;
      }
   }
}

/** updates minimum and maximum activity and maximum absolute value for coefficient addition */
static
void consdataUpdateAddCoef(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(consdata != NULL);

   /* update maximum absolute value */
   if( consdata->validmaxabsval )
   {
      Real absval;

      assert(consdata->maxabsval < SCIP_INVALID);

      absval = ABS(val);
      consdata->maxabsval = MAX(consdata->maxabsval, absval);
   }

   /* update pseudo, minimal and maximal activity */
   if( consdata->validactivities )
   {
      assert(consdata->pseudoactivity < SCIP_INVALID);
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);

      consdataUpdateChgLb(scip, consdata, var, 0.0, SCIPvarGetLbLocal(var), val);
      consdataUpdateChgUb(scip, consdata, var, 0.0, SCIPvarGetUbLocal(var), val);
   }
}

/** updates minimum and maximum activity for coefficient deletion, invalidates maximum absolute value if necessary */
static
void consdataUpdateDelCoef(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(consdata != NULL);

   /* invalidate maximum absolute value, if this coefficient was the maximum */
   if( consdata->validmaxabsval )
   {
      if( val == consdata->maxabsval ) /*lint !e777*/ /* check for equality without epsilon is correct here! */
      {
         consdata->validmaxabsval = FALSE;
         consdata->maxabsval = SCIP_INVALID;
      }
   }

   /* update pseudo, minimal and maximal activity */
   if( consdata->validactivities )
   {
      assert(consdata->pseudoactivity < SCIP_INVALID);
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);

      consdataUpdateChgLb(scip, consdata, var, SCIPvarGetLbLocal(var), 0.0, val);
      consdataUpdateChgUb(scip, consdata, var, SCIPvarGetUbLocal(var), 0.0, val);
   }
}

/** updates minimum and maximum activity for coefficient change, invalidates maximum absolute value if necessary */
static
void consdataUpdateChgCoef(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   VAR*             var,                /**< variable of constraint entry */
   Real             oldval,             /**< old coefficient of constraint entry */
   Real             newval              /**< new coefficient of constraint entry */
   )
{
   consdataUpdateDelCoef(scip, consdata, var, oldval);
   consdataUpdateAddCoef(scip, consdata, var, newval);
}

/** calculates maximum absolute value of coefficients */
static
void consdataCalcMaxAbsval(
   CONSDATA*        consdata            /**< linear constraint data */
   )
{
   Real absval;
   int i;

   assert(consdata != NULL);
   assert(!consdata->validmaxabsval);
   assert(consdata->maxabsval >= SCIP_INVALID);

   consdata->validmaxabsval = TRUE;
   consdata->maxabsval = 0.0;
   for( i = 0; i < consdata->nvars; ++i )
   {
      absval = consdata->vals[i];
      absval = ABS(absval);
      consdata->maxabsval = MAX(consdata->maxabsval, absval);
   }   
}

/** returns the maximum absolute value of all coefficients in the constraint */
static
Real consdataGetMaxAbsval(
   CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->validmaxabsval )
      consdataCalcMaxAbsval(consdata);
   assert(consdata->validmaxabsval);
   assert(consdata->maxabsval < SCIP_INVALID);

   return consdata->maxabsval;
}

/** calculates pseudo activity, and minimum and maximum activity for constraint;
 *  additionally recalculates maximum absolute value of coefficients
 */
static
void consdataCalcActivities(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   
   assert(consdata != NULL);
   assert(!consdata->validactivities);
   assert(consdata->pseudoactivity >= SCIP_INVALID);
   assert(consdata->minactivity >= SCIP_INVALID);
   assert(consdata->maxactivity >= SCIP_INVALID);
   
   consdata->validmaxabsval = TRUE;
   consdata->validactivities = TRUE;
   consdata->maxabsval = 0.0;
   consdata->pseudoactivity = 0.0;
   consdata->minactivity = 0.0;
   consdata->maxactivity = 0.0;
   consdata->minactivityinf = 0;
   consdata->maxactivityinf = 0;

   for( i = 0; i < consdata->nvars; ++i )
      consdataUpdateAddCoef(scip, consdata, consdata->vars[i], consdata->vals[i]);
}

/** gets activity bounds for constraint */
static
Real consdataGetPseudoActivity(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< linear constraint */
   )
{
   assert(consdata != NULL);

   if( !consdata->validactivities )
      consdataCalcActivities(scip, consdata);
   assert(consdata->pseudoactivity < SCIP_INVALID);
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);

   debugMessage("pseudo activity of linear constraint: %g\n", consdata->pseudoactivity);

   return consdata->pseudoactivity;
}

/** calculates the feasibility of the linear constraint for given solution */
static
Real consdataGetPseudoFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< linear constraint data */
   )
{
   Real activity;

   assert(consdata != NULL);

   activity = consdataGetPseudoActivity(scip, consdata);

   return MIN(consdata->rhs - activity, activity - consdata->lhs);
}

/** gets activity bounds for constraint */
static
void consdataGetActivityBounds(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint */
   Real*            minactivity,        /**< pointer to store the minimal activity */
   Real*            maxactivity         /**< pointer to store the maximal activity */
   )
{
   assert(consdata != NULL);
   assert(minactivity != NULL);
   assert(maxactivity != NULL);

   if( !consdata->validactivities )
      consdataCalcActivities(scip, consdata);
   assert(consdata->pseudoactivity < SCIP_INVALID);
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);

   if( consdata->minactivityinf > 0 )
      *minactivity = -SCIPinfinity(scip);
   else
      *minactivity = consdata->minactivity;
   if( consdata->maxactivityinf > 0 )
      *maxactivity = SCIPinfinity(scip);
   else
      *maxactivity = consdata->maxactivity;
}

/** gets activity bounds for constraint after setting variable to zero */
static
void consdataGetActivityResiduals(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint */
   VAR*             var,                /**< variable to calculate activity residual for */
   Real             val,                /**< coefficient value of variable in linear constraint */
   Real*            minresactivity,     /**< pointer to store the minimal residual activity */
   Real*            maxresactivity      /**< pointer to store the maximal residual activity */
   )
{
   Real lb;
   Real ub;
   
   assert(consdata != NULL);
   assert(var != NULL);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);

   /* get activity bounds of linear constraint */
   if( !consdata->validactivities )
      consdataCalcActivities(scip, consdata);
   assert(consdata->pseudoactivity < SCIP_INVALID);
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);
   assert(consdata->minactivityinf >= 0);
   assert(consdata->maxactivityinf >= 0);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(!SCIPisInfinity(scip, lb));
   assert(!SCIPisInfinity(scip, -ub));

   if( val > 0.0 )
   {
      if( SCIPisInfinity(scip, -lb) )
      {
         assert(consdata->minactivityinf >= 1);
         if( consdata->minactivityinf >= 2 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->minactivity;
      }
      else
      {
         if( consdata->minactivityinf >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->minactivity - val * lb;
      }
      if( SCIPisInfinity(scip, ub) )
      {
         assert(consdata->maxactivityinf >= 1);
         if( consdata->maxactivityinf >= 2 )
            *maxresactivity = +SCIPinfinity(scip);
         else
            *maxresactivity = consdata->maxactivity;
      }
      else
      {
         if( consdata->maxactivityinf >= 1 )
            *maxresactivity = +SCIPinfinity(scip);
         else
            *maxresactivity = consdata->maxactivity - val * ub;
      }
   }
   else
   {
      if( SCIPisInfinity(scip, ub) )
      {
         assert(consdata->minactivityinf >= 1);
         if( consdata->minactivityinf >= 2 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->minactivity;
      }
      else
      {
         if( consdata->minactivityinf >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->minactivity - val * ub;
      }
      if( SCIPisInfinity(scip, -lb) )
      {
         assert(consdata->maxactivityinf >= 1);
         if( consdata->maxactivityinf >= 2 )
            *maxresactivity = +SCIPinfinity(scip);
         else
            *maxresactivity = consdata->maxactivity;
      }
      else
      {
         if( consdata->maxactivityinf >= 1 )
            *maxresactivity = +SCIPinfinity(scip);
         else
            *maxresactivity = consdata->maxactivity - val * lb;
      }
   }
}

/** invalidates pseudo activity and activity bounds, such that they are recalculated in next get */
static
void consdataInvalidateActivities(
   CONSDATA*        consdata            /**< linear constraint */
   )
{
   assert(consdata != NULL);

   consdata->validactivities = FALSE;
   consdata->pseudoactivity = SCIP_INVALID;
   consdata->minactivity = SCIP_INVALID;
   consdata->maxactivity = SCIP_INVALID;
   consdata->minactivityinf = -1;
   consdata->maxactivityinf = -1;
}

/** calculates the activity of the linear constraint for given solution */
static
Real consdataGetActivity(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   SOL*             sol                 /**< solution to get activity for, NULL to current solution */
   )
{
   Real activity;
   Real infinity;

   assert(consdata != NULL);

   if( sol == NULL && !SCIPhasActNodeLP(scip) )
   {
      /* for performance reasons, the pseudo activity is updated with each bound change, so we don't have to
       * recalculate it
       */
      activity = consdataGetPseudoActivity(scip, consdata);
   }
   else
   {
      Real solval;
      int v;

      activity = 0.0;
      for( v = 0; v < consdata->nvars; ++v )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[v]);
         activity += consdata->vals[v] * solval;
      }

      debugMessage("activity of linear constraint: %g\n", activity);
   }

   infinity = SCIPinfinity(scip);
   activity = MAX(activity, -infinity);
   activity = MIN(activity, +infinity);

   return activity;
}

/** calculates the feasibility of the linear constraint for given solution */
static
Real consdataGetFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< linear constraint data */
   SOL*             sol                 /**< solution to get feasibility for, NULL to current solution */
   )
{
   Real activity;

   assert(consdata != NULL);

   activity = consdataGetActivity(scip, consdata, sol);

   return MIN(consdata->rhs - activity, activity - consdata->lhs);
}

#define BOUNDSCALETOL 1e-5
/** tightens bounds of a single variable due to activity bounds */
static
RETCODE consdataTightenVarBounds(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< constraint data */
   VAR*             var,                /**< variable to tighten bounds for */
   Real             val,                /**< coefficient value of variable in linear constraint */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*             nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   Real lb;
   Real ub;
   Real newlb;
   Real newub;
   Real minresactivity;
   Real maxresactivity;
   Real lhs;
   Real rhs;

   assert(consdata != NULL);
   assert(var != NULL);
   assert(!SCIPisZero(scip, val));
   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   *cutoff = FALSE;

   lhs = consdata->lhs;
   rhs = consdata->rhs;
   consdataGetActivityResiduals(scip, consdata, var, val, &minresactivity, &maxresactivity);
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));
   assert(!SCIPisInfinity(scip, minresactivity));
   assert(!SCIPisInfinity(scip, -maxresactivity));

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPisLE(scip, lb, ub));

   if( val > 0.0 )
   {
      /* check, if we can tighten the variable's bounds */
      if( !SCIPisInfinity(scip, -minresactivity) && !SCIPisInfinity(scip, rhs) )
      {
         newub = (rhs - minresactivity)/val;
         newub = SCIPceil(scip, newub/BOUNDSCALETOL) * BOUNDSCALETOL;
         if( SCIPisUbBetter(scip, newub, ub) )
         {
            /* tighten upper bound */
            debugMessage("linear constraint: tighten <%s>, old bds=[%.9f,%.9f], val=%g, resactivity=[%g,%g], sides=[%g,%g]\n",
               SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs);
            if( SCIPisFeasLT(scip, newub, lb) )
            {
               debugMessage("linear constraint: cutoff  <%s>, new bds=[%.9f,%.9f]\n", SCIPvarGetName(var), lb, newub);
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            else if( newub <= lb )
               newub = lb;  /* avoid infeasibilities in consequence of numerical inaccuracies */
            CHECK_OKAY( SCIPchgVarUb(scip, var, newub) );
            ub = SCIPvarGetUbLocal(var); /* get bound again, because it may be additionally modified due to integrality */
            assert(SCIPisFeasLE(scip, ub, newub));
            (*nchgbds)++;
            debugMessage("linear constraint: tighten <%s>, new bds=[%.9f,%.9f]\n", SCIPvarGetName(var), lb, ub);
         }
      }
      if( !SCIPisInfinity(scip, maxresactivity) && !SCIPisInfinity(scip, -lhs) )
      {
         newlb = (lhs - maxresactivity)/val;
         newlb = SCIPfloor(scip, newlb/BOUNDSCALETOL) * BOUNDSCALETOL;
         if( SCIPisLbBetter(scip, newlb, lb) )
         {
            /* tighten lower bound */
            debugMessage("linear constraint: tighten <%s>, old bds=[%.9f,%.9f], val=%g, resactivity=[%g,%g], sides=[%g,%g]\n",
               SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs);
            if( SCIPisFeasGT(scip, newlb, ub) )
            {
               debugMessage("linear constraint: cutoff  <%s>, new bds=[%.9f,%.9f]\n", SCIPvarGetName(var), newlb, ub);
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            else if( newlb >= ub )
               newlb = ub; /* avoid infeasibilities in consequence of numerical inaccuracies */
            CHECK_OKAY( SCIPchgVarLb(scip, var, newlb) );
            lb = SCIPvarGetLbLocal(var); /* get bound again, because it may be additionally modified due to integrality */
            assert(SCIPisFeasGE(scip, lb, newlb));
            (*nchgbds)++;
            debugMessage("linear constraint: tighten <%s>, new bds=[%.9f,%.9f]\n", SCIPvarGetName(var), lb, ub);
         }
      }
   }
   else
   {
      /* check, if we can tighten the variable's bounds */
      if( !SCIPisInfinity(scip, -minresactivity) && !SCIPisInfinity(scip, rhs) )
      {
         newlb = (rhs - minresactivity)/val;
         newlb = SCIPfloor(scip, newlb/BOUNDSCALETOL) * BOUNDSCALETOL;
         if( SCIPisLbBetter(scip, newlb, lb) )
         {
            /* tighten lower bound */
            debugMessage("linear constraint: tighten <%s>, old bds=[%.9f,%.9f], val=%g, resactivity=[%g,%g], sides=[%g,%g]\n",
               SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs);
            if( SCIPisFeasGT(scip, newlb, ub) )
            {
               debugMessage("linear constraint: cutoff  <%s>, new bds=[%.9f,%.9f]\n", SCIPvarGetName(var), newlb, ub);
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            else if( newlb >= ub )
               newlb = ub; /* avoid infeasibilities in consequence of numerical inaccuracies */
            CHECK_OKAY( SCIPchgVarLb(scip, var, newlb) );
            lb = SCIPvarGetLbLocal(var); /* get bound again, because it may be additionally modified due to integrality */
            assert(SCIPisFeasGE(scip, lb, newlb));
            (*nchgbds)++;
            debugMessage("linear constraint: tighten <%s>, new bds=[%.9f,%.9f]\n", SCIPvarGetName(var), lb, ub);
         }
      }
      if( !SCIPisInfinity(scip, maxresactivity) && !SCIPisInfinity(scip, -lhs) )
      {
         newub = (lhs - maxresactivity)/val;
         newub = SCIPceil(scip, newub/BOUNDSCALETOL) * BOUNDSCALETOL;
         if( SCIPisUbBetter(scip, newub, ub) )
         {
            /* tighten upper bound */
            debugMessage("linear constraint: tighten <%s>, old bds=[%.9f,%.9f], val=%g, resactivity=[%g,%g], sides=[%g,%g]\n",
               SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs);
            if( SCIPisFeasLT(scip, newub, lb) )
            {
               debugMessage("linear constraint: cutoff  <%s>, new bds=[%.9f,%.9f]\n", SCIPvarGetName(var), lb, newub);
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            else if( newub <= lb )
               newub = lb; /* avoid infeasibilities in consequence of numerical inaccuracies */
            CHECK_OKAY( SCIPchgVarUb(scip, var, newub) );
            ub = SCIPvarGetUbLocal(var); /* get bound again, because it may be additionally modified due to integrality */
            assert(SCIPisFeasLE(scip, ub, newub));
            (*nchgbds)++;
            debugMessage("linear constraint: tighten <%s>, new bds=[%.9f,%.9f]\n", SCIPvarGetName(var), lb, ub);
         }
      }
   }
   
   return SCIP_OKAY;
}

/** index comparison method of linear constraints: compares two indices of the variable set in the linear constraint */
static
DECL_SORTINDCOMP(consdataCmpVar)
{  /*lint --e{715}*/
   CONSDATA* consdata = (CONSDATA*)dataptr;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nvars);
   assert(0 <= ind2 && ind2 < consdata->nvars);
   
   return SCIPvarCmp(consdata->vars[ind1], consdata->vars[ind2]);
}

/** sorts linear constraint's variables */
static
RETCODE consdataSort(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->nvars == 0 )
      consdata->sorted = TRUE;
   else if( !consdata->sorted )
   {
      VAR* varv;
      EVENTDATA* eventdatav;
      Real valv;
      int* perm;
      int v;
      int i;
      int nexti;

      /* get temporary memory to store the sorted permutation */
      CHECK_OKAY( SCIPallocBufferArray(scip, &perm, consdata->nvars) );

      /* call bubble sort */
      SCIPbsort((void*)consdata, consdata->nvars, consdataCmpVar, perm);

      /* permute the variables in the linear constraint according to the resulting permutation */
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( perm[v] != v )
         {
            varv = consdata->vars[v];
            valv = consdata->vals[v];
            eventdatav = consdata->eventdatas[v];
            i = v;
            do
            {
               assert(0 <= perm[i] && perm[i] < consdata->nvars);
               assert(perm[i] != i);
               consdata->vars[i] = consdata->vars[perm[i]];
               consdata->vals[i] = consdata->vals[perm[i]];
               consdata->eventdatas[i] = consdata->eventdatas[perm[i]];
               consdata->eventdatas[i]->varpos = i;
               nexti = perm[i];
               perm[i] = i;
               i = nexti;
            }
            while( perm[i] != v );
            consdata->vars[i] = varv;
            consdata->vals[i] = valv;
            consdata->eventdatas[i] = eventdatav;
            consdata->eventdatas[i]->varpos = i;
            perm[i] = i;
         }
      }
      consdata->sorted = TRUE;

#ifdef DEBUG
      /* check sorting */
      for( v = 0; v < consdata->nvars; ++v )
      {
         assert(v == consdata->nvars-1 || SCIPvarCmp(consdata->vars[v], consdata->vars[v+1]) <= 0);
         assert(perm[v] == v);
         assert(consdata->eventdatas[v]->varpos == v);
      }
#endif

      /* free temporary memory */
      CHECK_OKAY( SCIPfreeBufferArray(scip, &perm) );
   }
   assert(consdata->sorted);

   return SCIP_OKAY;
}




/*
 * local linear constraint handler methods
 */

/** sets left hand side of linear constraint */
static
RETCODE chgLhs(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   Real             lhs                 /**< new left hand side */
   )
{
   CONSDATA* consdata;

   assert(!SCIPisInfinity(scip, lhs));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->vars != NULL && consdata->vals != NULL));
   assert(!SCIPisInfinity(scip, consdata->lhs));

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      if( SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, -lhs) )
      {
         VAR** vars;
         Real* vals;
         int v;
   
         /* the left hand side switched from -infinity to a non-infinite value -> forbid rounding */
         vars = consdata->vars;
         vals = consdata->vals;
         
         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));
            
            if( SCIPisPositive(scip, vals[v]) )
               SCIPvarLockDownCons(vars[v], cons);
            else
               SCIPvarLockUpCons(vars[v], cons);
         }
      }
      else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, -lhs) )
      {
         VAR** vars;
         Real* vals;
         int v;
   
         /* the left hand side switched from a non-infinte value to -infinity -> allow rounding */
         vars = consdata->vars;
         vals = consdata->vals;
         
         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));
            
            if( SCIPisPositive(scip, vals[v]) )
               SCIPvarUnlockDownCons(vars[v], cons);
            else
               SCIPvarUnlockUpCons(vars[v], cons);
         }
      }
   }

   /* set new left hand side */
   consdata->lhs = lhs;
   consdata->propagated = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;

   /* update the lhs of the LP row */
   if( consdata->row != NULL )
   {
      CHECK_OKAY( SCIPchgRowLhs(scip, consdata->row, lhs) );
   }

   return SCIP_OKAY;
}

/** sets right hand side of linear constraint */
static
RETCODE chgRhs(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   Real             rhs                 /**< new right hand side */
   )
{
   CONSDATA* consdata;

   assert(!SCIPisInfinity(scip, -rhs));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->vars != NULL && consdata->vals != NULL));
   assert(!SCIPisInfinity(scip, -consdata->rhs));

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));

      if( SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, rhs) )
      {
         VAR** vars;
         Real* vals;
         int v;
   
         /* the right hand side switched from infinity to a non-infinite value -> forbid rounding */
         vars = consdata->vars;
         vals = consdata->vals;
         
         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));
            
            if( SCIPisPositive(scip, vals[v]) )
               SCIPvarLockUpCons(vars[v], cons);
            else
               SCIPvarLockDownCons(vars[v], cons);
         }
      }
      else if( !SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, rhs) )
      {
         VAR** vars;
         Real* vals;
         int v;
   
         /* the right hand side switched from a non-infinte value to infinity -> allow rounding */
         vars = consdata->vars;
         vals = consdata->vals;
         
         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));
            
            if( SCIPisPositive(scip, vals[v]) )
               SCIPvarUnlockUpCons(vars[v], cons);
            else
               SCIPvarUnlockDownCons(vars[v], cons);
         }
      }
   }

   /* set new right hand side */
   consdata->rhs = rhs;
   consdata->propagated = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;

   /* update the rhs of the LP row */
   if( consdata->row != NULL )
   {
      CHECK_OKAY( SCIPchgRowRhs(scip, consdata->row, rhs) );
   }

   return SCIP_OKAY;
}

/** adds coefficient in linear constraint */
static
RETCODE addCoef(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   CONSDATA* consdata;
   Bool transformed;

   assert(var != NULL);
   assert(!SCIPisZero(scip, val));

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

   CHECK_OKAY( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1, transformed) );
   consdata->vars[consdata->nvars] = var;
   consdata->vals[consdata->nvars] = val;
   consdata->nvars++;

   /* if we are in transformed problem, the variable needs an additional event data */
   if( transformed )
   {
      CONSHDLR* conshdlr;
      CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* initialize eventdatas array */
      consdata->eventdatas[consdata->nvars-1] = NULL;

      /* catch bound change events of variable */
      CHECK_OKAY( consdataCatchEvent(scip, consdata, conshdlrdata->eventhdlr, consdata->nvars-1) );

      /* update minimum and maximum activities */
      consdataUpdateAddCoef(scip, consdata, var, val);
   }

   /* if necessary, update the rounding locks of variable */
   if( SCIPconsIsLocked(cons) )
   {
      assert(transformed);
      consdataLockRounding(scip, consdata, var, val, (int)SCIPconsIsLockedPos(cons), (int)SCIPconsIsLockedNeg(cons));
   }

   consdata->propagated = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   if( consdata->nvars == 1 )
      consdata->sorted = TRUE;
   else
      consdata->sorted = consdata->sorted
         && (SCIPvarCmp(consdata->vars[consdata->nvars-2], consdata->vars[consdata->nvars-1]) == -1);

   /* add the new coefficient to the LP row */
   if( consdata->row != NULL )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->row, var, val) );
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from linear constraint data */
static
RETCODE delCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int              pos                 /**< position of coefficient to delete */
   )
{
   CONSDATA* consdata;
   VAR* var;
   Real val;
   Bool transformed;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   var = consdata->vars[pos];
   val = consdata->vals[pos];
   assert(var != NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* if necessary, update the rounding locks of variable */
   if( SCIPconsIsLocked(cons) )
   {
      assert(transformed);
      consdataUnlockRounding(scip, consdata, var, val, (int)SCIPconsIsLockedPos(cons), (int)SCIPconsIsLockedNeg(cons));
   }

   /* if we are in transformed problem, delete the event data of the variable */
   if( transformed )
   {
      CONSHDLR* conshdlr;
      CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* update minimum and maximum activities */
      consdataUpdateDelCoef(scip, consdata, var, val);

      /* drop bound change events of variable */
      CHECK_OKAY( consdataDropEvent(scip, consdata, conshdlrdata->eventhdlr, pos) );
      assert(consdata->eventdatas[pos] == NULL);
   }

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->vals[pos] = consdata->vals[consdata->nvars-1];
   if( pos != consdata->nvars-1 )
   {
      if( transformed )
      {
         consdata->eventdatas[pos] = consdata->eventdatas[consdata->nvars-1];
         assert(consdata->eventdatas[pos] != NULL);
         consdata->eventdatas[pos]->varpos = pos;
      }
      consdata->sorted = FALSE;
   }
   consdata->nvars--;

   consdata->propagated = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;

   return SCIP_OKAY;
}

/** changes coefficient value at given position of linear constraint data */
static
RETCODE chgCoefPos(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int              pos,                /**< position of coefficient to delete */
   Real             newval              /**< new value of coefficient */
   )
{
   CONSDATA* consdata;
   VAR* var;
   Real val;

   assert(!SCIPisZero(scip, newval));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   var = consdata->vars[pos];
   val = consdata->vals[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   if( SCIPconsIsTransformed(cons) )
   {
      /* update minimum and maximum activities */
      consdataUpdateDelCoef(scip, consdata, var, val);
      consdataUpdateAddCoef(scip, consdata, var, newval);
   }

   /* if necessary, update the rounding locks of the variable */
   if( SCIPconsIsLocked(cons) && newval * val < 0.0 )
   {
      assert(SCIPconsIsTransformed(cons));
      consdataUnlockRounding(scip, consdata, var, val, (int)SCIPconsIsLockedPos(cons), (int)SCIPconsIsLockedNeg(cons));
      consdataLockRounding(scip, consdata, var, newval, (int)SCIPconsIsLockedPos(cons), (int)SCIPconsIsLockedNeg(cons));
   }

   /* change the value */
   consdata->vals[pos] = newval;

   consdata->propagated = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;

   return SCIP_OKAY;
}

/** scales a linear constraint with a constant scalar */
static
RETCODE scaleCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint to scale */
   Real             scalar              /**< value to scale constraint with */
   )
{
   CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   /* scale the coefficients */
   for( i = 0; i < consdata->nvars; ++i )
   {
      Real oldval;
      
      oldval = consdata->vals[i];
      consdata->vals[i] *= scalar;
      if( SCIPisIntegral(scip, consdata->vals[i]) )
         consdata->vals[i] = SCIPfloor(scip, consdata->vals[i]);
      if( SCIPisZero(scip, consdata->vals[i]) )
      {
         warningMessage("coefficient of variable <%s> in linear constraint <%s> scaled to zero (scale: %g)\n", 
            SCIPvarGetName(consdata->vars[i]), SCIPconsGetName(cons), scalar);
         consdata->vals[i] = oldval;
         CHECK_OKAY( delCoefPos(scip, cons, i) );
         --i;
      }
   }

   /* scale the sides */
   if( scalar < 0.0 )
   {
      Real lhs;
      lhs = consdata->lhs;
      consdata->lhs = -consdata->rhs;
      consdata->rhs = -lhs;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      consdata->lhs *= ABS(scalar);
      if( SCIPisIntegral(scip, consdata->lhs) )
         consdata->lhs = SCIPfloor(scip, consdata->lhs);
   }
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      consdata->rhs *= ABS(scalar);
      if( SCIPisIntegral(scip, consdata->rhs) )
         consdata->rhs = SCIPfloor(scip, consdata->rhs);
   }

   consdataInvalidateActivities(consdata);

   return SCIP_OKAY;
}

/** normalizes a linear constraint with the following rules:
 *  - multiplication with +1 or -1:
 *      Apply the following rules in the given order, until the sign of the factor is determined. Later rules only apply,
 *      if the current rule doesn't determine the sign):
 *        1. the right hand side must not be negative
 *        2. the right hand side must not be infinite
 *        3. the absolute value of the right hand side must be greater than that of the left hand side
 *        4. the number of positive coefficients must not be smaller than the number of negative coefficients
 *        5. multiply with +1
 *  - rationals to integrals
 *      Try to identify a rational representation of the fractional coefficients, and multiply all coefficients
 *      by the smallest common multiple of all denominators to get integral coefficients.
 *      Forbid large denominators due to numerical stability.
 *  - division by greatest common divisor
 *      If all coefficients are integral, divide them by the greatest common divisor.
 */
static
RETCODE normalizeCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< linear constraint to normalize */
   )
{
   CONSDATA* consdata;
   VAR** vars;
   Real* vals;
   Longint scm;
   Longint nominator;
   Longint denominator;
   Longint gcd;
   Longint maxmult;
   Real epsilon;
   Real feastol;
   Bool success;
   int nvars;
   int mult;
   int nposcoeffs;
   int nnegcoeffs;
   int i;

   /* we must not change a modifiable constraint in any way */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check, if the constraint is already normalized */
   if( consdata->normalized )
      return SCIP_OKAY;

   /* get coefficient arrays */
   vars = consdata->vars;
   vals = consdata->vals;
   nvars = consdata->nvars;
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   /* calculate the maximal multiplier for common divisor calculation:
    *   |p/q - val| < epsilon  and  q < feastol/epsilon  =>  |p - q*val| < feastol
    * which means, a value of feastol/epsilon should be used as maximal multiplier
    */
   epsilon = SCIPepsilon(scip);
   feastol = SCIPfeastol(scip);
   maxmult = (Longint)(feastol/epsilon + feastol);

   /*
    * multiplication with +1 or -1
    */
   mult = 0;
   
   /* 1. the right hand side must not be negative */
   if( SCIPisPositive(scip, consdata->lhs) )
      mult = +1;
   else if( SCIPisNegative(scip, consdata->rhs) )
      mult = -1;

   if( mult == 0 )
   {
      /* 2. the right hand side must not be infinite */
      if( SCIPisInfinity(scip, -consdata->lhs) )
         mult = +1;
      else if( SCIPisInfinity(scip, consdata->rhs) )
         mult = -1;
   }

   if( mult == 0 )
   {
      /* 3. the absolute value of the right hand side must be greater than that of the left hand side */
      if( SCIPisGT(scip, ABS(consdata->rhs), ABS(consdata->lhs)) )
         mult = +1;
      else if( SCIPisLT(scip, ABS(consdata->rhs), ABS(consdata->lhs)) )
         mult = -1;
   }
   
   if( mult == 0 )
   {
      /* 4. the number of positive coefficients must not be smaller than the number of negative coefficients */
      nposcoeffs = 0;
      nnegcoeffs = 0;
      for( i = 0; i < nvars; ++i )
      {
         if( vals[i] > 0.0 )
            nposcoeffs++;
         else
            nnegcoeffs++;
      }
      if( nposcoeffs > nnegcoeffs )
         mult = +1;
      else if( nposcoeffs < nnegcoeffs )
         mult = -1;
   }

   if( mult == 0 )
   {
      /* 5. multiply with +1 */
      mult = +1;
   }

   assert(mult == +1 || mult == -1);
   if( mult == -1 )
   {
      /* scale the constraint with -1 */
      debugMessage("multiply linear constraint with -1.0\n");
      debug(consdataPrint(scip, consdata, NULL));
      CHECK_OKAY( scaleCons(scip, cons, -1.0) );
   }

   /*
    * rationals to integrals
    */
   success = TRUE;
   scm = 1;
   for( i = 0; i < nvars && success && scm <= maxmult; ++i )
   {
      if( !SCIPisIntegral(scip, vals[i]) )
      {
         success = SCIPrealToRational(vals[i], epsilon, maxmult, &nominator, &denominator);
         if( success )
            scm = SCIPcalcSmaComMul(scm, denominator);
      }
   }
   assert(scm >= 1);
   success = success && (scm <= maxmult);
   if( success && scm != 1 )
   {
      /* scale the constraint with the smallest common multiple of all denominators */
      debugMessage("scale linear constraint with %lld to make coefficients integral\n", scm);
      debug(consdataPrint(scip, consdata, NULL));
      CHECK_OKAY( scaleCons(scip, cons, (Real)scm) );
   }

   /*
    * division by greatest common divisor
    */
   if( success && nvars >= 1 )
   {
      /* all coefficients are integral: divide them by their greatest common divisor */
      assert(SCIPisIntegral(scip, vals[0]));
      gcd = (Longint)(ABS(vals[0]) + feastol);
      assert(gcd >= 1);
      for( i = 1; i < nvars && gcd > 1; ++i )
      {
         assert(SCIPisIntegral(scip, vals[i]));
         gcd = SCIPcalcGreComDiv(gcd, (Longint)(ABS(vals[i]) + feastol));
      }

      if( gcd > 1 )
      {
         /* divide the constaint by the greatest common divisor of the coefficients */
         debugMessage("divide linear constraint by greatest common divisor %lld\n", gcd);
         debug(consdataPrint(scip, consdata, NULL));
         CHECK_OKAY( scaleCons(scip, cons, 1.0/(Real)gcd) );
      }
   }

   /* mark constraint to be normalized */
   consdata->normalized = TRUE;

   debugMessage("normalized constraint: ");
   debug(consdataPrint(scip, consdata, NULL));

   return SCIP_OKAY;
}

/** replaces multiple occurrences of a variable by a single coefficient */
static
RETCODE mergeMultiples(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< linear constraint */
   )
{
   CONSDATA* consdata;
   VAR* var;
   Real valsum;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* sort the constraint */
   CHECK_OKAY( consdataSort(scip, consdata) );
   
   /* go backwards through the constraint looking for multiple occurrences of the same variable;
    * backward direction is necessary, since consdataDelCoefPos() modifies the given position and
    * the subsequent ones
    */
   for( v = consdata->nvars-1; v >= 1; --v )
   {
      var = consdata->vars[v];
      if( consdata->vars[v-1] == var )
      {
         valsum = consdata->vals[v];
         do
         {
            CHECK_OKAY( delCoefPos(scip, cons, v) );
            --v;
            valsum += consdata->vals[v];
         }
         while( v >= 1 && consdata->vars[v-1] == var );

         /* modify the last existing occurrence of the variable */
         assert(consdata->vars[v] == var);
         if( SCIPisZero(scip, valsum) )
         {
            CHECK_OKAY( delCoefPos(scip, cons, v) );
         }
         else
         {
            CHECK_OKAY( chgCoefPos(scip, cons, v, valsum) );
         }
      }
   }

   return SCIP_OKAY;
}

/** replaces all fixed and aggregated variables by their non-fixed counterparts */
static
RETCODE applyFixings(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   Bool*            conschanged         /**< pointer to store TRUE, if changes were made to the constraint */
   )
{
   CONSDATA* consdata;
   VAR* var;
   Real val;
   Real fixedval;
   Real aggrconst;
   Bool cleanup;
   int v;

   assert(conschanged != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   cleanup = FALSE;
   v = 0;
   while( v < consdata->nvars )
   {
      var = consdata->vars[v];
      val = consdata->vals[v];
      assert(SCIPvarIsTransformed(var));

      switch( SCIPvarGetStatus(var) )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         errorMessage("original variable in transformed linear constraint\n");
         return SCIP_INVALIDDATA;

      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      case SCIP_VARSTATUS_MULTAGGR:
         ++v;
         break;

      case SCIP_VARSTATUS_FIXED:
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)));
         fixedval = SCIPvarGetLbGlobal(var);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            CHECK_OKAY( chgLhs(scip, cons, consdata->lhs - val * fixedval) );
         }
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            CHECK_OKAY( chgRhs(scip, cons, consdata->rhs - val * fixedval) );
         }
         CHECK_OKAY( delCoefPos(scip, cons, v) );
         *conschanged = TRUE;
         break;

      case SCIP_VARSTATUS_AGGREGATED:
         CHECK_OKAY( addCoef(scip, cons, SCIPvarGetAggrVar(var), val * SCIPvarGetAggrScalar(var)) );
         aggrconst = SCIPvarGetAggrConstant(var);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            CHECK_OKAY( chgLhs(scip, cons, consdata->lhs - val * aggrconst) );
         }
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            CHECK_OKAY( chgRhs(scip, cons, consdata->rhs - val * aggrconst) );
         }
         CHECK_OKAY( delCoefPos(scip, cons, v) );
         *conschanged = TRUE;
         cleanup = TRUE;
         break;

      case SCIP_VARSTATUS_NEGATED:
         CHECK_OKAY( addCoef(scip, cons, SCIPvarGetNegationVar(var), -val) );
         aggrconst = SCIPvarGetNegationConstant(var);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            CHECK_OKAY( chgLhs(scip, cons, consdata->lhs - val * aggrconst) );
         }
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            CHECK_OKAY( chgRhs(scip, cons, consdata->rhs - val * aggrconst) );
         }
         CHECK_OKAY( delCoefPos(scip, cons, v) );
         *conschanged = TRUE;
         cleanup = TRUE;
         break;

      default:
         errorMessage("unknown variable status\n");
         abort();
      }
   }

   debugMessage("after fixings: ");
   debug(consdataPrint(scip, consdata, NULL));

   /* if aggregated variables have been replaced, multiple entries of the same variable are possible and we have
    * to clean up the constraint
    */
   if( cleanup )
   {
      CHECK_OKAY( mergeMultiples(scip, cons) );
      
      debugMessage("after merging: ");
      debug(consdataPrint(scip, consdata, NULL));
   }

   return SCIP_OKAY;
}

/** tightens bounds of variables in constraint due to activity bounds */
static
RETCODE tightenBounds(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*             nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   CONSDATA* consdata;
   VAR** vars;
   Real* vals;
   int nvars;

   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   *cutoff = FALSE;

   /* we cannot tighten variables' bounds, if the constraint may be not complete */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   if( nvars > 0 )
   {
      int lastnchgbds;
      int lastsuccess;
      int v;
   
      vars = consdata->vars;
      vals = consdata->vals;
      assert(vars != NULL);
      assert(vals != NULL);
      lastsuccess = 0;
      v = 0;
      do
      {
         assert(0 <= v && v < nvars);
         lastnchgbds = *nchgbds;
         CHECK_OKAY( consdataTightenVarBounds(scip, consdata, vars[v], vals[v], cutoff, nchgbds) );
         if( *nchgbds > lastnchgbds )
            lastsuccess = v;
         v++;
         if( v == nvars )
            v = 0;
      }
      while( v != lastsuccess && !(*cutoff) );
   }

   return SCIP_OKAY;
}

/** checks linear constraint for feasibility of given solution or current solution */
static
RETCODE checkCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   SOL*             sol,                /**< solution to be checked, or NULL for current solution */
   Bool             checklprows,        /**< has linear constraint to be checked, if it is already in current LP? */
   Real*            violation,          /**< pointer to store the constraint's violation, or NULL */
   Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   CONSDATA* consdata;
   Real feasibility;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   debugMessage("checking linear constraint <%s>\n", SCIPconsGetName(cons));
   debug(consdataPrint(scip, consdata, NULL));

   *violated = FALSE;

   if( consdata->row != NULL )
   {
      if( !checklprows && SCIProwIsInLP(consdata->row) )
         return SCIP_OKAY;
      else if( sol == NULL && !SCIPhasActNodeLP(scip) )
         feasibility = consdataGetPseudoFeasibility(scip, consdata);
      else
         feasibility = SCIPgetRowSolFeasibility(scip, consdata->row, sol);
   }
   else
      feasibility = consdataGetFeasibility(scip, consdata, sol);
   
   debugMessage("  consdata feasibility=%g (lhs=%g, rhs=%g, row=%p, checklprows=%d, rowinlp=%d, sol=%p, hasactnodelp=%d)\n",
      feasibility, consdata->lhs, consdata->rhs, consdata->row, checklprows,
      consdata->row == NULL ? 0 : SCIProwIsInLP(consdata->row), sol, SCIPhasActNodeLP(scip));

   if( SCIPisFeasible(scip, feasibility) )
   {
      *violated = FALSE;
      CHECK_OKAY( SCIPincConsAge(scip, cons) );
   }
   else
   {
      *violated = TRUE;
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
   }

   if( violation != NULL )
      *violation = -feasibility;
   
   return SCIP_OKAY;
}

/** creates an LP row in a linear constraint data */
static
RETCODE createRow(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< linear constraint */
   )
{
   CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), consdata->lhs, consdata->rhs,
                  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   
   CHECK_OKAY( SCIPaddVarsToRow(scip, consdata->row, consdata->nvars, consdata->vars, consdata->vals) );

   return SCIP_OKAY;
}

/** adds linear constraint as cut to the LP */
static
RETCODE addCut(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   Real             violation           /**< absolute violation of the constraint */
   )
{
   CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   if( consdata->row == NULL )
   {
      /* convert consdata object into LP row */
      CHECK_OKAY( createRow(scip, cons) );
   }
   assert(consdata->row != NULL);
   assert(!SCIProwIsInLP(consdata->row));
   
   /* insert LP row as cut */
   CHECK_OKAY( SCIPaddCut(scip, consdata->row, 
                  violation/SCIProwGetNorm(consdata->row)/(SCIProwGetNNonz(consdata->row)+1)) );

   return SCIP_OKAY;
}

/** separates linear constraint: adds linear constraint as cut, if violated by current LP solution */
static
RETCODE separateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   RESULT*          result              /**< pointer to store result of separation */
   )
{
   Real violation;
   Bool violated;

   assert(cons != NULL);
   assert(result != NULL);

   CHECK_OKAY( checkCons(scip, cons, NULL, FALSE, &violation, &violated) );

   if( violated )
   {
      /* insert LP row as cut */
      CHECK_OKAY( addCut(scip, cons, violation) );
      *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}

/** propagation method for linear constraints */
static
RETCODE propagateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   Bool             tightenbounds,      /**< should the variable's bounds be tightened? */
   Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*             nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   CONSDATA* consdata;
   Real minactivity;
   Real maxactivity;

   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   debugMessage("propagating linear constraint <%s>\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *cutoff = FALSE;

   /* check, if constraint is already propagated */
   if( consdata->propagated )
      return SCIP_OKAY;

   /* we can only infer activity bounds of the linear constraint, if it is not modifiable */
   if( !SCIPconsIsModifiable(cons) )
   {
      /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
      CHECK_OKAY( SCIPincConsAge(scip, cons) );

      /* tighten the variable's bounds */
      if( tightenbounds )
      {
         int oldnchgbds;

         oldnchgbds = *nchgbds;
         CHECK_OKAY( tightenBounds(scip, cons, cutoff, nchgbds) );
         if( *nchgbds > oldnchgbds )
         {
            CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         }            
      }
      
      /* check constraint for infeasibility and redundancy */
      consdataGetActivityBounds(scip, consdata, &minactivity, &maxactivity);
      
      if( SCIPisFeasGT(scip, minactivity, consdata->rhs) || SCIPisFeasLT(scip, maxactivity, consdata->lhs) )
      {
         debugMessage("linear constraint <%s> is infeasible: activitybounds=[%g,%g], sides=[%g,%g]\n",
            SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *cutoff = TRUE;
      }
      else if( SCIPisGE(scip, minactivity, consdata->lhs) && SCIPisLE(scip, maxactivity, consdata->rhs) )
      {
         debugMessage("linear constraint <%s> is redundant: activitybounds=[%g,%g], sides=[%g,%g]\n",
            SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      }
   }

   /* mark constraint to be propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
DECL_CONSFREE(consFreeLinear)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitLinear NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitLinear NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreLinear NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
DECL_CONSEXITPRE(consExitpreLinear)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int c;

   /* delete all linear constraints that were upgraded to a more specific constraint type */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->upgraded )
      {
         CHECK_OKAY( SCIPdelCons(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolLinear NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
DECL_CONSEXITSOL(consExitsolLinear)
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
DECL_CONSDELETE(consDeleteLinear)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);
   
   /* free linear constraint */
   CHECK_OKAY( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransLinear)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   /*debugMessage("Trans method of linear constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPstage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* create linear constraint data for target constraint */
   CHECK_OKAY( consdataCreateTransformed(scip, &targetdata, conshdlrdata->eventhdlr,
                  sourcedata->nvars, sourcedata->vars, sourcedata->vals, sourcedata->lhs, sourcedata->rhs) );

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                  SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                  SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                  SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpLinear)
{  /*lint --e{715}*/
   int c;

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

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
DECL_CONSSEPA(consSepaLinear)
{  /*lint --e{715}*/
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Sepa method of linear constraints\n");*/

   *result = SCIP_DIDNOTFIND;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      /*debugMessage("separating linear constraint <%s>\n", SCIPconsGetName(conss[c]));*/
      CHECK_OKAY( separateCons(scip, conss[c], result) );
   }

   /* combine linear constraints to get more cuts */
   /**@todo further cuts of linear constraints */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpLinear)
{  /*lint --e{715}*/
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Enfolp method of linear constraints\n");*/

   /* check for violated constraints
    * LP is processed at current node -> we can add violated linear constraints to the LP
    */

   *result = SCIP_FEASIBLE;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], result) );
   }

   /* check all obsolete linear constraints for feasibility */
   for( c = nusefulconss; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      CHECK_OKAY( separateCons(scip, conss[c], result) );
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsLinear)
{  /*lint --e{715}*/
   Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Enfops method of linear constraints\n");*/

   /* if the solution is infeasible anyway due to objective value, skip the enforcement */
   if( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* check all linear constraints for feasibility */
   violated = FALSE;
   for( c = 0; c < nconss && !violated; ++c )
   {
      CHECK_OKAY( checkCons(scip, conss[c], NULL, TRUE, NULL, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckLinear)
{  /*lint --e{715}*/
   Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Check method of linear constraints\n");*/

   /* check all linear constraints for feasibility */
   violated = FALSE;
   for( c = 0; c < nconss && !violated; ++c )
   {
      CHECK_OKAY( checkCons(scip, conss[c], sol, checklprows, NULL, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
DECL_CONSPROP(consPropLinear)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   Bool tightenbounds;
   Bool cutoff;
   int nchgbds;
   int propfreq;
   int depth;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Prop method of linear constraints\n");*/

   /* check, if we want to tighten variable's bounds */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   propfreq = SCIPconshdlrGetPropFreq(conshdlr);
   depth = SCIPgetDepth(scip);
   tightenbounds = (conshdlrdata->tightenboundsfreq == 0 && depth == 0)
      || (conshdlrdata->tightenboundsfreq >= 1 && (depth % (propfreq * conshdlrdata->tightenboundsfreq) == 0));

   cutoff = FALSE;
   nchgbds = 0;

   /* process useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      CHECK_OKAY( propagateCons(scip, conss[c], tightenbounds, &cutoff, &nchgbds) );
   }

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}




/*
 * Presolving
 */

/* tightens left and right hand side of constraint due to integrality */
static
RETCODE tightenSides(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int*             nchgsides,          /**< pointer to count number of side changes */
   Bool*            conschanged         /**< pointer to store TRUE, if changes were made to the constraint */
   )
{
   CONSDATA* consdata;
   Bool integral;
   int i;

   assert(nchgsides != NULL);
   assert(conschanged != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPisIntegral(scip, consdata->lhs) || !SCIPisIntegral(scip, consdata->rhs) )
   {
      integral = TRUE;
      for( i = 0; i < consdata->nvars && integral; ++i )
      {
         integral = integral
            && SCIPisIntegral(scip, consdata->vals[i])
            && (SCIPvarGetType(consdata->vars[i]) != SCIP_VARTYPE_CONTINUOUS);
      }
      if( integral )
      {
         debugMessage("linear constraint <%s>: make sides integral: sides=[%g,%g]\n",
            SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
         if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisIntegral(scip, consdata->lhs) )
         {
            CHECK_OKAY( chgLhs(scip, cons, SCIPceil(scip, consdata->lhs)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
            *conschanged = TRUE;
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisIntegral(scip, consdata->rhs) )
         {
            CHECK_OKAY( chgRhs(scip, cons, SCIPfloor(scip, consdata->rhs)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
            *conschanged = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}

/** tightens coefficients of binary, integer, and implicit integer variables due to activity bounds in presolving:
 *  given an inequality  lhs <= a*x + ai*xi <= rhs, with a non-continouos variable  li <= xi <= ui
 *  let minact := min{a*x + ai*xi}, maxact := max{a*x + ai*xi}
 *  (i) ai >= 0:
 *      if  minact + ai >= lhs  and  maxact - ai <= rhs:
 *       - a deviation from the lower/upper bound of xi would make the left/right hand side redundant
 *       - ai, lhs and rhs can be changed to have the same redundancy effect and the same results for
 *         xi fixed to its bounds, but with a reduced ai and tightened sides to tighten the LP relaxation
 *       - change coefficients:
 *           ai'  := max(lhs - minact, maxact - rhs)
 *           lhs' := lhs - (ai - ai')*li
 *           rhs' := rhs - (ai - ai')*ui
 * (ii) ai < 0:
 *      if  minact - ai >= lhs  and  maxact + ai <= rhs:
 *       - a deviation from the upper/lower bound of xi would make the left/right hand side redundant
 *       - ai, lhs and rhs can be changed to have the same redundancy effect and the same results for
 *         xi fixed to its bounds, but with a reduced ai and tightened sides to tighten the LP relaxation
 *       - change coefficients:
 *           ai'  := min(rhs - maxact, minact - lhs)
 *           lhs' := lhs - (ai - ai')*ui
 *           rhs' := rhs - (ai - ai')*li
 */
static
RETCODE consdataTightenCoefs(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int*             nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*             nchgsides,          /**< pointer to count number of side changes */
   Bool*            conschanged         /**< pointer to store TRUE, if changes were made to the constraint */
   )
{
   CONSDATA* consdata;
   VAR* var;
   Real minactivity;
   Real maxactivity;
   Real val;
   Real newval;
   Real newlhs;
   Real newrhs;
   Real lb;
   Real ub;
   int i;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(conschanged != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get the minimal and maximal activity of the constraint */
   consdataGetActivityBounds(scip, consdata, &minactivity, &maxactivity);

   /* try to tighten each coefficient */
   for( i = 0; i < consdata->nvars; ++i )
   {
      var = consdata->vars[i];

      /* ignore continouos variables */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      /* get coefficient and variable's bounds */
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      val = consdata->vals[i];
      assert(!SCIPisZero(scip, val));

      /* check sign of coefficient */
      if( val >= 0.0 )
      {
         /* check, if a deviation from lower/upper bound would make lhs/rhs redundant */
         if( SCIPisGE(scip, minactivity + val, consdata->lhs) && SCIPisLE(scip, maxactivity - val, consdata->rhs) )
         {
            /* change coefficients:
             *   ai'  := max(lhs - minact, maxact - rhs)
             *   lhs' := lhs - (ai - ai')*li
             *   rhs' := rhs - (ai - ai')*ui
             */
            newval = MAX(consdata->lhs - minactivity, maxactivity - consdata->rhs);
            newlhs = consdata->lhs - (val - newval)*lb;
            newrhs = consdata->rhs - (val - newval)*ub;
            if( !SCIPisEQ(scip, newval, val) )
            {
               debugMessage("linear constraint <%s>: change coefficient %+g<%s> to %+g<%s>, act=[%g,%g], side=[%g,%g]\n",
                  SCIPconsGetName(cons), val, SCIPvarGetName(var), newval, SCIPvarGetName(var), 
                  minactivity, maxactivity, consdata->lhs, consdata->rhs);

               /* update the coefficient and the activity bounds */
               consdata->vals[i] = newval;
               consdataUpdateChgCoef(scip, consdata, var, val, newval);
               (*nchgcoefs)++;
               *conschanged = TRUE;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, &minactivity, &maxactivity);
            }
            if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisEQ(scip, newlhs, consdata->lhs) )
            {
               debugMessage("linear constraint <%s>: change lhs %g to %g\n", SCIPconsGetName(cons), consdata->lhs, newlhs);

               CHECK_OKAY( chgLhs(scip, cons, newlhs) );
               (*nchgsides)++;
               *conschanged = TRUE;
               assert(SCIPisEQ(scip, consdata->lhs, newlhs));
            }
            if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisEQ(scip, newrhs, consdata->rhs) )
            {
               debugMessage("linear constraint <%s>: change rhs %g to %g\n", SCIPconsGetName(cons), consdata->rhs, newrhs);

               CHECK_OKAY( chgRhs(scip, cons, newrhs) );
               (*nchgsides)++;
               *conschanged = TRUE;
               assert(SCIPisEQ(scip, consdata->rhs, newrhs));
            }
         }
      }
      else
      {
         /* check, if a deviation from lower/upper bound would make lhs/rhs redundant */
         if( SCIPisGE(scip, minactivity - val, consdata->lhs) && SCIPisLE(scip, maxactivity + val, consdata->rhs) )
         {
            /* change coefficients:
             *   ai'  := min(rhs - maxact, minact - lhs)
             *   lhs' := lhs - (ai - ai')*ui
             *   rhs' := rhs - (ai - ai')*li
             */
            newval = MIN(consdata->rhs - maxactivity, minactivity - consdata->lhs);
            newlhs = consdata->lhs - (val - newval)*ub;
            newrhs = consdata->rhs - (val - newval)*lb;
            if( !SCIPisEQ(scip, newval, val) )
            {
               debugMessage("linear constraint <%s>: change coefficient %+g<%s> to %+g<%s>, act=[%g,%g], side=[%g,%g]\n",
                  SCIPconsGetName(cons), val, SCIPvarGetName(var), newval, SCIPvarGetName(var), 
                  minactivity, maxactivity, consdata->lhs, consdata->rhs);

               /* update the coefficient and the activity bounds */
               consdata->vals[i] = newval;
               consdataUpdateChgCoef(scip, consdata, var, val, newval);
               (*nchgcoefs)++;
               *conschanged = TRUE;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, &minactivity, &maxactivity);
            }
            if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisEQ(scip, newlhs, consdata->lhs) )
            {
               debugMessage("linear constraint <%s>: change lhs %g to %g\n", SCIPconsGetName(cons), consdata->lhs, newlhs);

               CHECK_OKAY( chgLhs(scip, cons, newlhs) );
               (*nchgsides)++;
               *conschanged = TRUE;
               assert(SCIPisEQ(scip, consdata->lhs, newlhs));
            }
            if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisEQ(scip, newrhs, consdata->rhs) )
            {
               debugMessage("linear constraint <%s>: change rhs %g to %g\n", SCIPconsGetName(cons), consdata->rhs, newrhs);

               CHECK_OKAY( chgRhs(scip, cons, newrhs) );
               (*nchgsides)++;
               *conschanged = TRUE;
               assert(SCIPisEQ(scip, consdata->rhs, newrhs));
            }
         }
      }
   }
   
   return SCIP_OKAY;
}

/* processes equality with only one variable by fixing the variable and deleting the constraint */
static
RETCODE convertUnaryEquality(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int*             nfixedvars,         /**< pointer to count number of fixed variables */
   int*             ndelconss,          /**< pointer to count number of deleted constraints */
   RESULT*          result,             /**< pointer to store result for successful conversions */
   Bool*            consdeleted         /**< pointer to store TRUE, if constraint was deleted */
   )
{
   CONSDATA* consdata;
   VAR* var;
   Real val;
   Real fixval;
   Bool infeasible;
   Bool fixed;

   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);
   assert(result != NULL);
   assert(consdeleted != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 1);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   /* calculate the value to fix the variable to */
   var = consdata->vars[0];
   val = consdata->vals[0];
   assert(!SCIPisZero(scip, val));
   fixval = consdata->rhs/val;

   debugMessage("linear equality <%s>: fix <%s> == %g\n",
      SCIPconsGetName(cons), SCIPvarGetName(var), fixval);
   
   /* fix variable */
   CHECK_OKAY( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
   if( infeasible )
   {
      debugMessage(" -> infeasible fixing\n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }
   if( fixed )
      (*nfixedvars)++;

   /* disable constraint */
   CHECK_OKAY( SCIPdelCons(scip, cons) );
   if( !consdata->upgraded )
      (*ndelconss)++;
   *result = SCIP_SUCCESS;
   *consdeleted = TRUE;

   return SCIP_OKAY;
}

/* processes equality with exactly two variables by aggregating one of the variables and deleting the constraint */
static
RETCODE convertBinaryEquality(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int*             naggrvars,          /**< pointer to count number of aggregated variables */
   int*             ndelconss,          /**< pointer to count number of deleted constraints */
   RESULT*          result,             /**< pointer to store result for successful conversions */
   Bool*            consdeleted         /**< pointer to store TRUE, if constraint was deleted */
   )
{
   CONSDATA* consdata;
   Bool infeasible;
   Bool redundant;
   Bool aggregated;

   assert(naggrvars != NULL);
   assert(ndelconss != NULL);
   assert(result != NULL);
   assert(consdeleted != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 2);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   debugMessage("linear constraint <%s>: aggregate %g<%s> + %g<%s> == %g\n",
      SCIPconsGetName(cons), consdata->vals[0], SCIPvarGetName(consdata->vars[0]), 
      consdata->vals[1], SCIPvarGetName(consdata->vars[1]), consdata->rhs);
   
   /* aggregate the equality */
   CHECK_OKAY( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], consdata->vals[0], consdata->vals[1],
                  consdata->rhs, &infeasible, &redundant, &aggregated) );

   /* check for infeasibility of aggregation */
   if( infeasible )
   {
      debugMessage(" -> infeasible aggregation\n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }
   
   /* count the aggregation */
   if( aggregated )
   {
      (*naggrvars)++;
      *result = SCIP_SUCCESS;
   }

   /* delete the constraint, if it is redundant */
   if( redundant )
   {
      CHECK_OKAY( SCIPdelCons(scip, cons) );
      
      if( !consdata->upgraded )
         (*ndelconss)++;
      *result = SCIP_SUCCESS;
      *consdeleted = TRUE;
   }

   return SCIP_OKAY;
}

/* processes equality with more than two variables by multi-aggregating one of the variables and converting the equality
 * into an inequality
 */
static
RETCODE convertLongEquality(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int*             naggrvars,          /**< pointer to count number of aggregated variables */
   RESULT*          result,             /**< pointer to store result for successful conversions */
   Bool*            conschanged         /**< pointer to store TRUE, if changes were made to the constraint */
   )
{
   CONSDATA* consdata;
   VAR** vars;
   Real* vals;
   VAR* var;
   Real val;
   VARTYPE bestslacktype;
   VARTYPE slacktype;
   Real bestslackdomrng;
   Real slackdomrng;
   Bool integral;
   int bestslackpos;
   int v;

   assert(naggrvars != NULL);
   assert(result != NULL);
   assert(conschanged != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 2);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   /* look for a slack variable s to convert a*x + s == b into lhs <= a*x <= rhs */
   vars = consdata->vars;
   vals = consdata->vals;
   bestslackpos = -1;
   bestslacktype = SCIP_VARTYPE_BINARY;
   bestslackdomrng = 0.0;
   integral = TRUE;
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(vars != NULL);
      assert(vals != NULL);
      var = vars[v];
      val = vals[v];

      assert(SCIPvarGetNLocksDown(var) >= 1); /* because variable is locked in this equality */
      assert(SCIPvarGetNLocksUp(var) >= 1);

      slacktype = SCIPvarGetType(var);
      integral = integral && (slacktype != SCIP_VARTYPE_CONTINUOUS) && SCIPisIntegral(scip, val);

      /* check, if variable is already fixed or aggregated */
      if( !SCIPvarIsActive(var) )
         continue;

      /* check, if variable is used in other constraints than this one */
      if( SCIPvarGetNLocksDown(var) > 1 || SCIPvarGetNLocksUp(var) > 1 )
         continue;

      /* check, if variable can be used as a slack variable */
      if( slacktype == SCIP_VARTYPE_CONTINUOUS
         || slacktype == SCIP_VARTYPE_IMPLINT
         || (integral && SCIPisEQ(scip, ABS(val), 1.0))
          )
      {
         slackdomrng = SCIPvarGetUbGlobal(var) - SCIPvarGetLbGlobal(var);
         if( bestslackpos == -1
            || slacktype > bestslacktype
            || (slacktype == bestslacktype && slackdomrng > bestslackdomrng) )
         {
            bestslackpos = v;
            bestslacktype = slacktype;
            bestslackdomrng = slackdomrng;
         }
      }
   }

   /* if all coefficients and variables are integral, the right hand side must also be integral */
   if( integral && !SCIPisIntegral(scip, consdata->rhs) )
   {
      debugMessage("linear equality <%s> is integer infeasible:", SCIPconsGetName(cons));
      debug(consdataPrint(scip, consdata, NULL));
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* if the slack variable is of integer type, and the constraint itself may take fractional values,
    * we cannot aggregate the variable, because the integrality condition would get lost
    */
   if( bestslackpos >= 0
      && (bestslacktype == SCIP_VARTYPE_CONTINUOUS || bestslacktype == SCIP_VARTYPE_IMPLINT || integral) )
   {
      VAR* slackvar;
      Real* scalars;
      Real slackcoef;
      Real slackvarlb;
      Real slackvarub;
      Real aggrconst;
      Real newlhs;
      Real newrhs;
      Bool infeasible;

      /* we found a slack variable that only occurs in this equality:
       *   a_1*x_1 + ... + a_k*x_k + a'*s == rhs  ->  s == rhs - a_1/a'*x_1 - ... - a_k/a'*x_k
       */

      /* convert equality into inequality by deleting the slack variable:
       *  x + a*s == b, l <= s <= u   ->  b - a*u <= x <= b - a*l
       */
      slackvar = vars[bestslackpos];
      slackcoef = vals[bestslackpos];
      assert(!SCIPisZero(scip, slackcoef));
      aggrconst = consdata->rhs/slackcoef;
      slackvarlb = SCIPvarGetLbGlobal(slackvar);
      slackvarub = SCIPvarGetUbGlobal(slackvar);
      if( slackcoef > 0.0 )
      {
         if( SCIPisInfinity(scip, -slackvarlb) )
            newrhs = SCIPinfinity(scip);
         else
            newrhs = consdata->rhs - slackcoef * slackvarlb;
         if( SCIPisInfinity(scip, slackvarub) )
            newlhs = -SCIPinfinity(scip);
         else
            newlhs = consdata->lhs - slackcoef * slackvarub;
      }
      else
      {
         if( SCIPisInfinity(scip, -slackvarlb) )
            newlhs = -SCIPinfinity(scip);
         else
            newlhs = consdata->rhs - slackcoef * slackvarlb;
         if( SCIPisInfinity(scip, slackvarub) )
            newrhs = SCIPinfinity(scip);
         else
            newrhs = consdata->lhs - slackcoef * slackvarub;
      }
      assert(SCIPisLE(scip, newlhs, newrhs));
      CHECK_OKAY( chgLhs(scip, cons, newlhs) );
      CHECK_OKAY( chgRhs(scip, cons, newrhs) );
      CHECK_OKAY( delCoefPos(scip, cons, bestslackpos) );

      /* allocate temporary memory */
      CHECK_OKAY( SCIPallocBufferArray(scip, &scalars, consdata->nvars) );

      /* set up the multi-aggregation */
      debugMessage("linear constraint <%s>: multi-aggregate <%s> ==", SCIPconsGetName(cons), SCIPvarGetName(slackvar));
      for( v = 0; v < consdata->nvars; ++v )
      {
         scalars[v] = -consdata->vals[v]/slackcoef;
         debug(printf(" %+g<%s>", scalars[v], SCIPvarGetName(consdata->vars[v])));
      }
      debug(printf(" %+g, bounds of <%s>: [%g,%g]\n", aggrconst, SCIPvarGetName(slackvar), slackvarlb, slackvarub));

      /* perform the multi-aggregation */
      CHECK_OKAY( SCIPmultiaggregateVar(scip, slackvar, consdata->nvars, consdata->vars, scalars, aggrconst,
                     &infeasible) );

      /* free temporary memory */
      CHECK_OKAY( SCIPfreeBufferArray(scip, &scalars) );

      /* check for infeasible aggregation */
      if( infeasible )
      {
         debugMessage("linear constraint <%s>: infeasible multi-aggregation\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      (*naggrvars)++;
      *result = SCIP_SUCCESS;
      *conschanged = TRUE;
   }

   return SCIP_OKAY;
}

/* converts special equalities */
static
RETCODE convertEquality(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int*             nfixedvars,         /**< pointer to count number of fixed variables */
   int*             naggrvars,          /**< pointer to count number of aggregated variables */
   int*             ndelconss,          /**< pointer to count number of deleted constraints */
   RESULT*          result,             /**< pointer to store result for successful conversions */
   Bool*            conschanged,        /**< pointer to store TRUE, if changes were made to the constraint */
   Bool*            consdeleted         /**< pointer to store TRUE, if constraint was deleted */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* do nothing on inequalities */
   if( !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      return SCIP_OKAY;

   /* depending on the number of variables, call a special conversion method */
   if( consdata->nvars == 1 )
   {
      /* fix variable */
      CHECK_OKAY( convertUnaryEquality(scip, cons, nfixedvars, ndelconss, result, consdeleted) );
   }
   else if( consdata->nvars == 2 )
   {
      /* aggregate one of the variables */
      CHECK_OKAY( convertBinaryEquality(scip, cons, naggrvars, ndelconss, result, consdeleted) );
   }
   else
   {
      /* try to multi-aggregate one of the variables */
      CHECK_OKAY( convertLongEquality(scip, cons, naggrvars, result, conschanged) );
   }

   return SCIP_OKAY;
}

/** converts all variables with fixed domain into FIXED variables */
static
RETCODE fixVariables(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< linear constraint */
   int*             nfixedvars,         /**< pointer to count the total number of fixed variables */
   RESULT*          result,             /**< pointer to store the result of the variable fixing */
   Bool*            conschanged         /**< pointer to store TRUE, if changes were made to the constraint */
   )
{
   CONSDATA* consdata;
   VAR* var;
   VARSTATUS varstatus;
   Real lb;
   Real ub;
   Bool fixed;
   Bool fixingsfound;
   Bool infeasible;
   int v;

   assert(nfixedvars != NULL);
   assert(result != NULL);
   assert(*result != SCIP_CUTOFF);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   fixingsfound = FALSE;
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars != NULL);
      var = consdata->vars[v];
      varstatus = SCIPvarGetStatus(var);

      if( varstatus != SCIP_VARSTATUS_FIXED )
      {
         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);
         if( SCIPisEQ(scip, lb, ub) )
         {
            debugMessage("converting variable <%s> with fixed bounds [%g,%g] into fixed variable\n",
               SCIPvarGetName(var), lb, ub);
            CHECK_OKAY( SCIPfixVar(scip, var, lb, &infeasible, &fixed) );
            if( infeasible )
            {
               debugMessage(" -> infeasible fixing\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            if( fixed )
            {
               (*nfixedvars)++;
               *result = SCIP_SUCCESS;
               fixingsfound = TRUE;
            }
         }
      }
   }

   if( fixingsfound )
   {
      CHECK_OKAY( applyFixings(scip, cons, conschanged) );
      assert(*conschanged);
   }
   
   return SCIP_OKAY;
}

#define BINWEIGHT  1
#define INTWEIGHT  4
#define CONTWEIGHT 8

/** gets weight for variable in a "weighted number of variables" sum */
static
int getVarWeight(
   VAR*             var                 /**< variable to get weight for */
   )
{
   switch( SCIPvarGetType(var) )
   {
   case SCIP_VARTYPE_BINARY:
      return BINWEIGHT;
   case SCIP_VARTYPE_INTEGER:
   case SCIP_VARTYPE_IMPLINT:
      return INTWEIGHT;
   case SCIP_VARTYPE_CONTINUOUS:
      return CONTWEIGHT;
   default:
      errorMessage("invalid variable type\n");
      abort();
   }
}

/* tries to aggregate an (in)equality and an equality in order to decrease the number of variables in the (in)equality:
 *   cons0 := a * cons0 + b * cons1, 
 * where a = val1[v] and b = -val0[v] for common variable v which removes most variable weight;
 * for numerical stability, we will only accept integral a and b;
 * the variable weight is a weighted sum over all included variables, where each binary variable weighs BINWEIGHT, 
 * each integer or implicit integer variable weighs INTWEIGHT and each continuous variable weighs CONTWEIGHT
 */
static 
RETCODE aggregateConstraints(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons0,              /**< (in)equality to modify */
   CONS*            cons1,              /**< equality to use for aggregation of cons0 */
   int*             commonidx0,         /**< array with indices of variables in cons0, that appear also in cons1 */
   int*             commonidx1,         /**< array with indices of variables in cons1, that appear also in cons0 */
   int*             diffidx0minus1,     /**< array with indices of variables in cons0, that don't appear in cons1 */
   int*             diffidx1minus0,     /**< array with indices of variables in cons1, that don't appear in cons0 */
   int              nvarscommon,        /**< number of variables, that appear in both constraints */
   int              commonidxweight,    /**< variable weight sum of common variables */
   int              diffidx0minus1weight, /**< variable weight sum of variables in cons0, that don't appear in cons1 */
   int              diffidx1minus0weight, /**< variable weight sum of variables in cons1, that don't appear in cons0 */
   Real             maxaggrnormscale,   /**< maximal allowed relative gain in maximum norm for constraint aggregation */
   int*             nupgdconss,         /**< pointer to count the number of upgraded constraints */
   int*             nchgcoefs,          /**< pointer to count the number of changed coefficients */
   RESULT*          result,             /**< pointer to store the result of the aggregation */
   Bool*            aggregated          /**< pointer to store whether an aggregation was made */
   )
{
   CONSDATA* consdata0;
   CONSDATA* consdata1;
   Real a;
   Real b;
   Real aggrcoef;
   Real scalarsum;
   Real bestscalarsum;
   Bool betterscalarsum;
   int varweight;
   int nvars;
   int bestvarweight;
   int bestnvars;
   int bestv;
   int v;
   int i;

   assert(commonidx0 != NULL);
   assert(commonidx1 != NULL);
   assert(diffidx0minus1 != NULL);
   assert(diffidx1minus0 != NULL);
   assert(nvarscommon >= 1);
   assert(commonidxweight >= nvarscommon);
   assert(nchgcoefs != NULL);
   assert(result != NULL);
   assert(aggregated != NULL);

   assert(SCIPconsIsActive(cons0));
   assert(SCIPconsIsActive(cons1));
   
   debugMessage("try aggregation of <%s> and <%s>\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));

   /* cons0 is an (in)equality */
   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);
   assert(SCIPisLE(scip, consdata0->lhs, consdata0->rhs));
   assert(diffidx0minus1weight >= consdata0->nvars - nvarscommon);
   
   /* cons1 is an equality */
   consdata1 = SCIPconsGetData(cons1);
   assert(consdata1 != NULL);
   assert(consdata1->nvars >= 1);
   assert(SCIPisEQ(scip, consdata1->lhs, consdata1->rhs));
   assert(diffidx1minus0weight >= consdata1->nvars - nvarscommon);

   *aggregated = FALSE;

   /* search for the best common variable such that
    *   val1[var] * consdata0 - val0[var] * consdata1
    * has least weighted number of variables
    */
   bestvarweight = commonidxweight + diffidx0minus1weight;
   bestnvars = consdata0->nvars;
   bestv = -1;
   bestscalarsum = 0.0;
   for( v = 0; v < nvarscommon; ++v )
   {
      assert(consdata0->vars[commonidx0[v]] == consdata1->vars[commonidx1[v]]);
      a = consdata1->vals[commonidx1[v]];
      b = -consdata0->vals[commonidx0[v]];

      /* only try aggregation, if coefficients are integral (numerical stability) */
      if( SCIPisIntegral(scip, a) && SCIPisIntegral(scip, b) )
      {
         /* count the number of variables in the potential new constraint  a * consdata0 + b * consdata1 */
         varweight = diffidx0minus1weight + diffidx1minus0weight;
         nvars = consdata0->nvars + consdata1->nvars - 2*nvarscommon;
         scalarsum = ABS(a) + ABS(b);
         betterscalarsum = (scalarsum < bestscalarsum);
         for( i = 0; i < nvarscommon
                 && (varweight < bestvarweight || (varweight == bestvarweight && betterscalarsum)); ++i )
         {
            aggrcoef = a * consdata0->vals[commonidx0[i]] + b * consdata1->vals[commonidx1[i]];
            if( !SCIPisZero(scip, aggrcoef) )
            {
               varweight += getVarWeight(consdata0->vars[commonidx0[i]]);
               nvars++;
            }
         }
         if( varweight < bestvarweight || (varweight == bestvarweight && betterscalarsum) )
         {
            bestv = v;
            bestvarweight = varweight;
            bestnvars = nvars;
            bestscalarsum = scalarsum;
         }
      }
   }

   /* if better aggregation was found, create new constraint and delete old one */
   if( bestv != -1 )
   {
      CONS* newcons;
      CONSDATA* newconsdata;
      VAR** newvars;
      Real* newvals;
      Real newlhs;
      Real newrhs;
      int newnvars;

      /* choose multipliers such that the multiplier for the (in)equality cons0 is positive */
      if( consdata1->vals[commonidx1[bestv]] > 0.0 )
      {
         a = consdata1->vals[commonidx1[bestv]];
         b = -consdata0->vals[commonidx0[bestv]];
      }
      else
      {
         a = -consdata1->vals[commonidx1[bestv]];
         b = consdata0->vals[commonidx0[bestv]];
      }
      assert(SCIPisIntegral(scip, a));
      assert(SCIPisPositive(scip, a));
      assert(SCIPisIntegral(scip, b));
      assert(!SCIPisZero(scip, b));

      debugMessage("aggregate linear constraints <%s> := %g*<%s> + %g*<%s>  ->  nvars: %d -> %d, weight: %d -> %d\n",
         SCIPconsGetName(cons0), a, SCIPconsGetName(cons0), b, SCIPconsGetName(cons1),
         consdata0->nvars, bestnvars, commonidxweight + diffidx0minus1weight, bestvarweight);
      debugMessage("<%s>: ", SCIPconsGetName(cons0));
      debug(consdataPrint(scip, consdata0, NULL));
      debugMessage("<%s>: ", SCIPconsGetName(cons1));
      debug(consdataPrint(scip, consdata1, NULL));

      /* get temporary memory for creating the new linear constraint */
      CHECK_OKAY( SCIPallocBufferArray(scip, &newvars, bestnvars) );
      CHECK_OKAY( SCIPallocBufferArray(scip, &newvals, bestnvars) );

      /* calculate the common coefficients */
      newnvars = 0;
      for( i = 0; i < nvarscommon; ++i )
      {
         assert(0 <= commonidx0[i] && commonidx0[i] < consdata0->nvars);
         assert(0 <= commonidx1[i] && commonidx1[i] < consdata1->nvars);

         aggrcoef = a * consdata0->vals[commonidx0[i]] + b * consdata1->vals[commonidx1[i]];
         if( !SCIPisZero(scip, aggrcoef) )
         {
            assert(newnvars < bestnvars);
            newvars[newnvars] = consdata0->vars[commonidx0[i]];
            newvals[newnvars] = aggrcoef;
            newnvars++;
         }
      }

      /* calculate the coefficients appearing in cons0 but not in cons1 */
      for( i = 0; i < consdata0->nvars - nvarscommon; ++i )
      {
         assert(0 <= diffidx0minus1[i] && diffidx0minus1[i] < consdata0->nvars);

         aggrcoef = a * consdata0->vals[diffidx0minus1[i]];
         assert(!SCIPisZero(scip, aggrcoef));
         assert(newnvars < bestnvars);
         newvars[newnvars] = consdata0->vars[diffidx0minus1[i]];
         newvals[newnvars] = aggrcoef;
         newnvars++;
      }

      /* calculate the coefficients appearing in cons1 but not in cons0 */
      for( i = 0; i < consdata1->nvars - nvarscommon; ++i )
      {
         assert(0 <= diffidx1minus0[i] && diffidx1minus0[i] < consdata1->nvars);

         aggrcoef = b * consdata1->vals[diffidx1minus0[i]];
         assert(!SCIPisZero(scip, aggrcoef));
         assert(newnvars < bestnvars);
         newvars[newnvars] = consdata1->vars[diffidx1minus0[i]];
         newvals[newnvars] = aggrcoef;
         newnvars++;
      }
      assert(newnvars == bestnvars);

      /* calculate the new left and right hand side of the (in)equality */
      assert(!SCIPisInfinity(scip, -consdata1->lhs));
      assert(!SCIPisInfinity(scip, consdata1->rhs));
      if( SCIPisInfinity(scip, -consdata0->lhs) )
         newlhs = -SCIPinfinity(scip);
      else
         newlhs = a * consdata0->lhs + b * consdata1->lhs;
      if( SCIPisInfinity(scip, consdata0->rhs) )
         newrhs = SCIPinfinity(scip);
      else
         newrhs = a * consdata0->rhs + b * consdata1->rhs;

      /* create the new linear constraint */
      CHECK_OKAY( SCIPcreateConsLinear(scip, &newcons, SCIPconsGetName(cons0), newnvars, newvars, newvals, newlhs, newrhs,
                     SCIPconsIsInitial(cons0), SCIPconsIsSeparated(cons0), SCIPconsIsEnforced(cons0), 
                     SCIPconsIsChecked(cons0), SCIPconsIsPropagated(cons0),
                     SCIPconsIsLocal(cons0), SCIPconsIsModifiable(cons0), SCIPconsIsRemoveable(cons0)) );

      newconsdata = SCIPconsGetData(newcons);
      assert(newconsdata != NULL);

      /* copy the upgraded flag from the old cons0 to the new constraint */
      newconsdata->upgraded = consdata0->upgraded;

      /* normalize the new constraint */
      CHECK_OKAY( normalizeCons(scip, newcons) );

      /* check, if we really want to use the new constraint instead of the old one:
       * use the new one, if the maximum norm doesn't grow too much
       */
      if( consdataGetMaxAbsval(SCIPconsGetData(newcons)) <= maxaggrnormscale * consdataGetMaxAbsval(consdata0) )
      {
         CONS* upgdcons;

         debugMessage(" -> aggregated to <%s>: ", SCIPconsGetName(newcons));
         debug(consdataPrint(scip, SCIPconsGetData(newcons), NULL));

         /* update the statistics: we changed all coefficients */
         if( !consdata0->upgraded )
            (*nchgcoefs) += consdata0->nvars + consdata1->nvars - nvarscommon;
         *result = SCIP_SUCCESS;
         *aggregated = TRUE;
         
         /* delete the old constraint, and add the new linear constraint to the problem */
         CHECK_OKAY( SCIPdelCons(scip, cons0) );
         CHECK_OKAY( SCIPaddCons(scip, newcons) );

         /* try to upgrade the new constraint */
         CHECK_OKAY( SCIPupgradeConsLinear(scip, newcons, &upgdcons) );
         if( upgdcons != NULL )
         {
            debugMessage("upgraded linear constraint <%s> to type <%s>:\n", 
               SCIPconsGetName(newcons), SCIPconshdlrGetName(SCIPconsGetHdlr(upgdcons)));
            debug(consdataPrint(scip, newconsdata, NULL));

            /* add the upgraded constraint to the problem */
            CHECK_OKAY( SCIPaddCons(scip, upgdcons) );
            CHECK_OKAY( SCIPreleaseCons(scip, &upgdcons) );
            (*nupgdconss)++;
            
            /* mark the linear constraint being upgraded and to be removed after presolving;
             * don't delete it directly, because it may help to preprocess other linear constraints
             */
            assert(!newconsdata->upgraded);
            newconsdata->upgraded = TRUE;
         }
      }

      /* release the new constraint */
      CHECK_OKAY( SCIPreleaseCons(scip, &newcons) );

      /* free temporary memory */
      CHECK_OKAY( SCIPfreeBufferArray(scip, &newvals) );
      CHECK_OKAY( SCIPfreeBufferArray(scip, &newvars) );
   }

   return SCIP_OKAY;
}

/** compares constraint with all prior constraints for possible redundancy or aggregation,
 *  and removes or changes constraint accordingly
 */
static
RETCODE preprocessConstraintPairs(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           conss,              /**< constraint set */
   int              firstchange,        /**< first constraint that changed since last pair preprocessing round */
   int              chkind,             /**< index of constraint to check against all prior indices upto startind */
   Real             maxaggrnormscale,   /**< maximal allowed relative gain in maximum norm for constraint aggregation */
   int*             nfixedvars,         /**< pointer to count number of fixed variables */
   int*             naggrvars,          /**< pointer to count number of aggregated variables */
   int*             ndelconss,          /**< pointer to count number of deleted constraints */
   int*             nupgdconss,         /**< pointer to count the number of upgraded constraints */
   int*             nchgsides,          /**< pointer to count number of changed left/right hand sides */
   int*             nchgcoefs,          /**< pointer to count number of changed coefficients */
   RESULT*          result              /**< pointer to store result for successful conversions */
   )
{
   CONS* cons0;
   CONS* cons1;
   CONSDATA* consdata0;
   CONSDATA* consdata1;
   VAR* var;
   int* commonidx0;
   int* commonidx1;
   int* diffidx0minus1;
   int* diffidx1minus0;
   Real val0;
   Real val1;
   Bool cons0dominateslhs;
   Bool cons1dominateslhs;
   Bool cons0dominatesrhs;
   Bool cons1dominatesrhs;
   Bool cons0isequality;
   Bool cons1isequality;
   int diffidx1minus0size;
   int nvarscommon;
   int nvars0minus1;
   int nvars1minus0;
   int commonidxweight;
   int diffidx0minus1weight;
   int diffidx1minus0weight;
   int varcmp;
   int c;
   int v0;
   int v1;

   assert(conss != NULL);
   assert(firstchange <= chkind);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);
   assert(nchgsides != NULL);
   assert(result != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(SCIPconsIsActive(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);
   cons0isequality = SCIPisEQ(scip, consdata0->lhs, consdata0->rhs);

   /* sort the constraint */
   CHECK_OKAY( consdataSort(scip, consdata0) );

   /* get temporary memory for indices of common variables */
   CHECK_OKAY( SCIPallocBufferArray(scip, &commonidx0, consdata0->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &commonidx1, consdata0->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &diffidx0minus1, consdata0->nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &diffidx1minus0, consdata0->nvars) );
   diffidx1minus0size = consdata0->nvars;

   /* check constraint against all prior constraints */
   for( c = (consdata0->changed ? 0 : firstchange); c < chkind && *result != SCIP_CUTOFF && SCIPconsIsActive(cons0);
        ++c )
   {
      cons1 = conss[c];
      assert(cons1 != NULL);

      /* ignore inactive constraints */
      if( !SCIPconsIsActive(cons1) )
         continue;

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

      /* if both constraints didn't change since last pair processing, we can ignore the pair */
      if( !consdata0->changed && !consdata1->changed )
         continue;

      /* if both constraints are already upgraded, skip the pair */
      if( consdata0->upgraded && consdata1->upgraded )
         continue;

      assert(consdata1->nvars >= 1);

      /* sort the constraint */
      CHECK_OKAY( consdataSort(scip, consdata1) );

      cons1isequality = SCIPisEQ(scip, consdata1->lhs, consdata1->rhs);
      
      /* make sure, we have enough memory for the index set of V_1 \ V_0 */
      if( consdata1->nvars > diffidx1minus0size )
      {
         CHECK_OKAY( SCIPfreeBufferArray(scip, &diffidx1minus0) );
         CHECK_OKAY( SCIPallocBufferArray(scip, &diffidx1minus0, consdata1->nvars) );
         diffidx1minus0size = consdata1->nvars;
      }

      /* because both constraints are normalized, a "<="-row and a ">="-row cannot be redundant */
      if( SCIPisInfinity(scip, -consdata0->lhs) != SCIPisInfinity(scip, -consdata1->lhs)
         && SCIPisInfinity(scip, consdata0->rhs) != SCIPisInfinity(scip, consdata1->rhs) )
         continue;

      /* check consdata0 against consdata1:
       * - if lhs0 >= lhs1 and for each variable v and each solution value x_v val0[v]*x_v <= val1[v]*x_v,
       *   consdata0 dominates consdata1 w.r.t. left hand side
       * - if rhs0 <= rhs1 and for each variable v and each solution value x_v val0[v]*x_v >= val1[v]*x_v,
       *   consdata0 dominates consdata1 w.r.t. right hand side
       * - if both at least one constraint is an equality, count the weighted number of common variables W_c
       *   and the weighted number of variable in the difference sets W_0 = w(V_0 \ V_1), W_1 = w(V_1 \ V_0),
       *   where the weight of each variable depends on its type, such that aggregations in order to remove the
       *   number of continuous and integer variables are prefered:
       *   - if W_c > W_1, try to aggregate  consdata0 := a * consdata0 + b * consdata1  in order to decrease the 
       *     variable weight in consdata0, where a = +/- val1[v] and b = -/+ val0[v] for common v which leads to
       *     the smallest weight; for numerical stability, we will only accept integral a and b; the sign of a has
       *     to be positive to not switch the sense of the (in)equality cons0
       *   - if W_c > W_0, try to aggregate  consdata1 := a * consdata1 + b * consdata0  in order to decrease the 
       *     variable weight in consdata1, where a = +/- val0[v] and b = -/+ val1[v] for common v which leads to
       *     the smallest weight; for numerical stability, we will only accept integral a and b; the sign of a has
       *     to be positive to not switch the sense of the (in)equality cons1
       */

      /* check consdata0 against consdata1 for redundancy */
      cons0dominateslhs = SCIPisGE(scip, consdata0->lhs, consdata1->lhs);
      cons1dominateslhs = SCIPisGE(scip, consdata1->lhs, consdata0->lhs);
      cons0dominatesrhs = SCIPisLE(scip, consdata0->rhs, consdata1->rhs);
      cons1dominatesrhs = SCIPisLE(scip, consdata1->rhs, consdata0->rhs);
      nvarscommon = 0;
      commonidxweight = 0;
      nvars0minus1 = 0;
      diffidx0minus1weight = 0;
      nvars1minus0 = 0;
      diffidx1minus0weight = 0;
      v0 = 0;
      v1 = 0;
      while( (v0 < consdata0->nvars || v1 < consdata1->nvars)
         && (cons0dominateslhs || cons1dominateslhs || cons0dominatesrhs || cons1dominatesrhs
            || cons0isequality || cons1isequality) )
      {
         /* test, if variable appears in only one or in both constraints */
         if( v0 < consdata0->nvars && v1 < consdata1->nvars )
            varcmp = SCIPvarCmp(consdata0->vars[v0], consdata1->vars[v1]);
         else if( v0 < consdata0->nvars )
            varcmp = -1;
         else
            varcmp = +1;

         switch( varcmp )
         {
         case -1:
            /* variable doesn't appear in consdata1 */
            var = consdata0->vars[v0];
            val0 = consdata0->vals[v0];
            val1 = 0.0;
            diffidx0minus1[nvars0minus1] = v0;
            nvars0minus1++;
            diffidx0minus1weight += getVarWeight(var);
            v0++;
            break;

         case +1:
            /* variable doesn't appear in consdata0 */
            var = consdata1->vars[v1];
            val0 = 0.0;
            val1 = consdata1->vals[v1];
            diffidx1minus0[nvars1minus0] = v1;
            nvars1minus0++;
            diffidx1minus0weight += getVarWeight(var);
            v1++;
            break;

         case 0:
            /* variable appears in both constraints */
            assert(consdata0->vars[v0] == consdata1->vars[v1]);
            var = consdata0->vars[v0];
            val0 = consdata0->vals[v0];
            val1 = consdata1->vals[v1];
            commonidx0[nvarscommon] = v0;
            commonidx1[nvarscommon] = v1;
            nvarscommon++;
            commonidxweight += getVarWeight(var);
            v0++;
            v1++;
            break;

         default:
            errorMessage("invalid comparison result\n");
            abort();
         }
         assert(var != NULL);

         /* update domination criteria w.r.t. the coefficient and the variable's bounds */
         if( SCIPisGT(scip, val0, val1) )
         {
            if( SCIPisNegative(scip, SCIPvarGetLbGlobal(var)) )
            {
               cons0dominatesrhs = FALSE;
               cons1dominateslhs = FALSE;
            }
            if( SCIPisPositive(scip, SCIPvarGetUbGlobal(var)) )
            {
               cons0dominateslhs = FALSE;
               cons1dominatesrhs = FALSE;
            }
         }
         else if( SCIPisLT(scip, val0, val1) )
         {
            if( SCIPisNegative(scip, SCIPvarGetLbGlobal(var)) )
            {
               cons0dominateslhs = FALSE;
               cons1dominatesrhs = FALSE;
            }
            if( SCIPisPositive(scip, SCIPvarGetUbGlobal(var)) )
            {
               cons0dominatesrhs = FALSE;
               cons1dominateslhs = FALSE;
            }
         }
      }

      /* check for domination */
      if( cons1dominateslhs && !SCIPisInfinity(scip, -consdata0->lhs) )
      {
         /* left hand side is dominated by consdata1: delete left hand side of consdata0 */
         debugMessage("left hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         debug(consdataPrint(scip, consdata0, NULL));
         debug(consdataPrint(scip, consdata1, NULL));
         /* check for infeasibility */
         if( SCIPisFeasGT(scip, consdata1->lhs, consdata0->rhs) )
         {
            debugMessage("linear constraint <%s> is infeasible\n", SCIPconsGetName(cons0));
            *result = SCIP_CUTOFF;
            continue;
         }
         CHECK_OKAY( chgLhs(scip, cons0, -SCIPinfinity(scip)) );
         cons0isequality = FALSE;
         if( !consdata0->upgraded )
            (*nchgsides)++;
         *result = SCIP_SUCCESS;
      }
      else if( cons0dominateslhs && !SCIPisInfinity(scip, -consdata1->lhs) )
      {
         /* left hand side is dominated by consdata0: delete left hand side of consdata1 */
         debugMessage("left hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons1), SCIPconsGetName(cons0));
         debug(consdataPrint(scip, consdata1, NULL));
         debug(consdataPrint(scip, consdata0, NULL));
         /* check for infeasibility */
         if( SCIPisFeasGT(scip, consdata0->lhs, consdata1->rhs) )
         {
            debugMessage("linear constraint <%s> is infeasible\n", SCIPconsGetName(cons1));
            *result = SCIP_CUTOFF;
            continue;
         }
         CHECK_OKAY( chgLhs(scip, cons1, -SCIPinfinity(scip)) );
         cons1isequality = FALSE;
         if( !consdata1->upgraded )
            (*nchgsides)++;
         *result = SCIP_SUCCESS;
      }
      if( cons1dominatesrhs && !SCIPisInfinity(scip, consdata0->rhs) )
      {
         /* right hand side is dominated by consdata1: delete right hand side of consdata0 */
         debugMessage("right hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         debug(consdataPrint(scip, consdata0, NULL));
         debug(consdataPrint(scip, consdata1, NULL));
         /* check for infeasibility */
         if( SCIPisFeasLT(scip, consdata1->rhs, consdata0->lhs) )
         {
            debugMessage("linear constraint <%s> is infeasible\n", SCIPconsGetName(cons0));
            *result = SCIP_CUTOFF;
            continue;
         }
         CHECK_OKAY( chgRhs(scip, cons0, SCIPinfinity(scip)) );
         cons0isequality = FALSE;
         if( !consdata0->upgraded )
            (*nchgsides)++;
         *result = SCIP_SUCCESS;
      }
      else if( cons0dominatesrhs && !SCIPisInfinity(scip, consdata1->rhs) )
      {
         /* right hand side is dominated by consdata0: delete right hand side of consdata1 */
         debugMessage("right hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons1), SCIPconsGetName(cons0));
         debug(consdataPrint(scip, consdata1, NULL));
         debug(consdataPrint(scip, consdata0, NULL));
         /* check for infeasibility */
         if( SCIPisFeasLT(scip, consdata0->rhs, consdata1->lhs) )
         {
            debugMessage("linear constraint <%s> is infeasible\n", SCIPconsGetName(cons1));
            *result = SCIP_CUTOFF;
            continue;
         }
         CHECK_OKAY( chgRhs(scip, cons1, SCIPinfinity(scip)) );
         cons1isequality = FALSE;
         if( !consdata1->upgraded )
            (*nchgsides)++;
         *result = SCIP_SUCCESS;
      }

      /* check for now redundant constraints */
      if( SCIPisInfinity(scip, -consdata0->lhs) && SCIPisInfinity(scip, consdata0->rhs) )
      {
         /* consdata0 became redundant */
         debugMessage("linear constraint <%s> is redundant\n", SCIPconsGetName(cons0));
         CHECK_OKAY( SCIPdelCons(scip, cons0) );
         if( !consdata0->upgraded )
            (*ndelconss)++;
         *result = SCIP_SUCCESS;
         continue;
      }
      if( SCIPisInfinity(scip, -consdata1->lhs) && SCIPisInfinity(scip, consdata1->rhs) )
      {
         /* consdata1 became redundant */
         debugMessage("linear constraint <%s> is redundant\n", SCIPconsGetName(cons1));
         CHECK_OKAY( SCIPdelCons(scip, cons1) );
         if( !consdata1->upgraded )
            (*ndelconss)++;
         *result = SCIP_SUCCESS;
         continue;
      }

      /* check, if we want to aggregate an (in)equality with an equality:
       *   consdata0 := a * consdata0 + b * consdata1  or  consdata1 := a * consdata1 + b * consdata0
       */
      if( cons0isequality || cons1isequality )
      {
         Bool aggregated;

         assert(consdata0->nvars == nvarscommon + nvars0minus1);
         assert(consdata1->nvars == nvarscommon + nvars1minus0);

         aggregated = FALSE;
         if( cons1isequality && !consdata0->upgraded && commonidxweight > diffidx1minus0weight )
         {
            /* W_c > W_1: try to aggregate  consdata0 := a * consdata0 + b * consdata1 */
            CHECK_OKAY( aggregateConstraints(scip, cons0, cons1, commonidx0, commonidx1, diffidx0minus1, diffidx1minus0, 
                           nvarscommon, commonidxweight, diffidx0minus1weight, diffidx1minus0weight, maxaggrnormscale,
                           nupgdconss, nchgcoefs, result, &aggregated) );
         }
         if( !aggregated && cons0isequality && !consdata1->upgraded && commonidxweight > diffidx0minus1weight )
         {
            /* W_c > W_0: try to aggregate  consdata1 := a * consdata1 + b * consdata0 */
            CHECK_OKAY( aggregateConstraints(scip, cons1, cons0, commonidx1, commonidx0, diffidx1minus0, diffidx0minus1,
                           nvarscommon, commonidxweight, diffidx1minus0weight, diffidx0minus1weight, maxaggrnormscale,
                           nupgdconss, nchgcoefs, result, &aggregated) );
         }
      }
   }

   /* free temporary memory */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &diffidx1minus0) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &diffidx0minus1) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &commonidx1) );
   CHECK_OKAY( SCIPfreeBufferArray(scip, &commonidx0) );

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
DECL_CONSPRESOL(consPresolLinear)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;
   CONS* cons;
   CONS* upgdcons;
   CONSDATA* consdata;
   Real minactivity;
   Real maxactivity;
   Bool consdeleted;
   Bool conschanged;
   Bool cutoff;
   int oldnfixedvars;
   int oldnaggrvars;
   int firstchange;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Presol method of linear constraints\n");*/

   *result = SCIP_DIDNOTFIND;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* process single constraints */
   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   firstchange = -1;
   for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
   {
      cons = conss[c];
      assert(SCIPconsIsActive(cons));
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->propagated = FALSE;

      /* check, if constraint is already propagated/preprocessed */
      if( consdata->propagated )
         continue;

      /* mark constraint being propagated/preprocessed */
      assert(SCIPconsIsActive(cons));
      consdata->propagated = TRUE;

      /* if inequality is already upgraded, delete it now; we only want to keep upgraded inequalities one presolving round
       * to help detecting redundancy of other linear constraints, but we want to keep equalities until the end of
       * presolving, because they may be used for aggregating other constraints
       */
      if( consdata->upgraded && SCIPisLT(scip, consdata->lhs, consdata->rhs) )
      {
         CHECK_OKAY( SCIPdelCons(scip, cons) );
         continue;
      }

      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == -1 && consdata->changed )
         firstchange = c;

      debugMessage("presolving linear constraint <%s>: ", SCIPconsGetName(cons));
      debug(consdataPrint(scip, consdata, NULL));

      consdeleted = FALSE;
      conschanged = FALSE;

      /* incorporate fixings and aggregations in constraint */
      if( nnewfixedvars > 0 || nnewaggrvars > 0 || *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars )
      {
         CHECK_OKAY( applyFixings(scip, cons, &conschanged) );
      }

      /* we can only presolve linear constraints, that are not modifiable */
      if( SCIPconsIsModifiable(cons) )
         continue;

      /* normalize constraint */
      CHECK_OKAY( normalizeCons(scip, cons) );

      /* tighten left and right hand side due to integrality */
      CHECK_OKAY( tightenSides(scip, cons, nchgsides, &conschanged) );

      /* check bounds */
      if( SCIPisFeasGT(scip, consdata->lhs, consdata->rhs) )
      {
         debugMessage("linear constraint <%s> is infeasible: sides=[%g,%g]\n",
            SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
         *result = SCIP_CUTOFF;
         break;
      }

      /* tighten variable's bounds */
      CHECK_OKAY( tightenBounds(scip, cons, &cutoff, nchgbds) );
      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      /* check for fixed variables */
      CHECK_OKAY( fixVariables(scip, cons, nfixedvars, result, &conschanged) );
      if( *result == SCIP_CUTOFF )
         break;

      /* check, if constraint is empty */
      if( consdata->nvars == 0 )
      {
         if( SCIPisFeasPositive(scip, consdata->lhs) || SCIPisFeasNegative(scip, consdata->rhs) )
         {
            debugMessage("linear constraint <%s> is empty and infeasible: sides=[%g,%g]\n",
               SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
            *result = SCIP_CUTOFF;
            break;
         }
         else
         {
            debugMessage("linear constraint <%s> is empty and redundant: sides=[%g,%g]\n",
               SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
            CHECK_OKAY( SCIPdelCons(scip, cons) );
            assert(!SCIPconsIsActive(cons));
            if( !consdata->upgraded )
               (*ndelconss)++;
            *result = SCIP_SUCCESS;
            continue;
         }
      }

      /* check constraint for infeasibility and redundancy */
      consdataGetActivityBounds(scip, consdata, &minactivity, &maxactivity);
      if( SCIPisFeasGT(scip, minactivity, consdata->rhs) || SCIPisFeasLT(scip, maxactivity, consdata->lhs) )
      {
         debugMessage("linear constraint <%s> is infeasible: activitybounds=[%g,%g], sides=[%g,%g]\n",
            SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
         *result = SCIP_CUTOFF;
         break;
      }
      else if( SCIPisGE(scip, minactivity, consdata->lhs) && SCIPisLE(scip, maxactivity, consdata->rhs) )
      {
         debugMessage("linear constraint <%s> is redundant: activitybounds=[%g,%g], sides=[%g,%g]\n",
            SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
         CHECK_OKAY( SCIPdelCons(scip, cons) );
         if( !consdata->upgraded )
            (*ndelconss)++;
         *result = SCIP_SUCCESS;
         continue;
      }
      else if( SCIPisGE(scip, minactivity, consdata->lhs) && !SCIPisInfinity(scip, -consdata->lhs) )
      {
         debugMessage("linear constraint <%s> left hand side is redundant: activitybounds=[%g,%g], sides=[%g,%g]\n",
            SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
         CHECK_OKAY( chgLhs(scip, cons, -SCIPinfinity(scip)) );
         if( !consdata->upgraded )
            (*nchgsides)++;
         *result = SCIP_SUCCESS;
         conschanged = TRUE;
      }
      else if( SCIPisLE(scip, maxactivity, consdata->rhs) && !SCIPisInfinity(scip, consdata->rhs) )
      {
         debugMessage("linear constraint <%s> right hand side is redundant: activitybounds=[%g,%g], sides=[%g,%g]\n",
            SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
         CHECK_OKAY( chgRhs(scip, cons, SCIPinfinity(scip)) );
         if( !consdata->upgraded )
            (*nchgsides)++;
         *result = SCIP_SUCCESS;
         conschanged = TRUE;
      }

      /* reduce big-M coefficients, that make the constraint redundant if the variable is on a bound */
      CHECK_OKAY( consdataTightenCoefs(scip, cons, nchgcoefs, nchgsides, &conschanged) );

      /* convert special equalities */
      CHECK_OKAY( convertEquality(scip, cons, nfixedvars, naggrvars, ndelconss, result,
                     &conschanged, &consdeleted) );
      if( *result == SCIP_CUTOFF || consdeleted )
         continue;

      /* try to upgrade the constraint into a more specific constraint type */
      CHECK_OKAY( SCIPupgradeConsLinear(scip, cons, &upgdcons) );
      if( upgdcons != NULL )
      {
         debugMessage("upgraded linear constraint <%s> to type <%s>:\n", 
            SCIPconsGetName(cons), SCIPconshdlrGetName(SCIPconsGetHdlr(upgdcons)));
         debug(consdataPrint(scip, consdata, NULL));
         
         /* add the upgraded constraint to the problem */
         CHECK_OKAY( SCIPaddCons(scip, upgdcons) );
         CHECK_OKAY( SCIPreleaseCons(scip, &upgdcons) );
         (*nupgdconss)++;
         
         /* mark the linear constraint being upgraded and to be removed after presolving;
          * don't delete it directly, because it may help to preprocess other linear constraints
          */
         assert(!consdata->upgraded);
         consdata->upgraded = TRUE;
      }
   }

   /* process pairs of constraints: check them for redundancy and try to aggregate them */
   if( *result != SCIP_CUTOFF && firstchange != -1 )
   {
      for( c = firstchange; c < nconss; ++c )
      {
         if( SCIPconsIsActive(conss[c]) )
         {
            CHECK_OKAY( preprocessConstraintPairs(scip, conss, firstchange, c, conshdlrdata->maxaggrnormscale,
                           nfixedvars, naggrvars, ndelconss, nupgdconss, nchgsides, nchgcoefs, result) );
         }
      }
      for( c = firstchange; c < nconss; ++c )
      {
         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         consdata->changed = FALSE;
      }
   }

   /* modify the result code */
   if( *result == SCIP_REDUCEDDOM )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** conflict variable resolving method of constraint handler */
#define consRescvarLinear NULL


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockLinear)
{  /*lint --e{715}*/
   consdataLockAllRoundings(scip, SCIPconsGetData(cons), nlockspos, nlocksneg);

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockLinear)
{  /*lint --e{715}*/
   consdataUnlockAllRoundings(scip, SCIPconsGetData(cons), nunlockspos, nunlocksneg);

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveLinear NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveLinear NULL


/** constraint enabling notification method of constraint handler */
#define consEnableLinear NULL


/** constraint disabling notification method of constraint handler */
#define consDisableLinear NULL




/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecLinear)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   VAR* var;
   Real oldbound;
   Real newbound;
   int varpos;
   EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   /*debugMessage("Exec method of bound change event handler for linear constraints\n");*/

   consdata = eventdata->consdata;
   varpos = eventdata->varpos;
   assert(consdata != NULL);
   assert(0 <= varpos && varpos < consdata->nvars);

   eventtype = SCIPeventGetType(event);
   var = SCIPeventGetVar(event);
   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);
   assert(var != NULL);
   assert(consdata->vars[varpos] == var);

   /*debugMessage(" -> eventtype=0x%x, var=<%s>, oldbound=%g, newbound=%g => activity: [%g,%g]", 
     eventtype, SCIPvarGetName(consdatadata->vars[varpos]), oldbound, newbound, consdatadata->minactivity, 
     consdatadata->maxactivity);*/

   if( (eventtype & SCIP_EVENTTYPE_LBCHANGED) != 0 )
      consdataUpdateChgLb(scip, consdata, var, oldbound, newbound, consdata->vals[varpos]);
   else
   {
      assert((eventtype & SCIP_EVENTTYPE_UBCHANGED) != 0);
      consdataUpdateChgUb(scip, consdata, var, oldbound, newbound, consdata->vals[varpos]);
   }

   consdata->propagated = FALSE;

   /*debug(printf(" -> [%g,%g]\n", consdatadata->minactivity, consdatadata->maxactivity));*/

   return SCIP_OKAY;
}




/*
 * Callback methods of conflict handler
 */

static
DECL_CONFLICTEXEC(conflictExecLinear)
{  /*lint --e{715}*/
   CONS* cons;
   CONS* upgdcons;
   Real* vals;
   char consname[MAXSTRLEN];
   int v;

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

   /* create array of ones for storing the coefficients: x1 + ... + xk >= 1 */
   CHECK_OKAY( SCIPallocBufferArray(scip, &vals, nconflictvars) );
   for( v = 0; v < nconflictvars; ++v )
      vals[v] = 1.0;

   /* create a constraint out of the conflict set */
   sprintf(consname, "cf%d", SCIPgetNGlobalConss(scip));
   CHECK_OKAY( SCIPcreateConsLinear(scip, &cons, consname, nconflictvars, conflictvars, vals, 1.0, SCIPinfinity(scip),
                  FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE) );

   /* release the vals buffer */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &vals) );

   /** try to automatically convert a linear constraint into a more specific and more specialized constraint */
   CHECK_OKAY( SCIPupgradeConsLinear(scip, cons, &upgdcons) );
   if( upgdcons != NULL )
   {
      CHECK_OKAY( SCIPreleaseCons(scip, &cons) );
      cons = upgdcons;
   }

   /* add constraint to SCIP */
   CHECK_OKAY( SCIPaddConsNode(scip, node, cons) );
   CHECK_OKAY( SCIPreleaseCons(scip, &cons) );

   *result = SCIP_CONSADDED;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for linear constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrLinear(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
                  NULL, NULL, NULL,
                  NULL, eventExecLinear,
                  NULL) );

   /* create conflict handler for linear constraints */
   CHECK_OKAY( SCIPincludeConflicthdlr(scip, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
                  NULL, NULL, NULL, conflictExecLinear,
                  NULL) );

   /* create constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler in SCIP */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
                  consFreeLinear, consInitLinear, consExitLinear, 
                  consInitpreLinear, consExitpreLinear, consInitsolLinear, consExitsolLinear,
                  consDeleteLinear, consTransLinear, 
                  consInitlpLinear, consSepaLinear, consEnfolpLinear, consEnfopsLinear, consCheckLinear, 
                  consPropLinear, consPresolLinear, consRescvarLinear,
                  consLockLinear, consUnlockLinear,
                  consActiveLinear, consDeactiveLinear, 
                  consEnableLinear, consDisableLinear,
                  conshdlrdata) );

   /* add linear constraint handler parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "constraints/linear/tightenboundsfreq",
                  "multiplier on propagation frequency, how often the bounds are tightened (-1: never, 0: only at root)",
                  &conshdlrdata->tightenboundsfreq, DEFAULT_TIGHTENBOUNDSFREQ, -1, INT_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "constraints/linear/maxaggrnormscale",
                  "maximal allowed relative gain in maximum norm for constraint aggregation",
                  &conshdlrdata->maxaggrnormscale, DEFAULT_MAXAGGRNORMSCALE, 0.0, REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** includes a linear constraint update method into the linear constraint handler */
RETCODE SCIPincludeLinconsUpgrade(
   SCIP*            scip,               /**< SCIP data structure */
   DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   int              priority            /**< priority of upgrading method */
   )
{
   CONSHDLR* conshdlr;
   CONSHDLRDATA* conshdlrdata;
   LINCONSUPGRADE* linconsupgrade;

   assert(linconsupgd != NULL);

   /* find the linear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("linear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create a linear constraint upgrade data object */
   CHECK_OKAY( linconsupgradeCreate(scip, &linconsupgrade, linconsupgd, priority) );

   /* insert linear constraint update method into constraint handler data */
   CHECK_OKAY( conshdlrdataIncludeUpgrade(scip, conshdlrdata, linconsupgrade) );

   return SCIP_OKAY;
}

/** creates and captures a linear constraint */
RETCODE SCIPcreateConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of constraint */
   Real             rhs,                /**< right hand side of constraint */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable during node processing (subject to col generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   /* find the linear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("linear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create the constraint specific data */
   if( SCIPstage(scip) == SCIP_STAGE_PROBLEM )
   {
      /* create constraint in original problem */
      CHECK_OKAY( consdataCreate(scip, &consdata, nvars, vars, vals, lhs, rhs) );
   }
   else
   {
      CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* create constraint in transformed problem */
      CHECK_OKAY( consdataCreateTransformed(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, vals, lhs, rhs) );
   }
   assert(consdata != NULL);

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
                  local, modifiable, removeable) );

   return SCIP_OKAY;
}

/** adds coefficient in linear constraint */
RETCODE SCIPaddCoefLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(var != NULL);

   /*debugMessage("adding coefficient %g * <%s> to linear constraint <%s>\n", 
     val, SCIPvarGetName(var), SCIPconsGetName(cons));*/

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }
   
   CHECK_OKAY( addCoef(scip, cons, var, val) );

   return SCIP_OKAY;
}

/** gets left hand side of linear constraint */
Real SCIPgetLhsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   )
{
   CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear\n");
      abort();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets right hand side of linear constraint */
Real SCIPgetRhsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   )
{
   CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear\n");
      abort();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** changes left hand side of linear constraint */
RETCODE SCIPchgLhsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real             lhs                 /**< new left hand side */
   )
{
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }
   
   CHECK_OKAY( chgLhs(scip, cons, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of linear constraint */
RETCODE SCIPchgRhsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real             rhs                 /**< new right hand side */
   )
{
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }
   
   CHECK_OKAY( chgRhs(scip, cons, rhs) );

   return SCIP_OKAY;
}

/** gets the dual solution of the linear constraint in the current LP */
Real SCIPgetDualsolLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   )
{
   CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear\n");
      abort();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** tries to automatically convert a linear constraint into a more specific and more specialized constraint */
RETCODE SCIPupgradeConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< source constraint to try to convert */
   CONS**           upgdcons            /**< pointer to store upgraded constraint, or NULL if not successful */
   )
{
   CONSHDLR* conshdlr;
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;
   VAR* var;
   Real val;
   Real lb;
   Real ub;
   Real poscoeffsum;
   Real negcoeffsum;
   Bool integral;
   int nposbin;
   int nnegbin;
   int nposint;
   int nnegint;
   int nposimpl;
   int nnegimpl;
   int nposcont;
   int nnegcont;
   int ncoeffspone;
   int ncoeffsnone;
   int ncoeffspint;
   int ncoeffsnint;
   int ncoeffspfrac;
   int ncoeffsnfrac;
   int i;

   assert(upgdcons != NULL);

   *upgdcons = NULL;

   /* we cannot upgrade a modifiable linear constraint, since we don't know what additional coefficients to expect */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* get the constraint handler and check, if it's really a linear constraint */
   conshdlr = SCIPconsGetHdlr(cons);
   if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   /* get constraint handler data and constraint data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check, if the constraint was already upgraded and will be deleted anyway after preprocessing */
   if( consdata->upgraded )
      return SCIP_OKAY;

   /* check, if the constraint is already stored as LP row */
   if( consdata->row != NULL )
   {
      if( SCIProwIsInLP(consdata->row) )
      {
         errorMessage("cannot upgrade linear constraint that is already stored as row in the LP\n");
         return SCIP_INVALIDDATA;
      }
      else
      {
         CHECK_OKAY( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   /* normalize constraint */
   CHECK_OKAY( normalizeCons(scip, cons) );


   /* 
    * calculate some statistics on linear constraint
    */

   nposbin = 0;
   nnegbin = 0;
   nposint = 0;
   nnegint = 0;
   nposimpl = 0;
   nnegimpl = 0;
   nposcont = 0;
   nnegcont = 0;
   ncoeffspone = 0;
   ncoeffsnone = 0;
   ncoeffspint = 0;
   ncoeffsnint = 0;
   ncoeffspfrac = 0;
   ncoeffsnfrac = 0;
   integral = TRUE;
   poscoeffsum = 0.0;
   negcoeffsum = 0.0;
   for( i = 0; i < consdata->nvars; ++i )
   {
      var = consdata->vars[i];
      val = consdata->vals[i];
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      assert(!SCIPisZero(scip, val));

      switch( SCIPvarGetType(var) )
      {
      case SCIP_VARTYPE_BINARY:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral &= SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposbin++;
         else
            nnegbin++;
         break;
      case SCIP_VARTYPE_INTEGER:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral &= SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposint++;
         else
            nnegint++;
         break;
      case SCIP_VARTYPE_IMPLINT:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral &= SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposimpl++;
         else
            nnegimpl++;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         integral = integral && SCIPisEQ(scip, lb, ub) && SCIPisIntegral(scip, val * lb);
         if( val >= 0.0 )
            nposcont++;
         else
            nnegcont++;
         break;
      default:
         errorMessage("unknown variable type\n");
         return SCIP_INVALIDDATA;
      }
      if( SCIPisEQ(scip, val, 1.0) )
         ncoeffspone++;
      else if( SCIPisEQ(scip, val, -1.0) )
         ncoeffsnone++;
      else if( SCIPisIntegral(scip, val) )
      {
         if( SCIPisPositive(scip, val) )
            ncoeffspint++;
         else
            ncoeffsnint++;
      }
      else
      {
         if( SCIPisPositive(scip, val) )
            ncoeffspfrac++;
         else
            ncoeffsnfrac++;
      }
      if( SCIPisPositive(scip, val) )
         poscoeffsum += val;
      else
         negcoeffsum += val;
   }



   /*
    * call the upgrading methods
    */

   debugMessage("upgrading linear constraint <%s> (%d upgrade methods):\n", 
      SCIPconsGetName(cons), conshdlrdata->nlinconsupgrades);
   debugMessage(" +bin=%d -bin=%d +int=%d -int=%d +impl=%d -impl=%d +cont=%d -cont=%d +1=%d -1=%d +I=%d -I=%d +F=%d -F=%d possum=%g negsum=%g integral=%d\n",
      nposbin, nnegbin, nposint, nnegint, nposimpl, nnegimpl, nposcont, nnegcont,
      ncoeffspone, ncoeffsnone, ncoeffspint, ncoeffsnint, ncoeffspfrac, ncoeffsnfrac,
      poscoeffsum, negcoeffsum, integral);

   /* try all upgrading methods in priority order */
   for( i = 0; i < conshdlrdata->nlinconsupgrades && *upgdcons == NULL; ++i )
   {
      CHECK_OKAY( conshdlrdata->linconsupgrades[i]->linconsupgd(scip, cons, consdata->nvars, 
                     consdata->vars, consdata->vals, consdata->lhs, consdata->rhs, 
                     nposbin, nnegbin, nposint, nnegint, nposimpl, nnegimpl, nposcont, nnegcont,
                     ncoeffspone, ncoeffsnone, ncoeffspint, ncoeffsnint, ncoeffspfrac, ncoeffsnfrac,
                     poscoeffsum, negcoeffsum, integral,
                     upgdcons) );
   }

#ifdef DEBUG
   if( *upgdcons != NULL )
   {
      conshdlr = SCIPconsGetHdlr(*upgdcons);
      debugMessage(" -> upgraded to constraint type <%s>\n", SCIPconshdlrGetName(conshdlr));
   }
#endif

   return SCIP_OKAY;
}
