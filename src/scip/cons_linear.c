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

/**@file   cons_linear.c
 * @brief  constraint handler for linear constraints
 * @author Tobias Achterberg
 */

/** Linear constraints can be stored in two different ways: as LP rows with
 *  coefficients for the columns, or as LINCONSDATA objects with coefficients for
 *  the variables. The second way is needed to be able to create the constraints
 *  before the solution process started, and before the necessary columns exist.
 *  At the first moment, a linear constraint is separated, it gets converted into
 *  an LP row.
 *
 *  Linear constraints are separated with a high priority, because they are easy
 *  to separate. The cut pool is implemented by adding linear constraints to the
 *  root node, such that it is separated each time, the linear constraints are
 *  separated. A constraint handler, which generates cuts for the pool should have
 *  a lower separation priority than the linear constraint handler, and it should
 *  have a separation frequency that is a multiple of the frequency of the linear
 *  constraint handler. In this way, it can be avoided to separate the same cut
 *  twice, because if a separation run of the handler is always preceded by a
 *  separation of the linear constraints, the pooled cuts are always satisfied.
 *
 *  Linear constraints are enforced and checked with a very low priority. Checking
 *  of (many) linear constraints is much more involved than checking the solution
 *  values for integrality. Because we are separating the linear constraints quite
 *  often, it is only necessary to enforce them for integral solutions. A constraint
 *  handler which generates pool cuts in its enforcing method should have an
 *  enforcing priority smaller than that of the linear constraint handler to avoid
 *  regenerating cuts which are already in the cut pool.
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
#define CONSHDLR_SEPAFREQ             4
#define CONSHDLR_PROPFREQ             4
#define CONSHDLR_NEEDSCONS         TRUE /**< the constraint handler should only be called, if linear constraints exist */

#define TIGHTENBOUNDSFREQ             5 /**< multiplier on propagation frequency, how often the bounds are tightened */

#define EVENTHDLR_NAME         "linear"
#define EVENTHDLR_DESC         "bound change event handler for linear constraints"


/** externally stored linear constraint data */
struct LinConsData
{
   char*            name;               /**< name of linear constraint */
   VAR**            vars;               /**< variables of constraint entries */
   Real*            vals;               /**< coefficients of constraint entries */
   EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   Real             lhs;                /**< left hand side of row (for ranged rows) */
   Real             rhs;                /**< right hand side of row */
   Real             minactivity;        /**< minimal value w.r.t. the variable's bounds for the constraint's activity,
                                         *   ignoring the coefficients contributing with infinite value */
   Real             maxactivity;        /**< maximal value w.r.t. the variable's bounds for the constraint's activity,
                                         *   ignoring the coefficients contributing with infinite value */
   int              minactivityinf;     /**< number of coefficients contributing with infinite value to minactivity */
   int              maxactivityinf;     /**< number of coefficients contributing with infinite value to maxactivity */
   int              varssize;           /**< size of the vars- and vals-arrays */
   int              nvars;              /**< number of nonzeros in constraint */
   unsigned int     local:1;            /**< is linear constraint only valid locally? */
   unsigned int     modifiable:1;       /**< is data modifiable during node processing (subject to column generation)? */
   unsigned int     transformed:1;      /**< does the linear constraint data belongs to the transformed problem? */
   unsigned int     validactivitybds:1; /**< are the activity bounds minactivity/maxactivity valid? */
};
typedef struct LinConsData LINCONSDATA; /**< externally stored linear constraint data */

/** linear constraint */
struct LinCons
{
   union
   {
      LINCONSDATA*  linconsdata;        /**< external linear constraint data structure */
      ROW*          row;                /**< LP row, if constraint is already stored in LP row format */
   } data;
   unsigned int     islprow:1;          /**< is linear constraint already stored as LP row? */
};

/** constraint data for linear constraints */
struct ConsData
{
   LINCONS*         lincons;            /**< linear constraint */
};

/** event data for bound change event */
struct EventData
{
   LINCONSDATA*     linconsdata;        /**< linear constraint data to process the bound change for */
   int              varpos;             /**< position of variable in vars array */
};

/** constraint handler data */
struct ConsHdlrData
{
   LINCONSUPGRADE** linconsupgrades;    /**< linear constraint upgrade methods for specializing linear constraints */
   int              linconsupgradessize;/**< size of linconsupgrade array */
   int              nlinconsupgrades;   /**< number of linear constraint upgrade methods */
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
RETCODE linconsdataEnsureVarsSize(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(linconsdata != NULL);
   assert(linconsdata->nvars <= linconsdata->varssize);
   
   if( num > linconsdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      CHECK_OKAY( SCIPreallocBlockMemoryArray(scip, &linconsdata->vars, linconsdata->varssize, newsize) );
      CHECK_OKAY( SCIPreallocBlockMemoryArray(scip, &linconsdata->vals, linconsdata->varssize, newsize) );
      if( linconsdata->transformed )
      {
         CHECK_OKAY( SCIPreallocBlockMemoryArray(scip, &linconsdata->eventdatas, linconsdata->varssize, newsize) );
      }
      else
         assert(linconsdata->eventdatas == NULL);
      linconsdata->varssize = newsize;
   }
   assert(num <= linconsdata->varssize);

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
RETCODE linconsupgradeFree(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSUPGRADE** linconsupgrade      /**< pointer to the linear constraint upgrade */
   )
{
   assert(linconsupgrade != NULL);
   assert(*linconsupgrade != NULL);

   SCIPfreeMemory(scip, linconsupgrade);

   return SCIP_OKAY;
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
 * linconsdata local methods
 */

/** creates event data for variable at given position, and catches events */
static
RETCODE linconsdataCatchEvent(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(linconsdata != NULL);
   assert(scip != NULL);
   assert(0 <= pos && pos < linconsdata->nvars);
   assert(linconsdata->vars != NULL);
   assert(linconsdata->vars[pos] != NULL);
   assert(linconsdata->eventdatas != NULL);
   assert(linconsdata->eventdatas[pos] == NULL);

   CHECK_OKAY( SCIPallocBlockMemory(scip, &linconsdata->eventdatas[pos]) );
   linconsdata->eventdatas[pos]->linconsdata = linconsdata;
   linconsdata->eventdatas[pos]->varpos = pos;

   CHECK_OKAY( SCIPcatchVarEvent(scip, linconsdata->vars[pos], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, 
                  linconsdata->eventdatas[pos]) );

   return SCIP_OKAY;
}

/** deletes event data for variable at given position, and drops events */
static
RETCODE linconsdataDropEvent(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int              pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(linconsdata != NULL);
   assert(scip != NULL);
   assert(0 <= pos && pos < linconsdata->nvars);
   assert(linconsdata->vars[pos] != NULL);
   assert(linconsdata->eventdatas[pos] != NULL);
   assert(linconsdata->eventdatas[pos]->linconsdata == linconsdata);
   assert(linconsdata->eventdatas[pos]->varpos == pos);
   
   CHECK_OKAY( SCIPdropVarEvent(scip, linconsdata->vars[pos], eventhdlr, linconsdata->eventdatas[pos]) );

   SCIPfreeBlockMemory(scip, &linconsdata->eventdatas[pos]);

   return SCIP_OKAY;
}

/** creates a linear constraint data object of the original problem */
static
RETCODE linconsdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA**    linconsdata,        /**< pointer to linear constraint data object */
   const char*      name,               /**< name of linear constraint */
   int              nvars,              /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             local,              /**< is linear constraint only valid locally? */
   Bool             modifiable          /**< is data modifiable during node processing (subject to column generation)? */
   )
{
   int i;

   assert(linconsdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   if( SCIPisGT(scip, lhs, rhs) )
   {
      char s[255];
      errorMessage("left hand side of linear constraint greater than right hand side");
      sprintf(s, "  (lhs=%f, rhs=%f)", lhs, rhs);
      errorMessage(s);
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPallocBlockMemory(scip, linconsdata) );

   CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*linconsdata)->name, name, strlen(name)+1) );

   if( nvars > 0 )
   {
      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*linconsdata)->vars, vars, nvars) );
      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*linconsdata)->vals, vals, nvars) );
   }
   else
   {
      (*linconsdata)->vars = NULL;
      (*linconsdata)->vals = NULL;
   }
   (*linconsdata)->eventdatas = NULL;

   (*linconsdata)->lhs = lhs;
   (*linconsdata)->rhs = rhs;
   (*linconsdata)->minactivity = SCIP_INVALID;
   (*linconsdata)->maxactivity = SCIP_INVALID;
   (*linconsdata)->minactivityinf = -1;
   (*linconsdata)->maxactivityinf = -1;
   (*linconsdata)->varssize = nvars;
   (*linconsdata)->nvars = nvars;
   (*linconsdata)->local = local;
   (*linconsdata)->modifiable = modifiable;
   (*linconsdata)->transformed = FALSE;
   (*linconsdata)->validactivitybds = FALSE;

   return SCIP_OKAY;
}

/** creates a linear constraint data object of the transformed problem */
static
RETCODE linconsdataCreateTransformed(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA**    linconsdata,        /**< pointer to linear constraint data object */
   const char*      name,               /**< name of linear constraint */
   int              nvars,              /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             local,              /**< is linear constraint only valid locally? */
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   )
{
   EVENTHDLR* eventhdlr;
   int i;

   assert(linconsdata != NULL);

   /* create linear constraint data */
   CHECK_OKAY( linconsdataCreate(scip, linconsdata, name, nvars, vars, vals, lhs, rhs, local, modifiable) );

   /* allocate the additional needed eventdatas array */
   assert((*linconsdata)->eventdatas == NULL);
   CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &(*linconsdata)->eventdatas, (*linconsdata)->varssize) );
   (*linconsdata)->transformed = TRUE;

   /* get event handler for updating linear constraint activity bounds */
   eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
   if( eventhdlr == NULL )
   {
      errorMessage("event handler for linear constraints not found");
      return SCIP_PLUGINNOTFOUND;
   }

   /* transform the variables, catch events */
   for( i = 0; i < (*linconsdata)->nvars; ++i )
   {
      if( SCIPvarGetStatus((*linconsdata)->vars[i]) == SCIP_VARSTATUS_ORIGINAL )
      {
         (*linconsdata)->vars[i] = SCIPvarGetTransformed((*linconsdata)->vars[i]);
         assert((*linconsdata)->vars[i] != NULL);
      }
      assert(SCIPvarGetStatus((*linconsdata)->vars[i]) != SCIP_VARSTATUS_ORIGINAL);
      (*linconsdata)->eventdatas[i] = NULL;
      CHECK_OKAY( linconsdataCatchEvent(scip, *linconsdata, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** frees a linear constraint data object */
static
RETCODE linconsdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA**    linconsdata         /**< pointer to linear constraint data object */
   )
{
   assert(linconsdata != NULL);
   assert(*linconsdata != NULL);
   assert((*linconsdata)->varssize >= 0);

   /* drop events for included variables */
   if( (*linconsdata)->transformed )
   {
      EVENTHDLR* eventhdlr;
      int i;

      /* get event handler for updating linear constraint activity bounds */
      eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
      if( eventhdlr == NULL )
      {
         errorMessage("event handler for linear constraints not found");
         return SCIP_PLUGINNOTFOUND;
      }
      
      /* drop all update activity bounds events for constraint's variables */
      for( i = 0; i < (*linconsdata)->nvars; ++i )
      {
         CHECK_OKAY( linconsdataDropEvent(scip, *linconsdata, eventhdlr, i) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*linconsdata)->eventdatas, (*linconsdata)->varssize);
   }
   assert((*linconsdata)->eventdatas == NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*linconsdata)->vars, (*linconsdata)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*linconsdata)->vals, (*linconsdata)->varssize);
   SCIPfreeBlockMemory(scip, linconsdata);

   return SCIP_OKAY;
}

/** adds coefficient in linear constraint data object */
static
RETCODE linconsdataAddCoef(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(linconsdata != NULL);
   assert(scip != NULL);
   assert(var != NULL);

   if( linconsdata->transformed && SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
   {
      var = SCIPvarGetTransformed(var);
      assert(var != NULL);
   }

   assert(linconsdata->transformed ^ (SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL));

   CHECK_OKAY( linconsdataEnsureVarsSize(scip, linconsdata, linconsdata->nvars+1) );
   linconsdata->vars[linconsdata->nvars] = var;
   linconsdata->vals[linconsdata->nvars] = val;
   linconsdata->nvars++;

   if( linconsdata->transformed )
   {
      EVENTHDLR* eventhdlr;

      /* get event handler for updating linear constraint activity bounds */
      eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
      if( eventhdlr == NULL )
      {
         errorMessage("event handler for linear constraints not found");
         return SCIP_PLUGINNOTFOUND;
      }

      /* catch bound change events on variable */
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
      linconsdata->eventdatas[linconsdata->nvars] = NULL;
      CHECK_OKAY( linconsdataCatchEvent(scip, linconsdata, eventhdlr, linconsdata->nvars-1) );
   }

   if( !linconsdata->local )
   {
      if( SCIPisPositive(scip, val) )
      {
         if( !SCIPisInfinity(scip, -linconsdata->lhs) )
            SCIPvarForbidRoundDown(var);
         if( !SCIPisInfinity(scip, linconsdata->rhs) )
            SCIPvarForbidRoundUp(var);
      }
      else if( SCIPisNegative(scip, val) )
      {
         if( !SCIPisInfinity(scip, linconsdata->rhs) )
            SCIPvarForbidRoundDown(var);
         if( !SCIPisInfinity(scip, -linconsdata->lhs) )
            SCIPvarForbidRoundUp(var);
      }
   }

   return SCIP_OKAY;
}

/** calculates the activity of the linear constraint for given solution */
static
RETCODE linconsdataGetActivity(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   SOL*             sol,                /**< solution to get activity for, NULL to actual solution */
   Real*            activity            /**< pointer to store the activity */
   )
{
   Real solval;
   int v;

   assert(linconsdata != NULL);
   assert(linconsdata->transformed);
   assert(activity != NULL);

   *activity = 0.0;
   for( v = 0; v < linconsdata->nvars; ++v )
   {
      solval = SCIPgetSolVal(scip, sol, linconsdata->vars[v]);
      *activity += linconsdata->vals[v] * solval;
   }

   return SCIP_OKAY;
}

/** calculates the feasibility of the linear constraint for given solution */
static
RETCODE linconsdataGetFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   SOL*             sol,                /**< solution to get feasibility for, NULL to actual solution */
   Real*            feasibility         /**< pointer to store the feasibility */
   )
{
   Real activity;

   assert(linconsdata != NULL);
   assert(linconsdata->transformed);
   assert(feasibility != NULL);

   CHECK_OKAY( linconsdataGetActivity(scip, linconsdata, sol, &activity) );
   *feasibility = MIN(linconsdata->rhs - activity, activity - linconsdata->lhs);

   return SCIP_OKAY;
}

/** creates an LP row from a linear constraint data object */
static
RETCODE linconsdataToRow(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   ROW**            row                 /**< pointer to an LP row data object */
   )
{
   int v;

   assert(linconsdata != NULL);
   assert(linconsdata->transformed);
   assert(row != NULL);

   CHECK_OKAY( SCIPcreateRow(scip, row, linconsdata->name, 0, NULL, NULL, linconsdata->lhs, linconsdata->rhs,
                  linconsdata->local, linconsdata->modifiable) );
   
   for( v = 0; v < linconsdata->nvars; ++v )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, *row, linconsdata->vars[v], linconsdata->vals[v]) );
   }

   return SCIP_OKAY;
}

/** updates minimum and maximum activity for a change in lower bound */
static
void linconsdataUpdateChgLb(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   Real             oldlb,              /**< old lower bound of variable */
   Real             newlb,              /**< new lower bound of variable */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(linconsdata != NULL);
   assert(linconsdata->transformed);
   
   if( linconsdata->validactivitybds )
   {
      assert(linconsdata->minactivity < SCIP_INVALID);
      assert(linconsdata->maxactivity < SCIP_INVALID);
      assert(linconsdata->minactivityinf >= 0);
      assert(linconsdata->maxactivityinf >= 0);
      assert(!SCIPisInfinity(scip, oldlb));
      assert(!SCIPisInfinity(scip, newlb));
      
      if( val > 0.0 )
      {
         if( SCIPisInfinity(scip, -oldlb) )
         {
            assert(linconsdata->minactivityinf >= 1);
            linconsdata->minactivityinf--;
         }
         else
            linconsdata->minactivity -= val * oldlb;

         if( SCIPisInfinity(scip, -newlb) )
            linconsdata->minactivityinf++;
         else
            linconsdata->minactivity += val * newlb;
      }
      else
      {
         if( SCIPisInfinity(scip, -oldlb) )
         {
            assert(linconsdata->maxactivityinf >= 1);
            linconsdata->maxactivityinf--;
         }
         else
            linconsdata->maxactivity -= val * oldlb;

         if( SCIPisInfinity(scip, -newlb) )
            linconsdata->maxactivityinf++;
         else
            linconsdata->maxactivity += val * newlb;
      }
   }
}

/** updates minimum and maximum activity for a change in upper bound */
static
void linconsdataUpdateChgUb(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   Real             oldub,              /**< old upper bound of variable */
   Real             newub,              /**< new upper bound of variable */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(linconsdata != NULL);
   assert(linconsdata->transformed);

   if( linconsdata->validactivitybds )
   {
      assert(linconsdata->minactivity < SCIP_INVALID);
      assert(linconsdata->maxactivity < SCIP_INVALID);
      assert(!SCIPisInfinity(scip, -oldub));
      assert(!SCIPisInfinity(scip, -newub));

      if( val > 0.0 )
      {
         if( SCIPisInfinity(scip, oldub) )
         {
            assert(linconsdata->maxactivityinf >= 1);
            linconsdata->maxactivityinf--;
         }
         else
            linconsdata->maxactivity -= val * oldub;

         if( SCIPisInfinity(scip, newub) )
            linconsdata->maxactivityinf++;
         else
            linconsdata->maxactivity += val * newub;
      }
      else
      {
         if( SCIPisInfinity(scip, oldub) )
         {
            assert(linconsdata->minactivityinf >= 1);
            linconsdata->minactivityinf--;
         }
         else
            linconsdata->minactivity -= val * oldub;

         if( SCIPisInfinity(scip, newub) )
            linconsdata->minactivityinf++;
         else
            linconsdata->minactivity += val * newub;
      }
   }
}

/** updates minimum and maximum activity for variable addition */
static
void linconsdataUpdateAddVar(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(linconsdata != NULL);
   assert(linconsdata->transformed);

   if( linconsdata->validactivitybds )
   {
      assert(linconsdata->minactivity < SCIP_INVALID);
      assert(linconsdata->maxactivity < SCIP_INVALID);

      linconsdataUpdateChgLb(scip, linconsdata, 0.0, SCIPvarGetLb(var), val);
      linconsdataUpdateChgUb(scip, linconsdata, 0.0, SCIPvarGetUb(var), val);
   }
}

/** calculates minimum and maximum activity for constraint */
static
void linconsdataCalcActivityBounds(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata         /**< linear constraint data object */
   )
{
   int i;
   
   assert(linconsdata != NULL);
   assert(linconsdata->transformed);
   assert(!linconsdata->validactivitybds);
   assert(linconsdata->minactivity >= SCIP_INVALID);
   assert(linconsdata->maxactivity >= SCIP_INVALID);
   
   linconsdata->validactivitybds = TRUE;
   linconsdata->minactivity = 0.0;
   linconsdata->maxactivity = 0.0;
   linconsdata->minactivityinf = 0;
   linconsdata->maxactivityinf = 0;

   for( i = 0; i < linconsdata->nvars; ++ i )
      linconsdataUpdateAddVar(scip, linconsdata, linconsdata->vars[i], linconsdata->vals[i]);
}

/** forbids roundings of variables in constraint that may violate constraint */
static
RETCODE linconsdataForbidRounding(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata         /**< linear constraint data object */
   )
{
   VAR** vars;
   Real* vals;
   Bool lhsexists;
   Bool rhsexists;
   int v;
   
   assert(linconsdata != NULL);
   assert(linconsdata->nvars == 0 || (linconsdata->vars != NULL && linconsdata->vals != NULL));
   assert(!SCIPisInfinity(scip, linconsdata->lhs));
   assert(!SCIPisInfinity(scip, -linconsdata->rhs));

   lhsexists = !SCIPisInfinity(scip, -linconsdata->lhs);
   rhsexists = !SCIPisInfinity(scip, linconsdata->rhs);
   vars = linconsdata->vars;
   vals = linconsdata->vals;

   for( v = 0; v < linconsdata->nvars; ++v )
   {
      assert(vars[v] != NULL);

      if( SCIPisPositive(scip, vals[v]) )
      {
         if( lhsexists )
            SCIPvarForbidRoundDown(vars[v]);
         if( rhsexists )
            SCIPvarForbidRoundUp(vars[v]);
      }
      else
      {
         assert(SCIPisNegative(scip, vals[v]));
         if( lhsexists )
            SCIPvarForbidRoundUp(vars[v]);
         if( rhsexists )
            SCIPvarForbidRoundDown(vars[v]);
      }
   }
   
   return SCIP_OKAY;
}

/** allows roundings of variables in constraint that may violate constraint */
static
RETCODE linconsdataAllowRounding(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONSDATA*     linconsdata         /**< linear constraint data object */
   )
{
   VAR** vars;
   Real* vals;
   Bool lhsexists;
   Bool rhsexists;
   int v;
   
   assert(linconsdata != NULL);
   assert(linconsdata->nvars == 0 || (linconsdata->vars != NULL && linconsdata->vals != NULL));
   assert(!SCIPisInfinity(scip, linconsdata->lhs));
   assert(!SCIPisInfinity(scip, -linconsdata->rhs));

   lhsexists = !SCIPisInfinity(scip, -linconsdata->lhs);
   rhsexists = !SCIPisInfinity(scip, linconsdata->rhs);
   vars = linconsdata->vars;
   vals = linconsdata->vals;

   for( v = 0; v < linconsdata->nvars; ++v )
   {
      assert(vars[v] != NULL);

      if( SCIPisPositive(scip, vals[v]) )
      {
         if( lhsexists )
            SCIPvarAllowRoundDown(vars[v]);
         if( rhsexists )
            SCIPvarAllowRoundUp(vars[v]);
      }
      else
      {
         assert(SCIPisNegative(scip, vals[v]));
         if( lhsexists )
            SCIPvarAllowRoundUp(vars[v]);
         if( rhsexists )
            SCIPvarAllowRoundDown(vars[v]);
      }
   }

   return SCIP_OKAY;
}

/** prints linear constraint to file stream */
static
void linconsdataPrint(
   LINCONSDATA*     linconsdata,        /**< linear constraint data object */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(linconsdata != NULL);

   if( file == NULL )
      file = stdout;

   /* print left hand side */
   fprintf(file, "%+f <= ", linconsdata->lhs);

   /* print coefficients */
   if( linconsdata->nvars == 0 )
      fprintf(file, "0 ");
   for( v = 0; v < linconsdata->nvars; ++v )
   {
      assert(linconsdata->vars[v] != NULL);
      assert(linconsdata->vars[v]->name != NULL);
      fprintf(file, "%+f%s ", linconsdata->vals[v], linconsdata->vars[v]->name);
   }

   /* print right hand side */
   fprintf(file, "<= %+f\n", linconsdata->rhs);
}



/*
 * lincons methods
 */

/** creates a linear constraint object */
RETCODE SCIPlinconsCreate(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS**        lincons,            /**< pointer to store the linear constraint */
   const char*      name,               /**< name of linear constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            var,                /**< array with variables of constraint entries */
   Real*            val,                /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             local,              /**< is linear constraint only valid locally? */
   Bool             modifiable          /**< is constraint modifiable during node processing (sbj. to column generation)? */
   )
{
   assert(lincons != NULL);
   assert(len == 0 || var != NULL);
   assert(len == 0 || val != NULL);

   /* create the constraint specific data */
   CHECK_OKAY( SCIPallocBlockMemory(scip, lincons) );
   (*lincons)->islprow = FALSE;
   if( SCIPstage(scip) == SCIP_STAGE_PROBLEM )
   {
      if( local )
      {
         errorMessage("problem constraint cannot be local");
         return SCIP_INVALIDDATA;
      }

      /* create constraint in original problem */
      CHECK_OKAY( linconsdataCreate(scip, &(*lincons)->data.linconsdata, name, len, var, val, lhs, rhs,
                     local, modifiable) );
   }
   else
   {
      /* create constraint in transformed problem */
      CHECK_OKAY( linconsdataCreateTransformed(scip, &(*lincons)->data.linconsdata, name, len, var, val, lhs, rhs, 
                     local, modifiable) );
   }

   /* forbid rounding of variables */
   if( !local )
   {
      CHECK_OKAY( SCIPlinconsForbidRounding(scip, *lincons) );
   }

   return SCIP_OKAY;
}

/** frees a linear constraint object */
RETCODE SCIPlinconsFree(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS**        lincons             /**< pointer to store the linear constraint */
   )
{
   assert(lincons != NULL);
   assert(*lincons != NULL);

   if( (*lincons)->islprow )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*lincons)->data.row) );
   }
   else
   {
      CHECK_OKAY( linconsdataFree(scip, &(*lincons)->data.linconsdata) );
   }
   SCIPfreeBlockMemory(scip, lincons);

   return SCIP_OKAY;
}

/** adds coefficient in linear constraint */
RETCODE SCIPlinconsAddCoef(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(lincons != NULL);

   if( lincons->islprow )
   {
      assert(lincons->data.row != NULL);
      CHECK_OKAY( SCIPaddVarToRow(scip, lincons->data.row, var, val) );
   }
   else
   {
      assert(lincons->data.linconsdata != NULL);
      CHECK_OKAY( linconsdataAddCoef(scip, lincons->data.linconsdata, var, val) );
   }

   return SCIP_OKAY;
}

/** forbids roundings of variables in linear constraint that may violate the constraint */
RETCODE SCIPlinconsForbidRounding(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons             /**< linear constraint */
   )
{
   assert(lincons != NULL);

   if( lincons->islprow )
   {
      CHECK_OKAY( SCIPforbidRowRounding(scip, lincons->data.row) );
   }
   else
   {
      CHECK_OKAY( linconsdataForbidRounding(scip, lincons->data.linconsdata) );
   }

   return SCIP_OKAY;
}

/** allows roundings of variables in linear constraint that may violate the constraint */
RETCODE SCIPlinconsAllowRounding(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons             /**< linear constraint */
   )
{
   assert(lincons != NULL);

   if( lincons->islprow )
   {
      CHECK_OKAY( SCIPallowRowRounding(scip, lincons->data.row) );
   }
   else
   {
      CHECK_OKAY( linconsdataAllowRounding(scip, lincons->data.linconsdata) );
   }

   return SCIP_OKAY;
}

/** gets activity bounds for constraint */
RETCODE SCIPlinconsGetActivityBounds(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   Real*            minactivity,        /**< pointer to store the minimal activity */
   Real*            maxactivity         /**< pointer to store the maximal activity */
   )
{
   assert(lincons != NULL);
   assert(scip != NULL);
   assert(minactivity != NULL);
   assert(maxactivity != NULL);

   if( lincons->islprow )
   {
      CHECK_OKAY( SCIPgetRowActivityBounds(scip, lincons->data.row, minactivity, maxactivity) );
   }
   else
   {
      LINCONSDATA* linconsdata;

      linconsdata = lincons->data.linconsdata;
      assert(linconsdata != NULL);
      if( !linconsdata->validactivitybds )
         linconsdataCalcActivityBounds(scip, linconsdata);
      assert(linconsdata->minactivity < SCIP_INVALID);
      assert(linconsdata->maxactivity < SCIP_INVALID);

      if( linconsdata->minactivityinf > 0 )
         *minactivity = -SCIPinfinity(scip);
      else
         *minactivity = linconsdata->minactivity;
      if( linconsdata->maxactivityinf > 0 )
         *maxactivity = SCIPinfinity(scip);
      else
         *maxactivity = linconsdata->maxactivity;
   }

   return SCIP_OKAY;
}

/** invalidates activity bounds, such that they are recalculated in next get */
RETCODE SCIPlinconsInvalidActivityBounds(
   LINCONS*         lincons             /**< linear constraint */
   )
{
   assert(lincons != NULL);

   if( lincons->islprow )
   {
      CHECK_OKAY( SCIProwInvalidateActivityBounds(lincons->data.row) );
   }
   else
   {
      LINCONSDATA* linconsdata;

      linconsdata = lincons->data.linconsdata;
      assert(linconsdata != NULL);
      linconsdata->validactivitybds = FALSE;
      linconsdata->minactivity = SCIP_INVALID;
      linconsdata->maxactivity = SCIP_INVALID;
      linconsdata->minactivityinf = -1;
      linconsdata->maxactivityinf = -1;
   }

   return SCIP_OKAY;
}

/** gets activity bounds for constraint after setting variable to zero */
RETCODE SCIPlinconsGetActivityResiduals(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   VAR*             var,                /**< variable to calculate activity residual for */
   Real             val,                /**< coefficient value of variable in linear constraint */
   Real*            minresactivity,     /**< pointer to store the minimal residual activity */
   Real*            maxresactivity      /**< pointer to store the maximal residual activity */
   )
{
   assert(lincons != NULL);
   assert(scip != NULL);
   assert(var != NULL);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);

   if( lincons->islprow )
   {
      CHECK_OKAY( SCIPgetRowActivityResiduals(scip, lincons->data.row, var, val, minresactivity, maxresactivity) );
   }
   else
   {
      LINCONSDATA* linconsdata;
      Real lb;
      Real ub;

      linconsdata = lincons->data.linconsdata;
      assert(linconsdata != NULL);
      if( !linconsdata->validactivitybds )
         linconsdataCalcActivityBounds(scip, linconsdata);
      assert(linconsdata->minactivity < SCIP_INVALID);
      assert(linconsdata->maxactivity < SCIP_INVALID);
      assert(linconsdata->minactivityinf >= 0);
      assert(linconsdata->maxactivityinf >= 0);

      lb = SCIPvarGetLb(var);
      ub = SCIPvarGetUb(var);
      assert(!SCIPisInfinity(scip, lb));
      assert(!SCIPisInfinity(scip, -ub));

      if( val > 0.0 )
      {
         if( SCIPisInfinity(scip, -lb) )
         {
            assert(linconsdata->minactivityinf >= 1);
            if( linconsdata->minactivityinf >= 2 )
               *minresactivity = -SCIPinfinity(scip);
            else
               *minresactivity = linconsdata->minactivity;
         }
         else
         {
            if( linconsdata->minactivityinf >= 1 )
               *minresactivity = -SCIPinfinity(scip);
            else
               *minresactivity = linconsdata->minactivity - val * lb;
         }
         if( SCIPisInfinity(scip, ub) )
         {
            assert(linconsdata->maxactivityinf >= 1);
            if( linconsdata->maxactivityinf >= 2 )
               *maxresactivity = +SCIPinfinity(scip);
            else
               *maxresactivity = linconsdata->maxactivity;
         }
         else
         {
            if( linconsdata->maxactivityinf >= 1 )
               *maxresactivity = +SCIPinfinity(scip);
            else
               *maxresactivity = linconsdata->maxactivity - val * ub;
         }
      }
      else
      {
         if( SCIPisInfinity(scip, ub) )
         {
            assert(linconsdata->minactivityinf >= 1);
            if( linconsdata->minactivityinf >= 2 )
               *minresactivity = -SCIPinfinity(scip);
            else
               *minresactivity = linconsdata->minactivity;
         }
         else
         {
            if( linconsdata->minactivityinf >= 1 )
               *minresactivity = -SCIPinfinity(scip);
            else
               *minresactivity = linconsdata->minactivity - val * ub;
         }
         if( SCIPisInfinity(scip, -lb) )
         {
            assert(linconsdata->maxactivityinf >= 1);
            if( linconsdata->maxactivityinf >= 2 )
               *maxresactivity = +SCIPinfinity(scip);
            else
               *maxresactivity = linconsdata->maxactivity;
         }
         else
         {
            if( linconsdata->maxactivityinf >= 1 )
               *maxresactivity = +SCIPinfinity(scip);
            else
               *maxresactivity = linconsdata->maxactivity - val * lb;
         }
      }
   }

   return SCIP_OKAY;
}

/** gets left hand side of linear constraint */
RETCODE SCIPlinconsGetLhs(
   LINCONS*         lincons,            /**< linear constraint */
   Real*            lhs                 /**< pointer to store left hand side */
   )
{
   assert(lincons != NULL);
   assert(lhs != NULL);

   if( lincons->islprow )
   {
      assert(lincons->data.row != NULL);
      *lhs = SCIProwGetLhs(lincons->data.row);
   }
   else
   {
      assert(lincons->data.linconsdata != NULL);
      *lhs = lincons->data.linconsdata->lhs;
   }

   return SCIP_OKAY;
}

/** gets right hand side of linear constraint */
RETCODE SCIPlinconsGetRhs(
   LINCONS*         lincons,            /**< linear constraint */
   Real*            rhs                 /**< pointer to store right hand side */
   )
{
   assert(lincons != NULL);
   assert(rhs != NULL);

   if( lincons->islprow )
   {
      assert(lincons->data.row != NULL);
      *rhs = SCIProwGetRhs(lincons->data.row);
   }
   else
   {
      assert(lincons->data.linconsdata != NULL);
      *rhs = lincons->data.linconsdata->rhs;
   }

   return SCIP_OKAY;
}

/** sets left hand side of linear constraint */
RETCODE SCIPlinconsChgLhs(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   Real             lhs                 /**< new left hand side */
   )
{
   assert(lincons != NULL);

   if( lincons->islprow )
   {
      todoMessage("introduce side change methods for rows");
      errorMessage("cannot change sides of linear constraint already stored as LP row (not implemented yet)");
      return SCIP_INVALIDDATA;
   }
   else
   {
      LINCONSDATA* linconsdata;
      Bool updaterounding;

      linconsdata = lincons->data.linconsdata;
      assert(linconsdata != NULL);

      updaterounding = (!linconsdata->local && SCIPisInfinity(scip, -lhs) != SCIPisInfinity(scip, -linconsdata->lhs));
      if( updaterounding )
      {
         CHECK_OKAY( SCIPlinconsAllowRounding(scip, lincons) );
      }
      linconsdata->lhs = lhs;
      if( updaterounding )
      {
         CHECK_OKAY( SCIPlinconsForbidRounding(scip, lincons) );
      }
   }

   return SCIP_OKAY;
}

/** sets right hand side of linear constraint */
RETCODE SCIPlinconsChgRhs(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   Real             rhs                 /**< new right hand side */
   )
{
   assert(lincons != NULL);

   if( lincons->islprow )
   {
      todoMessage("introduce side change methods for rows");
      errorMessage("cannot change sides of linear constraint already stored as LP row (not implemented yet)");
      return SCIP_INVALIDDATA;
   }
   else
   {
      LINCONSDATA* linconsdata;
      Bool updaterounding;

      linconsdata = lincons->data.linconsdata;
      assert(linconsdata != NULL);

      updaterounding = (!linconsdata->local && SCIPisInfinity(scip, rhs) != SCIPisInfinity(scip, linconsdata->lhs));
      if( updaterounding )
      {
         CHECK_OKAY( SCIPlinconsAllowRounding(scip, lincons) );
      }
      linconsdata->rhs = rhs;
      if( updaterounding )
      {
         CHECK_OKAY( SCIPlinconsForbidRounding(scip, lincons) );
      }
   }

   return SCIP_OKAY;
}

/** prints linear constraint to file stream */
void SCIPlinconsPrint(
   LINCONS*         lincons,            /**< linear constraint */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   assert(lincons != NULL);

   if( lincons->islprow )
   {
      assert(lincons->data.row != NULL);
      SCIProwPrint(lincons->data.row, file);
   }
   else
   {
      assert(lincons->data.linconsdata != NULL);
      linconsdataPrint(lincons->data.linconsdata, file);
   }
}

/** tightens bounds of a single variable due to activity bounds */
static
RETCODE tightenVarBounds(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,           /**< constraint data */
   VAR*             var,                /**< variable to tighten bounds for */
   Real             val,                /**< coefficient value of variable in linear constraint */
   RESULT*          result,             /**< pointer to store SCIP_CUTOFF, if node is infeasible */
   Bool*            success             /**< pointer to store whether a bound was tightened */
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

   assert(var != NULL);
   assert(!SCIPisZero(scip, val));
   assert(result != NULL);
   assert(success != NULL);

   CHECK_OKAY( SCIPlinconsGetLhs(lincons, &lhs) );
   CHECK_OKAY( SCIPlinconsGetRhs(lincons, &rhs) );
   CHECK_OKAY( SCIPlinconsGetActivityResiduals(scip, lincons, var, val, &minresactivity, &maxresactivity) );
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));
   assert(!SCIPisInfinity(scip, minresactivity));
   assert(!SCIPisInfinity(scip, -maxresactivity));
   
   lb = SCIPvarGetLb(var);
   ub = SCIPvarGetUb(var);
   assert(SCIPisLE(scip, lb, ub));

   *success = FALSE;
   if( val > 0.0 )
   {
      /* check, if we can tighten the variable's bounds */
      if( !SCIPisInfinity(scip, -minresactivity) && !SCIPisInfinity(scip, rhs) )
      {
         newub = (rhs - minresactivity)/val;
         if( SCIPisSumLT(scip, newub, ub) )
         {
            /* tighten upper bound */
            debugMessage("linear constraint: tighten <%s>, old bds=[%f,%f], val=%g, resactivity=[%g,%g], sides=[%g,%g]\n",
               SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs);
            debugMessage("  -> newub=%f\n", newub);
            if( SCIPisSumLT(scip, newub, lb) )
            {
               debugMessage("linear constraint: cutoff  <%s>, new bds=[%f,%f]\n", SCIPvarGetName(var), lb, newub);
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            CHECK_OKAY( SCIPchgVarUb(scip, var, newub) );
            ub = SCIPvarGetUb(var); /* get bound again, because it may be additionally modified due to integrality */
            *success = TRUE;
            debugMessage("linear constraint: tighten <%s>, new bds=[%f,%f]\n", SCIPvarGetName(var), lb, ub);
         }
      }
      if( !SCIPisInfinity(scip, maxresactivity) && !SCIPisInfinity(scip, -lhs) )
      {
         newlb = (lhs - maxresactivity)/val;
         if( SCIPisSumGT(scip, newlb, lb) )
         {
            /* tighten lower bound */
            debugMessage("linear constraint: tighten <%s>, old bds=[%f,%f], val=%g, resactivity=[%g,%g], sides=[%g,%g]\n",
               SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs);
            debugMessage("  -> newlb=%f\n", newlb);
            if( SCIPisSumGT(scip, newlb, ub) )
            {
               debugMessage("linear constraint: cutoff  <%s>, new bds=[%f,%f]\n", SCIPvarGetName(var), newlb, ub);
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            CHECK_OKAY( SCIPchgVarLb(scip, var, newlb) );
            lb = SCIPvarGetLb(var); /* get bound again, because it may be additionally modified due to integrality */
            *success = TRUE;
            debugMessage("linear constraint: tighten <%s>, new bds=[%f,%f]\n", SCIPvarGetName(var), lb, ub);
         }
      }
   }
   else
   {
      /* check, if we can tighten the variable's bounds */
      if( !SCIPisInfinity(scip, -minresactivity) && !SCIPisInfinity(scip, rhs) )
      {
         newlb = (rhs - minresactivity)/val;
         if( SCIPisSumGT(scip, newlb, lb) )
         {
            /* tighten lower bound */
            debugMessage("linear constraint: tighten <%s>, old bds=[%f,%f], val=%g, resactivity=[%g,%g], sides=[%g,%g]\n",
               SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs);
            debugMessage("  -> newlb=%f\n", newlb);
            if( SCIPisSumGT(scip, newlb, ub) )
            {
               debugMessage("linear constraint: cutoff  <%s>, new bds=[%f,%f]\n", SCIPvarGetName(var), newlb, ub);
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            CHECK_OKAY( SCIPchgVarLb(scip, var, newlb) );
            lb = SCIPvarGetLb(var); /* get bound again, because it may be additionally modified due to integrality */
            *success = TRUE;
            debugMessage("linear constraint: tighten <%s>, new bds=[%f,%f]\n", SCIPvarGetName(var), lb, ub);
         }
      }
      if( !SCIPisInfinity(scip, maxresactivity) && !SCIPisInfinity(scip, -lhs) )
      {
         newub = (lhs - maxresactivity)/val;
         if( SCIPisSumLT(scip, newub, ub) )
         {
            /* tighten upper bound */
            debugMessage("linear constraint: tighten <%s>, old bds=[%f,%f], val=%g, resactivity=[%g,%g], sides=[%g,%g]\n",
               SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs);
            debugMessage("  -> newub=%f\n", newub);
            if( SCIPisSumLT(scip, newub, lb) )
            {
               debugMessage("linear constraint: cutoff  <%s>, new bds=[%f,%f]\n", SCIPvarGetName(var), lb, newub);
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            CHECK_OKAY( SCIPchgVarUb(scip, var, newub) );
            ub = SCIPvarGetUb(var); /* get bound again, because it may be additionally modified due to integrality */
            *success = TRUE;
            debugMessage("linear constraint: tighten <%s>, new bds=[%f,%f]\n", SCIPvarGetName(var), lb, ub);
         }
      }
   }
   
   return SCIP_OKAY;
}

/** tightens variable's bounds due to activity bounds */
RETCODE SCIPlinconsTightenBounds(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   RESULT*          result              /**< pointer to store SCIP_CUTOFF, if node is infeasible */
   )
{
   Bool success;
   int lastsuccess;

   assert(lincons != NULL);
   assert(scip != NULL);
   assert(result != NULL);
   assert(*result != SCIP_CUTOFF);

   if( lincons->islprow )
   {
      ROW* row;
      COL** cols;
      Real* vals;
      int ncols;
      int c;

      row = lincons->data.row;
      ncols = SCIProwGetNNonz(row);
      if( ncols > 0 )
      {
         cols = SCIProwGetCols(row);
         vals = SCIProwGetVals(row);
         assert(cols != NULL);
         assert(vals != NULL);
         lastsuccess = 0;
         c = 0;
         do
         {
            assert(0 <= c && c < ncols);
            assert(cols[c] != NULL);
            
            CHECK_OKAY( tightenVarBounds(scip, lincons, SCIPcolGetVar(cols[c]), vals[c], result, &success) );
            if( success )
            {
               *result = SCIP_REDUCEDDOM;
               lastsuccess = c;
            }
            c++;
            if( c == ncols )
               c = 0;
         }
         while( c != lastsuccess && *result != SCIP_CUTOFF );
      }
   }
   else
   {
      LINCONSDATA* linconsdata;
      VAR** vars;
      Real* vals;
      int nvars;
      int v;

      linconsdata = lincons->data.linconsdata;
      assert(linconsdata != NULL);
      nvars = linconsdata->nvars;
      if( nvars > 0 )
      {
         vars = linconsdata->vars;
         vals = linconsdata->vals;
         assert(vars != NULL);
         assert(vals != NULL);
         lastsuccess = 0;
         v = 0;
         do
         {
            assert(0 <= v && v < nvars);
            CHECK_OKAY( tightenVarBounds(scip, lincons, vars[v], vals[v], result, &success) );
            if( success )
            {
               *result = SCIP_REDUCEDDOM;
               lastsuccess = v;
            }
            v++;
            if( v == nvars )
               v = 0;
         }
         while( v != lastsuccess && *result != SCIP_CUTOFF );
      }
   }

   return SCIP_OKAY;
}

/** separates linear constraint: adds linear constraint as cut, if violated by current LP solution */
RETCODE SCIPlinconsSeparate(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   Bool*            violated            /**< pointer to store information, if linear constraint was violated */
   )
{
   assert(lincons != NULL);
   assert(violated != NULL);

   *violated = FALSE;

   if( lincons->islprow )
   {
      ROW* row = lincons->data.row;
      
      if( !SCIProwIsInLP(row) )
      {         
         Real feasibility;
         
         feasibility = SCIPgetRowFeasibility(scip, row);
         /*debugMessage("  row feasibility = %g\n", feasibility);*/
         if( !SCIPisFeasible(scip, feasibility) )
         {
            /* insert LP row as cut */
            CHECK_OKAY( SCIPaddCut(scip, row, -feasibility/SCIProwGetNorm(row)/(SCIProwGetNNonz(row)+1)) );
            *violated = TRUE;
         }
      }
   }
   else
   {
      LINCONSDATA* linconsdata = lincons->data.linconsdata;
      Real feasibility;
      
      /* check feasibility of linear constraint */
      CHECK_OKAY( linconsdataGetFeasibility(scip, linconsdata, NULL, &feasibility) );
      /*debugMessage("  linconsdata feasibility = %g\n", feasibility);*/
      if( !SCIPisFeasible(scip, feasibility) )
      {
         ROW* row;
         
         /* convert linconsdata data into LP row */
         CHECK_OKAY( linconsdataToRow(scip, linconsdata, &row) );
#ifndef NDEBUG
         {
            Real rowfeasibility;
            rowfeasibility = SCIPgetRowFeasibility(scip, row);
            assert(SCIPisSumEQ(scip, rowfeasibility, feasibility));
         }
#endif

         /* free the linconsdata data and convert lincons to point to the row */
         CHECK_OKAY( linconsdataFree(scip, &lincons->data.linconsdata) );
         lincons->islprow = TRUE;
         lincons->data.row = row;
         
         /* don't release the row, because we need it as data storage */
         
         /* insert LP row as cut */
         CHECK_OKAY( SCIPaddCut(scip, row, -feasibility/SCIProwGetNorm(row)/(SCIProwGetNNonz(row)+1)) );
         *violated = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** checks linear constraint for feasibility of given solution or actual LP/pseudo solution */
RETCODE SCIPlinconsCheck(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   SOL*             sol,                /**< solution to be checked, or NULL for actual solution */
   Bool             checklprows,        /**< has linear constraint to be checked, if it is already in current LP? */
   Bool*            violated            /**< pointer to store information, if linear constraint is violated */
   )
{
   assert(lincons != NULL);
   assert(violated != NULL);

   *violated = FALSE;

   if( lincons->islprow )
   {
      ROW* row = lincons->data.row;

      if( checklprows || !SCIProwIsInLP(row) )
      {         
         Real feasibility;
         
         feasibility = SCIPgetRowSolFeasibility(scip, row, sol);
         /*debugMessage("  row feasibility = %g\n", feasibility);*/
         *violated |= !SCIPisFeasible(scip, feasibility);
      }
   }
   else
   {
      LINCONSDATA* linconsdata = lincons->data.linconsdata;
      Real feasibility;
      
      /* check feasibility of linear constraint */
      CHECK_OKAY( linconsdataGetFeasibility(scip, linconsdata, sol, &feasibility) );
      /*debugMessage("  linconsdata feasibility = %g\n", feasibility);*/
      *violated |= !SCIPisFeasible(scip, feasibility);
   }

   return SCIP_OKAY;
}




/*
 * local linear constraint handler methods
 */

/** separates violated inequalities; called from Sepa and Enfo methods */
static
RETCODE separateConstraints(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< linear constraint handler */
   CONS**           conss,              /**< array of constraints to process */
   int              nconss,             /**< number of constraints to process */
   Bool*            found               /**< pointer to store information, if a violated constraint has been found */
   )
{
   CONSDATA* consdata;
   Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(nconss >= 0);
   assert(found != NULL);

   *found = FALSE;

   debugMessage("separating %d linear constraints at %p\n", nconss, conss);
   for( c = 0; c < nconss; ++c )
   {
      /*debugMessage("separating linear constraint <%s>\n", SCIPconsGetName(conss[c]));*/
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      CHECK_OKAY( SCIPlinconsSeparate(scip, consdata->lincons, &violated) );
      *found |= violated;

      /* update age of constraint */
      if( violated )
      {
         CHECK_OKAY( SCIPresetConsAge(scip, conss[c]) );
      }
      else
      {
         CHECK_OKAY( SCIPincConsAge(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** checks pseudo solution for violated inequalities */
static
RETCODE checkConstraints(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLR*        conshdlr,           /**< linear constraint handler */
   CONS**           conss,              /**< array of constraints to process */
   int              nconss,             /**< number of constraints to process */
   SOL*             sol,                /**< solution to be checked, NULL to check pseudo solution */
   Bool             checklprows,         /**< have current LP rows to be checked? */
   Bool*            found               /**< pointer to store information, if a violated constraint has been found */
   )
{
   CONSDATA* consdata;
   Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(nconss >= 0);
   assert(found != NULL);

   *found = FALSE;

   debugMessage("checking solution %p for %d linear constraints at %p\n", sol, nconss, conss);
   for( c = 0; c < nconss && !(*found); ++c )
   {
      /*debugMessage("checking linear constraint <%s>\n", SCIPconsGetName(conss[c]));*/
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      CHECK_OKAY( SCIPlinconsCheck(scip, consdata->lincons, sol, checklprows, &violated) );
      *found |= violated;

      /* update age of constraint */
      if( violated )
      {
         CHECK_OKAY( SCIPresetConsAge(scip, conss[c]) );
      }
      else
      {
         CHECK_OKAY( SCIPincConsAge(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

static
DECL_CONSFREE(consFreeLinear)
{
   CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

static
DECL_CONSDELETE(consDeleteLinear)
{
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* free linear constraint */
   CHECK_OKAY( SCIPlinconsFree(scip, &(*consdata)->lincons) );

   /* free constraint data object */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

static
DECL_CONSTRANS(consTransLinear)
{
   CONSDATA* sourcedata;
   CONSDATA* targetdata;
   LINCONSDATA* linconsdata;

   /*debugMessage("Trans method of linear constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(SCIPstage(scip) == SCIP_STAGE_INITSOLVE);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->lincons != NULL);
   assert(!sourcedata->lincons->islprow);  /* in original problem, there cannot be LP rows */
   assert(sourcedata->lincons->data.linconsdata != NULL);

   /* create constraint data for target constraint */
   CHECK_OKAY( SCIPallocBlockMemory(scip, &targetdata) );
   
   linconsdata = sourcedata->lincons->data.linconsdata;

   todoMessage("normalize linear constraints");

   CHECK_OKAY( SCIPlinconsCreate(scip, &targetdata->lincons, linconsdata->name, 
                  linconsdata->nvars, linconsdata->vars, linconsdata->vals, linconsdata->lhs, linconsdata->rhs,
                  linconsdata->local, linconsdata->modifiable) );

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                  SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
                  SCIPconsIsPropagated(sourcecons)) );

   /* try to upgrade target linear constraint into more specific constraint */
   CHECK_OKAY( SCIPupgradeConsLinear(scip, targetcons) );

   return SCIP_OKAY;
}

static
DECL_CONSSEPA(consSepaLinear)
{
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /*debugMessage("Sepa method of linear constraints\n");*/

   /* separate violated constraints; only process useful ones */
   CHECK_OKAY( separateConstraints(scip, conshdlr, conss, nusefulconss, &found) );

   if( found )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

static
DECL_CONSENFOLP(consEnfolpLinear)
{
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /*debugMessage("Enfolp method of linear constraints\n");*/

   /* check for violated constraints */

   /* LP is processed at current node -> we can add violated linear constraints to the LP */

   /* process useful constraints first */
   CHECK_OKAY( separateConstraints(scip, conshdlr, conss, nusefulconss, &found) );

   /* if no violation was found, process obsolete constraints */
   if( !found )
   {
      CHECK_OKAY( separateConstraints(scip, conshdlr, &conss[nusefulconss], nconss - nusefulconss, &found) );
   }
   
   if( found )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
DECL_CONSENFOPS(consEnfopsLinear)
{
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /*debugMessage("Enfops method of linear constraints\n");*/

   /* check for violated constraints */

   /* LP is not processed at current node -> we just have to check pseudo solution for feasibility */
   CHECK_OKAY( checkConstraints(scip, conshdlr, conss, nconss, NULL, TRUE, &found) );
   
   if( found )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
DECL_CONSCHECK(consCheckLinear)
{
   Bool found;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   CHECK_OKAY( checkConstraints(scip, conshdlr, conss, nconss, sol, checklprows, &found) );

   if( found )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

static
DECL_CONSPROP(consPropLinear)
{
   CONS* cons;
   CONSDATA* consdata;
   Real minactivity;
   Real maxactivity;
   Real lhs;
   Real rhs;
   Bool redundant;
   Bool tightenbounds;
   int propfreq;
   int actdepth;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /*debugMessage("Prop method of linear constraints\n");*/

   /* check, if we want to tighten variable's bounds */
   propfreq = SCIPconshdlrGetPropFreq(conshdlr);
   actdepth = SCIPgetActDepth(scip);
   tightenbounds = (actdepth % (propfreq * TIGHTENBOUNDSFREQ) == 0);

   /* process useful constraints */
   *result = SCIP_DIDNOTFIND;
   for( c = 0; c < nusefulconss && *result != SCIP_CUTOFF; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* tighten the variable's bounds */
      if( tightenbounds )
      {
         CHECK_OKAY( SCIPlinconsTightenBounds(scip, consdata->lincons, result) );
#ifndef NDEBUG
         {
            Real newminactivity;
            Real newmaxactivity;
            Real recalcminactivity;
            Real recalcmaxactivity;
            
            CHECK_OKAY( SCIPlinconsGetActivityBounds(scip, consdata->lincons, &newminactivity, &newmaxactivity) );
            CHECK_OKAY( SCIPlinconsInvalidActivityBounds(consdata->lincons) );
            CHECK_OKAY( SCIPlinconsGetActivityBounds(scip, consdata->lincons, &recalcminactivity, &recalcmaxactivity) );

            assert(SCIPisSumRelEQ(scip, newminactivity, recalcminactivity));
            assert(SCIPisSumRelEQ(scip, newmaxactivity, recalcmaxactivity));
         }
#endif
      }

      /* check constraint for infeasibility and redundancy */
      CHECK_OKAY( SCIPlinconsGetActivityBounds(scip, consdata->lincons, &minactivity, &maxactivity) );
      CHECK_OKAY( SCIPlinconsGetLhs(consdata->lincons, &lhs) );
      CHECK_OKAY( SCIPlinconsGetRhs(consdata->lincons, &rhs) );

      if( SCIPisGT(scip, minactivity, rhs) || SCIPisLT(scip, maxactivity, lhs) )
      {
         debugMessage("linear constraint <%s> is infeasible: activitybounds=[%g,%g], sides=[%g,%g]\n",
            SCIPconsGetName(cons), minactivity, maxactivity, lhs, rhs);
         CHECK_OKAY( SCIPresetConsAge(scip, cons) );
         *result = SCIP_CUTOFF;
      }
      else if( SCIPisGE(scip, minactivity, lhs) && SCIPisLE(scip, maxactivity, rhs) )
      {
         debugMessage("linear constraint <%s> is redundant: activitybounds=[%g,%g], sides=[%g,%g]\n",
            SCIPconsGetName(cons), minactivity, maxactivity, lhs, rhs);
         CHECK_OKAY( SCIPincConsAge(scip, cons) );
         CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      }
   }

   return SCIP_OKAY;
}



/*
 * Callback methods of event handler
 */

static
DECL_EVENTEXEC(eventExecLinear)
{
   LINCONSDATA* linconsdata;
   Real oldbound;
   Real newbound;
   int varpos;
   EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(event != NULL);

   /*debugMessage("Exec method of bound change event handler for linear constraints\n");*/

   linconsdata = eventdata->linconsdata;
   varpos = eventdata->varpos;
   assert(linconsdata != NULL);
   assert(0 <= varpos && varpos < linconsdata->nvars);

   CHECK_OKAY( SCIPeventGetType(event, &eventtype) );
   CHECK_OKAY( SCIPeventGetOldbound(event, &oldbound) );
   CHECK_OKAY( SCIPeventGetNewbound(event, &newbound) );
#ifndef NDEBUG
   {
      VAR* var;
      CHECK_OKAY( SCIPeventGetVar(event, &var) );
      assert(linconsdata->vars[varpos] == var);
   }
#endif

   /*debugMessage(" -> eventtype=0x%x, var=<%s>, oldbound=%g, newbound=%g => activity: [%g,%g]", 
     eventtype, SCIPvarGetName(linconsdata->vars[varpos]), oldbound, newbound, linconsdata->minactivity, 
     linconsdata->maxactivity);*/

   if( (eventtype & SCIP_EVENTTYPE_LBCHANGED) != 0 )
      linconsdataUpdateChgLb(scip, linconsdata, oldbound, newbound, linconsdata->vals[varpos]);
   else
   {
      assert((eventtype & SCIP_EVENTTYPE_UBCHANGED) != 0);
      linconsdataUpdateChgUb(scip, linconsdata, oldbound, newbound, linconsdata->vals[varpos]);
   }

   /*debug(printf(" -> [%g,%g]\n", linconsdata->minactivity, linconsdata->maxactivity));*/

   return SCIP_OKAY;
}



/*
 * constraint specific interface methods
 */

/** creates the handler for linear constraints and includes it in SCIP */
RETCODE SCIPincludeConsHdlrLinear(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
                  NULL, NULL, NULL,
                  NULL, eventExecLinear,
                  NULL) );

   /* create constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler in SCIP */
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ,
                  CONSHDLR_NEEDSCONS,
                  consFreeLinear, NULL, NULL,
                  consDeleteLinear, consTransLinear, 
                  consSepaLinear, consEnfolpLinear, consEnfopsLinear, consCheckLinear, consPropLinear,
                  NULL, NULL,
                  conshdlrdata) );

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
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("Linear constraint handler not found");
      return SCIP_INVALIDCALL;
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
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is linear constraint only valid locally? */
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   assert(scip != NULL);

   /* find the linear constraint handler */
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("Linear constraint handler not found");
      return SCIP_INVALIDCALL;
   }

   /* create the constraint specific data */
   CHECK_OKAY( SCIPallocBlockMemory(scip, &consdata) );
   CHECK_OKAY( SCIPlinconsCreate(scip, &consdata->lincons, name, nvars, vars, vals, lhs, rhs, local, modifiable) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, separate, enforce, check, propagate) );

   return SCIP_OKAY;
}

/** adds coefficient in linear constraint */
RETCODE SCIPaddCoefConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);
   assert(var != NULL);

   /*debugMessage("adding coefficient %g * <%s> to linear constraint <%s>\n", val, var->name, SCIPconsGetName(cons));*/

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   CHECK_OKAY( SCIPlinconsAddCoef(scip, consdata->lincons, var, val) );

   return SCIP_OKAY;
}

/** gets left hand side of linear constraint */
RETCODE SCIPgetLhsConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real*            lhs                 /**< pointer to store left hand side */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);
   assert(lhs != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   CHECK_OKAY( SCIPlinconsGetLhs(consdata->lincons, lhs) );

   return SCIP_OKAY;
}

/** gets right hand side of linear constraint */
RETCODE SCIPgetRhsConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real*            rhs                 /**< pointer to store right hand side */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);
   assert(rhs != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   CHECK_OKAY( SCIPlinconsGetRhs(consdata->lincons, rhs) );

   return SCIP_OKAY;
}

/** changes left hand side of linear constraint */
RETCODE SCIPchgLhsConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real             lhs                 /**< new left hand side */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   CHECK_OKAY( SCIPlinconsChgLhs(scip, consdata->lincons, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of linear constraint */
RETCODE SCIPchgRhsConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real             rhs                 /**< new right hand side */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(scip != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   CHECK_OKAY( SCIPlinconsChgRhs(scip, consdata->lincons, rhs) );

   return SCIP_OKAY;
}

/** tries to automatically convert a linear constraint into a more specific and more specialized constraint */
RETCODE SCIPupgradeConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons                /**< pointer to constraint to convert */
   )
{
   CONSHDLR* conshdlr;
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;
   LINCONSDATA* linconsdata;
   CONS* upgdcons;
   VAR* var;
   Real val;
   Real lb;
   Real ub;
   Bool integral;
   Bool upgraded;
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

   assert(scip != NULL);
   assert(cons != NULL);
   assert(*cons != NULL);

   conshdlr = SCIPconsGetHdlr(*cons);
   if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not linear");
      return SCIP_INVALIDDATA;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   consdata = SCIPconsGetData(*cons);
   assert(consdata != NULL);
   assert(consdata->lincons != NULL);

   if( consdata->lincons->islprow )
   {
      errorMessage("cannot upgrade linear constraint that is already stored as LP row");
      return SCIP_INVALIDDATA;
   }

   linconsdata = consdata->lincons->data.linconsdata;
   assert(linconsdata != NULL);

   /* we cannot upgrade a modifiable linear constraint, since we don't know what additional coefficients to expect */
   if( linconsdata->modifiable )
      return SCIP_OKAY;

   /* calculate some statistics on linear constraint */
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
   for( i = 0; i < linconsdata->nvars; ++i )
   {
      var = linconsdata->vars[i];
      val = linconsdata->vals[i];
      lb = SCIPvarGetLb(var);
      ub = SCIPvarGetUb(var);
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
      case SCIP_VARTYPE_CONTINOUS:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral &= (SCIPisIntegral(scip, val) && SCIPisEQ(scip, lb, ub));
         if( val >= 0.0 )
            nposcont++;
         else
            nnegcont++;
         break;
      default:
         errorMessage("unknown variable type");
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
   }

   debugMessage("upgrading linear constraint <%s> (%d upgrade methods):\n", 
      SCIPconsGetName(*cons), conshdlrdata->nlinconsupgrades);
   debugMessage(" +bin=%d -bin=%d +int=%d -int=%d +impl=%d -impl=%d +cont=%d -cont=%d +1=%d -1=%d +I=%d -I=%d +F=%d -F=%d integral=%d\n",
      nposbin, nnegbin, nposint, nnegint, nposimpl, nnegimpl, nposcont, nnegcont,
      ncoeffspone, ncoeffsnone, ncoeffspint, ncoeffsnint, ncoeffspfrac, ncoeffsnfrac, integral);

   /* try all upgrading methods in priority order */
   upgdcons = NULL;
   upgraded = FALSE;
   for( i = 0; i < conshdlrdata->nlinconsupgrades && !upgraded; ++i )
   {
      CHECK_OKAY( conshdlrdata->linconsupgrades[i]->linconsupgd(scip, *cons, linconsdata->nvars, 
                     linconsdata->vars, linconsdata->vals, linconsdata->lhs, linconsdata->rhs, linconsdata->local,
                     nposbin, nnegbin, nposint, nnegint, nposimpl, nnegimpl, nposcont, nnegcont,
                     ncoeffspone, ncoeffsnone, ncoeffspint, ncoeffsnint, ncoeffspfrac, ncoeffsnfrac, integral,
                     &upgdcons, &upgraded) );
      assert(upgraded ^ (upgdcons == NULL));
   }

   /* if upgrading was successful, release the old constraint and use the upgraded constraint instead */
   if( upgraded )
   {
      assert(upgdcons != NULL);
      CHECK_OKAY( SCIPreleaseCons(scip, cons) );
      *cons = upgdcons;
   }

   return SCIP_OKAY;
}
