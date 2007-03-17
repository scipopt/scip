/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
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
#pragma ident "@(#) $Id: cons_linear.c,v 1.232 2007/03/17 18:09:11 bzfpfend Exp $"

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

#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"


#define CONSHDLR_NAME          "linear"
#define CONSHDLR_DESC          "linear constraints of the form  lhs <= a^T x <= rhs"
#define CONSHDLR_SEPAPRIORITY   +100000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             5 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define EVENTHDLR_NAME         "linear"
#define EVENTHDLR_DESC         "bound change event handler for linear constraints"

#define CONFLICTHDLR_NAME      "linear"
#define CONFLICTHDLR_DESC      "conflict handler creating linear constraints"
#define CONFLICTHDLR_PRIORITY  -1000000

#define DEFAULT_TIGHTENBOUNDSFREQ     1 /**< multiplier on propagation frequency, how often the bounds are tightened */
#define DEFAULT_MAXROUNDS             5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     200 /**< maximal number of cuts separated per separation round in root node */
#define DEFAULT_MAXPRESOLPAIRROUNDS  -1 /**< maximal number of presolving rounds with pairwise constraint comparison
                                         *   (-1: no limit) */
#define DEFAULT_MAXAGGRNORMSCALE    0.0 /**< maximal allowed relative gain in maximum norm for constraint aggregation
                                         *   (0.0: disable constraint aggregation) */
#define DEFAULT_SEPARATEALL       FALSE /**< should all constraints be subject to cover cut generation instead of only
                                         *   the ones with non-zero dual value? */

#define KNAPSACKRELAX_MAXDELTA      0.1 /**< maximal allowed rounding distance for scaling in knapsack relaxation */
#define KNAPSACKRELAX_MAXDNOM    1000LL /**< maximal allowed denominator in knapsack rational relaxation */
#define KNAPSACKRELAX_MAXSCALE   1000.0 /**< maximal allowed scaling factor in knapsack rational relaxation */

#define MAXDNOM                 10000LL /**< maximal denominator for simple rational fixed values */
#define MAXSCALEDCOEF             1e+06 /**< maximal coefficient value after scaling */


/** constraint data for linear constraints */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of row (for ranged rows) */
   SCIP_Real             rhs;                /**< right hand side of row */
   SCIP_Real             maxabsval;          /**< maximum absolute value of all coefficients */
   SCIP_Real             pseudoactivity;     /**< pseudo activity value in current pseudo solution */
   SCIP_Real             minactivity;        /**< minimal value w.r.t. the variable's local bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             maxactivity;        /**< maximal value w.r.t. the variable's local bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             glbminactivity;     /**< minimal value w.r.t. the variable's global bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             glbmaxactivity;     /**< maximal value w.r.t. the variable's global bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Longint          possignature;       /**< bit signature of coefficients that may take a positive value */
   SCIP_Longint          negsignature;       /**< bit signature of coefficients that may take a negative value */
   SCIP_ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   SCIP_VAR**            vars;               /**< variables of constraint entries */
   SCIP_Real*            vals;               /**< coefficients of constraint entries */
   SCIP_EVENTDATA**      eventdatas;         /**< event datas for bound change events of the variables */
   int                   minactivityneginf;  /**< number of coefficients contributing with neg. infinite value to minactivity */
   int                   minactivityposinf;  /**< number of coefficients contributing with pos. infinite value to minactivity */
   int                   maxactivityneginf;  /**< number of coefficients contributing with neg. infinite value to maxactivity */
   int                   maxactivityposinf;  /**< number of coefficients contributing with pos. infinite value to maxactivity */
   int                   glbminactivityneginf;/**< number of coefficients contrib. with neg. infinite value to glbminactivity */
   int                   glbminactivityposinf;/**< number of coefficients contrib. with pos. infinite value to glbminactivity */
   int                   glbmaxactivityneginf;/**< number of coefficients contrib. with neg. infinite value to glbmaxactivity */
   int                   glbmaxactivityposinf;/**< number of coefficients contrib. with pos. infinite value to glbmaxactivity */
   int                   varssize;           /**< size of the vars- and vals-arrays */
   int                   nvars;              /**< number of nonzeros in constraint */
   unsigned int          validmaxabsval:1;   /**< is the maximum absolute value valid? */
   unsigned int          validactivities:1;  /**< are the pseudo activity and activity bounds (local and global) valid? */
   unsigned int          propagated:1;       /**< is constraint already propagated? */
   unsigned int          boundstightened:1;  /**< is constraint already propagated with bound tightening? */
   unsigned int          presolved:1;        /**< is constraint already presolved? */
   unsigned int          removedfixings:1;   /**< are all fixed variables removed from the constraint? */
   unsigned int          validsignature:1;   /**< is the bit signature valid? */
   unsigned int          changed:1;          /**< was constraint changed since last aggregation round in preprocessing? */
   unsigned int          normalized:1;       /**< is the constraint in normalized form? */
   unsigned int          upgradetried:1;     /**< was the constraint already tried to be upgraded? */
   unsigned int          upgraded:1;         /**< is the constraint upgraded and will it be removed after preprocessing? */
   unsigned int          sorted:1;           /**< are the constraint's variables sorted? */
   unsigned int          merged:1;           /**< are the constraint's equal variables already merged? */
   unsigned int          cliquesadded:1;     /**< were the cliques of the constraint already extracted? */
};

/** event data for bound change event */
struct SCIP_EventData
{
   SCIP_CONSDATA*        consdata;           /**< linear constraint data to process the bound change for */
   int                   varpos;             /**< position of variable in vars array */
   int                   filterpos;          /**< position of event in variable's event filter */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   SCIP_LINCONSUPGRADE** linconsupgrades;    /**< linear constraint upgrade methods for specializing linear constraints */
   SCIP_Real             maxaggrnormscale;   /**< maximal allowed relative gain in maximum norm for constraint aggregation
                                              *   (0.0: disable constraint aggregation) */
   int                   linconsupgradessize;/**< size of linconsupgrade array */
   int                   nlinconsupgrades;   /**< number of linear constraint upgrade methods */
   int                   tightenboundsfreq;  /**< multiplier on propagation frequency, how often the bounds are tightened */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cuts separated per separation round in root node */
   int                   maxpresolpairrounds;/**< maximal number of presolving rounds with pairwise constraint comparison
                                              *   (-1: no limit) */
   SCIP_Bool             separateall;        /**< should all constraints be subject to cover cut generation instead of only
                                              *   the ones with non-zero dual value? */
};

/** linear constraint update method */
struct SCIP_LinConsUpgrade
{
   SCIP_DECL_LINCONSUPGD((*linconsupgd));    /**< method to call for upgrading linear constraint */
   int                   priority;           /**< priority of upgrading method */
};




/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_1_RHS        = 1,                /**< activity residuals of all other variables tighten bounds of single
                                              *   variable due to the right hand side of the inequality */
   PROPRULE_1_LHS        = 2,                /**< activity residuals of all other variables tighten bounds of single
                                              *   variable due to the left hand side of the inequality */
   PROPRULE_INVALID      = 0                 /**< propagation was applied without a specific propagation rule */
};
typedef enum Proprule PROPRULE;

/** inference information */
struct InferInfo
{
   union
   {
      struct
      {
         unsigned int    proprule:8;         /**< propagation rule that was applied */
         unsigned int    pos:24;             /**< variable position, the propagation rule was applied at */
      } asbits;
      int                asint;              /**< inference information as a single int value */
   } val;
};
typedef struct InferInfo INFERINFO;

/** convers an integer into an inference information */
static
INFERINFO intToInferInfo(
   int                   i                   /**< integer to convert */
   )
{
   INFERINFO inferinfo;

   inferinfo.val.asint = i;

   return inferinfo;
}

/** converts an inference information into an int */
static
int inferInfoToInt(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asint;
}

/** returns the propagation rule stored in the inference information */
static
int inferInfoGetProprule(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asbits.proprule;
}

/** returns the position stored in the inference information */
static
int inferInfoGetPos(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asbits.pos;
}

/** constructs an inference information out of a propagation rule and a position number */
static
INFERINFO getInferInfo(
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   int                   pos                 /**< variable position, the propagation rule was applied at */
   )
{
   INFERINFO inferinfo;

   inferinfo.val.asbits.proprule = proprule; /*lint !e641*/
   inferinfo.val.asbits.pos = pos; /*lint !e732*/

   return inferinfo;
}

/** constructs an inference information out of a propagation rule and a position number, returns info as int */
static
int getInferInt(
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   int                   pos                 /**< variable position, the propagation rule was applied at */
   )
{
   return inferInfoToInt(getInferInfo(proprule, pos));
}




/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that linconsupgrades array can store at least num entries */
static
SCIP_RETCODE conshdlrdataEnsureLinconsupgradesSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< linear constraint handler data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->nlinconsupgrades <= conshdlrdata->linconsupgradessize);

   if( num > conshdlrdata->linconsupgradessize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->linconsupgrades, newsize) );
      conshdlrdata->linconsupgradessize = newsize;
   }
   assert(num <= conshdlrdata->linconsupgradessize);

   return SCIP_OKAY;
}

/** ensures, that vars and vals arrays can store at least num entries */
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
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vals, consdata->varssize, newsize) );
      if( consdata->eventdatas != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->eventdatas, consdata->varssize, newsize) );
      }
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
SCIP_RETCODE linconsupgradeCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LINCONSUPGRADE** linconsupgrade,     /**< pointer to store the linear constraint upgrade */
   SCIP_DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   int                   priority            /**< priority of upgrading method */
   )
{
   assert(linconsupgrade != NULL);
   assert(linconsupgd != NULL);

   SCIP_CALL( SCIPallocMemory(scip, linconsupgrade) );
   (*linconsupgrade)->linconsupgd = linconsupgd;
   (*linconsupgrade)->priority = priority;

   return SCIP_OKAY;
}

/** frees a linear constraint upgrade data object */
static
void linconsupgradeFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LINCONSUPGRADE** linconsupgrade      /**< pointer to the linear constraint upgrade */
   )
{
   assert(linconsupgrade != NULL);
   assert(*linconsupgrade != NULL);

   SCIPfreeMemory(scip, linconsupgrade);
}

/** creates constaint handler data for linear constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );
   (*conshdlrdata)->linconsupgrades = NULL;
   (*conshdlrdata)->linconsupgradessize = 0;
   (*conshdlrdata)->nlinconsupgrades = 0;

   /* get event handler for updating linear constraint activity bounds */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for linear constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** frees constraint handler data for linear constraint handler */
static
void conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
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
SCIP_RETCODE conshdlrdataIncludeUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_LINCONSUPGRADE*  linconsupgrade      /**< linear constraint upgrade method */
   )
{
   int i;

   assert(conshdlrdata != NULL);
   assert(linconsupgrade != NULL);

   SCIP_CALL( conshdlrdataEnsureLinconsupgradesSize(scip, conshdlrdata, conshdlrdata->nlinconsupgrades+1) );

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

/** installs rounding locks for the given variable associated to the given coefficient in the linear constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));

   if( SCIPisPositive(scip, val) )
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip, consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable associated to the given coefficient in the linear constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));

   if( SCIPisPositive(scip, val) )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip, consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** creates event data for variable at given position, and catches events */
static
SCIP_RETCODE consdataCatchEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
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

   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata->eventdatas[pos]) );
   consdata->eventdatas[pos]->consdata = consdata;
   consdata->eventdatas[pos]->varpos = pos;

   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[pos],
         SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARUNLOCKED | SCIP_EVENTTYPE_GBDCHANGED,
         eventhdlr, consdata->eventdatas[pos], &consdata->eventdatas[pos]->filterpos) );

   return SCIP_OKAY;
}

/** deletes event data for variable at given position, and drops events */
static
SCIP_RETCODE consdataDropEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars[pos] != NULL);
   assert(consdata->eventdatas != NULL);
   assert(consdata->eventdatas[pos] != NULL);
   assert(consdata->eventdatas[pos]->consdata == consdata);
   assert(consdata->eventdatas[pos]->varpos == pos);

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos],
         SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARUNLOCKED | SCIP_EVENTTYPE_GBDCHANGED,
         eventhdlr, consdata->eventdatas[pos], consdata->eventdatas[pos]->filterpos) );

   SCIPfreeBlockMemory(scip, &consdata->eventdatas[pos]);

   return SCIP_OKAY;
}

/** catches bound change events for all variables in transformed linear constraint */
static
SCIP_RETCODE consdataCatchAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);
   assert(consdata->eventdatas == NULL);

   /* allocate eventdatas array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eventdatas, consdata->varssize) );
   assert(consdata->eventdatas != NULL);
   BMSclearMemoryArray(consdata->eventdatas, consdata->nvars);

   /* catch event for every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( consdataCatchEvent(scip, consdata, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for all variables in transformed linear constraint */
static
SCIP_RETCODE consdataDropAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int i;

   assert(consdata != NULL);
   assert(consdata->eventdatas != NULL);

   /* drop event of every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( consdataDropEvent(scip, consdata, eventhdlr, i) );
   }

   /* free eventdatas array */
   SCIPfreeBlockMemoryArray(scip, &consdata->eventdatas, consdata->varssize);
   assert(consdata->eventdatas == NULL);

   return SCIP_OKAY;
}

/** returns whether we are in a stage, where the variable events should be catched */
static
SCIP_Bool needEvents(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return (SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED && SCIPgetStage(scip) < SCIP_STAGE_FREETRANS);
}

/** creates a linear constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs                 /**< right hand side of row */
   )
{
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   if( SCIPisGT(scip, lhs, rhs) )
   {
      SCIPerrorMessage("left hand side of linear constraint greater than right hand side\n");
      SCIPerrorMessage(" -> lhs=%f, rhs=%f", lhs, rhs);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   if( nvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vals, vals, nvars) );
#ifndef NDEBUG
      {
         int v;
         for( v = 0; v < nvars; ++v )
            assert(vars[v] != NULL);
      }
#endif
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
   (*consdata)->minactivityneginf = -1;
   (*consdata)->minactivityposinf = -1;
   (*consdata)->maxactivityneginf = -1;
   (*consdata)->maxactivityposinf = -1;
   (*consdata)->glbminactivity = SCIP_INVALID;
   (*consdata)->glbmaxactivity = SCIP_INVALID;
   (*consdata)->glbminactivityneginf = -1;
   (*consdata)->glbminactivityposinf = -1;
   (*consdata)->glbmaxactivityneginf = -1;
   (*consdata)->glbmaxactivityposinf = -1;
   (*consdata)->possignature = 0;
   (*consdata)->negsignature = 0;
   (*consdata)->varssize = nvars;
   (*consdata)->nvars = nvars;
   (*consdata)->validmaxabsval = FALSE;
   (*consdata)->validactivities = FALSE;
   (*consdata)->propagated = FALSE;
   (*consdata)->boundstightened = FALSE;
   (*consdata)->presolved = FALSE;
   (*consdata)->removedfixings = FALSE;
   (*consdata)->validsignature = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->normalized = FALSE;
   (*consdata)->upgradetried = FALSE;
   (*consdata)->upgraded = FALSE;
   (*consdata)->sorted = (nvars <= 1);
   (*consdata)->merged = (nvars <= 1);
   (*consdata)->cliquesadded = FALSE;

   if( SCIPisTransformed(scip) )
   {
      /* get transformed variables */
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

      /* catch bound change events of variables */
      if( needEvents(scip) )
      {
         SCIP_CALL( consdataCatchAllEvents(scip, *consdata, eventhdlr) );
         assert((*consdata)->eventdatas != NULL);
      }
   }

   return SCIP_OKAY;
}

/** frees a linear constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to linear constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->varssize >= 0);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* free event datas */
   if( (*consdata)->eventdatas != NULL )
   {
      /* drop bound change events of variables */
      SCIP_CALL( consdataDropAllEvents(scip, *consdata, eventhdlr) );
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
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%g <= ", consdata->lhs);

   /* print coefficients */
   if( consdata->nvars == 0 )
      SCIPinfoMessage(scip, file, "0 ");
   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      SCIPinfoMessage(scip, file, "%+g<%s> ", consdata->vals[v], SCIPvarGetName(consdata->vars[v]));
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "== %g\n", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, "<= %g\n", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, ">= %g\n", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]\n");
}

/** updates minimum and maximum activity for a change in lower bound */
static
void consdataUpdateChgLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable that has been changed */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb,              /**< new lower bound of variable */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      assert(consdata->pseudoactivity < SCIP_INVALID);
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(consdata->minactivityneginf >= 0);
      assert(consdata->minactivityposinf >= 0);
      assert(consdata->maxactivityneginf >= 0);
      assert(consdata->maxactivityposinf >= 0);

      if( SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_LOWER )
         consdata->pseudoactivity += val * (newlb - oldlb);

      if( val > 0.0 )
      {
         if( SCIPisInfinity(scip, oldlb) )
         {
            assert(consdata->minactivityposinf >= 1);
            consdata->minactivityposinf--;
         }
         else if( SCIPisInfinity(scip, -oldlb) )
         {
            assert(consdata->minactivityneginf >= 1);
            consdata->minactivityneginf--;
         }
         else
            consdata->minactivity -= val * oldlb;

         if( SCIPisInfinity(scip, newlb) )
            consdata->minactivityposinf++;
         else if( SCIPisInfinity(scip, -newlb) )
            consdata->minactivityneginf++;
         else
            consdata->minactivity += val * newlb;
      }
      else
      {
         if( SCIPisInfinity(scip, oldlb) )
         {
            assert(consdata->maxactivityneginf >= 1);
            consdata->maxactivityneginf--;
         }
         else if( SCIPisInfinity(scip, -oldlb) )
         {
            assert(consdata->maxactivityposinf >= 1);
            consdata->maxactivityposinf--;
         }
         else
            consdata->maxactivity -= val * oldlb;

         if( SCIPisInfinity(scip, newlb) )
            consdata->maxactivityneginf++;
         else if( SCIPisInfinity(scip, -newlb) )
            consdata->maxactivityposinf++;
         else
            consdata->maxactivity += val * newlb;
      }
   }
}

/** updates minimum and maximum activity for a change in upper bound */
static
void consdataUpdateChgUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable that has been changed */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub,              /**< new upper bound of variable */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      assert(consdata->pseudoactivity < SCIP_INVALID);
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(consdata->minactivityneginf >= 0);
      assert(consdata->minactivityposinf >= 0);
      assert(consdata->maxactivityneginf >= 0);
      assert(consdata->maxactivityposinf >= 0);

      if( SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_UPPER )
         consdata->pseudoactivity += val * (newub - oldub);

      if( val > 0.0 )
      {
         if( SCIPisInfinity(scip, oldub) )
         {
            assert(consdata->maxactivityposinf >= 1);
            consdata->maxactivityposinf--;
         }
         else if( SCIPisInfinity(scip, -oldub) )
         {
            assert(consdata->maxactivityneginf >= 1);
            consdata->maxactivityneginf--;
         }
         else
            consdata->maxactivity -= val * oldub;

         if( SCIPisInfinity(scip, newub) )
            consdata->maxactivityposinf++;
         else if( SCIPisInfinity(scip, -newub) )
            consdata->maxactivityneginf++;
         else
            consdata->maxactivity += val * newub;
      }
      else
      {
         if( SCIPisInfinity(scip, oldub) )
         {
            assert(consdata->minactivityneginf >= 1);
            consdata->minactivityneginf--;
         }
         else if( SCIPisInfinity(scip, -oldub) )
         {
            assert(consdata->minactivityposinf >= 1);
            consdata->minactivityposinf--;
         }
         else
            consdata->minactivity -= val * oldub;

         if( SCIPisInfinity(scip, newub) )
            consdata->minactivityneginf++;
         else if( SCIPisInfinity(scip, -newub) )
            consdata->minactivityposinf++;
         else
            consdata->minactivity += val * newub;
      }
   }
}

/** updates minimum and maximum global activity for a change in global lower bound */
static
void consdataUpdateChgGlbLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb,              /**< new lower bound of variable */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      assert(consdata->glbminactivity < SCIP_INVALID);
      assert(consdata->glbmaxactivity < SCIP_INVALID);
      assert(consdata->glbminactivityneginf >= 0);
      assert(consdata->glbminactivityposinf >= 0);
      assert(consdata->glbmaxactivityneginf >= 0);
      assert(consdata->glbmaxactivityposinf >= 0);

      if( val > 0.0 )
      {
         if( SCIPisInfinity(scip, oldlb) )
         {
            assert(consdata->glbminactivityposinf >= 1);
            consdata->glbminactivityposinf--;
         }
         else if( SCIPisInfinity(scip, -oldlb) )
         {
            assert(consdata->glbminactivityneginf >= 1);
            consdata->glbminactivityneginf--;
         }
         else
            consdata->glbminactivity -= val * oldlb;

         if( SCIPisInfinity(scip, newlb) )
            consdata->glbminactivityposinf++;
         else if( SCIPisInfinity(scip, -newlb) )
            consdata->glbminactivityneginf++;
         else
            consdata->glbminactivity += val * newlb;
      }
      else
      {
         if( SCIPisInfinity(scip, oldlb) )
         {
            assert(consdata->glbmaxactivityneginf >= 1);
            consdata->glbmaxactivityneginf--;
         }
         else if( SCIPisInfinity(scip, -oldlb) )
         {
            assert(consdata->glbmaxactivityposinf >= 1);
            consdata->glbmaxactivityposinf--;
         }
         else
            consdata->glbmaxactivity -= val * oldlb;

         if( SCIPisInfinity(scip, newlb) )
            consdata->glbmaxactivityneginf++;
         else if( SCIPisInfinity(scip, -newlb) )
            consdata->glbmaxactivityposinf++;
         else
            consdata->glbmaxactivity += val * newlb;
      }
   }
}

/** updates minimum and maximum global activity for a change in global upper bound */
static
void consdataUpdateChgGlbUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub,              /**< new upper bound of variable */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      assert(consdata->glbminactivity < SCIP_INVALID);
      assert(consdata->glbmaxactivity < SCIP_INVALID);
      assert(consdata->glbminactivityneginf >= 0);
      assert(consdata->glbminactivityposinf >= 0);
      assert(consdata->glbmaxactivityneginf >= 0);
      assert(consdata->glbmaxactivityposinf >= 0);

      if( val > 0.0 )
      {
         if( SCIPisInfinity(scip, oldub) )
         {
            assert(consdata->glbmaxactivityposinf >= 1);
            consdata->glbmaxactivityposinf--;
         }
         else if( SCIPisInfinity(scip, -oldub) )
         {
            assert(consdata->glbmaxactivityneginf >= 1);
            consdata->glbmaxactivityneginf--;
         }
         else
            consdata->glbmaxactivity -= val * oldub;

         if( SCIPisInfinity(scip, newub) )
            consdata->glbmaxactivityposinf++;
         else if( SCIPisInfinity(scip, -newub) )
            consdata->glbmaxactivityneginf++;
         else
            consdata->glbmaxactivity += val * newub;
      }
      else
      {
         if( SCIPisInfinity(scip, oldub) )
         {
            assert(consdata->glbminactivityneginf >= 1);
            consdata->glbminactivityneginf--;
         }
         else if( SCIPisInfinity(scip, -oldub) )
         {
            assert(consdata->glbminactivityposinf >= 1);
            consdata->glbminactivityposinf--;
         }
         else
            consdata->glbminactivity -= val * oldub;

         if( SCIPisInfinity(scip, newub) )
            consdata->glbminactivityneginf++;
         else if( SCIPisInfinity(scip, -newub) )
            consdata->glbminactivityposinf++;
         else
            consdata->glbminactivity += val * newub;
      }
   }
}

/** updates minimum and maximum activity and maximum absolute value for coefficient addition */
static
void consdataUpdateAddCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(consdata != NULL);

   /* update maximum absolute value */
   if( consdata->validmaxabsval )
   {
      SCIP_Real absval;

      assert(consdata->maxabsval < SCIP_INVALID);

      absval = REALABS(val);
      consdata->maxabsval = MAX(consdata->maxabsval, absval);
   }

   /* update pseudo, minimal and maximal activity */
   if( consdata->validactivities )
   {
      assert(consdata->pseudoactivity < SCIP_INVALID);
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(consdata->glbminactivity < SCIP_INVALID);
      assert(consdata->glbmaxactivity < SCIP_INVALID);

      consdataUpdateChgLb(scip, consdata, var, 0.0, SCIPvarGetLbLocal(var), val);
      consdataUpdateChgUb(scip, consdata, var, 0.0, SCIPvarGetUbLocal(var), val);
      consdataUpdateChgGlbLb(scip, consdata, 0.0, SCIPvarGetLbGlobal(var), val);
      consdataUpdateChgGlbUb(scip, consdata, 0.0, SCIPvarGetUbGlobal(var), val);
   }
}

/** updates minimum and maximum activity for coefficient deletion, invalidates maximum absolute value if necessary */
static
void consdataUpdateDelCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
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
      assert(consdata->glbminactivity < SCIP_INVALID);
      assert(consdata->glbmaxactivity < SCIP_INVALID);

      consdataUpdateChgLb(scip, consdata, var, SCIPvarGetLbLocal(var), 0.0, val);
      consdataUpdateChgUb(scip, consdata, var, SCIPvarGetUbLocal(var), 0.0, val);
      consdataUpdateChgGlbLb(scip, consdata, SCIPvarGetLbGlobal(var), 0.0, val);
      consdataUpdateChgGlbUb(scip, consdata, SCIPvarGetUbGlobal(var), 0.0, val);
   }
}

/** updates minimum and maximum activity for coefficient change, invalidates maximum absolute value if necessary */
static
void consdataUpdateChgCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             oldval,             /**< old coefficient of constraint entry */
   SCIP_Real             newval              /**< new coefficient of constraint entry */
   )
{
   consdataUpdateDelCoef(scip, consdata, var, oldval);
   consdataUpdateAddCoef(scip, consdata, var, newval);
}

/** calculates maximum absolute value of coefficients */
static
void consdataCalcMaxAbsval(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   SCIP_Real absval;
   int i;

   assert(consdata != NULL);
   assert(!consdata->validmaxabsval);
   assert(consdata->maxabsval >= SCIP_INVALID);

   consdata->validmaxabsval = TRUE;
   consdata->maxabsval = 0.0;
   for( i = 0; i < consdata->nvars; ++i )
   {
      absval = consdata->vals[i];
      absval = REALABS(absval);
      consdata->maxabsval = MAX(consdata->maxabsval, absval);
   }
}

/** returns the maximum absolute value of all coefficients in the constraint */
static
SCIP_Real consdataGetMaxAbsval(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->validmaxabsval )
      consdataCalcMaxAbsval(consdata);
   assert(consdata->validmaxabsval);
   assert(consdata->maxabsval < SCIP_INVALID);

   return consdata->maxabsval;
}

/** calculates pseudo activity, and minimum and maximum local and global activity for constraint;
 *  additionally recalculates maximum absolute value of coefficients
 */
static
void consdataCalcActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;

   assert(consdata != NULL);
   assert(!consdata->validactivities);
   assert(consdata->pseudoactivity >= SCIP_INVALID);
   assert(consdata->minactivity >= SCIP_INVALID);
   assert(consdata->maxactivity >= SCIP_INVALID);
   assert(consdata->glbminactivity >= SCIP_INVALID);
   assert(consdata->glbmaxactivity >= SCIP_INVALID);

   consdata->validmaxabsval = TRUE;
   consdata->validactivities = TRUE;
   consdata->maxabsval = 0.0;
   consdata->pseudoactivity = 0.0;
   consdata->minactivity = 0.0;
   consdata->maxactivity = 0.0;
   consdata->minactivityneginf = 0;
   consdata->minactivityposinf = 0;
   consdata->maxactivityneginf = 0;
   consdata->maxactivityposinf = 0;
   consdata->glbminactivity = 0.0;
   consdata->glbmaxactivity = 0.0;
   consdata->glbminactivityneginf = 0;
   consdata->glbminactivityposinf = 0;
   consdata->glbmaxactivityneginf = 0;
   consdata->glbmaxactivityposinf = 0;

   for( i = 0; i < consdata->nvars; ++i )
      consdataUpdateAddCoef(scip, consdata, consdata->vars[i], consdata->vals[i]);
}

/** gets activity bounds for constraint */
static
SCIP_Real consdataGetPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint */
   )
{
   assert(consdata != NULL);

   if( !consdata->validactivities )
      consdataCalcActivities(scip, consdata);
   assert(consdata->pseudoactivity < SCIP_INVALID);
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);
   assert(consdata->glbminactivity < SCIP_INVALID);
   assert(consdata->glbmaxactivity < SCIP_INVALID);

   SCIPdebugMessage("pseudo activity of linear constraint: %g\n", consdata->pseudoactivity);

   return consdata->pseudoactivity;
}

#if 0
/** calculates the feasibility of the linear constraint for given solution */
static
SCIP_Real consdataGetPseudoFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   SCIP_Real activity;

   assert(consdata != NULL);

   activity = consdataGetPseudoActivity(scip, consdata);

   return MIN(consdata->rhs - activity, activity - consdata->lhs);
}
#endif

/** gets activity bounds for constraint */
static
void consdataGetActivityBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_Real*            minactivity,        /**< pointer to store the minimal activity */
   SCIP_Real*            maxactivity         /**< pointer to store the maximal activity */
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
   assert(consdata->minactivityneginf >= 0);
   assert(consdata->minactivityposinf >= 0);
   assert(consdata->maxactivityneginf >= 0);
   assert(consdata->maxactivityposinf >= 0);

   if( consdata->minactivityposinf > 0 )
      *minactivity = SCIPinfinity(scip);
   else if( consdata->minactivityneginf > 0 )
      *minactivity = -SCIPinfinity(scip);
   else
      *minactivity = consdata->minactivity;

   if( consdata->maxactivityneginf > 0 )
      *maxactivity = -SCIPinfinity(scip);
   else if( consdata->maxactivityposinf > 0 )
      *maxactivity = SCIPinfinity(scip);
   else
      *maxactivity = consdata->maxactivity;
}

/** gets activity bounds for constraint after setting variable to zero */
static
void consdataGetActivityResiduals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_VAR*             var,                /**< variable to calculate activity residual for */
   SCIP_Real             val,                /**< coefficient value of variable in linear constraint */
   SCIP_Real*            minresactivity,     /**< pointer to store the minimal residual activity */
   SCIP_Real*            maxresactivity      /**< pointer to store the maximal residual activity */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

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
   assert(consdata->minactivityneginf >= 0);
   assert(consdata->minactivityposinf >= 0);
   assert(consdata->maxactivityneginf >= 0);
   assert(consdata->maxactivityposinf >= 0);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   if( val > 0.0 )
   {
      if( SCIPisInfinity(scip, lb) )
      {
         assert(consdata->minactivityposinf >= 1);
         if( consdata->minactivityposinf >= 2 )
            *minresactivity = SCIPinfinity(scip);
         else if( consdata->minactivityneginf >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->minactivity;
      }
      else if( SCIPisInfinity(scip, -lb) )
      {
         assert(consdata->minactivityneginf >= 1);
         if( consdata->minactivityposinf >= 1 )
            *minresactivity = SCIPinfinity(scip);
         else if( consdata->minactivityneginf >= 2 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->minactivity;
      }
      else
      {
         if( consdata->minactivityposinf >= 1 )
            *minresactivity = SCIPinfinity(scip);
         if( consdata->minactivityneginf >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->minactivity - val * lb;
      }

      if( SCIPisInfinity(scip, -ub) )
      {
         assert(consdata->maxactivityneginf >= 1);
         if( consdata->maxactivityneginf >= 2 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->maxactivityposinf >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->maxactivity;
      }
      else if( SCIPisInfinity(scip, ub) )
      {
         assert(consdata->maxactivityposinf >= 1);
         if( consdata->maxactivityneginf >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->maxactivityposinf >= 2 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->maxactivity;
      }
      else
      {
         if( consdata->maxactivityneginf >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->maxactivityposinf >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->maxactivity - val * ub;
      }
   }
   else
   {
      if( SCIPisInfinity(scip, -ub) )
      {
         assert(consdata->minactivityposinf >= 1);
         if( consdata->minactivityposinf >= 2 )
            *minresactivity = SCIPinfinity(scip);
         else if( consdata->minactivityneginf >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->minactivity;
      }
      else if( SCIPisInfinity(scip, ub) )
      {
         assert(consdata->minactivityneginf >= 1);
         if( consdata->minactivityposinf >= 1 )
            *minresactivity = SCIPinfinity(scip);
         else if( consdata->minactivityneginf >= 2 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->minactivity;
      }
      else
      {
         if( consdata->minactivityposinf >= 1 )
            *minresactivity = SCIPinfinity(scip);
         else if( consdata->minactivityneginf >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->minactivity - val * ub;
      }

      if( SCIPisInfinity(scip, lb) )
      {
         assert(consdata->maxactivityneginf >= 1);
         if( consdata->maxactivityneginf >= 2 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->maxactivityposinf >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->maxactivity;
      }
      else if( SCIPisInfinity(scip, -lb) )
      {
         assert(consdata->maxactivityposinf >= 1);
         if( consdata->maxactivityneginf >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->maxactivityposinf >= 2 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->maxactivity;
      }
      else
      {
         if( consdata->maxactivityneginf >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->maxactivityposinf >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->maxactivity - val * lb;
      }
   }
}

/** gets global activity bounds for constraint */
static
void consdataGetGlbActivityBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_Real*            glbminactivity,     /**< pointer to store the minimal activity */
   SCIP_Real*            glbmaxactivity      /**< pointer to store the maximal activity */
   )
{
   assert(consdata != NULL);
   assert(glbminactivity != NULL);
   assert(glbmaxactivity != NULL);

   if( !consdata->validactivities )
      consdataCalcActivities(scip, consdata);
   assert(consdata->glbminactivity < SCIP_INVALID);
   assert(consdata->glbmaxactivity < SCIP_INVALID);
   assert(consdata->glbminactivityneginf >= 0);
   assert(consdata->glbminactivityposinf >= 0);
   assert(consdata->glbmaxactivityneginf >= 0);
   assert(consdata->glbmaxactivityposinf >= 0);

   if( consdata->glbminactivityposinf > 0 )
      *glbminactivity = SCIPinfinity(scip);
   else if( consdata->glbminactivityneginf > 0 )
      *glbminactivity = -SCIPinfinity(scip);
   else
      *glbminactivity = consdata->glbminactivity;

   if( consdata->glbmaxactivityneginf > 0 )
      *glbmaxactivity = -SCIPinfinity(scip);
   else if( consdata->glbmaxactivityposinf > 0 )
      *glbmaxactivity = SCIPinfinity(scip);
   else
      *glbmaxactivity = consdata->glbmaxactivity;
}

/** gets global activity bounds for constraint after setting variable to zero */
static
void consdataGetGlbActivityResiduals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_VAR*             var,                /**< variable to calculate activity residual for */
   SCIP_Real             val,                /**< coefficient value of variable in linear constraint */
   SCIP_Real*            minresactivity,     /**< pointer to store the minimal residual activity */
   SCIP_Real*            maxresactivity      /**< pointer to store the maximal residual activity */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(consdata != NULL);
   assert(var != NULL);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);

   /* get activity bounds of linear constraint */
   if( !consdata->validactivities )
      consdataCalcActivities(scip, consdata);
   assert(consdata->glbminactivity < SCIP_INVALID);
   assert(consdata->glbmaxactivity < SCIP_INVALID);
   assert(consdata->glbminactivityneginf >= 0);
   assert(consdata->glbminactivityposinf >= 0);
   assert(consdata->glbmaxactivityneginf >= 0);
   assert(consdata->glbmaxactivityposinf >= 0);

   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);

   if( val > 0.0 )
   {
      if( SCIPisInfinity(scip, lb) )
      {
         assert(consdata->glbminactivityposinf >= 1);
         if( consdata->glbminactivityposinf >= 2 )
            *minresactivity = SCIPinfinity(scip);
         else if( consdata->glbminactivityneginf >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->glbminactivity;
      }
      else if( SCIPisInfinity(scip, -lb) )
      {
         assert(consdata->glbminactivityneginf >= 1);
         if( consdata->glbminactivityposinf >= 1 )
            *minresactivity = SCIPinfinity(scip);
         else if( consdata->glbminactivityneginf >= 2 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->glbminactivity;
      }
      else
      {
         if( consdata->glbminactivityposinf >= 1 )
            *minresactivity = SCIPinfinity(scip);
         if( consdata->glbminactivityneginf >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->glbminactivity - val * lb;
      }

      if( SCIPisInfinity(scip, -ub) )
      {
         assert(consdata->glbmaxactivityneginf >= 1);
         if( consdata->glbmaxactivityneginf >= 2 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->glbmaxactivityposinf >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->glbmaxactivity;
      }
      else if( SCIPisInfinity(scip, ub) )
      {
         assert(consdata->glbmaxactivityposinf >= 1);
         if( consdata->glbmaxactivityneginf >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->glbmaxactivityposinf >= 2 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->glbmaxactivity;
      }
      else
      {
         if( consdata->glbmaxactivityneginf >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->glbmaxactivityposinf >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->glbmaxactivity - val * ub;
      }
   }
   else
   {
      if( SCIPisInfinity(scip, -ub) )
      {
         assert(consdata->glbminactivityposinf >= 1);
         if( consdata->glbminactivityposinf >= 2 )
            *minresactivity = SCIPinfinity(scip);
         else if( consdata->glbminactivityneginf >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->glbminactivity;
      }
      else if( SCIPisInfinity(scip, ub) )
      {
         assert(consdata->glbminactivityneginf >= 1);
         if( consdata->glbminactivityposinf >= 1 )
            *minresactivity = SCIPinfinity(scip);
         else if( consdata->glbminactivityneginf >= 2 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->glbminactivity;
      }
      else
      {
         if( consdata->glbminactivityposinf >= 1 )
            *minresactivity = SCIPinfinity(scip);
         else if( consdata->glbminactivityneginf >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = consdata->glbminactivity - val * ub;
      }

      if( SCIPisInfinity(scip, lb) )
      {
         assert(consdata->glbmaxactivityneginf >= 1);
         if( consdata->glbmaxactivityneginf >= 2 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->glbmaxactivityposinf >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->glbmaxactivity;
      }
      else if( SCIPisInfinity(scip, -lb) )
      {
         assert(consdata->glbmaxactivityposinf >= 1);
         if( consdata->glbmaxactivityneginf >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->glbmaxactivityposinf >= 2 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->glbmaxactivity;
      }
      else
      {
         if( consdata->glbmaxactivityneginf >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( consdata->glbmaxactivityposinf >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = consdata->glbmaxactivity - val * lb;
      }
   }
}

/** invalidates pseudo activity and activity bounds, such that they are recalculated in next get */
static
void consdataInvalidateActivities(
   SCIP_CONSDATA*        consdata            /**< linear constraint */
   )
{
   assert(consdata != NULL);

   consdata->validactivities = FALSE;
   consdata->pseudoactivity = SCIP_INVALID;
   consdata->minactivity = SCIP_INVALID;
   consdata->maxactivity = SCIP_INVALID;
   consdata->minactivityneginf = -1;
   consdata->minactivityposinf = -1;
   consdata->maxactivityneginf = -1;
   consdata->maxactivityposinf = -1;
   consdata->glbminactivity = SCIP_INVALID;
   consdata->glbmaxactivity = SCIP_INVALID;
   consdata->glbminactivityneginf = -1;
   consdata->glbminactivityposinf = -1;
   consdata->glbmaxactivityneginf = -1;
   consdata->glbmaxactivityposinf = -1;
}

/** calculates the activity of the linear constraint for given solution */
static
SCIP_Real consdataGetActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_SOL*             sol                 /**< solution to get activity for, NULL to current solution */
   )
{
   SCIP_Real activity;
   SCIP_Real scipinf;

   assert(consdata != NULL);

   if( sol == NULL && !SCIPhasCurrentNodeLP(scip) )
   {
      /* for performance reasons, the pseudo activity is updated with each bound change, so we don't have to
       * recalculate it
       */
      activity = consdataGetPseudoActivity(scip, consdata);
   }
   else
   {
      SCIP_Real solval;
      int v;

      activity = 0.0;
      for( v = 0; v < consdata->nvars; ++v )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[v]);
         activity += consdata->vals[v] * solval;
      }

      SCIPdebugMessage("activity of linear constraint: %g\n", activity);
   }

   scipinf = SCIPinfinity(scip);
   activity = MAX(activity, -scipinf);
   activity = MIN(activity, +scipinf);

   return activity;
}

/** calculates the feasibility of the linear constraint for given solution */
static
SCIP_Real consdataGetFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_SOL*             sol                 /**< solution to get feasibility for, NULL to current solution */
   )
{
   SCIP_Real activity;

   assert(consdata != NULL);

   activity = consdataGetActivity(scip, consdata, sol);

   return MIN(consdata->rhs - activity, activity - consdata->lhs);
}

/** returns the signature bitmask for the given variable */
static
SCIP_Longint getVarSignature(
   SCIP_VAR*             var                 /**< variable */
   )
{
   int sigidx;

   sigidx = SCIPvarGetIndex(var) % (int)(8*sizeof(SCIP_Longint));
   return ((SCIP_Longint)1) << sigidx;
}

/** updates bit signatures after adding a single coefficient */
static
void consdataUpdateSignatures(
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   int                   pos                 /**< position of coefficient to update signatures for */
   )
{
   SCIP_Longint varsignature;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real val;

   assert(consdata != NULL);
   assert(consdata->validsignature);

   varsignature = getVarSignature(consdata->vars[pos]);
   lb = SCIPvarGetLbGlobal(consdata->vars[pos]);
   ub = SCIPvarGetUbGlobal(consdata->vars[pos]);
   val = consdata->vals[pos];
   if( (val > 0.0 && ub > 0.0) || (val < 0.0 && lb < 0.0) )
      consdata->possignature |= varsignature;
   if( (val > 0.0 && lb < 0.0) || (val < 0.0 && ub > 0.0) )
      consdata->negsignature |= varsignature;
}

/** calculates the bit signatures of the given constraint data */
static
void consdataCalcSignatures(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->validsignature )
   {
      int i;

      consdata->validsignature = TRUE;
      consdata->possignature = 0;
      consdata->negsignature = 0;
      for( i = 0; i < consdata->nvars; ++i )
         consdataUpdateSignatures(consdata, i);
   }
}

/** index comparison method of linear constraints: compares two indices of the variable set in the linear constraint */
static
SCIP_DECL_SORTINDCOMP(consdataCompVar)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;
   int cmp;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nvars);
   assert(0 <= ind2 && ind2 < consdata->nvars);
   assert(SCIP_VARTYPE_BINARY < SCIP_VARTYPE_INTEGER);
   assert(SCIP_VARTYPE_INTEGER < SCIP_VARTYPE_IMPLINT);
   assert(SCIP_VARTYPE_IMPLINT < SCIP_VARTYPE_CONTINUOUS);

   cmp = (int)SCIPvarGetType(consdata->vars[ind1]) - (int)SCIPvarGetType(consdata->vars[ind2]);
   if( cmp == 0 )
      return SCIPvarCompare(consdata->vars[ind1], consdata->vars[ind2]);
   else
      return cmp;
}

/** sorts linear constraint's variables by binaries, integers, implicit integers, and continuous variables,
 *  and sorts the variables of the same type by non-decreasing variable index
 */
static
SCIP_RETCODE consdataSort(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->nvars == 0 )
      consdata->sorted = TRUE;
   else if( !consdata->sorted )
   {
      SCIP_VAR* varv;
      SCIP_EVENTDATA* eventdatav;
      SCIP_Real valv;
      int* perm;
      int v;
      int i;
      int nexti;

      /* get temporary memory to store the sorted permutation */
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, consdata->nvars) );

      /* call bubble sort */
      SCIPbsort((void*)consdata, consdata->nvars, consdataCompVar, perm);

      /* permute the variables in the linear constraint according to the resulting permutation */
      eventdatav = NULL;
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( perm[v] != v )
         {
            varv = consdata->vars[v];
            valv = consdata->vals[v];
            if( consdata->eventdatas != NULL )
               eventdatav = consdata->eventdatas[v];
            i = v;
            do
            {
               assert(0 <= perm[i] && perm[i] < consdata->nvars);
               assert(perm[i] != i);
               consdata->vars[i] = consdata->vars[perm[i]];
               consdata->vals[i] = consdata->vals[perm[i]];
               if( consdata->eventdatas != NULL )
               {
                  consdata->eventdatas[i] = consdata->eventdatas[perm[i]];
                  consdata->eventdatas[i]->varpos = i;
               }
               nexti = perm[i];
               perm[i] = i;
               i = nexti;
            }
            while( perm[i] != v );
            consdata->vars[i] = varv;
            consdata->vals[i] = valv;
            if( consdata->eventdatas != NULL )
            {
               consdata->eventdatas[i] = eventdatav;
               consdata->eventdatas[i]->varpos = i;
            }
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
         assert(consdata->eventdatas == NULL || consdata->eventdatas[v]->varpos == v);
      }
#endif

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &perm);
   }
   assert(consdata->sorted);

   return SCIP_OKAY;
}




/*
 * local linear constraint handler methods
 */

/** sets left hand side of linear constraint */
static
SCIP_RETCODE chgLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(!SCIPisInfinity(scip, lhs));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->vars != NULL && consdata->vals != NULL));
   assert(!SCIPisInfinity(scip, consdata->lhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->lhs, lhs) )
      return SCIP_OKAY;

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      if( SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, -lhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the left hand side switched from -infinity to a non-infinite value -> install rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
         }
      }
      else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, -lhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the left hand side switched from a non-infinte value to -infinity -> remove rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
         }
      }
   }

   /* set new left hand side */
   consdata->lhs = lhs;
   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;

   /* update the lhs of the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPchgRowLhs(scip, consdata->row, lhs) );
   }

   return SCIP_OKAY;
}

/** sets right hand side of linear constraint */
static
SCIP_RETCODE chgRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(!SCIPisInfinity(scip, -rhs));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->vars != NULL && consdata->vals != NULL));
   assert(!SCIPisInfinity(scip, -consdata->rhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->rhs, rhs) )
      return SCIP_OKAY;

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));

      if( SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, rhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the right hand side switched from infinity to a non-infinite value -> install rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
            else
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
         }
      }
      else if( !SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, rhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the right hand side switched from a non-infinte value to infinity -> remove rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
            else
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
         }
      }
   }

   /* set new right hand side */
   consdata->rhs = rhs;
   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;

   /* update the rhs of the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPchgRowRhs(scip, consdata->row, rhs) );
   }

   return SCIP_OKAY;
}

/** adds coefficient in linear constraint */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(var != NULL);
   assert(!SCIPisZero(scip, val));

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

   SCIP_CALL( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1) );
   consdata->vars[consdata->nvars] = var;
   consdata->vals[consdata->nvars] = val;
   consdata->nvars++;

   /* if we are in transformed problem, the variable needs an additional event data */
   if( transformed )
   {
      if( consdata->eventdatas != NULL )
      {
         SCIP_CONSHDLR* conshdlr;
         SCIP_CONSHDLRDATA* conshdlrdata;

         /* get event handler */
         conshdlr = SCIPconsGetHdlr(cons);
         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);
         assert(conshdlrdata->eventhdlr != NULL);

         /* initialize eventdatas array */
         consdata->eventdatas[consdata->nvars-1] = NULL;

         /* catch bound change events of variable */
         SCIP_CALL( consdataCatchEvent(scip, consdata, conshdlrdata->eventhdlr, consdata->nvars-1) );
      }

      /* update minimum and maximum activities */
      consdataUpdateAddCoef(scip, consdata, var, val);
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockRounding(scip, cons, var, val) );

   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->removedfixings = consdata->removedfixings || !SCIPvarIsActive(var);
   if( consdata->validsignature )
      consdataUpdateSignatures(consdata, consdata->nvars-1);
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;
   if( consdata->nvars == 1 )
   {
      consdata->sorted = TRUE;
      consdata->merged = TRUE;
   }
   else
   {
      consdata->sorted = consdata->sorted
         && (SCIPvarCompare(consdata->vars[consdata->nvars-2], consdata->vars[consdata->nvars-1]) == -1);
      consdata->merged = FALSE;
   }

   /* add the new coefficient to the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, val) );
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from linear constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   var = consdata->vars[pos];
   val = consdata->vals[pos];
   assert(var != NULL);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockRounding(scip, cons, var, val) );

   /* if we are in transformed problem, delete the event data of the variable */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* update minimum and maximum activities */
      consdataUpdateDelCoef(scip, consdata, var, val);

      /* drop bound change events of variable */
      if( consdata->eventdatas != NULL )
      {
         SCIP_CALL( consdataDropEvent(scip, consdata, conshdlrdata->eventhdlr, pos) );
         assert(consdata->eventdatas[pos] == NULL);
      }
   }

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->vals[pos] = consdata->vals[consdata->nvars-1];
   if( pos != consdata->nvars-1 )
   {
      if( consdata->eventdatas != NULL )
      {
         consdata->eventdatas[pos] = consdata->eventdatas[consdata->nvars-1];
         assert(consdata->eventdatas[pos] != NULL);
         consdata->eventdatas[pos]->varpos = pos;
      }
      consdata->sorted = FALSE;
   }
   consdata->nvars--;

   /* if no more variables are left, the activies should be recalculated (to give exactly 0.0) */
   if( consdata->nvars == 0 )
      consdataInvalidateActivities(consdata);

   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->validsignature = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;

   /* delete coefficient from the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, -val) );
   }

   return SCIP_OKAY;
}

/** changes coefficient value at given position of linear constraint data */
static
SCIP_RETCODE chgCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos,                /**< position of coefficient to delete */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;

   assert(!SCIPisZero(scip, newval));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(!SCIPisZero(scip, newval));

   var = consdata->vars[pos];
   val = consdata->vals[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   if( SCIPconsIsTransformed(cons) )
   {
      /* update minimum and maximum activities */
      consdataUpdateChgCoef(scip, consdata, var, val, newval);
   }

   /* if necessary, update the rounding locks of the variable */
   if( SCIPconsIsLocked(cons) && newval * val < 0.0 )
   {
      assert(SCIPconsIsTransformed(cons));

      /* remove rounding locks for variable with old coefficient */
      SCIP_CALL( unlockRounding(scip, cons, var, val) );

      /* install rounding locks for variable with new coefficient */
      SCIP_CALL( lockRounding(scip, cons, var, newval) );
   }

   /* change the value */
   consdata->vals[pos] = newval;

   consdata->propagated = FALSE;
   consdata->boundstightened = FALSE;
   consdata->presolved = FALSE;
   consdata->validsignature = consdata->validsignature && (newval * val > 0.0);
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;

   return SCIP_OKAY;
}

/** scales a linear constraint with a constant scalar */
static
SCIP_RETCODE scaleCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint to scale */
   SCIP_Real             scalar              /**< value to scale constraint with */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real newval;
   SCIP_Real absscalar;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);
   assert(!SCIPisEQ(scip, scalar, 1.0));

   /* scale the coefficients */
   for( i = 0; i < consdata->nvars; ++i )
   {
      newval = scalar * consdata->vals[i];
      if( SCIPisScalingIntegral(scip, consdata->vals[i], scalar) )
         newval = SCIPfeasFloor(scip, newval);
      if( SCIPisZero(scip, newval) )
      {
         SCIPwarningMessage("coefficient %g of variable <%s> in linear constraint <%s> scaled to zero (scalar: %g)\n",
            consdata->vals[i], SCIPvarGetName(consdata->vars[i]), SCIPconsGetName(cons), scalar);
         SCIP_CALL( delCoefPos(scip, cons, i) );
         --i;
      }
      else
         consdata->vals[i] = newval;
   }

   /* scale the sides */
   if( scalar < 0.0 )
   {
      SCIP_Real lhs;

      lhs = consdata->lhs;
      consdata->lhs = -consdata->rhs;
      consdata->rhs = -lhs;
   }
   absscalar = REALABS(scalar);
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      newval = absscalar * consdata->lhs;
      if( SCIPisScalingIntegral(scip, consdata->lhs, absscalar) )
         consdata->lhs = SCIPfeasFloor(scip, newval);
      else
         consdata->lhs = newval;
   }
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      newval = absscalar * consdata->rhs;
      if( SCIPisScalingIntegral(scip, consdata->rhs, absscalar) )
         consdata->rhs = SCIPfeasCeil(scip, newval);
      else
         consdata->rhs = newval;
   }

   consdataInvalidateActivities(consdata);
   consdata->cliquesadded = FALSE;

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
SCIP_RETCODE normalizeCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint to normalize */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Longint scm;
   SCIP_Longint nominator;
   SCIP_Longint denominator;
   SCIP_Longint gcd;
   SCIP_Longint maxmult;
   SCIP_Real epsilon;
   SCIP_Real feastol;
   SCIP_Real maxabsval;
   SCIP_Bool success;
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
    * which means, a value of feastol/epsilon should be used as maximal multiplier;
    * additionally, we don't want to scale the constraint if this would lead to too
    * large coefficients
    */
   epsilon = SCIPepsilon(scip);
   feastol = SCIPfeastol(scip);
   maxmult = (SCIP_Longint)(feastol/epsilon + feastol);
   maxabsval = consdataGetMaxAbsval(consdata);
   maxmult = MIN(maxmult, (SCIP_Longint)(MAXSCALEDCOEF/MAX(maxabsval, 1.0)));

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
      if( SCIPisGT(scip, REALABS(consdata->rhs), REALABS(consdata->lhs)) )
         mult = +1;
      else if( SCIPisLT(scip, REALABS(consdata->rhs), REALABS(consdata->lhs)) )
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
      SCIPdebugMessage("multiply linear constraint with -1.0\n");
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));
      SCIP_CALL( scaleCons(scip, cons, -1.0) );
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
         success = SCIPrealToRational(vals[i], -epsilon, epsilon, maxmult, &nominator, &denominator);
         if( success )
            scm = SCIPcalcSmaComMul(scm, denominator);
      }
   }
   assert(scm >= 1);
   success = success && (scm <= maxmult);
   if( success && scm != 1 )
   {
      /* scale the constraint with the smallest common multiple of all denominators */
      SCIPdebugMessage("scale linear constraint with %"SCIP_LONGINT_FORMAT" to make coefficients integral\n", scm);
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));
      SCIP_CALL( scaleCons(scip, cons, (SCIP_Real)scm) );
   }

   /*
    * division by greatest common divisor
    */
   if( success && nvars >= 1 )
   {
      /* all coefficients are integral: divide them by their greatest common divisor */
      assert(SCIPisIntegral(scip, vals[0]));
      gcd = (SCIP_Longint)(REALABS(vals[0]) + feastol);
      assert(gcd >= 1);
      for( i = 1; i < nvars && gcd > 1; ++i )
      {
         assert(SCIPisIntegral(scip, vals[i]));
         gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[i]) + feastol));
      }

      if( gcd > 1 )
      {
         /* divide the constaint by the greatest common divisor of the coefficients */
         SCIPdebugMessage("divide linear constraint by greatest common divisor %"SCIP_LONGINT_FORMAT"\n", gcd);
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));
         SCIP_CALL( scaleCons(scip, cons, 1.0/(SCIP_Real)gcd) );
      }
   }

   /* mark constraint to be normalized */
   consdata->normalized = TRUE;

   SCIPdebugMessage("normalized constraint:\n");
   SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));

   return SCIP_OKAY;
}

/** replaces multiple occurrences of a variable by a single coefficient */
static
SCIP_RETCODE mergeMultiples(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real valsum;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->merged )
      return SCIP_OKAY;

   /* sort the constraint */
   SCIP_CALL( consdataSort(scip, consdata) );

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
            SCIP_CALL( delCoefPos(scip, cons, v) );
            --v;
            valsum += consdata->vals[v];
         }
         while( v >= 1 && consdata->vars[v-1] == var );

         /* modify the last existing occurrence of the variable */
         assert(consdata->vars[v] == var);
         if( SCIPisZero(scip, valsum) )
         {
            SCIP_CALL( delCoefPos(scip, cons, v) );
         }
         else
         {
            SCIP_CALL( chgCoefPos(scip, cons, v, valsum) );
         }
      }
   }

   consdata->merged = TRUE;

   return SCIP_OKAY;
}

/** replaces all fixed and aggregated variables by their non-fixed counterparts */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_VAR** aggrvars;
   SCIP_Real val;
   SCIP_Real* aggrscalars;
   SCIP_Real fixedval;
   SCIP_Real aggrconst;
   int v;
   int naggrvars;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !consdata->removedfixings )
   {
      SCIPdebugMessage("applying fixings:\n");
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));

      v = 0;
      while( v < consdata->nvars )
      {
         var = consdata->vars[v];
         val = consdata->vals[v];
         assert(SCIPvarIsTransformed(var));

         switch( SCIPvarGetStatus(var) )
         {
         case SCIP_VARSTATUS_ORIGINAL:
            SCIPerrorMessage("original variable in transformed linear constraint\n");
            return SCIP_INVALIDDATA;

         case SCIP_VARSTATUS_LOOSE:
         case SCIP_VARSTATUS_COLUMN:
            ++v;
            break;

         case SCIP_VARSTATUS_FIXED:
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)));
            fixedval = SCIPvarGetLbGlobal(var);
            if( !SCIPisInfinity(scip, -consdata->lhs) )
            {
               SCIP_CALL( chgLhs(scip, cons, consdata->lhs - val * fixedval) );
            }
            if( !SCIPisInfinity(scip, consdata->rhs) )
            {
               SCIP_CALL( chgRhs(scip, cons, consdata->rhs - val * fixedval) );
            }
            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         case SCIP_VARSTATUS_AGGREGATED:
            SCIP_CALL( addCoef(scip, cons, SCIPvarGetAggrVar(var), val * SCIPvarGetAggrScalar(var)) );
            aggrconst = SCIPvarGetAggrConstant(var);
            if( !SCIPisInfinity(scip, -consdata->lhs) )
            {
               SCIP_CALL( chgLhs(scip, cons, consdata->lhs - val * aggrconst) );
            }
            if( !SCIPisInfinity(scip, consdata->rhs) )
            {
               SCIP_CALL( chgRhs(scip, cons, consdata->rhs - val * aggrconst) );
            }
            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         case SCIP_VARSTATUS_MULTAGGR:
            naggrvars = SCIPvarGetMultaggrNVars(var);
            aggrvars = SCIPvarGetMultaggrVars(var);
            aggrscalars = SCIPvarGetMultaggrScalars(var);
            for( i = 0; i < naggrvars; ++i )
            {
               SCIP_CALL( addCoef(scip, cons, aggrvars[i], val * aggrscalars[i]) );
            }
            aggrconst = SCIPvarGetMultaggrConstant(var);
            if( !SCIPisInfinity(scip, -consdata->lhs) )
            {
               SCIP_CALL( chgLhs(scip, cons, consdata->lhs - val * aggrconst) );
            }
            if( !SCIPisInfinity(scip, consdata->rhs) )
            {
               SCIP_CALL( chgRhs(scip, cons, consdata->rhs - val * aggrconst) );
            }
            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         case SCIP_VARSTATUS_NEGATED:
            SCIP_CALL( addCoef(scip, cons, SCIPvarGetNegationVar(var), -val) );
            aggrconst = SCIPvarGetNegationConstant(var);
            if( !SCIPisInfinity(scip, -consdata->lhs) )
            {
               SCIP_CALL( chgLhs(scip, cons, consdata->lhs - val * aggrconst) );
            }
            if( !SCIPisInfinity(scip, consdata->rhs) )
            {
               SCIP_CALL( chgRhs(scip, cons, consdata->rhs - val * aggrconst) );
            }
            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         default:
            SCIPerrorMessage("unknown variable status\n");
            SCIPABORT();
         }
      }
      consdata->removedfixings = TRUE;

      SCIPdebugMessage("after fixings:\n");
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));

      /* if aggregated variables have been replaced, multiple entries of the same variable are possible and we have
       * to clean up the constraint
       */
      SCIP_CALL( mergeMultiples(scip, cons) );

      SCIPdebugMessage("after merging:\n");
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));
   }
   assert(consdata->removedfixings);

#ifndef NDEBUG
   /* check, if all fixings are applied */
   for( v = 0; v < consdata->nvars; ++v )
      assert(SCIPvarIsActive(consdata->vars[v]) || SCIPvarGetStatus(consdata->vars[v]) == SCIP_VARSTATUS_MULTAGGR);
#endif

   return SCIP_OKAY;
}

/** for each variable in the linear constraint, except the inferred variable, adds one bound to the conflict analysis'
 *  candidate store (bound depends on sign of coefficient and whether the left or right hand side was the reason for the
 *  inference variable's bound change); the conflict analysis can be initialized with the linear constraint being the
 *  conflict detecting constraint by using NULL as inferred variable
 */
static
SCIP_RETCODE addConflictBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   int                   inferpos,           /**< position of the inferred variable in the vars array */
   SCIP_Bool             reasonisrhs         /**< is the right hand side responsible for the bound change? */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   vars = consdata->vars;
   vals = consdata->vals;
   nvars = consdata->nvars;
   assert(vars != NULL);
   assert(vals != NULL);
   assert(-1 <= inferpos && inferpos < nvars);
   assert((infervar == NULL) == (inferpos == -1));
   assert(inferpos == -1 || vars[inferpos] == infervar);

   /* for each variable, add the bound to the conflict queue, that is responsible for the minimal or maximal
    * residual value, depending on whether the left or right hand side is responsible for the bound change:
    *  - if the right hand side is the reason, the minimal residual activity is responsible
    *  - if the left hand side is the reason, the maximal residual activity is responsible
    */

   /* if the variable is integral we only need to add reason bounds until the propagation could be applied */
   if( infervar == NULL || SCIPvarIsIntegral(infervar) )
   {
      SCIP_Real minresactivity;
      SCIP_Real maxresactivity;

      /* calculate the minimal and maximal global activity of all other variables involved in the constraint */
      if( infervar != NULL )
         consdataGetGlbActivityResiduals(scip, consdata, infervar, vals[inferpos], &minresactivity, &maxresactivity);
      else
         consdataGetGlbActivityBounds(scip, consdata, &minresactivity, &maxresactivity);

      /* we can only do something clever, if the residual activity is finite */
      if( (reasonisrhs && !SCIPisInfinity(scip, -minresactivity))
         || (!reasonisrhs && !SCIPisInfinity(scip, maxresactivity)) )
      {
         SCIP_Real rescap;

         /* calculate the residual capacity that would be left, if the variable would be set to one more / one less
          * than its inferred bound
          */
         rescap = (reasonisrhs ? consdata->rhs - minresactivity : consdata->lhs - maxresactivity);
         if( infervar != NULL )
         {
            if( reasonisrhs == (vals[inferpos] > 0.0) )
               rescap -= vals[inferpos] * (SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) + 1.0);
            else
               rescap -= vals[inferpos] * (SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) - 1.0);
         }

         /* now add bounds as reasons until the residual capacity is exceeded */
         for( i = 0; i < nvars; ++i )
         {
            /* zero coefficients and the infered variable can be ignored */
            if( vars[i] == infervar || SCIPisZero(scip, vals[i]) )
               continue;

            /* check if the residual capacity is exceeded */
            if( (reasonisrhs && SCIPisFeasNegative(scip, rescap))
               || (!reasonisrhs && SCIPisFeasPositive(scip, rescap)) )
               break;

            /* update the residual capacity due to the local bound of this variable */
            if( reasonisrhs == (vals[i] > 0.0) )
            {
               /* rhs is reason and coeff is positive, or lhs is reason and coeff is negative -> lower bound */
               SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
               rescap -= vals[i] * (SCIPvarGetLbAtIndex(vars[i], bdchgidx, FALSE) - SCIPvarGetLbGlobal(vars[i]));
            }
            else
            {
               /* lhs is reason and coeff is positive, or rhs is reason and coeff is negative -> upper bound */
               SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
               rescap -= vals[i] * (SCIPvarGetUbAtIndex(vars[i], bdchgidx, FALSE) - SCIPvarGetUbGlobal(vars[i]));
            }
         }

         return SCIP_OKAY;
      }
   }

   /* for a bound change on a continuous variable, all locally changed bounds are responsible */
   for( i = 0; i < nvars; ++i )
   {
      /* zero coefficients and the infered variable can be ignored */
      if( vars[i] == infervar || SCIPisZero(scip, vals[i]) )
         continue;

      if( reasonisrhs == (vals[i] > 0.0) )
      {
         /* rhs is reason and coeff is positive, or lhs is reason and coeff is negative -> lower bound is responsible */
         SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
      }
      else
      {
         /* lhs is reason and coeff is positive, or rhs is reason and coeff is negative -> upper bound is responsible */
         SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
      }
   }

   return SCIP_OKAY;
}

/** resolves a propagation on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) activity residuals of all other variables tighten bounds of single variable
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   INFERINFO             inferinfo,          /**< inference information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   int inferpos;

   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   vars = consdata->vars;
   vals = consdata->vals;
   nvars = consdata->nvars;
   assert(vars != NULL);
   assert(vals != NULL);

   /* get the position of the inferred variable in the vars array */
   inferpos = inferInfoGetPos(inferinfo);
   if( inferpos >= nvars || vars[inferpos] != infervar )
   {
      /* find inference variable in constraint */
      /**@todo use a binary search here; the variables can be sorted by variable index */
      for( inferpos = 0; inferpos < nvars && vars[inferpos] != infervar; ++inferpos )
      {}
   }
   assert(inferpos < nvars);
   assert(vars[inferpos] == infervar);
   assert(!SCIPisZero(scip, vals[inferpos]));

   switch( inferInfoGetProprule(inferinfo) )
   {
   case PROPRULE_1_RHS:
      /* the bound of the variable was tightened, because the minimal or maximal residual activity of the linear
       * constraint (only taking the other variables into account) didn't leave enough space for a larger
       * domain in order to not exceed the right hand side of the inequality
       */
      assert((vals[inferpos] > 0.0) == (boundtype == SCIP_BOUNDTYPE_UPPER));
      SCIP_CALL( addConflictBounds(scip, cons, infervar, bdchgidx, inferpos, TRUE) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_1_LHS:
      /* the bound of the variable was tightened, because the minimal or maximal residual activity of the linear
       * constraint (only taking the other variables into account) didn't leave enough space for a larger
       * domain in order to not fall below the left hand side of the inequality
       */
      assert((vals[inferpos] > 0.0) == (boundtype == SCIP_BOUNDTYPE_LOWER));
      SCIP_CALL( addConflictBounds(scip, cons, infervar, bdchgidx, inferpos, FALSE) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in linear constraint <%s> at position %d for %s bound of variable <%s>\n",
         inferInfoGetProprule(inferinfo), SCIPconsGetName(cons), inferInfoGetPos(inferinfo),
         boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", SCIPvarGetName(infervar));
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** analyzes conflicting bounds on given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< conflict detecting constraint */
   SCIP_Bool             reasonisrhs         /**< is the right hand side responsible for the conflict? */
   )
{
   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* initialize conflict analysis */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );

   /* add the conflicting bound for each variable of infeasible constraint to conflict candidate queue */
   SCIP_CALL( addConflictBounds(scip, cons, NULL, NULL, -1, reasonisrhs) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

#define BOUNDSCALETOL 1e-05
/** tightens bounds of a single variable due to activity bounds */
static
SCIP_RETCODE tightenVarBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos,                /**< position of the variable in the vars array */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   /* we cannot tighten variables' bounds, if the constraint may be not complete */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   *cutoff = FALSE;

   var = consdata->vars[pos];

   /* we cannot tighten bounds of multi-aggregated variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIP_OKAY;

   val = consdata->vals[pos];
   lhs = consdata->lhs;
   rhs = consdata->rhs;
   consdataGetActivityResiduals(scip, consdata, var, val, &minresactivity, &maxresactivity);
   assert(var != NULL);
   assert(!SCIPisZero(scip, val));
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
         SCIP_Real newub;

         newub = (rhs - minresactivity)/val;
         if( (SCIPvarIsIntegral(var) && SCIPisFeasLT(scip, newub, ub))
            || SCIPisUbBetter(scip, newub, lb, ub) )
         {
            /* adjust bound for numerical reasons */
            newub = SCIPfeasCeil(scip, newub/BOUNDSCALETOL) * BOUNDSCALETOL;

            /* tighten upper bound */
            SCIPdebugMessage("linear constraint <%s>: tighten <%s>, old bds=[%g,%g], val=%g, resactivity=[%g,%g], sides=[%g,%g] -> newub=%g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newub);
            SCIP_CALL( SCIPinferVarUbCons(scip, var, newub, cons, getInferInt(PROPRULE_1_RHS, pos),
                  &infeasible, &tightened) );
            if( infeasible )
            {
               SCIPdebugMessage("linear constraint <%s>: cutoff  <%s>, new bds=[%.9f,%.9f]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

               /* analyze conflict */
               SCIP_CALL( analyzeConflict(scip, cons, TRUE) );

               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( tightened )
            {
               ub = SCIPvarGetUbLocal(var); /* get bound again: it may be additionally modified due to integrality */
               assert(SCIPisFeasLE(scip, ub, newub));
               (*nchgbds)++;
               SCIPdebugMessage("linear constraint <%s>: tighten <%s>, new bds=[%.9f,%.9f]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }
      if( !SCIPisInfinity(scip, maxresactivity) && !SCIPisInfinity(scip, -lhs) )
      {
         SCIP_Real newlb;

         newlb = (lhs - maxresactivity)/val;
         if( (SCIPvarIsIntegral(var) && SCIPisFeasGT(scip, newlb, lb))
            || SCIPisLbBetter(scip, newlb, lb, ub) )
         {
            /* adjust bound for numerical reasons */
            newlb = SCIPfeasFloor(scip, newlb/BOUNDSCALETOL) * BOUNDSCALETOL;

            /* tighten lower bound */
            SCIPdebugMessage("linear constraint <%s>: tighten <%s>, old bds=[%g,%g], val=%g, resactivity=[%g,%g], sides=[%g,%g] -> newlb=%g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newlb);
            SCIP_CALL( SCIPinferVarLbCons(scip, var, newlb, cons, getInferInt(PROPRULE_1_LHS, pos),
                  &infeasible, &tightened) );
            if( infeasible )
            {
               SCIPdebugMessage("linear constraint <%s>: cutoff  <%s>, new bds=[%.9f,%.9f]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

               /* analyze conflict */
               SCIP_CALL( analyzeConflict(scip, cons, FALSE) );

               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( tightened )
            {
               lb = SCIPvarGetLbLocal(var); /* get bound again: it may be additionally modified due to integrality */
               assert(SCIPisFeasGE(scip, lb, newlb));
               (*nchgbds)++;
               SCIPdebugMessage("linear constraint <%s>: tighten <%s>, new bds=[%.9f,%.9f]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }
   }
   else
   {
      /* check, if we can tighten the variable's bounds */
      if( !SCIPisInfinity(scip, -minresactivity) && !SCIPisInfinity(scip, rhs) )
      {
         SCIP_Real newlb;

         newlb = (rhs - minresactivity)/val;
         if( SCIPisLbBetter(scip, newlb, lb, ub) )
         {
            /* adjust bound for numerical reasons */
            newlb = SCIPfeasFloor(scip, newlb/BOUNDSCALETOL) * BOUNDSCALETOL;

            /* tighten lower bound */
            SCIPdebugMessage("linear constraint <%s>: tighten <%s>, old bds=[%g,%g], val=%g, resactivity=[%g,%g], sides=[%g,%g] -> newlb=%g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newlb);
            SCIP_CALL( SCIPinferVarLbCons(scip, var, newlb, cons, getInferInt(PROPRULE_1_RHS, pos),
                  &infeasible, &tightened) );
            if( infeasible )
            {
               SCIPdebugMessage("linear constraint <%s>: cutoff  <%s>, new bds=[%.9f,%.9f]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

               /* analyze conflict */
               SCIP_CALL( analyzeConflict(scip, cons, TRUE) );

               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( tightened )
            {
               lb = SCIPvarGetLbLocal(var); /* get bound again: it may be additionally modified due to integrality */
               assert(SCIPisFeasGE(scip, lb, newlb));
               (*nchgbds)++;
               SCIPdebugMessage("linear constraint <%s>: tighten <%s>, new bds=[%.9f,%.9f]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }
      if( !SCIPisInfinity(scip, maxresactivity) && !SCIPisInfinity(scip, -lhs) )
      {
         SCIP_Real newub;

         newub = (lhs - maxresactivity)/val;
         if( SCIPisUbBetter(scip, newub, lb, ub) )
         {
            /* adjust bound for numerical reasons */
            newub = SCIPfeasCeil(scip, newub/BOUNDSCALETOL) * BOUNDSCALETOL;

            /* tighten upper bound */
            SCIPdebugMessage("linear constraint <%s>: tighten <%s>, old bds=[%g,%g], val=%g, resactivity=[%g,%g], sides=[%g,%g], newub=%g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newub);
            SCIP_CALL( SCIPinferVarUbCons(scip, var, newub, cons, getInferInt(PROPRULE_1_LHS, pos),
                  &infeasible, &tightened) );
            if( infeasible )
            {
               SCIPdebugMessage("linear constraint <%s>: cutoff  <%s>, new bds=[%.9f,%.9f]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

               /* analyze conflict */
               SCIP_CALL( analyzeConflict(scip, cons, FALSE) );

               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( tightened )
            {
               ub = SCIPvarGetUbLocal(var); /* get bound again: it may be additionally modified due to integrality */
               assert(SCIPisFeasLE(scip, ub, newub));
               (*nchgbds)++;
               SCIPdebugMessage("linear constraint <%s>: tighten <%s>, new bds=[%.9f,%.9f]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }
   }

   return SCIP_OKAY;
}

#define MAXTIGHTENROUNDS 10
/** tightens bounds of variables in constraint due to activity bounds */
static
SCIP_RETCODE tightenBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   int nvars;
   int nrounds;
   int lastchange;

   assert(nchgbds != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* we cannot tighten variables' bounds, if the constraint may be not complete */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* as long as the bounds might be tightened again, try to tighten them; abort after a maximal number of rounds */
   lastchange = -1;
   for( nrounds = 0; !consdata->boundstightened && nrounds < MAXTIGHTENROUNDS; ++nrounds )
   {
      int v;

      /* mark the constraint to have the variables' bounds tightened */
      consdata->boundstightened = TRUE;

      /* try to tighten the bounds of each variable in the constraint */
      for( v = 0; v < nvars && v != lastchange && !(*cutoff); ++v )
      {
         int oldnchgbds;

         oldnchgbds = *nchgbds;
         SCIP_CALL( tightenVarBounds(scip, cons, v, cutoff, nchgbds) );
         if( *nchgbds > oldnchgbds )
            lastchange = v;
      }
   }

   return SCIP_OKAY;
}

/** checks linear constraint for feasibility of given solution or current solution */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_SOL*             sol,                /**< solution to be checked, or NULL for current solution */
   SCIP_Bool             checklprows,        /**< has linear constraint to be checked, if it is already in current LP? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real activity;

   assert(violated != NULL);

   SCIPdebugMessage("checking linear constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   if( consdata->row != NULL )
   {
      if( !checklprows && SCIProwIsInLP(consdata->row) )
         return SCIP_OKAY;
      else if( sol == NULL && !SCIPhasCurrentNodeLP(scip) )
         activity = consdataGetPseudoActivity(scip, consdata);
      else
         activity = SCIPgetRowSolActivity(scip, consdata->row, sol);
   }
   else
      activity = consdataGetActivity(scip, consdata, sol);

   SCIPdebugMessage("  consdata activity=%g (lhs=%g, rhs=%g, row=%p, checklprows=%d, rowinlp=%d, sol=%p, hascurrentnodelp=%d)\n",
      activity, consdata->lhs, consdata->rhs, consdata->row, checklprows,
      consdata->row == NULL ? 0 : SCIProwIsInLP(consdata->row), sol,
      consdata->row == NULL ? FALSE : SCIPhasCurrentNodeLP(scip));

   if( SCIPisFeasLT(scip, activity, consdata->lhs) || SCIPisFeasGT(scip, activity, consdata->rhs) )
   {
      *violated = TRUE;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }
   else
   {
      *violated = FALSE;
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** creates an LP row in a linear constraint data */
static
SCIP_RETCODE createRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), consdata->lhs, consdata->rhs,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   SCIP_CALL( SCIPaddVarsToRow(scip, consdata->row, consdata->nvars, consdata->vars, consdata->vals) );

   return SCIP_OKAY;
}

/** adds linear constraint as cut to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_SOL*             sol                 /**< primal CIP solution, NULL for current LP solution */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      /* convert consdata object into LP row */
      SCIP_CALL( createRow(scip, cons) );
   }
   assert(consdata->row != NULL);

   /* insert LP row as cut */
   if( !SCIProwIsInLP(consdata->row) )
   {
      SCIPdebugMessage("adding relaxation of linear constraint <%s>: ", SCIPconsGetName(cons));
      SCIPdebug( SCIProwPrint(consdata->row, NULL) );
      SCIP_CALL( SCIPaddCut(scip, sol, consdata->row, FALSE) );
   }

   return SCIP_OKAY;
}

/* separates relaxed knapsack constraint */
static
SCIP_RETCODE separateRelaxedKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   nknapvars,          /**< number of variables in the continuous knapsack constraint */
   SCIP_VAR**            knapvars,           /**< variables in the continuous knapsack constraint */
   SCIP_Real*            knapvals,           /**< coefficientce of the variables in the continuous knapsack constraint */
   SCIP_Real             valscale,           /**< -1.0 if lhs of row is used as rhs of c. k. constraint, +1.0 otherwise */
   SCIP_Real             rhs,                /**< right hand side of the continuous knapsack constraint */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_VAR** binvars;
   SCIP_VAR** consvars;
   SCIP_Real* binvals;
   SCIP_Longint* consvals;
   SCIP_Longint maxact;
   SCIP_Real intscalar;
   SCIP_Bool success;
   int nbinvars;
   int nconsvars;
   int i;

   assert(nknapvars > 0);
   assert(knapvars != NULL);

   SCIPdebugMessage("separate linear constraint <%s> relaxed to knapsack\n", SCIPconsGetName(cons));
   SCIPdebug(SCIPprintCons(scip, cons, NULL));

   SCIP_CALL( SCIPgetVarsData(scip, &binvars, NULL, &nbinvars, NULL, NULL, NULL) );

   if( nbinvars == 0 )
      return SCIP_OKAY;

   /* set up data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvals, nbinvars) );
   BMSclearMemoryArray(binvals, nbinvars);

   /* relax continuous knapsack constraint:
    * 1. make all variables binary:
    *    if x_j is continuous or integer variable substitute:
    *      - a_j < 0: x_j = lb  or  x_j = b*z + d with variable lower bound b*z + d with binary variable z
    *      - a_j > 0: x_j = ub  or  x_j = b*z + d with variable upper bound b*z + d with binary variable z
    * 2. convert coefficients of all variables to positive integers:
    *      - scale all coefficients a_j to a~_j integral
    *      - substitute  x~_j = 1 - x_j if a~_j < 0
    */

   /* replace integer and continuous variables with binary variables */
   for( i = 0; i < nknapvars; i++ )
   {
      SCIP_VAR* var;

      var = knapvars[i];

      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      {
         assert(0 <= SCIPvarGetProbindex(var) && SCIPvarGetProbindex(var) < nbinvars);
         binvals[SCIPvarGetProbindex(var)] += valscale * knapvals[i];
         SCIPdebugMessage(" -> binary variable %+g<%s>(%g)\n", 
            valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var));
      }
      else if( valscale * knapvals[i] > 0.0 )
      {
         SCIP_VAR** zvlb;
         SCIP_Real* bvlb;
         SCIP_Real* dvlb;
         SCIP_Real bestlbsol;
         int bestlbtype;
         int nvlb;
         int j;

         /* a_j > 0: substitution with lb or vlb */
         nvlb = SCIPvarGetNVlbs(var);
         zvlb = SCIPvarGetVlbVars(var);
         bvlb = SCIPvarGetVlbCoefs(var);
         dvlb = SCIPvarGetVlbConstants(var);

         /* search for lb or vlb with maximal bound value */
         bestlbsol = SCIPvarGetLbGlobal(var);
         bestlbtype = -1;
         for( j = 0; j < nvlb; j++ )
         {
            /* use only vlb with binary variable z */
            if( SCIPvarGetType(zvlb[j]) == SCIP_VARTYPE_BINARY && SCIPvarIsActive(zvlb[j]) )
            {
               SCIP_Real vlbsol;

               assert(0 <= SCIPvarGetProbindex(zvlb[j]) && SCIPvarGetProbindex(zvlb[j]) < nbinvars);
               vlbsol = bvlb[j] * SCIPgetSolVal(scip, sol, zvlb[j]) + dvlb[j];
               if( SCIPisGE(scip, vlbsol, bestlbsol) )
               {
                  bestlbsol = vlbsol;
                  bestlbtype = j;
               }
            }
         }

         /* if no lb or vlb with binary variable was found, we have to abort */
         if( SCIPisInfinity(scip, -bestlbsol) )
            goto TERMINATE;

         if( bestlbtype == -1 )
         {
            rhs -= valscale * knapvals[i] * bestlbsol;
            SCIPdebugMessage(" -> non-binary variable %+g<%s>(%g) replaced with lower bound %g (rhs=%g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetLbGlobal(var), rhs);
         }
         else
         {
            assert(0 <= SCIPvarGetProbindex(zvlb[bestlbtype]) && SCIPvarGetProbindex(zvlb[bestlbtype]) < nbinvars);
            rhs -= valscale * knapvals[i] * dvlb[bestlbtype];
            binvals[SCIPvarGetProbindex(zvlb[bestlbtype])] += valscale * knapvals[i] * bvlb[bestlbtype];
            SCIPdebugMessage(" -> non-binary variable %+g<%s>(%g) replaced with variable lower bound %+g<%s>(%g) %+g (rhs=%g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var),
               bvlb[bestlbtype], SCIPvarGetName(zvlb[bestlbtype]),
               SCIPgetSolVal(scip, sol, zvlb[bestlbtype]), dvlb[bestlbtype], rhs);
         }
      }
      else
      {
         SCIP_VAR** zvub;
         SCIP_Real* bvub;
         SCIP_Real* dvub;
         SCIP_Real bestubsol;
         int bestubtype;
         int nvub;
         int j;

         assert(valscale * knapvals[i] < 0.0);

         /* a_j < 0: substitution with ub or vub */
         nvub = SCIPvarGetNVubs(var);
         zvub = SCIPvarGetVubVars(var);
         bvub = SCIPvarGetVubCoefs(var);
         dvub = SCIPvarGetVubConstants(var);

         /* search for ub or vub with minimal bound value */
         bestubsol = SCIPvarGetUbGlobal(var);
         bestubtype = -1;
         for( j = 0; j < nvub; j++ )
         {
            /* use only vub with active binary variable z */
            if( SCIPvarGetType(zvub[j]) == SCIP_VARTYPE_BINARY && SCIPvarIsActive(zvub[j])  )
            {
               SCIP_Real vubsol;

               assert(0 <= SCIPvarGetProbindex(zvub[j]) && SCIPvarGetProbindex(zvub[j]) < nbinvars);
               vubsol = bvub[j] * SCIPgetSolVal(scip, sol, zvub[j]) + dvub[j];
               if( SCIPisLE(scip, vubsol, bestubsol) )
               {
                  bestubsol = vubsol;
                  bestubtype = j;
               }
            }
         }

         /* if no ub or vub with binary variable was found, we have to abort */
         if( SCIPisInfinity(scip, bestubsol) )
            goto TERMINATE;

         if( bestubtype == -1 )
         {
            rhs -= valscale * knapvals[i] * bestubsol;
            SCIPdebugMessage(" -> non-binary variable %+g<%s>(%g) replaced with upper bound %g (rhs=%g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetUbGlobal(var), rhs);
         }
         else
         {
            assert(0 <= SCIPvarGetProbindex(zvub[bestubtype]) && SCIPvarGetProbindex(zvub[bestubtype]) < nbinvars);
            rhs -= valscale * knapvals[i] * dvub[bestubtype];
            binvals[SCIPvarGetProbindex(zvub[bestubtype])] += valscale * knapvals[i] * bvub[bestubtype];
            SCIPdebugMessage(" -> non-binary variable %+g<%s>(%g) replaced with variable upper bound %+g<%s>(%g) %+g (rhs=%g)\n",
               valscale * knapvals[i], SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var),
               bvub[bestubtype], SCIPvarGetName(zvub[bestubtype]),
               SCIPgetSolVal(scip, sol, zvub[bestubtype]), dvub[bestubtype], rhs);
         }
      }
   }

   /* convert coefficents of all (now binary) variables to positive integers:
    *   - make all coefficients integral
    *   - make all coefficients positive (substitute negated variable)
    */
   nconsvars = 0;

   /* calculate scalar which makes all coefficients integral */
   SCIP_CALL( SCIPcalcIntegralScalar(binvals, nbinvars, -SCIPepsilon(scip), KNAPSACKRELAX_MAXDELTA,
         KNAPSACKRELAX_MAXDNOM, KNAPSACKRELAX_MAXSCALE, &intscalar, &success) );
   SCIPdebugMessage(" -> intscalar = %g\n", intscalar);

   /* if coefficients can not be made integral, we have to use a scalar of 1.0 and only round fractional coefficients down */
   if( !success )
      intscalar = 1.0;

   /* make all coefficients integral and positive:
    *  - scale a~_j = a_j * intscalar
    *  - substitute x~_j = 1 - x_j if a~_j < 0
    */
   rhs = rhs*intscalar;

   SCIPdebugMessage(" -> rhs = %g\n", rhs);
   maxact = 0;
   for( i = 0; i < nbinvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_Longint val;

      val = (SCIP_Longint)SCIPfloor(scip, binvals[i]*intscalar);
      if( val == 0 )
         continue;

      if( val > 0 )
      {
         var = binvars[i];
         SCIPdebugMessage(" -> positive scaled binary variable %+"SCIP_LONGINT_FORMAT"<%s> (unscaled %g): not changed (rhs=%g)\n",
            val, SCIPvarGetName(var), binvals[i], rhs);
      }
      else
      {
         assert(val < 0);

         SCIP_CALL( SCIPgetNegatedVar(scip, binvars[i], &var) );
         val = -val;
         rhs += val;
         SCIPdebugMessage(" -> negative scaled binary variable %+"SCIP_LONGINT_FORMAT"<%s> (unscaled %g): substituted by (1 - <%s>) (rhs=%g)\n",
            -val, SCIPvarGetName(binvars[i]), binvals[i], SCIPvarGetName(var), rhs);
      }

      maxact += val;
      consvals[nconsvars] = val;
      consvars[nconsvars] = var;
      nconsvars++;
   }

   if( nconsvars > 0 )
   {
      SCIP_Longint capacity;

      assert(consvars != NULL);
      assert(consvals != NULL);
      capacity = (SCIP_Longint)SCIPfeasFloor(scip, rhs);
      if( maxact > capacity )
      {
#ifdef SCIP_DEBUG
         SCIP_Real act;

         SCIPdebugMessage(" -> linear constraint <%s> relaxed to knapsack:", SCIPconsGetName(cons));
         act = 0.0;
         for( i = 0; i < nconsvars; ++i )
         {
            SCIPdebugPrintf(" %+"SCIP_LONGINT_FORMAT"<%s>(%g)", consvals[i], SCIPvarGetName(consvars[i]),
               SCIPgetSolVal(scip, sol, consvars[i]));
            act += consvals[i] * SCIPgetSolVal(scip, sol, consvars[i]);
         }
         SCIPdebugPrintf(" <= %"SCIP_LONGINT_FORMAT" (%g) [act: %g, max: %"SCIP_LONGINT_FORMAT"]\n",
            capacity, rhs, act, maxact);
#endif

         /* separate lifted cut from relaxed knapsack constraint */
         SCIP_CALL( SCIPseparateKnapsackCover(scip, cons, consvars, nconsvars, consvals, capacity, sol, -1, ncuts) );
      }
   }

 TERMINATE:
   /* free data structures */
   SCIPfreeBufferArray(scip, &binvals);
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}

/** separates linear constraint: adds linear constraint as cut, if violated by given solution */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool             separateall,        /**< should all constraints be subject to cover cut generation instead of only
                                              *   the ones with non-zero dual value? */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool violated;
   int oldncuts;

   consdata = SCIPconsGetData(cons);
   assert(ncuts != NULL);
   assert(consdata != NULL);
   assert(cons != NULL);

   oldncuts = *ncuts;

   SCIP_CALL( checkCons(scip, cons, sol, (sol != NULL), &violated) );

   if( violated )
   {
      /* insert LP row as cut */
      SCIP_CALL( addRelaxation(scip, cons, sol) );
      (*ncuts)++;
   }
   else if( !SCIPconsIsModifiable(cons) )
   {
      /* relax linear constraint into knapsack constraint and separate lifted cardinality cuts */
      if( !separateall && sol == NULL )
      {
         /* we only want to call the knapsack cover separator for rows that have a non-zero dual solution */
         if( consdata->row != NULL && SCIProwIsInLP(consdata->row) )
         {
            SCIP_Real dualsol;

            dualsol = SCIProwGetDualsol(consdata->row);
            if( SCIPisFeasNegative(scip, dualsol) )
            {
               if( !SCIPisInfinity(scip, consdata->rhs) )
               {
                  SCIP_CALL( separateRelaxedKnapsack(scip, cons, consdata->nvars, consdata->vars,
                        consdata->vals, +1.0, consdata->rhs, sol, ncuts) );
               }
            }
            else if( SCIPisFeasPositive(scip, dualsol) )
            {
               if( !SCIPisInfinity(scip, -consdata->lhs) )
               {
                  SCIP_CALL( separateRelaxedKnapsack(scip, cons, consdata->nvars, consdata->vars,
                        consdata->vals, -1.0, -consdata->lhs, sol, ncuts) );
               }
            }
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIP_CALL( separateRelaxedKnapsack(scip, cons, consdata->nvars, consdata->vars,
                  consdata->vals, +1.0, consdata->rhs, sol, ncuts) );
         }
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_CALL( separateRelaxedKnapsack(scip, cons, consdata->nvars, consdata->vars,
                  consdata->vals, -1.0, -consdata->lhs, sol, ncuts) );
         }
      }
   }

   if( *ncuts > oldncuts )
   {
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** propagation method for linear constraints */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool             tightenbounds,      /**< should the variable's bounds be tightened? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;

   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   /*SCIPdebugMessage("propagating linear constraint <%s>\n", SCIPconsGetName(cons));*/

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *cutoff = FALSE;

   /* check, if constraint is already propagated */
   if( consdata->propagated && (!tightenbounds || consdata->boundstightened) )
      return SCIP_OKAY;

   /* mark constraint to be propagated */
   consdata->propagated = TRUE;

   /* we can only infer activity bounds of the linear constraint, if it is not modifiable */
   if( !SCIPconsIsModifiable(cons) )
   {
      /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
      if( !SCIPinRepropagation(scip) )
      {
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }

      /* tighten the variable's bounds */
      if( tightenbounds )
      {
         int oldnchgbds;

         oldnchgbds = *nchgbds;
         SCIP_CALL( tightenBounds(scip, cons, cutoff, nchgbds) );
         if( *nchgbds > oldnchgbds )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
      }

      /* check constraint for infeasibility and redundancy */
      if( !(*cutoff) )
      {
         consdataGetActivityBounds(scip, consdata, &minactivity, &maxactivity);

         if( SCIPisFeasGT(scip, minactivity, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is infeasible (rhs): activitybounds=[%g,%g], sides=[%g,%g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);

            /* analyze conflict */
            SCIP_CALL( analyzeConflict(scip, cons, TRUE) );

            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( SCIPisFeasLT(scip, maxactivity, consdata->lhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is infeasible (lhs): activitybounds=[%g,%g], sides=[%g,%g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);

            /* analyze conflict */
            SCIP_CALL( analyzeConflict(scip, cons, FALSE) );

            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( SCIPisGE(scip, minactivity, consdata->lhs) && SCIPisLE(scip, maxactivity, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is redundant: activitybounds=[%g,%g], sides=[%g,%g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         }
      }
   }

   return SCIP_OKAY;
}




/*
 * Presolving methods
 */

/** extracts cliques of the constraint and adds them to SCIP */
static
SCIP_RETCODE extractCliques(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool lhsclique;
   SCIP_Bool rhsclique;
   int i;
   int nposcoefs;
   int nnegcoefs;

   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check if we already added the cliques of the constraints */
   if( consdata->cliquesadded )
      return SCIP_OKAY;

   consdata->cliquesadded = TRUE;

   /* sort variables by variable type */
   SCIP_CALL( consdataSort(scip, consdata) );

   /* currently, we only check whether the constraint is a set packing / partitioning constraint */
   /**@todo extract more cliques from linear constraints */
   if( SCIPvarGetType(consdata->vars[consdata->nvars-1]) != SCIP_VARTYPE_BINARY )
      return SCIP_OKAY;

   /* all variables are binary: check, if the coefficients are +1 or -1, and if the right hand side is equal
    * to 1 - number of negative coefficients, or if the left hand side is equal to number of positive coefficients - 1
    */
   nposcoefs = 0;
   nnegcoefs = 0;
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( SCIPisEQ(scip, consdata->vals[i], +1.0) )
         nposcoefs++;
      else if( SCIPisEQ(scip, consdata->vals[i], -1.0) )
         nnegcoefs++;
      else
         return SCIP_OKAY;
   }
   lhsclique = SCIPisEQ(scip, consdata->lhs, (SCIP_Real)nposcoefs - 1.0);
   rhsclique = SCIPisEQ(scip, consdata->rhs, 1.0 - (SCIP_Real)nnegcoefs);
   if( lhsclique || rhsclique )
   {
      SCIP_Bool* values;
      SCIP_Bool infeasible;
      int nbdchgs;

      SCIPdebugMessage("linear constraint <%s>: adding clique with %d vars (%d pos, %d neg)\n",
         SCIPconsGetName(cons), consdata->nvars, nposcoefs, nnegcoefs);
      SCIP_CALL( SCIPallocBufferArray(scip, &values, consdata->nvars) );
      for( i = 0; i < consdata->nvars; ++i )
         values[i] = (rhsclique == (consdata->vals[i] > 0.0));
      SCIP_CALL( SCIPaddClique(scip, consdata->vars, values, consdata->nvars, &infeasible, &nbdchgs) );
      if( infeasible )
         *cutoff = TRUE;
      *nchgbds += nbdchgs;
      SCIPfreeBufferArray(scip, &values);
   }

   return SCIP_OKAY;
}

/** tightens left and right hand side of constraint due to integrality */
static
SCIP_RETCODE tightenSides(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool integral;
   int i;

   assert(nchgsides != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPisIntegral(scip, consdata->lhs) || !SCIPisIntegral(scip, consdata->rhs) )
   {
      integral = TRUE;
      for( i = 0; i < consdata->nvars && integral; ++i )
      {
         integral = SCIPisIntegral(scip, consdata->vals[i])
                    && (SCIPvarGetType(consdata->vars[i]) != SCIP_VARTYPE_CONTINUOUS);
      }
      if( integral )
      {
         SCIPdebugMessage("linear constraint <%s>: make sides integral: sides=[%g,%g]\n",
            SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
         if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisIntegral(scip, consdata->lhs) )
         {
            SCIP_CALL( chgLhs(scip, cons, SCIPfeasCeil(scip, consdata->lhs)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisIntegral(scip, consdata->rhs) )
         {
            SCIP_CALL( chgRhs(scip, cons, SCIPfeasFloor(scip, consdata->rhs)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
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
SCIP_RETCODE consdataTightenCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Real val;
   SCIP_Real newval;
   SCIP_Real newlhs;
   SCIP_Real newrhs;
   SCIP_Real lb;
   SCIP_Real ub;
   int i;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);

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
               SCIPdebugMessage("linear constraint <%s>: change coefficient %+g<%s> to %+g<%s>, act=[%g,%g], side=[%g,%g]\n",
                  SCIPconsGetName(cons), val, SCIPvarGetName(var), newval, SCIPvarGetName(var),
                  minactivity, maxactivity, consdata->lhs, consdata->rhs);

               /* update the coefficient and the activity bounds */
               if( SCIPisZero(scip, newval) )
               {
                  SCIP_CALL( delCoefPos(scip, cons, i) );
                  i--;
               }
               else
               {
                  SCIP_CALL( chgCoefPos(scip, cons, i, newval) );
               }
               (*nchgcoefs)++;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, &minactivity, &maxactivity);
            }
            if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisEQ(scip, newlhs, consdata->lhs) )
            {
               SCIPdebugMessage("linear constraint <%s>: change lhs %g to %g\n", SCIPconsGetName(cons), consdata->lhs, newlhs);

               SCIP_CALL( chgLhs(scip, cons, newlhs) );
               (*nchgsides)++;
               assert(SCIPisEQ(scip, consdata->lhs, newlhs));
            }
            if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisEQ(scip, newrhs, consdata->rhs) )
            {
               SCIPdebugMessage("linear constraint <%s>: change rhs %g to %g\n", SCIPconsGetName(cons), consdata->rhs, newrhs);

               SCIP_CALL( chgRhs(scip, cons, newrhs) );
               (*nchgsides)++;
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
               SCIPdebugMessage("linear constraint <%s>: change coefficient %+g<%s> to %+g<%s>, act=[%g,%g], side=[%g,%g]\n",
                  SCIPconsGetName(cons), val, SCIPvarGetName(var), newval, SCIPvarGetName(var),
                  minactivity, maxactivity, consdata->lhs, consdata->rhs);

               /* update the coefficient and the activity bounds */
               if( SCIPisZero(scip, newval) )
               {
                  SCIP_CALL( delCoefPos(scip, cons, i) );
                  i--;
               }
               else
               {
                  SCIP_CALL( chgCoefPos(scip, cons, i, newval) );
               }
               (*nchgcoefs)++;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, &minactivity, &maxactivity);
            }
            if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisEQ(scip, newlhs, consdata->lhs) )
            {
               SCIPdebugMessage("linear constraint <%s>: change lhs %g to %g\n", SCIPconsGetName(cons), consdata->lhs, newlhs);

               SCIP_CALL( chgLhs(scip, cons, newlhs) );
               (*nchgsides)++;
               assert(SCIPisEQ(scip, consdata->lhs, newlhs));
            }
            if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisEQ(scip, newrhs, consdata->rhs) )
            {
               SCIPdebugMessage("linear constraint <%s>: change rhs %g to %g\n", SCIPconsGetName(cons), consdata->rhs, newrhs);

               SCIP_CALL( chgRhs(scip, cons, newrhs) );
               (*nchgsides)++;
               assert(SCIPisEQ(scip, consdata->rhs, newrhs));
            }
         }
      }
   }

   return SCIP_OKAY;
}

/* processes equality with only one variable by fixing the variable and deleting the constraint */
static
SCIP_RETCODE convertUnaryEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real fixval;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 1);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   /* calculate the value to fix the variable to */
   var = consdata->vars[0];
   val = consdata->vals[0];
   assert(!SCIPisZero(scip, val));
   fixval = SCIPselectSimpleValue(consdata->lhs/val - SCIPepsilon(scip), consdata->rhs/val + SCIPepsilon(scip), MAXDNOM);
   SCIPdebugMessage("linear equality <%s>: fix <%s> == %g\n",
      SCIPconsGetName(cons), SCIPvarGetName(var), fixval);

   /* fix variable */
   SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
   if( infeasible )
   {
      SCIPdebugMessage(" -> infeasible fixing\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }
   if( fixed )
      (*nfixedvars)++;

   /* disable constraint */
   SCIP_CALL( SCIPdelCons(scip, cons) );
   if( !consdata->upgraded )
      (*ndelconss)++;

   return SCIP_OKAY;
}

/* processes equality with exactly two variables by aggregating one of the variables and deleting the constraint */
static
SCIP_RETCODE convertBinaryEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;

   assert(cutoff != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 2);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   SCIPdebugMessage("linear constraint <%s>: aggregate %g<%s> + %g<%s> == %g\n",
      SCIPconsGetName(cons), consdata->vals[0], SCIPvarGetName(consdata->vars[0]),
      consdata->vals[1], SCIPvarGetName(consdata->vars[1]), consdata->rhs);

   /* aggregate the equality */
   SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], consdata->vals[0], consdata->vals[1],
         consdata->rhs, &infeasible, &redundant, &aggregated) );

   /* check for infeasibility of aggregation */
   if( infeasible )
   {
      SCIPdebugMessage(" -> infeasible aggregation\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* count the aggregation */
   if( aggregated )
      (*naggrvars)++;

   /* delete the constraint, if it is redundant */
   if( redundant )
   {
      SCIP_CALL( SCIPdelCons(scip, cons) );

      if( !consdata->upgraded )
         (*ndelconss)++;
   }

   return SCIP_OKAY;
}

/* processes equality with more than two variables by multi-aggregating one of the variables and converting the equality
 * into an inequality; if multi-aggregation is not possible, tries to identify one continuous or integer variable that is
 * implicitly integral by this constraint
 */
static
SCIP_RETCODE convertLongEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  naggrvars           /**< pointer to count number of aggregated variables */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_VARTYPE bestslacktype;
   SCIP_VARTYPE slacktype;
   SCIP_Real bestslackdomrng;
   SCIP_Bool coefszeroone;
   SCIP_Bool coefsintegral;
   SCIP_Bool varsintegral;
   int bestslackpos;
   int ncontvars;
   int contvarpos;
   int nintvars;
   int intvarpos;
   int v;

   assert(cutoff != NULL);
   assert(naggrvars != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 2);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   SCIPdebugMessage("linear constraint <%s>: try to multi-aggregate equality\n", SCIPconsGetName(cons));

   /* look for a slack variable s to convert a*x + s == b into lhs <= a*x <= rhs */
   vars = consdata->vars;
   vals = consdata->vals;
   bestslackpos = -1;
   bestslacktype = SCIP_VARTYPE_BINARY;
   bestslackdomrng = 0.0;
   coefszeroone = TRUE;
   coefsintegral = TRUE;
   varsintegral = TRUE;
   ncontvars = 0;
   contvarpos = -1;
   nintvars = 0;
   intvarpos = -1;
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Real absval;
      SCIP_Real varlb;
      SCIP_Real varub;
      SCIP_Bool iscont;
      int maxnlocks;

      assert(vars != NULL);
      assert(vals != NULL);

      var = vars[v];
      assert(!SCIPconsIsChecked(cons) || SCIPvarGetNLocksDown(var) >= 1); /* because variable is locked in this equality */
      assert(!SCIPconsIsChecked(cons) || SCIPvarGetNLocksUp(var) >= 1);
      varlb = SCIPvarGetLbGlobal(var);
      varub = SCIPvarGetUbGlobal(var);

      val = vals[v];
      absval = REALABS(val);
      assert(SCIPisPositive(scip, absval));

      slacktype = SCIPvarGetType(var);
      coefszeroone = coefszeroone && SCIPisEQ(scip, absval, 1.0);
      coefsintegral = coefsintegral && SCIPisIntegral(scip, val);
      varsintegral = varsintegral && (slacktype != SCIP_VARTYPE_CONTINUOUS);
      iscont = (slacktype == SCIP_VARTYPE_CONTINUOUS || slacktype == SCIP_VARTYPE_IMPLINT);

      /* update candidates for continuous -> implint and integer -> implint conversion */
      if( iscont )
      {
         ncontvars++;
         contvarpos = v;
      }
      else if( slacktype == SCIP_VARTYPE_INTEGER )
      {
         nintvars++;
         intvarpos = v;
      }

      /* check, if variable is already fixed or aggregated */
      if( !SCIPvarIsActive(var) )
         continue;

      /* check, if variable is used in other constraints */
      maxnlocks = SCIPconsIsChecked(cons) ? 1 : 0;
      if( SCIPvarGetNLocksDown(var) > maxnlocks || SCIPvarGetNLocksUp(var) > maxnlocks )
         continue;

      /* check, if variable can be used as a slack variable */
      if( iscont || (coefsintegral && varsintegral && SCIPisEQ(scip, absval, 1.0)) )
      {
         SCIP_Real slackdomrng;

         slackdomrng = (varub - varlb)*absval;
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
   if( coefsintegral && varsintegral && !SCIPisFeasIntegral(scip, consdata->rhs) )
   {
      SCIPdebugMessage("linear equality <%s> is integer infeasible\n", SCIPconsGetName(cons));
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if the slack variable is of integer type, and the constraint itself may take fractional values,
    * we cannot aggregate the variable, because the integrality condition would get lost
    */
   if( bestslackpos >= 0
      && (bestslacktype == SCIP_VARTYPE_CONTINUOUS || bestslacktype == SCIP_VARTYPE_IMPLINT
         || (coefsintegral && varsintegral)) )
   {
      SCIP_VAR* slackvar;
      SCIP_Real* scalars;
      SCIP_Real slackcoef;
      SCIP_Real slackvarlb;
      SCIP_Real slackvarub;
      SCIP_Real aggrconst;
      SCIP_Real newlhs;
      SCIP_Real newrhs;
      SCIP_Bool infeasible;
      SCIP_Bool aggregated;

      /* we found a slack variable that only occurs in at most one other constraint:
       *   a_1*x_1 + ... + a_k*x_k + a'*s == rhs  ->  s == rhs - a_1/a'*x_1 - ... - a_k/a'*x_k
       */
      assert(bestslackpos < consdata->nvars);

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
      SCIP_CALL( chgLhs(scip, cons, newlhs) );
      SCIP_CALL( chgRhs(scip, cons, newrhs) );
      SCIP_CALL( delCoefPos(scip, cons, bestslackpos) );

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &scalars, consdata->nvars) );

      /* set up the multi-aggregation */
      SCIPdebugMessage("linear constraint <%s>: multi-aggregate <%s> ==", SCIPconsGetName(cons), SCIPvarGetName(slackvar));
      for( v = 0; v < consdata->nvars; ++v )
      {
         scalars[v] = -consdata->vals[v]/slackcoef;
         SCIPdebugPrintf(" %+g<%s>", scalars[v], SCIPvarGetName(vars[v]));
      }
      SCIPdebugPrintf(" %+g, bounds of <%s>: [%g,%g]\n", aggrconst, SCIPvarGetName(slackvar), slackvarlb, slackvarub);

      /* perform the multi-aggregation */
      SCIP_CALL( SCIPmultiaggregateVar(scip, slackvar, consdata->nvars, vars, scalars, aggrconst,
            &infeasible, &aggregated) );
      assert(aggregated);

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &scalars);

      /* check for infeasible aggregation */
      if( infeasible )
      {
         SCIPdebugMessage("linear constraint <%s>: infeasible multi-aggregation\n", SCIPconsGetName(cons));
         *cutoff = TRUE;
         return SCIP_OKAY;
      }

      (*naggrvars)++;
   }
   else if( ncontvars == 1 )
   {
      SCIP_VAR* var;

      assert(0 <= contvarpos && contvarpos < consdata->nvars);
      var = vars[contvarpos];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);

      if( coefsintegral
         && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
         && SCIPisEQ(scip, REALABS(vals[contvarpos]), 1.0)
         && SCIPisFeasIntegral(scip, consdata->rhs) )
      {
         /* convert the continuous variable with coefficient 1.0 into an implicit integer variable */
         SCIPdebugMessage("linear constraint <%s>: converting continuous variable <%s> to implicit integer variable\n",
            SCIPconsGetName(cons), SCIPvarGetName(var));
         SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_IMPLINT) );
      }
   }
   else if( ncontvars == 0 && nintvars == 1 && !coefszeroone )
   {
      SCIP_VAR* var;

      /* this seems to help for rococo instances, but does not for rout (where all coefficients are +/- 1.0)
       *  -> we don't convert integers into implints if the row is a 0/1-row
       */
      assert(varsintegral);
      assert(0 <= intvarpos && intvarpos < consdata->nvars);
      var = vars[intvarpos];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

      if( coefsintegral
         && SCIPisEQ(scip, REALABS(vals[intvarpos]), 1.0)
         && SCIPisFeasIntegral(scip, consdata->rhs) )
      {
         /* convert the integer variable with coefficient 1.0 into an implicit integer variable */
         SCIPdebugMessage("linear constraint <%s>: converting integer variable <%s> to implicit integer variable\n",
            SCIPconsGetName(cons), SCIPvarGetName(var));
         SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_IMPLINT) );
      }
   }

   return SCIP_OKAY;
}

/* converts special equalities */
static
SCIP_RETCODE convertEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* do nothing on inequalities */
   if( !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      return SCIP_OKAY;

   /* depending on the number of variables, call a special conversion method */
   if( consdata->nvars == 1 )
   {
      /* fix variable */
      SCIP_CALL( convertUnaryEquality(scip, cons, cutoff, nfixedvars, ndelconss) );
   }
   else if( consdata->nvars == 2 )
   {
      /* aggregate one of the variables */
      SCIP_CALL( convertBinaryEquality(scip, cons, cutoff, naggrvars, ndelconss) );
   }
   else
   {
      /* try to multi-aggregate one of the variables */
      SCIP_CALL( convertLongEquality(scip, cons, cutoff, naggrvars) );
   }

   return SCIP_OKAY;
}

/** returns whether the linear sum of all variables/coefficients except the given one divided by the given value is always
 *  integral
 */
static
SCIP_Bool consdataIsResidualIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   int                   pos,                /**< position of variable to be left out */
   SCIP_Real             val                 /**< value to divide the coefficients by */
   )
{
   int v;
   
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   for( v = 0; v < consdata->nvars; ++v )
   {
      if( v != pos && (!SCIPvarIsIntegral(consdata->vars[v]) || !SCIPisIntegral(scip, consdata->vals[v]/val)) )
         return FALSE;
   }

   return TRUE;
}

/* applies dual presolving for variables that are locked only once in a direction, and this locking is due to a
 * linear inequality
 */
static
SCIP_RETCODE dualPresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool lhsexists;
   SCIP_Bool rhsexists;
   SCIP_Bool bestisint;
   SCIP_Bool bestislhs;
   int bestpos;
   int i;

   assert(cutoff != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);

   /* only process checked constraints (for which the locks are increased);
    * otherwise we would have to check for variables with nlocks == 0, and these are already processed by the
    * dualfix presolver
    */
   if( !SCIPconsIsChecked(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   lhsexists = !SCIPisInfinity(scip, -consdata->lhs);
   rhsexists = !SCIPisInfinity(scip, consdata->rhs);

   /* search for a single-locked variable which can be multi-aggregated; if a valid continuous variable was found, we
    * can use it safely for aggregation and break the search loop
    */
   bestpos = -1;
   bestisint = TRUE;
   bestislhs = FALSE;
   for( i = 0; i < consdata->nvars && bestisint; ++i )
   {
      SCIP_VAR* var;
      SCIP_Bool isint;
      SCIP_Real val;
      SCIP_Real obj;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool agglhs;
      SCIP_Bool aggrhs;

      var = consdata->vars[i];
      isint = (SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

      /* if we already found a candidate, skip integers */
      if( bestpos >= 0 && isint )
         continue;

      /* better do not multi-aggregate binary variables, since most plugins rely on their binary variables to be either
       * active, fixed, or single-aggregated with another binary variable
       */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY && consdata->nvars > 2 )
         continue;

      val = consdata->vals[i];
      obj = SCIPvarGetObj(var);
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      /* lhs <= a_0 * x_0 + a_1 * x_1 + ... + a_{n-1} * x_{n-1} <= rhs
       *
       * a_i >= 0, c_i >= 0, lhs exists, nlocksdown(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its lower bound
       *  - fix x_i to the smallest value for this constraint: x_i := lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *
       * a_i <= 0, c_i <= 0, lhs exists, nlocksup(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its upper bound
       *  - fix x_i to the largest value for this constraint: x_i := lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *
       * a_i >= 0, c_i <= 0, rhs exists, nlocksup(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its upper bound
       *  - fix x_i to the largest value for this constraint: x_i := rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *
       * a_i <= 0, c_i >= 0, rhs exists, nlocksdown(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its lower bound
       *  - fix x_i to the smallest value for this constraint: x_i := rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *
       * but: all this is only applicable, if the aggregated value is inside x_i's bounds for all possible values
       *      of all x_j
       */
      agglhs = lhsexists
         && ((val > 0.0 && SCIPvarGetNLocksDown(var) == 1 && !SCIPisNegative(scip, obj))
            || (val < 0.0 && SCIPvarGetNLocksUp(var) == 1 && !SCIPisPositive(scip, obj)));
      aggrhs = rhsexists
         && ((val > 0.0 && SCIPvarGetNLocksUp(var) == 1 && !SCIPisPositive(scip, obj))
            || (val < 0.0 && SCIPvarGetNLocksDown(var) == 1 && !SCIPisNegative(scip, obj)));
      if( agglhs || aggrhs )
      {
         SCIP_Real minresactivity;
         SCIP_Real maxresactivity;
         SCIP_Real minval;
         SCIP_Real maxval;

         /* calculate bounds for \sum_{j \neq i} a_j * x_j */
         consdataGetActivityResiduals(scip, consdata, var, val, &minresactivity, &maxresactivity);
         assert(minresactivity <= maxresactivity);

         if( agglhs )
         {
            /* check if lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i */
            if( val > 0.0 )
            {
               minval = (consdata->lhs - maxresactivity)/val;
               maxval = (consdata->lhs - minresactivity)/val;
            }
            else
            {
               minval = (consdata->lhs - minresactivity)/val;
               maxval = (consdata->lhs - maxresactivity)/val;
            }
            if( SCIPisFeasGE(scip, minval, lb) && SCIPisFeasLE(scip, maxval, ub) )
            {
               /* if the variable is integer, we have to check whether the integrality condition would always be satisfied
                * in the multi-aggregation
                */
               if( !isint || (SCIPisIntegral(scip, consdata->lhs/val) && consdataIsResidualIntegral(scip, consdata, i, val)) )
               {
                  bestpos = i;
                  bestisint = isint;
                  bestislhs = TRUE;
                  continue; /* no need to also look at the right hand side */
               }
            }
         }

         if( aggrhs )
         {
            /* check if rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i */
            if( val > 0.0 )
            {
               minval = (consdata->rhs - maxresactivity)/val;
               maxval = (consdata->rhs - minresactivity)/val;
            }
            else
            {
               minval = (consdata->rhs - minresactivity)/val;
               maxval = (consdata->rhs - maxresactivity)/val;
            }
            if( SCIPisFeasGE(scip, minval, lb) && SCIPisFeasLE(scip, maxval, ub) )
            {
               /* if the variable is integer, we have to check whether the integrality condition would always be satisfied
                * in the multi-aggregation
                */
               if( !isint || (SCIPisIntegral(scip, consdata->rhs/val) && consdataIsResidualIntegral(scip, consdata, i, val)) )
               {
                  bestpos = i;
                  bestisint = isint;
                  bestislhs = FALSE;
               }
            }
         }
      }
   }

   if( bestpos >= 0 )
   {
      SCIP_VAR** aggrvars;
      SCIP_Real* aggrcoefs;
      SCIP_Real aggrconst;
      SCIP_VAR* bestvar;
      SCIP_Real bestval;
      int naggrs;
      int j;
      SCIP_Bool infeasible;
      SCIP_Bool aggregated;

      assert(!bestislhs || lhsexists);
      assert(bestislhs || rhsexists);

      bestvar = consdata->vars[bestpos];
      bestval = consdata->vals[bestpos];
      assert(bestisint ==
         (SCIPvarGetType(bestvar) == SCIP_VARTYPE_BINARY || SCIPvarGetType(bestvar) == SCIP_VARTYPE_INTEGER));

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &aggrvars, consdata->nvars-1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &aggrcoefs, consdata->nvars-1) );
            
      /* set up the multi-aggregation */
      SCIPdebug(SCIPprintCons(scip, cons, NULL));
      SCIPdebugMessage("linear constraint <%s> (dual): multi-aggregate <%s> ==", SCIPconsGetName(cons), SCIPvarGetName(bestvar));
      naggrs = 0;
      for( j = 0; j < consdata->nvars; ++j )
      {
         if( j != bestpos )
         {
            aggrvars[naggrs] = consdata->vars[j];
            aggrcoefs[naggrs] = -consdata->vals[j]/consdata->vals[bestpos];
            SCIPdebugPrintf(" %+g<%s>", aggrcoefs[naggrs], SCIPvarGetName(aggrvars[naggrs]));
            assert(!bestisint || SCIPisIntegral(scip, aggrcoefs[naggrs]));
            naggrs++;
         }
      }
      aggrconst = (bestislhs ? consdata->lhs/bestval : consdata->rhs/bestval);
      SCIPdebugPrintf(" %+g, bounds of <%s>: [%g,%g]\n", aggrconst, SCIPvarGetName(bestvar),
         SCIPvarGetLbGlobal(bestvar), SCIPvarGetUbGlobal(bestvar));
      assert(!bestisint || SCIPisIntegral(scip, aggrconst));
      assert(naggrs == consdata->nvars-1);

      /* perform the multi-aggregation */
      SCIP_CALL( SCIPmultiaggregateVar(scip, bestvar, naggrs, aggrvars, aggrcoefs, aggrconst, &infeasible, &aggregated) );
       
      /* free temporary memory */
      SCIPfreeBufferArray(scip, &aggrcoefs);
      SCIPfreeBufferArray(scip, &aggrvars);
            
      /* check for infeasible aggregation */
      if( infeasible )
      {
         SCIPdebugMessage("linear constraint <%s>: infeasible multi-aggregation\n", SCIPconsGetName(cons));
         *cutoff = TRUE;
         return SCIP_OKAY;
      }

      /* delete the constraint, if the aggregation was successful */
      if( aggregated )
      {
         SCIP_CALL( SCIPdelCons(scip, cons) );
         
         if( !consdata->upgraded )
            (*ndelconss)++;
         (*naggrvars)++;
      }
   }

   return SCIP_OKAY;
}

/** converts all variables with fixed domain into FIXED variables */
static
SCIP_RETCODE fixVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars          /**< pointer to count the total number of fixed variables */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_VARSTATUS varstatus;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool fixed;
   SCIP_Bool infeasible;
   int v;

   assert(cutoff != NULL);
   assert(nfixedvars != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

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
            SCIP_Real fixval;

            fixval = SCIPselectSimpleValue(lb - SCIPepsilon(scip), ub + SCIPepsilon(scip), MAXDNOM);
            SCIPdebugMessage("converting variable <%s> with fixed bounds [%g,%g] into fixed variable fixed at %g\n",
               SCIPvarGetName(var), lb, ub, fixval);
            SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
            if( infeasible )
            {
               SCIPdebugMessage(" -> infeasible fixing\n");
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( fixed )
               (*nfixedvars)++;
         }
      }
   }

   SCIP_CALL( applyFixings(scip, cons) );
   assert(consdata->removedfixings);

   return SCIP_OKAY;
}

#define BINWEIGHT  1
#define INTWEIGHT  4
#define CONTWEIGHT 8

/** gets weight for variable in a "weighted number of variables" sum */
static
int getVarWeight(
   SCIP_VAR*             var                 /**< variable to get weight for */
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
      SCIPerrorMessage("invalid variable type\n");
      SCIPABORT();
      return 0; /*lint !e527*/
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
SCIP_RETCODE aggregateConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< (in)equality to modify */
   SCIP_CONS*            cons1,              /**< equality to use for aggregation of cons0 */
   int*                  commonidx0,         /**< array with indices of variables in cons0, that appear also in cons1 */
   int*                  commonidx1,         /**< array with indices of variables in cons1, that appear also in cons0 */
   int*                  diffidx0minus1,     /**< array with indices of variables in cons0, that don't appear in cons1 */
   int*                  diffidx1minus0,     /**< array with indices of variables in cons1, that don't appear in cons0 */
   int                   nvarscommon,        /**< number of variables, that appear in both constraints */
   int                   commonidxweight,    /**< variable weight sum of common variables */
   int                   diffidx0minus1weight, /**< variable weight sum of variables in cons0, that don't appear in cons1 */
   int                   diffidx1minus0weight, /**< variable weight sum of variables in cons1, that don't appear in cons0 */
   SCIP_Real             maxaggrnormscale,   /**< maximal allowed relative gain in maximum norm for constraint aggregation */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   SCIP_Bool*            aggregated          /**< pointer to store whether an aggregation was made */
   )
{
   SCIP_CONSDATA* consdata0;
   SCIP_CONSDATA* consdata1;
   SCIP_Real a;
   SCIP_Real b;
   SCIP_Real aggrcoef;
   SCIP_Real scalarsum;
   SCIP_Real bestscalarsum;
   SCIP_Bool betterscalarsum;
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
   assert(aggregated != NULL);

   assert(SCIPconsIsActive(cons0));
   assert(SCIPconsIsActive(cons1));

   SCIPdebugMessage("try aggregation of <%s> and <%s>\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));

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
         scalarsum = REALABS(a) + REALABS(b);
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
      SCIP_CONS* newcons;
      SCIP_CONSDATA* newconsdata;
      SCIP_VAR** newvars;
      SCIP_Real* newvals;
      SCIP_Real newlhs;
      SCIP_Real newrhs;
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

      SCIPdebugMessage("aggregate linear constraints <%s> := %g*<%s> + %g*<%s>  ->  nvars: %d -> %d, weight: %d -> %d\n",
         SCIPconsGetName(cons0), a, SCIPconsGetName(cons0), b, SCIPconsGetName(cons1),
         consdata0->nvars, bestnvars, commonidxweight + diffidx0minus1weight, bestvarweight);
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ));
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ));

      /* get temporary memory for creating the new linear constraint */
      SCIP_CALL( SCIPallocBufferArray(scip, &newvars, bestnvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newvals, bestnvars) );

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
      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, SCIPconsGetName(cons0), newnvars, newvars, newvals, newlhs, newrhs,
            SCIPconsIsInitial(cons0), SCIPconsIsSeparated(cons0), SCIPconsIsEnforced(cons0),
            SCIPconsIsChecked(cons0), SCIPconsIsPropagated(cons0),
            SCIPconsIsLocal(cons0), SCIPconsIsModifiable(cons0),
            SCIPconsIsDynamic(cons0), SCIPconsIsRemovable(cons0), SCIPconsIsStickingAtNode(cons0)) );

      newconsdata = SCIPconsGetData(newcons);
      assert(newconsdata != NULL);

      /* copy the upgraded flag from the old cons0 to the new constraint */
      newconsdata->upgraded = consdata0->upgraded;

      /* normalize the new constraint */
      SCIP_CALL( normalizeCons(scip, newcons) );

      /* check, if we really want to use the new constraint instead of the old one:
       * use the new one, if the maximum norm doesn't grow too much
       */
      if( consdataGetMaxAbsval(SCIPconsGetData(newcons)) <= maxaggrnormscale * consdataGetMaxAbsval(consdata0) )
      {
         SCIPdebugMessage(" -> aggregated to <%s>\n", SCIPconsGetName(newcons));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, newcons, NULL) ));

         /* update the statistics: we changed all coefficients */
         if( !consdata0->upgraded )
            (*nchgcoefs) += consdata0->nvars + consdata1->nvars - nvarscommon;
         *aggregated = TRUE;

         /* delete the old constraint, and add the new linear constraint to the problem */
         SCIP_CALL( SCIPdelCons(scip, cons0) );
         SCIP_CALL( SCIPaddCons(scip, newcons) );
      }

      /* release the new constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &newvals);
      SCIPfreeBufferArray(scip, &newvars);
   }

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
   SCIP_Real             maxaggrnormscale,   /**< maximal allowed relative gain in maximum norm for constraint aggregation */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgsides,          /**< pointer to count number of changed left/right hand sides */
   int*                  nchgcoefs           /**< pointer to count number of changed coefficients */
   )
{
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   int* commonidx0;
   int* commonidx1;
   int* diffidx0minus1;
   int* diffidx1minus0;
   SCIP_Longint possignature0;
   SCIP_Longint negsignature0;
   SCIP_Bool cons0changed;
   SCIP_Bool cons0isequality;
   int diffidx1minus0size;
   int c;

   assert(conss != NULL);
   assert(firstchange <= chkind);
   assert(cutoff != NULL);
   assert(ndelconss != NULL);
   assert(nchgsides != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);
   cons0isequality = SCIPisEQ(scip, consdata0->lhs, consdata0->rhs);

   /* sort the constraint */
   SCIP_CALL( consdataSort(scip, consdata0) );

   /* calculate bit signatures of cons0 for potentially positive and negative coefficients */
   consdataCalcSignatures(consdata0);
   possignature0 = consdata0->possignature;
   negsignature0 = consdata0->negsignature;

   /* get temporary memory for indices of common variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &commonidx0, consdata0->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &commonidx1, consdata0->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &diffidx0minus1, consdata0->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &diffidx1minus0, consdata0->nvars) );
   diffidx1minus0size = consdata0->nvars;

   /* check constraint against all prior constraints */
   cons0changed = consdata0->changed;
   consdata0->changed = FALSE;
   for( c = (cons0changed ? 0 : firstchange); c < chkind && !(*cutoff) && SCIPconsIsActive(cons0); ++c )
   {
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata1;
      SCIP_Longint possignature1;
      SCIP_Longint negsignature1;
      SCIP_Bool cons0dominateslhs;
      SCIP_Bool cons1dominateslhs;
      SCIP_Bool cons0dominatesrhs;
      SCIP_Bool cons1dominatesrhs;
      SCIP_Bool cons1isequality;
      SCIP_Bool coefsequal;
      SCIP_Bool coefsnegated;
      SCIP_Bool tryaggregation;
      int nvarscommon;
      int nvars0minus1;
      int nvars1minus0;
      int commonidxweight;
      int diffidx0minus1weight;
      int diffidx1minus0weight;
      int v0;
      int v1;

      cons1 = conss[c];

      /* ignore inactive and modifiable constraints */
      if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
         continue;

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

#if 0
      SCIPdebugMessage("preprocess linear constraint pair <%s>[chgd:%d, upgd:%d] and <%s>[chgd:%d, upgd:%d]\n",
         SCIPconsGetName(cons0), cons0changed, consdata0->upgraded,
         SCIPconsGetName(cons1), consdata1->changed, consdata1->upgraded);
#endif

      /* if both constraints didn't change since last pair processing, we can ignore the pair */
      if( !cons0changed && !consdata1->changed )
         continue;

      /* if both constraints are already upgraded, skip the pair */
      if( consdata0->upgraded && consdata1->upgraded )
         continue;

      assert(consdata1->nvars >= 1);

      /* sort the constraint */
      SCIP_CALL( consdataSort(scip, consdata1) );

      /* calculate bit signatures of cons1 for potentially positive and negative coefficients */
      consdataCalcSignatures(consdata1);
      possignature1 = consdata1->possignature;
      negsignature1 = consdata1->negsignature;

      /* the signatures give a quick test to check for domination and equality of coefficients */
      coefsequal = (possignature0 == possignature1) && (negsignature0 == negsignature1);
      coefsnegated = (possignature0 == negsignature1) && (negsignature0 == possignature1);
      cons0dominateslhs = SCIPisGE(scip, consdata0->lhs, consdata1->lhs)
         && ((possignature0 | possignature1) == possignature1)  /* possignature0 <= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature0); /* negsignature0 >= negsignature1 (as bit vector) */
      cons1dominateslhs = SCIPisGE(scip, consdata1->lhs, consdata0->lhs)
         && ((possignature0 | possignature1) == possignature0)  /* possignature0 >= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature1); /* negsignature0 <= negsignature1 (as bit vector) */
      cons0dominatesrhs = SCIPisLE(scip, consdata0->rhs, consdata1->rhs)
         && ((possignature0 | possignature1) == possignature0)  /* possignature0 >= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature1); /* negsignature0 <= negsignature1 (as bit vector) */
      cons1dominatesrhs = SCIPisLE(scip, consdata1->rhs, consdata0->rhs)
         && ((possignature0 | possignature1) == possignature1)  /* possignature0 <= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature0); /* negsignature0 >= negsignature1 (as bit vector) */
      cons1isequality = SCIPisEQ(scip, consdata1->lhs, consdata1->rhs);
      tryaggregation = (cons0isequality || cons1isequality) && (maxaggrnormscale > 0.0);
      if( !cons0dominateslhs && !cons1dominateslhs && !cons0dominatesrhs && !cons1dominatesrhs
         && !coefsequal && !coefsnegated && !tryaggregation )
         continue;

      /* make sure, we have enough memory for the index set of V_1 \ V_0 */
      if( tryaggregation && consdata1->nvars > diffidx1minus0size )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &diffidx1minus0, consdata1->nvars) );
         diffidx1minus0size = consdata1->nvars;
      }

      /* check consdata0 against consdata1:
       * - if lhs0 >= lhs1 and for each variable v and each solution value x_v val0[v]*x_v <= val1[v]*x_v,
       *   consdata0 dominates consdata1 w.r.t. left hand side
       * - if rhs0 <= rhs1 and for each variable v and each solution value x_v val0[v]*x_v >= val1[v]*x_v,
       *   consdata0 dominates consdata1 w.r.t. right hand side
       * - if val0[v] == -val1[v] for all variables v, the two inequalities can be replaced by a single
       *   ranged row (or equality)
       * - if at least one constraint is an equality, count the weighted number of common variables W_c
       *   and the weighted number of variable in the difference sets W_0 = w(V_0 \ V_1), W_1 = w(V_1 \ V_0),
       *   where the weight of each variable depends on its type, such that aggregations in order to remove the
       *   number of continuous and integer variables are preferred:
       *   - if W_c > W_1, try to aggregate  consdata0 := a * consdata0 + b * consdata1  in order to decrease the
       *     variable weight in consdata0, where a = +/- val1[v] and b = -/+ val0[v] for common v which leads to
       *     the smallest weight; for numerical stability, we will only accept integral a and b; the sign of a has
       *     to be positive to not switch the sense of the (in)equality cons0
       *   - if W_c > W_0, try to aggregate  consdata1 := a * consdata1 + b * consdata0  in order to decrease the
       *     variable weight in consdata1, where a = +/- val0[v] and b = -/+ val1[v] for common v which leads to
       *     the smallest weight; for numerical stability, we will only accept integral a and b; the sign of a has
       *     to be positive to not switch the sense of the (in)equality cons1
       */

      /* check consdata0 against consdata1 for redundancy, or ranged row accumulation */
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
            || coefsequal || coefsnegated || tryaggregation) )
      {
         SCIP_VAR* var;
         SCIP_Real val0;
         SCIP_Real val1;
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
            var = consdata0->vars[v0];
            val0 = consdata0->vals[v0];
            val1 = 0.0;
            if( tryaggregation )
            {
               diffidx0minus1[nvars0minus1] = v0;
               nvars0minus1++;
               diffidx0minus1weight += getVarWeight(var);
            }
            v0++;
            coefsequal = FALSE;
            coefsnegated = FALSE;
            break;

         case +1:
            /* variable doesn't appear in consdata0 */
            var = consdata1->vars[v1];
            val0 = 0.0;
            val1 = consdata1->vals[v1];
            if( tryaggregation )
            {
               diffidx1minus0[nvars1minus0] = v1;
               nvars1minus0++;
               diffidx1minus0weight += getVarWeight(var);
            }
            v1++;
            coefsequal = FALSE;
            coefsnegated = FALSE;
            break;

         case 0:
            /* variable appears in both constraints */
            assert(consdata0->vars[v0] == consdata1->vars[v1]);
            var = consdata0->vars[v0];
            val0 = consdata0->vals[v0];
            val1 = consdata1->vals[v1];
            if( tryaggregation )
            {
               commonidx0[nvarscommon] = v0;
               commonidx1[nvarscommon] = v1;
               nvarscommon++;
               commonidxweight += getVarWeight(var);
            }
            v0++;
            v1++;
            coefsequal = coefsequal && (SCIPisEQ(scip, val0, val1));
            coefsnegated = coefsnegated && (SCIPisEQ(scip, val0, -val1));
            break;

         default:
            SCIPerrorMessage("invalid comparison result\n");
            var = NULL;
            val0 = 0.0;
            val1 = 0.0;
            SCIPABORT();
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

      /* check for domination: remove dominated sides, but don't touch equalities as long as they are not totally
       * redundant
       */
      if( cons1dominateslhs && (!cons0isequality || cons1dominatesrhs || SCIPisInfinity(scip, consdata0->rhs) ) )
      {
         /* left hand side is dominated by consdata1: delete left hand side of consdata0 */
         SCIPdebugMessage("left hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ));

         /* check for infeasibility */
         if( SCIPisFeasGT(scip, consdata1->lhs, consdata0->rhs) )
         {
            SCIPdebugMessage("linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant left hand side */
         if( !SCIPisInfinity(scip, -consdata0->lhs) )
         {
            SCIP_CALL( chgLhs(scip, cons0, -SCIPinfinity(scip)) );
            cons0isequality = FALSE;
            if( !consdata0->upgraded )
               (*nchgsides)++;
         }
      }
      else if( cons0dominateslhs && (!cons1isequality || cons0dominatesrhs || SCIPisInfinity(scip, consdata1->rhs)) )
      {
         /* left hand side is dominated by consdata0: delete left hand side of consdata1 */
         SCIPdebugMessage("left hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons1), SCIPconsGetName(cons0));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ));

         /* check for infeasibility */
         if( SCIPisFeasGT(scip, consdata0->lhs, consdata1->rhs) )
         {
            SCIPdebugMessage("linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant left hand side */
         if( !SCIPisInfinity(scip, -consdata1->lhs) )
         {
            SCIP_CALL( chgLhs(scip, cons1, -SCIPinfinity(scip)) );
            cons1isequality = FALSE;
            if( !consdata1->upgraded )
               (*nchgsides)++;
         }
      }
      if( cons1dominatesrhs && (!cons0isequality || cons1dominateslhs || SCIPisInfinity(scip, -consdata0->lhs)) )
      {
         /* right hand side is dominated by consdata1: delete right hand side of consdata0 */
         SCIPdebugMessage("right hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ));

         /* check for infeasibility */
         if( SCIPisFeasLT(scip, consdata1->rhs, consdata0->lhs) )
         {
            SCIPdebugMessage("linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant right hand side */
         if( !SCIPisInfinity(scip, consdata0->rhs) )
         {
            SCIP_CALL( chgRhs(scip, cons0, SCIPinfinity(scip)) );
            cons0isequality = FALSE;
            if( !consdata0->upgraded )
               (*nchgsides)++;
         }
      }
      else if( cons0dominatesrhs && (!cons1isequality || cons0dominateslhs || SCIPisInfinity(scip, -consdata1->lhs)) )
      {
         /* right hand side is dominated by consdata0: delete right hand side of consdata1 */
         SCIPdebugMessage("right hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons1), SCIPconsGetName(cons0));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ));

         /* check for infeasibility */
         if( SCIPisFeasLT(scip, consdata0->rhs, consdata1->lhs) )
         {
            SCIPdebugMessage("linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant right hand side */
         if( !SCIPisInfinity(scip, consdata1->rhs) )
         {
            SCIP_CALL( chgRhs(scip, cons1, SCIPinfinity(scip)) );
            cons1isequality = FALSE;
            if( !consdata1->upgraded )
               (*nchgsides)++;
         }
      }

      /* check for now redundant constraints */
      if( SCIPisInfinity(scip, -consdata0->lhs) && SCIPisInfinity(scip, consdata0->rhs) )
      {
         /* consdata0 became redundant */
         SCIPdebugMessage("linear constraint <%s> is redundant\n", SCIPconsGetName(cons0));
         SCIP_CALL( SCIPdelCons(scip, cons0) );
         if( !consdata0->upgraded )
            (*ndelconss)++;
         continue;
      }
      if( SCIPisInfinity(scip, -consdata1->lhs) && SCIPisInfinity(scip, consdata1->rhs) )
      {
         /* consdata1 became redundant */
         SCIPdebugMessage("linear constraint <%s> is redundant\n", SCIPconsGetName(cons1));
         SCIP_CALL( SCIPdelCons(scip, cons1) );
         if( !consdata1->upgraded )
            (*ndelconss)++;
         continue;
      }

      /* check for disaggregated ranged rows */
      if( coefsequal || coefsnegated )
      {
         SCIP_Real lhs;
         SCIP_Real rhs;

         /* the coefficients in both rows are either equal or negated: create a new constraint with same coeffients and
          * best left and right hand sides; delete the old constraints afterwards
          */
         SCIPdebugMessage("aggregate linear constraints <%s> and <%s> with %s coefficients into single ranged row\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1), coefsequal ? "equal" : "negated");
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ));
         SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ));

         if( coefsequal )
         {
            /* the coefficients of both rows are equal */
            lhs = MAX(consdata0->lhs, consdata1->lhs);
            rhs = MIN(consdata0->rhs, consdata1->rhs);
         }
         else
         {
            /* the coefficients of both rows are negations */
            lhs = MAX(consdata0->lhs, -consdata1->rhs);
            rhs = MIN(consdata0->rhs, -consdata1->lhs);
         }
         if( SCIPisFeasLT(scip, rhs, lhs) )
         {
            SCIPdebugMessage("aggregated linear constraint <%s> is infeasible\n", SCIPconsGetName(cons0));
            *cutoff = TRUE;
            break;
         }

         /* update the sides of cons0 */
         SCIP_CALL( chgLhs(scip, cons0, lhs) );
         SCIP_CALL( chgRhs(scip, cons0, rhs) );

         /* update the flags of cons0 */
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

         /* delete cons1 */
         SCIP_CALL( SCIPdelCons(scip, cons1) );
         (*ndelconss)++;
         continue;
      }

      /* check, if we want to aggregate an (in)equality with an equality:
       *   consdata0 := a * consdata0 + b * consdata1  or  consdata1 := a * consdata1 + b * consdata0
       */
      if( tryaggregation )
      {
         SCIP_Bool aggregated;

         assert(consdata0->nvars == nvarscommon + nvars0minus1);
         assert(consdata1->nvars == nvarscommon + nvars1minus0);

         aggregated = FALSE;
         if( cons1isequality && !consdata0->upgraded && commonidxweight > diffidx1minus0weight )
         {
            /* W_c > W_1: try to aggregate  consdata0 := a * consdata0 + b * consdata1 */
            SCIP_CALL( aggregateConstraints(scip, cons0, cons1, commonidx0, commonidx1, diffidx0minus1, diffidx1minus0,
                  nvarscommon, commonidxweight, diffidx0minus1weight, diffidx1minus0weight, maxaggrnormscale,
                  nchgcoefs, &aggregated) );
         }
         if( !aggregated && cons0isequality && !consdata1->upgraded && commonidxweight > diffidx0minus1weight )
         {
            /* W_c > W_0: try to aggregate  consdata1 := a * consdata1 + b * consdata0 */
            SCIP_CALL( aggregateConstraints(scip, cons1, cons0, commonidx1, commonidx0, diffidx1minus0, diffidx0minus1,
                  nvarscommon, commonidxweight, diffidx1minus0weight, diffidx0minus1weight, maxaggrnormscale,
                  nchgcoefs, &aggregated) );
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &diffidx1minus0);
   SCIPfreeBufferArray(scip, &diffidx0minus1);
   SCIPfreeBufferArray(scip, &commonidx1);
   SCIPfreeBufferArray(scip, &commonidx0);

   return SCIP_OKAY;
}

/** applies full dual presolving on variables that only appear in linear constraints */
static
SCIP_RETCODE fullDualPresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints */
   int*                  nchgbds             /**< pointer to count the number of bound changes */
   )
{
   SCIP_Real* redlb;
   SCIP_Real* redub;
   int* nlocksdown;
   int* nlocksup;
   SCIP_VAR** vars;
   int nvars;
   int nbinvars;
   int v;
   int c;

   /* we calculate redundancy bounds with the following meaning:
    *   redlb[v] == k : if x_v >= k, we can always round x_v down to x_v == k without violating any constraint
    *   redub[v] == k : if x_v <= k, we can always round x_v up to x_v == k without violating any constraint
    * then:
    *   c_v >= 0 : x_v <= redlb[v] is feasible due to optimality
    *   c_v <= 0 : x_v >= redub[v] is feasible due to optimality
    */

   assert(nconss == 0 || conss != NULL);
   assert(nchgbds != NULL);

   /* get active variables */
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   /* if the problem is a pure binary program, nothing can be achieved by full dual presolve */
   nbinvars = SCIPgetNBinVars(scip);
   if( nbinvars == nvars )
      return SCIP_OKAY;

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &redlb, nvars - nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redub, nvars - nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksdown, nvars - nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksup, nvars - nbinvars) );

   /* initialize redundancy bounds */
   for( v = nbinvars; v < nvars; ++v )
   {
      redlb[v - nbinvars] = SCIPvarGetLbGlobal(vars[v]);
      redub[v - nbinvars] = SCIPvarGetUbGlobal(vars[v]);
   }
   BMSclearMemoryArray(nlocksdown, nvars - nbinvars);
   BMSclearMemoryArray(nlocksup, nvars - nbinvars);

   /* scan all constraints */
   for( c = 0; c < nconss; ++c )
   {
      /* we only need to consider constraints that have been locked (i.e., checked constraints or constraints that are
       * part of checked disjunctions)
       */
      if( SCIPconsIsLocked(conss[c]) )
      {
         SCIP_CONSDATA* consdata;
         SCIP_Bool lhsexists;
         SCIP_Bool rhsexists;
         int nlockspos;
         int i;

         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);

         /* get number of times the constraint was locked */
         nlockspos = SCIPconsGetNLocksPos(conss[c]);

         /* we do not want to include constraints with locked negation (this would be to weird) */
         if( SCIPconsGetNLocksNeg(conss[c]) > 0 )
            continue;

         /* check for existing sides */
         lhsexists = !SCIPisInfinity(scip, -consdata->lhs);
         rhsexists = !SCIPisInfinity(scip, consdata->rhs);

         /* count locks and update redundancy bounds */
         for( i = 0; i < consdata->nvars; ++i )
         {
            SCIP_VAR* var;
            SCIP_Real val;
            SCIP_Real minresactivity;
            SCIP_Real maxresactivity;
            SCIP_Real newredlb;
            SCIP_Real newredub;
            int arrayindex;

            /* we do not need to process binary variables */
            var = consdata->vars[i];
            if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
               continue;

            /* calculate residual activity bounds if variable would be fixed to zero */
            val = consdata->vals[i];
            consdataGetGlbActivityResiduals(scip, consdata, var, val, &minresactivity, &maxresactivity);

            arrayindex = SCIPvarGetProbindex(consdata->vars[i]) - nbinvars;
            assert(0 <= arrayindex && arrayindex < nvars - nbinvars); /* variable should be active due to applyFixings() */

            newredlb = redlb[arrayindex];
            newredub = redub[arrayindex];
            if( val > 0.0 )
            {
               if( lhsexists )
               {
                  /* lhs <= d*x + a*y, d > 0  ->  redundant in y if  x >= (lhs - min{a*y})/d */
                  nlocksdown[arrayindex] += nlockspos;
                  newredlb = (SCIPisInfinity(scip, -minresactivity) ? SCIPinfinity(scip)
                     : (consdata->lhs - minresactivity)/val);
               }
               if( rhsexists )
               {
                  /* d*x + a*y <= rhs, d > 0  ->  redundant in y if  x <= (rhs - max{a*y})/d */
                  nlocksup[arrayindex] += nlockspos;
                  newredub = (SCIPisInfinity(scip, maxresactivity) ? -SCIPinfinity(scip)
                     : (consdata->rhs - maxresactivity)/val);
               }
            }
            else
            {
               if( lhsexists )
               {
                  /* lhs <= d*x + a*y, d < 0  ->  redundant in y if  x <= (lhs - min{a*y})/d */
                  nlocksup[arrayindex] += nlockspos;
                  newredub = (SCIPisInfinity(scip, -minresactivity) ? -SCIPinfinity(scip)
                     : (consdata->lhs - minresactivity)/val);
               }
               if( rhsexists )
               {
                  /* d*x + a*y <= rhs, d < 0  ->  redundant in y if  x >= (rhs - max{a*y})/d */
                  nlocksdown[arrayindex] += nlockspos;
                  newredlb = (SCIPisInfinity(scip, maxresactivity) ? SCIPinfinity(scip)
                     : (consdata->rhs - maxresactivity)/val);
               }
            }

            /* if the variable is integer, we have to round the value to the next integral value */
            if( SCIPvarIsIntegral(var) )
            {
               if( !SCIPisInfinity(scip, newredlb) )
                  newredlb = SCIPceil(scip, newredlb);
               if( !SCIPisInfinity(scip, -newredub) )
                  newredub = SCIPfloor(scip, newredub);
            }

            /* update redundancy bounds */
            redlb[arrayindex] = MAX(redlb[arrayindex], newredlb);
            redub[arrayindex] = MIN(redub[arrayindex], newredub);
         }
      }
   }

   /* check if any bounds can be tightened due to optimality */
   for( v = nbinvars; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real obj;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      assert(SCIPvarGetProbindex(vars[v]) == v);
      assert(SCIPvarGetType(vars[v]) != SCIP_VARTYPE_BINARY);
      assert(SCIPvarGetNLocksDown(vars[v]) >= nlocksdown[v - nbinvars]);
      assert(SCIPvarGetNLocksUp(vars[v]) >= nlocksup[v - nbinvars]);

      var = vars[v];
      obj = SCIPvarGetObj(var);
      if( obj >= 0.0 )
      {
         /* making the variable as small as possible does not increase the objective:
          * check if all down locks of the variables are due to linear constraints
          */
         if( SCIPvarGetNLocksDown(var) == nlocksdown[v - nbinvars] && redlb[v - nbinvars] < SCIPvarGetUbGlobal(var) )
         {
            /* if x_v >= redlb[v], we can always round x_v down to x_v == redlb[v] without violating any constraint
             *  -> tighten upper bound to x_v <= redlb[v]
             */
            SCIPdebugMessage("variable <%s> only locked down in linear constraints: dual presolve <%s>[%g,%g] <= %g\n",
               SCIPvarGetName(var), SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
               redlb[v - nbinvars]);
            SCIP_CALL( SCIPtightenVarUb(scip, var, redlb[v - nbinvars], &infeasible, &tightened) );
            assert(!infeasible);
            if( tightened )
               (*nchgbds)++;
         }
      }
      if( obj <= 0.0 )
      {
         /* making the variable as large as possible does not increase the objective:
          * check if all up locks of the variables are due to linear constraints
          */
         if( SCIPvarGetNLocksUp(var) == nlocksup[v - nbinvars] && redub[v - nbinvars] > SCIPvarGetLbGlobal(var) )
         {
            /* if x_v <= redub[v], we can always round x_v up to x_v == redub[v] without violating any constraint
             *  -> tighten lower bound to x_v >= redub[v]
             */
            SCIPdebugMessage("variable <%s> only locked up in linear constraints: dual presolve <%s>[%g,%g] >= %g\n",
               SCIPvarGetName(var), SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
               redub[v - nbinvars]);
            SCIP_CALL( SCIPtightenVarLb(scip, var, redub[v - nbinvars], &infeasible, &tightened) );
            assert(!infeasible);
            if( tightened )
               (*nchgbds)++;
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &nlocksup);
   SCIPfreeBufferArray(scip, &nlocksdown);
   SCIPfreeBufferArray(scip, &redub);
   SCIPfreeBufferArray(scip, &redlb);

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

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
static
SCIP_DECL_CONSINIT(consInitLinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* catch events for the constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->eventdatas == NULL);

      /* catch all events */
      SCIP_CALL( consdataCatchAllEvents(scip, consdata, conshdlrdata->eventhdlr) );
      assert(consdata->eventdatas != NULL);
   }

   return SCIP_OKAY;

}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitLinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* drop events for the constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->eventdatas != NULL )
      {
         /* drop all events */
         SCIP_CALL( consdataDropAllEvents(scip, consdata, conshdlrdata->eventhdlr) );
         assert(consdata->eventdatas == NULL);
      }
   }

   return SCIP_OKAY;

}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreLinear NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreLinear)
{  /*lint --e{715}*/
   int c;

   /* delete all linear constraints that were upgraded to a more specific constraint type;
    * make sure, only active variables remain in the remaining constraints
    */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->upgraded )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
      }
      else
      {
         SCIP_CALL( applyFixings(scip, conss[c]) );
      }
   }

   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolLinear NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolLinear)
{  /*lint --e{715}*/
   int c;

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   /* if this is a restart, convert cutpool rows into linear constraints */
   if( restart )
   {
      SCIP_CUT** cuts;
      int ncuts;
      int ncutsadded;

      ncutsadded = 0;
      cuts = SCIPgetPoolCuts(scip);
      ncuts = SCIPgetNPoolCuts(scip);
      for( c = 0; c < ncuts; ++c )
      {
         SCIP_ROW* row;

         row = SCIPcutGetRow(cuts[c]);
         assert(!SCIProwIsLocal(row));
         assert(!SCIProwIsModifiable(row));
         if( SCIPcutGetAge(cuts[c]) == 0 && SCIProwIsInLP(row) )
         {
            char name[SCIP_MAXSTRLEN];
            SCIP_CONS* cons;
            SCIP_COL** cols;
            SCIP_VAR** vars;
            int ncols;
            int i;

            /* create a linear constraint out of the cut */
            cols = SCIProwGetCols(row);
            ncols = SCIProwGetNNonz(row);
            
            SCIP_CALL( SCIPallocBufferArray(scip, &vars, ncols) );
            for( i = 0; i < ncols; ++i )
               vars[i] = SCIPcolGetVar(cols[i]);
            
            sprintf(name, "%s_%d", SCIProwGetName(row), SCIPgetNRuns(scip));
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ncols, vars, SCIProwGetVals(row),
                  SCIProwGetLhs(row) - SCIProwGetConstant(row), SCIProwGetRhs(row) - SCIProwGetConstant(row),
                  TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            SCIPfreeBufferArray(scip, &vars);

            ncutsadded++;
         }
      }

      if( ncutsadded > 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "(restart) converted %d cuts from the global cut pool into linear constraints\n\n", ncutsadded);
      }
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* free linear constraint */
   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   /*debugMessage("Trans method of linear constraints\n");*/

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
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
   SCIP_CALL( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr,
         sourcedata->nvars, sourcedata->vars, sourcedata->vals, sourcedata->lhs, sourcedata->rhs) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpLinear)
{  /*lint --e{715}*/
   int c;

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   for( c = 0; c < nconss; ++c )
   {
      assert(SCIPconsIsInitial(conss[c]));
      SCIP_CALL( addRelaxation(scip, conss[c], NULL) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;
   int depth;
   int nrounds;
   int maxsepacuts;
   int ncuts;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   /*debugMessage("Sepa method of linear constraints\n");*/

   *result = SCIP_DIDNOTRUN;

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && conshdlrdata->maxroundsroot >= 0 && nrounds >= conshdlrdata->maxroundsroot)
      || (depth > 0 && conshdlrdata->maxrounds >= 0 && nrounds >= conshdlrdata->maxrounds) )
      return SCIP_OKAY;

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss && ncuts < maxsepacuts; ++c )
   {
      /*debugMessage("separating linear constraint <%s>\n", SCIPconsGetName(conss[c]));*/
      SCIP_CALL( separateCons(scip, conss[c], NULL, conshdlrdata->separateall, &ncuts) );
   }

   /* adjust return value */
   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* combine linear constraints to get more cuts */
   /**@todo further cuts of linear constraints */

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;
   int depth;
   int nrounds;
   int maxsepacuts;
   int ncuts;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   /*debugMessage("Sepa method of linear constraints\n");*/

   *result = SCIP_DIDNOTRUN;

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && conshdlrdata->maxroundsroot >= 0 && nrounds >= conshdlrdata->maxroundsroot)
      || (depth > 0 && conshdlrdata->maxrounds >= 0 && nrounds >= conshdlrdata->maxrounds) )
      return SCIP_OKAY;

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss && ncuts < maxsepacuts; ++c )
   {
      /*debugMessage("separating linear constraint <%s>\n", SCIPconsGetName(conss[c]));*/
      SCIP_CALL( separateCons(scip, conss[c], sol, conshdlrdata->separateall, &ncuts) );
   }

   /* adjust return value */
   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* combine linear constraints to get more cuts */
   /**@todo further cuts of linear constraints */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /*debugMessage("Enfolp method of linear constraints\n");*/

   /* check for violated constraints
    * LP is processed at current node -> we can add violated linear constraints to the SCIP_LP
    */
   *result = SCIP_FEASIBLE;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], NULL, FALSE, &violated) );

      if( violated )
      {
         /* insert LP row as cut */
         SCIP_CALL( addRelaxation(scip, conss[c], NULL) );
         *result = SCIP_SEPARATED;
      }
   }

   /* check all obsolete linear constraints for feasibility */
   for( c = nusefulconss; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], NULL, FALSE, &violated) );

      if( violated )
      {
         /* insert LP row as cut */
         SCIP_CALL( addRelaxation(scip, conss[c], NULL) );
         *result = SCIP_SEPARATED;
      }
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsLinear)
{  /*lint --e{715}*/
   SCIP_Bool violated;
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
      SCIP_CALL( checkCons(scip, conss[c], NULL, TRUE, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckLinear)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Check method of linear constraints\n");*/

   /* check all linear constraints for feasibility */
   violated = FALSE;
   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], sol, checklprows, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropLinear)
{  /*lint --e{715}*/
   SCIP_Bool tightenbounds;
   SCIP_Bool cutoff;
   int nchgbds;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Prop method of linear constraints\n");*/

   /* check, if we want to tighten variable's bounds (in probing, we always want to tighten the bounds) */
   if( SCIPinProbing(scip) )
      tightenbounds = TRUE;
   else
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int depth;
      int propfreq;
      int tightenboundsfreq;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      depth = SCIPgetDepth(scip);
      propfreq = SCIPconshdlrGetPropFreq(conshdlr);
      tightenboundsfreq = propfreq * conshdlrdata->tightenboundsfreq;
      tightenbounds = (conshdlrdata->tightenboundsfreq >= 0)
         && ((tightenboundsfreq == 0 && depth == 0) || (tightenboundsfreq >= 1 && (depth % tightenboundsfreq == 0)));
   }

   cutoff = FALSE;
   nchgbds = 0;

   /* process useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( propagateCons(scip, conss[c], tightenbounds, &cutoff, &nchgbds) );
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


#define MAXCONSPRESOLROUNDS 10
/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Bool cutoff;
   SCIP_Bool delay;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnupgdconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int firstchange;
   int firstupgradetry;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMessage("Presol method of linear constraints\n");*/

   /* remember old preprocessing counters */
   cutoff = FALSE;
   delay = FALSE;
   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* process single constraints */
   firstchange = INT_MAX;
   firstupgradetry = INT_MAX;
   for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
   {
      int npresolrounds;

      cons = conss[c];
      assert(SCIPconsIsActive(cons));
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* constraint should not be already presolved in the initial round */
      assert(SCIPgetNRuns(scip) > 0 || nrounds > 0 || !consdata->propagated);
      assert(SCIPgetNRuns(scip) > 0 || nrounds > 0 || !consdata->boundstightened);
      assert(SCIPgetNRuns(scip) > 0 || nrounds > 0 || !consdata->presolved);
      assert(consdata->propagated || !consdata->presolved);

      /* incorporate fixings and aggregations in constraint */
      SCIP_CALL( applyFixings(scip, cons) );
      assert(consdata->removedfixings);

      /* we can only presolve linear constraints, that are not modifiable */
      if( SCIPconsIsModifiable(cons) )
         continue;

      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == INT_MAX && consdata->changed )
         firstchange = c;

      /* remember the first constraint that was not yet tried to be upgraded, to begin the next upgrading round with */
      if( firstupgradetry == INT_MAX && !consdata->upgradetried )
         firstupgradetry = c;

      /* check, if constraint is already preprocessed */
      if( consdata->presolved )
         continue;

      assert(SCIPconsIsActive(cons));

      SCIPdebugMessage("presolving linear constraint <%s>\n", SCIPconsGetName(cons));
      SCIPdebug(SCIP_CALL( SCIPprintCons(scip, cons, NULL) ));

      /* apply presolving as long as possible on the single constraint (however, abort after a certain number of rounds
       * to avoid nearly infinite cycling due to very small bound changes)
       */
      npresolrounds = 0;
      while( !consdata->presolved && npresolrounds < MAXCONSPRESOLROUNDS )
      {
         assert(!cutoff);

         npresolrounds++;

         /* mark constraint being presolved and propagated */
         consdata->presolved = TRUE;
         consdata->propagated = TRUE;

         /* normalize constraint */
         SCIP_CALL( normalizeCons(scip, cons) );

         /* tighten left and right hand side due to integrality */
         SCIP_CALL( tightenSides(scip, cons, nchgsides) );

         /* check bounds */
         if( SCIPisFeasGT(scip, consdata->lhs, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is infeasible: sides=[%g,%g]\n",
               SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
            cutoff = TRUE;
            break;
         }

         /* tighten variable's bounds */
         SCIP_CALL( tightenBounds(scip, cons, &cutoff, nchgbds) );
         if( cutoff )
            break;

         /* check for fixed variables */
         SCIP_CALL( fixVariables(scip, cons, &cutoff, nfixedvars) );
         if( cutoff )
            break;

         /* check constraint for infeasibility and redundancy */
         consdataGetActivityBounds(scip, consdata, &minactivity, &maxactivity);
         if( SCIPisFeasGT(scip, minactivity, consdata->rhs) || SCIPisFeasLT(scip, maxactivity, consdata->lhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is infeasible: activitybounds=[%g,%g], sides=[%g,%g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            cutoff = TRUE;
            break;
         }
         else if( SCIPisFeasGE(scip, minactivity, consdata->lhs) && SCIPisFeasLE(scip, maxactivity, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s> is redundant: activitybounds=[%g,%g], sides=[%g,%g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( SCIPdelCons(scip, cons) );
            assert(!SCIPconsIsActive(cons));
            if( !consdata->upgraded )
               (*ndelconss)++;
            break;
         }
         else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasGE(scip, minactivity, consdata->lhs) )
         {
            SCIPdebugMessage("linear constraint <%s> left hand side is redundant: activitybounds=[%g,%g], sides=[%g,%g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( chgLhs(scip, cons, -SCIPinfinity(scip)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         else if( !SCIPisInfinity(scip, consdata->rhs) && SCIPisFeasLE(scip, maxactivity, consdata->rhs) )
         {
            SCIPdebugMessage("linear constraint <%s> right hand side is redundant: activitybounds=[%g,%g], sides=[%g,%g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( chgRhs(scip, cons, SCIPinfinity(scip)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         assert(consdata->nvars >= 1); /* otherwise, it should be redundant or infeasible */

         /* extract cliques from constraint */
         SCIP_CALL( extractCliques(scip, cons, &cutoff, nchgbds) );
         if( cutoff )
            break;

         /* reduce big-M coefficients, that make the constraint redundant if the variable is on a bound */
         SCIP_CALL( consdataTightenCoefs(scip, cons, nchgcoefs, nchgsides) );
      }

      /* convert special equalities */
      if( !cutoff && SCIPconsIsActive(cons) )
      {
         assert(consdata->propagated);
         assert(consdata->boundstightened);
         assert(consdata->presolved);

         SCIP_CALL( convertEquality(scip, cons, &cutoff, nfixedvars, naggrvars, ndelconss) );
      }

      /* apply dual presolving for variables that appear in only one constraint */
      if( !cutoff && SCIPconsIsActive(cons) )
      {
         SCIP_CALL( dualPresolve(scip, cons, &cutoff, nfixedvars, naggrvars, ndelconss) );
      }


      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == INT_MAX && consdata->changed )
         firstchange = c;

      /* remember the first constraint that was not yet tried to be upgraded, to begin the next upgrading round with */
      if( firstupgradetry == INT_MAX && !consdata->upgradetried )
         firstupgradetry = c;
   }

   /* process pairs of constraints: check them for redundancy and try to aggregate them;
    * only apply this expensive procedure, if the single constraint preprocessing did not find any reductions
    * (otherwise, we delay the presolving to be called again next time)
    */
   if( !cutoff
      && *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars && *nchgbds == oldnchgbds && *ndelconss == oldndelconss
      && *nupgdconss == oldnupgdconss && *nchgcoefs == oldnchgcoefs && *nchgsides == oldnchgsides )
   {
      if( conshdlrdata->maxpresolpairrounds == -1 || nrounds < conshdlrdata->maxpresolpairrounds )
      {
         for( c = firstchange; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
         {
            if( SCIPconsIsActive(conss[c]) && !SCIPconsIsModifiable(conss[c]) )
            {
               SCIP_CALL( preprocessConstraintPairs(scip, conss, firstchange, c, conshdlrdata->maxaggrnormscale,
                     &cutoff, ndelconss, nchgsides, nchgcoefs) );
            }
         }
      }
   }
   else
      delay = TRUE;

   /* before upgrading, check whether we can apply some additional dual presolving, because a variable only appears
    * in linear constraints and we therefore have full information about it
    */
   if( !cutoff && firstupgradetry < nconss
      && *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars && *nchgbds == oldnchgbds && *ndelconss == oldndelconss
      && *nupgdconss == oldnupgdconss && *nchgcoefs == oldnchgcoefs && *nchgsides == oldnchgsides
       )
   {
      SCIP_CALL( fullDualPresolve(scip, conss, nconss, nchgbds) );
   }

   /* try to upgrade constraints into a more specific constraint type;
    * only upgrade constraints, if no reductions were found in this round (otherwise, the linear constraint handler
    * may find additional reductions before giving control away to other (less intelligent?) constraint handlers)
    */
   if( !cutoff
      && *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars && *nchgbds == oldnchgbds && *ndelconss == oldndelconss
      && *nupgdconss == oldnupgdconss && *nchgcoefs == oldnchgcoefs && *nchgsides == oldnchgsides
       )
   {
      for( c = firstupgradetry; c < nconss; ++c )
      {
         cons = conss[c];

         /* don't upgrade modifiable constraints */
         if( SCIPconsIsModifiable(cons) )
            continue;

         consdata = SCIPconsGetData(cons);
         assert(consdata != NULL);

         /* only upgrade completely presolved constraints, that changed since the last upgrading call */
         if( consdata->upgradetried )
            continue;
         if( !consdata->presolved )
         {
            delay = TRUE;
            continue;
         }

         consdata->upgradetried = TRUE;
         if( SCIPconsIsActive(cons) )
         {
            SCIP_CONS* upgdcons;

            SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &upgdcons) );
            if( upgdcons != NULL )
            {
               /* add the upgraded constraint to the problem */
               SCIP_CALL( SCIPaddCons(scip, upgdcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &upgdcons) );
               (*nupgdconss)++;

               /* mark the linear constraint being upgraded and to be removed after presolving;
                * don't delete it directly, because it may help to preprocess other linear constraints
                */
               assert(!consdata->upgraded);
               consdata->upgraded = TRUE;

               /* delete upgraded inequalities immediately;
                * delete upgraded equalities, if we don't need it anymore for aggregation and redundancy checking
                */
               if( SCIPisLT(scip, consdata->lhs, consdata->rhs)
                  || (conshdlrdata->maxpresolpairrounds != -1 && nrounds > conshdlrdata->maxpresolpairrounds)
                  || (conshdlrdata->maxaggrnormscale == 0.0) )
               {
                  SCIP_CALL( SCIPdelCons(scip, cons) );
               }
            }
         }
      }
   }
   else
      delay = TRUE;

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( delay )
      *result = SCIP_DELAYED;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nchgbds > oldnchgbds || *ndelconss > oldndelconss
      || *nupgdconss > oldnupgdconss || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropLinear)
{  /*lint --e{715}*/
   SCIP_CALL( resolvePropagation(scip, cons, infervar, intToInferInfo(inferinfo), boundtype, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   haslhs = !SCIPisInfinity(scip, -consdata->lhs);
   hasrhs = !SCIPisInfinity(scip, consdata->rhs);

   /* update rounding locks of every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( SCIPisPositive(scip, consdata->vals[i]) )
      {
         if( haslhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos, nlocksneg) );
         }
         if( hasrhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( haslhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
         }
         if( hasrhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos, nlocksneg) );
         }
      }
   }

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

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintLinear)
{  /*lint --e{715}*/
   consdataPrint(scip, SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}




/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   consdata = eventdata->consdata;
   assert(consdata != NULL);

   eventtype = SCIPeventGetType(event);
   var = SCIPeventGetVar(event);

   if( (eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) != 0 )
   {
      SCIP_Real oldbound;
      SCIP_Real newbound;
      SCIP_Real val;
      int varpos;

      varpos = eventdata->varpos;
      assert(0 <= varpos && varpos < consdata->nvars);
      oldbound = SCIPeventGetOldbound(event);
      newbound = SCIPeventGetNewbound(event);
      assert(var != NULL);
      assert(consdata->vars[varpos] == var);
      val = consdata->vals[varpos];

      /* update the activity values */
      if( (eventtype & SCIP_EVENTTYPE_LBCHANGED) != 0 )
         consdataUpdateChgLb(scip, consdata, var, oldbound, newbound, val);
      else
      {
         assert((eventtype & SCIP_EVENTTYPE_UBCHANGED) != 0);
         consdataUpdateChgUb(scip, consdata, var, oldbound, newbound, val);
      }

      consdata->presolved = FALSE;

      /* bound change can turn the constraint infeasible or redundant only if it was a tightening */
      if( (eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED) != 0 )
         consdata->propagated = FALSE;

      /* check whether bound tightening might now be successful (if the current bound was relaxed, it might be
       * that it can be tightened again)
       */
      if( consdata->boundstightened )
      {
         switch( eventtype )
         {
         case SCIP_EVENTTYPE_LBTIGHTENED:
         case SCIP_EVENTTYPE_UBRELAXED:
            consdata->boundstightened = (val > 0.0 && SCIPisInfinity(scip, consdata->rhs))
               || (val < 0.0 && SCIPisInfinity(scip, -consdata->lhs));
            break;
         case SCIP_EVENTTYPE_LBRELAXED:
         case SCIP_EVENTTYPE_UBTIGHTENED:
            consdata->boundstightened = (val > 0.0 && SCIPisInfinity(scip, -consdata->lhs))
               || (val < 0.0 && SCIPisInfinity(scip, consdata->rhs));
            break;
         default:
            SCIPerrorMessage("invalid event type %d\n", eventtype);
            return SCIP_INVALIDDATA;
         }
      }
   }

   if( (eventtype & SCIP_EVENTTYPE_VARFIXED) != 0 )
   {
      /* we want to remove the fixed variable */
      consdata->presolved = FALSE;
      consdata->removedfixings = FALSE;
   }

   if( (eventtype & SCIP_EVENTTYPE_VARUNLOCKED) != 0 )
   {
      /* there is only one lock left: we may multiaggregate the variable as slack of an equation */
      assert(SCIPvarGetNLocksDown(var) <= 1);
      assert(SCIPvarGetNLocksUp(var) <= 1);
      consdata->presolved = FALSE;
   }

   if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
   {
      SCIP_Real oldbound;
      SCIP_Real newbound;
      SCIP_Real val;
      int varpos;

      varpos = eventdata->varpos;
      assert(0 <= varpos && varpos < consdata->nvars);
      oldbound = SCIPeventGetOldbound(event);
      newbound = SCIPeventGetNewbound(event);
      assert(var != NULL);
      assert(consdata->vars[varpos] == var);
      val = consdata->vals[varpos];

      /* update the activity values */
      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
         consdataUpdateChgGlbLb(scip, consdata, oldbound, newbound, val);
      else
      {
         assert((eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0);
         consdataUpdateChgGlbUb(scip, consdata, oldbound, newbound, val);
      }
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of conflict handler
 */

static
SCIP_DECL_CONFLICTEXEC(conflictExecLinear)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real lhs;
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

   /* create array of variables and coefficients: sum_{i \in P} x_i - sum_{i \in N} x_i >= 1 - |N| */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nbdchginfos) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nbdchginfos) );
   lhs = 1.0;
   for( i = 0; i < nbdchginfos; ++i )
   {
      assert(bdchginfos != NULL);

      vars[i] = SCIPbdchginfoGetVar(bdchginfos[i]);

      /* we can only treat binary variables */
      /**@todo extend linear conflict constraints to some non-binary cases */
      if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_BINARY )
         break;

      /* check whether the variable is fixed to zero (P) or one (N) in the conflict set */
      if( SCIPbdchginfoGetNewbound(bdchginfos[i]) < 0.5 )
         vals[i] = 1.0;
      else
      {
         vals[i] = -1.0;
         lhs -= 1.0;
      }
   }

   if( i == nbdchginfos )
   {
      SCIP_CONS* cons;
      SCIP_CONS* upgdcons;
      char consname[SCIP_MAXSTRLEN];

      /* create a constraint out of the conflict set */
      sprintf(consname, "cf%"SCIP_LONGINT_FORMAT, SCIPgetNConflictConssApplied(scip));
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, nbdchginfos, vars, vals, lhs, SCIPinfinity(scip),
            FALSE, TRUE, FALSE, FALSE, TRUE, local, FALSE, dynamic, removable, FALSE) );

      /** try to automatically convert a linear constraint into a more specific and more specialized constraint */
      SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &upgdcons) );
      if( upgdcons != NULL )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         cons = upgdcons;
      }
      
      /* add constraint to SCIP */
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, validnode) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      
      *result = SCIP_CONSADDED;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for linear constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrLinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecLinear,
         NULL) );

   /* create conflict handler for linear constraints */
   SCIP_CALL( SCIPincludeConflicthdlr(scip, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
         NULL, NULL, NULL, NULL, NULL, conflictExecLinear,
         NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler in SCIP */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeLinear, consInitLinear, consExitLinear,
         consInitpreLinear, consExitpreLinear, consInitsolLinear, consExitsolLinear,
         consDeleteLinear, consTransLinear, consInitlpLinear,
         consSepalpLinear, consSepasolLinear, consEnfolpLinear, consEnfopsLinear, consCheckLinear,
         consPropLinear, consPresolLinear, consRespropLinear, consLockLinear,
         consActiveLinear, consDeactiveLinear,
         consEnableLinear, consDisableLinear,
         consPrintLinear,
         conshdlrdata) );

   /* add linear constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/tightenboundsfreq",
         "multiplier on propagation frequency, how often the bounds are tightened (-1: never, 0: only at root)",
         &conshdlrdata->tightenboundsfreq, DEFAULT_TIGHTENBOUNDSFREQ, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &conshdlrdata->maxrounds, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/maxroundsroot",
         "maximal number of separation rounds per node in the root node (-1: unlimited)",
         &conshdlrdata->maxroundsroot, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/maxsepacuts",
         "maximal number of cuts separated per separation round",
         &conshdlrdata->maxsepacuts, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/maxsepacutsroot",
         "maximal number of cuts separated per separation round in the root node",
         &conshdlrdata->maxsepacutsroot, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/linear/maxpresolpairrounds",
         "maximal number of presolving rounds with pairwise constraint comparison (-1: no limit)",
         &conshdlrdata->maxpresolpairrounds, DEFAULT_MAXPRESOLPAIRROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/linear/maxaggrnormscale",
         "maximal allowed relative gain in maximum norm for constraint aggregation (0.0: disable constraint aggregation)",
         &conshdlrdata->maxaggrnormscale, DEFAULT_MAXAGGRNORMSCALE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/linear/separateall",
         "should all constraints be subject to cover cut generation instead of only the ones with non-zero dual value?",
         &conshdlrdata->separateall, DEFAULT_SEPARATEALL, NULL, NULL) );

   return SCIP_OKAY;
}

/** includes a linear constraint update method into the linear constraint handler */
SCIP_RETCODE SCIPincludeLinconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   int                   priority            /**< priority of upgrading method */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_LINCONSUPGRADE* linconsupgrade;

   assert(linconsupgd != NULL);

   /* find the linear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("linear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create a linear constraint upgrade data object */
   SCIP_CALL( linconsupgradeCreate(scip, &linconsupgrade, linconsupgd, priority) );

   /* insert linear constraint update method into constraint handler data */
   SCIP_CALL( conshdlrdataIncludeUpgrade(scip, conshdlrdata, linconsupgrade) );

   return SCIP_OKAY;
}

/** creates and captures a linear constraint */
SCIP_RETCODE SCIPcreateConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries; coefs must not be zero! */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable during node processing (subject to col generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode      /**< should the node always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the linear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("linear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, vals, lhs, rhs) );
   assert(consdata != NULL);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** adds coefficient to linear constraint (if it is not zero) */
SCIP_RETCODE SCIPaddCoefLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   if( !SCIPisZero(scip, val) )
   {
      SCIP_CALL( addCoef(scip, cons, var, val) );
   }

   return SCIP_OKAY;
}

/** gets left hand side of linear constraint */
SCIP_Real SCIPgetLhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets right hand side of linear constraint */
SCIP_Real SCIPgetRhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** changes left hand side of linear constraint */
SCIP_RETCODE SCIPchgLhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( chgLhs(scip, cons, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of linear constraint */
SCIP_RETCODE SCIPchgRhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( chgRhs(scip, cons, rhs) );

   return SCIP_OKAY;
}

/** gets the number of variables in the linear constraint */
int SCIPgetNVarsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets the array of variables in the linear constraint; the user must not modify this array! */
SCIP_VAR** SCIPgetVarsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets the array of coefficient values in the linear constraint; the user must not modify this array! */
SCIP_Real* SCIPgetValsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vals;
}

/** gets the activity of the linear constraint in the given solution */
SCIP_Real SCIPgetActivityLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIPgetRowSolActivity(scip, consdata->row, sol);
   else
      return consdataGetActivity(scip, consdata, sol);
}

/** gets the feasibility of the linear constraint in the given solution */
SCIP_Real SCIPgetFeasibilityLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIPgetRowSolFeasibility(scip, consdata->row, sol);
   else
      return consdataGetFeasibility(scip, consdata, sol);
}

/** gets the dual solution of the linear constraint in the current LP */
SCIP_Real SCIPgetDualsolLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets the dual farkas value of the linear constraint in the current infeasible LP */
SCIP_Real SCIPgetDualfarkasLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}

/** returns the linear relaxation of the given linear constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->row;
}

/** tries to automatically convert a linear constraint into a more specific and more specialized constraint */
SCIP_RETCODE SCIPupgradeConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_CONS**           upgdcons            /**< pointer to store upgraded constraint, or NULL if not successful */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real poscoeffsum;
   SCIP_Real negcoeffsum;
   SCIP_Bool integral;
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
      SCIPerrorMessage("constraint is not linear\n");
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
         SCIPerrorMessage("cannot upgrade linear constraint that is already stored as row in the LP\n");
         return SCIP_INVALIDDATA;
      }
      else
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   /* normalize constraint */
   SCIP_CALL( normalizeCons(scip, cons) );


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
            integral = integral && SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposbin++;
         else
            nnegbin++;
         break;
      case SCIP_VARTYPE_INTEGER:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposint++;
         else
            nnegint++;
         break;
      case SCIP_VARTYPE_IMPLINT:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, val);
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
         SCIPerrorMessage("unknown variable type\n");
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

   SCIPdebugMessage("upgrading linear constraint <%s> (%d upgrade methods):\n",
      SCIPconsGetName(cons), conshdlrdata->nlinconsupgrades);
   SCIPdebugMessage(" +bin=%d -bin=%d +int=%d -int=%d +impl=%d -impl=%d +cont=%d -cont=%d +1=%d -1=%d +I=%d -I=%d +F=%d -F=%d possum=%g negsum=%g integral=%d\n",
      nposbin, nnegbin, nposint, nnegint, nposimpl, nnegimpl, nposcont, nnegcont,
      ncoeffspone, ncoeffsnone, ncoeffspint, ncoeffsnint, ncoeffspfrac, ncoeffsnfrac,
      poscoeffsum, negcoeffsum, integral);

   /* try all upgrading methods in priority order */
   for( i = 0; i < conshdlrdata->nlinconsupgrades && *upgdcons == NULL; ++i )
   {
      SCIP_CALL( conshdlrdata->linconsupgrades[i]->linconsupgd(scip, cons, consdata->nvars,
            consdata->vars, consdata->vals, consdata->lhs, consdata->rhs,
            nposbin, nnegbin, nposint, nnegint, nposimpl, nnegimpl, nposcont, nnegcont,
            ncoeffspone, ncoeffsnone, ncoeffspint, ncoeffsnint, ncoeffspfrac, ncoeffsnfrac,
            poscoeffsum, negcoeffsum, integral,
            upgdcons) );
   }

#ifdef SCIP_DEBUG
   if( *upgdcons != NULL )
   {
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      SCIPdebugMessage(" -> upgraded to constraint type <%s>\n", SCIPconshdlrGetName(SCIPconsGetHdlr(*upgdcons)));
      SCIP_CALL( SCIPprintCons(scip, *upgdcons, NULL) );
   }
#endif

   return SCIP_OKAY;
}
