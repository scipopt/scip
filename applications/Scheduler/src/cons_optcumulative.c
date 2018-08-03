/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_optcumulative.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for cumulative constraints with optional activities
 * @author Chris Beck
 * @author Stefan Heinz
 *
 * Given a set of jobs \f$J\f$. Each job~\f$j\f$ has a binary variables \f$x_j\f$ which is one if this job is scheduled
 * on that machine (otherwise it is zero), an integer start time variables \f$S_j\f$, a processing time \f$p_j\f$, and a
 * demands \f$d_j\f$. Besides that an integer resource capacity \f$C\f$.
 *
 * The optcumulative enfoces the cumulative conditions for those jobs which are assigned to that machine. Let \f$J'\f$
 * be the subset of jobs assigned to that optcumulative constraint, then the cumulative constraint ensures that for
 * each point in time \f$t\f$ \f$\sum_{j\in J': S_j \leq t < S_j + p_j} d_j \leq C\f$ holds.
 *
 *
 * Propagation:
 *
 *
 * LP Relaxation:
 *
 * - let est(J) the earliest start time of all jobs of set \f$J\f$ and lct(J) the latest completion time for all jobs of
 *   set \f$J\f$, then the following linear constraint has to hold
 *   \f$\sum_{j\in J} p_j \cdot d_j \leq (lct(J) - est(J)) \cdot C\f$
 *
 */

/*
 * @todo Find subsets \f$J'\f$ of jobs which are together not schedulable and create knapsack constraint
 *       \f$\sum_{j\in J'} p_j \cdot d_j \leq (lct(J') - est(J')) \cdot C\f$
 * @todo Use a rectangle relaxation to determine if jobs which run in a certain interval can be packed feasible. this
 *       relaxation ignores the actual start and end time of a job.
 * @todo Adjsut relaxation after jobs are removed during search
 *
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "cons_optcumulative.h"

#include "scip/cons_cumulative.h"
#include "scip/cons_knapsack.h"
#include "scip/scipdefplugins.h"

/**@name Constraint handler properties
 *
 * @{
 */

/* constraint handler properties */
#define CONSHDLR_NAME          "optcumulative"
#define CONSHDLR_DESC          "constraint handler for cumulative constraints with optional activities"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -2060000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -3100000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_MEDIUM

/**@} */

/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_BINVARS_NAME  "optcumulativebinvars"
#define EVENTHDLR_BINVARS_DESC  "bound change event handler for binary variables of optcumulative constraints"

#define EVENTHDLR_INTVARS_NAME  "optcumulativeintvars"
#define EVENTHDLR_INTVARS_DESC  "bound change event handler for integer variables of optcumulative constraints"

/**@} */

/**@name Default parameter values
 *
 * @{
 */

#define DEFAULT_ROWRELAX          FALSE /**< add linear relaxation as LP row (otherwise a knapsack constraint is created)? */
#define DEFAULT_CONFLICTANALYSIS   TRUE /**< participate in conflict analysis?" */
#define DEFAULT_INTERVALRELAX      TRUE /**< create a relaxation for each start and end time point interval */

/**@} */


/*
 * Data structures
 */

/** constraint data for optcumulative constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< array of variable representing the start time of each job */
   SCIP_VAR**            binvars;            /**< array of variable representing if the job has to be processed on this machine */
   SCIP_Bool*            downlocks;          /**< array to store if the start time variable has a down lock */
   SCIP_Bool*            uplocks;            /**< array to store if the start time variable has an up lock */
   SCIP_ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   SCIP_CONS*            cons;               /**< knapsack relaxation, if created */
   int*                  demands;            /**< array containing corresponding demands */
   int*                  durations;          /**< array containing corresponding durations */
   int                   nvars;              /**< number of variables */
   int                   varssize;           /**< number of available slots in variable arrays */
   int                   capacity;           /**< available cumulative capacity */

   int                   hmin;               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax;               /**< right bound of time axis to be considered  (not including hmax) */

   int                   nglbfixedzeros;     /**< number of binary variable globally fixed to zero */
   int                   nglbfixedones;      /**< number of binary variable globally fixed to one */
   int                   nfixedzeros;        /**< number of binary variable fixed to zero */
   int                   nfixedones;         /**< number of binary variable fixed to one */
   int                   est;                /**< used earliest start time for the relaxation */
   int                   lct;                /**< used latest completion time for the relaxation */
   unsigned int          propagated:1;       /**< is constraint already propagated? */
   unsigned int          relaxadded:1;       /**< was relaxation added? */
   unsigned int          triedsolving:1;     /**< bool to store if it was tried to solve the cumulative sub-problem */
   unsigned int          normalized:1;       /**< is the constraint normalized */
   unsigned int          triedredundant:1;   /**< bool to store if the redundancy check was applied */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlrbinvars;   /**< event handler for bound change events on binary variables */
   SCIP_EVENTHDLR*       eventhdlrintvars;   /**< event handler for bound change events on integer variables */
   SCIP_HEUR*            heurtrysol;         /**< trysol heuristic */
   SCIP_Bool             rowrelax;           /**< add linear relaxation as LP row (otherwise a knapsack constraint is created)? */
   SCIP_Bool             conflictanalysis;   /**< participate in conflict analysis? */
   SCIP_Bool             intervalrelax;      /**< create a relaxation for each start and end time point interval */
};

/**@name Debug Methods
 *
 */

#ifndef NDEBUG
/** check constraint state (nglbfixedones and nglbfixedzeros) */
static
void checkCounters(
   SCIP_CONSDATA*        consdata            /**< optcumulative constraint data */
   )
{
   int nglbfixedones;
   int nglbfixedzeors;
   int nfixedones;
   int nfixedzeors;
   int v;

   nglbfixedones = 0;
   nglbfixedzeors = 0;
   nfixedones = 0;
   nfixedzeors = 0;

   for( v = 0; v < consdata->nvars; ++v )
   {
      if( SCIPvarGetLbGlobal(consdata->binvars[v]) > 0.5 )
         nglbfixedones++;

      if( SCIPvarGetUbGlobal(consdata->binvars[v]) < 0.5 )
         nglbfixedzeors++;

      if( SCIPvarGetLbLocal(consdata->binvars[v]) > 0.5 )
         nfixedones++;

      if( SCIPvarGetUbLocal(consdata->binvars[v]) < 0.5 )
         nfixedzeors++;
   }

   assert(nglbfixedones == consdata->nglbfixedones);
   assert(nglbfixedzeors == consdata->nglbfixedzeros);
   assert(nfixedones == consdata->nfixedones);
   assert(nfixedzeors == consdata->nfixedzeros);
}
#else
#define checkCounters(x) /* */
#endif

/**@} */

/**@name Miscellaneous Methods
 *
 * @{
 */

#ifndef NDEBUG
/** converts the given double bound which is integral to an int; in optimized mode the function gets inlined for
 *  performance; in debug mode we check some additional conditions
 */
static
int convertBoundToInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bound               /**< double bound to convert */
   )
{
   assert(SCIPisIntegral(scip, bound));
   assert(SCIPisEQ(scip, bound, (SCIP_Real)(int)(bound + 0.5)));

   return (int)(bound + 0.5);
}
#else
#define convertBoundToInt(x, y) ((int)((y) + 0.5))
#endif

/**@} */

/**@name Constraint data methods
 *
 * @{
 */

/** creates constraint data of optcumulative constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to consdata */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< array of integer variables */
   SCIP_VAR**            binvars,            /**< array of variable representing if the job has to be processed on this machine */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   SCIP_Bool             check               /**< is the corresponding constraint a check constraint */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(vars != NULL || nvars > 0);
   assert(binvars != NULL || nvars > 0);
   assert(demands != NULL);
   assert(durations != NULL);
   assert(capacity >= 0);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->capacity = capacity;
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->hmin = 0;
   (*consdata)->hmax = INT_MAX;
   (*consdata)->nglbfixedzeros = 0;
   (*consdata)->est = -1;
   (*consdata)->lct = INT_MAX;
   (*consdata)->row = NULL;
   (*consdata)->cons = NULL;
   (*consdata)->nglbfixedzeros = 0;
   (*consdata)->nglbfixedones = 0;
   (*consdata)->nfixedzeros = 0;
   (*consdata)->nfixedones = 0;
   (*consdata)->propagated = FALSE;
   (*consdata)->relaxadded = FALSE;
   (*consdata)->triedsolving = FALSE;
   (*consdata)->normalized = FALSE;
   (*consdata)->triedredundant = FALSE;

   if( nvars > 0 )
   {
      int v;

      assert(vars != NULL); /* for flexelint */
      assert(binvars != NULL); /* for flexelint */

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->binvars, binvars, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->downlocks, demands, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->uplocks, demands, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->demands, demands, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->durations, durations, nvars) );

      /* initialize locking arrays */
      for( v = 0; v < nvars; ++v )
      {
         /* the locks are only used if the contraint is a check constraint */
         (*consdata)->downlocks[v] = check;
         (*consdata)->uplocks[v] = check;
      }

      /* transform variables, if they are not yet transformed */
      if( SCIPisTransformed(scip) )
      {
         SCIPdebugMessage("get tranformed variables and constraints\n");

         /* get transformed variables and do NOT captures these */
         SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
         SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->binvars, (*consdata)->binvars) );

         for( v = 0; v < nvars; ++v )
         {
            SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->vars[v]) );
            SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->binvars[v]) );
         }
      }
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->binvars = NULL;
      (*consdata)->downlocks = NULL;
      (*consdata)->uplocks = NULL;
      (*consdata)->demands = NULL;
      (*consdata)->durations = NULL;
   }

   return SCIP_OKAY;
}


/** frees a optcumulative constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to linear constraint data */
   )
{
   int varssize;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* release the row */
   if( (*consdata)->cons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*consdata)->cons) );
   }

   varssize = (*consdata)->varssize;

   if( varssize > 0 )
   {
      /* free arrays */
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->durations, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->demands, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->uplocks, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->downlocks, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->binvars, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, varssize);
   }

   /* free memory */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints optcumulative constraint to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< optcumulative constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);

   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      if( v > 0 )
         SCIPinfoMessage(scip, file, ", ");

      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->vars[v], FALSE) );

      SCIPinfoMessage(scip, file, "[%g,%g](%d,%d)", SCIPvarGetLbLocal(consdata->vars[v]),
         SCIPvarGetUbLocal(consdata->vars[v]), consdata->durations[v], consdata->demands[v]);

      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->binvars[v], FALSE) );

   }
   SCIPinfoMessage(scip, file, " [%d,%d]<= %d", consdata->hmin, consdata->hmax, consdata->capacity);

   return SCIP_OKAY;
}

/**@} */

/**@name Constraint handler data
 *
 * Method used to create and free the constraint handler data when including and removing the cumulative constraint
 * handler.
 *
 * @{
 */

/** creates constaint handler data for set partitioning / packing / covering constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata,       /**< pointer to store the constraint handler data */
   SCIP_EVENTHDLR*       eventhdlrbinvars,   /**< used event handler for tracing bound changes on binary variables */
   SCIP_EVENTHDLR*       eventhdlrintvars    /**< used event handler for tracing bound changes on integer variables */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(eventhdlrbinvars != NULL);
   assert(eventhdlrintvars != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, conshdlrdata) );

   (*conshdlrdata)->eventhdlrbinvars = eventhdlrbinvars;
   (*conshdlrdata)->eventhdlrintvars = eventhdlrintvars;
   (*conshdlrdata)->heurtrysol = NULL;

   return SCIP_OKAY;
}

/** frees constraint handler data for set partitioning / packing / covering constraint handler */
static
SCIP_RETCODE conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/**@} */

/** removes rounding locks for the given variable in the given optcumulative constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< optcumulative constraint */
   SCIP_VAR*             binvar,             /**< decision variable */
   SCIP_VAR*             var,                /**< start time variable */
   SCIP_Bool             downlock,           /**< has the integer start time variable a down lock */
   SCIP_Bool             uplock              /**< has the integer start time variable an up lock */
   )
{
   /* rounding up may violate the constraint */
   SCIP_CALL( SCIPunlockVarCons(scip, binvar, cons, FALSE, TRUE) );

   /* rounding in both directions may violate the constraint */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, downlock, uplock) );

   return SCIP_OKAY;
}

/** catches events for binary variable at given position */
static
SCIP_RETCODE catchEventBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR* binvar;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->binvars != NULL);

   binvar = consdata->binvars[pos];
   assert(binvar != NULL);

   /* we are catching the following events for the binary variables:
    *
    * - SCIP_EVENTTYPE_BOUNDRELAXED: This allows for counting locally fixed variables to one or zero
    * - SCIP_EVENTTYPE_GBDCHANGED: This allows to check if the optcumulative can be converted into an cumulative
    *   constraint
    * - SCIP_EVENTTYPE_BOUNDRELAXED: This allows us to detect the moment when we can retry to solve a local cumulative
    *   constraint again
    */
   eventtype = SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_BOUNDRELAXED;

   /* catch bound change events on variable */
   SCIP_CALL( SCIPcatchVarEvent(scip, binvar, eventtype, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

   /* update the globally fixed variables counter for this variable */
   if( SCIPvarGetUbGlobal(binvar) < 0.5)
      consdata->nglbfixedzeros++;
   else if( SCIPvarGetLbGlobal(binvar) > 0.5 )
      consdata->nglbfixedones++;

   /* update the locally fixed variables counter for this variable */
   if( SCIPvarGetUbLocal(binvar) < 0.5)
      consdata->nfixedzeros++;
   else if( SCIPvarGetLbLocal(binvar) > 0.5 )
      consdata->nfixedones++;

   assert(consdata->nglbfixedzeros + consdata->nglbfixedones <= consdata->nvars);
   assert(consdata->nfixedzeros + consdata->nfixedones <= consdata->nvars);

   return SCIP_OKAY;
}

/** drops events for binary variable at given position */
static
SCIP_RETCODE dropEventBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR* binvar;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->binvars != NULL);

   binvar = consdata->binvars[pos];
   assert(binvar != NULL);

   eventtype = SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_BOUNDRELAXED;

   /* drop events on variable */
   SCIP_CALL( SCIPdropVarEvent(scip, binvar, eventtype, eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

   /* update the globally fixed variables counter for this variable */
   if( SCIPvarGetUbGlobal(binvar) < 0.5)
      consdata->nglbfixedzeros--;
   else if( SCIPvarGetLbGlobal(binvar) > 0.5 )
      consdata->nglbfixedones--;

   /* update the locally fixed variables counter for this variable */
   if( SCIPvarGetUbLocal(binvar) < 0.5)
      consdata->nfixedzeros--;
   else if( SCIPvarGetLbLocal(binvar) > 0.5 )
      consdata->nfixedones--;

   assert(consdata->nglbfixedzeros >= 0);
   assert(consdata->nglbfixedones >= 0);
   assert(consdata->nfixedzeros >= 0);
   assert(consdata->nfixedones >= 0);

   return SCIP_OKAY;
}

/** catches events for integer variable at given position */
static
SCIP_RETCODE catchEventIntvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR* var;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars != NULL);

   var = consdata->vars[pos];
   assert(var != NULL);

   /* we are catching the following events for the integer variables:
    *
    * - SCIP_EVENTTYPE_GBDCHANGED: This allows to check if the optcumulative can be converted into an cumulative
    *   constraint
    */
   eventtype = SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_GBDCHANGED;

   /* catch bound change events on variable */
   SCIP_CALL( SCIPcatchVarEvent(scip, var, eventtype, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

   return SCIP_OKAY;
}

/** drops events for integer variable at given position */
static
SCIP_RETCODE dropEventIntvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR* var;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars != NULL);

   var = consdata->vars[pos];
   assert(var != NULL);

   eventtype = SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_GBDCHANGED;

   /* drop events on variable */
   SCIP_CALL( SCIPdropVarEvent(scip, var, eventtype, eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

   return SCIP_OKAY;
}

/** catches bound change events for all variables in transformed optcumulative constraint */
static
SCIP_RETCODE catchAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   SCIP_EVENTHDLR*       eventhdlrbinvars,   /**< event handler to call for the event processing on binary variables */
   SCIP_EVENTHDLR*       eventhdlrintvars    /**< event handler to call for the event processing on integer variables */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nglbfixedones == 0);
   assert(consdata->nglbfixedzeros == 0);
   assert(consdata->nfixedones == 0);
   assert(consdata->nfixedzeros == 0);

   /* catch event for every single variable */
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_CALL( catchEventBinvar(scip, cons, eventhdlrbinvars, v) );

      SCIP_CALL( catchEventIntvar(scip, cons, eventhdlrintvars, v) );
   }

   /* (debug) check if the counter of the constraint are correct */
   checkCounters(consdata);

   return SCIP_OKAY;
}

/** drops bound change events for all variables in transformed optcumulative constraint */
static
SCIP_RETCODE dropAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< set partitioning / packing / covering constraint */
   SCIP_EVENTHDLR*       eventhdlrbinvars,   /**< event handler to call for the event processing on binary variables */
   SCIP_EVENTHDLR*       eventhdlrintvars    /**< event handler to call for the event processing on integer variables */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* drop event of every single variable */
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_CALL( dropEventBinvar(scip, cons, eventhdlrbinvars, v) );

      SCIP_CALL( dropEventIntvar(scip, cons, eventhdlrintvars, v) );
   }

   /* check that the internal constraint state is rested */
   assert(consdata->nglbfixedones == 0);
   assert(consdata->nglbfixedzeros == 0);
   assert(consdata->nfixedones == 0);
   assert(consdata->nfixedzeros == 0);

   return SCIP_OKAY;
}

/** initialize the sorted event point arrays */
static
void createSortedEventpoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int*                  starttimes,         /**< array to store sorted start events */
   int*                  endtimes,           /**< array to store sorted end events */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int*                  endindices,         /**< permutation with rspect to the end times */
   SCIP_Bool             local               /**< shall local bounds be used */
   )
{
   SCIP_VAR* var;
   int nvars;
   int j;

   nvars = consdata->nvars;

   /* assign variables, start and endpoints to arrays */
   for ( j = 0; j < nvars; ++j )
   {
      var = consdata->vars[j];
      if( local )
         starttimes[j] = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      else
         starttimes[j] = convertBoundToInt(scip, SCIPvarGetLbGlobal(var));

      startindices[j] = j;

      if( local )
         endtimes[j] = convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + consdata->durations[j];
      else
         endtimes[j] = convertBoundToInt(scip, SCIPvarGetUbGlobal(var)) + consdata->durations[j];

      endindices[j] = j;
   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues (and sort the indices in the same way) */
   SCIPsortIntInt(starttimes, startindices, nvars);
   SCIPsortIntInt(endtimes, endindices, nvars);
}

/** computes the maximum energy for all variables which correspond to jobs which start between the given start time and
 *  end time
 *
 *  @return Maximum energy for the given time window
 */
static
SCIP_Longint computeMaxEnergy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< optcumulative constraint data */
   int                   starttime,          /**< start time */
   int                   endtime             /**< end time */
   )
{
   SCIP_VAR* var;
   SCIP_Longint maxenergy;
   int v;

   assert(starttime < endtime);
   maxenergy = 0LL;

   for( v = 0; v < consdata->nvars; ++v )
   {
      var = consdata->vars[v];

      /* collect jobs which run between the start and end time */
      if( convertBoundToInt(scip, SCIPvarGetUbGlobal(var)) + consdata->durations[v] <= endtime
         && convertBoundToInt(scip, SCIPvarGetLbGlobal(var)) >= starttime)
      {
         maxenergy += (SCIP_Longint)(consdata->durations[v] * consdata->demands[v]); /*lint !e647*/
      }
   }

   return maxenergy;
}

/** collects all variables which correspond to jobs which start between the given start time and end time */
static
SCIP_RETCODE collectVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< optcumulative constraint data */
   SCIP_VAR**            vars,               /**< array to store the variables */
   SCIP_Longint*         weights,            /**< array to store the weights */
   int*                  nvars,              /**< pointer to store the number of collected variables */
   int                   starttime,          /**< start time */
   int                   endtime             /**< end time */
   )
{
   SCIP_VAR* var;
   int v;

   assert(starttime < endtime);
   (*nvars) = 0;

   for( v = 0; v < consdata->nvars; ++v )
   {
      var = consdata->vars[v];

      /* collect jobs which run between the start and end time */
      if( convertBoundToInt(scip, SCIPvarGetUbGlobal(var)) + consdata->durations[v] <= endtime
         && convertBoundToInt(scip, SCIPvarGetLbGlobal(var)) >= starttime)
      {
         vars[*nvars] = consdata->binvars[v];
         weights[*nvars] = (SCIP_Longint)(consdata->durations[v] * consdata->demands[v]); /*lint !e647*/
         (*nvars)++;
      }
   }

   return SCIP_OKAY;
}

/** remove row which have a tightness which is smaller or equal to the given one
 *
 *  @return The number of remaining rows
 */
static
int removeRedundantRows(
   SCIP_Longint*         rowtightness,       /**< array containing the tightness for the previously selected rows */
   int*                  startidxs,          /**< array containing for each row the index for the start event */
   int                   nrows,              /**< current number of rows */
   SCIP_Longint          tightness           /**< tightness to use to detect redundant rows */
   )
{
   int keptrows;
   int j;

   keptrows = 0;

   for( j = 0; j < nrows; ++j )
   {
      rowtightness[keptrows] = rowtightness[j];
      startidxs[keptrows] = startidxs[j];

      /* only keep this row if the tightness is better as the (current) given one */
      if( rowtightness[j] > tightness )
         keptrows++;
   }

   return keptrows;
}

/** depending on the parameters setting a row or an knapsack constraint is created */
static
SCIP_RETCODE createRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   const char*           name,               /**< name of the row */
   SCIP_VAR**            vars,               /**< array of variable representing if the job has to be processed on this machine */
   SCIP_Longint*         weights,            /**< start time variables of the activities which are assigned */
   int                   nvars,              /**< number of variables */
   SCIP_Longint          capacity,           /**< available cumulative capacity */
   SCIP_Bool             local,              /**< create local row */
   SCIP_Bool*            rowadded,           /**< pointer to store if a row was added */
   SCIP_Bool*            consadded,          /**< pointer to store if a constraint was added */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff occurred */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *cutoff = FALSE;
   if( conshdlrdata->rowrelax || SCIPgetDepth(scip) > 0 )
   {
      SCIP_ROW* row;
      int v;

      /* create empty row */
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, name, -SCIPinfinity(scip), (SCIP_Real)capacity, local, FALSE, FALSE) );

      /* w.r.t. performance we cache the row extension and flush them in the end */
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      for( v = 0; v < nvars; ++v )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[v], (SCIP_Real)weights[v]) );
      }

      /* w.r.t. performance we flush the row extension in the end */
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );

      assert(!SCIProwIsInLP(row));

      if( SCIPgetDepth(scip) == 0 || SCIPisCutEfficacious(scip, NULL, row) )
      {
         SCIPdebug( SCIPprintRow(scip, row, NULL) );
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
         (*rowadded) = TRUE;
      }

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }
   else
   {
      SCIP_CONS* cons;

      /* create knapsack constraint */
      SCIP_CALL( SCIPcreateConsKnapsack(scip, &cons, name, nvars, vars, weights, capacity,
            FALSE, TRUE, TRUE, FALSE, TRUE, local, FALSE, FALSE, TRUE, FALSE) );

      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );

      /* add and releasse knapsack constraint */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      (*consadded) = TRUE;
   }

   return SCIP_OKAY;
}

/** adds linear relaxation as cut to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< optcumulative constraint */
   SCIP_Bool*            rowadded,           /**< pointer to store if a row was added */
   SCIP_Bool*            consadded,          /**< pointer to store if a constraint was added */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff occurred */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert( cutoff != NULL );

   *cutoff = FALSE;
   if( consdata->relaxadded )
      return SCIP_OKAY;

   SCIPdebugMessage("add relaxation for optcumulative constraint <%s>\n", SCIPconsGetName(cons));

   if( conshdlrdata->intervalrelax )
   {
      SCIP_Longint** rowtightness;
      int** startidxs;
      int* nrows;
      int* starttimes;
      int* endtimes;
      int* startindices;
      int* endindices;
      int starttime;
      int endtime;
      int i;
      int j;

      SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &startindices, consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &endindices, consdata->nvars) );

      SCIP_CALL( SCIPallocBufferArray(scip, &nrows, consdata->nvars) );
      BMSclearMemoryArray(nrows, consdata->nvars);
      SCIP_CALL( SCIPallocBufferArray(scip, &rowtightness, consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &startidxs, consdata->nvars) );
      for( j = 0; j < consdata->nvars; ++j )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &rowtightness[j], consdata->nvars) ); /*lint !e866*/
         SCIP_CALL( SCIPallocBufferArray(scip, &startidxs[j], consdata->nvars) ); /*lint !e866*/
      }

      createSortedEventpoints(scip, consdata, starttimes, endtimes, startindices, endindices, TRUE);

      starttime = -INT_MAX;

      /* check each startpoint of a job whether the capacity is kept or not */
      for( j = 0; j < consdata->nvars; ++j )
      {
         SCIP_Longint besttightness;

         assert(starttime <= starttimes[j]);

         /* if we hit the same start time again we skip the loop */
         if( starttime == starttimes[j])
            continue;

         starttime = starttimes[j];
         endtime = -INT_MAX;
         besttightness = 0LL;

         for( i = 0; i < consdata->nvars; ++i )
         {
            SCIP_Longint energy;
            SCIP_Longint maxenergy;
            SCIP_Longint tightness;

            assert(endtime <= endtimes[i]);

            /* if we hit the same end time again we skip the loop */
            if( endtime == endtimes[i] )
               continue;

            endtime = endtimes[i];

            /* skip all end times which are smaller than the start time */
            if( endtime <= starttime )
               continue;

            maxenergy = computeMaxEnergy(scip, consdata, starttime, endtime);

            energy = (endtime - starttime) * consdata->capacity; /*lint !e647*/
            tightness = maxenergy - energy;

            /* check if the linear constraint is not trivially redundant */
            if( tightness > besttightness )
            {
               besttightness = tightness;

               nrows[i] = removeRedundantRows(rowtightness[i], startidxs[i], nrows[i], tightness);

               /* add row information */
               rowtightness[i][nrows[i]] = tightness;
               startidxs[i][nrows[i]] = j;
               nrows[i]++;
            }
         }
      }

      for( j = consdata->nvars-1; j >= 0 && ! (*cutoff); --j )
      {
         for( i = 0; i < nrows[j] && ! (*cutoff); ++i )
         {
            SCIP_VAR** vars;
            SCIP_Longint* weights;
            SCIP_Longint energy;
            char name[SCIP_MAXSTRLEN];
            int nvars;

            SCIP_CALL( SCIPallocBufferArray(scip, &vars, consdata->nvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &weights, consdata->nvars) );

            starttime = starttimes[startidxs[j][i]];
            endtime = endtimes[j];

            energy = (endtime - starttime) * consdata->capacity; /*lint !e647*/

            SCIP_CALL( collectVars(scip, consdata, vars, weights, &nvars, starttime, endtime) );

            SCIPdebugMessage("create linear relaxation for <%s> time interval [%d,%d] <= %"SCIP_LONGINT_FORMAT" (tightness %"SCIP_LONGINT_FORMAT")\n",
               SCIPconsGetName(cons), starttime, endtime, energy, rowtightness[j][i]);

            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s[%d,%d]", SCIPconsGetName(cons), starttime, endtime);
            SCIP_CALL( createRow(scip, conshdlr, name, vars, weights, nvars, energy, TRUE, rowadded, consadded, cutoff) );

            SCIPfreeBufferArray(scip, &weights);
            SCIPfreeBufferArray(scip, &vars);
         }
      }

      /* free buffers */
      for( j = consdata->nvars-1; j >= 0; --j )
      {
         SCIPfreeBufferArray(scip, &startidxs[j]);
         SCIPfreeBufferArray(scip, &rowtightness[j]);
      }
      SCIPfreeBufferArray(scip, &startidxs);
      SCIPfreeBufferArray(scip, &rowtightness);
      SCIPfreeBufferArray(scip, &nrows);

      SCIPfreeBufferArray(scip, &endindices);
      SCIPfreeBufferArray(scip, &endtimes);
      SCIPfreeBufferArray(scip, &startindices);
      SCIPfreeBufferArray(scip, &starttimes);
   }
   else
   {
      SCIP_VAR** vars;
      SCIP_Longint* weights;
      SCIP_Longint maxenergy;
      SCIP_Longint energy;
      int* durations;
      int* demands;
      int est;
      int lct;
      int nvars;
      int v;

      nvars = consdata->nvars;
      vars = consdata->vars;
      durations = consdata->durations;
      demands = consdata->demands;
      maxenergy = 0LL;

      SCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );

      est = INT_MAX;
      lct = 0;

      for( v = 0; v < nvars; ++v )
      {
         weights[v] = (SCIP_Longint)(durations[v] * demands[v]); /*lint !e647*/
         maxenergy += weights[v];

         /* adjust earlier start time */
         est = MIN(est, convertBoundToInt(scip, SCIPvarGetLbLocal(vars[v]))); /*lint !e666*/

         /* adjust latest completion */
         lct = MAX(lct, convertBoundToInt(scip, SCIPvarGetUbLocal(vars[v]) + durations[v])); /*lint !e666*/
      }

      energy = (lct - est) * consdata->capacity; /*lint !e647*/

      if( maxenergy > energy )
      {
         char name[SCIP_MAXSTRLEN];

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s[%d,%d]", SCIPconsGetName(cons), est, lct);

         SCIPdebugMessage("create linear relaxation for <%s> (nvars %d) time interval [%d,%d] <= %"SCIP_LONGINT_FORMAT"\n",
            SCIPconsGetName(cons), nvars, est, lct, energy);

         SCIP_CALL( createRow(scip, conshdlr, name, consdata->binvars, weights, nvars, energy, TRUE, rowadded, consadded, cutoff) );
      }

      /* free buffer */
      SCIPfreeBufferArray(scip, &weights);
   }

   consdata->relaxadded = TRUE;

#if 0
   if( !conshdlrdata->rowrelax )
   {
      SCIP_CALL( SCIPrestartSolve(scip) );
   }
#endif

   return SCIP_OKAY;
}

/** collect all activities which are locally (that means in the current branch and bound node) assigned to that
 *  machine
 */
static
void collectActivities(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR**            binvars,            /**< array of variable representing if the job has to be processed on this machine */
   SCIP_VAR**            vars,               /**< start time variables of the activities which are assigned */
   int*                  durations,          /**< durations of the activities */
   int*                  demands,            /**< demands of the activities */
   int*                  nfixedones,         /**< pointer to store number of activities assigned to that machine */
   int*                  nfixedzeros,        /**< pointer to store number of binary variables fixed to zeor */
   SCIP_Bool*            auxiliary           /**< pointer to store if the integer start time variables of the assigned
                                              *   activities are auxiliary variables; that is the case if the optcumulative
                                              *   choice constraints is the only one having locks on these variables */
   )
{
   int v;

   /* collect all jobs which have to be processed */
   (*auxiliary) = TRUE;
   (*nfixedones) = 0;
   (*nfixedzeros) = 0;

   for( v = 0; v < consdata->nvars; ++v )
   {
      if( SCIPvarGetLbLocal(consdata->binvars[v]) > 0.5 )
      {
         /* binary variable is fixed one */

         SCIPdebugMessage("collect variable <%s>[%g,%g](%d)\n",
            SCIPvarGetName(consdata->vars[v]), SCIPvarGetLbLocal(consdata->vars[v]), SCIPvarGetUbGlobal(consdata->vars[v]), consdata->durations[v]);

         binvars[*nfixedones] = consdata->binvars[v];
         vars[*nfixedones] = consdata->vars[v];
         durations[*nfixedones] = consdata->durations[v];
         demands[*nfixedones] = consdata->demands[v];

         (*nfixedones)++;

         /* check the locks on the integer start time variable to determine if its a auxiliary variable (only locked by
          * this constraint)
          */
         if( SCIPvarGetNLocksDown(consdata->vars[v]) > (int)consdata->downlocks[v]
            || SCIPvarGetNLocksUp(consdata->vars[v]) > (int)consdata->uplocks[v] )
         {
            (*auxiliary) = FALSE;
         }
      }
      else if( SCIPvarGetUbLocal(consdata->binvars[v]) < 0.5 )
         (*nfixedzeros)++;
   }

   assert(consdata->nfixedzeros == *nfixedzeros);
   assert(consdata->nfixedones == *nfixedones);
}

/** collect all activities which are assigned to that machine in the given solution */
static
void collectSolActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   SCIP_VAR**            binvars,            /**< array of variable representing if the job has to be processed on this machine */
   SCIP_VAR**            vars,               /**< start time variables of the activities which are assigned */
   int*                  durations,          /**< durations of the activities */
   int*                  demands,            /**< demands of the activities */
   int*                  nvars,              /**< pointer to store number of activities assigned to that machine */
   int*                  nfixedones,         /**< pointer to store number of binary variables locally fixed to one */
   int*                  nfixedzeros,        /**< pointer to store number of binary variables locally fixed to zero */
   SCIP_Bool*            auxiliary           /**< pointer to store if the integer start time variables of the assigned
                                              *   activities are auxiliary variables; that is the case if the machine
                                              *   choice constraints is the only one haveing locks on these variables */
   )
{
   int v;

   (*nvars) = 0;
   (*nfixedones) = 0;
   (*nfixedzeros) = 0;
   (*auxiliary) = TRUE;

   /* collect all jobs which have to be processed */
   for( v = 0; v < consdata->nvars; ++v )
   {
      if( SCIPgetSolVal(scip, sol, consdata->binvars[v]) > 0.5 )
      {
         SCIPdebugMessage("collect variable <%s>\n", SCIPvarGetName(consdata->vars[v]));
         binvars[*nvars] = consdata->binvars[v];
         vars[*nvars] = consdata->vars[v];
         durations[*nvars] = consdata->durations[v];
         demands[*nvars] = consdata->demands[v];
         (*nvars)++;

         /* check the locks on the integer start time variable to determine if its a auxiliary variable */
         if( SCIPvarGetNLocksDown(consdata->vars[v]) > (int)consdata->downlocks[v]
            || SCIPvarGetNLocksUp(consdata->vars[v]) > (int)consdata->uplocks[v]
            )
            (*auxiliary) = FALSE;
      }

      if( SCIPvarGetLbLocal(consdata->binvars[v]) > 0.5 )
         nfixedones++;
      else if( SCIPvarGetUbLocal(consdata->binvars[v]) < 0.5 )
         nfixedzeros++;
   }
}

/** solves given cumulative condition as independent sub problem
 *
 *  @note The time and memory limit of the SCIP environment in transferred to sub solver
 *
 *  @note If the problem was solved to the earliest start times (ests) and latest start times (lsts) array contain the
 *        solution values; If the problem was not solved these two arrays contain the global bounds at the time the sub
 *        solver was interrupted.
 */
static
SCIP_RETCODE solveCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_Bool             local,              /**< use local bounds, otherwise global */
   SCIP_Real*            ests,               /**< array to store the earlier start time for each job */
   SCIP_Real*            lsts,               /**< array to store the latest start time for each job */
   SCIP_Longint          maxnodes,           /**< maximum number of branch-and-bound nodes to solve the single cumulative constraint  (-1: no limit) */
   SCIP_Bool*            solved,             /**< pointer to store if the problem is solved (to optimality) */
   SCIP_Bool*            infeasible,         /**< pointer to store if the problem is infeasible */
   SCIP_Bool*            unbounded,          /**< pointer to store if the problem is unbounded */
   SCIP_Bool*            error               /**< pointer to store if an error occurred */
   )
{
   SCIP_Real* objvals;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   int v;

   SCIP_CALL( SCIPallocBufferArray(scip, &objvals, nvars) );

   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;

      var = vars[v];
      assert(var != NULL);

      if( local )
      {
         ests[v] = SCIPvarGetLbLocal(var);
         lsts[v] = SCIPvarGetUbLocal(var);
      }
      else
      {
         ests[v] = SCIPvarGetLbGlobal(var);
         lsts[v] = SCIPvarGetUbGlobal(var);
      }

      objvals[v] = SCIPvarGetObj(var);
   }

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( !SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   SCIP_CALL( SCIPsolveCumulative(scip, nvars, ests, lsts, objvals,  durations, demands,
         capacity, hmin, hmax, timelimit, memorylimit, maxnodes,
         solved, infeasible, unbounded, error) );

   SCIPfreeBufferArray(scip, &objvals);

   return SCIP_OKAY;
}


/** create a logicor constraint which ensures that the jobs related to binary variables are not assigned in the same
 *  time to this optional cumulative constraint
 */
static
SCIP_RETCODE createConflictCons(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of conflict constraint */
   SCIP_VAR**            binvars,            /**< array of binary variables */
   int                   nvars               /**< number of variables */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* negatedvar;
   int v;

   /* one of the jobs cannot be processed on that resource */
   SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, 0, NULL,
         TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarGetLbGlobal(binvars[v]) > 0.5 )
         continue;

      SCIP_CALL( SCIPgetNegatedVar(scip, binvars[v], &negatedvar) );

      SCIP_CALL( SCIPaddCoefLogicor(scip, cons, negatedvar) );
   }

   /* add and release to constraint */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/** check of the given constraint is redundant */
static
SCIP_RETCODE checkRedundancy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< optcumulative constraint which collapsed to a cumulative constraint locally */
   int*                  ndelconss,          /**< pointer to store the number of deleted constraints */
   SCIP_Bool*            redundant           /**< pointer to store if the constraint is redundant */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool solved;
   SCIP_Bool infeasible;
   SCIP_Bool unbounded;
   SCIP_Bool error;
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(!SCIPinProbing(scip));

   (*redundant) = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nglbfixedzeros == 0);

   if( consdata->triedredundant )
      return SCIP_OKAY;

   consdata->triedredundant = TRUE;

   nvars = consdata->nvars;

   /* check the locks on the integer start time variable to determine if its a auxiliary variable */
   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarGetNLocksDown(consdata->vars[v]) > (int)consdata->downlocks[v]
         || SCIPvarGetNLocksUp(consdata->vars[v]) > (int)consdata->uplocks[v]
         )
         return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) );

   /* solve the cumulative condition separately */
   SCIP_CALL( solveCumulative(scip, nvars, consdata->vars, consdata->durations, consdata->demands,
         consdata->capacity, consdata->hmin, consdata->hmax, FALSE,
         lbs, ubs, 2000LL, &solved, &infeasible, &unbounded, &error) );
   assert(!unbounded);

   if( !error )
   {
      if( infeasible )
      {
         SCIP_VAR** binvars;
         SCIP_VAR** vars;
         int* durations;
         int* demands;
         SCIP_Real* weights;

         SCIP_CALL( SCIPallocBufferArray(scip, &binvars, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &durations, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &demands, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );

         for( v = 0; v < nvars; ++v )
         {
            SCIP_VAR* var;
            int est;
            int lst;

            var = consdata->vars[v];
            assert(var != NULL);

            est = convertBoundToInt(scip, SCIPvarGetLbGlobal(var));
            lst = convertBoundToInt(scip, SCIPvarGetUbGlobal(var));

            if( consdata->demands[v] == 0.0 || consdata->durations[v] == 0.0 )
               return SCIP_ERROR;

            weights[v] = (lst - est) / (consdata->demands[v] * consdata->durations[v]); /*lint !e653*/

            binvars[v] = consdata->binvars[v];
            vars[v] = var;
            durations[v] = consdata->durations[v];
            demands[v] = consdata->demands[v];
         }
         SCIPsortRealPtrPtrIntInt(weights, (void*)binvars, (void*)vars, durations, demands, nvars);

         while( nvars > 1 )
         {
            SCIP_CALL( solveCumulative(scip, nvars-1, vars, consdata->durations, consdata->demands, consdata->capacity, consdata->hmin, consdata->hmax, TRUE,
                  lbs, ubs, 2000LL, &solved, &infeasible, &unbounded, &error) );

            if( !infeasible )
               break;

            nvars--;
         }

         SCIP_CALL( createConflictCons(scip, SCIPconsGetName(cons), binvars,  nvars) );

         SCIPfreeBufferArray(scip, &weights);
         SCIPfreeBufferArray(scip, &demands);
         SCIPfreeBufferArray(scip, &durations);
         SCIPfreeBufferArray(scip, &vars);
         SCIPfreeBufferArray(scip, &binvars);
      }
      else if( solved )
      {
         for( v = 0; v < nvars; ++v )
         {
            SCIP_VAR* var;

            /* check if variable is fixed */
            assert(lbs[v] + 0.5 > ubs[v]);

            var = consdata->vars[v];
            assert(var != NULL);

            if( SCIPvarGetLbGlobal(var) + 0.5 < lbs[v] )
            {
               SCIP_CALL( SCIPchgVarLbGlobal(scip, var, lbs[v]) );
            }

            if( SCIPvarGetUbGlobal(var) - 0.5 > lbs[v] )
            {
               SCIP_CALL( SCIPchgVarUbGlobal(scip, var, lbs[v]) );
            }
         }

         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         (*ndelconss)++;
         (*redundant) = TRUE;
      }
   }

   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);

   return SCIP_OKAY;
}

/** solve the cumulative sub problem */
static
SCIP_RETCODE solveSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< optcumulative constraint which collapsed to a cumulative constraint locally */
   SCIP_Bool             conflictanalysis,   /**< should conflict analysis be called for infeasible subproblems */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR**            binvars,            /**< array of variable representing if the job has to be processed on this machine */
   SCIP_VAR**            vars,               /**< start time variables of the activities which are assigned */
   int*                  durations,          /**< durations of the activities */
   int*                  demands,            /**< demands of the activities */
   int                   nvars,              /**< number of activities assigned to that machine */
   int*                  nfixedvars,         /**< pointer to store the numbver of fixed variables */
   int*                  nchgbds,            /**< pointer to store the number of changed bounds */
   int*                  ndelconss,          /**< pointer to store the number of deleted constraints */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is violated */
   )
{
   SCIP_Bool unbounded;
   SCIP_Bool solved;
   SCIP_Bool error;
   SCIP_Real* lbs;
   SCIP_Real* ubs;

   assert(scip != NULL);
   assert(!SCIPinProbing(scip));

   /* if we already tried solving this subproblem we do not do it again */
   if( consdata->triedsolving )
      return SCIP_OKAY;

   consdata->triedsolving = TRUE;

   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) );

   /* solve the cumulative condition separately */
   SCIP_CALL( solveCumulative(scip, nvars, vars, durations, demands, consdata->capacity, consdata->hmin, consdata->hmax, TRUE,
         lbs, ubs, 2000LL, &solved, cutoff, &unbounded, &error) );
   assert(!unbounded);

   if( !error )
   {
      if( *cutoff && conflictanalysis )
      {
         SCIP_Real* weights;
         SCIP_Bool infeasible;
         int v;

         SCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );

         for( v = 0; v < nvars; ++v )
         {
            int est;
            int lst;

            est = convertBoundToInt(scip, SCIPvarGetLbLocal(vars[v]));
            lst = convertBoundToInt(scip, SCIPvarGetUbLocal(vars[v]));

            if( demands[v] == 0.0 || durations[v] == 0.0 )
               return SCIP_ERROR;

            weights[v] = (lst - est) / (demands[v] * durations[v]); /*lint !e653*/
         }
         SCIPsortRealPtrPtrIntInt(weights, (void*)binvars, (void*)vars, durations, demands, nvars);

         SCIPfreeBufferArray(scip, &weights);

         while( nvars > 1 )
         {
            SCIP_CALL( solveCumulative(scip, nvars-1, vars, durations, demands, consdata->capacity, consdata->hmin, consdata->hmax, TRUE,
                  lbs, ubs, 2000LL, &solved, &infeasible, &unbounded, &error) );

            if( !infeasible )
               break;
            nvars--;
         }

         /**@todo try to shrink the initial explanation */

         SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

         for( v = 0; v < nvars; ++v )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, binvars[v]) );

            /* we have to add the lower and upper bounds of of the start time variable to have a valid reason */
            SCIP_CALL( SCIPaddConflictLb(scip, vars[v], NULL) );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[v], NULL) );
         }

         /* perform conflict analysis */
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }
      else
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_Bool allfixed;
         int v;

         allfixed = TRUE;

         for( v = 0; v < nvars; ++v )
         {
            /* check if variable is fixed */
            if( lbs[v] + 0.5 > ubs[v] )
            {
               SCIP_CALL( SCIPfixVar(scip, vars[v], lbs[v], &infeasible, &tightened) );
               assert(!infeasible);

               if( tightened )
               {
                  (*nfixedvars)++;
                  consdata->triedsolving = FALSE;
               }
            }
            else
            {
               SCIP_CALL( SCIPtightenVarLb(scip, vars[v], lbs[v], TRUE, &infeasible, &tightened) );
               assert(!infeasible);

               if( tightened )
               {
                  (*nchgbds)++;
                  consdata->triedsolving = FALSE;
               }

               SCIP_CALL( SCIPtightenVarUb(scip, vars[v], ubs[v], TRUE, &infeasible, &tightened) );
               assert(!infeasible);

               if( tightened )
               {
                  (*nchgbds)++;
                  consdata->triedsolving = FALSE;
               }

               allfixed = FALSE;
            }
         }

         /* if all variables are fixed, remove the optcumulative constraint since it is redundant */
         if( allfixed )
         {
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
            (*ndelconss)++;
         }
      }
   }

   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);

   return SCIP_OKAY;
}

/** check if the given constrait is valid; checks each starting point of a job whether the remaining capacity is at
 *  least zero or not. If not (*violated) is set to TRUE
 */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   SCIP_Bool*            violated,           /**< pointer to store if the constraint is violated */
   SCIP_Bool             printreason         /**< should the reason for the violation be printed? */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_VAR** vars;
   SCIP_Bool auxiliary;
   int* demands;
   int* durations;
   int nfixedones;
   int nfixedzeros;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("check optcumulative constraints <%s>\n", SCIPconsGetName(cons));

   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &durations, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, consdata->nvars) );

   /* collect information of all activities which are assigned to that machine in the given solution */
   collectSolActivities(scip, consdata, sol, binvars, vars, durations, demands, &nvars, &nfixedones, &nfixedzeros, &auxiliary);

   if( nvars > 0 )
   {
      /* check the cumulative condition */
      SCIP_CALL( SCIPcheckCumulativeCondition(scip, sol, nvars, vars,
            durations, demands, consdata->capacity, consdata->hmin, consdata->hmax, violated, cons, printreason) );
   }

   /* free all buffers */
   SCIPfreeBufferArray(scip, &demands);
   SCIPfreeBufferArray(scip, &durations);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &binvars);

   return SCIP_OKAY;
}

/** check if the given constrait is valid; checks each starting point of a job whether the remaining capacity is at
 *  least zero or not. If not (*violated) is set to TRUE
 */
static
SCIP_RETCODE enfopsCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   SCIP_SOL*             trysol,             /**< primal solution to construct, or NULL */
   SCIP_Bool*            violated,           /**< pointer to store if the constraint is violated/infeasible */
   SCIP_Bool*            consadded,          /**< pointer to store if a constraint was added */
   SCIP_Bool*            solfeasible         /**< pointer to store if the constraint solution is potentially feasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_VAR** vars;
   SCIP_Bool auxiliary;
   int* demands;
   int* durations;
   int nfixedones;
   int nfixedzeros;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("enforce optcumulative constraints <%s>\n", SCIPconsGetName(cons));

   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &durations, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, consdata->nvars) );

   /* collect information of all activities which are assigned to that machine in the given solution */
   collectSolActivities(scip, consdata, NULL, binvars, vars, durations, demands, &nvars, &nfixedones, &nfixedzeros, &auxiliary);

   (*violated) = FALSE;

   if( nvars > 0 )
   {
      /* check the cumulative condition */
      SCIP_CALL( SCIPcheckCumulativeCondition(scip, NULL, nvars, vars,
            durations, demands, consdata->capacity, consdata->hmin, consdata->hmax, violated, cons, FALSE) );

      if( *violated && auxiliary && !consdata->triedsolving )
      {
         SCIP_Real* lbs;
         SCIP_Real* ubs;
         SCIP_Bool infeasible;
         SCIP_Bool unbounded;
         SCIP_Bool error;
         SCIP_Bool solved;

         if( nfixedones == nvars )
            consdata->triedsolving = TRUE;

         SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) );

         /* solve the cumulative condition separately */
         SCIP_CALL( solveCumulative(scip, nvars, vars, durations, demands, consdata->capacity, consdata->hmin, consdata->hmax,
               FALSE, lbs, ubs, 1000LL, &solved, &infeasible, &unbounded, &error) );
         assert(!unbounded);

         if( !error )
         {
            if( infeasible )
            {

#ifdef SCIP_DISABLED_CODE
               SCIP_Real* weights;
               int v;

               SCIP_CALL( SCIPallocBufferArray(scip, &weights, nvars) );

               for( v = 0; v < nvars; ++v )
               {
                  int est;
                  int lst;

                  est = convertBoundToInt(scip, SCIPvarGetLbGlobal(vars[v]));
                  lst = convertBoundToInt(scip, SCIPvarGetUbGlobal(vars[v]));
                  weights[v] = (lst - est) / (consdata->demands[v] * consdata->durations[v]);
               }
               SCIPsortRealPtrPtrIntInt(weights, (void*)binvars, (void*)vars, durations, demands, nvars);

               SCIPfreeBufferArray(scip, &weights);

               while( nvars > 1 && !SCIPisStopped(scip) )
               {
                  SCIP_CALL( solveCumulative(scip, nvars-1, vars, durations, demands, consdata->capacity, consdata->hmin, consdata->hmax,
                        FALSE, lbs, ubs, 1000LL, &solved, &infeasible, &unbounded, &error) );

                  if( !infeasible )
                     break;

                  nvars--;
               }
#endif

               /* create and adds a conflict constraint (logicor constraint) */
               SCIP_CALL( createConflictCons(scip, SCIPconsGetName(cons), binvars,  nvars) );

               (*solfeasible) = FALSE;
               (*consadded) = TRUE;
            }
            else if( solved && *solfeasible && trysol != NULL )
            {
               int v;

               for(v = 0; v < nvars; ++v )
               {
                  SCIP_CALL( SCIPsetSolVal(scip, trysol, vars[v], lbs[v]) );
               }
            }
            else
               (*solfeasible) = FALSE;
         }

         SCIPfreeBufferArray(scip, &ubs);
         SCIPfreeBufferArray(scip, &lbs);
      }
   }

   /* free all buffers */
   SCIPfreeBufferArray(scip, &demands);
   SCIPfreeBufferArray(scip, &durations);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &binvars);

   return SCIP_OKAY;
}

#if 0
/** enforce the LP or pseudo solution */
static
SCIP_RETCODE enfoCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   SCIP_Bool*            violated,           /**< pointer to store if the constraint is violated */
   SCIP_Bool*            rowadded            /**< pointer to store if a row was added */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_VAR** vars;
   int* demands;
   int* durations;
   SCIP_Bool auxiliary;
   SCIP_Bool cutoff;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violated != NULL);

   SCIPdebugMessage("check optcumulative constraints <%s>\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &durations, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, consdata->nvars) );

   /* collect information of all activities which are assigned to that machine in the given solution */
   collectSolActivities(scip, consdata, NULL, binvars, vars, durations, demands, &nvars, &auxiliary);

   if( nvars > 0 )
   {
      /* check the cumulative condition */
      SCIP_CALL( SCIPcheckCumulativeCondition(scip, NULL, nvars, vars,
            durations, demands, consdata->capacity, consdata->hmin, consdata->hmax, violated, cons, FALSE) );

      if( *violated )
      {
#if 0
         /* create row */
         SCIP_CALL( createRow(scip, SCIPconsGetName(cons), binvars, vars, durations, demands, nvars,
               consdata->capacity, TRUE, &cutoff) );
#endif
         /* reset constraint age since it successfully detected infeasibility */
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      else
      {
         /* increase constraint age since it did not detected infeasibility */
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }
   }

   /* free all buffers */
   SCIPfreeBufferArray(scip, &demands);
   SCIPfreeBufferArray(scip, &durations);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &binvars);

   return SCIP_OKAY;
}
#endif

/** upgrade constraints to an cumulative constraint */
static
SCIP_RETCODE upgradeCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   int*                  ndelconss,          /**< pointer to store the number of deleted constraints */
   int*                  nupgdconss,         /**< pointer to store the number of upgrade constraints */
   SCIP_Bool*            mustpropagate       /**< pointer to store if the constraints has to be propagated */
   )
{
   SCIP_CONSDATA* consdata;
   int nvars;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* (debug) check if the counter of the constraint are correct */
   checkCounters(consdata);

   if( nvars == 0 && consdata->nfixedzeros == nvars )
   {
      SCIPdebugMessage("delete optcumulative constraint <%s> since it contains no jobs\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPdelCons(scip, cons) );
      (*ndelconss)++;
      (*mustpropagate) = FALSE;
   }
   else if( nvars == 1 )
   {
      SCIPdebugMessage("delete optcumulative constraint <%s> since it contains only one jobs\n", SCIPconsGetName(cons));

      if( consdata->capacity < consdata->demands[0] )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;

         SCIP_CALL( SCIPfixVar(scip, consdata->binvars[0], 0.0, &infeasible, &tightened) );
         assert(!infeasible);
         assert(tightened);
      }

      SCIP_CALL( SCIPdelCons(scip, cons) );
      (*ndelconss)++;
      (*mustpropagate) = FALSE;
   }
   else if( consdata->nglbfixedones == nvars )
   {
      SCIP_CONS* cumulativecons;
      char name[SCIP_MAXSTRLEN];

      SCIPdebugMessage("upgrade optcumulative constraint <%s> to cumulative constraint\n", SCIPconsGetName(cons));

      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_cumulative", SCIPconsGetName(cons));

      SCIP_CALL( SCIPcreateConsCumulative(scip, &cumulativecons, name, consdata->nvars, consdata->vars, consdata->durations, consdata->demands, consdata->capacity,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPsetHminCumulative(scip, cumulativecons, consdata->hmin) );
      SCIP_CALL( SCIPsetHmaxCumulative(scip, cumulativecons, consdata->hmax) );
      SCIP_CALL( SCIPaddCons(scip, cumulativecons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cumulativecons) );

      assert(!SCIPconsIsDeleted(cons));
      SCIP_CALL( SCIPdelCons(scip, cons) );

      (*nupgdconss)++;
      (*mustpropagate) = FALSE;
   }
   else if( consdata->nfixedones + consdata->nfixedzeros == nvars && consdata->nfixedones > 0 )
   {
      SCIP_CONS* cumulativecons;

      SCIP_VAR** binvars;
      SCIP_VAR** vars;
      int* durations;
      int* demands;
      int nfixedzeros;
      int nfixedones;

      SCIP_Bool auxiliary;

      char name[SCIP_MAXSTRLEN];

      SCIPdebugMessage("upgrade optcumulative constraint <%s> to cumulative constraint (locally)\n", SCIPconsGetName(cons));

      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_cumulative", SCIPconsGetName(cons));

      SCIP_CALL( SCIPallocBufferArray(scip, &vars, consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &binvars, consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &demands, consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &durations, consdata->nvars) );

      /* collect all activities which are locally assigned to that machine */
      collectActivities(consdata, binvars, vars, durations, demands, &nfixedones, &nfixedzeros, &auxiliary);

      SCIP_CALL( SCIPcreateConsCumulative(scip, &cumulativecons, name, nfixedones, vars, durations, demands, consdata->capacity,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            TRUE, SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPsetHminCumulative(scip, cumulativecons, consdata->hmin) );
      SCIP_CALL( SCIPsetHmaxCumulative(scip, cumulativecons, consdata->hmax) );
      SCIP_CALL( SCIPaddConsLocal(scip, cumulativecons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cumulativecons) );

      /* free all buffers */
      SCIPfreeBufferArray(scip, &durations);
      SCIPfreeBufferArray(scip, &demands);
      SCIPfreeBufferArray(scip, &binvars);
      SCIPfreeBufferArray(scip, &vars);

      assert(!SCIPconsIsDeleted(cons));
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );

      (*nupgdconss)++;
      (*mustpropagate) = FALSE;
   }
   else
      assert(consdata->nvars > 1);

   return SCIP_OKAY;
}

/** since the binary variable is fixed to zero, depending in the objective coefficient of the integer variable and the
 *  rounding locks, we might can fix the integer variable
 */
static
SCIP_RETCODE fixIntegerVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< integer variable to fix */
   SCIP_Bool             downlock,           /**< does the variable has down lock given by the optcumulative constraint */
   SCIP_Bool             uplock,             /**< does the variable has up lock given by the optcumulative constraint */
   int*                  nchgbds             /**< pointer to store the number changed variable bounds */
   )
{
   SCIP_Real objval;
   SCIP_Real fixvalue;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   objval = SCIPvarGetObj(var);
   fixvalue = SCIP_INVALID;

   /* if SCIP is in probing mode or during repropagation we cannot perform this dual reductions since this dual
    * reduction would end in an implication which can lead to cutoff the optimal solution
    */
   if( SCIPinProbing(scip) || SCIPinRepropagation(scip) )
      return SCIP_OKAY;

   assert(SCIPvarGetNLocksDown(var) >= (int)downlock);
   assert(SCIPvarGetNLocksUp(var) >= (int)uplock);

   if( SCIPisZero(scip, objval) )
   {
      /* the integer start time variable has a zero objective value; if only the optcumulative constraint
       * handler has a problem with rounding it down or up, then this issue is obsolete since binary
       * variable is fixed zero; therefore, rounding the integer down or up is a feasible dual reduction
       */
      if( SCIPvarGetNLocksDown(var) == (int)downlock )
         fixvalue = SCIPvarGetLbLocal(var);
      else if( SCIPvarGetNLocksUp(var) == (int)uplock )
         fixvalue = SCIPvarGetUbLocal(var);
      else
         return SCIP_OKAY;
   }
   else if( SCIPisNegative(scip, objval) && SCIPvarGetNLocksUp(var) == (int)uplock )
   {
      /* the integer start time variable has a negative objective value and only the optcumulative constraint
       * handler has a problem with rounding it up; since the binary variable is fixed the rounding up
       * issue is obsolete; there rounding it to the upper bound is the best thing we can do
       */
      fixvalue = SCIPvarGetUbLocal(var);
   }
   else if( SCIPisPositive(scip, objval) && SCIPvarGetNLocksDown(var) == (int)downlock )
   {
      /* the integer start time variable has a positive objective value and only the optcumulative
       * constraint handler has a problem with rounding it down; since the binary variable is fixed the
       * rounding down issue is obsolete; there rounding it to the lower bound is the best thing we can do
       */
      fixvalue = SCIPvarGetLbLocal(var);
   }
   else
      return SCIP_OKAY;

   /* the integer start time variable has a positive objective value and only the optcumulative
    * constraint handler has a problem with rounding it down; since the binary variable is fixed the
    * rounding down issue is obsolete; there rounding it to the lower bound is the best thing we can do
    */
   assert(fixvalue < SCIP_INVALID);
   SCIP_CALL( SCIPfixVar(scip, var, fixvalue, &infeasible, &tightened) );
   assert(!infeasible);

   if( tightened )
      (*nchgbds)++;

   return SCIP_OKAY;
}

/** deletes coefficient at given position from constraint data */
static
SCIP_RETCODE consdataDeletePos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< cumulative constraint data */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   assert(consdata != NULL);
   assert(pos < consdata->nvars);

   /* remove the rounding locks for the deleted variable */
   SCIP_CALL( unlockRounding(scip, cons, consdata->binvars[pos],
         consdata->vars[pos], consdata->downlocks[pos], consdata->uplocks[pos]) );

   consdata->downlocks[pos] = FALSE;
   consdata->uplocks[pos] = FALSE;

   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlrbinvars != NULL);
      assert(conshdlrdata->eventhdlrintvars != NULL);

      /* drop bound change events of variable */
      SCIP_CALL( dropEventBinvar(scip, cons, conshdlrdata->eventhdlrbinvars, pos) );
      SCIP_CALL( dropEventIntvar(scip, cons, conshdlrdata->eventhdlrintvars, pos) );
   }

   SCIPdebugMessage("remove variable <%s> from optcumulative constraint <%s>\n",
      SCIPvarGetName(consdata->binvars[pos]), SCIPconsGetName(cons));

   if( pos != consdata->nvars - 1 )
   {
      consdata->binvars[pos] = consdata->binvars[consdata->nvars-1];
      consdata->vars[pos] = consdata->vars[consdata->nvars-1];
      consdata->demands[pos] = consdata->demands[consdata->nvars-1];
      consdata->durations[pos] = consdata->durations[consdata->nvars-1];
      consdata->downlocks[pos] = consdata->downlocks[consdata->nvars-1];
      consdata->uplocks[pos] = consdata->uplocks[consdata->nvars-1];
   }

   consdata->nvars--;

   /* (debug) check if the counter of the constraint are correct */
   checkCounters(consdata);

   consdata->relaxadded = FALSE;
   consdata->normalized = FALSE;

   return SCIP_OKAY;
}

/** remove all jobs for which the binary variable is globally fixed to zero */
static
SCIP_RETCODE applyZeroFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   int*                  nchgcoefs,          /**< pointer to store the number changed coefficients */
   int*                  nchgbds             /**< pointer to store the number changed variable bounds */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( v = consdata->nvars-1; v >= 0 && consdata->nglbfixedzeros > 0; --v )
   {
      assert(consdata->binvars[v] != NULL);
      if( SCIPvarGetUbGlobal(consdata->binvars[v]) < 0.5 )
      {
         SCIPdebugMessage("variable <%s> is globally fixed to zero\n", SCIPvarGetName(consdata->binvars[v]));

         /* fix integer start time variable if possible */
         if( SCIPconsIsChecked(cons) )
         {
            SCIP_CALL( fixIntegerVariable(scip, consdata->vars[v], consdata->downlocks[v], consdata->uplocks[v], nchgbds) );
         }

         /* remove the job */
         SCIP_CALL( consdataDeletePos(scip, consdata, cons, v) );
         (*nchgcoefs)++;

         /* mark constraint to be checked for redundancy */
         consdata->triedredundant = TRUE;
      }
   }

   /* (debug) check if the counter of the constraint are correct */
   checkCounters(consdata);

   /* check that all variables fixed to zero are removed */
   assert(consdata->nglbfixedzeros == 0);

   return SCIP_OKAY;
}

/** remove jobs which have a duration or demand of zero (zero energy) or lay outside the efficient horizon [hmin, hmax);
 *  this is done in the SCIP_DECL_CONSINITPRE() callback
 */
static
SCIP_RETCODE removeIrrelevantJobs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to propagate */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int demand;
   int duration;
   int hmin;
   int hmax;
   int est;
   int lct;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   hmin = consdata->hmin;
   hmax = consdata->hmax;

   SCIPdebugMessage("check for irrelevant jobs within cumulative constraint <%s>[%d,%d)\n",
      SCIPconsGetName(cons), hmin, hmax);

   for( j = consdata->nvars-1; j >= 0; --j )
   {
      var = consdata->vars[j];
      demand = consdata->demands[j];
      duration = consdata->durations[j];

      /* earliest completion time (ect) and latest start time (lst) */
      est = convertBoundToInt(scip, SCIPvarGetLbGlobal(var));
      lct = convertBoundToInt(scip, SCIPvarGetUbGlobal(var)) + duration;

      if( demand == 0 || duration == 0 )
      {
         /* jobs with zero demand or zero duration can be removed */
         SCIPdebugMessage("  remove variable <%s> due to zero %s\n",
            SCIPvarGetName(var), demand == 0 ? "demand" : "duration");

         /* remove variable form constraint */
         SCIP_CALL( consdataDeletePos(scip, consdata, cons, j) );
      }
      else if( est >= hmax || lct <= hmin )
      {
         SCIPdebugMessage("  remove variable <%s>[%d,%d] with duration <%d>\n",
            SCIPvarGetName(var), est, lct - duration, duration);

         /* delete variable at the given position */
         SCIP_CALL( consdataDeletePos(scip, consdata, cons, j) );
      }
   }

   return SCIP_OKAY;
}

/** presolve cumulative condition w.r.t. effective horizon by detecting irrelevant variables */
static
SCIP_RETCODE presolveCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   int*                  nchgcoefs,          /**< pointer to store the number of changed coefficients */
   int*                  nchgsides,          /**< pointer to store the number of changed sides */
   SCIP_Bool*            cutoff              /**< buffer to store whether a cutoff is detected */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool* irrelevants;
   int nvars;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   assert(nvars > 1);

   SCIP_CALL( SCIPallocBufferArray(scip, &irrelevants, nvars) );
   BMSclearMemoryArray(irrelevants, nvars);

   /* use presolving of cumulative constraint handler to process cumulative condition */
   SCIP_CALL( SCIPpresolveCumulativeCondition(scip, nvars, consdata->vars, consdata->durations,
         consdata->hmin, consdata->hmax, consdata->downlocks, consdata->uplocks, cons,
         irrelevants, nfixedvars, nchgsides, cutoff) );

   /* remove all variable which are irrelevant; note we have to iterate backwards do to the functionality of of
    * consdataDeletePos()
    */
   for( v = nvars-1; v >= 0; --v )
   {
      SCIP_VAR* var;
      int ect;
      int lst;

      if( !irrelevants[v] )
         continue;

      var = consdata->vars[v];
      assert(var != NULL);

      ect = convertBoundToInt(scip, SCIPvarGetLbGlobal(var)) + consdata->durations[v];
      lst = convertBoundToInt(scip, SCIPvarGetUbGlobal(var));

      /* check if the jobs runs completely during the effective horizon */
      if( lst <= consdata->hmin && ect >= consdata->hmax )
      {
         assert(!consdata->downlocks[v]);
         assert(!consdata->uplocks[v]);

         if( consdata->capacity < consdata->demands[v] )
         {
            SCIP_Bool infeasible;
            SCIP_Bool tightened;

            SCIP_CALL( SCIPfixVar(scip, consdata->binvars[0], 0.0, &infeasible, &tightened) );
            assert(!infeasible);
            assert(tightened);
            (*nfixedvars)++;

            consdata->capacity -= consdata->demands[v];

            SCIP_CALL( consdataDeletePos(scip, consdata, cons, v) );
            (*nchgcoefs)++;
         }
      }
      else
      {
         SCIP_CALL( consdataDeletePos(scip, consdata, cons, v) );
         (*nchgcoefs)++;
      }
   }

   SCIPdebugMessage("constraint <%s>[%d,%d) <= %d has %d variables left\n", SCIPconsGetName(cons),
      consdata->hmin, consdata->hmax, consdata->capacity, nvars);

   SCIPfreeBufferArray(scip, &irrelevants);

   return SCIP_OKAY;
}

/** create an an set partitioning constraint */
static
SCIP_RETCODE createSetPackingCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_VAR*             var2                /**< second variable */
   )
{
   SCIP_CONS* cons;

   SCIP_CALL( SCIPcreateConsBasicSetpack(scip, &cons, "implication", 0, NULL) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( SCIPaddCoefSetppc(scip, cons, var1) );
   SCIP_CALL( SCIPaddCoefSetppc(scip, cons, var2) );
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/** create variable bound constraint */
static
SCIP_RETCODE createVarboundCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             binvar,             /**< binary variable x */
   SCIP_VAR*             intvar,             /**< integer variable y */
   int                   bound,              /**< variable bound */
   SCIP_Bool             lower               /**< variable lower bound? (Otherwise upper bound) */
   )
{
   SCIP_CONS* cons;
   SCIP_Real coef;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(scip != NULL);

   if( lower )
   {
      lhs = SCIPvarGetLbGlobal(intvar);
      rhs = SCIPinfinity(scip);
      coef = lhs - bound;
   }
   else
   {
      lhs = -SCIPinfinity(scip);
      rhs = SCIPvarGetUbGlobal(intvar);
      coef = rhs - bound;
   }

   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, "implication", intvar, binvar, coef, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/** create bound disjunction constraint */
static
SCIP_RETCODE createBounddisjunctionCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             binvar,             /**< binary variable x */
   SCIP_VAR*             intvar,             /**< integer variable y */
   int                   lb,                 /**< lower bound */
   int                   ub                  /**< lower bound */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, 3) );

   /* intvar >= ub */
   vars[0] = intvar;
   boundtypes[0] = SCIP_BOUNDTYPE_LOWER;
   bounds[0] = ub;

   /* intvar <= lb */
   vars[1] = intvar;
   boundtypes[1] = SCIP_BOUNDTYPE_UPPER;
   bounds[1] = lb;

   /* binvar <= 0.0 */
   vars[2] = binvar;
   boundtypes[2] = SCIP_BOUNDTYPE_LOWER;
   bounds[2] = 0.0;

   SCIP_CALL( SCIPcreateConsBasicBounddisjunction(scip, &cons, "implication", 3, vars, boundtypes, bounds) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &boundtypes);
   SCIPfreeBufferArray(scip, &bounds);

   return SCIP_OKAY;
}

/** detect implication */
static
SCIP_RETCODE detectImplications(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< optcumulative constraint */
   int*                  nchgcoefs,          /**< pointer to store the number of changed coefficients */
   int*                  naddconss           /**< pointer to store the number of added constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_VAR** vars;
   int* durations;
   int hmin;
   int hmax;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   binvars = consdata->binvars;
   durations = consdata->durations;

   hmin = consdata->hmin;
   hmax = consdata->hmax;
   assert(hmin < hmax);

   SCIPdebugMessage("search for implications <%s>[%d,%d) <= %d\n", SCIPconsGetName(cons), hmin, hmax, consdata->capacity);

   /* we loop backwards since we are deleting variable out of the constraint */
   for( v = consdata->nvars-1; v >= 0; --v )
   {
      SCIP_VAR* var;
      int start;
      int end;

      var = vars[v];
      assert(var != NULL);

      /* skip start time variables which are not globally fixed */
      if( SCIPvarGetLbGlobal(var) + 0.5 < SCIPvarGetUbGlobal(var) )
         continue;

      /* adjust the code for resources with capacity larger than one ??????????????? */
      if( consdata->demands[v] < consdata->capacity )
         continue;

      start = convertBoundToInt(scip, SCIPvarGetLbGlobal(var));
      assert(start < hmax);

      end = start + durations[v];
      assert(end > hmin);

      SCIPdebugMessage("candidate <%s> (start %d, end %d, demand %d)\n", SCIPvarGetName(var), start, end, consdata->demands[v]);

      if( start <= hmin && end >= hmax )
      {
         int j;

         /* job runs during the complete time horizon */
         for( j = 0; j < consdata->nvars; ++j )
         {
            SCIP_VAR* implvar;
            int est;
            int ect;
            int lst;

            if( j == v )
               continue;

            implvar = vars[j];
            assert(implvar != NULL);

            est = convertBoundToInt(scip, SCIPvarGetLbGlobal(implvar));
            ect = est + durations[j];
            lst = convertBoundToInt(scip, SCIPvarGetUbGlobal(implvar));

            SCIPdebugMessage("variable <%s>[%d,%d] (duration %d, demand %d)\n", SCIPvarGetName(implvar), est, lst, durations[j], consdata->demands[j]);

            /* check if the job will overlap with effective horizon, hence, only one of the two jobs can be scheduled on
             * that machine
             */
            if( ect > hmin && lst < hmax )
            {
               SCIP_CALL( createSetPackingCons(scip, binvars[v], binvars[j]) );
               (*naddconss)++;
            }
            else if( lst < hmax )
            {
               SCIP_CALL( createVarboundCons(scip, binvars[v], implvar, hmin - durations[j], FALSE) );
               (*naddconss)++;
            }
            else if( ect > hmin )
            {
               SCIP_CALL( createVarboundCons(scip, binvars[v], implvar, hmax, TRUE) );
               (*naddconss)++;
            }
            else
            {
               SCIP_CALL( createBounddisjunctionCons(scip, binvars[v], implvar, hmin - durations[j], hmax) );
               (*naddconss)++;
            }
         }
      }
      else if( start <= hmin )
      {
         int j;

         assert(end > hmin);

         /* job overlaps with hmin */
         for( j = 0; j < consdata->nvars; ++j )
         {
            SCIP_VAR* implvar;
            int est;
            int ect;
            int lst;

            if( j == v )
               continue;

            implvar = vars[j];
            assert(implvar != NULL);

            est = convertBoundToInt(scip, SCIPvarGetLbGlobal(implvar));
            ect = est + durations[j];
            lst = convertBoundToInt(scip, SCIPvarGetUbGlobal(implvar));

            SCIPdebugMessage("variable <%s>[%d,%d] (duration %d, demand %d)\n", SCIPvarGetName(implvar), est, lst, durations[j], consdata->demands[j]);

            if( lst < ect && hmin < ect && lst < end )
            {
               /* job j has a core which overlaps with job v within the effective horizon, hence, both jobs cannot run
                * at same time on that machine
                */
               SCIP_CALL( createSetPackingCons(scip, binvars[v], binvars[j]) );
               (*naddconss)++;
            }
            else if( end > lst )
            {
               SCIP_CALL( createSetPackingCons(scip, binvars[v], binvars[j]) );
               (*naddconss)++;
            }
            else if( est < end )
            {
               SCIP_CALL( createVarboundCons(scip, binvars[v], implvar, end, TRUE) );
               (*naddconss)++;
            }
         }
      }
      else if( end >= hmax )
      {
         int j;

         assert(start < hmax);

         /* job overlaps with hmax; that means if the job is scheduled on that machine all other jobs have to finish
          * before that job starts
          */
         for( j = 0; j < consdata->nvars; ++j )
         {
            SCIP_VAR* implvar;
            int ect;
            int lst;
            int lct;

            if( j == v )
               continue;

            implvar = vars[j];
            assert(implvar != NULL);

            ect = convertBoundToInt(scip, SCIPvarGetLbGlobal(implvar)) + durations[j];
            lst = convertBoundToInt(scip, SCIPvarGetUbGlobal(implvar));
            lct = lst + durations[j];

            SCIPdebugMessage("variable <%s>[%d,%d] (duration %d, demand %d)\n", SCIPvarGetName(implvar), ect - durations[j], lst, durations[j], consdata->demands[j]);

            if( lst < ect && start < ect && lst < hmax )
            {
               /* job j has a core which overlaps with job v within the effective horizon, hence, both jobs cannot run
                * at same time on that machine
                */
               SCIP_CALL( createSetPackingCons(scip, binvars[v], binvars[j]) );
               (*naddconss)++;
            }
            else if( start < ect )
            {
               SCIP_CALL( createSetPackingCons(scip, binvars[v], binvars[j]) );
               (*naddconss)++;
            }
            else if( lct > start )
            {
               /* job j potentially finishes to late, hence, if job v runs on that machine we can bound the start time
                * variable of job j form above
                */
               SCIP_CALL( createVarboundCons(scip, binvars[v], implvar, start - durations[j], FALSE) );
               (*naddconss)++;
            }
         }
      }
      else
         continue;

      SCIP_CALL( consdataDeletePos(scip, consdata, cons, v) );
      (*nchgcoefs)++;
   }

   return SCIP_OKAY;
}

/** propgates given constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   SCIP_Bool             conflictanalysis,   /**< should conflict analysis be called for infeasible subproblems */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   int*                  nchgbds,            /**< pointer to store the number changed variable bounds */
   int*                  ndelconss,          /**< pointer to store the number of deleted constraints */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff (infeasibility) was detected */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   SCIP_VAR** vars;
   SCIP_Bool auxiliary;
   int* durations;
   int* demands;
   int nfixedones;
   int nfixedzeros;
   int v;

   assert(cutoff != NULL);
   assert(*cutoff == FALSE);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 1);

   /* (debug) check if the counter of the constraint are correct */
   checkCounters(consdata);

   if( consdata->propagated && (consdata->nfixedones + consdata->nfixedzeros < consdata->nvars || consdata->triedsolving) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &durations, consdata->nvars) );

   /* collect all activities which are locally assigned to that machine */
   collectActivities(consdata, binvars, vars, durations, demands, &nfixedones, &nfixedzeros, &auxiliary);

   /* if more than one variable is assigned to that machine propagate the cumulative condition */
   if( !consdata->propagated && nfixedones > 1 )
   {
      SCIP_Bool* explanation;
      SCIP_Bool initialized;

      initialized = FALSE;

      SCIP_CALL( SCIPallocBufferArray(scip, &explanation, nfixedones) );
      BMSclearMemoryArray(explanation, nfixedones);

      /* propagate cumulative condition */
      SCIP_CALL( SCIPpropCumulativeCondition(scip, SCIP_PRESOLTIMING_ALWAYS, nfixedones, vars,
            durations, demands, consdata->capacity, consdata->hmin, consdata->hmax, cons, nchgbds, &initialized, explanation, cutoff) );

      /* in case of a conflict we have to extend the initial reason before the conflict analysis starts */
      if( initialized  && conflictanalysis )
      {
         assert(*cutoff == TRUE);

         for( v = 0; v < nfixedones; ++v )
         {
            if( explanation[v] )
            {
               SCIP_CALL( SCIPaddConflictBinvar(scip, binvars[v]) );
            }
         }

         /* perform conflict analysis */
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }

      SCIPfreeBufferArray(scip, &explanation);
   }
   assert(consdata->nvars > 1);

   /* if we are still feasible we can try to perform dual reductions; Note that we have to avoid dual reductions during
    * probing since these dual reductions can lead to wrong implications; the same hold in case of repropagating
    */
   if( !(*cutoff) && !SCIPinProbing(scip) && !SCIPinRepropagation(scip) )
   {
      if( nfixedzeros + nfixedones == consdata->nvars )
      {
         /* all binary variables are fixed */

         if( auxiliary )
         {
            /* we have an independent subproblems since all binary variables are fixed and the integer start time
             * variables belonging to the binary variables which are fixed to one are only locked by this constraint
             */
            SCIP_CALL( solveSubproblem(scip, cons, conflictanalysis, consdata, binvars, vars, durations, demands,
                  nfixedones, nfixedvars, nchgbds, ndelconss, cutoff) );
         }
      }
      else if( !consdata->propagated && nfixedones < consdata->nvars )
      {
         SCIP_PROFILE* profile;
         int hmin;
         int est;
         int lct;
         int pos;

         /* create empty resource profile with infinity resource capacity */
         SCIP_CALL( SCIPprofileCreate(&profile, INT_MAX) );

         /* create worst case resource profile */
         SCIP_CALL( SCIPcreateWorstCaseProfile(scip, profile, nfixedones, vars, durations, demands) );

         hmin = SCIPcomputeHmin(scip, profile, consdata->capacity);

         if( hmin < INT_MAX )
         {
            /* check if the not selected variables can be discard from the machine */
            for( v = 0; v < consdata->nvars && !(*cutoff) && !SCIPisStopped(scip) ; ++v )
            {
               SCIP_VAR* binvar;
               SCIP_VAR* var;

               binvar = consdata->binvars[v];
               assert(binvar != NULL);

               var = consdata->vars[v];
               assert(var != NULL);

               /* check if the binary choice variable is not fixed yet */
               if( SCIPvarGetLbLocal(binvar) + 0.5 < SCIPvarGetUbLocal(binvar) )
               {
                  SCIP_Real lb;
                  SCIP_Real ub;
                  SCIP_Bool infeasible;

                  assert(SCIPvarGetLbLocal(binvar) < 0.5);
                  assert(SCIPvarGetUbLocal(binvar) > 0.5);

                  est = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
                  lct = convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + consdata->durations[v];

                  SCIP_CALL( SCIPprofileInsertCore(profile, est, lct, consdata->demands[v], &pos, &infeasible) );
                  assert(!infeasible);
                  assert(pos == -1);

                  hmin = SCIPcomputeHmin(scip, profile, consdata->capacity);

                  SCIP_CALL( SCIPprofileDeleteCore(profile, est, lct, consdata->demands[v]) );

                  if( hmin == INT_MAX )
                     continue;

                  /* start probing mode */
                  SCIPdebugMessage("start probing\n");
                  SCIP_CALL( SCIPstartProbing(scip) );

                  SCIP_CALL( SCIPnewProbingNode(scip) );

                  SCIPdebugMessage("  fix variables <%s>[%g,%g] to 1.0\n",
                     SCIPvarGetName(binvar), SCIPvarGetLbLocal(binvar), SCIPvarGetUbLocal(binvar));

                  SCIP_CALL( SCIPfixVarProbing(scip, binvar, 1.0) );

                  SCIPdebugMessage("  run propagation\n");
                  SCIP_CALL( SCIPpropagateProbing(scip, 0, &infeasible, NULL) );

                  lb = SCIPvarGetLbLocal(var);
                  ub = SCIPvarGetUbLocal(var);

                  /* end probing mode */
                  SCIP_CALL( SCIPendProbing(scip) );
                  SCIPdebugMessage("end probing\n");

                  if( infeasible )
                  {
                     SCIP_Bool tightened;

                     /* propagation detected infeasibility, therefore, job cannot be processed by that machine */
                     SCIPdebugMessage("  probing detect infeasibility\n");
                     SCIPdebugMessage("  fix variable <%s> to 0.0\n", SCIPvarGetName(binvar));

                     /* since this bound change is dual reduction we have to avoid that this bound change is analyzed
                      * during the conflict analysis; otherwise all optimal solution might be removed: therefore, we
                      * SCIPtightenVarUb instead of SCIPinferBinvarCons()
                      */
                     SCIP_CALL( SCIPtightenVarUb(scip, binvar, 0.0, FALSE, &infeasible, &tightened) );
                     if( infeasible )
                        (*cutoff) = TRUE;
                     else if( tightened )
                     {
                        (*nchgbds)++;

                        /* fix integer start time variable if possible (before calling that method we have to leave the
                         * probing mode)
                         */
                        if( SCIPconsIsChecked(cons) )
                        {
                           SCIP_CALL( fixIntegerVariable(scip, var, consdata->downlocks[v], consdata->uplocks[v], nchgbds) );
                        }
                     }
                  }
                  else
                  {
                     SCIP_Bool tightened;

                     /* probing was feasible, therefore, we can adjust the bounds of the start time variable for that job */
                     SCIPdebugMessage("  probing stayed feasible\n");

                     assert(SCIPvarGetNLocksUp(var) >= (int)consdata->uplocks[v]);
                     if( SCIPvarGetNLocksUp(var) == (int)consdata->uplocks[v] )
                     {
                        SCIPdebugMessage("  variable <%s> change lower bound from <%g> to <%g>\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), lb);

                        /* for this bound change there is no inference information needed since no other constraint can
                         * use this bound change to reason something
                         */
                        SCIP_CALL( SCIPtightenVarLb(scip, var, lb, FALSE, &infeasible, &tightened) );
                        assert(!infeasible);

                        if( tightened )
                           (*nchgbds)++;
                     }

                     assert(SCIPvarGetNLocksDown(var) >= (int)consdata->downlocks[v]);
                     if( SCIPvarGetNLocksDown(var) == (int)consdata->downlocks[v] )
                     {
                        SCIPdebugMessage("  variable <%s> change upper bound from <%g> to <%g>\n", SCIPvarGetName(var), SCIPvarGetUbLocal(var), ub);

                        /* for this boound change there is no inference information needed since no other constraint can
                         * use this bound change to reason something
                         */
                        SCIP_CALL( SCIPtightenVarUb(scip, var, ub, FALSE, &infeasible, &tightened) );
                        assert(!infeasible);

                        if( tightened )
                           (*nchgbds)++;
                     }
                  }
               }
               else if( SCIPvarGetUbLocal(binvar) < 0.5 && SCIPconsIsChecked(cons) )
               {
                  /* if the binary choice variable is fixed to zero we can try to perform a dual reductions */
                  SCIP_CALL( fixIntegerVariable(scip, var, consdata->downlocks[v], consdata->uplocks[v], nchgbds) );
               }
            }
         }

         /* free worst case profile */
         SCIPprofileFree(&profile);
      }
   }

   /* mark constraint to be propagated */
   if( !SCIPinProbing(scip) )
      consdata->propagated = TRUE;

   /* free all buffers */
   SCIPfreeBufferArray(scip, &durations);
   SCIPfreeBufferArray(scip, &demands);
   SCIPfreeBufferArray(scip, &binvars);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyOptcumulative)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrOptcumulative(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitOptcumulative NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitOptcumulative NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      /* remove jobs which have a duration or demand of zero (zero energy) or lay outside the effective horizon [hmin,
       * hmax)
       */
      SCIP_CALL( removeIrrelevantJobs(scip, conss[c]) );
   }

   /* find trysol heuristic */
   if( conshdlrdata->heurtrysol == NULL )
   {
      conshdlrdata->heurtrysol = SCIPfindHeur(scip, "trysol");
   }

   return SCIP_OKAY;
}

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreOptcumulative NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolOptcumulative NULL

/** constraint enforcing method of constraint handler for relaxation solutions */
#define consEnforelaxOptcomulative NULL

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolOptcumulative)
{  /*lint --e{715}*/
   int c;

   assert(scip != NULL);

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

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL );
   assert(*consdata != NULL );

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlrbinvars != NULL);
   assert(conshdlrdata->eventhdlrintvars != NULL);

   /* if constraint belongs to transformed problem space, drop bound change events on variables */
   if( (*consdata)->nvars > 0 && SCIPvarIsTransformed((*consdata)->vars[0]) )
   {
      SCIP_CALL( dropAllEvents(scip, cons, conshdlrdata->eventhdlrbinvars, conshdlrdata->eventhdlrintvars) );
   }

   /* free optcumulative constraint data */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlrbinvars != NULL);
   assert(conshdlrdata->eventhdlrintvars != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);

   SCIPdebugMessage("transform optcumulative constraint <%s>\n", SCIPconsGetName(sourcecons));

   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->nvars, sourcedata->vars, sourcedata->binvars,
         sourcedata->durations, sourcedata->demands, sourcedata->capacity, SCIPconsIsChecked(sourcecons)) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   assert(targetdata->nglbfixedones == 0);
   assert(targetdata->nglbfixedzeros == 0);
   assert(targetdata->nfixedones == 0);
   assert(targetdata->nfixedzeros == 0);

   /* catch bound change events of variables */
   SCIP_CALL( catchAllEvents(scip, *targetcons, conshdlrdata->eventhdlrbinvars, conshdlrdata->eventhdlrintvars) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool rowadded;
   SCIP_Bool consadded;
   SCIP_Bool cutoff;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   rowadded = FALSE;
   consadded = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(SCIPconsIsInitial(conss[c]));
      SCIP_CALL( addRelaxation(scip, conshdlr, conshdlrdata, conss[c], &rowadded, &consadded, &cutoff) );
      /* ignore cutoff value */
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpOptcumulative)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool rowadded;
   SCIP_Bool consadded;
   SCIP_Bool cutoff;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   rowadded = FALSE;
   consadded = FALSE;
   cutoff = FALSE;

   for( c = 0; c < nconss && ! cutoff; ++c )
   {
      SCIP_CALL( addRelaxation(scip, conshdlr, conshdlrdata, conss[c], &rowadded, &consadded, &cutoff) );
   }

   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if( consadded )
      *result = SCIP_CONSADDED;
   else if( rowadded )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}/*lint !e715*/


/** separation method of constraint handler for arbitrary primal solutions */
#define consSepasolOptcumulative NULL


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_SOL* trysol;
   SCIP_Bool violated;
   SCIP_Bool consviolated;
   SCIP_Bool consadded;
   SCIP_Bool solfeasible;
   int c;

   SCIPdebugMessage("method: enforce LP solution (nconss %d)\n", nconss);

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   violated = FALSE;
   consviolated = FALSE;
   consadded = FALSE;
   solfeasible = TRUE;
   trysol = NULL;

   /* create pseudo solution */
   if( conshdlrdata->heurtrysol != NULL )
   {
      SCIP_CALL( SCIPcreateCurrentSol(scip, &trysol, NULL) );
   }

   /* check all constraints even if one is dectected be violated */
   for( c = 0; c < nconss && (!violated || solfeasible); ++c )
   {
      SCIP_CALL( enfopsCons(scip, conss[c], trysol, &consviolated, &consadded, &solfeasible) );
      violated = violated || consviolated;
   }

   /* add a potentially feasible solution was constructed we pass it to the heuristic try sol */
   if( solfeasible && violated && trysol != NULL )
   {
#ifdef SCIP_DEBUG
      FILE* file;
      file = fopen("build.sol", "w");

      if( file != NULL )
      {
         SCIP_CALL( SCIPprintSol(scip, trysol, file, FALSE) );
         fclose(file);
      }
#endif

      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->heurtrysol, trysol) );
   }

   SCIP_CALL( SCIPfreeSol(scip, &trysol) );

   if( consadded )
      *result = SCIP_CONSADDED;
   else if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_SOL* trysol;
   SCIP_Bool violated;
   SCIP_Bool consadded;
   SCIP_Bool solfeasible;
   int c;

   SCIPdebugMessage("method: enforce pseudo solution\n");

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   violated = FALSE;
   consadded = FALSE;
   solfeasible = TRUE;
   trysol = NULL;

   /* create pseudo solution */
   if( conshdlrdata->heurtrysol != NULL )
   {
      SCIP_CALL( SCIPcreateCurrentSol(scip, &trysol, NULL) );
   }

   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( enfopsCons(scip, conss[c], trysol, &violated, &consadded, &solfeasible) );
   }

   /* add a potentially feasible solution was constructed we pass it to the heuristic try sol */
   if( solfeasible && violated && trysol != NULL )
   {
      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->heurtrysol, trysol) );
   }

   SCIP_CALL( SCIPfreeSol(scip, &trysol) );

   if( consadded )
      *result = SCIP_CONSADDED;
   else if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckOptcumulative)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   violated = FALSE;

   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], sol, &violated, printreason) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_Bool cutoff;
   int nfixedvars;
   int nupgdconss;
   int ndelconss;
   int nchgcoefs;
   int nchgbds;
   int c;

   assert(scip != NULL);
   assert(nconss > 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   nfixedvars = 0;
   nupgdconss = 0;
   ndelconss = 0;
   nchgcoefs = 0;
   nchgbds = 0;
   cutoff = FALSE;

   SCIPdebugMessage("propagate %d optcumulative constraints (probing: %u)\n", nusefulconss, SCIPinProbing(scip));

   /* first propagate only the useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_Bool mustpropagate;
      int oldnchgcoefs;
      int oldnchgbds;

      cons = conss[c];
      mustpropagate = TRUE;
      oldnchgcoefs = nchgcoefs;
      oldnchgbds = nchgbds;

      /* it might be that the constraint is already deleted which can be case if SCIP is in probing mode */
      if( SCIPconsIsDeleted(cons) )
      {
         assert(SCIPinProbing(scip));
         continue;
      }

      /* try to upgrade optcumulative to cumulative constraint which is possible if all remaining binary variables are
       * fixed to one; in case the constraint has no variable left it is removed
       */
      if( !SCIPinProbing(scip) )
      {
         SCIP_Bool redundant;

         /* remove all jobs for which the binary variable is globally fixed to zero */
         SCIP_CALL( applyZeroFixings(scip, cons, &nchgcoefs, &nchgbds) );

         SCIP_CALL( checkRedundancy(scip, cons, &ndelconss, &redundant) );

         if( redundant )
            continue;

         SCIP_CALL( upgradeCons(scip, cons, &ndelconss, &nupgdconss, &mustpropagate) );
      }

      if( mustpropagate )
      {
         SCIP_CALL( propagateCons(scip, cons, conshdlrdata->conflictanalysis, &nfixedvars, &nchgbds, &ndelconss, &cutoff) );
      }

      /* update the age of the constraint w.r.t. success of the propagation rule */
      if( oldnchgbds < nchgbds || oldnchgcoefs < nchgcoefs )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      else
      {
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }
   }

   if( cutoff )
   {
      SCIPdebugMessage("propagation detected a cutoff\n");
      *result = SCIP_CUTOFF;
   }
   else if( nfixedvars > 0 || nchgbds > 0 || nupgdconss > 0 )
   {
      SCIPdebugMessage("propagation detected %d bound changes\n", nchgbds);
      *result = SCIP_REDUCEDDOM;
   }
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_Bool cutoff;
   SCIP_Bool mustpropagate;
   int oldnchgbds;
   int oldndelconss;
   int oldnupgdconss;
   int oldnfixedvars;
   int c;

   assert(scip != NULL);
   assert(nconss > 0);
   assert(!SCIPinProbing(scip));

   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;
   oldnfixedvars = *nfixedvars;
   cutoff = FALSE;

   SCIPdebugMessage("presolve %d optcumulative constraints\n", nconss);

   for( c = 0; c < nconss && !cutoff; ++c )
   {
      SCIP_CONSDATA* consdata;

      cons = conss[c];
      mustpropagate = TRUE;

      /* remove all jobs for which the binary variable is globally fixed to zero */
      SCIP_CALL( applyZeroFixings(scip, cons, nchgcoefs, nchgbds) );

      /* try to upgrade optcumulative to cumulative constraint which is possible if all remaining binary variables are
       * fixed to one; in case the constraint has no or one variable left it is removed
       */
      SCIP_CALL( upgradeCons(scip, cons, ndelconss, nupgdconss, &mustpropagate) );

      if( mustpropagate )
      {
         int nvars;
         int hmin;
         int hmax;
         int split;

         consdata = SCIPconsGetData(cons);
         assert(consdata != NULL);

         nvars = consdata->nvars;
         assert(nvars > 1);

         if( !consdata->normalized )
         {
            /* divide demands and capacity by their greatest common divisor */
            SCIP_CALL( SCIPnormalizeCumulativeCondition(scip, nvars, consdata->vars, consdata->durations,
                  consdata->demands, &consdata->capacity, nchgcoefs, nchgsides) );
            consdata->normalized = TRUE;
         }

         /* propagate the constaint */
         SCIP_CALL( propagateCons(scip, cons,  FALSE, nfixedvars, nchgbds, ndelconss, &cutoff) );

         /* if a cutoff was detected we are done */
         if( cutoff )
            break;

         /* check if the optimal cumulative constraint can be decomposed */
         SCIP_CALL( SCIPsplitCumulativeCondition(scip, nvars, consdata->vars, consdata->durations,
               consdata->demands, consdata->capacity, &hmin, &hmax, &split) );

         /* check if this time point improves the effective horizon */
         if( consdata->hmin < hmin )
         {
            SCIPdebugMessage("cumulative constraint <%s> adjust hmin <%d> -> <%d>\n", SCIPconsGetName(cons), consdata->hmin, hmin);

            consdata->hmin = hmin;
            (*nchgsides)++;
         }

         /* check if this time point improves the effective horizon */
         if( consdata->hmax > hmax )
         {
            SCIPdebugMessage("cumulative constraint <%s> adjust hmax <%d> -> <%d>\n", SCIPconsGetName(cons), consdata->hmax,  hmax);
            consdata->hmax = hmax;
            (*nchgsides)++;
         }

         /* check if the constraint is redundant */
         if( consdata->hmax <= consdata->hmin )
         {
            SCIPdebugMessage("constraint <%s> is redundant since hmax(%d) <= hmin(%d)\n",
               SCIPconsGetName(cons), consdata->hmax, consdata->hmin);

            SCIP_CALL( SCIPdelCons(scip, cons) );
            (*ndelconss)++;

            continue;
         }

         /* check if the cumulative constraint can be decomposed */
         if( consdata->hmin < split && split < consdata->hmax )
         {
            SCIP_CONS* splitcons;
            SCIP_CONSDATA* splitconsdata;
            char name[SCIP_MAXSTRLEN];

            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "(%s)'", SCIPconsGetName(cons));

            SCIPdebugMessage("split optcumulative constraint <%s>[%d,%d) with %d jobs at time point %d\n",
               SCIPconsGetName(cons), consdata->hmin, consdata->hmax, nvars, split);

            SCIP_CALL( SCIPcreateConsOptcumulative(scip, &splitcons, name, nvars, consdata->vars, consdata->binvars,
                  consdata->durations, consdata->demands, consdata->capacity,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
                  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

            splitconsdata = SCIPconsGetData(splitcons);
            assert(splitconsdata != NULL);

            /* adjust the effective time horizon of the new constraint */
            splitconsdata->hmin = split;
            splitconsdata->hmax = consdata->hmax;

            assert(split < consdata->hmax);

            /* add and release new cumulative constraint */
            SCIP_CALL( SCIPaddCons(scip, splitcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &splitcons) );

            /* adjust the effective time horizon of the constraint */
            consdata->hmax = split;

            assert(consdata->hmin < consdata->hmax);

            (*naddconss)++;
         }

         /* presolve cumulative condition w.r.t. effective horizon by detecting irrelevant variables */
         SCIP_CALL( presolveCumulativeCondition(scip, cons, nfixedvars, nchgcoefs, nchgsides, &cutoff) );

         /* detect implications */
         SCIP_CALL( detectImplications(scip, cons, nchgcoefs, naddconss) );

         /* try to upgrade optcumulative to cumulative constraint which is possible if all remaining binary variables
          * are fixed to one; in case the constraint has no variable left it is removed
          */
         assert(!SCIPinProbing(scip));
         SCIP_CALL( upgradeCons(scip, cons, ndelconss, nupgdconss, &mustpropagate) );
      }
   }

   if( cutoff )
   {
      SCIPdebugMessage("presolving detected a cutoff\n");
      *result = SCIP_CUTOFF;
   }
   else if( oldnfixedvars < *nfixedvars || oldnchgbds < *nchgbds || oldnupgdconss < *nupgdconss || oldndelconss < *ndelconss )
   {
      SCIPdebugMessage("presolving detected %d bound changes\n", *nchgbds - oldnchgbds);
      *result = SCIP_SUCCESS;
   }
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_VAR** binvars;
   int* durations;
   int* demands;
   SCIP_Bool choicevar;
   int nvars;
   int v;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check if the constraint handler wants to participate in the conflict analysis */
   if( !conshdlrdata->conflictanalysis )
   {
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("resolve propagate of optcumulative constraints <%s>\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &durations, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, consdata->nvars) );

   nvars = 0;
   choicevar = FALSE;

   /* collect all activities which are were locally assigned to that machine before the bound change was made */
   for( v = 0; v < consdata->nvars; ++v )
   {
      if( SCIPvarGetLbAtIndex(consdata->binvars[v], bdchgidx, FALSE) > 0.5 )
      {
         vars[nvars] = consdata->vars[v];
         binvars[nvars] = consdata->binvars[v];
         durations[nvars] = consdata->durations[v];
         demands[nvars] = consdata->demands[v];
         nvars++;
      }
      else if( consdata->binvars[v] == infervar )
         choicevar = TRUE;
   }

   assert(nvars > 0);

   if( choicevar )
   {
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( SCIPvarGetLbAtIndex(consdata->binvars[v], bdchgidx, FALSE) > 0.5 )
         {
            SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->binvars[v]) );

            SCIP_CALL( SCIPaddConflictLb(scip, consdata->vars[v], bdchgidx) );
            SCIP_CALL( SCIPaddConflictUb(scip, consdata->vars[v], bdchgidx) );
         }
         else if( consdata->binvars[v] == infervar )
         {
            SCIP_CALL( SCIPaddConflictLb(scip, consdata->vars[v], bdchgidx) );
            SCIP_CALL( SCIPaddConflictUb(scip, consdata->vars[v], bdchgidx) );
         }
      }

      *result = SCIP_SUCCESS;
   }
   else
   {
      SCIP_Bool* explanation;

      SCIP_CALL( SCIPallocBufferArray(scip, &explanation, nvars) );
      BMSclearMemoryArray(explanation, nvars);

      /* resolve propagate of cumulative condition */
      SCIP_CALL( SCIPrespropCumulativeCondition(scip, nvars, vars, durations, demands, consdata->capacity, consdata->hmin, consdata->hmax,
            infervar, inferinfo, boundtype, bdchgidx, relaxedbd, explanation, result) );

      /* if the cumulative constraint handler successfully create an explanation for the propagate we extend this
       * explanation with the required choice variables
       */
      if( *result == SCIP_SUCCESS )
      {
         for( v = 0; v < nvars; ++v )
         {
            if( explanation[v] )
            {
               /* add the lower bounds of the choice variables as part of the initial reason */
               SCIP_CALL( SCIPaddConflictBinvar(scip, binvars[v]) );
            }
         }
      }

      SCIPfreeBufferArray(scip, &explanation);
   }

   /* free all buffers */
   SCIPfreeBufferArray(scip, &demands);
   SCIPfreeBufferArray(scip, &durations);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &binvars);

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   assert(vars != NULL);

   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      if( consdata->downlocks[v] && consdata->uplocks[v] )
      {
         /* the integer start variable should not get rounded in both direction  */
         SCIP_CALL( SCIPaddVarLocksType(scip, vars[v], SCIP_LOCKTYPE_MODEL, nlockspos + nlocksneg, nlockspos + nlocksneg) );
      }
      else if( consdata->downlocks[v]  )
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, vars[v], SCIP_LOCKTYPE_MODEL, nlockspos, nlocksneg) );
      }
      else if( consdata->uplocks[v] )
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, vars[v], SCIP_LOCKTYPE_MODEL, nlocksneg, nlockspos) );
      }

      /* the binary decision variable should not get rounded up; rounding down does not influence the feasibility */
      assert(consdata->binvars[v] != NULL);
      SCIP_CALL( SCIPaddVarLocksType(scip, consdata->binvars[v], SCIP_LOCKTYPE_MODEL, nlocksneg, nlockspos) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveOptcumulative NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveOptcumulative NULL


/** constraint enabling notification method of constraint handler */
#define consEnableOptcumulative NULL


/** constraint disabling notification method of constraint handler */
#define consDisableOptcumulative NULL

/** variable deletion method of constraint handler */
#define consDelvarsOptcumulative NULL

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintOptcumulative)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyOptcumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourceconsdata;
   SCIP_VAR** sourcebinvars;
   SCIP_VAR** sourcevars;
   SCIP_VAR** binvars;
   SCIP_VAR** vars;
   SCIP_Bool success;
   const char* consname;

   int nvars;
   int v;

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert(sourceconsdata != NULL);

   /* get variables of the source constraint */
   sourcebinvars = sourceconsdata->binvars;
   sourcevars = sourceconsdata->vars;
   nvars = sourceconsdata->nvars;

   (*valid) = TRUE;

   if( nvars == 0 )
      return SCIP_OKAY;

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   success = TRUE;

   for( v = 0; v < nvars && success; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcebinvars[v], &binvars[v], varmap, consmap, global, &success) );
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &vars[v], varmap, consmap, global, &success) );
   }

   if( success )
   {
      if( name != NULL )
         consname = name;
      else
         consname = SCIPconsGetName(sourcecons);

      /* copy the logic using the linear constraint copy method */
      SCIP_CALL( SCIPcreateConsOptcumulative(scip, cons, consname, nvars, vars, binvars,
            sourceconsdata->durations, sourceconsdata->demands, sourceconsdata->capacity,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   }
   else
      (*valid) = FALSE;

   /* free buffer array */
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &binvars);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
#define consParseOptcumulative NULL


/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecOptcumulativeBinvars)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_BINVARS_NAME) == 0);
   assert(event != NULL);

   /* collect event information */
   consdata = (SCIP_CONSDATA*)eventdata;
   eventtype = SCIPeventGetType(event);

   switch( eventtype )
   {
   case SCIP_EVENTTYPE_GLBCHANGED:
      consdata->nglbfixedones++;
      break;
   case SCIP_EVENTTYPE_GUBCHANGED:
      consdata->nglbfixedzeros++;
      break;
   case SCIP_EVENTTYPE_LBTIGHTENED:
      consdata->nfixedones++;
      consdata->propagated = FALSE;
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      consdata->nfixedzeros++;
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      consdata->nfixedones--;
      consdata->triedsolving = FALSE;
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      consdata->nfixedzeros--;
      consdata->triedsolving = FALSE;

      if( !SCIPinProbing(scip) )
         consdata->propagated = FALSE;
      break;
   default:
      SCIPerrorMessage("invalid event type %x\n", eventtype);
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXEC(eventExecOptcumulativeIntvars)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_INTVARS_NAME) == 0);
   assert(event != NULL);

   /* collect event information */
   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   /* a bound of a start time variable was tightened; therefore we mark to constraint to create a new local linear
    * relaxation
    */
   if( consdata->nfixedzeros + consdata->nfixedones < consdata->nvars )
      consdata->relaxadded = FALSE;

   if( !SCIPinProbing(scip) )
      consdata->propagated = FALSE;

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for optcumulative constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrOptcumulative(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EVENTHDLR* eventhdlrbinvars;
   SCIP_EVENTHDLR* eventhdlrintvars;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlrbinvars, EVENTHDLR_BINVARS_NAME, EVENTHDLR_BINVARS_DESC,
         eventExecOptcumulativeBinvars, NULL) );

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlrintvars, EVENTHDLR_INTVARS_NAME, EVENTHDLR_INTVARS_DESC,
         eventExecOptcumulativeIntvars, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlrbinvars, eventhdlrintvars) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
         conshdlrCopyOptcumulative,
         consFreeOptcumulative, consInitOptcumulative, consExitOptcumulative,
         consInitpreOptcumulative, consExitpreOptcumulative, consInitsolOptcumulative, consExitsolOptcumulative,
         consDeleteOptcumulative, consTransOptcumulative, consInitlpOptcumulative,
         consSepalpOptcumulative, consSepasolOptcumulative, consEnfolpOptcumulative, consEnforelaxOptcomulative, consEnfopsOptcumulative, consCheckOptcumulative,
         consPropOptcumulative, consPresolOptcumulative, consRespropOptcumulative, consLockOptcumulative,
         consActiveOptcumulative, consDeactiveOptcumulative,
         consEnableOptcumulative, consDisableOptcumulative,
         consDelvarsOptcumulative, consPrintOptcumulative, consCopyOptcumulative, consParseOptcumulative,
         NULL, NULL, NULL,
         conshdlrdata) );

   /* add optcumulative constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/rowrelax",
         "add linear relaxation as LP row (otherwise a knapsack constraint is created)?",
         &conshdlrdata->rowrelax, FALSE, DEFAULT_ROWRELAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/conflictanalysis",
         "participate in conflict analysis?",
         &conshdlrdata->conflictanalysis, FALSE, DEFAULT_CONFLICTANALYSIS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/intervalrelax",
         "create a relaxation for each start and end time point interval",
         &conshdlrdata->intervalrelax, FALSE, DEFAULT_INTERVALRELAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a optcumulative constraint */
SCIP_RETCODE SCIPcreateConsOptcumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   SCIP_VAR**            binvars,            /**< array of variable representing if the job has to be processed on this machine */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsOptcumulative() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the optcumulative constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("optcumulative constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* the optcumulative constraint handler currently does not support modifiable constraints */
   assert(modifiable == FALSE);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, nvars, vars, binvars, durations, demands, capacity, check) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlrbinvars != NULL);
      assert(conshdlrdata->eventhdlrintvars != NULL);

      assert(consdata->nglbfixedones == 0);
      assert(consdata->nglbfixedzeros == 0);

      /* catch bound change events of variables */
      SCIP_CALL( catchAllEvents(scip, *cons, conshdlrdata->eventhdlrbinvars, conshdlrdata->eventhdlrintvars) );
   }

   return SCIP_OKAY;
}
