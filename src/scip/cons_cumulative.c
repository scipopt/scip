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

/**@file   cons_cumulative.c
 * @brief  constraint handler for cumulative constraints
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Jens Schulz
 *
 * Given:
 * - a set of jobs, represented by their integer start time variables \f$S_j\f$, their array of processing times \f$p_j\f$ and of
 *   their demands \f$d_j\f$.
 * - an integer resource capacity \f$C\f$
 *
 * The cumulative constraint ensures that for each point in time \f$t\f$ \f$\sum_{j: S_j \leq t < S_j + p_j} d_j \leq C\f$ holds.
 *
 * Separation:
 * - can be done using binary start time model, see Pritskers, Watters and Wolfe
 * - or by just separating relatively weak cuts on the start time variables
 *
 * Propagation:
 * - time tabling, Klein & Scholl (1999)
 * - Edge-finding from Petr Vilim, adjusted and simplified for dynamic repropagation
 *   (2009)
 * - energetic reasoning, see Baptiste, Le Pape, Nuijten (2001)
 *
 *
 * TODOS:
 * - normalize demands and capacities of each cumulative constraint
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_cumulative.h"
#include "scip/cons_linking.h"
#include "scip/cons_knapsack.h"
#include "scip/scipdefplugins.h"

/**@name Constraint handler properties
 *
 * @{
 */

/* constraint handler properties */
#define CONSHDLR_NAME          "cumulative"
#define CONSHDLR_DESC          "cumulative constraint handler"
#define CONSHDLR_SEPAPRIORITY   2100000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -2040000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -3030000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING   SCIP_PROPTIMING_BEFORELP

/**@} */

/**@name Default parameter values
 *
 * @{
 */

/* default parameter values */

/* separation */
#define DEFAULT_USEBINVARS             FALSE /**< should the binary representation be used? */
#define DEFAULT_LOCALCUTS              FALSE /**< should cuts be added only locally? */
#define DEFAULT_USECOVERCUTS            TRUE /**< should covering cuts be added? */
#define DEFAULT_CUTSASCONSS             TRUE /**< should the cuts be created as knapsack constraints? */
#define DEFAULT_SEPAOLD                 TRUE /**< shall old sepa algo be applied? */

/* propagation */
#define DEFAULT_CORETIMES               TRUE /**< should core-times be propagated (time tabling)? */
#define DEFAULT_OVERLOAD                TRUE /**< should edge finding be used to detect an overload? */

/* presolving */
#define DEFAULT_DUALPRESOLVE            TRUE /**< should dual presolving be applied? */
#define DEFAULT_COEFTIGHTENING         FALSE /**< should coeffisient tightening be applied? */
#define DEFAULT_NORMALIZE               TRUE /**< should demands and capacity be normalized? */
#define DEFAULT_MAXNODES             10000LL /**< number of branch-and-bound nodes to solve an independent cumulative constraint  (-1: no limit) */

/* enforcement */
#define DEFAULT_FILLBRANCHCANDS        FALSE /**< should branching candidates be added to storage? */

/**@} */

/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_NAME         "cumulative"
#define EVENTHDLR_DESC         "bound change event handler for cumulative constraints"

/**@} */

/*
 * Data structures
 */

/** constraint data for cumulative constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< array of variable representing the start time of each job */
   SCIP_Bool*            downlocks;          /**< array to store if the variable has a down lock */
   SCIP_Bool*            uplocks;            /**< array to store if the variable has a down lock */
   SCIP_CONS**           linkingconss;       /**< array of linking constraints for the integer variables */
   SCIP_ROW**            demandrows;         /**< array of rows of linear relaxation of this problem */
   SCIP_ROW**            scoverrows;         /**< array of rows of small cover cuts of this problem */
   SCIP_ROW**            bcoverrows;         /**< array of rows of big cover cuts of this problem */
   int*                  demands;            /**< array containing corresponding demands */
   int*                  durations;          /**< array containing corresponding durations */
   SCIP_Real             resstrength1;       /**< stores the resource strength 1*/
   SCIP_Real             resstrength2;       /**< stores the resource strength 2 */
   SCIP_Real             cumfactor1;         /**< stroes the cumulativeness of the constraint */
   SCIP_Real             disjfactor1;        /**< stores the disjunctiveness of the constraint */
   SCIP_Real             disjfactor2;        /**< stores the disjunctiveness of the constraint */
   SCIP_Real             estimatedstrength;
   int                   nvars;              /**< number of variables */
   int                   varssize;           /**< size of the arrays */
   int                   ndemandrows;        /**< number of rows of cumulative constrint for linear relaxation */
   int                   demandrowssize;     /**< size of array rows of demand rows */
   int                   nscoverrows;        /**< number of rows of small cover cuts */
   int                   scoverrowssize;     /**< size of array of small cover cuts */
   int                   nbcoverrows;        /**< number of rows of big cover cuts */
   int                   bcoverrowssize;     /**< size of array of big cover cuts */
   int                   capacity;           /**< available cumulative capacity */

   int                   hmin;               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax;               /**< right bound of time axis to be considered  (not including hmax) */

   unsigned int          normalized:1;       /**< is the constraint normalized */
   unsigned int          covercuts:1;        /**< cover cuts are created? */
   unsigned int          propagated:1;       /**< is constraint propagted */

#ifdef SCIP_STATISTIC
   int                   maxpeak;
   int                   nirrelevantjobs;
   int                   nalwaysruns;
   int                   ndualfixs;
   int                   nremovedlocks;
   int                   ndecomps;
#endif
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */

   SCIP_Bool             usebinvars;         /**< should the binary variables be used? */
   SCIP_Bool             cutsasconss;        /**< should the cumulative constraint create cuts as knapsack constraints? */
   SCIP_Bool             coretimes;          /**< should core-times be propagated (time tabling)? */
   SCIP_Bool             overload;           /**< should edge finding be used to detect an overload? */
   SCIP_Bool             localcuts;          /**< should cuts be added only locally? */
   SCIP_Bool             usecovercuts;       /**< should covering cuts be added? */
   SCIP_Bool             sepaold;            /**< shall old sepa algo be applied? */


   SCIP_Bool             fillbranchcands;    /**< should branching candidates be added to storage? */

   SCIP_Bool             dualpresolve;       /**< should dual presolving be applied? */
   SCIP_Bool             coeftightening;     /**< should coeffisient tightening be applied? */
   SCIP_Bool             normalize;          /**< should demands and capacity be normalized? */

   SCIP_Longint          maxnodes;           /**< number of branch-and-bound nodes to solve an independent cumulative constraint  (-1: no limit) */

#ifdef SCIP_STATISTIC
   int                   nirrelevantjobs;
   int                   nalwaysruns;
   int                   ndualfixs;
   int                   nremovedlocks;
   int                   ndecomps;
   int                   nallconsdualfixs;
#endif
};

/**@name Inference Information Methods
 *
 *  An inference information can be passed with each domain reduction to SCIP. This information is passed back to the
 *  constraint handler if the corresponding bound change has to be explained. It can be used to store information which
 *  help to construct a reason/explanation for a bound change. The inference information is limited to size of integer.
 *
 *  In case of the cumulative constraint handler we store the used propagation algorithms for that particular bound
 *  change and the earliest start and latest completion time of all jobs in the conflict set.
 *
 * @{
 */

/** Propagation rules */
enum Proprule
{
   PROPRULE_1_CORETIMES          = 1,        /**< core-time propagator */
   PROPRULE_2_EDGEFINDING        = 2,        /**< edge-finder */
   PROPRULE_3_ENERGETICREASONING = 3         /**< energetic reasoning */
};
typedef enum Proprule PROPRULE;

/** inference information */
struct InferInfo
{
   union
   {
      /** struct to use the inference information */
      struct
      {
         unsigned int    proprule:2;         /**< propagation rule that was applied */
         unsigned int    data1:15;           /**< data field one */
         unsigned int    data2:15;           /**< data field two */
      } asbits;
      int                asint;              /**< inference information as a single int value */
   } val;
};
typedef struct InferInfo INFERINFO;

/** converts an integer into an inference information */
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
PROPRULE inferInfoGetProprule(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return (PROPRULE) inferinfo.val.asbits.proprule;
}

/** returns data field one of the inference information */
static
int inferInfoGetData1(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asbits.data1;
}

#if 0
/** returns data field two of the inference information */
static
int inferInfoGetData2(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asbits.data2;
}
#endif


/** constructs an inference information out of a propagation rule, an earliest start and a latest completion time */
static
INFERINFO getInferInfo(
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   int                   data1,              /**< data field one */
   int                   data2               /**< data field two */
   )
{
   INFERINFO inferinfo;

   inferinfo.val.asbits.proprule = proprule; /*lint !e641*/
   inferinfo.val.asbits.data1 = data1; /*lint !e732*/
   inferinfo.val.asbits.data2 = data2; /*lint !e732*/

   return inferinfo;
}

/**@} */

/*
 * Local methods
 */

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
   assert(SCIPisFeasIntegral(scip, bound));
   assert(SCIPisFeasEQ(scip, bound, (SCIP_Real)(int)(bound + 0.5)));

   return (int)(bound + 0.5);
}
#else
#define convertBoundToInt(x, y) ((int)((y) + 0.5))
#endif

/** returns the implied earliest start time */
static
SCIP_RETCODE computeImpliedEst(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the implied est should be returned */
   SCIP_HASHMAP*         addedvars,          /**< hash map containig the variable which are already added */
   int*                  est                 /**< pointer to store the implied earliest start time */
   )
{
#if 0
   SCIP_VAR** vbdvars;
   SCIP_VAR* vbdvar;
   SCIP_Real* vbdcoefs;
   SCIP_Real* vbdconsts;
   void* image;
   int nvbdvars;
   int v;
#endif

   (*est) = convertBoundToInt(scip, SCIPvarGetLbLocal(var));

#if 0
   /* the code contains a bug; we need to check if an implication forces that the jobs do not run in parallel */

   /* check if ????? */
   nvbdvars = SCIPvarGetNVlbs(var);
   vbdvars = SCIPvarGetVlbVars(var);
   vbdcoefs = SCIPvarGetVlbCoefs(var);
   vbdconsts = SCIPvarGetVlbConstants(var);

   for( v = 0; v < nvbdvars; ++v )
   {
      vbdvar = vbdvars[v];
      assert(vbdvar != NULL);

      image = SCIPhashmapGetImage(addedvars, (void*)vbdvar);

      if( image != NULL && SCIPisEQ(scip, vbdcoefs[v], 1.0 ) )
      {
         int duration;
         int vbdconst;

         duration = (int)(size_t)image;
         vbdconst = convertBoundToInt(scip, vbdconsts[v]);

         SCIPdebugMessage("check implication <%s>[%g,%g] >= <%s>[%g,%g] + <%g>\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
            SCIPvarGetName(vbdvar), SCIPvarGetLbLocal(vbdvar), SCIPvarGetUbLocal(vbdvar), vbdconsts[v]);

         if( duration >= vbdconst )
         {
            int impliedest;

            impliedest =  convertBoundToInt(scip, SCIPvarGetUbLocal(vbdvar)) + duration;

            if( (*est) < impliedest )
            {
               (*est) = impliedest;

               SCIP_CALL( SCIPhashmapRemove(addedvars, (void*)vbdvar) );
            }
         }
      }
   }
#endif

   return SCIP_OKAY;
}

/** returns the implied latest completion time */
static
SCIP_RETCODE computeImpliedLct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the implied est should be returned */
   int                   duration,           /**< duration of the given job */
   SCIP_HASHMAP*         addedvars,          /**< hash map containig the variable which are already added */
   int*                  lct                 /**< pointer to store the implied latest completion time */
   )
{
#if 0
   SCIP_VAR** vbdvars;
   SCIP_VAR* vbdvar;
   SCIP_Real* vbdcoefs;
   SCIP_Real* vbdconsts;
   int nvbdvars;
   int v;
#endif

   (*lct) = convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + duration;

#if 0
   /* the code contains a bug; we need to check if an implication forces that the jobs do not run in parallel */

   /* check if ????? */
   nvbdvars = SCIPvarGetNVubs(var);
   vbdvars = SCIPvarGetVubVars(var);
   vbdcoefs = SCIPvarGetVubCoefs(var);
   vbdconsts = SCIPvarGetVubConstants(var);

   for( v = 0; v < nvbdvars; ++v )
   {
      vbdvar = vbdvars[v];
      assert(vbdvar != NULL);

      if( SCIPhashmapExists(addedvars, (void*)vbdvar) &&  SCIPisEQ(scip, vbdcoefs[v], 1.0 ) )
      {
         int vbdconst;

         vbdconst = convertBoundToInt(scip, -vbdconsts[v]);

         SCIPdebugMessage("check implication <%s>[%g,%g] <= <%s>[%g,%g] + <%g>\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
            SCIPvarGetName(vbdvar), SCIPvarGetLbLocal(vbdvar), SCIPvarGetUbLocal(vbdvar), vbdconsts[v]);

         if( duration >= -vbdconst )
         {
            int impliedlct;

            impliedlct = convertBoundToInt(scip, SCIPvarGetLbLocal(vbdvar));

            if( (*lct) > impliedlct )
            {
               (*lct) = impliedlct;

               SCIP_CALL( SCIPhashmapRemove(addedvars, (void*)vbdvar) );
            }
         }
      }
   }
#endif

   return SCIP_OKAY;
}

/** creates the worst case resource profile, that is, all jobs are inserted with the earliest start and latest
 *  completion time
 */
static
SCIP_RETCODE createWorstCaseProfile(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROFILE*         profile,            /**< resource profile */
   SCIP_CONSDATA*        consdata            /**< constraint data to use for filling the worst case profile */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_HASHMAP* addedvars;
   int* demands;
   int* perm;
   int duration;
   int nvars;
   int impliedest;
   int est;
   int impliedlct;
   int lct;
   int v;

   nvars = consdata->nvars;
   vars = consdata->vars;

   assert(SCIPgetDepth(scip) <= 0);

   /* create hash map for variables which are added, mapping to their duration */
   SCIP_CALL( SCIPhashmapCreate(&addedvars, SCIPblkmem(scip), SCIPcalcHashtableSize(nvars)) );

   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, nvars) );

   /* sort variables  w.r.t. job demands */
   for( v = 0; v < nvars; ++v )
   {
      demands[v] = consdata->demands[v];
      perm[v] = v;
   }
   SCIPsortDownIntInt(demands, perm, nvars);

   /* add each job with its earliest start and latest completion time into the resource profile */
   for( v = 0; v < nvars; ++v )
   {
      int idx;

      idx = perm[v];
      assert(idx >= 0 && idx < nvars),

         var = vars[idx];
      assert(var != NULL);

      duration = consdata->durations[idx];
      assert(duration > 0);

      est = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      SCIP_CALL( computeImpliedEst(scip, var, addedvars, &impliedest) );

      lct = convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + duration;
      SCIP_CALL( computeImpliedLct(scip, var, duration, addedvars, &impliedlct) );

      if( impliedest < impliedlct )
      {
         SCIP_Bool infeasible;
         int pos;

         SCIP_CALL( SCIPprofileInsertCore(profile, impliedest, impliedlct, demands[v], &pos, &infeasible) );
         assert(!infeasible);
         assert(pos == -1);
      }

      if( est == impliedest && lct == impliedlct )
      {
         SCIP_CALL( SCIPhashmapInsert(addedvars, (void*)var, (void*)(size_t)duration) );
      }
   }

   SCIPfreeBufferArray(scip, &demands);
   SCIPfreeBufferArray(scip, &perm);

   SCIPhashmapFree(&addedvars);

   return SCIP_OKAY;
}

/** collects all necessary binary variables to represent the jobs which can be active at time point of interest */
static
SCIP_RETCODE collectBinaryVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR***           vars,               /**< pointer to the array to store the binary variables */
   int**                 coefs,              /**< pointer to store the coefficients */
   int*                  nvars,              /**< number if collect binary variables */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int                   curtime,            /**< current point in time */
   int                   nstarted,           /**< number of jobs that start before the curtime or at curtime */
   int                   nfinished           /**< number of jobs that finished before curtime or at curtime */
   )
{
   SCIP_VAR** binvars;
   SCIP_VAR* var;
   int nbinvars;
   int nrowvars;
   int startindex;
   int endtime;
   int duration;
   int demand;
   int varidx;
   int offset;
   int minub;
   int size;

   size = 10;
   nrowvars = 0;
   startindex = nstarted - 1;

   SCIP_CALL( SCIPallocBufferArray(scip, vars, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, coefs, size) );

   /* search for the (nstarted - nfinished) jobs which are active at curtime */
   while( nstarted - nfinished > nrowvars )
   {
      /* collect job information */
      varidx = startindices[startindex];
      assert(varidx >= 0 && varidx < consdata->nvars);

      var = consdata->vars[varidx];
      duration = consdata->durations[varidx];
      demand = consdata->demands[varidx];
      assert(var != NULL);

      endtime = convertBoundToInt(scip, SCIPvarGetUbGlobal(var)) + duration;

      /* check the end time of this job is larger than the curtime; in this case the job is still running */
      if( endtime > curtime )
      {
         int tau;  /* counter from curtime - duration + 1 to curtime */

         /* check if the linking constraints exists */
         assert(SCIPexistsConsLinking(scip, var));
         assert(SCIPgetConsLinking(scip, var) != NULL);
         assert(SCIPgetConsLinking(scip, var) == consdata->linkingconss[varidx]);

         /* collect linking constraint information */
         SCIP_CALL( SCIPgetBinvarsLinking(scip, consdata->linkingconss[varidx], &binvars, &nbinvars) );
         offset = SCIPgetOffsetLinking(scip, consdata->linkingconss[varidx]);

         minub = MIN(curtime, endtime - duration);

         for (tau = MAX(curtime - duration + 1, offset); tau <= minub; ++tau )
         {
            assert(tau >= offset && tau < nbinvars + offset);
            assert(binvars[tau-offset] != NULL);

            /* ensure array proper array size */
            if( size == *nvars )
            {
               size *= 2;
               SCIP_CALL( SCIPreallocBufferArray(scip, vars, size) );
               SCIP_CALL( SCIPreallocBufferArray(scip, coefs, size) );
            }

            (*vars)[*nvars] = binvars[tau-offset];
            (*coefs)[*nvars] = demand;
            (*nvars)++;
         }
         nrowvars++;
      }

      startindex--;
   }

   return SCIP_OKAY;
}

/** collect all integer variable which belong to jobs which can run at the point of interest */
static
SCIP_RETCODE collectIntVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR***           activevars,         /**< jobs that are currently running */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int                   curtime,            /**< current point in time */
   int                   nstarted,           /**< number of jobs that start before the curtime or at curtime */
   int                   nfinished,          /**< number of jobs that finished before curtime or at curtime */
   SCIP_Bool             lower,              /**< shall cuts be created due to lower or upper bounds? */
   int*                  lhs                 /**< lhs for the new row sum of lbs + minoffset */
   )
{
   SCIP_VAR* var;
   int startindex;
   int endtime;
   int duration;
   int starttime;

   int varidx;
   int sumofstarts;
   int mindelta;
   int counter;

   counter = 0;
   sumofstarts = 0;

   mindelta = INT_MAX;

   startindex = nstarted - 1;

   /* search for the (nstarted - nfinished) jobs which are active at curtime */
   while( nstarted - nfinished > counter )
   {
      assert(startindex >= 0);

      /* collect job information */
      varidx = startindices[startindex];
      assert(varidx >= 0 && varidx < consdata->nvars);

      var = consdata->vars[varidx];
      duration = consdata->durations[varidx];
      assert(duration > 0);
      assert(var != NULL);

      if( lower )
         starttime = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      else
         starttime = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

      endtime = starttime + duration;

      /* check the end time of this job is larger than the curtime; in this case the job is still running */
      if( endtime > curtime )
      {
         (*activevars)[counter] = var;
         sumofstarts += starttime;
         mindelta = MIN(mindelta, endtime - curtime); /* this amount of schifting holds for lb and ub */
         counter++;
      }

      startindex--;
   }

   assert(mindelta > 0);
   *lhs = lower ? sumofstarts + mindelta : sumofstarts - mindelta;

   return SCIP_OKAY;
}

/** initialize the sorted event point arrays */
static
void createSortedEventpoints(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations per start time variable */
   int*                  starttimes,         /**< array to store sorted start events */
   int*                  endtimes,           /**< array to store sorted end events */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int*                  endindices,         /**< permutation with rspect to the end times */
   SCIP_Bool             local               /**< shall local bounds be used */
   )
{
   SCIP_VAR* var;
   int j;

   assert(vars != NULL || nvars == 0);

   /* assign variables, start and endpoints to arrays */
   for ( j = 0; j < nvars; ++j )
   {
      assert(vars != NULL);

      var = vars[j];
      assert(var != NULL);

      if( local )
         starttimes[j] = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      else
         starttimes[j] = convertBoundToInt(scip, SCIPvarGetLbGlobal(var));

      startindices[j] = j;

      if( local )
         endtimes[j] = convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + durations[j];
      else
         endtimes[j] = convertBoundToInt(scip, SCIPvarGetUbGlobal(var)) + durations[j];

      endindices[j] = j;

   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues (and sort the indices in the same way) */
   SCIPsortIntInt(starttimes, startindices, j);
   SCIPsortIntInt(endtimes, endindices, j);
}

/** initialize the sorted event point arrays w.r.t. the given primal solutions */
static
void createSortedEventpointsSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations per start time variable */
   int*                  starttimes,         /**< array to store sorted start events */
   int*                  endtimes,           /**< array to store sorted end events */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int*                  endindices         /**< permutation with rspect to the end times */
   )
{
   SCIP_VAR* var;
   int j;

   assert(vars != NULL || nvars == 0);

   /* assign variables, start and endpoints to arrays */
   for ( j = 0; j < nvars; ++j )
   {
      assert(vars != NULL);

      var = vars[j];
      assert(var != NULL);

      starttimes[j] = convertBoundToInt(scip, SCIPgetSolVal(scip, sol, var));
      startindices[j] = j;

      endtimes[j] = convertBoundToInt(scip, SCIPgetSolVal(scip, sol, var)) + durations[j];
      endindices[j] = j;

   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues (and sort the indices in the same way) */
   SCIPsortIntInt(starttimes, startindices, j);
   SCIPsortIntInt(endtimes, endindices, j);
}

/** initialize the sorted event point arrays
 *
 * @todo Check the separation process!
 */
static
void createSelectedSortedEventpointsSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   int*                  starttimes,         /**< array to store sorted start events */
   int*                  endtimes,           /**< array to store sorted end events */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int*                  endindices,         /**< permutation with rspect to the end times */
   int*                  nvars,              /**< number of variables that are integral */
   SCIP_Bool             lower               /**< shall the constraints be derived for lower or upper bounds? */
   )
{
   SCIP_VAR* var;
   int tmpnvars;
   int j;

   tmpnvars = consdata->nvars;
   *nvars = 0;

   /* assign variables, start and endpoints to arrays */
   for ( j = 0; j < tmpnvars; ++j )
   {
      var = consdata->vars[j];
      assert(var != NULL);

      if( lower )
      {
         /* only consider jobs that are at their lower or upper bound */
         if( !SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, var))
            || !SCIPisFeasEQ(scip, SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var)) )
            continue;

         if( consdata->durations[j] == 0 || consdata->demands[j] == 0 )
            continue;

         starttimes[*nvars] = convertBoundToInt(scip, SCIPgetSolVal(scip, sol, var));
         startindices[*nvars] = j;

         endtimes[*nvars] =  starttimes[*nvars] + consdata->durations[j];
         endindices[*nvars] = j;

         (*nvars) = *nvars + 1;

         SCIPdebugMessage("lower bounds are considered:\n");
         SCIPdebugMessage("%d: job[%d] starttime %d, endtime = %d, demand = %d\n ", *nvars-1,
            startindices[*nvars-1], starttimes[*nvars-1], starttimes[*nvars-1] + consdata->durations[startindices[*nvars-1]],
            consdata->demands[startindices[*nvars-1]]);
      }
      else
      {
         if( !SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, var))
            || !SCIPisFeasEQ(scip, SCIPgetSolVal(scip, sol, var), SCIPvarGetUbLocal(var)) )
            continue;

         starttimes[*nvars] = convertBoundToInt(scip, SCIPgetSolVal(scip, sol, var));
         startindices[*nvars] = j;

         endtimes[*nvars] =  starttimes[*nvars] + consdata->durations[j];
         endindices[*nvars] = j;

         (*nvars) = *nvars + 1;

         SCIPdebugMessage("upper bounds are considered:\n");
         SCIPdebugMessage("%d: job[%d] starttime %d, endtime = %d, demand = %d\n ", *nvars-1,
            startindices[*nvars-1], starttimes[*nvars-1], starttimes[*nvars-1] + consdata->durations[startindices[*nvars-1]],
            consdata->demands[startindices[*nvars-1]]);
      }
   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues (and sort the indices in the same way) */
   SCIPsortIntInt(starttimes, startindices, *nvars);
   SCIPsortIntInt(endtimes, endindices, *nvars);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("sorted output %d\n", *nvars);

   for ( j = 0; j < *nvars; ++j )
   {
      SCIPdebugMessage("%d: job[%d] starttime %d, endtime = %d, demand = %d\n", j,
         startindices[j], starttimes[j], starttimes[j] + consdata->durations[startindices[j]],
         consdata->demands[startindices[j]]);
   }

   for ( j = 0; j < *nvars; ++j )
   {
      SCIPdebugMessage("%d: job[%d] endtime %d,  demand = %d\n", j, endindices[j], endtimes[j],
         consdata->demands[endindices[j]]);
   }
   SCIPdebugMessage("capacity = %d\n", consdata->capacity);
#endif
}

#ifdef SCIP_STATISTIC
/** this method checks for relevant intervals for energetic reasoning */
static
SCIP_RETCODE computeRelevantEnergyIntervals(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   int**                 timepoints,         /**< array to store relevant points in time */
   SCIP_Real**           cumulativedemands,  /**< array to store the estimated cumulative demand for each point in time */
   int*                  ntimepoints,        /**< pointer to store the number of timepoints */
   int*                  maxdemand,          /**< pointer to store maximum over all demands */
   SCIP_Real*            minfreecapacity     /**< pointer to store the minimum free capacity */
   )
{
   int* starttimes;         /* stores when each job is starting */
   int* endtimes;           /* stores when each job ends */
   int* startindices;       /* we will sort the startsolvalues, thus we need to know wich index of a job it corresponds to */
   int* endindices;         /* we will sort the endsolvalues, thus we need to know wich index of a job it corresponds to */

   SCIP_Real totaldemand;
   int curtime;             /* point in time which we are just checking */
   int endindex;            /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int j;

   assert( scip != NULL );
   assert(durations != NULL);
   assert(demands != NULL);
   assert(capacity >= 0);

   /* if no activities are associated with this cumulative then this constraint is redundant */
   if( nvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   /* create event point arrays */
   createSortedEventpoints(scip, nvars, vars, durations, starttimes, endtimes, startindices, endindices, TRUE);

   endindex = 0;
   totaldemand = 0.0;

   *ntimepoints = 0;
   (*timepoints)[0] = starttimes[0];
   (*cumulativedemands)[0] = 0;
   *maxdemand = 0;

   /* check each startpoint of a job whether the capacity is kept or not */
   for( j = 0; j < nvars; ++j )
   {
      int lct;
      int idx;

      curtime = starttimes[j];

      if( curtime >= hmax )
         break;

      /* free all capacity usages of jobs the are no longer running */
      while( endindex < nvars && endtimes[endindex] <= curtime )
      {
         int est;

         if( (*timepoints)[*ntimepoints] < endtimes[endindex] )
         {
            (*ntimepoints)++;
            (*timepoints)[*ntimepoints] = endtimes[endindex];
            (*cumulativedemands)[*ntimepoints] = 0;
         }

         idx = endindices[endindex];
         est = convertBoundToInt(scip, SCIPvarGetLbLocal(vars[idx]));
         totaldemand -= (SCIP_Real) demands[idx] * durations[idx] / (endtimes[endindex] - est);
         endindex++;

         (*cumulativedemands)[*ntimepoints] = totaldemand;
      }

      idx = startindices[j];
      lct = convertBoundToInt(scip, SCIPvarGetUbLocal(vars[idx]) + durations[idx]);
      totaldemand += (SCIP_Real) demands[idx] * durations[idx] / (lct - starttimes[j]);

      if( (*timepoints)[*ntimepoints] < curtime )
      {
         (*ntimepoints)++;
         (*timepoints)[*ntimepoints] = curtime;
         (*cumulativedemands)[*ntimepoints] = 0;
      }

      (*cumulativedemands)[*ntimepoints] = totaldemand;

      /* add the relative capacity requirements for all job which start at the curtime */
      while( j+1 < nvars && starttimes[j+1] == curtime )
      {
         ++j;
         idx = startindices[j];
         lct = convertBoundToInt(scip, SCIPvarGetUbLocal(vars[idx]) + durations[idx]);
         totaldemand += (SCIP_Real) demands[idx] * durations[idx] / (lct - starttimes[j]);

         (*cumulativedemands)[*ntimepoints] = totaldemand;
      }
   } /*lint --e{850}*/

   /* free all capacity usages of jobs that are no longer running */
   while( endindex < nvars/* && endtimes[endindex] < hmax*/)
   {
      int est;
      int idx;

      if( (*timepoints)[*ntimepoints] < endtimes[endindex] )
      {
         (*ntimepoints)++;
         (*timepoints)[*ntimepoints] = endtimes[endindex];
         (*cumulativedemands)[*ntimepoints] = 0;
      }

      idx = endindices[endindex];
      est = convertBoundToInt(scip, SCIPvarGetLbLocal(vars[idx]));
      totaldemand -= (SCIP_Real) demands[idx] * durations[idx] / (endtimes[endindex] - est);
      (*cumulativedemands)[*ntimepoints] = totaldemand;

      ++endindex;
   }

   (*ntimepoints)++;
   /* compute minimum free capacity */
   (*minfreecapacity) = INT_MAX;
   for( j = 0; j < *ntimepoints; ++j )
   {
      if( (*timepoints)[j] >= hmin && (*timepoints)[j] < hmax )
         *minfreecapacity = MIN( *minfreecapacity, (SCIP_Real)capacity - (*cumulativedemands)[j] );
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endtimes);
   SCIPfreeBufferArray(scip, &starttimes);

   return SCIP_OKAY;
}

/** evaluates the cumulativeness and disjointness factor of a cumulative constraint */
static
SCIP_RETCODE evaluateCumulativeness(
   SCIP*                 scip,               /**< pointer to scip */
   SCIP_CONS*            cons                /**< cumulative constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int nvars;
   int v;
   int capacity;

   /* output values: */
   SCIP_Real disjfactor2; /* (peak-capacity)/capacity * (large demands/nvars_t) */
   SCIP_Real cumfactor1;
   SCIP_Real resstrength1; /* overall strength */
   SCIP_Real resstrength2; /* timepoint wise maximum */

   /* helpful variables: */
   SCIP_Real globalpeak;
   SCIP_Real globalmaxdemand;

   /* get constraint data structure */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   capacity = consdata->capacity;
   globalpeak = 0.0;
   globalmaxdemand = 0.0;

   disjfactor2 = 0.0;
   cumfactor1 = 0.0;
   resstrength2 = 0.0;

   /* check each starting time (==each job, but inefficient) */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Real peak;
      SCIP_Real maxdemand;
      SCIP_Real deltademand;
      int ndemands;
      int nlarge;

      int timepoint;
      int j;
      timepoint = convertBoundToInt(scip, SCIPvarGetLbLocal(consdata->vars[v]));
      peak = consdata->demands[v];
      ndemands = 1;
      maxdemand = 0;
      nlarge = 0;

      if( consdata->demands[v] > capacity / 3 )
         nlarge++;

      for( j = 0; j < nvars; ++j )
      {
         int lb;

         if( j == v )
            continue;

         maxdemand = 0.0;
         lb = convertBoundToInt(scip, SCIPvarGetLbLocal(consdata->vars[j]));

         if( lb <= timepoint && lb + consdata->durations[j] > timepoint )
         {
            peak += consdata->demands[j];
            ndemands++;

            if( consdata->demands[j] > consdata->capacity / 3 )
               nlarge++;
         }
      }

      deltademand = (SCIP_Real)peak / (SCIP_Real)ndemands;
      globalpeak = MAX(globalpeak, peak);
      globalmaxdemand = MAX(globalmaxdemand, maxdemand);

      if( peak > capacity )
      {
         disjfactor2 = MAX( disjfactor2, (peak-(SCIP_Real)capacity)/peak * (nlarge/(SCIP_Real)ndemands) );
         cumfactor1 = MAX( cumfactor1, (peak-capacity)/peak * (capacity-deltademand)/(SCIP_Real)capacity );
         resstrength2 = MAX(resstrength2, (capacity-maxdemand)/(peak-maxdemand) );
      }
   }

   resstrength1 = (capacity-globalmaxdemand) / (globalpeak-globalmaxdemand);

   consdata->maxpeak = convertBoundToInt(scip, globalpeak);
   consdata->disjfactor2 = disjfactor2;
   consdata->cumfactor1 = cumfactor1;
   consdata->resstrength2 = resstrength2;
   consdata->resstrength1 = resstrength1;

   /* get estimated res strength */
   {
      int* timepoints;
      SCIP_Real* estimateddemands;
      int ntimepoints;
      int maxdemand;
      SCIP_Real minfreecapacity;

      SCIP_CALL( SCIPallocBufferArray(scip, &timepoints, 2*nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &estimateddemands, 2*nvars) );

      ntimepoints = 0;
      minfreecapacity = INT_MAX;


      SCIP_CALL( computeRelevantEnergyIntervals(scip, nvars, consdata->vars,
            consdata->durations, consdata->demands,
            capacity, consdata->hmin, consdata->hmax, &timepoints, &estimateddemands,
            &ntimepoints, &maxdemand, &minfreecapacity) );

      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &timepoints);
      SCIPfreeBufferArray(scip, &estimateddemands);

      consdata->estimatedstrength = (SCIP_Real)(capacity - minfreecapacity) / (SCIP_Real) capacity;
   }

   SCIPstatisticMessage("cumulative constraint<%s>: DISJ1=%g, DISJ2=%g, CUM=%g, RS1 = %g, RS2 = %g, EST = %g\n",
      SCIPconsGetName(cons), consdata->disjfactor1, disjfactor2, cumfactor1, resstrength1, resstrength2,
      consdata->estimatedstrength);

   return SCIP_OKAY;
}
#endif

/**@} */

/**@name Constraint handler data
 *
 * Method used to create and free the constraint handler data when including and removing the cumulative constraint
 * handler.
 *
 * @{
 */

/** creates constaint handler data for cumulative constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   /* create precedence constraint handler data */
   assert(conshdlrdata != NULL);
   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   /* get event handler for updating linear constraint activity bounds */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for linear constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

#ifdef SCIP_STATISTIC
   (*conshdlrdata)->nirrelevantjobs = 0;
   (*conshdlrdata)->nalwaysruns = 0;
   (*conshdlrdata)->ndualfixs = 0;
   (*conshdlrdata)->nremovedlocks = 0;
   (*conshdlrdata)->ndecomps = 0;
   (*conshdlrdata)->nallconsdualfixs = 0;
#endif

   return SCIP_OKAY;
}

/** frees constraint handler data for logic or constraint handler */
static
void conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeMemory(scip, conshdlrdata);
}

/**@} */


/**@name Constraint data methods
 *
 * @{
 */

/** catches bound change events for all variables in transformed cumulative constraint */
static
SCIP_RETCODE consdataCatchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< cumulative constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);

   /* catch event for every single variable */
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[v],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );
   }

   return SCIP_OKAY;
}

/** drops events for variable at given position */
static
SCIP_RETCODE consdataDropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< cumulative constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars[pos] != NULL);

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos],
         SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

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
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);

   /* drop event of every single variable */
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_CALL( consdataDropEvents(scip, consdata, eventhdlr, v) );
   }

   return SCIP_OKAY;
}

/** creates constraint data of cumulative constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to consdata */
   SCIP_VAR**            vars,               /**< array of integer variables */
   SCIP_CONS**           linkingconss,       /**< array of linking constraints for the integer variables, or NULL */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   nvars,              /**< number of variables */
   int                   capacity,           /**< available cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax                /**< right bound of time axis to be considered (not including hmax) */
   )
{
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(vars != NULL || nvars > 0);
   assert(demands != NULL);
   assert(durations != NULL);
   assert(capacity >= 0);
   assert(hmin >= 0);
   assert(hmin < hmax);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->hmin = hmin;
   (*consdata)->hmax = hmax;

   (*consdata)->capacity = capacity;
   (*consdata)->demandrows = NULL;
   (*consdata)->demandrowssize = 0;
   (*consdata)->ndemandrows = 0;
   (*consdata)->scoverrows = NULL;
   (*consdata)->nscoverrows = 0;
   (*consdata)->scoverrowssize = 0;
   (*consdata)->bcoverrows = NULL;
   (*consdata)->nbcoverrows = 0;
   (*consdata)->bcoverrowssize = 0;
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->covercuts = FALSE;
   (*consdata)->normalized = FALSE;

#ifdef SCIP_STATISTIC
   (*consdata)->nirrelevantjobs = 0;
   (*consdata)->nalwaysruns = 0;
   (*consdata)->ndualfixs = 0;
   (*consdata)->nremovedlocks = 0;
   (*consdata)->ndecomps = 0;
#endif

   if( nvars > 0 )
   {
      assert(vars != NULL); /* for flexelint */

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->demands, demands, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->durations, durations, nvars) );
      (*consdata)->linkingconss = NULL;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->downlocks, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->uplocks, nvars) );

      /* initialize locking arrays */
      for( v = 0; v < nvars; ++v )
      {
         (*consdata)->downlocks[v] = TRUE;
         (*consdata)->uplocks[v] = TRUE;
      }

      if( linkingconss != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->linkingconss, linkingconss, nvars) );
      }

      /* transform variables, if they are not yet transformed */
      if( SCIPisTransformed(scip) )
      {
         SCIPdebugMessage("get tranformed variables and constraints\n");

         /* get transformed variables and do NOT captures these */
         SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

         if( linkingconss != NULL )
         {
            /* get transformed constraints and captures these */
            SCIP_CALL( SCIPtransformConss(scip, (*consdata)->nvars, (*consdata)->linkingconss, (*consdata)->linkingconss) );

            for( v = 0; v < nvars; ++v )
               assert(SCIPgetConsLinking(scip, (*consdata)->vars[v]) == (*consdata)->linkingconss[v]);
         }
      }
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->downlocks = NULL;
      (*consdata)->uplocks = NULL;
      (*consdata)->demands = NULL;
      (*consdata)->durations = NULL;
      (*consdata)->linkingconss = NULL;
   }

   /* initialize values for running propagation algorithms efficiently */
   (*consdata)->resstrength1 = -1.0;
   (*consdata)->resstrength2 = -1.0;
   (*consdata)->cumfactor1 = -1.0;
   (*consdata)->disjfactor1 = -1.0;
   (*consdata)->disjfactor2 = -1.0;
   (*consdata)->estimatedstrength = -1.0;

   SCIPstatistic( (*consdata)->maxpeak = -1 );

   return SCIP_OKAY;
}

/** releases LP rows of constraint data and frees rows array */
static
SCIP_RETCODE consdataFreeRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< constraint data */
   )
{
   int r;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   for( r = 0; r < (*consdata)->ndemandrows; ++r )
   {
      assert((*consdata)->demandrows[r] != NULL);
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->demandrows[r]) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->demandrows, (*consdata)->demandrowssize);

   (*consdata)->ndemandrows = 0;
   (*consdata)->demandrowssize = 0;

   /* free rows of cover cuts */
   for( r = 0; r < (*consdata)->nscoverrows; ++r )
   {
      assert((*consdata)->scoverrows[r] != NULL);
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->scoverrows[r]) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->scoverrows, (*consdata)->scoverrowssize);

   (*consdata)->nscoverrows = 0;
   (*consdata)->scoverrowssize = 0;

   for( r = 0; r < (*consdata)->nbcoverrows; ++r )
   {
      assert((*consdata)->bcoverrows[r] != NULL);
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->bcoverrows[r]) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->bcoverrows, (*consdata)->bcoverrowssize);

   (*consdata)->nbcoverrows = 0;
   (*consdata)->bcoverrowssize = 0;

   (*consdata)->covercuts = FALSE;

   return SCIP_OKAY;
}

/** frees a cumulative constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to linear constraint data */
   )
{
   int varssize;
   int nvars;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   nvars =  (*consdata)->nvars;
   varssize = (*consdata)->varssize;

   if( varssize > 0 )
   {
      int v;

      /* release and free the rows */
      SCIP_CALL( consdataFreeRows(scip, consdata) );

      /* release the linking constraints if they were generated */
      if( (*consdata)->linkingconss != NULL )
      {
         for( v = nvars-1; v >= 0; --v )
         {
            assert((*consdata)->linkingconss[v] != NULL );
            SCIP_CALL( SCIPreleaseCons(scip, &(*consdata)->linkingconss[v]) );
         }

         SCIPfreeBlockMemoryArray(scip, &(*consdata)->linkingconss, varssize);
      }

      /* free arrays */
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->downlocks, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->uplocks, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->durations, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->demands, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, varssize);
   }

   /* free memory */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints cumulative constraint to file stream */
static
void consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< cumulative constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);

   /* print coefficients */
   SCIPinfoMessage( scip, file, "cumulative(");

   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      if( v > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "<%s>(%d)[%d]", SCIPvarGetName(consdata->vars[v]),
         consdata->durations[v], consdata->demands[v]);
   }
   SCIPinfoMessage(scip, file, ")[%d,%d) <= %d", consdata->hmin, consdata->hmax, consdata->capacity);
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
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   assert(pos < consdata->nvars);

   SCIPdebugMessage("cumulative constraint <%s>: remove variable <%s>\n",
      SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[pos]));

   /* remove the rounding locks for the deleted variable */
   SCIP_CALL( SCIPunlockVarCons(scip, consdata->vars[pos], cons, consdata->downlocks[pos], consdata->uplocks[pos]) );

   consdata->downlocks[pos] = FALSE;
   consdata->uplocks[pos] = FALSE;

   if( consdata->linkingconss != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &consdata->linkingconss[pos]) );
   }

   /* get event handler */
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* drop events */
   SCIP_CALL( consdataDropEvents(scip, consdata, conshdlrdata->eventhdlr, pos) );

   SCIPdebugMessage("remove variable <%s>[%g,%g] from cumulative constraint <%s>\n",
      SCIPvarGetName(consdata->vars[pos]), SCIPvarGetLbGlobal(consdata->vars[pos]), SCIPvarGetUbGlobal(consdata->vars[pos]), SCIPconsGetName(cons));


   /* in case the we did not remove the variable in the last slot of the arrays we move the current last to this
    * position
    */
   if( pos != consdata->nvars - 1 )
   {
      consdata->vars[pos] = consdata->vars[consdata->nvars-1];
      consdata->downlocks[pos] = consdata->downlocks[consdata->nvars-1];
      consdata->uplocks[pos] = consdata->uplocks[consdata->nvars-1];
      consdata->demands[pos] = consdata->demands[consdata->nvars-1];
      consdata->durations[pos] = consdata->durations[consdata->nvars-1];

      if( consdata->linkingconss != NULL )
      {
         consdata->linkingconss[pos]= consdata->linkingconss[consdata->nvars-1];
      }
   }

   consdata->nvars--;
   consdata->normalized = FALSE;

   return SCIP_OKAY;
}

/** collect linking constraints for each integer variable */
static
SCIP_RETCODE consdataCollectLinkingCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< pointer to consdata */
   )
{
   int nvars;
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   assert(nvars > 0);
   assert(consdata->linkingconss == NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->linkingconss, consdata->varssize) );

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CONS* cons;
      SCIP_VAR* var;

      var = consdata->vars[v];
      assert(var != NULL);

      SCIPdebugMessage("linking constraint (%d of %d) for variable <%s>\n", v+1, nvars, SCIPvarGetName(var));

      /* create linking constraint if it does not exist yet */
      if( !SCIPexistsConsLinking(scip, var) )
      {
         char name[SCIP_MAXSTRLEN];

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "link(%s)", SCIPvarGetName(var));

         /** creates and captures an linking constraint */
         SCIP_CALL( SCIPcreateConsLinking(scip, &cons, name, var, NULL, 0, 0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE /*TRUE*/, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         consdata->linkingconss[v] = cons;
      }
      else
      {
         consdata->linkingconss[v] = SCIPgetConsLinking(scip, var);
         SCIP_CALL( SCIPcaptureCons(scip, consdata->linkingconss[v]) );
      }

      assert(SCIPexistsConsLinking(scip, var));
      assert(consdata->linkingconss[v] != NULL);
      assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(consdata->linkingconss[v])), "linking") == 0 );
      assert(SCIPgetConsLinking(scip, var) == consdata->linkingconss[v]);
   }

   return SCIP_OKAY;
}

/**@} */

/**@name Sorting methods
 *
 * @{
 */

/** comparison method for two variables w.r.t. the lower bounds (earliest start time) */
static
SCIP_DECL_SORTPTRCOMP(compVarsEst)
{
   int est1;
   int est2;

   est1 = (int)(SCIPvarGetLbLocal((SCIP_VAR*)elem1) + 0.5);
   est2 = (int)(SCIPvarGetLbLocal((SCIP_VAR*)elem2) + 0.5);

   return (est1 - est2);
}

/** check if the variables are sorted in a non-increasing w.r.t. the earliest start time */
#ifndef NDEBUG
static
void checkSortVariablesEst(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   nvars               /**< number of variables to check */
   )
{
   int i;
   int lb1;
   int lb2;

   for( i = 0; i < nvars-1; ++i )
   {
      lb1 = convertBoundToInt(scip, SCIPvarGetLbGlobal(consdata->vars[i]));
      lb2 = convertBoundToInt(scip, SCIPvarGetLbGlobal(consdata->vars[i+1]));
      assert(lb1 >= lb2);
   }
}
#endif

/** sorts variables in non-increasing order w.r.t. their earliest start times */
static
SCIP_RETCODE sortVariablesEst(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< cumulative constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(consdata->nvars == 0 || consdata->durations != NULL);

   /* sort array w.r.t. earliest start time */
   SCIPsortDownPtrIntIntBoolBool((void**)consdata->vars, consdata->durations, consdata->demands, consdata->downlocks, consdata->uplocks,
      compVarsEst, consdata->nvars);

#ifndef NDEBUG
   checkSortVariablesEst(scip, consdata, consdata->nvars);
#endif

   return SCIP_OKAY;
}

/** check if the variables are sorted in a non-decreasing w.r.t. the latest completion time */
#ifndef NDEBUG
static
void checkSortVariablesLct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   nvars               /**< number of variables to check */
   )
{
   int i;
   int ub1;
   int ub2;

   for( i = 0; i < nvars-1; ++i )
   {
      ub1 = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[i])) + consdata->durations[i];
      ub2 = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[i+1]))  + consdata->durations[i+1];
      assert(ub1 <= ub2);
   }
}
#endif

/** sorts variables in non-decreasing order w.r.t. their latest completion time */
static
SCIP_RETCODE sortVariablesLct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< cumulative constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int* lcts;
   int nvars;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(consdata->nvars == 0 || consdata->durations != NULL);

   nvars = consdata->nvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &lcts, nvars) );

   for( v = 0; v < nvars; ++v )
      lcts[v] = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[v])) + consdata->durations[v];

   /* sort of three joint arrays of ints/ints/Ptr, sorted by first array in non-decreasing order via sort template */
   SCIPsortIntPtrIntIntBoolBool(lcts, (void**)consdata->vars, consdata->durations, consdata->demands, consdata->downlocks, consdata->uplocks, nvars);

#ifndef NDEBUG
   checkSortVariablesLct(scip, consdata, consdata->nvars);
#endif

   SCIPfreeBufferArray(scip, &lcts);

   return SCIP_OKAY;
}

/**@} */

/**@name Check methods
 *
 * @{
 */

/** check for the given starting time variables with their demands and durations if the cumulative conditions for the
 *  given solution is satisfied
 */
static
SCIP_RETCODE checkCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_Bool*            violated,           /**< pointer to store if the cumulative condition is violated */
   SCIP_CONS*            cons,               /**< constraint which is checked */
   SCIP_Bool             printreason         /**< should the reason for the violation be printed? */
   )
{
   int* startsolvalues;       /* stores when each job is starting */
   int* endsolvalues;         /* stores when each job ends */
   int* startindices;         /* we will sort the startsolvalues, thus we need to know which index of a job it corresponds to */
   int* endindices;           /* we will sort the endsolvalues, thus we need to know which index of a job it corresponds to */

   int freecapacity;
   int curtime;            /* point in time which we are just checking */
   int endindex;           /* index of endsolvalues with: endsolvalues[endindex] > curtime */
   int j;

   assert(scip != NULL);
   assert(violated != NULL);

   (*violated) = FALSE;

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);
   assert(demands != NULL);
   assert(durations != NULL);

   /* compute time points where we have to check whether capacity constraint is infeasible or not */
   SCIP_CALL( SCIPallocBufferArray(scip, &startsolvalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endsolvalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   /* assign variables, start and endpoints to arrays */
   for ( j = 0; j < nvars; ++j )
   {
      /* the constraint of the cumulative constraint handler should be called after the integrality check */
      assert(SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[j])));

      startsolvalues[j] = convertBoundToInt(scip, SCIPgetSolVal(scip, sol, vars[j]));
      startindices[j] = j;

      endsolvalues[j] = startsolvalues[j] + durations[j];
      endindices[j] = j;
   }

   /* sort the arrays not-decreasing according to start solution values and end solution values (and sort the
    * corresponding indices in the same way)
    */
   SCIPsortIntInt(startsolvalues, startindices, nvars);
   SCIPsortIntInt(endsolvalues, endindices, nvars);

   endindex = 0;
   freecapacity = capacity;

   /* check each start point of a job whether the capacity is kept or not */
   for( j = 0; j < nvars; ++j )
   {
      /* only check intervals [hmin,hmax) */
      curtime = startsolvalues[j];

      if( curtime >= hmax )
         break;

      /* subtract all capacity needed up to this point */
      freecapacity -= demands[startindices[j]];
      while( j+1 < nvars && startsolvalues[j+1] == curtime )
      {
         j++;
         freecapacity -= demands[startindices[j]];
      }

      /* free all capacity usages of jobs that are no longer running */
      while( endindex < nvars && curtime >= endsolvalues[endindex] )
      {
         freecapacity += demands[endindices[endindex]];
         ++endindex;
      }
      assert(freecapacity <= capacity);

      /* check freecapacity to be smaller than zero */
      if( freecapacity < 0 && curtime >= hmin )
      {
         SCIPdebugMessage("freecapacity = %3d \n", freecapacity);
         (*violated) = TRUE;

         if( printreason )
         {
            int i;

            /* first state the violated constraints */
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );

            /* second state the reason */
            SCIPinfoMessage(scip, NULL,
               ";\nviolation: at time point %d available capacity = %d, needed capacity = %d\n",
               curtime, capacity, capacity - freecapacity);

            for( i = 0; i <= j; ++i )
            {
               if( startsolvalues[i] + durations[startindices[i]] > curtime )
               {
                  SCIPinfoMessage(scip, NULL, "activity %s, start = %i, duration = %d, demand = %d \n",
                     SCIPvarGetName(vars[startindices[i]]), startsolvalues[i], durations[startindices[i]],
                     demands[startindices[i]]);
               }
            }
         }
         break;
      }
   } /*lint --e{850}*/

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endsolvalues);
   SCIPfreeBufferArray(scip, &startsolvalues);

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

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violated != NULL);

   SCIPdebugMessage("check cumulative constraints <%s>\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check the cumulative condition */
   SCIP_CALL( checkCumulativeCondition(scip, sol, consdata->nvars, consdata->vars,
         consdata->durations, consdata->demands, consdata->capacity, consdata->hmin, consdata->hmax,
         violated, cons, printreason) );

   return SCIP_OKAY;
}

/**@} */

/**@name Conflict analysis
 *
 * @{
 */

/** relax lower bound of give variable as long as the given inference lower bound is reached and add that bound change
 *  to the conflict set
 */
static
SCIP_RETCODE relaxVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the upper bound should be relaxed */
   int                   inferlb,            /**< smallest lower bound which still leads to the propagation */
   SCIP_BDCHGIDX*        bdchgidx            /**< the index of the bound change, representing the point of time where the change took place */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_Real conflictlb;
   SCIP_Real relaxedlb;
   SCIP_Real newlb;
   int nbdchgs;

   /* get number of bound changes */
   nbdchgs = SCIPvarGetNBdchgInfosLb(var);

   /* get current conflict lower bound */
   conflictlb = SCIPgetConflictVarLb(scip, var);
   relaxedlb = SCIPgetConflictVarRelaxedLb(scip, var);

   assert(nbdchgs > 0 || SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetLbGlobal(var)));
   assert(nbdchgs == 0 || SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPbdchginfoGetNewbound(SCIPvarGetBdchgInfoLb(var, nbdchgs-1))));

   /* first we move in the lower bound change array to given bdchgidx */
   while( nbdchgs > 0 && SCIPbdchgidxIsEarlier(bdchgidx, SCIPbdchginfoGetIdx(SCIPvarGetBdchgInfoLb(var, nbdchgs-1))) )
      nbdchgs--;

   SCIPdebugMessage("variable <%s>[%.15g,%.15g]: nbdchgs %d try to relax lower bound to %d (conflict lower bound %g)\n",
      SCIPvarGetName(var), SCIPvarGetLbAtIndex(var, bdchgidx, FALSE), SCIPvarGetUbAtIndex(var, bdchgidx, FALSE), nbdchgs, inferlb, conflictlb);

   newlb = SCIPvarGetLbAtIndex(var, bdchgidx, FALSE);

   /* if the current conflict lower bound is worse or equal (greater than or equal) than the new lowerbound we can
    * return since the a stronger bound is already part of the conflict (widening does not help)
    */
   if( conflictlb - newlb > 0.5 )
      return SCIP_OKAY;

   /* try to relax lower bound; we can stop if we reached the current conflict lower bound */
   while( nbdchgs > 0 && newlb > conflictlb )
   {
      bdchginfo = SCIPvarGetBdchgInfoLb(var, nbdchgs-1);

      SCIPdebugMessage("lower bound change %d oldbd=%.15g, newbd=%.15g, depth=%d, pos=%d, redundant=%u\n",
         nbdchgs, SCIPbdchginfoGetOldbound(bdchginfo), SCIPbdchginfoGetNewbound(bdchginfo),
         SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo), SCIPbdchginfoIsRedundant(bdchginfo));

      /* check if the old upper bound is sufficient for propagation (creationg a core) */
      if( inferlb - SCIPbdchginfoGetOldbound(bdchginfo) > 0.5 )
         break;

      SCIPdebugMessage("***** relaxed lower bound of inference variable <%s> from <%g> to <%g>\n",
         SCIPvarGetName(var), SCIPbdchginfoGetNewbound(bdchginfo), SCIPbdchginfoGetOldbound(bdchginfo));

      bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
      newlb = SCIPbdchginfoGetOldbound(bdchginfo);
      nbdchgs--;
   }

   /* adjust relaxed lower bound w.r.t. already known one */
   relaxedlb = MAX(relaxedlb, (SCIP_Real)inferlb);

   /* if the nbdchgs is zero then the local bound matches the global bound, therefore bdchgidx equal to NULL represents
    * the right time point and SCIP finds out that this bound is redundant since it is global
    */
   SCIPdebugMessage("add lower bound of bound change info %d to conflict set\n", nbdchgs);
   SCIP_CALL( SCIPaddConflictRelaxedLb(scip, var, bdchgidx, relaxedlb) );

   return SCIP_OKAY;
}

/** relax upper bound of give variable as long as the given inference upper bound is reached and add that bound change
 *  to the conflict set
 */
static
SCIP_RETCODE relaxVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the upper bound should be relaxed */
   int                   inferub,            /**< largest upper bound which still leads to the propagation */
   SCIP_BDCHGIDX*        bdchgidx            /**< the index of the bound change, representing the point of time where the change took place */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_Real conflictub;
   SCIP_Real relaxedub;
   SCIP_Real newub;
   int nbdchgs;

   assert(SCIPvarIsActive(var));

   /* get number of bound changes */
   nbdchgs = SCIPvarGetNBdchgInfosUb(var);

   /* get current conflict lower bound */
   conflictub = SCIPgetConflictVarUb(scip, var);
   relaxedub = SCIPgetConflictVarRelaxedUb(scip, var);

   assert(nbdchgs > 0 || SCIPisEQ(scip, SCIPvarGetUbLocal(var), SCIPvarGetUbGlobal(var)));
   assert(nbdchgs == 0 || SCIPisEQ(scip, SCIPvarGetUbLocal(var), SCIPbdchginfoGetNewbound(SCIPvarGetBdchgInfoUb(var, nbdchgs-1))));

   /* first we move in the upper bound change array to given bdchgidx */
   while( nbdchgs > 0 && SCIPbdchgidxIsEarlier(bdchgidx, SCIPbdchginfoGetIdx(SCIPvarGetBdchgInfoUb(var, nbdchgs-1))) )
      nbdchgs--;

   SCIPdebugMessage("variable <%s>[%.15g,%.15g]: nbdchgs %d try to relax upper bound to %d (conflict upper bound %g)\n",
      SCIPvarGetName(var), SCIPvarGetLbAtIndex(var, bdchgidx, FALSE), SCIPvarGetUbAtIndex(var, bdchgidx, FALSE), nbdchgs, inferub, conflictub);

   newub = SCIPvarGetUbAtIndex(var, bdchgidx, FALSE);

   /* if the current conflict upper bound is worse or equal (smaller than  or equal) than the new upper bound we can return
    * since the a stronger bound is already part of the conflict (widening does not help)
    */
   if( newub - conflictub > 0.5 )
      return SCIP_OKAY;

   /* try to relax upperbound; we can stop if we reached the current conflict upper bound */
   while( nbdchgs > 0 && newub < conflictub )
   {
      bdchginfo = SCIPvarGetBdchgInfoUb(var, nbdchgs-1);

      SCIPdebugMessage("upper bound change %d oldbd=%.15g, newbd=%.15g, depth=%d, pos=%d, redundant=%u\n",
         nbdchgs, SCIPbdchginfoGetOldbound(bdchginfo), SCIPbdchginfoGetNewbound(bdchginfo),
         SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo), SCIPbdchginfoIsRedundant(bdchginfo));

      /* check if the old upper bound is sufficient for propagation (creationg a core) */
      if( SCIPbdchginfoGetOldbound(bdchginfo) - inferub > 0.5 )
         break;

      SCIPdebugMessage("***** relaxed lower bound of inference variable <%s> from <%g> to <%g>\n",
         SCIPvarGetName(var), SCIPbdchginfoGetNewbound(bdchginfo), SCIPbdchginfoGetOldbound(bdchginfo));

      bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
      newub = SCIPbdchginfoGetOldbound(bdchginfo);
      nbdchgs--;
   }

   /* adjust relaxed upper bound w.r.t. already known one */
   relaxedub = MIN(relaxedub, (SCIP_Real)inferub);

   /* if the nbdchgs is zero then the local bound matches the global bound, therefore bdchgidx equal to NULL represents
    * the right time point and SCIP finds out that this bound is redundant since it is global
    */
   SCIPdebugMessage("add lower bound of bound change info %d to conflict set\n", nbdchgs);
   SCIP_CALL( SCIPaddConflictRelaxedUb(scip, var, bdchgidx, relaxedub) );

   return SCIP_OKAY;
}
/** resolves the propagation of the core time algorithm */
static
SCIP_RETCODE resolvePropagationCoretimes(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             infervar,           /**< inference variable */
   int                   inferdemand,        /**< demand of the inference variable */
   int                   inferpeak,          /**< time point which causes the propagation */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Bool*            explanation         /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   )
{
   SCIP_VAR* var;
   SCIP_Bool* reported;
   int duration;
   int ect;
   int lst;
   int j;

   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPinProbing(scip));

   SCIPdebugMessage("variable <%s>: (demand %d) resolve propagation of core time algorithm (peak %d)\n",
      SCIPvarGetName(infervar), inferdemand, inferpeak);

   assert(nvars > 0);

   capacity -= inferdemand;

   SCIP_CALL( SCIPallocBufferArray(scip, &reported, nvars) );
   BMSclearMemoryArray(reported, nvars);

   /* first we loop over all variable adjust the capacity with those which provide a global core at the inference peak
    * and those where the current conflict bounds provide a core at the inference peak
    */
   for( j = 0; j < nvars && capacity >= 0; ++j )
   {
      var = vars[j];
      assert(var != NULL);

      /* skip inference variable */
      if( var == infervar )
         continue;

      duration = durations[j];
      assert(duration > 0);

      /* compute cores of jobs; if core overlaps interval of inference variable add this job to the array */
      assert(SCIPisFeasEQ(scip, SCIPvarGetUbAtIndex(var, bdchgidx, TRUE), SCIPvarGetUbAtIndex(var, bdchgidx, FALSE)));
      assert(SCIPisFeasIntegral(scip, SCIPvarGetUbAtIndex(var, bdchgidx, TRUE)));
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbAtIndex(var, bdchgidx, TRUE), SCIPvarGetLbAtIndex(var, bdchgidx, FALSE)));
      assert(SCIPisFeasIntegral(scip, SCIPvarGetLbAtIndex(var, bdchgidx, TRUE)));

      SCIPdebugMessage("variable <%s>: glb=[%g,%g] conflict=[%g,%g] (duration %d, demand %d)\n",
         SCIPvarGetName(var),  SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
         SCIPgetConflictVarLb(scip, var), SCIPgetConflictVarUb(scip, var), duration, demands[j]);

      /* collect the conflict bound core (the conflict bounds are those bounds which are already part of the conflict)
       * hence these bound are already reported by other resolve propation steps. In case a bound (lower or upper) is
       * not part of the conflict yet we get the global bounds back.
       */
      ect = convertBoundToInt(scip, SCIPgetConflictVarRelaxedLb(scip, var)) + duration;
      lst = convertBoundToInt(scip, SCIPgetConflictVarRelaxedUb(scip, var));

      /* check if the inference peak is part of the global/conflict bound core; if so we decreasing the capacity by the
       * demand of that job without adding anything to the explanation
       */
      if( inferpeak < ect && lst <= inferpeak )
      {
         capacity -= demands[j];
         reported[j] = TRUE;

         if( explanation != NULL )
            explanation[j] = TRUE;
      }
   }

   /* collect all cores of the variables which lay in the considered time window except the inference variable */
   for( j = 0; j < nvars && capacity >= 0; ++j )
   {
      var = vars[j];
      assert(var != NULL);

      /* skip inference variable */
      if( var == infervar || reported[j] )
         continue;

      duration = durations[j];
      assert(duration > 0);

      /* compute cores of jobs; if core overlaps interval of inference variable add this job to the array */
      assert(SCIPisFeasEQ(scip, SCIPvarGetUbAtIndex(var, bdchgidx, TRUE), SCIPvarGetUbAtIndex(var, bdchgidx, FALSE)));
      assert(SCIPisFeasIntegral(scip, SCIPvarGetUbAtIndex(var, bdchgidx, TRUE)));
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbAtIndex(var, bdchgidx, TRUE), SCIPvarGetLbAtIndex(var, bdchgidx, FALSE)));
      assert(SCIPisFeasIntegral(scip, SCIPvarGetLbAtIndex(var, bdchgidx, TRUE)));

      /* collect local core information */
      ect = convertBoundToInt(scip, SCIPvarGetLbAtIndex(var, bdchgidx, FALSE)) + duration;
      lst = convertBoundToInt(scip, SCIPvarGetUbAtIndex(var, bdchgidx, FALSE));

      SCIPdebugMessage("variable <%s>: loc=[%g,%g] glb=[%g,%g] (duration %d, demand %d)\n",
         SCIPvarGetName(var), SCIPvarGetLbAtIndex(var, bdchgidx, FALSE), SCIPvarGetUbAtIndex(var, bdchgidx, FALSE),
         SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), duration, demands[j]);

      /* check if the inference peak is part of the core */
      if( inferpeak < ect && lst <= inferpeak )
      {
         int aggrpeak;

         /* check the current status of the variable in */
         SCIPdebugMessage("variable <%s> (durations %d demands %d) loc=[%g,%g] glb=[%g,%g] conflict=[%g,%g]\n",
            SCIPvarGetName(var), duration, demands[j],
            SCIPvarGetLbAtIndex (var, bdchgidx, FALSE), SCIPvarGetUbAtIndex(var, bdchgidx, FALSE),
            SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
            SCIPgetConflictVarLb(scip, var), SCIPgetConflictVarUb(scip, var));

         aggrpeak = inferpeak;

         if( !SCIPvarIsActive(var) )
         {
            SCIP_Real scalar;
            SCIP_Real constant;

            scalar = 1.0;
            constant = 0.0;

            SCIP_CALL( SCIPgetProbvarSum(scip, &var, &scalar, &constant) );
            assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN);
            assert(!SCIPisZero(scip, scalar));

            /* compute inference peak w.r.t. aggregation */
            aggrpeak = convertBoundToInt(scip, (aggrpeak - constant) / scalar);
         }

         SCIP_CALL( relaxVarLb(scip, var, aggrpeak - duration + 1, bdchgidx) );
         SCIP_CALL( relaxVarUb(scip, var, aggrpeak, bdchgidx) );

         capacity -= demands[j];

         if( explanation != NULL )
            explanation[j] = TRUE;
      }
   }

   assert(capacity < 0);

   SCIPfreeBufferArray(scip, &reported);

   return SCIP_OKAY;
}

/** resolve propagation w.r.t. the cumulative condition */
static
SCIP_RETCODE respropCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             infervar,           /**< the conflict variable whose bound change has to be resolved */
   INFERINFO             inferinfo,          /**< the user information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Bool*            explanation,        /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   int inferdemand;
   int inferduration;
   int inferpos;

   *result = SCIP_DIDNOTFIND;

   /* get the position of the inferred variable in the vars array */
   inferpos = inferInfoGetData1(inferinfo);
   if( inferpos >= nvars || vars[inferpos] != infervar )
   {
      /* find inference variable in constraint */
      for( inferpos = 0; inferpos < nvars && vars[inferpos] != infervar; ++inferpos )
      {}
   }
   assert(inferpos < nvars);
   assert(vars[inferpos] == infervar);

   inferdemand = demands[inferpos];
   inferduration = durations[inferpos];

   SCIPdebugMessage("variable <%s>: duration %d, demand %d\n", SCIPvarGetName(infervar), inferduration, inferdemand);

   switch( inferInfoGetProprule(inferinfo) )
   {
   case PROPRULE_1_CORETIMES:
   {
      int inferpeak;

      if( boundtype == SCIP_BOUNDTYPE_UPPER )
      {
         /* we propagated the latest start time (upper bound) step wise with a step length of at most the duration of
          *  the inference variable
          */
         assert(SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE) - SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < inferduration + 0.5);

         SCIPdebugMessage("variable <%s>: upper bound changed from %g to %g\n",
            SCIPvarGetName(infervar), SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE),
            SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE));

         /* get the inference peak that the time point which lead to the that propagtion */
         inferpeak = convertBoundToInt(scip, SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE)) + inferduration;

         /* old upper bound of variable itself is part of the explanation */
         SCIP_CALL( SCIPaddConflictUb(scip, infervar, bdchgidx) );
      }
      else
      {
         assert(boundtype == SCIP_BOUNDTYPE_LOWER);

         SCIPdebugMessage("variable <%s>: lower bound changed from %g to %g\n",
            SCIPvarGetName(infervar), SCIPvarGetLbAtIndex(infervar, bdchgidx, FALSE),
            SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE));

         /* get the time interval where the job could not be scheduled */
         inferpeak = convertBoundToInt(scip, SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE)) - 1;

         /* old lower bound of variable itself is part of the explanation */
         SCIP_CALL( SCIPaddConflictLb(scip, infervar, bdchgidx) );
      }

      SCIP_CALL( resolvePropagationCoretimes(scip, nvars, vars, durations, demands, capacity,
            infervar, inferdemand, inferpeak, bdchgidx, explanation) );

      (*result) = SCIP_SUCCESS;

      break;
   }

   case PROPRULE_2_EDGEFINDING:
      break;

   case PROPRULE_3_ENERGETICREASONING:
      break;
   }

   return SCIP_OKAY;
}

/**@} */


/**@name Enforcement methods
 *
 * @{
 */

/** apply all fixings which are given by the alternative bounds */
static
SCIP_RETCODE applyAlternativeBoundsBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of active variables */
   int                   nvars,              /**< number of active variables */
   int*                  alternativelbs,     /**< alternative lower bounds */
   int*                  alternativeubs,     /**< alternative lower bounds */
   int*                  downlocks,          /**< number of constraints with down lock participating by the computation */
   int*                  uplocks,            /**< number of constraints with up lock participating by the computation */
   SCIP_Bool*            branched            /**< pointer to store if a branching was applied */
   )
{
   int v;

   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;

      var = vars[v];
      assert(var != NULL);


      if( SCIPvarGetNLocksDown(var) == downlocks[v] && SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_LOWER )
      {
         int ub;

         ub = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

         if( alternativelbs[v] <= ub )
         {
            SCIP_CALL( SCIPbranchVarHole(scip, var, SCIPvarGetLbLocal(var), (SCIP_Real)alternativelbs[v], NULL, NULL) );
            (*branched) = TRUE;

            SCIPdebugMessage("variable <%s> branched domain hole (%g,%d)\n", SCIPvarGetName(var),
               SCIPvarGetLbLocal(var), alternativelbs[v]);

            return SCIP_OKAY;
         }
      }

      if( SCIPvarGetNLocksUp(var) == uplocks[v] && SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_UPPER )
      {
         int lb;

         lb = convertBoundToInt(scip, SCIPvarGetLbLocal(var));

         if( alternativeubs[v] >= lb )
         {
            SCIP_CALL( SCIPbranchVarHole(scip, var, (SCIP_Real)alternativeubs[v], SCIPvarGetUbLocal(var), NULL, NULL) );
            (*branched) = TRUE;

            SCIPdebugMessage("variable <%s> branched domain hole (%d,%g)\n", SCIPvarGetName(var),
               alternativeubs[v],  SCIPvarGetUbLocal(var));

            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** remove the capacity requirments for all job which start at the curtime */
static
void subtractStartingJobDemands(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   curtime,            /**< current point in time */
   int*                  starttimes,         /**< array of start times */
   int*                  startindices,       /**< permutation with respect to the start times */
   int*                  freecapacity,       /**< pointer to store the resulting free capacity */
   int*                  idx,                /**< pointer to index in start time array */
   int                   nvars               /**< number of vars in array of starttimes and startindices */
   )
{

#if defined SCIP_DEBUG && !defined NDEBUG
   int oldidx;
   oldidx = *idx;
#endif

   assert(idx != NULL);
   assert(starttimes != NULL);
   assert(starttimes != NULL);
   assert(freecapacity != NULL);
   assert(starttimes[*idx] == curtime);
   assert(consdata->demands != NULL);
   assert(freecapacity != idx);

   /* subtract all capacity needed up to this point */
   (*freecapacity) -= consdata->demands[startindices[*idx]];

   while( (*idx)+1 < nvars && starttimes[(*idx)+1] == curtime )
   {
      ++(*idx);
      (*freecapacity) -= consdata->demands[startindices[(*idx)]];
      assert(freecapacity != idx);
   }
#ifdef SCIP_DEBUG
   assert(oldidx <= *idx);
#endif
}

/** add the capacity requirments for all job which end at the curtime */
static
void addEndingJobDemands(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   curtime,            /**< current point in time */
   int*                  endtimes,           /**< array of end times */
   int*                  endindices,         /**< permutation with rspect to the end times */
   int*                  freecapacity,       /**< pointer to store the resulting free capacity */
   int*                  idx,                /**< pointer to index in end time array */
   int                   nvars               /**< number of vars in array of starttimes and startindices */
   )
{
#if defined SCIP_DEBUG && !defined NDEBUG
   int oldidx;
   oldidx = *idx;
#endif

   /* free all capacity usages of jobs the are no longer running */
   while( endtimes[*idx] <= curtime && *idx < nvars)
   {
      (*freecapacity) += consdata->demands[endindices[*idx]];
      ++(*idx);
   }

#ifdef SCIP_DEBUG
   assert(oldidx <= *idx);
#endif
}

/** computes a point in time when the capacity is exceeded returns hmax if this does not happen */
static
SCIP_RETCODE computePeak(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint handler data */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int*                  timepoint           /**< pointer to store the time point of the peak */
   )
{
   int* starttimes;         /* stores when each job is starting */
   int* endtimes;           /* stores when each job ends */
   int* startindices;       /* we will sort the startsolvalues, thus we need to know wich index of a job it corresponds to */
   int* endindices;         /* we will sort the endsolvalues, thus we need to know wich index of a job it corresponds to */

   int nvars;               /* number of activities for this constraint */
   int freecapacity;        /* remaining capacity */
   int curtime;             /* point in time which we are just checking */
   int endindex;            /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int hmin;
   int hmax;

   int j;

   assert(consdata != NULL);

   nvars = consdata->nvars;
   assert(nvars > 0);

   *timepoint = consdata->hmax;

   assert(consdata->vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   /* create event point arrays */
   createSortedEventpointsSol(scip, sol, consdata->nvars, consdata->vars, consdata->durations,
      starttimes, endtimes, startindices, endindices);

   endindex = 0;
   freecapacity = consdata->capacity;
   hmin = consdata->hmin;
   hmax = consdata->hmax;

   /* check each startpoint of a job whether the capacity is kept or not */
   for( j = 0; j < nvars; ++j )
   {
      curtime = starttimes[j];
      SCIPdebugMessage("look at %d-th job with start %d\n", j, curtime);

      if( curtime >= hmax )
         break;

      /* remove the capacity requirments for all job which start at the curtime */
      subtractStartingJobDemands(consdata, curtime, starttimes, startindices, &freecapacity, &j, nvars);

      /* add the capacity requirments for all job which end at the curtime */
      addEndingJobDemands(consdata, curtime, endtimes, endindices, &freecapacity, &endindex, nvars);

      assert(freecapacity <= consdata->capacity);
      assert(endindex <= nvars);

      /* endindex - points to the next job which will finish */
      /* j - points to the last job that has been released */

      /* if free capacity is smaller than zero, then add branching candidates */
      if( freecapacity < 0 && curtime >= hmin )
      {
         *timepoint = curtime;
         break;
      }
   } /*lint --e{850}*/

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endtimes);
   SCIPfreeBufferArray(scip, &starttimes);

   return SCIP_OKAY;
}

/** checks all cumulative constraints for infeasibility and add branching candidates to storage */
static
SCIP_RETCODE collectBranchingCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to be processed */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int*                  nbranchcands        /**< pointer to store the number of branching variables */
   )
{
   SCIP_HASHTABLE* collectedvars;
   int c;

   assert(scip != NULL);
   assert(conss != NULL);

   /* create a hash table */
   SCIP_CALL( SCIPhashtableCreate(&collectedvars, SCIPblkmem(scip), SCIPcalcHashtableSize(SCIPgetNVars(scip)),
         SCIPvarGetHashkey, SCIPvarIsHashkeyEq, SCIPvarGetHashkeyVal, NULL) );

   assert(scip != NULL);
   assert(conss != NULL);

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;

      int curtime;
      int j;

      cons = conss[c];
      assert(cons != NULL);

      if( !SCIPconsIsActive(cons) )
         continue;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* get point in time when capacity is exceeded */
      SCIP_CALL( computePeak(scip, consdata, sol, &curtime) );

      if( curtime < consdata->hmin || curtime >= consdata->hmax )
         continue;

      /* report all variables that are running at that point in time */
      for( j = 0; j < consdata->nvars; ++j )
      {
         SCIP_VAR* var;
         int lb;
         int ub;

         var = consdata->vars[j];
         assert(var != NULL);

         /* check if the variable was already added */
         if( SCIPhashtableExists(collectedvars, (void*)var) )
            continue;

         lb = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
         ub = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

         if( lb <= curtime && ub + consdata->durations[j] > curtime && lb < ub  )
         {
            SCIP_Real solval;
            SCIP_Real score;

            solval = SCIPgetSolVal(scip, sol, var);
            score = MIN(solval - lb, ub - solval) / ((SCIP_Real)ub-lb);

            SCIPdebugMessage("add var <%s> to branch cand storage\n", SCIPvarGetName(var));
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, score, lb + (ub - lb) / 2.0 + 0.2) );
            (*nbranchcands)++;

            SCIP_CALL( SCIPhashtableInsert(collectedvars, var) );
         }
      }
   }

   SCIPhashtableFree(&collectedvars);

   SCIPdebugMessage("found %d branching candidates\n", *nbranchcands);

   return SCIP_OKAY;
}

/** enforcement pseudo or LP solution */
static
SCIP_RETCODE enforceSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to be processed */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             branch,             /**< should branching candidates be collected */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   if( branch )
   {
      int nbranchcands;

      nbranchcands = 0;
      SCIP_CALL( collectBranchingCands(scip, conss, nconss, NULL, &nbranchcands) );

      if( nbranchcands > 0 )
         (*result) = SCIP_INFEASIBLE;
   }
   else
   {
      SCIP_Bool violated;
      int c;

      violated = FALSE;

      /* first check if a constraints is violated */
      for( c = 0; c < nconss && !violated; ++c )
      {
         SCIP_CONS* cons;

         cons = conss[c];
         assert(cons != NULL);

         SCIP_CALL( checkCons(scip, cons, NULL, &violated, FALSE) );
      }

      if( violated )
         (*result) = SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}

/**@} */

/**@name Propagation
 *
 * @{
 */

/** check if cumulative constraint is independently of all other constraints */
static
SCIP_Bool isConsIndependently(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< cumulative constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Bool* downlocks;
   SCIP_Bool* uplocks;
   int nvars;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   vars = consdata->vars;
   downlocks = consdata->downlocks;
   uplocks = consdata->uplocks;

   /* check if the cumulative constraint has the only locks on the involved variables */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;

      var = vars[v];
      assert(var != NULL);

      if( SCIPvarGetNLocksDown(var) > (int)downlocks[v] || SCIPvarGetNLocksUp(var) > (int)uplocks[v] )
         return FALSE;
   }

   return TRUE;
}

/** in case the cumulative constraint is independent of every else, solve the cumulative problem and apply the fixings
 *  (dual reductions)
 */
static
SCIP_RETCODE solveIndependentCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   SCIP_Longint          maxnodes,           /**< number of branch-and-bound nodes to solve an independent cumulative constraint  (-1: no limit) */
   int*                  nchgbds,            /**< pointer to store the number changed variable bounds */
   int*                  nfixedvars,         /**< pointer to count number of fixings */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints  */
   SCIP_Bool*            cutoff,             /**< pointer to store if the constraint is infeasible */
   SCIP_Bool*            unbounded           /**< pointer to store if the constraint is unbounded */
   )
{
   SCIP* subscip;
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP* varmapfw;
   SCIP_CONS* targetcons;
   SCIP_VAR** vars;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Bool succeed;
   char probname[SCIP_MAXSTRLEN];
   int nvars;
   int v;

   assert(!SCIPconsIsModifiable(cons));
   assert(SCIPgetNConss(scip) > 0);

   /* if the cumulative constraint is the only constraint do nothing */
   if( SCIPgetNConss(scip) == 1 )
      return SCIP_OKAY;

   /* constraints for which the check flag is set to FALSE, did not contribute to the lock numbers; therefore, we cannot
    * use the locks to decide for a dual reduction using this constraint;
    */
   if( !SCIPconsIsChecked(cons) )
      return SCIP_OKAY;

   /* check if constraint is independently */
   if( !isConsIndependently(scip, cons) )
      return SCIP_OKAY;

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

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( timelimit <= 0.0 || memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   vars = consdata->vars;

   SCIPdebugMessage("the cumulative constraint <%s> is independent from rest of the problem\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

   /* copy all plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* get name of the original problem and add the string "_cumulative" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_cumulative", SCIPgetProbName(scip));

   /* create the subproblem */
   SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* copy cumulative constraint */
   SCIP_CALL( SCIPgetConsCopy(scip, subscip, cons, &targetcons, SCIPconsGetHdlr(cons), varmapfw, NULL, NULL,
         FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, &succeed) );

   if( succeed )
   {
      /* add constraint to subscip */
      SCIP_CALL( SCIPaddCons(subscip, targetcons) );

      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

      /* set limits for the subproblem */
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", maxnodes) );
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

      /* forbid recursive call of heuristics and separators solving subMIPs */
      SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

      /* solve single cumulative constraint by branch and bound */
      SCIP_CALL( SCIPsolve(subscip) );

      /* evaluated solution status */
      switch( SCIPgetStatus(subscip) )
      {
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_INFEASIBLE:
         *cutoff = TRUE;
         break;
      case SCIP_STATUS_UNBOUNDED:
         *unbounded = TRUE;
         break;
      case SCIP_STATUS_OPTIMAL:
      {
         /* copy optimal as dual reduction into the original SCIP instance */
         SCIP_SOL* sol;

         sol = SCIPgetBestSol(subscip);

         for( v = 0; v < nvars; ++v )
         {
            SCIP_VAR* subvar;
            SCIP_VAR* var;
            SCIP_Real fixval;
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            var = vars[v];

            subvar = (SCIP_VAR*)SCIPhashmapGetImage(varmapfw, var);
            fixval =  SCIPgetSolVal(subscip, sol, subvar);

            SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
            assert(!infeasible);

            if( fixed )
               (*nfixedvars)++;
         }

         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss)++;

         break;
      }
      case SCIP_STATUS_NODELIMIT:
      case SCIP_STATUS_TIMELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      {
         SCIP_VAR* var;
         SCIP_Bool infeasible;
         SCIP_Bool tightened;

         /* transfer the bound changes */
         for( v = 0; v < nvars; ++v )
         {
            var = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[v]);

            SCIP_CALL( SCIPtightenVarLb(scip, vars[v], SCIPvarGetLbGlobal(var), TRUE, &infeasible, &tightened) );
            assert(!infeasible);

            if( tightened )
               (*nchgbds)++;

            SCIP_CALL( SCIPtightenVarUb(scip, vars[v], SCIPvarGetUbGlobal(var), TRUE, &infeasible, &tightened) );
            assert(!infeasible);

            if( tightened )
               (*nchgbds)++;

         }
         break;
      }
      case SCIP_STATUS_UNKNOWN:
      case SCIP_STATUS_USERINTERRUPT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_BESTSOLLIMIT:
      case SCIP_STATUS_TOTALNODELIMIT:
         SCIPerrorMessage("invalid status code <%d>\n", SCIPgetStatus(scip));
         return SCIP_INVALIDDATA;
      }
   }

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/** start conflict analysis to analysis the core insertion which is infeasible */
static
SCIP_RETCODE analyseInfeasibelCoreInsertion(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             infervar,           /**< start time variable which lead to the infeasibilty */
   int                   inferduration,      /**< duration of the start time variable */
   int                   inferdemand,        /**< demand of the start time variable */
   int                   inferpeak,          /**< profile preak which causes the infeasibilty */
   SCIP_Bool*            initialized,        /**< pointer to store if the conflict analysis was initialized */
   SCIP_Bool*            explanation         /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   )
{
   SCIPdebugMessage("detected infeasibility due to adding a core to the core resource profile\n");
   SCIPdebugMessage("variable <%s>[%g,%g] (demand %d)\n", SCIPvarGetName(infervar),
      SCIPvarGetLbLocal(infervar), SCIPvarGetUbLocal(infervar), inferdemand);

   /* initialize conflict analysis if conflict analysis is applicable */
   if( SCIPisConflictAnalysisApplicable(scip) )
   {
      SCIP_CALL( SCIPinitConflictAnalysis(scip) );

      SCIPdebugMessage("add lower and upper bounds of variable <%s>\n", SCIPvarGetName(infervar));

      SCIP_CALL( resolvePropagationCoretimes(scip, nvars, vars, durations, demands, capacity,
            infervar, inferdemand, inferpeak, NULL, explanation) );

      /* add both bound of the inference variable since these biuld the core which we could not inserted */
      SCIP_CALL( relaxVarLb(scip, infervar, inferpeak - inferduration + 1, NULL) );
      SCIP_CALL( relaxVarUb(scip, infervar, inferpeak, NULL) );

      *initialized = TRUE;
   }

   return SCIP_OKAY;
}

/** We are using the core resource profile which contains all core except the one of the start time variable which we
 *  want to propagate, to incease the earliest start time. This we are doing in steps of length at most the duration of
 *  the job. The reason for that is, that this makes it later easier to resolve this propagation during the conflict
 *  analysis
 */
static
SCIP_RETCODE coretimesUpdateLb(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which is propagated */
   SCIP_PROFILE*         profile,            /**< resource profile */
   int                   pos,                /**< position of the variable to propagate */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            explanation,        /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   SCIP_Bool*            infeasible          /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_VAR* var;
   int ntimepoints;
   int duration;
   int demand;
   int peak;
   int newlb;
   int est;
   int lst;

   var = vars[pos];
   assert(var != NULL);

   duration = durations[pos];
   assert(duration > 0);

   demand = demands[pos];
   assert(demand > 0);

   est = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
   lst = convertBoundToInt(scip, SCIPvarGetUbLocal(var));
   ntimepoints = SCIPprofileGetNTimepoints(profile);

   /* first we find left position of earliest start time (lower bound) in resource profile; this position gives us the
    * load which we have at the earliest start time (lower bound)
    */
   (void) SCIPprofileFindLeft(profile, est, &pos);

   SCIPdebugMessage("propagate earliest start time (lower bound) (pos %d)\n", pos);

   /* we now trying to move the earliest start time in steps of at most "duration" length */
   do
   {
      INFERINFO inferinfo;
      SCIP_Bool tightened;
      int ect;

#ifndef NDEBUG
      {
         /* in debug mode we check that we adjust the search position correctly */
         int tmppos;

         (void)SCIPprofileFindLeft(profile, est, &tmppos);
         assert(pos == tmppos);
      }
#endif
      ect = est + duration;
      peak = -1;

      /* we search for a peak within the core profile which conflicts with the demand of the start time variable; we
       * want a peak which is closest to the earliest completion time
       */
      do
      {
         /* check if the profile load conflicts with the demand of the start time variable */
         if( SCIPprofileGetLoad(profile, pos) + demand > capacity )
            peak = pos;

         pos++;
      }
      while( pos < ntimepoints && SCIPprofileGetTime(profile, pos) < ect );

      /* if we found no peak that means current the job could be scheduled at its earliest start time without
       * conflicting to the core resource profile
       */
      if( peak == -1 )
         break;

      /* the peak position gives us a time point where the start time variable is in conflict with the resource
       * profile. That means we have to move it to the next time point in the resource profile but at most to the
       * earliest completion time (the remaining move will done in the next loop)
       */
      newlb = SCIPprofileGetTime(profile, peak+1);
      newlb = MIN(newlb, ect);

      /* if the earliest start time is greater than the lst we detected an infeasibilty */
      if( newlb > lst )
      {
         SCIPdebugMessage("variable <%s>: cannot be scheduled\n", SCIPvarGetName(var));

         /* use conflict analysis to analysis the core insertion which was infeasible */
         SCIP_CALL( analyseInfeasibelCoreInsertion(scip, nvars, vars, durations, demands, capacity,
               var, duration, demand, newlb-1, initialized, explanation) );

         *infeasible = TRUE;

         if( explanation != NULL )
            explanation[pos] = TRUE;

         break;
      }

      /* construct the inference information which we are using with the conflict analysis to resolve that particular
       * bound change
       */
      inferinfo = getInferInfo(PROPRULE_1_CORETIMES, pos, 0);

      /* perform the bound lower bound change */
      SCIP_CALL( SCIPinferVarLbCons(scip, var, (SCIP_Real)newlb, cons, inferInfoToInt(inferinfo), TRUE, infeasible, &tightened) );
      assert(tightened);
      assert(!(*infeasible));

      SCIPdebugMessage("variable <%s> new lower bound <%d> -> <%d>\n", SCIPvarGetName(var), est, newlb);
      (*nchgbds)++;

      /* adjust the earliest start time */
      est = newlb;

      /* adjust the search position for the resource profile for the next step */
      if( est == SCIPprofileGetTime(profile, peak+1) )
         pos = peak + 1;
      else
         pos = peak;
   }
   while( est < lst );

   return SCIP_OKAY;
}

/** We are using the core resource profile which contains all core except the one of the start time variable which we
 *  want to propagate, to decrease the latest start time. This we are doing in steps of length at most the duration of
 *  the job. The reason for that is, that this makes it later easier to resolve this propagation during the conflict
 *  analysis
 */
static
SCIP_RETCODE coretimesUpdateUb(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which is propagated */
   SCIP_PROFILE*         profile,            /**< resource profile */
   int                   pos,                /**< position of the variable to propagate */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            explanation,        /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   SCIP_Bool*            infeasible          /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_VAR* var;
   int ntimepoints;
   int duration;
   int demand;
   int newub;
   int peak;
   int est;
   int lst;
   int lct;

   var = vars[pos];
   assert(var != NULL);

   duration = durations[pos];
   assert(duration > 0);

   demand = demands[pos];
   assert(demand > 0);

   est = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
   lst = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

   /* in case the start time variable is fixed do nothing */
   if( est == lst )
      return SCIP_OKAY;

   ntimepoints = SCIPprofileGetNTimepoints(profile);


   lct = lst + duration;

   /* first we find left position of latest completion time minus 1 (upper bound + duration) in resource profile; That
    * is the last time point where the job would run if schedule it at its latest start time (upper bound). This
    * position gives us the load which we have at the latest completion time minus one
    */
   (void) SCIPprofileFindLeft(profile, lct - 1, &pos);

   SCIPdebugMessage("propagate upper bound (pos %d)\n", pos);
   SCIPdebug( SCIPprofilePrint(profile, SCIPgetMessagehdlr(scip), NULL) );

   if( pos == ntimepoints-1 && SCIPprofileGetTime(profile, pos) == lst )
      return SCIP_OKAY;

   /* we now trying to move the latest start time in steps of at most "duration" length */
   do
   {
      INFERINFO inferinfo;
      SCIP_Bool tightened;

      peak = -1;

#ifndef NDEBUG
      {
         /* in debug mode we check that we adjust the search position correctly */
         int tmppos;

         (void)SCIPprofileFindLeft(profile, lct - 1, &tmppos);
         assert(pos == tmppos);
      }
#endif

      /* we search for a peak within the core profile which conflicts with the demand of the start time variable; we
       * want a peak which is closest to the latest start time
       */
      do
      {
         if( SCIPprofileGetLoad(profile, pos) + demand > capacity )
            peak = pos;

         pos--;
      }
      while( pos >= 0 && SCIPprofileGetTime(profile, pos+1) > lst);

      /* if we found no peak that means current the job could be scheduled at its latest strar time without
       * conflicting to the core resource profile
       */
      if( peak == -1 )
         break;

      /* the peak position gives us a time point where the start time variable is in conflict with the resource
       * profile. That means the job has be done until that point. Hence that gives us the latest completion
       * time. Note that that we want to move the bound by at most the duration length (the remaining move we are
       * doing in the next loop)
       */
      newub = SCIPprofileGetTime(profile, peak);
      newub = MAX(newub, lst) - duration;
      assert(newub >= est);

      /* construct the inference information which we are using with the conflict analysis to resolve that particular
       * bound change
       */
      inferinfo = getInferInfo(PROPRULE_1_CORETIMES, pos, 0);

      /* perform the bound upper bound change */
      SCIP_CALL( SCIPinferVarUbCons(scip, var, (SCIP_Real)newub, cons, inferInfoToInt(inferinfo), TRUE, infeasible, &tightened) );
      assert(tightened);
      assert(!(*infeasible));

      SCIPdebugMessage("variable <%s>: new upper bound <%d> -> <%d>\n", SCIPvarGetName(var), lst, newub);
      (*nchgbds)++;

      /* adjust the latest start and completion time */
      lst = newub;
      lct = lst + duration;

      /* adjust the search position for the resource profile for the next step */
      if( SCIPprofileGetTime(profile, peak) == lct )
         pos = peak - 1;
      else
         pos = peak;
   }
   while( est < lst );

   return SCIP_OKAY;
}

/** a cumulative condition is not satisfied if its capacity is exceeded at a time where jobs cannot be shifted (core)
 *  anymore we build up a cumulative profile of all cores of jobs and try to improve bounds of all jobs
 */
static
SCIP_RETCODE propagateCoretimes(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_CONS*            cons,               /**< constraint which is propagated (needed to SCIPinferVar**Cons()) */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            explanation,        /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_PROFILE* profile;
   SCIP_Bool* cores;
   SCIP_Bool* fixeds;
   SCIP_Bool* ignorejobs;
   int* starts;
   int* ends;

   SCIP_Bool infeasible;

   int ncores;
   int j;

   assert(scip != NULL);
   assert(nvars > 0);
   assert(cons != NULL);
   assert(cutoff !=  NULL);
   assert(*cutoff ==  FALSE);
   assert(*initialized ==  FALSE);

   SCIPdebugMessage("check/propagate cores of cumulative condition of constraint <%s>\n", SCIPconsGetName(cons));

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &cores, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixeds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ignorejobs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &starts, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ends, nvars) );

   infeasible = FALSE;
   ncores = 0;

   /* create an empty resource profile for profiling the cores of the jobs */
   SCIP_CALL( SCIPprofileCreate(&profile, capacity) );

   /* insert all cores */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      int duration;
      int demand;
      int est;
      int lst;

      var = vars[j];
      assert(var != NULL);
      assert(SCIPisFeasIntegral(scip, SCIPvarGetLbLocal(var)));
      assert(SCIPisFeasIntegral(scip, SCIPvarGetUbLocal(var)));

      duration = durations[j];
      assert(duration > 0);

      demand =  demands[j];
      assert(demand > 0);

      /* collect earliest and latest start time */
      est = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      lst = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

      /* check if the job runs completly outside of the effective horizon [hmin, hmax); if so skip it */
      if( lst + duration <= hmin || est >= hmax )
      {
         ignorejobs[j] = TRUE;
         continue;
      }

      ignorejobs[j] = FALSE;

      /* compute core intervall */
      starts[j] = lst;
      ends[j] = est + duration;

      /* check if a core exists */
      if( starts[j] < ends[j] )
      {
         int pos;

         SCIPdebugMessage("variable <%s>[%d,%d] (duration %d, demand %d): add core [%d,%d)\n",
            SCIPvarGetName(var), est, lst, duration, demand, starts[j], ends[j]);

         /* insert the core into core resource profile */
         SCIP_CALL( SCIPprofileInsertCore(profile, starts[j], ends[j], demand, &pos, &infeasible) );

         /* in case the insertion of the core leads to an infeasibility; start the conflict analysis */
         if( infeasible )
         {
            assert(starts[j] <= SCIPprofileGetTime(profile, pos));
            assert(ends[j] > SCIPprofileGetTime(profile, pos));

            /* use conflict analysis to analysis the core insertion which was infeasible */
            SCIP_CALL( analyseInfeasibelCoreInsertion(scip, nvars, vars, durations, demands, capacity,
                  var, duration, demand, SCIPprofileGetTime(profile, pos), initialized, explanation) );

            *cutoff = TRUE;

            if( explanation != NULL )
               explanation[pos] = TRUE;

            break;
         }

         /* remenber that the job has a core */
         cores[j] = TRUE;
         ncores++;

         /* check if the start time variables is already fixed; in that case we can ignore the job */
         if( est == lst )
            ignorejobs[j] = TRUE;
      }
      else
         cores[j] = FALSE;
   }

   if( !(*cutoff) && ncores > 0 )
   {
      /* start checking each job whether the bounds can be improved */
      for( j = 0; j < nvars; ++j )
      {
         SCIP_VAR* var;
         int demand;
         int duration;
         int est;
         int lst;

         /* check if the job can be ignored */
         if( ignorejobs[j] )
            continue;

         var = vars[j];
         assert(var != NULL);

         duration = durations[j];
         assert(duration > 0);

         demand = demands[j];
         assert(demand > 0);

         /* if the job has a core, remove it first */
         if( cores[j] )
         {
            SCIPdebugMessage("variable <%s>[%g,%g] (duration %d, demand %d): remove core [%d,%d)\n",
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), duration, demand, starts[j], ends[j]);

            SCIP_CALL( SCIPprofileDeleteCore(profile, starts[j], ends[j], demand) );
         }

         /* first try to updates the earliest start time */
         SCIP_CALL( coretimesUpdateLb(scip, nvars, vars, durations, demands, capacity, cons,
               profile, j, nchgbds, initialized, explanation, cutoff) );

         if( *cutoff )
            break;

         SCIP_CALL( coretimesUpdateUb(scip, nvars, vars, durations, demands, capacity, cons,
               profile, j, nchgbds, explanation, cutoff) );

         if( *cutoff )
            break;

         /* collect earliest and latest start time */
         est = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
         lst = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

         /* compute core intervall */
         starts[j] = lst;
         ends[j] = est + duration;

         /* after updating the bound we might have a new core */
         if( starts[j] < ends[j] )
         {
            int pos;

            SCIPdebugMessage("variable <%s>[%d,%d] (duration %d, demand %d): remove core [%d,%d)\n",
               SCIPvarGetName(var), est, lst, duration, demand, starts[j], ends[j]);

            SCIP_CALL( SCIPprofileInsertCore(profile, starts[j], ends[j], demand, &pos, &infeasible) );

            if( infeasible )
            {
               /* use conflict analysis to analysis the core insertion which was infeasible */
               SCIP_CALL( analyseInfeasibelCoreInsertion(scip, nvars, vars, durations, demands, capacity,
                     var, duration, demand, SCIPprofileGetTime(profile, pos), initialized, explanation) );

               *cutoff = TRUE;

               if( explanation != NULL )
                  explanation[pos] = TRUE;

               break;
            }

            cores[j] = TRUE;
         }
      }
   }

   /* free resource profile */
   SCIPprofileFree(&profile);

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &ends);
   SCIPfreeBufferArray(scip, &starts);
   SCIPfreeBufferArray(scip, &ignorejobs);
   SCIPfreeBufferArray(scip, &fixeds);
   SCIPfreeBufferArray(scip, &cores);

   return SCIP_OKAY;
}


struct SCIP_Envelop
{
   SCIP_VAR*             var;                /**< start time variable of the job */
   SCIP_Real             key;                /**< key which is to insert the corresponding search node */
   int                   envelop;            /**< envelop of ?????????????? */
   int                   energy;             /**< energy of the job */
};
typedef struct SCIP_Envelop SCIP_ENVELOP;

/** update node data structure strating form the given node along the path to the root node */
static
SCIP_RETCODE updateEnvelop(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BSTNODE*         node,               /**< search node which inserted */
   SCIP_ENVELOP**        nodedatas,          /**< array with search node data */
   int*                  nnodedatas          /**< pointer to store the number of search node data structure */
   )
{
   SCIP_BSTNODE* left;
   SCIP_BSTNODE* right;
   SCIP_ENVELOP* nodedata;
   SCIP_ENVELOP* leftdata;
   SCIP_ENVELOP* rightdata;

   node = SCIPbstnodeGetParent(node);

   while( node != NULL )
   {
      /* get node data */
      nodedata = (SCIP_ENVELOP*)SCIPbstnodeGetData(node);

      /* if the node data is NULL, we have an internal node which was created due to inserting a new search node */
      if( nodedata == NULL )
      {
         SCIP_CALL( SCIPallocBuffer(scip, &nodedata) );
         nodedata->var = NULL;
         nodedata->key = SCIP_INVALID;
         nodedata->envelop = 0;
         nodedata->energy = 0;

         /* attach the node data to the search node */
         SCIPbstnodeSetData(node, (void*)nodedata);

         /* store node data to be able to delete them latter */
         nodedatas[*nnodedatas] = nodedata;
         (*nnodedatas)++;
      }
      assert(nodedata != NULL);

      /* collect node data from left node */
      left = SCIPbstnodeGetLeftchild(node);
      assert(left != NULL);
      leftdata = (SCIP_ENVELOP*)SCIPbstnodeGetData(left);
      assert(leftdata != NULL);

      /* collect node data from right node */
      right = SCIPbstnodeGetRightchild(node);
      assert(right != NULL);
      rightdata = (SCIP_ENVELOP*)SCIPbstnodeGetData(right);
      assert(rightdata != NULL);

      /* update envelop and energy */
      nodedata->envelop = MAX(leftdata->envelop + rightdata->energy, rightdata->envelop);
      nodedata->energy = leftdata->energy + rightdata->energy;

      /* go to parent */
      node = SCIPbstnodeGetParent(node);
   }

   return SCIP_OKAY;
}

/** insert mmethod for theta trees */
static
SCIP_DECL_BSTINSERT(thetatreeInsert)
{
   (*inserted) = TRUE;

   /* if the tree is empty the node will be the root node */
   if( SCIPbstIsEmpty(tree) )
   {
      SCIPbstSetRoot(tree, node);
   }
   else
   {
      SCIP_BSTNODE* leaf;
      SCIP_BSTNODE* newnode;
      SCIP_BSTNODE* parent;

      leaf = SCIPbstGetRoot(tree);
      assert(leaf != NULL);

      /* find the position to insert the node */
      while( !SCIPbstnodeIsLeaf(leaf) )
      {
         if( SCIPbstComp(tree, node, leaf) <= 0 )
            leaf = SCIPbstnodeGetLeftchild(leaf);
         else
            leaf = SCIPbstnodeGetRightchild(leaf);
      }
      assert(leaf != NULL);
      assert(leaf != node);

      /* create a new node */
      SCIP_CALL( SCIPbstnodeCreate(tree, &newnode, NULL, NULL) );
      assert(newnode != NULL);

      parent = SCIPbstnodeGetParent(leaf);

      if( parent != NULL )
      {
         SCIPbstnodeSetParent(newnode, parent);

         /* check if the node is the left child */
         if( SCIPbstnodeGetLeftchild(parent) == leaf )
         {
            SCIPbstnodeSetLeftchild(parent, newnode);
         }
         else
         {
            SCIPbstnodeSetRightchild(parent, newnode);
         }
      }
      else
         SCIPbstSetRoot(tree, newnode);

      if( SCIPbstComp(tree, node, leaf) <= 0 )
      {
         /* node is on the left */
         SCIPbstnodeSetLeftchild(newnode, node);
         SCIPbstnodeSetRightchild(newnode, leaf);
         SCIPbstnodeSetKey(newnode, SCIPbstnodeGetKey(node));
      }
      else
      {
         /* leaf is on the left */
         SCIPbstnodeSetLeftchild(newnode, leaf);
         SCIPbstnodeSetRightchild(newnode, node);
         SCIPbstnodeSetKey(newnode, SCIPbstnodeGetKey(leaf));
      }

      SCIPbstnodeSetParent(leaf, newnode);
      SCIPbstnodeSetParent(node, newnode);
   }

   return SCIP_OKAY;
}

/** compare the keys of the search nodes belonging to a theta tree */
static
SCIP_DECL_SORTPTRCOMP(thetatreeComp)
{
   SCIP_Real key1;
   SCIP_Real key2;

   key1 = ((SCIP_ENVELOP*)elem1)->key;
   key2 = ((SCIP_ENVELOP*)elem2)->key;

   if( key1 < key2 )
      return -1;
   else if( key1 > key2 )
      return 1;

   return 0;
}

/** checks whether the instance is infeasible due to a overload within a certain time frame
 *
 *  @note The algorithm is based on the paper: Petr Vilim, "Max Energy Filtering Algorithm for Discrete Cumulative
 *        Resources". In: Willem Jan van Hoeve and John N. Hooker (Eds.), Integration of AI and OR Techniques in
 *        Constraint Programming for Combinatorial Optimization Problems (CPAIOR 2009), LNCS 5547, pp 294--308
 */
static
SCIP_RETCODE checkOverload(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_CONS*            cons,               /**< constraint which is propagated */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            explanation,        /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_ENVELOP** nodedatas;
   SCIP_BST* bsttree;

   int* ests;
   int* lcts;
   int* perm;

   int nnodedatas;
   int ncands;

   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(initialized != NULL);
   assert(cutoff != NULL);

   SCIPdebugMessage("check overload of cumulative condition of constraint <%s> (capacity %d)\n", SCIPconsGetName(cons), capacity);

   SCIP_CALL( SCIPallocBufferArray(scip, &ests, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lcts, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodedatas, 2*nvars) );

   ncands = 0;

   SCIP_CALL( SCIPbstCreate(&bsttree, SCIPblkmem(scip), thetatreeInsert, NULL, thetatreeComp) );

   /* collect earliest and latest completion times and ignore jobs which do not run completion within the effective
    * horizon
    */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      int est;
      int lct;
      int energy;

      var = vars[j];
      assert(var != NULL);

      est = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      lct = convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + durations[j];

      /* ignore jobs which do not lie completely in the effective horizon [hmin,hmax) */
      if( est < hmin || lct > hmax )
      {
         nodedatas[j] = NULL;
         continue;
      }

      energy = demands[j] * durations[j];
      assert(energy > 0);

      /* create search node data */
      SCIP_CALL( SCIPallocBuffer(scip, &nodedatas[j]) );

      /* initialize search node data */
      nodedatas[j]->var = var;

      /* adjust earliest start time to make it unique in case several jobs have the same earliest start time */
      nodedatas[j]->key = est + j / (2.0 * nvars);

      /* the envelop is the energy of the job plus the total amount of energy which is available in the time period
       * before that job can start, that is [0,est). The envelop is later used to compare the energy consumption of a
       * particular time interval [a,b] against the time interval [0,b].
       */
      nodedatas[j]->envelop = capacity * est + energy;
      nodedatas[j]->energy = energy;

      ests[ncands] = est;
      lcts[ncands] = lct;
      perm[ncands] = j;
      ncands++;
   }

   nnodedatas = nvars;

   /* sort the latest completion times */
   SCIPsortIntIntInt(lcts, ests, perm, ncands);

   /* iterate over all jobs in non-decreasing order of their latest completion times */
   for( j = 0; j < ncands && !(*cutoff); ++j )
   {
      SCIP_VAR* var;
      SCIP_BSTNODE* root;
      SCIP_BSTNODE* node;
      SCIP_ENVELOP* data;
      SCIP_Bool inserted;
      int idx;

      idx = perm[j];

      /* create search node */
      SCIP_CALL( SCIPbstnodeCreate(bsttree, &node, (void*)nodedatas[idx], (void*)nodedatas[idx]) );

      /* insert new search node */
      SCIP_CALL( SCIPbstInsert(bsttree, node, &inserted) );
      assert(inserted);

      /* update envelop */
      SCIP_CALL( updateEnvelop(scip, node, nodedatas, &nnodedatas) );
      assert(nnodedatas <= 2*nvars);

      root = SCIPbstGetRoot(bsttree);
      data = SCIPbstnodeGetData(root);

      /* check for overload */
      if( data->envelop > capacity * lcts[j] )
      {
         (*cutoff) = TRUE;

         /* check if the conflict analysis is applicable */
         if( SCIPisConflictAnalysisApplicable(scip) )
         {
            int reportedenergy;
            int energy;
            int est;
            int lct;
            int c;

            /* collect the earliest and latest completion time of the last job which lead to the overload detection */
            est = ests[j];
            lct = lcts[j];

            /* adjust number of candidates to the once we added to the search tree */
            ncands = j+1;

            /* sort the start time variables which were added to search tree w.r.t. earliest start time */
            SCIPsortDownIntIntInt(ests, lcts, perm, ncands);

            energy = (lct - est) * capacity;
            reportedenergy = 0;
            c = 0;

            /* collect the energy of those jobs which run within the time frame [est,lct) of the job which led to
             * infeasibility
             */
            for( c = 0; c < ncands && ests[c] >= est && reportedenergy <= energy; ++c )
            {
               assert(lcts[c] <= lct);

               idx = perm[c];
               var = vars[idx];

               /* collect energy which is contributed by the start time variable */
               reportedenergy += demands[idx] * durations[idx];

               SCIPdebugMessage("variable <%s>: loc=[%g,%g] glb=[%g,%g] (duration %d, demand %d)\n",
                  SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
                  SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), durations[idx], demands[idx]);
            }

            /* continue with remaining candidates and adjust the earliest start time */
            for( ; c < ncands && reportedenergy <= energy; ++c )
            {
               assert(lcts[c] <= lct);

               idx = perm[c];
               var = vars[idx];

               reportedenergy += demands[idx] * durations[idx];

               SCIPdebugMessage("variable <%s>: loc=[%g,%g] glb=[%g,%g] (duration %d, demand %d)\n",
                  SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
                  SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), durations[idx], demands[idx]);

               /* adjust energy */
               est = ests[c];
               assert(energy <= (lct - est) * capacity);
               energy = (lct - est) * capacity;
            }
            assert(reportedenergy > energy);

            /* shrink number of candidates to once we need to report */
            ncands = c;

            /* initialize conflict analysis */
            SCIP_CALL( SCIPinitConflictAnalysis(scip) );

            /* report the variables and relax their bounds to final time interval [est,lct) which was been detected to
             * be overloaded
             */
            for( c = 0; c < ncands; ++c )
            {
               idx = perm[c];
               var = vars[idx];
               assert(var != NULL);

               SCIP_CALL( relaxVarLb(scip, var, est, NULL) );
               SCIP_CALL( relaxVarUb(scip, var, lct - durations[idx], NULL) );

               if( explanation != NULL )
                  explanation[idx] = TRUE;
            }

            (*initialized) = TRUE;
         }
         else
            assert((SCIPgetDepth(scip) == 0 || SCIPgetStage(scip) != SCIP_STAGE_SOLVING) && !SCIPinProbing(scip));
      }
   }

   /* free the search nodes data */
   for( j = nnodedatas - 1; j >= 0; --j )
   {
      SCIPfreeBufferNull(scip, &nodedatas[j]);
   }

   /* free theta tree */
   SCIPbstFree(&bsttree);

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &nodedatas);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &lcts);
   SCIPfreeBufferArray(scip, &ests);

   return SCIP_OKAY;
}
/** checks if the constraint is redundant; that is if its capacity can never be exceeded; therefore we check with
 *  respect to the lower and upper bounds of the integer variables the maximum capacity usage for all event points
 */
static
SCIP_RETCODE consCheckRedundancy(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_Bool*            redundant           /**< pointer to store whether this constraint is redundant */
   )
{

   SCIP_VAR* var;
   int* starttimes;              /* stores when each job is starting */
   int* endtimes;                /* stores when each job ends */
   int* startindices;            /* we will sort the startsolvalues, thus we need to know wich index of a job it corresponds to */
   int* endindices;              /* we will sort the endsolvalues, thus we need to know wich index of a job it corresponds to */

   int lb;
   int ub;
   int freecapacity;             /* remaining capacity */
   int curtime;                  /* point in time which we are just checking */
   int endindex;                 /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int j;

   assert(scip != NULL);
   assert(redundant != NULL);

   (*redundant) = TRUE;

   /* if no activities are associated with this cumulative then this constraint is redundant */
   if( nvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   /* assign variables, start and endpoints to arrays */
   for( j = 0; j < nvars; ++j )
   {
      var = vars[j];
      assert(var != NULL);

      lb = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      ub = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

      starttimes[j] = MAX(lb, hmin);
      startindices[j] = j;

      endtimes[j] =  ub + durations[j];
      endindices[j] = j;
   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues (and sort the indices in the same way) */
   SCIPsortIntInt(starttimes, startindices, nvars);
   SCIPsortIntInt(endtimes, endindices, nvars);

   endindex = 0;
   freecapacity = capacity;

   /* check each start point of a job whether the capacity is violated or not */
   for( j = 0; j < nvars; ++j )
   {
      curtime = starttimes[j];

      /* stop checking, if time point is above hmax */
      if( curtime >= hmax )
         break;

      /* subtract all capacity needed up to this point */
      freecapacity -= demands[startindices[j]];
      while( j+1 < nvars && starttimes[j+1] == curtime )
      {
         ++j;
         freecapacity -= demands[startindices[j]];
      }

      /* free all capacity usages of jobs the are no longer running */
      while( endtimes[endindex] <= curtime )
      {
         freecapacity += demands[endindices[endindex]];
         ++endindex;
      }
      assert(freecapacity <= capacity);

      /* check freecapacity to be smaller than zero */
      if( freecapacity < 0 && curtime >= hmin )
      {
         (*redundant) = FALSE;
         break;
      }
   } /*lint --e{850}*/

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endtimes);
   SCIPfreeBufferArray(scip, &starttimes);

   return SCIP_OKAY;
}

/** propagate the cumulative condition */
static
SCIP_RETCODE propagateCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_CONS*            cons,               /**< constraint which is propagated (needed to SCIPinferVar**Cons()) */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            redundant,          /**< pointer to store if the constraint is redundant */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            explanation,        /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   assert(nchgbds != NULL);
   assert(initialized != NULL);
   assert(cutoff != NULL);
   assert((*cutoff) == FALSE);

   /**@todo avoid always sorting the variable array */

   /* check if the constraint is redundant */
   SCIP_CALL( consCheckRedundancy(scip, nvars, vars, durations, demands, capacity, hmin, hmax, redundant) );

   if( *redundant )
      return SCIP_OKAY;

   /* propagate the job cores until nothing else can be detected */
   if( conshdlrdata->coretimes )
   {
      SCIP_CALL( propagateCoretimes(scip, nvars, vars, durations, demands, capacity, hmin, hmax, cons,
            nchgbds, initialized, explanation, cutoff) );

      if( *cutoff )
         return SCIP_OKAY;
   }

   if( conshdlrdata->overload )
   {
      /* check for overload, which may result in a cutoff */
      SCIP_CALL( checkOverload(scip, nvars, vars, durations, demands, capacity, hmin, hmax,
            cons, initialized, explanation, cutoff) );

      if( *cutoff )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** propagate the cumulative constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to propagate */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   int*                  ndelconss,          /**< pointer to store the number of deleted constraints */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool initialized;
   SCIP_Bool redundant;
   int oldnchgbds;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   oldnchgbds = *nchgbds;
   initialized = FALSE;
   redundant = FALSE;

   /* if the constraint marked to be propagated, do nothing */
   if( consdata->propagated )
      return SCIP_OKAY;

   SCIP_CALL( propagateCumulativeCondition(scip, conshdlrdata,
         consdata->nvars, consdata->vars, consdata->durations, consdata->demands, consdata->capacity,
         consdata->hmin, consdata->hmax, cons,
         nchgbds, &redundant, &initialized, NULL, cutoff) );

   if( redundant )
   {
      SCIPdebugMessage("%s deletes cumulative constraint <%s> since it is redundant\n",
         SCIPgetDepth(scip) == 0 ? "globally" : "locally", SCIPconsGetName(cons));

      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      (*ndelconss)++;
   }
   else
   {
      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && initialized && !SCIPisStopped(scip) )
      {
         /* run conflict analysis since it was initialized */
         assert(*cutoff == TRUE);
         SCIPdebugMessage("start conflict analysis\n");
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }

      /* if successful, reset age of constraint */
      if( *cutoff || *nchgbds > oldnchgbds )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   /* mark the constraint to be propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** computes w.r.t. the given worst case resource profile the first time point where the given capacity can be violated */
static
int computeHmin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROFILE*         profile,            /**< worst case resource profile */
   int                   capacity            /**< capacity to check */
   )
{
   int* timepoints;
   int* loads;
   int ntimepoints;
   int t;

   ntimepoints = SCIPprofileGetNTimepoints(profile);
   timepoints = SCIPprofileGetTimepoints(profile);
   loads = SCIPprofileGetLoads(profile);

   /* find first time point which potentially violates the capacity restriction */
   for( t = 0; t < ntimepoints - 1; ++t )
   {
      /* check if the time point exceed w.r.t. worst case profile the capacity */
      if( loads[t] > capacity )
      {
         assert(t == 0 || loads[t-1] <= capacity);
         return timepoints[t];
      }
   }

   return INT_MAX;
}

/** computes w.r.t. the given worst case resource profile the first time point where the given capacity is satisfied for sure */
static
int computeHmax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROFILE*         profile,            /**< worst case profile */
   int                   capacity            /**< capacity to check */
   )
{
   int* timepoints;
   int* loads;
   int ntimepoints;
   int t;

   ntimepoints = SCIPprofileGetNTimepoints(profile);
   timepoints = SCIPprofileGetTimepoints(profile);
   loads = SCIPprofileGetLoads(profile);

   /* find last time point which potentially violates the capacity restriction */
   for( t = ntimepoints - 1; t >= 0; --t )
   {
      /* check if at time point t the worst case resource profile  exceeds the capacity */
      if( loads[t] > capacity )
      {
         assert(t == ntimepoints-1 || loads[t+1] <= capacity);
         return timepoints[t+1];
      }
   }

   return INT_MIN;
}
/** For each variable we compute an alternative lower and upper bounds. That is, if the variable is not fixed to its
 *  lower or upper bound the next reasonable lower or upper bound would be this alternative bound (implying that certain
 *  values are not of interest). An alternative bound for a particular is only valied if the cumulative constarints are
 *  the only one locking this variable in the corresponding direction.
 */
static
SCIP_RETCODE computeAlternativeBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of cumulative constraint constraints */
   int                   nconss,             /**< number of cumulative constraints */
   SCIP_Bool             local,              /**< use local bounds effective horizon? */
   int*                  alternativelbs,     /**< alternative lower bounds */
   int*                  alternativeubs,     /**< alternative lower bounds */
   int*                  downlocks,          /**< number of constraints with down lock participating by the computation */
   int*                  uplocks             /**< number of constraints with up lock participating by the computation */
   )
{
   int c;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      SCIP_VAR* var;
      int nvars;
      int hmin;
      int hmax;
      int v;

      cons = conss[c];
      assert(cons != NULL);

      /* ignore constraints which are already deletet and those which are not check constraints */
      if( SCIPconsIsDeleted(cons) || !SCIPconsIsChecked(cons) )
         continue;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      assert(consdata->nvars > 1);

      /* compute the hmin and hmax */
      if( local )
      {
         SCIP_PROFILE* profile;

         /* create empty resource profile with infinity resource capacity */
         SCIP_CALL( SCIPprofileCreate(&profile, INT_MAX) );

         SCIP_CALL( createWorstCaseProfile(scip, profile, consdata) );

         hmin = computeHmin(scip, profile, consdata->capacity);
         hmax = computeHmax(scip, profile, consdata->capacity);

         /* free worst case profile */
         SCIPprofileFree(&profile);
      }
      else
      {
         hmin = consdata->hmin;
         hmax = consdata->hmax;
      }

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      nvars = consdata->nvars;

      for( v = 0; v < nvars; ++v )
      {
         int constant;
         int idx;

         var = consdata->vars[v];
         assert(var != NULL);

         /* ignore variable locally fixed variables */
         if( SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var) < 0.5 )
            continue;

         /* ignore multi-aggregated variables */
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
            continue;

         if( !SCIPvarIsActive(var) )
         {
            SCIP_VAR* aggrvar;

            assert(SCIPvarGetUbLocal(var)  - SCIPvarGetLbLocal(var) < 1.5);
            assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED);
            assert(SCIPisEQ(scip, SCIPvarGetAggrScalar(var), 1.0));

            aggrvar = SCIPvarGetAggrVar(var);
            constant = convertBoundToInt(scip, SCIPvarGetAggrConstant(var));
            assert(SCIPvarIsActive(aggrvar));
            assert(SCIPvarGetType(aggrvar) == SCIP_VARTYPE_BINARY);

            idx = SCIPvarGetProbindex(aggrvar);
         }
         else
         {
            idx = SCIPvarGetProbindex(var);
            constant = 0;
         }
         assert(idx >= 0);

         /* first check lower bound fixing */
         if( consdata->downlocks[v] )
         {
            int ect;
            int est;

            /* the variable has a down locked */
            est = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
            ect = est + consdata->durations[v];

            if( ect <= hmin || hmin >= hmax )
               downlocks[idx]++;
            else if( est < hmin && alternativelbs[idx] >= hmin + 1 - constant )
            {
               alternativelbs[idx] = hmin + 1 - constant;
               downlocks[idx]++;
            }
         }

         /* second check upper bound fixing */
         if( consdata->uplocks[v] )
         {
            int lct;
            int lst;

            /* the variable has a up lock locked */
            lst = convertBoundToInt(scip, SCIPvarGetUbLocal(var));
            lct = lst + consdata->durations[v];

            if( lst >= hmax || hmin >= hmax  )
               uplocks[idx]++;
            else if( lct > hmax && alternativeubs[idx] <= hmax - 1 - constant )
            {
               alternativeubs[idx] = hmax - 1 - constant;
               uplocks[idx]++;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** apply all fixings which are given by the alternative bounds */
static
SCIP_RETCODE applyAlternativeBoundsFixing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of active variables */
   int                   nvars,              /**< number of active variables */
   int*                  alternativelbs,     /**< alternative lower bounds */
   int*                  alternativeubs,     /**< alternative lower bounds */
   int*                  downlocks,          /**< number of constraints with down lock participating by the computation */
   int*                  uplocks,            /**< number of constraints with up lock participating by the computation */
   int*                  nfixedvars          /**< pointer to store the number of fixed variables */
   )
{
   int v;

   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;
      int ub;
      int lb;

      var = vars[v];
      assert(var != NULL);

      lb = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      ub = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

      /* ignore fixed variables */
      if( ub - lb < 0.5 )
         continue;

      if( SCIPvarGetNLocksDown(var) == downlocks[v] && SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_LOWER )
      {

         if( alternativelbs[v] > ub )
         {
            SCIP_CALL( SCIPfixVar(scip, var, SCIPvarGetLbLocal(var), &infeasible, &fixed) );
            assert(!infeasible);
            assert(fixed);

            (*nfixedvars)++;
         }
      }

      if( SCIPvarGetNLocksUp(var) == uplocks[v] && SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_UPPER )
      {
         if( alternativeubs[v] < lb )
         {
            SCIP_CALL( SCIPfixVar(scip, var, SCIPvarGetUbLocal(var), &infeasible, &fixed) );
            assert(!infeasible);
            assert(fixed);

            (*nfixedvars)++;
         }
      }
   }

   return SCIP_OKAY;
}

/** propagate all constraints together */
static
SCIP_RETCODE propagateAllConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS**           conss,              /**< all cumulative constraint */
   int                   nconss,             /**< number of cumulative constraints */
   SCIP_Bool             local,              /**< use local bounds effective horizon? */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   SCIP_Bool*            branched            /**< pointer to store if a branching was applied, or NULL to avoid branching */
   )
{
   SCIP_VAR** vars;
   int* downlocks;
   int* uplocks;
   int* alternativelbs;
   int* alternativeubs;
   int oldnfixedvars;
   int nvars;
   int v;

   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   nvars = SCIPgetNVars(scip);
   oldnfixedvars = *nfixedvars;

   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, SCIPgetVars(scip), nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &downlocks, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &uplocks, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &alternativelbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &alternativeubs, nvars) );

   /* initialize arrays */
   for( v = 0; v < nvars; ++v )
   {
      downlocks[v] = 0;
      uplocks[v] = 0;
      alternativelbs[v] = INT_MAX;
      alternativeubs[v] = INT_MIN;
   }

   SCIPstatistic( conshdlrdata->nallconsdualfixs -= *nfixedvars );

   /* compute alternative bounds */
   SCIP_CALL( computeAlternativeBounds(scip, conss, nconss, local, alternativelbs, alternativeubs, downlocks, uplocks) );

   /* apply fixing which result of the alternative bounds directly */
   SCIP_CALL( applyAlternativeBoundsFixing(scip, vars, nvars, alternativelbs, alternativeubs, downlocks, uplocks, nfixedvars) );

   SCIPstatistic( conshdlrdata->nallconsdualfixs += *nfixedvars );

   if( oldnfixedvars == *nfixedvars && branched != NULL )
   {
      SCIP_CALL( applyAlternativeBoundsBranching(scip, vars, nvars, alternativelbs, alternativeubs, downlocks, uplocks, branched) );
   }

   /* free all buffers */
   SCIPfreeBufferArray(scip, &alternativeubs);
   SCIPfreeBufferArray(scip, &alternativelbs);
   SCIPfreeBufferArray(scip, &uplocks);
   SCIPfreeBufferArray(scip, &downlocks);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/**@} */

/**@name Linear relaxations
 *
 * @{
 */

/** creates covering cuts for jobs violating resource constraints */
static
SCIP_RETCODE createCoverCutsTimepoint(
   SCIP*            scip,                 /**< SCIP data structure */
   SCIP_CONS*       cons,                 /**< constraint to be checked */
   int*             startvalues,          /**< upper bounds on finishing time per job for activities from 0,..., nactivities -1 */
   int              time                  /**< at this point in time covering constraints are valid */
   )
{
   SCIP_VAR** binvars;    /* binary variables of some integer variable */
   SCIP_CONSDATA* consdata;
   SCIP_ROW* row;
   int* flexibleids;
   int* demands;

   char rowname[SCIP_MAXSTRLEN];

   int remainingcap;
   int smallcoversize;    /* size of a small cover */
   int bigcoversize;    /* size of a big cover */
   int nbinvars;
   int offset;
   int nvars;

   int nflexible;
   int D;            /* demand of all jobs up to a certain index */
   int j;
   int i;

   assert(cons != NULL);

   /* get constraint data structure */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL );

   nvars = consdata->nvars;

   /* sort jobs according to demands */
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &flexibleids, nvars) );

   nflexible = 0;
   remainingcap = consdata->capacity;

   /* get all jobs intersecting point 'time' with their bounds */
   for( j = 0; j < nvars; ++j )
   {
      int ub;

      ub = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[j]));

      /* only add jobs to array if they intersect with point 'time' */
      if( startvalues[j] <= time && ub + consdata->durations[j] > time )
      {
         /* if job is fixed, capacity has to be decreased */
         if( startvalues[j] == ub )
         {
            remainingcap -= consdata->demands[j];
         }
         else
         {
            demands[nflexible] = consdata->demands[j];
            flexibleids[nflexible] = j;
            ++nflexible;
         }
      }
   }
   assert(remainingcap >= 0);

   /* sort demands and job ids */
   SCIPsortIntInt(demands, flexibleids, nflexible);

   /*
    * version 1:
    * D_j := sum_i=0,...,j  d_i, finde j maximal, so dass D_j <= remainingcap
    * erzeuge cover constraint
    *
    */

   /* find maximum number of jobs that can run in parallel (-->coversize = j) */
   D = 0;
   j = 0;

   while( j < nflexible && D <= remainingcap )
   {
      D += demands[j];
      ++j;
   }

   /* j jobs form a conflict, set coversize to 'j - 1' */
   bigcoversize = j-1;
   assert(D > remainingcap);
   assert(bigcoversize < nflexible);

   /* - create a row for all jobs and their binary variables.
    * - at most coversize many binary variables of jobs can be set to one
    */

   /* construct row name */
   (void)SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "capacity_coverbig_%d", time);
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), rowname, -SCIPinfinity(scip), (SCIP_Real)bigcoversize,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for( j = 0; j < nflexible; ++j )
   {
      int idx;
      int end;
      int start;
      int lb;
      int ub;

      idx = flexibleids[j];

      /* get and add binvars into var array */
      SCIP_CALL( SCIPgetBinvarsLinking(scip, consdata->linkingconss[idx], &binvars, &nbinvars) );
      assert(nbinvars != 0);
      offset = SCIPgetOffsetLinking(scip, consdata->linkingconss[idx]);

      lb = convertBoundToInt(scip, SCIPvarGetLbLocal(consdata->vars[idx]));
      ub = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[idx]));
      /* compute start and finishing time */
      start = MAX(lb, time + 1 - consdata->durations[idx]) - offset;
      end =  MIN(time, ub) + 1 - offset;

      /* add all neccessary binary variables */
      for( i = start; i < end; ++i )
      {
         assert(i >= 0);
         assert(i < nbinvars);
         assert(binvars[i] != NULL);
         SCIP_CALL( SCIPaddVarToRow(scip, row, binvars[i], 1.0) );
      }
   }

   /* insert and release row */
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   if( consdata->bcoverrowssize == 0 )
   {
      consdata->bcoverrowssize = 10;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->bcoverrows, consdata->bcoverrowssize) );
   }
   if( consdata->nbcoverrows == consdata->bcoverrowssize )
   {
      consdata->bcoverrowssize *= 2;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->bcoverrows, consdata->nbcoverrows, consdata->bcoverrowssize) );
   }

   consdata->bcoverrows[consdata->nbcoverrows] = row;
   consdata->nbcoverrows++;

   /*
    * version 2:
    * D_j := sum_i=j,...,0  d_i, finde j minimal, so dass D_j <= remainingcap
    * erzeuge cover constraint und fuege alle jobs i hinzu, mit d_i = d_largest
    */
   /* find maximum number of jobs that can run in parallel (= coversize -1) */
   D = 0;
   j = nflexible -1;
   while( D <= remainingcap )
   {
      assert(j >= 0);
      D += demands[j];
      --j;
   }

   smallcoversize = nflexible - (j + 1) - 1;
   while( j > 0 && demands[j] == demands[nflexible-1] )
      --j;

   assert(smallcoversize < nflexible);

   if( smallcoversize != 1 || smallcoversize != nflexible - (j + 1) - 1 )
   {
      /* construct row name */
      (void)SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "capacity_coversmall_%d", time);
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), rowname, -SCIPinfinity(scip), (SCIP_Real)smallcoversize,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      /* filter binary variables for each unfixed job */
      for( j = j + 1; j < nflexible; ++j )
      {
         int idx;
         int end;
         int start;
         int lb;
         int ub;

         idx = flexibleids[j];

         /* get and add binvars into var array */
         SCIP_CALL( SCIPgetBinvarsLinking(scip, consdata->linkingconss[idx], &binvars, &nbinvars) );
         assert(nbinvars != 0);
         offset = SCIPgetOffsetLinking(scip, consdata->linkingconss[idx]);

         lb = convertBoundToInt(scip, SCIPvarGetLbLocal(consdata->vars[idx]));
         ub = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[idx]));
         /* compute start and finishing time */
         start = MAX(lb, time + 1 - consdata->durations[idx]) - offset;
         end =  MIN(time, ub) + 1 - offset;

         /* add  all neccessary binary variables */
         for( i = start; i < end; ++i )
         {
            assert(i >= 0);
            assert(i < nbinvars);
            assert(binvars[i] != NULL);
            SCIP_CALL( SCIPaddVarToRow(scip, row, binvars[i], 1.0) );
         }
      }

      /* insert and release row */
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
      if( consdata->scoverrowssize == 0 )
      {
         consdata->scoverrowssize = 10;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->scoverrows, consdata->scoverrowssize) );
      }
      if( consdata->nscoverrows == consdata->scoverrowssize )
      {
         consdata->scoverrowssize *= 2;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->scoverrows, consdata->nscoverrows, consdata->scoverrowssize) );
      }

      consdata->scoverrows[consdata->nscoverrows] = row;
      consdata->nscoverrows++;
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &flexibleids);
   SCIPfreeBufferArray(scip, &demands);

   return SCIP_OKAY;
}

/** method to construct cover cuts for all points in time */
static
SCIP_RETCODE createCoverCuts(
   SCIP*            scip,                      /**< SCIP data structure */
   SCIP_CONS*       cons                       /**< constraint to be separated */
   )
{
   SCIP_CONSDATA* consdata;

   int* startvalues;        /* stores when each job is starting */
   int* endvalues;          /* stores when each job ends */
   int* startvaluessorted;  /* stores when each job is starting */
   int* endvaluessorted;    /* stores when each job ends */
   int* startindices;     /* we sort the startvalues, so we need to know wich index of a job it corresponds to */
   int* endindices;       /* we sort the endvalues, so we need to know wich index of a job it corresponds to */

   int nvars;               /* number of jobs for this constraint */
   int freecapacity;        /* remaining capacity */
   int curtime;             /* point in time which we are just checking */
   int endidx;              /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int hmin;
   int hmax;

   int j;
   int t;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* if no activities are associated with this resource then this constraint is redundant */
   if( consdata->vars == NULL )
      return SCIP_OKAY;

   nvars = consdata->nvars;
   hmin = consdata->hmin;
   hmax = consdata->hmax;

   SCIP_CALL( SCIPallocBufferArray(scip, &startvalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endvalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startvaluessorted, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endvaluessorted, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   /* assign start and endpoints to arrays */
   for ( j = 0; j < nvars; ++j )
   {
      startvalues[j] = convertBoundToInt(scip, SCIPvarGetLbLocal(consdata->vars[j]));
      startvaluessorted[j] = startvalues[j];

      endvalues[j] = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[j])) + consdata->durations[j];
      endvaluessorted[j] = endvalues[j];

      startindices[j] = j;
      endindices[j] = j;
   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues
    * (and sort the indices in the same way) */
   SCIPsortIntInt(startvaluessorted, startindices, nvars);
   SCIPsortIntInt(endvaluessorted, endindices, nvars);

   endidx = 0;
   freecapacity = consdata->capacity;

   /* check each startpoint of a job whether the capacity is kept or not */
   for( j = 0; j < nvars; ++j )
   {
      curtime = startvaluessorted[j];
      if( curtime >= hmax )
         break;

      /* subtract all capacity needed up to this point */
      freecapacity -= consdata->demands[startindices[j]];

      while( j+1 < nvars && startvaluessorted[j+1] == curtime )
      {
         ++j;
         freecapacity -= consdata->demands[startindices[j]];
      }

      /* free all capacity usages of jobs the are no longer running */
      while( endidx < nvars && curtime >= endvaluessorted[endidx] )
      {
         freecapacity += consdata->demands[endindices[endidx]];
         ++endidx;
      }

      assert(freecapacity <= consdata->capacity);
      assert(endidx <= nvars);

      /* --> endindex - points to the next job which will finish
       *     j        - points to the last job that has been released
       */


      /* check freecapacity to be smaller than zero
       * then we will add cover constraints to the MIP
       */
      if( freecapacity < 0 && curtime >= hmin )
      {
         int nextprofilechange;

         /* we can create covering constraints for each pint in time in interval [curtime; nextprofilechange[ */
         if( j < nvars-1 )
            nextprofilechange = MIN( startvaluessorted[j+1], endvaluessorted[endidx] );
         else
            nextprofilechange = endvaluessorted[endidx];

         nextprofilechange = MIN(nextprofilechange, hmax);

         for( t = curtime; t < nextprofilechange; ++t )
         {
            SCIPdebugMessage("add cover constraint for time %d\n", curtime);

            /* create covering constraint */
            SCIP_CALL( createCoverCutsTimepoint(scip, cons, startvalues, t)  );

         }
      } /* end if freecapacity > 0 */

   } /*lint --e{850}*/

   consdata->covercuts = TRUE;

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endvaluessorted);
   SCIPfreeBufferArray(scip, &startvaluessorted);
   SCIPfreeBufferArray(scip, &endvalues);
   SCIPfreeBufferArray(scip, &startvalues);

   return SCIP_OKAY;
}

/** this method creates a row for time point curtime which insures the capacity restriction of the cumulative
 *  constraint
 */
static
SCIP_RETCODE createCapacityRestriction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int                   curtime,            /**< current point in time */
   int                   nstarted,           /**< number of jobs that start before the curtime or at curtime */
   int                   nfinished,          /**< number of jobs that finished before curtime or at curtime */
   SCIP_Bool             cutsasconss         /**< should the cumulative constraint create the cuts as constraints? */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   int* coefs;
   int nbinvars;
   char name[SCIP_MAXSTRLEN];
   int capacity;
   int b;

   assert(nstarted > nfinished);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);

   capacity = consdata->capacity;
   assert(capacity > 0);

   nbinvars = 0;
   SCIP_CALL( collectBinaryVars(scip, consdata, &binvars, &coefs, &nbinvars, startindices, curtime, nstarted, nfinished) );

   /* construct row name */
   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d[%d]", SCIPconsGetName(cons), nstarted-1, curtime);

   if( cutsasconss )
   {
      SCIP_CONS* lincons;

      /* create linear constraint for the linking between the binary variables and the integer variable */
      SCIP_CALL( SCIPcreateConsKnapsack(scip, &lincons, name, 0, NULL, NULL, (SCIP_Longint)(capacity),
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );

      for( b = 0; b < nbinvars; ++b )
      {
         SCIP_CALL( SCIPaddCoefKnapsack(scip, lincons, binvars[b], (SCIP_Longint)coefs[b]) );
      }

      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   }
   else
   {
      SCIP_ROW* row;

      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), (SCIP_Real)capacity, FALSE, FALSE, SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      for( b = 0; b < nbinvars; ++b )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, binvars[b], (SCIP_Real)coefs[b]) );
      }

      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
      SCIPdebug( SCIP_CALL(SCIPprintRow(scip, row, NULL)) );

      if( consdata->demandrowssize == 0 )
      {
         consdata->demandrowssize = 10;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->demandrows, consdata->demandrowssize) );
      }
      if( consdata->ndemandrows == consdata->demandrowssize )
      {
         consdata->demandrowssize *= 2;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->demandrows, consdata->ndemandrows, consdata->demandrowssize) );
      }

      consdata->demandrows[consdata->ndemandrows] = row;
      consdata->ndemandrows++;
   }

   SCIPfreeBufferArrayNull(scip, &binvars);
   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** this method checks how many cumulatives can run at most at one time if this is greater than the capacity it creates
 *  row
 */
static
SCIP_RETCODE consCapacityConstraintsFinder(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   SCIP_Bool             cutsasconss         /**< should the cumulative constraint create the cuts as constraints? */
   )
{
   SCIP_CONSDATA* consdata;

   int* starttimes;         /* stores when each job is starting */
   int* endtimes;           /* stores when each job ends */
   int* startindices;       /* we will sort the startsolvalues, thus we need to know wich index of a job it corresponds to */
   int* endindices;         /* we will sort the endsolvalues, thus we need to know wich index of a job it corresponds to */

   int nvars;               /* number of activities for this constraint */
   int freecapacity;        /* remaining capacity */
   int curtime;             /* point in time which we are just checking */
   int endindex;            /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int hmin;
   int hmax;

   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* if no activities are associated with this cumulative then this constraint is redundant */
   if( nvars == 0 )
      return SCIP_OKAY;

   assert(consdata->vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   SCIPdebugMessage("create sorted event points for cumulative constraint <%s> with %d jobs\n",
      SCIPconsGetName(cons), nvars);

   /* create event point arrays */
   createSortedEventpoints(scip, nvars, consdata->vars, consdata->durations,
      starttimes, endtimes, startindices, endindices, FALSE);

   endindex = 0;
   freecapacity = consdata->capacity;
   hmin = consdata->hmin;
   hmax = consdata->hmax;

   /* check each startpoint of a job whether the capacity is kept or not */
   for( j = 0; j < nvars; ++j )
   {
      curtime = starttimes[j];
      SCIPdebugMessage("look at %d-th job with start %d\n", j, curtime);

      if( curtime >= hmax )
         break;

      /* remove the capacity requirments for all job which start at the curtime */
      subtractStartingJobDemands(consdata, curtime, starttimes, startindices, &freecapacity, &j, nvars);

      /* add the capacity requirments for all job which end at the curtime */
      addEndingJobDemands(consdata, curtime, endtimes, endindices, &freecapacity, &endindex, nvars);

      assert(freecapacity <= consdata->capacity);
      assert(endindex <= nvars);

      /* endindex - points to the next job which will finish */
      /* j - points to the last job that has been released */

      /* if free capacity is smaller than zero, then add rows to the LP */
      if( freecapacity < 0 && curtime >= hmin )
      {
         int nextstarttime;
         int t;

         /* step forward until next job is released and see whether capacity constraint is met or not */
         if( j < nvars-1 )
            nextstarttime = starttimes[j+1];
         else
            nextstarttime = endtimes[nvars-1];

         nextstarttime = MIN(nextstarttime, hmax);

         /* create capacity restriction row for current event point */
         SCIP_CALL( createCapacityRestriction(scip, cons, startindices, curtime, j+1, endindex, cutsasconss) );

         /* create for all points in time between the current event point and next start event point a row if the free
          * capacity is still smaller than zero  */
         for( t = curtime+1 ; t < nextstarttime; ++t )
         {
            /* add the capacity requirments for all job which end at the curtime */
            addEndingJobDemands(consdata, t, endtimes, endindices, &freecapacity, &endindex, nvars);

            if( freecapacity < 0 )
            {
               /* add constraint */
               SCIPdebugMessage("add capacity constraint at time %d\n", t);

               /* create capacity restriction row */
               SCIP_CALL( createCapacityRestriction(scip, cons, startindices, t, j+1, endindex, cutsasconss) );
            }
            else
               break;
         }
      }
   } /*lint --e{850}*/

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endtimes);
   SCIPfreeBufferArray(scip, &starttimes);

   return SCIP_OKAY;
}

/** creates LP rows corresponding to cumulative constraint; therefore, check each point in time if the maximal needed
 *  capacity is larger than the capacity of the cumulative constraint
 *  - for each necessary point in time:
 *
 *    sum_j sum_t demand_j * x_{j,t} <= capacity
 *
 *    where x(j,t) is the binary variables of job j at time t
 */
static
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   SCIP_Bool             cutsasconss         /**< should the cumulative constraint create the cuts as constraints? */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->demandrows == NULL);
   assert(consdata->ndemandrows == 0);

   /* collect the linking constraints */
   if( consdata->linkingconss == NULL )
   {
      SCIP_CALL( consdataCollectLinkingCons(scip, consdata) );
   }

   SCIP_CALL( consCapacityConstraintsFinder(scip, cons, cutsasconss) );

   /* switch of separation for the cumulative constraint if linear constraints are add as cuts */
   if( cutsasconss )
   {
      if( SCIPconsIsInitial(cons) )
      {
         SCIP_CALL( SCIPsetConsInitial(scip, cons, FALSE) );
      }
      if( SCIPconsIsSeparated(cons) )
      {
         SCIP_CALL( SCIPsetConsSeparated(scip, cons, FALSE) );
      }
      if( SCIPconsIsEnforced(cons) )
      {
         SCIP_CALL( SCIPsetConsEnforced(scip, cons, FALSE) );
      }
   }

   return SCIP_OKAY;
}

/** adds linear relaxation of cumulative constraint to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   SCIP_Bool             cutsasconss         /**< should the cumulative constraint create the cuts as constraints? */
   )
{
   SCIP_CONSDATA* consdata;
   int r;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->demandrows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons, cutsasconss) );
   }

   for( r = 0; r < consdata->ndemandrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->demandrows[r]) )
      {
         assert(consdata->demandrows[r] != NULL);
         SCIP_CALL( SCIPaddCut(scip, NULL, consdata->demandrows[r], FALSE) );
      }
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
SCIP_RETCODE separateConsBinaryRepresentation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint to be separated */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            separated           /**< pointer to store TRUE, if a cut was found */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int ncuts;
   int r;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("separate cumulative constraint <%s>\n", SCIPconsGetName(cons));

   if( consdata->demandrows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons, FALSE) );
   }

   ncuts = 0;

   /* check each row that is not contained in LP */
   for( r = 0; r < consdata->ndemandrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->demandrows[r]) )
      {
         SCIP_Real feasibility;

         if( sol != NULL )
            feasibility = SCIPgetRowSolFeasibility(scip, consdata->demandrows[r], sol);
         else
            feasibility = SCIPgetRowLPFeasibility(scip, consdata->demandrows[r]);

         if( SCIPisFeasNegative(scip, feasibility) )
         {
            SCIP_CALL( SCIPaddCut(scip, sol,  consdata->demandrows[r], FALSE) );
            ncuts++;
         }
      }
   }

   if( ncuts > 0 )
   {
      SCIPdebugMessage("cumulative constraint <%s> separated %d cuts\n",  SCIPconsGetName(cons), ncuts);

      /* if successful, reset age of constraint */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      (*separated) = TRUE;
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
SCIP_RETCODE separateCoverCutsCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint to be separated */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            separated           /**< pointer to store TRUE, if a cut was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_ROW* row;
   SCIP_Real minfeasibility;
   int r;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("separate cumulative constraint <%s>\n", SCIPconsGetName(cons));

   /* collect the linking constraints */
   if( consdata->linkingconss == NULL )
   {
      SCIP_CALL( consdataCollectLinkingCons(scip, consdata) );
   }

   if( !consdata->covercuts )
   {
      SCIP_CALL( createCoverCuts(scip, cons) );
   }

   row = NULL;
   minfeasibility = SCIPinfinity(scip);

   /* check each row of small covers that is not contained in LP */
   for( r = 0; r < consdata->nscoverrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->scoverrows[r]) )
      {
         SCIP_Real feasibility;

         assert(consdata->scoverrows[r] != NULL);
         if( sol != NULL )
            feasibility = SCIPgetRowSolFeasibility(scip, consdata->scoverrows[r], sol);
         else
            feasibility = SCIPgetRowLPFeasibility(scip, consdata->scoverrows[r]);

         if( minfeasibility > feasibility )
         {
            minfeasibility = feasibility;
            row =  consdata->scoverrows[r];
         }
      }
   }

   if( SCIPisFeasNegative(scip, minfeasibility) )
   {
      SCIPdebugMessage("cumulative constraint <%s> separated 1 cover cut with feasibility %g\n",
         SCIPconsGetName(cons), minfeasibility);

      assert(row != NULL);
      SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );

      /* if successful, reset age of constraint */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      (*separated) = TRUE;
   }

   minfeasibility = SCIPinfinity(scip);
   row = NULL;

   /* check each row of small covers that is not contained in LP */
   for( r = 0; r < consdata->nbcoverrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->bcoverrows[r]) )
      {
         SCIP_Real feasibility;

         assert(consdata->bcoverrows[r] != NULL);
         if( sol != NULL )
            feasibility = SCIPgetRowSolFeasibility(scip, consdata->bcoverrows[r], sol);
         else
            feasibility = SCIPgetRowLPFeasibility(scip, consdata->bcoverrows[r]);

         if( minfeasibility > feasibility )
         {
            minfeasibility = feasibility;
            row =  consdata->bcoverrows[r];
         }
      }
   }

   if( SCIPisFeasNegative(scip, minfeasibility) )
   {
      SCIPdebugMessage("cumulative constraint <%s> separated 1 cover cut with feasibility %g\n",
         SCIPconsGetName(cons), minfeasibility);

      assert(row != NULL);
      SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );

      /* if successful, reset age of constraint */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      (*separated) = TRUE;
   }

   return SCIP_OKAY;
}

/** this method creats a row for time point curtime which ensures the capacity restriction of the cumulative constraint */
static
SCIP_RETCODE createCapacityRestrictionIntvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int                   curtime,            /**< current point in time */
   int                   nstarted,           /**< number of jobs that start before the curtime or at curtime */
   int                   nfinished,          /**< number of jobs that finished before curtime or at curtime */
   SCIP_Bool             lower               /**< shall cuts be created due to lower or upper bounds? */
   )
{
   SCIP_CONSDATA* consdata;
   char name[SCIP_MAXSTRLEN];
   int lhs; /* left hand side of constraint */

   SCIP_VAR** activevars;
   SCIP_ROW* row;

   int v;

   assert(nstarted > nfinished);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);


   SCIP_CALL( SCIPallocBufferArray(scip, &activevars, nstarted-nfinished) );

   SCIP_CALL( collectIntVars(scip, consdata, &activevars, startindices, curtime, nstarted, nfinished, lower, &lhs ) );

   if( lower )
   {
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "lower(%d)", curtime);

      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, (SCIP_Real) lhs, SCIPinfinity(scip),  TRUE, FALSE, SCIPconsIsRemovable(cons)) );
   }
   else
   {
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "upper(%d)", curtime);
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), (SCIP_Real) lhs, TRUE, FALSE, SCIPconsIsRemovable(cons)) );
   }

   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for( v = 0; v < nstarted - nfinished; ++v )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, row, activevars[v], 1.) );
   }

   SCIP_CALL( SCIPflushRowExtensions(scip, row) );
   SCIPdebug( SCIP_CALL(SCIPprintRow(scip, row, NULL)) );

   SCIP_CALL( SCIPaddCut(scip, sol, row, TRUE) );

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   /* free buffers */
   SCIPfreeBufferArrayNull(scip, &activevars);

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
SCIP_RETCODE separateConsOnIntegerVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint to be separated */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool             lower,              /**< shall cuts be created according to lower bounds? */
   SCIP_Bool*            separated           /**< pointer to store TRUE, if a cut was found */
   )
{

   SCIP_CONSDATA* consdata;

   int* starttimes;         /* stores when each job is starting */
   int* endtimes;           /* stores when each job ends */
   int* startindices;       /* we will sort the startsolvalues, thus we need to know wich index of a job it corresponds to */
   int* endindices;         /* we will sort the endsolvalues, thus we need to know wich index of a job it corresponds to */

   int nvars;               /* number of activities for this constraint */
   int freecapacity;        /* remaining capacity */
   int curtime;             /* point in time which we are just checking */
   int endindex;            /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int hmin;
   int hmax;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* if no activities are associated with this cumulative then this constraint is redundant */
   if( nvars <= 1 )
      return SCIP_OKAY;

   assert(consdata->vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   SCIPdebugMessage("create sorted event points for cumulative constraint <%s> with %d jobs\n",
      SCIPconsGetName(cons), nvars);

   /* create event point arrays */
   createSelectedSortedEventpointsSol(scip, consdata, sol, starttimes, endtimes, startindices, endindices, &nvars, lower);

   /* now nvars might be smaller than before! */

   endindex = 0;
   freecapacity = consdata->capacity;
   hmin = consdata->hmin;
   hmax = consdata->hmax;

   /* check each startpoint of a job whether the capacity is kept or not */
   for( j = 0; j < nvars; ++j )
   {
      curtime = starttimes[j];

      if( curtime >= hmax )
         break;

      /* remove the capacity requirements for all job which start at the curtime */
      subtractStartingJobDemands(consdata, curtime, starttimes, startindices, &freecapacity, &j, nvars);

      /* add the capacity requirments for all job which end at the curtime */
      addEndingJobDemands(consdata, curtime, endtimes, endindices, &freecapacity, &endindex, nvars);

      assert(freecapacity <= consdata->capacity);
      assert(endindex <= nvars);

      /* endindex - points to the next job which will finish */
      /* j - points to the last job that has been released */

      /* if free capacity is smaller than zero, then add rows to the LP */
      if( freecapacity < 0 && curtime >= hmin)
      {
         /* create capacity restriction row for current event point */
         SCIP_CALL( createCapacityRestrictionIntvars(scip, cons, sol, startindices, curtime, j+1, endindex, lower) );
         *separated = TRUE;
      }
   } /*lint --e{850}*/

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endtimes);
   SCIPfreeBufferArray(scip, &starttimes);

   return SCIP_OKAY;
}

/**@} */


/**@name Presolving
 *
 * @{
 */

#ifndef NDEBUG
/** returns TRUE if all demands are smaller than the capacity of the cumulative constraint and if the total demand is
 *  correct
 */
static
SCIP_Bool checkDemands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to be checked */
   )
{
   SCIP_CONSDATA* consdata;
   int capacity;
   int nvars;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* if no activities are associated with this cumulative then this constraint is not infeasible, return */
   if( nvars == 0 )
      return TRUE;

   assert(consdata->vars != NULL);
   capacity = consdata->capacity;

   /* check each activity: if demand is larger than capacity the problem is infeasible */
   for ( j = 0; j < nvars; ++j )
   {
      if( consdata->demands[j] > capacity )
         return FALSE;
   }

   return TRUE;
}
#endif

/** delete constraint if it consists of at most one job */
static
SCIP_RETCODE deleteTrivilCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to propagate */
   int*                  ndelconss,          /**< pointer to store the number of deleted constraints */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nvars == 0 )
   {
      SCIPdebugMessage("delete cumulative constraints <%s>\n", SCIPconsGetName(cons));

      SCIP_CALL( SCIPdelCons(scip, cons) );
      (*ndelconss)++;
   }
   else if( consdata->nvars == 1 )
   {
      if( consdata->demands[0] > consdata->capacity )
         (*cutoff) = TRUE;
      else
      {
         SCIPdebugMessage("delete cumulative constraints <%s>\n", SCIPconsGetName(cons));

         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss)++;
      }
   }

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

/** adjust bounds of over sizeed job (the demand is larger than the capacity) */
static
SCIP_RETCODE adjustOversizedJobBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   pos,                /**< position of job in the consdata */
   int*                  nchgbds,            /**< pointer to store the number of changed bounds */
   int*                  naddconss,          /**< pointer to store the number of added constraints */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_VAR* var;
   SCIP_Bool tightened;
   int duration;
   int ect;
   int lst;

   assert(scip != NULL);

   /* zero energy jobs should be removed already */
   assert(consdata->durations[pos] > 0);
   assert(consdata->demands[pos] > 0);

   var = consdata->vars[pos];
   assert(var != NULL);
   duration =  consdata->durations[pos];

   /* jobs with a demand greater than the the capacity have to moved outside the time interval [hmin,hmax) */
   SCIPdebugMessage("  variable <%s>: demand <%d> is larger than the capacity <%d>\n",
      SCIPvarGetName(var), consdata->demands[pos], consdata->capacity);

   /* earliest completion time (ect) and latest start time (lst) */
   ect = convertBoundToInt(scip, SCIPvarGetLbGlobal(var)) + duration;
   lst = convertBoundToInt(scip, SCIPvarGetUbGlobal(var));

   if( ect > consdata->hmin && lst < consdata->hmax )
   {
      /* the job will at least run partly in the time interval [hmin,hmax) this means the problem is infeasible */
      *cutoff = TRUE;
   }
   else if( lst < consdata->hmax )
   {
      /* move the latest start time of this job in such a way that it finishes before or at hmin */
      SCIP_CALL( SCIPtightenVarUb(scip, var, (SCIP_Real)(consdata->hmin - duration), TRUE, cutoff, &tightened) );
      assert(tightened);
      assert(!(*cutoff));
      (*nchgbds)++;
   }
   else if( ect > consdata->hmin )
   {
      /* move the earliest start time of this job in such a way that it starts after or at hmax */
      SCIP_CALL( SCIPtightenVarLb(scip, var, (SCIP_Real)(consdata->hmax), TRUE, cutoff, &tightened) );
      assert(tightened);
      assert(!(*cutoff));
      (*nchgbds)++;
   }
   else
   {
      /* this job can run before or after the time interval [hmin,hmax) thus we create a bound disjunction
       * constraint to ensure that it does not overlap with the time interval [hmin,hmax); that is:
       *
       * (var <= hmin - duration) /\ (var >= hamx)
       */
      SCIP_CONS* cons;

      SCIP_VAR* vartuple[2];
      SCIP_BOUNDTYPE boundtypetuple[2];
      SCIP_Real boundtuple[2];

      char name[SCIP_MAXSTRLEN];
      int leftbound;
      int rightbound;

      leftbound = consdata->hmin - duration;
      rightbound = consdata->hmax;

      /* allocate temporary memory for arrays */
      vartuple[0] = var;
      vartuple[1] = var;
      boundtuple[0] = (SCIP_Real)leftbound;
      boundtuple[1] = (SCIP_Real)rightbound;
      boundtypetuple[0] = SCIP_BOUNDTYPE_UPPER;
      boundtypetuple[1] = SCIP_BOUNDTYPE_LOWER;

      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s<=%d or %s >= %d",
         SCIPvarGetName(var), leftbound, SCIPvarGetName(var), rightbound);

      /* creat bounddisjunction constraint */
      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, 2, vartuple, boundtypetuple, boundtuple,
            TRUE, FALSE, TRUE, TRUE /*check*/, TRUE/*prop*/, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIPdebugPrintCons(scip, cons, NULL);

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      (*naddconss)++;
   }

   return SCIP_OKAY;
}

/** try to removed over sizeed jobs (the demand is larger than the capacity) */
static
SCIP_RETCODE removeOversizedJobs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int*                  nchgbds,            /**< pointer to store the number of changed bounds */
   int*                  nchgcoefs,          /**< pointer to store the number of changed coefficient */
   int*                  naddconss,          /**< pointer to store the number of added constraints */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_CONSDATA* consdata;
   int capacity;
   int j;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* if a cutoff was already detected just return */
   if( *cutoff )
      return SCIP_OKAY;

   if( consdata->nvars == 0 )
      return SCIP_OKAY;

   capacity = consdata->capacity;

   for( j = consdata->nvars-1; j >= 0 && !(*cutoff);  --j )
   {
      if( consdata->demands[j] > capacity )
      {
         SCIP_CALL( adjustOversizedJobBounds(scip, consdata, j, nchgbds, naddconss, cutoff) );

         /* remove variable form constraint */
         SCIP_CALL( consdataDeletePos(scip, consdata, cons, j) );
         (*nchgcoefs)++;
      }
   }

   SCIPdebugMessage("cumulative constraint <%s> has %d jobs left, cutoff %u\n", SCIPconsGetName(cons), consdata->nvars, *cutoff);

   return SCIP_OKAY;
}

/** fix integer variable to upper bound if the rounding locks and the object coefficient are in favor of that */
static
SCIP_RETCODE fixIntegerVariableUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< integer variable to fix */
   SCIP_Bool             uplock,             /**< has thet start time variable a up lock */
   int*                  nfixedvars          /**< pointer to store the number fixed variables */
   )
{
   SCIP_Real objval;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   /* if SCIP is in probing mode we cannot perform this dual reductions since this dual reduction would end in an
    * implication which can lead to an infeasible cutoff (optimal solution)
    */
   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* rounding the variable to the upper bound is only a feasible dual reduction if the cumulative constraint
    * handler is the only one locking that variable up
    */
   assert(uplock == TRUE || uplock == FALSE);
   assert((int)TRUE == 1);
   assert((int)FALSE == 0);

   if( SCIPvarGetNLocksUp(var) > (int)(uplock) )
      return SCIP_OKAY;

   objval = SCIPvarGetObj(var);

   /* rounding the integer variable up is only a valid dual reduction if the object coefficient is zero or negative
    * (the transformed problem is always a minimization problem)
    */
   if( SCIPisPositive(scip, objval) )
      return SCIP_OKAY;

   SCIPdebugMessage("try fixing variable <%s>[%g,%g] to upper bound %g\n", SCIPvarGetName(var),
      SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), SCIPvarGetUbLocal(var));

   SCIP_CALL( SCIPfixVar(scip, var, SCIPvarGetUbLocal(var), &infeasible, &tightened) );
   assert(!infeasible);

   if( tightened )
   {
      SCIPdebugMessage("fix variable <%s> to upper bound %g\n", SCIPvarGetName(var), SCIPvarGetUbLocal(var));
      (*nfixedvars)++;
   }

   return SCIP_OKAY;
}

/** fix integer variable to lower bound if the rounding locks and the object coefficient are in favor of that */
static
SCIP_RETCODE fixIntegerVariableLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< integer variable to fix */
   SCIP_Bool             downlock,           /**< has the variable a down lock */
   int*                  nfixedvars          /**< pointer to store the number fixed variables */
   )
{
   SCIP_Real objval;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   /* if SCIP is in probing mode we cannot perform this dual reductions since this dual reduction would end in an
    * implication which can lead to cutoff the optimal solution
    */
   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* rounding the variable to the lower bound is only a feasible dual reduction if the cumulative constraint
    * handler is the only one locking that variable down
    */
   assert(downlock == TRUE || downlock == FALSE);
   assert((int)TRUE == 1);
   assert((int)FALSE == 0);

   if( SCIPvarGetNLocksDown(var) > (int)(downlock) )
      return SCIP_OKAY;

   objval = SCIPvarGetObj(var);

   /* rounding the integer variable down is only a valid dual reduction if the object coefficient is zero or positive
    * (the transformed problem is always a minimization problem)
    */
   if( SCIPisNegative(scip, objval) )
      return SCIP_OKAY;


   SCIP_CALL( SCIPfixVar(scip, var, SCIPvarGetLbLocal(var), &infeasible, &tightened) );
   assert(!infeasible);

   if( tightened )
   {
      SCIPdebugMessage("fix variable <%s> to lower bound %g\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var));
      (*nfixedvars)++;
   }

   return SCIP_OKAY;
}

/** divides demands by their greatest common divisor and divides capacity by the same value, rounding down the result;
 *  in case the the smallest demands add up to more than the capacity we reductions all demands to one as well as the
 *  capacity since in that case none of the jobs can run in parallel
 */
static
SCIP_RETCODE normalizeDemands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Longint gcd;
   int capacity;
   int mindemand1;
   int mindemand2;
   int nvars;
   int v;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(!SCIPconsIsModifiable(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   capacity = consdata->capacity;
   nvars = consdata->nvars;

   if( consdata->normalized || capacity == 1 || nvars <= 1 )
      return SCIP_OKAY;

   /**@todo sort items w.r.t. the demands, because we can stop earlier if the smaller weights are evaluated first */

   gcd = (SCIP_Longint)consdata->demands[nvars-1];
   mindemand1 = MIN(consdata->demands[nvars-1], consdata->demands[nvars-2]);
   mindemand2 = MAX(consdata->demands[nvars-1], consdata->demands[nvars-2]);

   for( v = nvars-2; v >= 0 && (gcd >= 2 || mindemand1 + mindemand2 > capacity); --v )
   {
      assert(mindemand1 <= mindemand2);
      gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)consdata->demands[v]);

      if( mindemand1 > consdata->demands[v] )
      {
         mindemand2 = mindemand1;
         mindemand1 = consdata->demands[v];
      }
      else if( mindemand2 > consdata->demands[v] )
         mindemand2 = consdata->demands[v];
   }

   if( mindemand1 + mindemand2 > capacity )
   {
      SCIPdebugMessage("update cumulative constraint <%s> (%d + %d > %d) to unary cumulative constraint\n", SCIPconsGetName(cons), mindemand1, mindemand2, capacity);

      for( v = 0; v < nvars; ++v )
         consdata->demands[v] = 1;

      consdata->capacity = 1;

      (*nchgcoefs) += nvars;
      (*nchgsides)++;
   }
   else if( gcd >= 2 )
   {
      SCIPdebugMessage("cumulative constraint <%s>: dividing demands by %"SCIP_LONGINT_FORMAT"\n", SCIPconsGetName(cons), gcd);

      for( v = 0; v < nvars; ++v )
         consdata->demands[v] /= gcd;

      consdata->capacity /= gcd;

      (*nchgcoefs) += nvars;
      (*nchgsides)++;
   }

   consdata->normalized = TRUE;

   return SCIP_OKAY;
}


/** computes the effective horizon and checks if the constraint can be decompesd */
static
SCIP_RETCODE computeEffectiveHorizon(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   int*                  ndelconss,          /**< pointer to store the number of deleted constraints */
   int*                  naddconss,          /**< pointer to store the number of added constraints */
   int*                  nchgsides           /**< pointer to store the number of changed sides */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_PROFILE* profile;
   int capacity;
   int hmin;
   int hmax;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nvars <= 1 )
      return SCIP_OKAY;

   capacity = consdata->capacity;

   /* create empty resource profile with infinity resource capacity */
   SCIP_CALL( SCIPprofileCreate(&profile, INT_MAX) );

   /* create worst case resource profile */
   SCIP_CALL( createWorstCaseProfile(scip, profile, consdata) );

   /* print resource profile in if SCIP_DEBUG is defined */
   SCIPdebug( SCIPprofilePrint(profile, SCIPgetMessagehdlr(scip), NULL) );

   /* computes the first time point where the resource capacity can be violated */
   hmin = computeHmin(scip, profile, capacity);

   /* check if this time point improves the effective horizon */
   if( consdata->hmin < hmin )
   {
      SCIPdebugMessage("cumulative constraint <%s> adjust hmin <%d> -> <%d>\n", SCIPconsGetName(cons), consdata->hmin, hmin);

      consdata->hmin = hmin;
      (*nchgsides)++;
   }

   /* computes the first time point where the resource capacity is satisfied for sure */
   hmax = computeHmax(scip, profile, capacity);

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
   }
   else if( !SCIPinRepropagation(scip) )
   {
      int* timepoints;
      int* loads;
      int ntimepoints;
      int t;

      /* If SCIP is repropagating the root node, it is not possible to decompose the constraints. This is the case since
       * the conflict analysis stores the constraint pointer for bound changes made by this constraint. These pointer
       * are used during the resolve propagation phase to explain bound changes. If we would decompose certain jobs into
       * a new cumulative constraint, the "old" pointer is not valid. More precise, the "old" constraint is not able to
       * explain the certain "old" bound changes
       */

      /* search for time points */
      ntimepoints = SCIPprofileGetNTimepoints(profile);
      timepoints = SCIPprofileGetTimepoints(profile);
      loads = SCIPprofileGetLoads(profile);

      /* check if there exist a time point within the effective horizon [hmin,hmax) such that the capacity is not exceed w.r.t. worst case profile */
      for( t = 0; t < ntimepoints; ++t )
      {
         /* ignore all time points before the effective horizon */
         if( timepoints[t] <= consdata->hmin )
            continue;

         /* ignore all time points after the effective horizon */
         if( timepoints[t] >= consdata->hmax )
            break;

         /* check if the current time point does not exceed the capacity w.r.t. worst case resource profile; if so we
          * can split the cumulative constraint into two cumulative constraints
          */
         if( loads[t] <= capacity )
         {
            SCIP_CONS* splitcons;
            char name[SCIP_MAXSTRLEN];

            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "(%s)'", SCIPconsGetName(cons));

            SCIPdebugMessage("split cumulative constraint <%s>[%d,%d) with %d jobs at time point %d\n",
               SCIPconsGetName(cons), consdata->hmin, consdata->hmax, consdata->nvars, timepoints[t]);

            SCIP_CALL( SCIPcreateConsCumulative(scip, &splitcons, name, consdata->nvars, consdata->vars,
                  consdata->durations, consdata->demands, capacity,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
                  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

            /* adjust the effective time horizon of the new constraint */
            SCIP_CALL( SCIPsetHminCumulative(scip, splitcons, timepoints[t]) );
            SCIP_CALL( SCIPsetHmaxCumulative(scip, splitcons, consdata->hmax) );

            assert(timepoints[t] < consdata->hmax);

            /* add and release new cumulative constraint */
            SCIP_CALL( SCIPaddCons(scip, splitcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &splitcons) );

            /* adjust the effective time horizon of the constraint */
            consdata->hmax = timepoints[t];

            assert(consdata->hmin < consdata->hmax);

            SCIPstatistic( consdata->ndecomps++ );
            (*naddconss)++;
            break;
         }
      }
   }

   /* free worst case profile */
   SCIPprofileFree(&profile);

   return SCIP_OKAY;
}

/** presolve constraint w.r.t. the earlier start times (est)
 *
 *  (1) Let the variables array be non-increasing ordered by there earlier start time and j1 and j2 the last two
 *      variables of w.r.t. that ordering. If the latest completion time (lct) of job j1 is less than or equal to the
 *      earlier start time (est) of job j2, then job j1 can be scheduled at any time of the feasible time window without
 *      interfering with other jobs
 *      => j1 can be removed from the cumulative constraint
 *
 *  (2) if the earliest completion time (ect) of job j1 is less than or equal to the earliest start time (est) of job j2
 *       and fixing the start time variable of job j1 to the lower bound is a feasible dual reduction
 *       => j1 can be removed from the cumulative constraint and fixed to its earliest start time
 */
static
SCIP_RETCODE presolveConsEst(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   int*                  nchgcoefs,          /**< pointer to store the number of changed coefficients */
   int*                  nchgsides,          /**< pointer to store the number of changed sides */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* cand;
   SCIP_Bool downlock;
   SCIP_Bool uplock;
   int duration;
   int nvars;
   int hmin;
   int est;
   int lst;
   int ect;
   int lct;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* check if a cutoff was already detected */
   if( *cutoff )
      return SCIP_OKAY;

   if( nvars <= 1 )
      return SCIP_OKAY;

   hmin = consdata->hmin;
   assert(hmin < consdata->hmax);

   SCIPdebugMessage("presolve constraint <%s>[%d,%d) <= %d w.r.t. earlier start time\n", SCIPconsGetName(cons), hmin, consdata->hmax, consdata->capacity);

   /* make sure, the variables are sorted by non-increasing order of the earliest start times */
   SCIP_CALL( sortVariablesEst(scip, cons) );

   for( v = nvars-1; v >= 0 && nvars > 1; --v )
   {
#ifndef NDEBUG
      checkSortVariablesEst(scip, consdata, v);
#endif

      cand = consdata->vars[v];
      duration = consdata->durations[v];
      downlock = consdata->downlocks[v];
      uplock = consdata->uplocks[v];

      /* collect latest completion time (lct) and earliest completion time (ect) from last job */
      est = convertBoundToInt(scip, SCIPvarGetLbGlobal(cand));
      ect = est + duration;
      lst = convertBoundToInt(scip, SCIPvarGetUbGlobal(cand));
      lct = lst + duration;

      /* ???????????? second condition ?????? */
      if( lct <= hmin || (v > 0 && v == nvars-1 && lct <= convertBoundToInt(scip, SCIPvarGetLbGlobal(consdata->vars[v-1]))) )
      {
         SCIPdebugMessage("  remove variable <%s>[%g,%g] with duration <%d> is irrelevant\n",
            SCIPvarGetName(cand), SCIPvarGetLbGlobal(cand), SCIPvarGetUbGlobal(cand), duration);

         /* delete variable at the given position */
         SCIP_CALL( consdataDeletePos(scip, consdata, cons, v) );
         (*nchgcoefs)++;

         /* adjust nvars after deleting the variable */
         nvars = consdata->nvars;

         SCIPstatistic( consdata->nirrelevantjobs++ );
      }
      else if( ect <= hmin )
      {
         /* job can be removed from the constraint only if the integer start time variable can be fixed to its lower
          * bound;
          */

         /* fix integer start time variable if possible to it lower bound */
         SCIP_CALL( fixIntegerVariableLb(scip, cand, downlock, nfixedvars) );

         if( SCIPvarGetLbGlobal(cand) + 0.5 > SCIPvarGetUbGlobal(cand) )
         {
            SCIPdebugMessage("  remove variable <%s>[%d,%d] with duration <%d> is irrelevant due to dual reductions wrt EST\n",
               SCIPvarGetName(cand), est, lct, duration);

            SCIP_CALL( consdataDeletePos(scip, consdata, cons, v) );
            (*nchgcoefs)++;

            /* adjust nvars after deleting the variable */
            nvars = consdata->nvars;
            SCIPstatistic( consdata->ndualfixs++ );
         }
      }
      else if( lst <= hmin )
      {
         if( !uplock )
         {
            /* the variables has no up lock and we can also remove the down lock;
             * => lst <= hmin and ect >= hmax
             * => remove job and reduce capacity by the demand of that job
             */

            assert(ect >= consdata->hmax);

            /* if the capacity is smaller than the demand we detected an infeasibility */
            if( consdata->capacity < consdata->demands[v] )
            {
               (*cutoff) = TRUE;
               break;
            }

            consdata->capacity -= consdata->demands[v];
            assert(consdata->capacity >= 0);

            (*nchgsides)++;

            SCIPdebugMessage("  remove variables <%s>[%d,%d] with duration <%d> due to no uplocks, new capacity = %d\n",
               SCIPvarGetName(cand), est, lst, duration, consdata->capacity);

            SCIP_CALL( consdataDeletePos(scip, consdata, cons, v) );
            (*nchgcoefs)++;

            /* adjust nvars after deleting the variable */
            nvars = consdata->nvars;
            SCIPstatistic( consdata->nalwaysruns++ );
         }
         else if( downlock )
         {
            SCIPdebugMessage("  remove down lock of variable <%s>[%g,%g] with duration <%d>\n",
               SCIPvarGetName(cand), SCIPvarGetLbGlobal(cand), SCIPvarGetUbGlobal(cand), consdata->durations[v]);

            SCIP_CALL( SCIPunlockVarCons(scip, cand, cons, TRUE, FALSE) );
            consdata->downlocks[v] = FALSE;
            (*nchgsides)++;

            SCIPstatistic( consdata->nremovedlocks++ );
         }
      }
      else if ( est >= hmin )
         break;
   }

   return SCIP_OKAY;
}

/** presolve constraint w.r.t. the latest completion time (lct)
 *
 *  (1) Let the variables array be non-decreasing ordered by there latest completion time and j1 and j2 the last two
 *      variables w.r.t. that ordering. If the latest completion time (lct) of job j1 is less than or equal to the
 *      earlier start time (est) of job j2, then job j2 can be scheduled at any time of the feasible time window without
 *      interfering with other jobs
 *       => j2 can be removed from the cumulative constraint
 *
 *  (2) if the latest start time (lst) of job j1 is greater than or equal to the latest completion time (lct) of job j2
 *      and fixing the start time variable of job j1 to the upper bound is a feasible dual reduction
 *      => j1 can be removed from the cumulative constraint and fixed to its latest start time
 */
static
SCIP_RETCODE presolveConsLct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   int*                  nchgcoefs,          /**< pointer to store the number of changed coefficients */
   int*                  nchgsides,          /**< pointer to store the number of changed sides */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* cand;
   SCIP_Bool downlock;
   SCIP_Bool uplock;
   int duration;
   int nvars;
   int hmin;
   int hmax;
   int est;
   int lst;
   int ect;
   int lct;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check if a cutoff was already detected */
   if( *cutoff )
      return SCIP_OKAY;

   nvars = consdata->nvars;

   if( nvars <= 1 )
      return SCIP_OKAY;

   hmin = consdata->hmin;
   hmax = consdata->hmax;
   assert(hmin < hmax);

   SCIPdebugMessage("presolve constraint <%s>[%d,%d) <= %d w.r.t. latest completion time\n", SCIPconsGetName(cons), hmin, hmax, consdata->capacity);

   /* make sure, the variables are sorted by non-increasing order of the latest completion time */
   SCIP_CALL( sortVariablesLct(scip, cons) );

   for( v = nvars-1; v >= 0 && nvars > 1; --v )
   {
#ifndef NDEBUG
      checkSortVariablesLct(scip, consdata, v);
#endif

      cand = consdata->vars[v];
      duration = consdata->durations[v];
      downlock = consdata->downlocks[v];
      uplock = consdata->uplocks[v];

      /* collect latest completion time (lct) and earliest completion time (ect) from last job */
      est = convertBoundToInt(scip, SCIPvarGetLbGlobal(cand));
      ect = est + duration;
      lst = convertBoundToInt(scip, SCIPvarGetUbGlobal(cand));
      lct = lst + duration;

      if( est >= hmax || (v > 0 && v == nvars-1 && est >= convertBoundToInt(scip, SCIPvarGetUbGlobal(consdata->vars[v-1])) + consdata->durations[v-1]) )
      {
         SCIPdebugMessage("  remove variable <%s>[%g,%g] with duration <%d> is irrelevant\n",
            SCIPvarGetName(cand), SCIPvarGetLbGlobal(cand), SCIPvarGetUbGlobal(cand), duration);

         /* delete variable at the given position */
         SCIP_CALL( consdataDeletePos(scip, consdata, cons, v) );
         (*nchgcoefs)++;

         /* adjust nvars after deleting the variable */
         nvars = consdata->nvars;

         SCIPstatistic( consdata->nirrelevantjobs++ );
      }
      else if( lst >= hmax )
      {
         /* job can be removed from the constraint only if the integer start time variable can be fixed to its upper
          * bound
          */

         /* fix integer start time variable if possible to its upper bound */
         SCIP_CALL( fixIntegerVariableUb(scip, cand, uplock, nfixedvars) );

         if( SCIPvarGetLbGlobal(cand) + 0.5 > SCIPvarGetUbGlobal(cand) )
         {
            SCIPdebugMessage("  remove variable <%s>[%d,%d] with duration <%d> is irrelevant due to dual reductions wrt LCT\n",
               SCIPvarGetName(cand), est, lst, duration);

            SCIP_CALL( consdataDeletePos(scip, consdata, cons, v) );
            (*nchgcoefs)++;

            /* adjust nvars after deleting the variable */
            nvars = consdata->nvars;
            SCIPstatistic( consdata->ndualfixs++ );
         }
      }
      else if( ect >= hmax )
      {
         if( !downlock )
         {
            /* the variables has no down lock and we can also remove the up lock;
             * => lst <= hmin and ect >= hmax
             * => remove job and reduce capacity by the demand of that job
             */

            assert(lst <= consdata->hmin);

            /* if the capacity is smaller than the demand we detected an infeasibility */
            if( consdata->capacity < consdata->demands[v] )
            {
               (*cutoff) = TRUE;
               break;
            }

            consdata->capacity -= consdata->demands[v];
            assert(consdata->capacity >= 0);

            (*nchgsides)++;

            SCIPdebugMessage("  remove variables <%s>[%d,%d] with duration <%d> due to no downlocks, new capacity = %d\n",
               SCIPvarGetName(cand), est, lst, duration, consdata->capacity);

            SCIP_CALL( consdataDeletePos(scip, consdata, cons, v) );
            (*nchgcoefs)++;

            /* adjust nvars after deleting the variable */
            nvars = consdata->nvars;

            SCIPstatistic( consdata->nalwaysruns++ );
         }
         else if( uplock )
         {
            SCIPdebugMessage("  remove up lock of variable <%s>[%g,%g] with duration <%d>\n",
               SCIPvarGetName(cand), SCIPvarGetLbGlobal(cand), SCIPvarGetUbGlobal(cand), consdata->durations[v]);

            SCIP_CALL( SCIPunlockVarCons(scip, cand, cons, FALSE, TRUE) );
            consdata->uplocks[v] = FALSE;
            (*nchgsides)++;

            SCIPstatistic( consdata->nremovedlocks++ );
         }
      }
      else if( lct <= hmax )
         break;
   }

   return SCIP_OKAY;
}

/** stores all demands which are smaller than the capacity of those jobs that are running at 'curtime' */
static
SCIP_RETCODE collectDemands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int                   curtime,            /**< current point in time */
   int                   nstarted,           /**< number of jobs that start before the curtime or at curtime */
   int                   nfinished,          /**< number of jobs that finished before curtime or at curtime */
   SCIP_Longint**        demands,            /**< pointer to array storing the demands */
   int*                  ndemands            /**< pointer to store the number of different demands */
   )
{
   int startindex;
   int ncountedvars;

   assert(demands != NULL);
   assert(ndemands != NULL);

   ncountedvars = 0;
   startindex = nstarted - 1;

   *ndemands = 0;

   /* search for the (nstarted - nfinished) jobs which are active at curtime */
   while( nstarted - nfinished > ncountedvars )
   {
      SCIP_VAR* var;
      int endtime;
      int varidx;

      /* collect job information */
      varidx = startindices[startindex];
      assert(varidx >= 0 && varidx < consdata->nvars);

      var = consdata->vars[varidx];
      assert(var != NULL);

      endtime = convertBoundToInt(scip, SCIPvarGetUbGlobal(var)) + consdata->durations[varidx];

      /* check the end time of this job is larger than the curtime; in this case the job is still running */
      if( endtime > curtime )
      {
         if( consdata->demands[varidx] < consdata->capacity )
         {
            (*demands)[*ndemands] = consdata->demands[varidx];
            (*ndemands)++;
         }
         ncountedvars++;
      }

      startindex--;
   }

   return SCIP_OKAY;
}

/** this method creates a row for time point curtime which insures the capacity restriction of the cumulative
 *  constraint
 */
static
SCIP_RETCODE getHighestCapacityUsage(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   int*                  startindices,       /**< permutation with rspect to the start times */
   int                   curtime,            /**< current point in time */
   int                   nstarted,           /**< number of jobs that start before the curtime or at curtime */
   int                   nfinished,          /**< number of jobs that finished before curtime or at curtime */
   int*                  bestcapacity        /**< pointer to store the maximum possible capacity usage */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Longint* demands;
   SCIP_Real* profits;
   int* items;
   int ndemands;
   SCIP_Bool success;
   SCIP_Real solval;
   int j;
   assert(nstarted > nfinished);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);
   assert(consdata->capacity > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &demands, consdata->nvars) );
   ndemands = 0;

   /* get demand array to initialize knapsack problem */
   SCIP_CALL( collectDemands(scip, consdata, startindices, curtime, nstarted, nfinished, &demands, &ndemands) );

   /* create array for profits */
   SCIP_CALL( SCIPallocBufferArray(scip, &profits, ndemands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &items, ndemands) );
   for( j = 0; j < ndemands; ++j )
   {
      profits[j] = (SCIP_Real) demands[j];
      items[j] = j;/* this is only a dummy value*/
   }

   /* solve knapsack problem and get maximum capacity usage <= capacity */
   SCIP_CALL( SCIPsolveKnapsackExactly(scip, ndemands, demands, profits, (SCIP_Longint)consdata->capacity,
         items, NULL, NULL, NULL, NULL, &solval, &success) );

   assert(SCIPisFeasIntegral(scip, solval));

   /* store result */
   *bestcapacity = convertBoundToInt(scip, solval);

   SCIPfreeBufferArray(scip, &items);
   SCIPfreeBufferArray(scip, &profits);
   SCIPfreeBufferArray(scip, &demands);

   return SCIP_OKAY;
}

/** try to tighten the capacity
 *  -- using DP for knapsack, we find the maximum possible capacity usage
 *  -- neglects hmin and hmax, such that it is also able to check solutions globally
 */
static
SCIP_RETCODE tightenCapacity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides           /**< pointer to store the number of changed sides */
   )
{
   SCIP_CONSDATA* consdata;
   int* starttimes;         /* stores when each job is starting */
   int* endtimes;           /* stores when each job ends */
   int* startindices;       /* we will sort the startsolvalues, thus we need to know wich index of a job it corresponds to */
   int* endindices;         /* we will sort the endsolvalues, thus we need to know wich index of a job it corresponds to */

   int nvars;               /* number of activities for this constraint */
   int freecapacity;        /* remaining capacity */
   int curtime;             /* point in time which we are just checking */
   int endindex;            /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int bestcapacity;

   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nchgsides != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* if no activities are associated with this cumulative or the capacity is 1, then this constraint is redundant */
   if( nvars <= 1 || consdata->capacity <= 1 )
      return SCIP_OKAY;

   assert(consdata->vars != NULL);

   SCIPdebugMessage("try to tighten capacity for cumulative constraint <%s> with capacity %d\n",
      SCIPconsGetName(cons), consdata->capacity);

   SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   /* create event point arrays */
   createSortedEventpoints(scip, nvars, consdata->vars, consdata->durations,
      starttimes, endtimes, startindices, endindices, FALSE);

   bestcapacity = 1;
   endindex = 0;
   freecapacity = consdata->capacity;

   /* check each startpoint of a job whether the capacity is kept or not */
   for( j = 0; j < nvars && bestcapacity < consdata->capacity; ++j )
   {
      curtime = starttimes[j];
      SCIPdebugMessage("look at %d-th job with start %d\n", j, curtime);

      /* remove the capacity requirments for all job which start at the curtime */
      subtractStartingJobDemands(consdata, curtime, starttimes, startindices, &freecapacity, &j, nvars);

      /* add the capacity requirments for all job which end at the curtime */
      addEndingJobDemands(consdata, curtime, endtimes, endindices, &freecapacity, &endindex, nvars);

      assert(freecapacity <= consdata->capacity);
      assert(endindex <= nvars);

      /* endindex - points to the next job which will finish */
      /* j - points to the last job that has been released */

      /* check point in time when capacity is exceeded (here, a knapsack problem must be solved) */
      if( freecapacity < 0 )
      {
         int newcapacity;

         newcapacity = 1;

         /* get best possible upper bound on capacity usage */
         SCIP_CALL( getHighestCapacityUsage(scip, cons, startindices, curtime, j+1, endindex, &newcapacity) );

         /* update bestcapacity */
         bestcapacity = MAX(bestcapacity, newcapacity);
         SCIPdebugMessage("after highest cap usage: bestcapacity = %d\n", bestcapacity);
      }

      /* also those points in time, where the capacity limit is not exceeded, must be taken into account */
      if( freecapacity > 0 && freecapacity != consdata->capacity )
      {
         bestcapacity = MAX(bestcapacity, consdata->capacity - freecapacity);
         SCIPdebugMessage("after peak < cap: bestcapacity = %d\n", bestcapacity);
      }

      /* capacity cannot be decreased if the demand sum over more than one job equals the capacity */
      if( freecapacity == 0 && consdata->demands[startindices[j]] < consdata->capacity)
      {
         /* if demands[startindices[j]] == cap then exactly that job is running */
         SCIPdebugMessage("--> cannot decrease capacity since sum equals capacity\n");
         bestcapacity = consdata->capacity;
         break;
      }
   }  /*lint --e{850}*/

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endtimes);
   SCIPfreeBufferArray(scip, &starttimes);

   /* check whether capacity can be tightened and whether demands need to be adjusted */
   if( bestcapacity < consdata->capacity )
   {
      int oldnchgcoefs;

      oldnchgcoefs = *nchgcoefs;

      SCIPdebugMessage("+-+-+-+-+-+ --> CHANGE capacity of cons<%s> from %d to %d",
         SCIPconsGetName(cons), consdata->capacity, bestcapacity);

      for( j = 0; j < nvars; ++j )
      {
         if( consdata->demands[j] == consdata->capacity )
         {
            consdata->demands[j] = bestcapacity;
            (*nchgcoefs)++;
         }
      }

      consdata->capacity = bestcapacity;
      (*nchgsides)++;

      SCIPdebugPrintf("; changed additionally %d coefficients\n", (*nchgcoefs)-oldnchgcoefs);

   }

   return SCIP_OKAY;
}

/** tries to change coefficients:
 *        demand_j < cap && all other parallel jobs in conflict
 *         ==> set demand_j := cap
 */
static
SCIP_RETCODE tightenCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   int*                  nchgcoefs           /**< pointer to count total number of changed coefficients */
   )
{
   SCIP_CONSDATA* consdata;
   int nvars;
   int j;
   int oldnchgcoefs;
   int mindemand;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nchgcoefs != NULL);

   /* get constraint data for some parameter testings only! */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   oldnchgcoefs = *nchgcoefs;

   if( nvars <= 0 )
      return SCIP_OKAY;

   /* PRE1:
    * check all jobs j whether: r_j + r_min > capacity holds
    * if so: adjust r_j to capacity
    */
   mindemand = consdata->demands[0];
   for( j = 0; j < nvars; ++j )
   {
      mindemand = MIN(mindemand, consdata->demands[j]);
   }

   /*check each job */
   for( j = 0; j < nvars; ++j )
   {
      if( mindemand + consdata->demands[j] > consdata->capacity  && consdata->demands[j] < consdata->capacity )
      {
         SCIPdebugMessage("+-+-+-+-+-+change demand of var<%s> from %d to capacity %d\n", SCIPvarGetName(consdata->vars[j]),
            consdata->demands[j], consdata->capacity);
         consdata->demands[j] = consdata->capacity;
         (*nchgcoefs)++;
      }
   }

   /* PRE2:
    * check for each job (with d_j < cap)
    * whether it is disjunctive to all others over the time horizon
    */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_Bool chgcoef;
      int est_j;
      int lct_j;
      int i;

      assert(consdata->demands[j] <= consdata->capacity);

      if( consdata->demands[j] == consdata->capacity )
         continue;

      chgcoef = TRUE;

      est_j = convertBoundToInt(scip, SCIPvarGetLbLocal(consdata->vars[j]));
      lct_j = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[j])) + consdata->durations[j];

      for( i = 0; i < nvars; ++i )
      {
         int est_i;
         int lct_i;

         if( i == j )
            continue;

         est_i = convertBoundToInt(scip, SCIPvarGetLbLocal(consdata->vars[i]));
         lct_i = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[i])) + consdata->durations[i];

         if( est_i >= lct_j || est_j >= lct_i )
            continue;

         if( consdata->demands[j] + consdata->demands[i] <= consdata->capacity )
         {
            chgcoef = FALSE;
            break;
         }
      }

      if( chgcoef )
      {
         SCIPdebugMessage("+-+-+-+-+-+change demand of var<%s> from %d to capacity %d\n", SCIPvarGetName(consdata->vars[j]),
            consdata->demands[j], consdata->capacity);
         consdata->demands[j] = consdata->capacity;
         (*nchgcoefs)++;
      }

   }

   if( (*nchgcoefs) > oldnchgcoefs )
   {
      SCIPdebugMessage("+-+-+-+-+-+changed %d coefficients of variables of cumulative constraint<%s>\n",
         (*nchgcoefs) - oldnchgcoefs, SCIPconsGetName(cons));
   }

   return SCIP_OKAY;
}

#if 0
/** try to reformulate constraint by replacing certain jobs */
static
SCIP_RETCODE reformulateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   int*                  naggrvars           /**< pointer to store the number of aggregated variables */
   )
{
   SCIP_CONSDATA* consdata;
   int hmin;
   int hmax;
   int nvars;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(cons != NULL);

   nvars = consdata->nvars;
   assert(nvars > 1);

   hmin = consdata->hmin;
   hmax = consdata->hmax;
   assert(hmin < hmax);

   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      int duration;
      int est;
      int ect;
      int lst;
      int lct;

      var = consdata->vars[v];
      assert(var != NULL);

      duration = consdata->durations[v];

      est = convertBoundToInt(scip, SCIPvarGetLbGlobal(var));
      ect = est + duration;
      lst = convertBoundToInt(scip, SCIPvarGetUbGlobal(var));
      lct = lst + duration;

      /* jobs for which the core [lst,ect) contains [hmin,hmax) should be removed already */
      assert(lst > hmin || ect < hmax);

      if( lst <= hmin && est < hmin - lct + MIN(hmin, ect) )
      {
         SCIP_VAR* aggrvar;
         char name[SCIP_MAXSTRLEN];
         SCIP_Bool infeasible;
         SCIP_Bool redundant;
         SCIP_Bool aggregated;
         int shift;

         shift = est - (hmin - lct + MIN(hmin, ect));
         assert(shift > 0);
         lst = hmin;
         duration = hmin - lct;

         SCIPdebugMessage("replace variable <%s>[%g,%g] by [%d,%d]\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), est + shift, lst);

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_aggr", SCIPvarGetName(var));
         SCIP_CALL( SCIPcreateVar(scip, &aggrvar, name, (SCIP_Real)(est+shift), (SCIP_Real)lst, 0.0, SCIPvarGetType(var),
               SCIPvarIsInitial(var), SCIPvarIsRemovable(var), NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, var) );
         SCIP_CALL( SCIPaggregateVars(scip, var, aggrvar, 1.0, -1.0, (SCIP_Real)shift, &infeasible, &redundant, &aggregated) );

         assert(!infeasible);
         assert(!redundant);
         assert(aggregated);

         /* replace variable */
         consdata->durations[v] = duration;
         consdata->vars[v] = aggrvar;

         /* remove and add locks */
         SCIP_CALL( SCIPunlockVarCons(scip, var, cons, consdata->downlocks[v], consdata->uplocks[v]) );
         SCIP_CALL( SCIPlockVarCons(scip, var, cons, consdata->downlocks[v], consdata->uplocks[v]) );

         SCIP_CALL( SCIPreleaseVar(scip, &aggrvar) );

         (*naggrvars)++;
      }
   }

   return SCIP_OKAY;
}
#endif

/** presolve given constraint */
static
SCIP_RETCODE presolveCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   int*                  nchgbds,            /**< pointer to store the number of changed bounds */
   int*                  ndelconss,          /**< pointer to store the number of deleted constraints */
   int*                  naddconss,          /**< pointer to store the number of added constraints */
   int*                  nchgcoefs,          /**< pointer to store the number of changed coefficients */
   int*                  nchgsides,          /**< pointer to store the number of changed sides */
   SCIP_Bool*            cutoff,             /**< pointer to store if a cutoff was detected */
   SCIP_Bool*            unbounded           /**< pointer to store if the problem is unbounded */
   )
{
#ifdef SCIP_STATISTIC
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
#endif

   /* over sized jobs should be removed */
   assert(checkDemands(scip, cons));
   assert(!SCIPconsIsDeleted(cons));

   if( conshdlrdata->dualpresolve )
   {
      SCIPstatistic( conshdlrdata->ndecomps -= consdata->ndecomps );

      /* in case the cumulative constraint is independent of every else, solve the cumulative problem and apply the
       * fixings (dual reductions)
       */
      SCIP_CALL( solveIndependentCons(scip, cons, conshdlrdata->maxnodes, nchgbds, nfixedvars, ndelconss, cutoff, unbounded) );

      if( *cutoff || *unbounded )
         return SCIP_OKAY;

      /* computes the effective horizon and checks if the constraint can be decompesd */
      SCIP_CALL( computeEffectiveHorizon(scip, cons, ndelconss, naddconss, nchgsides) );

      SCIPstatistic( conshdlrdata->ndecomps += consdata->ndecomps );

      if( SCIPconsIsDeleted(cons) )
         return SCIP_OKAY;

      SCIPstatistic( conshdlrdata->nirrelevantjobs -= consdata->nirrelevantjobs );
      SCIPstatistic( conshdlrdata->nremovedlocks -= consdata->nremovedlocks );
      SCIPstatistic( conshdlrdata->ndualfixs -= consdata->ndualfixs );
      SCIPstatistic( conshdlrdata->nalwaysruns -= consdata->nalwaysruns );

      /* presolve constraint form the earlier start time point of view */
      SCIP_CALL( presolveConsEst(scip, cons, nfixedvars, nchgcoefs, nchgsides, cutoff) );

      /* presolve constraint form the latest completion time point of view */
      SCIP_CALL( presolveConsLct(scip, cons, nfixedvars, nchgcoefs, nchgsides, cutoff) );

      SCIPstatistic( conshdlrdata->nirrelevantjobs += consdata->nirrelevantjobs );
      SCIPstatistic( conshdlrdata->nremovedlocks += consdata->nremovedlocks );
      SCIPstatistic( conshdlrdata->ndualfixs += consdata->ndualfixs );
      SCIPstatistic( conshdlrdata->nalwaysruns += consdata->nalwaysruns );

      /* remove jobs which have a duration or demand of zero ????????????? */
      SCIP_CALL( removeOversizedJobs(scip, cons, nchgbds, nchgcoefs, naddconss, cutoff) );

      if( *cutoff || SCIPconsIsDeleted(cons) )
         return SCIP_OKAY;
   }

   if( conshdlrdata->normalize )
   {
      /* divide demands by their greatest common divisor */
      SCIP_CALL( normalizeDemands(scip, cons, nchgcoefs, nchgsides) );
   }

   /* delete constraint with one job */
   SCIP_CALL( deleteTrivilCons(scip, cons, ndelconss, cutoff) );

   if( *cutoff || SCIPconsIsDeleted(cons) )
      return SCIP_OKAY;

   if( conshdlrdata->coeftightening )
   {
      /* try to tighten the capacity */
      SCIP_CALL( tightenCapacity(scip, cons, nchgcoefs, nchgsides) );

      /* try to tighten the coefficients */
      SCIP_CALL( tightenCoefs(scip, cons, nchgcoefs) );
   }

   assert(checkDemands(scip, cons) || *cutoff);

#if 0
   SCIP_CALL( reformulateCons(scip, cons, naggrvars) );
#endif

   return SCIP_OKAY;
}

/**@} */


/**@name Callback methods of constraint handler
 *
 * @{
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyCumulative)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrCumulative(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreCumulative)
{  /*lint --e{715}*/
   int c;

   for( c = 0; c < nconss; ++c )
   {
      /* remove jobs which have a duration or demand of zero (zero energy) or lay outside the effective horizon [hmin,
       * hmax)
       */
      SCIP_CALL( removeIrrelevantJobs(scip, conss[c]) );
   }

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#ifdef SCIP_STATISTIC
static
SCIP_DECL_CONSEXITPRE(consExitpreCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( evaluateCumulativeness(scip, conss[c]) );

      SCIP_CALL( SCIPconsdataVisualize(scip, conss[c]) );
   }

   SCIPstatisticMessage("@33  irrelevant %d\n", conshdlrdata->nirrelevantjobs);
   SCIPstatisticMessage("@44  dual %d\n", conshdlrdata->ndualfixs);
   SCIPstatisticMessage("@55  locks %d\n", conshdlrdata->nremovedlocks);
   SCIPstatisticMessage("@66  decomp %d\n", conshdlrdata->ndecomps);
   SCIPstatisticMessage("@77  allconsdual %d\n", conshdlrdata->nallconsdualfixs);
   SCIPstatisticMessage("@88  alwaysruns %d\n", conshdlrdata->nalwaysruns);

   return SCIP_OKAY;
}
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* free rows */
      SCIP_CALL( consdataFreeRows(scip, &consdata) );
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteCumulative)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL );
   assert(*consdata != NULL );

   /* if constraint belongs to transformed problem space, drop bound change events on variables */
   if( (*consdata)->nvars > 0 && SCIPvarIsTransformed((*consdata)->vars[0]) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( consdataDropAllEvents(scip, *consdata, conshdlrdata->eventhdlr) );
   }

   /* free cumulative constraint data */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->demandrows == NULL);

   SCIPdebugMessage("transform cumulative constraint <%s>\n", SCIPconsGetName(sourcecons));

   /* get event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->vars, sourcedata->linkingconss,
         sourcedata->durations, sourcedata->demands, sourcedata->nvars, sourcedata->capacity,
         sourcedata->hmin, sourcedata->hmax) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* catch bound change events of variables */
   SCIP_CALL( consdataCatchEvents(scip, targetdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpCumulative)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("initialize LP relaxation for %d cumulative constraints\n", nconss);

   if( conshdlrdata->usebinvars )
   {
      /* add rows to LP */
      for( c = 0; c < nconss; ++c )
      {
         assert(SCIPconsIsInitial(conss[c]));
         SCIP_CALL( addRelaxation(scip, conss[c], conshdlrdata->cutsasconss) );

         if( conshdlrdata->cutsasconss )
         {
            SCIP_CALL( SCIPrestartSolve(scip) );
         }
      }
   }

   /**@todo if we want to use only the integer variables; only these will be in cuts
    *       create some initial cuts, currently these are only separated */

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpCumulative)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool reducedom;
   SCIP_Bool separated;
   int c;

   SCIPdebugMessage("consSepalpCumulative\n");

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("separating %d/%d cumulative constraints\n", nusefulconss, nconss);

   cutoff = FALSE;
   reducedom = FALSE;
   separated = FALSE;
   (*result) = SCIP_DIDNOTRUN;

   if( !conshdlrdata->localcuts && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   (*result) = SCIP_DIDNOTFIND;

   if( conshdlrdata->usebinvars )
   {
      /* check all useful cumulative constraints for feasibility  */
      for( c = 0; c < nusefulconss && !reducedom && !cutoff; ++c )
      {
         SCIP_CALL( separateConsBinaryRepresentation(scip, conss[c], NULL, &separated) );
      }

      if( !cutoff && !reducedom && conshdlrdata->usecovercuts )
      {
         for( c = 0; c < nusefulconss; ++c )
         {
            SCIP_CALL( separateCoverCutsCons(scip, conss[c], NULL, &separated) );
         }
      }
   }

   if( conshdlrdata->sepaold )
   {
      /* separate cuts containing only integer variables */
      for( c = 0; c < nusefulconss; ++c )
      {
         SCIP_CALL( separateConsOnIntegerVariables(scip, conss[c], NULL, TRUE, &separated) );
         SCIP_CALL( separateConsOnIntegerVariables(scip, conss[c], NULL, FALSE, &separated) );
      }
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reducedom )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool reducedom;
   SCIP_Bool separated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( !conshdlrdata->localcuts && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("separating %d/%d cumulative constraints\n", nusefulconss, nconss);

   cutoff = FALSE;
   reducedom = FALSE;
   separated = FALSE;
   (*result) = SCIP_DIDNOTFIND;

   if( conshdlrdata->usebinvars )
   {
      /* check all useful cumulative constraints for feasibility  */
      for( c = 0; c < nusefulconss && !cutoff && !reducedom; ++c )
      {
         SCIP_CALL( separateConsBinaryRepresentation(scip, conss[c], NULL, &separated) );
      }

      if( !cutoff && !reducedom && conshdlrdata->usecovercuts )
      {
         for( c = 0; c < nusefulconss; ++c )
         {
            SCIP_CALL( separateCoverCutsCons(scip, conss[c], sol, &separated) );
         }
      }
   }
   if( conshdlrdata->sepaold )
   {
      /* separate cuts containing only integer variables */
      for( c = 0; c < nusefulconss; ++c )
      {
         SCIP_CALL( separateConsOnIntegerVariables(scip, conss[c], NULL, TRUE, &separated) );
         SCIP_CALL( separateConsOnIntegerVariables(scip, conss[c], NULL, FALSE, &separated) );
      }
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reducedom )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   if( solinfeasible )
   {
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("LP enforcing %d useful resource constraints of %d constraints\n", nusefulconss, nconss);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   (*result) = SCIP_FEASIBLE;

   if( conshdlrdata->usebinvars )
   {
      SCIP_Bool separated;
      int c;

      separated = FALSE;

      /* first check if a constraints is violated */
      for( c = 0; c < nusefulconss; ++c )
      {
         SCIP_CONS* cons;
         SCIP_Bool violated;

         cons = conss[c];
         assert(cons != NULL);

         SCIP_CALL( checkCons(scip, cons, NULL, &violated, FALSE) );

         if( !violated )
            continue;

         SCIP_CALL( separateConsBinaryRepresentation(scip, cons, NULL, &separated) );
      }

      for( ; c < nconss && !separated; ++c )
      {
         SCIP_CONS* cons;
         SCIP_Bool violated;

         cons = conss[c];
         assert(cons != NULL);

         SCIP_CALL( checkCons(scip, cons, NULL, &violated, FALSE) );

         if( !violated )
            continue;

         SCIP_CALL( separateConsBinaryRepresentation(scip, cons, NULL, &separated) );
      }

      if( separated )
         (*result) = SCIP_SEPARATED;
   }
   else
   {
      SCIP_CALL( enforceSolution(scip, conss, nconss, conshdlrdata->fillbranchcands, result) );
   }

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIPdebugMessage("method: enforce pseudo solution\n");

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   if( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   (*result) = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( enforceSolution(scip, conss, nconss, conshdlrdata->fillbranchcands, result) );

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckCumulative)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   violated = FALSE;

   SCIPdebugMessage("check %d cumulative constraints\n", nconss);

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
SCIP_DECL_CONSPROP(consPropCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   int nchgbds;
   int ndelconss;
   int c;

   SCIPdebugMessage("propagate %d of %d useful cumulative constraints\n", nusefulconss, nconss);

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   nchgbds = 0;
   ndelconss = 0;
   cutoff = FALSE;
   (*result) = SCIP_DIDNOTRUN;

   /* propgate all useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CONS* cons;

      cons = conss[c];
      assert(cons != NULL);

      if( SCIPgetDepth(scip) == 0 )
      {
         SCIP_CALL( presolveCons(scip, cons, conshdlrdata, &nchgbds, &nchgbds, &ndelconss, &nchgbds, &nchgbds, &nchgbds, &cutoff, &cutoff) );
      }

      if( cutoff )
         break;

      if( SCIPconsIsDeleted(cons) )
         continue;

      SCIP_CALL( propagateCons(scip, cons, conshdlrdata, &nchgbds, &ndelconss, &cutoff) );
   }

   if( !cutoff && nchgbds == 0 )
   {
      /* propgate all other constraints */
      for( c = nusefulconss; c < nconss && !cutoff; ++c )
      {
         SCIP_CALL( propagateCons(scip, conss[c], conshdlrdata, &nchgbds, &ndelconss, &cutoff) );
      }
   }

#if 0
   if( !cutoff && conshdlrdata->dualpresolve )
   {
      SCIP_CALL( propagateAllConss(scip, conshdlrdata, conss, nconss, TRUE, &nchgbds, NULL) );
   }
#endif

   if( cutoff )
   {
      SCIPdebugMessage("detected infeasible\n");
      *result = SCIP_CUTOFF;
   }
   else if( nchgbds > 0 )
   {
      SCIPdebugMessage("delete (locally) %d constraints and changed %d variable bounds\n", ndelconss, nchgbds);
      *result = SCIP_REDUCEDDOM;
   }
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_Bool cutoff;
   SCIP_Bool unbounded;
   int oldnfixedvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnaddconss;
   int oldnupgdconss;
   int oldnchgsides;
   int oldnchgcoefs;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("presolve cumulative constraints\n");

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTRUN;

   oldnfixedvars = *nfixedvars;
   oldnchgbds = *nchgbds;
   oldnchgsides = *nchgsides;
   oldnchgcoefs = *nchgcoefs;
   oldnupgdconss = *nupgdconss;
   oldndelconss = *ndelconss;
   oldnaddconss = *naddconss;
   cutoff = FALSE;
   unbounded = FALSE;

   /* process constraints */
   for( c = 0; c < nconss && !cutoff; ++c )
   {
      cons = conss[c];

      SCIP_CALL( presolveCons(scip, cons, conshdlrdata, nfixedvars, nchgbds, ndelconss, naddconss, nchgcoefs, nchgsides, &cutoff, &unbounded) );

      if( cutoff || unbounded )
         break;

      if( SCIPconsIsDeleted(cons) )
         continue;

      /* propagate cumulative constraint */
      SCIP_CALL( propagateCons(scip, cons, conshdlrdata, nchgbds, ndelconss, &cutoff) );
      assert(checkDemands(scip, cons) || cutoff);
   }

   if( !cutoff && !unbounded && conshdlrdata->dualpresolve )
   {
      SCIP_CALL( propagateAllConss(scip, conshdlrdata, conss, nconss, FALSE, nfixedvars, NULL) );
   }

   SCIPdebugMessage("delete %d constraints and changed %d variable bounds (cutoff %u)\n",
      *ndelconss - oldndelconss, *nchgbds - oldnchgbds, cutoff);

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( unbounded )
      *result = SCIP_UNBOUNDED;
   else if( *nchgbds > oldnchgbds || *nfixedvars > oldnfixedvars || *nchgsides > oldnchgsides
      || *nchgcoefs > oldnchgcoefs  || *nupgdconss > oldnupgdconss || *ndelconss > oldndelconss || *naddconss > oldnaddconss )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(infervar != NULL);
   assert(bdchgidx != NULL);

   /* process constraint */
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("resolve propagation: variable <%s>, cumulative constraint <%s> (capacity %d, propagation %d)\n",
      SCIPvarGetName(infervar), SCIPconsGetName(cons), consdata->capacity, inferInfoGetProprule(intToInferInfo(inferinfo)));

   SCIP_CALL( respropCumulativeCondition(scip, consdata->nvars, consdata->vars,
         consdata->durations, consdata->demands, consdata->capacity,
         infervar, intToInferInfo(inferinfo), boundtype, bdchgidx, NULL, result) );

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int v;

   SCIPdebugMessage("lock cumulative constraint <%s> with nlockspos = %d, nlocksneg = %d\n", SCIPconsGetName(cons), nlockspos, nlocksneg);

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   assert(vars != NULL);

   for( v = 0; v < consdata->nvars; ++v )
   {
      if( consdata->downlocks[v] && consdata->uplocks[v] )
      {
         /* the integer start variable should not get rounded in both direction  */
         SCIP_CALL( SCIPaddVarLocks(scip, vars[v], nlockspos + nlocksneg, nlockspos + nlocksneg) );
      }
      else if( consdata->downlocks[v]  )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, vars[v], nlockspos, nlocksneg) );
      }
      else if( consdata->uplocks[v] )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, vars[v], nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintCumulative)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdataPrint(scip, SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourceconsdata;
   SCIP_VAR** sourcevars;
   SCIP_VAR** vars;
   const char* consname;

   int nvars;
   int v;

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert(sourceconsdata != NULL);

   /* get variables of the source constraint */
   nvars = sourceconsdata->nvars;
   sourcevars = sourceconsdata->vars;

   (*valid) = TRUE;

   if( nvars == 0 )
      return SCIP_OKAY;

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   for( v = 0; v < nvars && *valid; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &vars[v], varmap, consmap, global, valid) );
      assert(!(*valid) || vars[v] != NULL);
   }

   /* only create the target constraint, if all variables could be copied */
   if( *valid )
   {
      if( name != NULL )
         consname = name;
      else
         consname = SCIPconsGetName(sourcecons);

      /* copy the logic using the linear constraint copy method */
      SCIP_CALL( SCIPcreateConsCumulative(scip, cons, consname, nvars, vars,
            sourceconsdata->durations, sourceconsdata->demands, sourceconsdata->capacity,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

      /* adjust left side if the time axis if needed */
      if( sourceconsdata->hmin > 0 )
      {
         SCIP_CALL( SCIPsetHminCumulative(scip, *cons, sourceconsdata->hmin) );
      }

      /* adjust right  side if the time axis if needed */
      if( sourceconsdata->hmax < INT_MAX )
      {
         SCIP_CALL( SCIPsetHmaxCumulative(scip, *cons, sourceconsdata->hmax) );
      }
   }

   /* free buffer array */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseCumulative)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real value;
   char strvalue[SCIP_MAXSTRLEN];
   char* endptr;
   int* demands;
   int* durations;
   int capacity;
   int duration;
   int demand;
   int hmin;
   int hmax;
   int varssize;
   int nvars;

   SCIPdebugMessage("parse <%s> as cumulative constraint\n", str);

   /* cutoff "cumulative" form the constraint string */
   SCIPstrCopySection(str, 'c', '(', strvalue, SCIP_MAXSTRLEN, &endptr);
   str = endptr;

   varssize = 100;
   nvars = 0;

   /* allocate buffer array for variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &durations, varssize) );

   do
   {
      SCIP_CALL( SCIPparseVarName(scip, str, &var, &endptr) );

      if( var != NULL )
      {
         str = endptr;

         SCIPstrCopySection(str, '(', ')', strvalue, SCIP_MAXSTRLEN, &endptr);
         duration = atoi(strvalue);
         str = endptr;

         SCIPstrCopySection(str, '[', ']', strvalue, SCIP_MAXSTRLEN, &endptr);
         demand = atoi(strvalue);
         str = endptr;

         SCIPdebugMessage("parse job <%s>, duration %d, demand %d\n", SCIPvarGetName(var), duration, demand);

         vars[nvars] = var;
         demands[nvars] = demand;
         durations[nvars] = duration;
         nvars++;
      }
   }
   while( var != NULL );

   /* parse effective time window */
   SCIPstrCopySection(str, '[', ',', strvalue, SCIP_MAXSTRLEN, &endptr);
   hmin = atoi(strvalue);
   str = endptr;

   if( SCIPstrToRealValue(str, &value, &endptr) )
   {
      hmax = (int)(value);
      str = endptr;

      /* parse capacity */
      SCIPstrCopySection(str, ')', '=', strvalue, SCIP_MAXSTRLEN, &endptr);
      str = endptr;
      if( SCIPstrToRealValue(str, &value, &endptr) )
      {
         capacity = (int)value;

         /* create cumulative constraint */
         SCIP_CALL( SCIPcreateConsCumulative(scip, cons, name, nvars, vars, durations, demands, capacity,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

         SCIP_CALL( SCIPsetHminCumulative(scip, *cons, hmin) );
         SCIP_CALL( SCIPsetHmaxCumulative(scip, *cons, hmax) );

         (*success) = TRUE;
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &durations);
   SCIPfreeBufferArray(scip, &demands);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsCumulative)
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
SCIP_DECL_CONSGETNVARS(consGetNVarsCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods of event handler
 *
 * @{
 */


/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   /* mark the constraint to be not propagated */
   consdata->propagated = FALSE;

   return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 * @{
 */

/*
 * constraint specific interface methods
 */

/** creates the handler for cumulative constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrCumulative(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecCumulative, NULL) );

   /* create cumulative constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpCumulative, consEnfopsCumulative, consCheckCumulative, consLockCumulative,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyCumulative, consCopyCumulative) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteCumulative) );
#ifdef SCIP_STATISTIC
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreCumulative) );
#endif
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolCumulative) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeCumulative) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsCumulative) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsCumulative) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreCumulative) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpCumulative) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseCumulative) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolCumulative, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintCumulative) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropCumulative, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropCumulative) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpCumulative, consSepasolCumulative, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransCumulative) );

   /* add cumulative constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/coretimes", "should core-times be propagated (time tabling)?",
         &conshdlrdata->coretimes, FALSE, DEFAULT_CORETIMES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/overload", "should edge finding be used to detect an overload?",
         &conshdlrdata->overload, FALSE, DEFAULT_OVERLOAD, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usebinvars", "should the binary representation be used?",
         &conshdlrdata->usebinvars, FALSE, DEFAULT_USEBINVARS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/localcuts", "should cuts be added only locally?",
         &conshdlrdata->localcuts, FALSE, DEFAULT_LOCALCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usecovercuts", "should covering cuts be added every node?",
         &conshdlrdata->usecovercuts, FALSE, DEFAULT_USECOVERCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/cutsasconss",
         "should the cumulative constraint create cuts as knapsack constraints?",
         &conshdlrdata->cutsasconss, FALSE, DEFAULT_CUTSASCONSS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/sepaold",
         "shall old sepa algo be applied?",
         &conshdlrdata->sepaold, FALSE, DEFAULT_SEPAOLD, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/fillbranchcands", "should branching candidates be added to storage?",
         &conshdlrdata->fillbranchcands, FALSE, DEFAULT_FILLBRANCHCANDS, NULL, NULL) );

   /* presolving parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/dualpresolve", "should dual presolving be applied?",
         &conshdlrdata->dualpresolve, FALSE, DEFAULT_DUALPRESOLVE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/coeftightening", "should coefficient tightening be applied?",
         &conshdlrdata->coeftightening, FALSE, DEFAULT_COEFTIGHTENING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/normalize", "should demands and capacity be normalized?",
         &conshdlrdata->normalize, FALSE, DEFAULT_NORMALIZE, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "constraints/"CONSHDLR_NAME"/maxnodes",
         "number of branch-and-bound nodes to solve an independent cumulative constraint (-1: no limit)?",
         &conshdlrdata->maxnodes, FALSE, DEFAULT_MAXNODES, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a cumulative constraint */
SCIP_RETCODE SCIPcreateConsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
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
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);

   /* find the cumulative constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage(""CONSHDLR_NAME" constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIPdebugMessage("create cumulative constraint <%s> with %d jobs\n", name, nvars);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, vars, NULL, durations, demands, nvars, capacity, 0, INT_MAX) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata,
         initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );


   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* catch bound change events of variables */
      SCIP_CALL( consdataCatchEvents(scip, consdata, conshdlrdata->eventhdlr) );
   }

   return SCIP_OKAY;
}

/** creates and captures a cumulative constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsCumulative(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsCumulative() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity            /**< available cumulative capacity */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsCumulative(scip, cons, name, nvars, vars, durations, demands, capacity,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** set the left bound of the time axis to be considered (including hmin) */
SCIP_RETCODE SCIPsetHminCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   int                   hmin                /**< left bound of time axis to be considered */
   )
{
   SCIP_CONSDATA* consdata;
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      return SCIP_INVALIDCALL;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(hmin >= 0);
   assert(hmin <= consdata->hmax);

   consdata->hmin = hmin;

   return SCIP_OKAY;
}

/** returns the left bound of the time axis to be considered */
int SCIPgetHminCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->hmin;
}

/** set the right bound of the time axis to be considered (not including hmax) */
SCIP_RETCODE SCIPsetHmaxCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   int                   hmax                /**< right bound of time axis to be considered */
   )
{
   SCIP_CONSDATA* consdata;
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
      return SCIP_INVALIDCALL; /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(hmax >= consdata->hmin);

   consdata->hmax = hmax;

   return SCIP_OKAY;
}

/** returns the right bound of the time axis to be considered */
int SCIPgetHmaxCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->hmax;
}

/** returns the activities of the cumulative constraint */
SCIP_VAR** SCIPgetVarsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** returns the activities of the cumulative constraint */
int SCIPgetNVarsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** returns the capacity of the cumulative constraint */
int SCIPgetCapacityCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->capacity;
}

/** returns the durations of the cumulative constraint */
int* SCIPgetDurationsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->durations;
}

/** returns the demands of the cumulative constraint */
int* SCIPgetDemandsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->demands;
}

/** check for the given starting time variables with their demands and durations if the cumulative conditions for the
 *  given solution is satisfied
 */
SCIP_RETCODE SCIPcheckCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_Bool*            violated,           /**< pointer to store if the cumulative condition is violated */
   SCIP_CONS*            cons,               /**< constraint which is checked */
   SCIP_Bool             printreason         /**< should the reason for the violation be printed? */
   )
{
   assert(scip != NULL);
   assert(violated != NULL);

   SCIP_CALL( checkCumulativeCondition(scip, sol, nvars, vars, durations, demands, capacity, hmin, hmax,
         violated, cons, printreason) );

   return SCIP_OKAY;
}

/** propagate the given cumulative condition */
SCIP_RETCODE SCIPpropCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_CONS*            cons,               /**< constraint which gets propagated */
   int*                  nchgbds,            /**< pointer to store the number of variable bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            explanation,        /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   SCIP_Bool*            cutoff              /**< pointer to store if the cumulative condition is violated */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool redundant;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(initialized != NULL);
   assert(*initialized == FALSE);
   assert(cutoff != NULL);
   assert(*cutoff == FALSE);

   /* find the cumulative constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage(""CONSHDLR_NAME" constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   redundant = FALSE;

   SCIP_CALL( propagateCumulativeCondition(scip, conshdlrdata,
         nvars, vars, durations, demands, capacity,  hmin, hmax, cons,
         nchgbds, &redundant, initialized, explanation, cutoff) );

   return SCIP_OKAY;
}

/** resolve propagation w.r.t. the cumulative condition */
SCIP_RETCODE SCIPrespropCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             infervar,           /**< the conflict variable whose bound change has to be resolved */
   int                   inferinfo,          /**< the user information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Bool*            explanation,        /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   SCIP_CALL( respropCumulativeCondition(scip, nvars, vars, durations, demands, capacity,
         infervar, intToInferInfo(inferinfo), boundtype, bdchgidx, explanation, result) );

   return SCIP_OKAY;
}

/** this method visualizes the cumulative structure in GML format */
SCIP_RETCODE SCIPconsdataVisualize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< cumulative constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_HASHTABLE* vars;
   FILE* file;
   SCIP_VAR* var;
   char filename[SCIP_MAXSTRLEN];
   int nvars;
   int v;

   /* open file */
   (void)SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s.gml", SCIPconsGetName(cons));
   file = fopen(filename, "w");

   /* check if the file was open */
   if( file == NULL )
   {
      SCIPerrorMessage("cannot create file <%s> for writing\n", filename);
      SCIPprintSysError(filename);
      return SCIP_FILECREATEERROR;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   SCIP_CALL( SCIPhashtableCreate(&vars, SCIPblkmem(scip), SCIPcalcHashtableSize(nvars),
         SCIPvarGetHashkey, SCIPvarIsHashkeyEq, SCIPvarGetHashkeyVal, NULL) );

   /* create opening of the GML format */
   SCIPgmlWriteOpening(file,  TRUE);

   for( v = 0; v < nvars; ++v )
   {
      char color[SCIP_MAXSTRLEN];

      var = consdata->vars[v];
      assert(var != NULL);

      SCIP_CALL( SCIPhashtableInsert(vars, (void*)var) );

      if( SCIPvarGetUbGlobal(var) - SCIPvarGetLbGlobal(var) < 0.5 )
         (void)SCIPsnprintf(color, SCIP_MAXSTRLEN, "%s", "#0000ff");
      else if( !consdata->downlocks[v] || !consdata->uplocks[v] )
         (void)SCIPsnprintf(color, SCIP_MAXSTRLEN, "%s", "#00ff00");
      else
         (void)SCIPsnprintf(color, SCIP_MAXSTRLEN, "%s", "#ff0000");

      SCIPgmlWriteNode(file, (unsigned int)(size_t)var, SCIPvarGetName(var), "rectangle", color, NULL);
   }

   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR** vbdvars;
      int nvbdvars;
      int b;

      var = consdata->vars[v];
      assert(var != NULL);

      vbdvars = SCIPvarGetVlbVars(var);
      nvbdvars = SCIPvarGetNVlbs(var);

      for( b = 0; b < nvbdvars; ++b )
      {
         if( SCIPhashtableExists(vars, (void*)vbdvars[b]) )
         {
            SCIPgmlWriteArc(file, (unsigned int)(size_t)vbdvars[b], (unsigned int)(size_t)var, NULL, NULL);
         }
      }

#if 0
      vbdvars = SCIPvarGetVubVars(var);
      nvbdvars = SCIPvarGetNVubs(var);

      for( b = 0; b < nvbdvars; ++b )
      {
         if( SCIPhashtableExists(vars, vbdvars[b]) )
         {
            SCIPgmlWriteArc(file, (unsigned int)(size_t)var, (unsigned int)(size_t)vbdvars[b], NULL, NULL);
         }
      }
#endif
   }

   /* create closing of the GML format */
   SCIPgmlWriteCosing(file);

   /* close file */
   fclose(file);

   SCIPhashtableFree(&vars);

   return SCIP_OKAY;
}

/**@} */
