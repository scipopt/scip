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

/**@file   probdata_rpa.c
 * @brief  Problem data for ringpacking problem
 * @author Benjamin Mueller
 *
 * This file handles the main problem data used in that project. For more details see \ref RINGPACKING_PROBLEMDATA page.
 *
 * @page RINGPACKING_PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the ringpacking problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access most data.
 *
 * The function SCIPprobdataCreate(), which is called in the \ref reader_bpa.c "reader plugin" after the input file was
 * parsed, initializes the problem data structure. Afterwards, the problem is setup in SCIPprobdataSetupProblem. For this,
 * it enumerates all dominating circular patterns, selects a set of initial rectangular patterns and creates the
 * corresponding variables and constraints. Note that the pattern constraints have to have the
 * <code>modifiable</code>-flag set to TRUE. This is necessary to tell the solver that these constraints are not
 * completed yet. This means, during the search new variables/patterns might be added. The solver needs this information
 * because certain reductions are not allowed.
 *
 * A list of all interface methods can be found in probdata_binpacking.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "probdata_rpa.h"
#include "pricer_rpa.h"

#include <string.h>
#include <math.h>

/* properties of the ringpacking statistics table */
#define TABLE_NAME_RPA                       "ringpacking"
#define TABLE_DESC_RPA                       "ringpacking statistics"
#define TABLE_POSITION_RPA                   12500                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_RPA             SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

#ifndef M_PI
#define M_PI           3.141592653589793238462643
#endif

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the ringpacking, all variables which are created, and all
 * constraints.
 */
struct SCIP_ProbData
{
   int*                  demands;            /**< array of demands */
   SCIP_Real*            rints;              /**< internal radii of each ring */
   SCIP_Real*            rexts;              /**< external radii of each ring */
   int                   ntypes;             /**< number of different types */

   SCIP_Real             width;              /**< height of each rectangle */
   SCIP_Real             height;             /**< width of each rectangle */

   SCIP_CONS**           patternconss;       /**< pattern constraints for each type */

   /* circular pattern data */
   SCIP_PATTERN**        cpatterns;          /**< array containing all circular patterns */
   SCIP_VAR**            cvars;              /**< variables corresponding to circular patterns */
   int                   ncpatterns;         /**< total number of circular patterns */
   int                   cpatternsize;       /**< size of cpatterns and cvars array */

   /* rectangular pattern data */
   SCIP_PATTERN**        rpatterns;          /**< array containing all rectangular patterns */
   SCIP_VAR**            rvars;              /**< variables corresponding to rectangular patterns */
   int                   nrpatterns;         /**< total number of rectangular patterns */
   int                   rpatternsize;       /**< size of rpatterns and rvars array */

   /* variables for statistics */
   int                   ncppatternsunknownbeg;/**< number of unknown circular patterns after enumeration step */
   SCIP_Real             enumtime;           /**< time spend for enumerating circular patterns */
   SCIP_Bool             isdualinvalid;      /**< whether the following reported dual bounds are valid */
   SCIP_Real             dualbound;          /**< valid dual bound for RCPP instance */

   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */


   /* parameters */
   SCIP_Real             nlptilimsoft;       /**< soft time limit for verification NLP */
   SCIP_Real             heurtilimsoft;      /**< soft time limit for verification heuristic */
   SCIP_Real             totaltilimsoft;     /**< soft time limit for enumerating circular patterns */
   SCIP_Longint          nlpnodelimsoft;     /**< soft node limit for verification NLP */
   int                   heuriterlimsoft;    /**< soft iteration limit for verification heuristic */
};


/**@name Local methods
 *
 * @{
 */

/** auxiliary function to create problem data;
 *
 * @note captures patterns and corresponding variables
 */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   SCIP_CONS**           patternconss,       /**< pattern constraints */
   SCIP_PATTERN**        cpatterns,          /**< circular patterns */
   SCIP_VAR**            cvars,              /**< variables corresponding to circular patterns */
   int                   ncpatterns,         /**< total number of circular patterns */
   SCIP_PATTERN**        rpatterns,          /**< rectangular patterns */
   SCIP_VAR**            rvars,              /**< variables corresponding to rectangular patterns */
   int                   nrpatterns,         /**< total number of rectangular patterns */
   int*                  demands,            /**< array containing the demands */
   SCIP_Real*            rints,              /**< interal radii of each ring */
   SCIP_Real*            rexts,              /**< external radii of each ring */
   int                   ntypes,             /**< number of different types */
   SCIP_Real             width,              /**< width of each rectangle */
   SCIP_Real             height              /**< height of each rectangle */
   )
{
   int i;

   assert(probdata != NULL);
   assert(demands != NULL);
   assert(rints != NULL);
   assert(rexts != NULL);
   assert(ntypes > 0);
   assert(height > 0.0);
   assert(width >= height);

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, probdata) );
   BMSclearMemory(*probdata);

   if( ncpatterns > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->cvars, cvars, ncpatterns) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->cpatterns, cpatterns, ncpatterns) );
      (*probdata)->ncpatterns = ncpatterns;
      (*probdata)->cpatternsize = ncpatterns;

      /* capture circular patterns */
      for( i = 0; i < ncpatterns; ++i )
         SCIPpatternCapture(cpatterns[i]);
   }

   if( nrpatterns > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->rvars, rvars, nrpatterns) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->rpatterns, rpatterns, nrpatterns) );
      (*probdata)->nrpatterns = nrpatterns;
      (*probdata)->rpatternsize = nrpatterns;

      /* capture rectangular patterns */
      for( i = 0; i < nrpatterns; ++i )
         SCIPpatternCapture(rpatterns[i]);
   }

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->demands, demands, ntypes) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->rints, rints, ntypes) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->rexts, rexts, ntypes) );

   /* copy pattern constraints if available, otherwise allocate enough memory */
   if( patternconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->patternconss, patternconss, ntypes) );
   }
   else
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->patternconss, ntypes) );
      BMSclearMemoryArray((*probdata)->patternconss, ntypes);
   }

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &(*probdata)->randnumgen, 0, TRUE) );

   (*probdata)->ntypes = ntypes;
   (*probdata)->width = width;
   (*probdata)->height = height;

   return SCIP_OKAY;
}

/** auxiliary function to free problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to release the probdata */
   )
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);
   assert(*probdata != NULL);

   /* free random number generator */
   if( (*probdata)->randnumgen != NULL )
   {
      SCIPfreeRandom(scip, &(*probdata)->randnumgen);
   }

   /* release pattern constraints */
   if( (*probdata)->patternconss != NULL )
   {
      for( i = 0; i < SCIPprobdataGetNTypes(*probdata); ++i )
      {
         if( (*probdata)->patternconss[i] != NULL )
         {
            SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->patternconss[i]));
         }
      }
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->patternconss, (*probdata)->ntypes);
   }

   /* release circular patterns */
   for( i = 0; i < (*probdata)->ncpatterns; ++i )
   {
      assert((*probdata)->cpatterns[i] != NULL);

      if( (*probdata)->cvars[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->cvars[i]) );
      }

      SCIPpatternRelease(scip, &(*probdata)->cpatterns[i]);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->cpatterns, (*probdata)->cpatternsize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->cvars, (*probdata)->cpatternsize);

   /* release rectangular patterns */
   for( i = 0; i < (*probdata)->nrpatterns; ++i )
   {
      assert((*probdata)->rpatterns[i] != NULL);
      assert((*probdata)->rvars[i] != NULL);

      if( (*probdata)->rvars[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->rvars[i]) );
      }

      SCIPpatternRelease(scip, &(*probdata)->rpatterns[i]);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->rpatterns, (*probdata)->rpatternsize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->rvars, (*probdata)->rpatternsize);

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->rexts, (*probdata)->ntypes);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->rints, (*probdata)->ntypes);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->demands, (*probdata)->ntypes);

   SCIPfreeBlockMemory(scip, probdata);
   SCIP_CALL( SCIPsetProbData(scip, NULL) );

   return SCIP_OKAY;
}

/** counts the number of circular patterns with a given packable status */
static
int getNCPatterns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PACKABLE         status              /**< packable status */
   )
{
   int count = 0;
   int p;

   assert(probdata != NULL);

   for( p = 0; p < probdata->ncpatterns; ++p )
   {
      if( SCIPpatternGetPackableStatus(probdata->cpatterns[p]) == status )
         ++count;
   }

   return count;
}

/** ensures a minimum size of the pattern and variable arrays */
static
SCIP_RETCODE ensureSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERNTYPE      type,               /**< pattern type */
   int                   size                /**< required size */
   )
{
   int newsize;

   assert(probdata != NULL);
   assert(size > 0);

   if( type == SCIP_PATTERNTYPE_CIRCULAR && size > probdata->cpatternsize )
   {
      newsize = MAX(size, 2 * probdata->cpatternsize);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->cpatterns, probdata->cpatternsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->cvars, probdata->cpatternsize, newsize) );
      probdata->cpatternsize = newsize;
   }
   else if( type == SCIP_PATTERNTYPE_RECTANGULAR && size > probdata->rpatternsize )
   {
      newsize = MAX(size, 2 * probdata->rpatternsize);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->rpatterns, probdata->rpatternsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->rvars, probdata->rpatternsize, newsize) );
      probdata->rpatternsize = newsize;
   }

   return SCIP_OKAY;
}

/** create variables for all existing circular and rectangular patterns */
static
SCIP_RETCODE createPatternVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   SCIP_VAR* var;
   SCIP_PATTERN* pattern;
   char name[SCIP_MAXSTRLEN];
   int k;

   assert(probdata != NULL);
   assert(probdata->ncpatterns > 0);
   assert(probdata->nrpatterns > 0);

   /* create variables for circular patterns */
   for( k = 0; k < probdata->ncpatterns; ++k )
   {
      SCIP_Real ub;
      int type;
      int i;

      pattern = probdata->cpatterns[k];
      assert(pattern != NULL);

      type = SCIPpatternGetCircleType(pattern);
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", type);
      ub = (SCIP_Real)SCIPprobdataGetDemands(probdata)[type];

      /* create variable name */
      for( i = 0; i < SCIPpatternGetNElemens(pattern); ++i )
      {
         char strtmp[SCIP_MAXSTRLEN];
         int elemtype = SCIPpatternGetElementType(pattern, i);
         (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", elemtype);
         (void) strcat(name, strtmp);
      }

      /* create variable */
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, ub, 0.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, var) );

      /* store variables in problem data */
      probdata->cvars[k] = var;
   }

   /* create variables for rectangular patterns */
   for( k = 0; k < probdata->nrpatterns; ++k )
   {
      int i;

      pattern = probdata->rpatterns[k];
      assert(pattern != NULL);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "r");

      /* create variable name */
      for( i = 0; i < SCIPpatternGetNElemens(pattern); ++i )
      {
         char strtmp[SCIP_MAXSTRLEN];
         int elemtype = SCIPpatternGetElementType(pattern, i);
         (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", elemtype);
         (void) strcat(name, strtmp);
      }

      /* create variable */
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, var) );

      /* store variables in problem data */
      probdata->rvars[k] = var;
   }

   return SCIP_OKAY;
}

/** upper bound on the number of circles of a single type that fit into a circular pattern of a given type */
static
int maxCircles(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   type,               /**< type of the circular pattern */
   int                   elemtype            /**< type of element to pack */
   )
{
   SCIP_Real _rint;
   SCIP_Real rext;
   SCIP_Real rintscaled;
   int demand;
   int n;

   assert(type >= 0 && type < SCIPprobdataGetNTypes(probdata));
   assert(elemtype >= 0 && elemtype < SCIPprobdataGetNTypes(probdata));

   _rint = SCIPprobdataGetRints(probdata)[type];
   rext = SCIPprobdataGetRexts(probdata)[elemtype];
   demand = SCIPprobdataGetDemands(probdata)[elemtype];

   /* volume-bsaed bound */
   n = MIN(demand, (int) SCIPceil(scip, SQR(_rint) / SQR(rext)));

   if( n <= 1 )
      return 1;

   /* use proven bounds on the density */
   rintscaled = _rint / rext;
   assert(rintscaled >= 1.0);

   if( SCIPisLT(scip, rintscaled, 2.0) )
      return MIN(1, n);
   else if( SCIPisLT(scip, rintscaled, 2.1547005383792515) )
      return MIN(2, n);
   else if( SCIPisLT(scip, rintscaled, 2.414213562373095) )
      return MIN(3, n);
   else if( SCIPisLT(scip, rintscaled, 2.7013016167040798) )
      return MIN(4, n);
   else if( SCIPisLT(scip, rintscaled, 3.0) )
      return MIN(5, n);
   else if( SCIPisLT(scip, rintscaled, 3.3047648709624866) )
      return MIN(7, n); /* note that here is a jump and 7 is correct */
   else if( SCIPisLT(scip, rintscaled, 3.613125929752753) )
      return MIN(8, n);

   return n;
}

/** helper function to compare two patterns; returns
 *
 *   -1 if p dominates q
 *   +1 if q dominates p
 *    0 otherwise
 */
static
int isPatternDominating(
   SCIP_PATTERN*         p,                  /**< pattern */
   SCIP_PATTERN*         q,                  /**< pattern */
   int*                  count,              /**< array for counting elements of patterns */
   int                   ntypes              /**< total number of types */
   )
{
   SCIP_Bool pdomq;
   SCIP_Bool qdomp;
   int i;

   /* patterns can only dominate each other if they have the same type */
   if(SCIPpatternGetCircleType(p) != SCIPpatternGetCircleType(q) )
      return 0;

   /* reset count array */
   BMSclearMemoryArray(count, ntypes);

   /* increase array entry for each element in p */
   for( i = 0; i < SCIPpatternGetNElemens(p); ++i )
   {
      int t = SCIPpatternGetElementType(p, i);
      count[t] += 1;
   }

   /* decrease array entry for each element in q */
   for( i = 0; i < SCIPpatternGetNElemens(q); ++i )
   {
      int t = SCIPpatternGetElementType(q, i);
      count[t] -= 1;
   }

   pdomq = TRUE;
   qdomp = TRUE;

   for( i = 0; i < ntypes && (pdomq || qdomp); ++i )
   {
      if( count[i] < 0 )
         pdomq = FALSE;
      else if( count[i] > 0 )
         qdomp = FALSE;
   }

   if( pdomq && (SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_YES || SCIPpatternGetPackableStatus(q) == SCIP_PACKABLE_UNKNOWN) )
      return -1;
   else if( qdomp && (SCIPpatternGetPackableStatus(q) == SCIP_PACKABLE_YES || SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_UNKNOWN) )
      return 1;
   return 0;
}

/** filter dominated patterns */
static
SCIP_RETCODE filterPatterns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   SCIP_PATTERN** cpatterns;
   SCIP_Bool* deleted;
   int* count;
   int ncpatterns;
   int i;

   SCIP_CALL( SCIPallocBufferArray(scip, &count, SCIPprobdataGetNTypes(probdata)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cpatterns, probdata->ncpatterns) );
   SCIP_CALL( SCIPallocBufferArray(scip, &deleted, probdata->ncpatterns) );
   BMSclearMemoryArray(deleted, probdata->ncpatterns);

   for( i = 0; i < probdata->ncpatterns - 1; ++i )
   {
      SCIP_PATTERN* p = probdata->cpatterns[i];
      int j;

      if( deleted[i] )
         continue;

      for( j = i + 1; j < probdata->ncpatterns; ++j )
      {
         SCIP_PATTERN* q = probdata->cpatterns[j];
         int res;

         if( deleted[j] )
            continue;

         res = isPatternDominating(p, q, count, SCIPprobdataGetNTypes(probdata));

         /* p dominates q */
         if( res == -1 )
            deleted[j] = TRUE;
         else if( res == 1 ) /* q dominates p */
            deleted[i] = TRUE;
      }
   }

   /* remove filtered patterns */
   ncpatterns = 0;
   for( i = 0; i < probdata->ncpatterns; ++i )
   {
      if( deleted[i] )
      {
         SCIPpatternRelease(scip, &probdata->cpatterns[i]);
      }
      else
      {
         cpatterns[ncpatterns] = probdata->cpatterns[i];
         ++ncpatterns;
      }
   }
   assert(ncpatterns > 0);

   BMScopyMemoryArray(probdata->cpatterns, cpatterns, ncpatterns);
   probdata->ncpatterns = ncpatterns;

   /* free memory */
   SCIPfreeBufferArray(scip, &deleted);
   SCIPfreeBufferArray(scip, &cpatterns);
   SCIPfreeBufferArray(scip, &count);

   return SCIP_OKAY;
}

/** enumerates all circular patterns for a given type */
static
SCIP_RETCODE enumeratePatterns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern,            /**< pattern (passed for performance reasons) */
   int*                  ms,                 /**< maximum number of elements for each type (passed for performance reasons) */
   int*                  nselected,          /**< number of selected elements for each type (passed for performance reasons) */
   SCIP_Real             nlptilim,           /**< time limit for each NLP verification */
   SCIP_Real             heurtilim,          /**< time limit for each call of the heuristics */
   SCIP_Longint          nlpnodelim,         /**< node limit for each NLP verification */
   int                   heuriterlim,        /**< iteration limit for each call of the heuristics */
   SCIP_Real*            timeleft            /**< pointer to update the remaining time for the enumeration */
   )
{
   SCIP_Real* rexts;
   SCIP_Real* _rints;
   SCIP_Real maxvolume;
   SCIP_Real volume;
   int ntypes;
   int type;
   int lasttype;

   assert(ms != NULL);
   assert(pattern != NULL);
   assert(timeleft != NULL);

   type = SCIPpatternGetCircleType(pattern);
   assert(type >= 0 && type < SCIPprobdataGetNTypes(probdata));

   /* get problem data */
   rexts = SCIPprobdataGetRexts(probdata);
   _rints = SCIPprobdataGetRints(probdata);
   ntypes = SCIPprobdataGetNTypes(probdata);
   lasttype = ntypes -1;
   volume = 0.0;
   maxvolume = SQR(_rints[SCIPpatternGetCircleType(pattern)]) * M_PI; /*lint !e666*/

   /* main loop */
   while( TRUE )
   {
      SCIP_Real timelim;
      int t = lasttype;

      /* reset packable status */
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_UNKNOWN);

      SCIPdebugMsg(scip, "volume = %g <= %g\n", volume, maxvolume);

      {
         int j;
         SCIPdebugMsg(scip, "verify c%d", type);

         for( j = 0; j < SCIPpatternGetNElemens(pattern); ++j )
            SCIPdebugMsgPrint(scip, "_%d", SCIPpatternGetElementType(pattern, j));
         SCIPdebugMsgPrint(scip, "\n");
      }

      /* check volume */
      if( SCIPisLE(scip, volume, maxvolume) )
      {
         /*
          * try to verify with heuristic
          */

         /* compute time limit */
         timelim = MIN(heurtilim, *timeleft);

         /* verify pattern */
         *timeleft += SCIPgetTotalTime(scip);
         SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, timelim, heuriterlim) );
         *timeleft -= SCIPgetTotalTime(scip);

         /*
          * try to verify with NLP
          */
         if( SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN )
         {
            /* compute time limit */
            timelim = MIN(*timeleft, nlptilim);

            /* verify pattern */
            *timeleft += SCIPgetTotalTime(scip);
            SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, timelim, nlpnodelim) );
            *timeleft -= SCIPgetTotalTime(scip);
         }

         /* pattern is not packable -> don't add more elements */
         if( SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_NO )
         {
            SCIPpatternRemoveLastElements(pattern, nselected[t]);
            volume -= SQR(rexts[t]) * M_PI * nselected[t];
            nselected[t] = 0;
            --t;
         }
         /* otherwise add the pattern (and hope for filtering) */
         else
         {
            SCIP_CALL( SCIPprobdataAddVar(scip, probdata, pattern, NULL) );
         }
      }

      /* update selection */
      while( t > type && (nselected[t] == ms[t] || SCIPisGT(scip, volume, maxvolume)) )
      {
         SCIPpatternRemoveLastElements(pattern, nselected[t]);
         volume -= SQR(rexts[t]) * M_PI * nselected[t];
         nselected[t] = 0;
         t--;
      }

      /* check termination criterion */
      if( t == type )
         break;

      /* add element of type i to the pattern */
      assert(nselected[t] < ms[t]);
      ++(nselected[t]);
      volume += SQR(rexts[t]) * M_PI;
      SCIP_CALL( SCIPpatternAddElement(pattern, t, SCIP_INVALID, SCIP_INVALID) );
   }

   assert(SCIPpatternGetNElemens(pattern) == 0);
   assert(SCIPisZero(scip, volume));

   return SCIP_OKAY;
}

/** auxiliary function to setup the master problem */
static
SCIP_RETCODE setupProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real* rexts;
   SCIP_Real* rints;
   int* demands;
   SCIP_Real dualbound;
   SCIP_Real minrext;
   SCIP_Real volume;
   int ntypes;
   int p;
   int t;

   assert(probdata != NULL);
   assert(SCIPprobdataGetNTypes(probdata) > 0);

   /* set objective sense; tell SCIP that the objective will be always integral */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   SCIP_CALL( SCIPsetObjIntegral(scip) );

   /* get problem data */
   ntypes = SCIPprobdataGetNTypes(probdata);
   rexts = SCIPprobdataGetRexts(probdata);
   rints = SCIPprobdataGetRints(probdata);
   demands = SCIPprobdataGetDemands(probdata);

   /* compute all non-dominated circular patterns */
   probdata->enumtime -= SCIPgetTotalTime(scip);
   SCIP_CALL( SCIPprobdataEnumeratePatterns(scip, probdata, probdata->nlptilimsoft, probdata->heurtilimsoft,
      probdata->totaltilimsoft, probdata->nlpnodelimsoft, probdata->heuriterlimsoft) );
   probdata->enumtime += SCIPgetTotalTime(scip);
   probdata->ncppatternsunknownbeg = getNCPatterns(scip, probdata, SCIP_PACKABLE_UNKNOWN);

   SCIPinfoMessage(scip, NULL, "+++++++++++++ starting with |CP|=%d\n", probdata->ncpatterns,
      probdata->nrpatterns);

   /* create initial rectangular patterns */
   for( t = 0; t < ntypes; ++t )
   {
      SCIP_PATTERN* pattern;

      /* create a pattern containing a single circle of type t; set position of the circle to the left-bottom */
      SCIP_CALL( SCIPpatternCreateRectangular(scip, &pattern) );
      SCIP_CALL( SCIPpatternAddElement(pattern, t, rexts[t], rexts[t]) );
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);

      /* add and release pattern */
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, pattern, NULL) );
      SCIPpatternRelease(scip, &pattern);
   }

   /* create variables for all existing patterns */
   SCIP_CALL( createPatternVars(scip, probdata) );

   /* create demand constraints */
   for( t = 0; t < ntypes; ++t )
   {
      SCIP_CONS* cons;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "demand_%d", t);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, (SCIP_Real)demands[t], SCIPinfinity(scip) ) );

      for( p = 0; p < probdata->ncpatterns; ++p )
      {
         SCIP_PATTERN* pattern;
         SCIP_VAR* var;

         pattern = probdata->cpatterns[p];
         assert(pattern != NULL);
         assert(SCIPpatternGetPatternType(pattern) == SCIP_PATTERNTYPE_CIRCULAR);

         var = probdata->cvars[p];
         assert(var != NULL);

         /* add coefficient to the pattern if the pattern is of type t */
         if(SCIPpatternGetCircleType(pattern) == t )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, 1.0) );
         }
      }

      /* add and release constraint */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* create pattern constraints */
   for( t = 0; t < ntypes; ++t )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "patterncons_%d", t);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->patternconss[t], name, 0, NULL, NULL, 0.0,
            SCIPinfinity(scip) ) );

      /* declare constraint modifiable for adding variables during pricing */
      SCIP_CALL( SCIPsetConsModifiable(scip, probdata->patternconss[t], TRUE) );
      SCIP_CALL( SCIPaddCons(scip, probdata->patternconss[t]) );
   }

   /* add coefficients for circular patterns */
   for( p = 0; p < probdata->ncpatterns; ++p )
   {
      SCIP_PATTERN* pattern = probdata->cpatterns[p];
      SCIP_VAR* var = probdata->cvars[p];
      int type;

      assert(pattern != NULL);
      assert(SCIPpatternGetPatternType(pattern) == SCIP_PATTERNTYPE_CIRCULAR);
      assert(var != NULL);

      type = SCIPpatternGetCircleType(pattern);
      assert(type >= 0 && type < ntypes);

      /* - z_C */
      SCIP_CALL( SCIPaddCoefLinear(scip, probdata->patternconss[type], var, -1.0) );

      for( t = 0; t < ntypes; ++t )
      {
         int nelems = SCIPpatternCountElements(pattern, t);

         if( nelems > 0 )
         {
            /* + P_t z_C */
            SCIP_CALL( SCIPaddCoefLinear(scip, probdata->patternconss[t], var, (SCIP_Real)nelems) );
         }
      }
   }

   /* add coefficients for rectangular patterns */
   for( p = 0; p < probdata->nrpatterns; ++p )
   {
      SCIP_PATTERN* pattern = probdata->rpatterns[p];
      SCIP_VAR* var = probdata->rvars[p];

      assert(pattern != NULL);
      assert(SCIPpatternGetPatternType(pattern) == SCIP_PATTERNTYPE_RECTANGULAR);
      assert(var != NULL);

      for( t = 0; t < ntypes; ++t )
      {
         int nelems = SCIPpatternCountElements(pattern, t);

         if( nelems > 0 )
         {
            /* + P_t z_P */
            SCIP_CALL( SCIPaddCoefLinear(scip, probdata->patternconss[t], var, (SCIP_Real)nelems) );
         }
      }
   }

   /* compute an initial dual bound by considering the volume of all rings */
   minrext = rexts[ntypes-1];
   volume = 0.0;
   for( t = 0; t < ntypes; ++t )
   {
      SCIP_Real vol;

      /* consider ring as circle if there is no ring with a smaller radius than than inner one */
      if( SCIPisFeasLT(scip, rints[t], minrext) )
         vol = M_PI * SQR(rexts[t]);
      else
         vol = M_PI * (SQR(rexts[t]) - SQR(rints[t]));

      volume += vol * demands[t];
   }

   /* update initial dual bound */
   dualbound = SCIPfeasCeil(scip, volume / (SCIPprobdataGetWidth(probdata) * SCIPprobdataGetHeight(probdata)));
   SCIP_CALL( SCIPupdateLocalDualbound(scip, dualbound) );
   SCIPinfoMessage(scip, NULL, "+++++++++++++ volume-based bound = ceil(%g / %g) = %g\n", volume,
      SCIPprobdataGetWidth(probdata) * SCIPprobdataGetHeight(probdata), dualbound);
   SCIPprobdataUpdateDualbound(scip, probdata, dualbound);

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputRpa)
{ /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_Real* rexts;
   SCIP_Real* rints;
   SCIP_Real dualbound;
   SCIP_Real maxrint;
   SCIP_Real minrext;
   int* demands;
   int ntypes;
   int nrings;
   int t;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   ntypes = SCIPprobdataGetNTypes(probdata);
   demands = SCIPprobdataGetDemands(probdata);
   rexts = SCIPprobdataGetRexts(probdata);
   rints = SCIPprobdataGetRints(probdata);
   nrings = 0;
   maxrint = 0.0;
   minrext = SCIPinfinity(scip);

   /* use global dual bound if it is still valid */
   if( !probdata->isdualinvalid )
   {
      assert(SCIPisGE(scip, SCIPgetDualbound(scip), probdata->dualbound));
      dualbound = SCIPgetDualbound(scip);
   }
   else
      dualbound = probdata->dualbound;

   /* count the number of rings */
   for( t = 0; t < ntypes; ++t )
   {
      nrings += demands[t];
      maxrint = MAX(maxrint, rints[t]);
      minrext = MIN(minrext, rexts[t]);
   }

   SCIPinfoMessage(scip, file, "Ringpacking        : %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
      "dual", "ntypes", "nrings", "width", "height", "CP", "CP_unk", "CP_unk_end" ,"CP_infeas", "RP", "enumtime", "radiiratio");

   SCIPinfoMessage(scip, file, "  %-17s:", "");
   SCIPinfoMessage(scip, file, " %10.2f", dualbound);
   SCIPinfoMessage(scip, file, " %10d", ntypes);
   SCIPinfoMessage(scip, file, " %10d", nrings);
   SCIPinfoMessage(scip, file, " %10.2f", SCIPprobdataGetWidth(probdata));
   SCIPinfoMessage(scip, file, " %10.2f", SCIPprobdataGetHeight(probdata));
   SCIPinfoMessage(scip, file, " %10d", probdata->ncpatterns);
   SCIPinfoMessage(scip, file, " %10d", probdata->ncppatternsunknownbeg);
   SCIPinfoMessage(scip, file, " %10d", getNCPatterns(scip, probdata, SCIP_PACKABLE_UNKNOWN));
   SCIPinfoMessage(scip, file, " %10d", getNCPatterns(scip, probdata, SCIP_PACKABLE_NO));
   SCIPinfoMessage(scip, file, " %10d", probdata->nrpatterns);
   SCIPinfoMessage(scip, file, " %10.2f", probdata->enumtime);
   SCIPinfoMessage(scip, file, " %10.1f", maxrint / minrext);
   SCIPinfoMessage(scip, file, "\n");

   return SCIP_OKAY;
}

/** auxiliary function to update the best known candidate */
static
void updateBestCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            xs,                 /**< x-coordinates of packed elements */
   SCIP_Real*            ys,                 /**< y-coordinates of packed elements */
   SCIP_Real*            rexts,              /**< radii of packed elements */
   SCIP_Real             rext,               /**< radii of element that should be packed */
   SCIP_Real             rbounding,          /**< inner radius of bounding circle (ignored for rectangular patterns) */
   SCIP_Real             wbounding,          /**< width of bounding rectangular (ignored for circular patterns) */
   SCIP_Real             hbounding,          /**< height of bounding rectangular (ignored for circular patterns) */
   SCIP_Real             rmax,               /**< maximum radius of elements in the pattern */
   SCIP_PATTERNTYPE      patterntype,        /**< pattern type */
   SCIP_Bool*            ispacked,           /**< array indicating which elements are already packed */
   int*                  elements,           /**< the order of the elements in the pattern */
   int                   nelements,          /**< the total number of elements */
   SCIP_Real*            bestx,              /**< buffer to update best x-coordinate */
   SCIP_Real*            besty,              /**< buffer to update best y-coordinate */
   SCIP_Real             x,                  /**< x-coordinate of a candidate point */
   SCIP_Real             y,                  /**< y-coordinate of a candidate point */
   int                   ncalls              /**< total number of calls of the packing heuristic */
   )
{
   SCIP_Real threshold;
   SCIP_Bool isoverthreshold;
   int i;

   /* candidate is not valid -> skip */
   if( x == SCIP_INVALID || y == SCIP_INVALID ) /*lint !e777*/
      return;

   /* check whether there is an intersection with the boundary */
   if( patterntype == SCIP_PATTERNTYPE_CIRCULAR )
   {
      if( SCIPisGT(scip, x*x + y*y, SQR(rbounding - rext)) )
         return;
   }
   else
   {
      if( SCIPisLT(scip, x, rext) || SCIPisGT(scip, x, wbounding - rext)
         || SCIPisLT(scip, y, rext) || SCIPisGT(scip, y, hbounding - rext) )
         return;
   }

   /* check whether circle intersects other circles */
   for( i = 0; i < nelements; ++i )
   {
      SCIP_Real dist;

      /* only consider packed elements */
      if( !ispacked[i] )
         continue;

      dist = SQR(x - xs[i]) + SQR(y - ys[i]);

      /* check if the distance between mid points is smaller than the sum of the radii */
      if( SCIPisLT(scip, dist, SQR(rext + rexts[elements[i]])) )
         return;
   }

   threshold = (patterntype == SCIP_PATTERNTYPE_RECTANGULAR ? wbounding - 2.0*rmax - rext : rbounding - 2.0*rmax - rext);
   isoverthreshold = (ncalls % 2) == 1 && SCIPisGT(scip, *bestx, threshold) && SCIPisGT(scip, x, threshold);

   /* check whether the candidate is better than the best known candidate */
   if( *bestx == SCIP_INVALID || *besty == SCIP_INVALID
      || ((!isoverthreshold || SCIPisEQ(scip, y, *besty)) && SCIPisLT(scip, x, *bestx)) /*lint !e777*/
      || ((isoverthreshold || SCIPisEQ(scip, x, *bestx)) && SCIPisLT(scip, y, *besty)) ) /*lint !e777*/
   {
      *bestx = x;
      *besty = y;
   }
}

/** auxiliary function for computing a candidate position between a circle and the outer ring */
static
void computePosRingCircle(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  elements,           /**< types of elements that have been packed */
   int                   nelements,          /**< the total number of elements */
   SCIP_Real*            rexts,              /**< external radii */
   SCIP_Real*            xs,                 /**< x-coordinate of circle */
   SCIP_Real*            ys,                 /**< y-coordinate of circle */
   int                   pos,                /**< position of element in the elements array */
   SCIP_Bool*            ispacked,           /**< array indicating whether an element has been packed already */
   SCIP_Real             rmax,               /**< maximum radius of elements in the pattern */
   SCIP_Real             rbound,             /**< radius of bounding circle */
   SCIP_Real*            bestx,              /**< pointer to store the best x-coordinate */
   SCIP_Real*            besty,              /**< pointer to store the best y-coordinate */
   int                   ncalls              /**< total number of calls of the packing heuristic */
   )
{
   int i;

   /* consider already packed patterns */
   for( i = 0; i < nelements; ++i )
   {
      SCIP_Real alpha, a, b, c, h, u, v, n1, n2;

      /* only consider packed elements */
      if( !ispacked[i] )
         continue;

      c = SQRT(xs[i]*xs[i] + ys[i]*ys[i]);

      /* inner ring is too far away from boundary or both rings can not fit */
      if( !SCIPisGE(scip, c + rexts[elements[i]] + 2.0*rexts[elements[pos]], rbound)
         || SCIPisGT(scip, rexts[elements[pos]] + rexts[elements[i]], rbound) )
         continue;

      a = rexts[elements[pos]] + rexts[elements[i]];
      b = rbound - rexts[elements[pos]];

      /* if a ring is in the center than there are infinitely many solutions; take an arbitrary point */
      if( SCIPisZero(scip, c) )
      {
         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], rbound, -1.0, -1.0, rmax, SCIP_PATTERNTYPE_CIRCULAR,
            ispacked, elements, nelements, bestx, besty, -rbound + rexts[elements[pos]], 0.0, ncalls);
         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], rbound, -1.0, -1.0, rmax, SCIP_PATTERNTYPE_CIRCULAR,
            ispacked, elements, nelements, bestx, besty, +rbound - rexts[elements[pos]], 0.0, ncalls);
      }
      else
      {
         assert(c != 0.0);
         alpha = (c*c - b*b + a*a) / (2*c);

         if( a*a >= alpha*alpha )
         {
            h = SQRT(a*a - alpha*alpha);
            u = (c - alpha) * xs[i] / c;
            v = (c - alpha) * ys[i] / c;

            n1 = SCIPisZero(scip, v) ? 0.0 : h * (v / SQRT(v*v + u*u));
            n2 = SCIPisZero(scip, u) ? 0.0 : h * (-u / SQRT(v*v + u*u));

            updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], rbound, -1.0, -1.0, rmax, SCIP_PATTERNTYPE_CIRCULAR,
               ispacked, elements, nelements, bestx, besty, u + n1, v + n2, ncalls);
            updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], rbound, -1.0, -1.0, rmax, SCIP_PATTERNTYPE_CIRCULAR,
               ispacked, elements, nelements, bestx, besty, u - n1, v - n2, ncalls);
         }
      }
   }
}

/** auxiliary function for computing trivial candidate positions */
static
void computePosTrivial(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  elements,           /**< types of elements that have been packed */
   int                   nelements,          /**< the total number of elements */
   SCIP_Real*            rexts,              /**< external radii */
   SCIP_Real*            xs,                 /**< x-coordinate of circle */
   SCIP_Real*            ys,                 /**< y-coordinate of circle */
   int                   pos,                /**< position of element in the elements array */
   SCIP_Bool*            ispacked,           /**< array indicating whether an element has been packed already */
   SCIP_Real             rmax,               /**< maximum radius of elements in the pattern */
   SCIP_Real             rbound,             /**< radius of bounding circle */
   SCIP_Real             width,              /**< width of the rectangle */
   SCIP_Real             height,             /**< height of the rectangle */
   SCIP_PATTERNTYPE      patterntype,        /**< the pattern type (rectangular or circular) */
   SCIP_Real*            bestx,              /**< pointer to store the best x-coordinate */
   SCIP_Real*            besty,              /**< pointer to store the best y-coordinate */
   int                   ncalls              /**< total number of calls of the packing heuristic */
   )
{
   SCIP_Real rext = rexts[elements[pos]];
   int i;

   if( patterntype == SCIP_PATTERNTYPE_CIRCULAR )
   {
      SCIP_Real xcands[4] = {-rbound + rext, +rbound - rext, 0.0, 0.0};
      SCIP_Real ycands[4] = {0.0, 0.0, -rbound + rext, +rbound - rext};

      for( i = 0; i < 4; ++i )
         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], rbound, width, height, rmax, patterntype,
            ispacked, elements, nelements, bestx, besty, xcands[i], ycands[i], ncalls);
   }
   else
   {
      SCIP_Real xcands[4] = {rext, width - rext, rext, width - rext};
      SCIP_Real ycands[4] = {rext, rext, height - rext, height - rext};

      for( i = 0; i < 4; ++i )
         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], rbound, width, height, rmax, patterntype,
            ispacked, elements, nelements, bestx, besty, xcands[i], ycands[i], ncalls);
   }
}

/** auxiliary function for computing a candidate position between a circle and the rectangle */
static
void computePosRectangleCircle(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  elements,           /**< types of elements that have been packed */
   int                   nelements,          /**< the total number of elements */
   SCIP_Real*            rexts,              /**< external radii */
   SCIP_Real*            xs,                 /**< x-coordinate of circle */
   SCIP_Real*            ys,                 /**< y-coordinate of circle */
   int                   pos,                /**< position of element in the elements array */
   SCIP_Bool*            ispacked,           /**< array indicating whether an element has been packed already */
   SCIP_Real             rmax,               /**< maximum radius of elements in the pattern */
   SCIP_Real             width,              /**< width of the rectangle */
   SCIP_Real             height,             /**< height of the rectangle */
   SCIP_Real*            bestx,              /**< pointer to store the best x-coordinate */
   SCIP_Real*            besty,              /**< pointer to store the best y-coordinate */
   int                   ncalls              /**< total number of calls of the packing heuristic */
   )
{
   SCIP_Real rext;
   int i;

   rext = rexts[elements[pos]];

   for( i = 0; i < nelements; ++i )
   {
      SCIP_Real xfix[2] = {rext, width - rext};
      SCIP_Real yfix[2] = {rext, height - rext};
      SCIP_Real Ri;
      int k;

      if( !ispacked[i] )
         continue;

      Ri = rexts[elements[i]];

      /* fix x */
      for( k = 0; k < 2; ++k )
      {
         SCIP_Real alpha = SQR(rext + Ri) - SQR(xfix[k] - xs[i]);

         if( alpha < 0.0 )
            continue;

         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], -1.0, width, height, rmax,
            SCIP_PATTERNTYPE_RECTANGULAR, ispacked, elements, nelements, bestx, besty, xfix[k], ys[i] + SQRT(alpha), ncalls);

         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], -1.0, width, height, rmax,
            SCIP_PATTERNTYPE_RECTANGULAR, ispacked, elements, nelements, bestx, besty, xfix[k], ys[i] - SQRT(alpha), ncalls);
      }

      /* fix y */
      for( k = 0; k < 2; ++k )
      {
         SCIP_Real alpha = SQR(rext + Ri) - SQR(yfix[k] - ys[i]);

         if( alpha < 0.0 )
            continue;

         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], -1.0, width, height, rmax,
            SCIP_PATTERNTYPE_RECTANGULAR, ispacked, elements, nelements, bestx, besty, xs[i] + SQRT(alpha), yfix[k], ncalls);

         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], -1.0, width, height, rmax,
            SCIP_PATTERNTYPE_RECTANGULAR, ispacked, elements, nelements, bestx, besty, xs[i] - SQRT(alpha), yfix[k], ncalls);
      }
   }
}

/** auxiliary function for computing a candidate position between two circles */
static
void computePosCircleCircle(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  elements,           /**< types of elements that have been packed */
   int                   nelements,          /**< the total number of elements */
   SCIP_Real*            rexts,              /**< external radii */
   SCIP_Real*            xs,                 /**< x-coordinate of circle */
   SCIP_Real*            ys,                 /**< y-coordinate of circle */
   int                   pos,                /**< position of element in the elements array */
   SCIP_Bool*            ispacked,           /**< array indicating whether an element has been packed already */
   SCIP_Real             rmax,               /**< maximum radius of elements in the pattern */
   SCIP_Real             rbound,             /**< radius of bounding circle */
   SCIP_Real             width,              /**< width of the rectangle */
   SCIP_Real             height,             /**< height of the rectangle */
   SCIP_PATTERNTYPE      patterntype,        /**< the pattern type (rectangular or circular) */
   SCIP_Real*            bestx,              /**< pointer to store the best x-coordinate */
   SCIP_Real*            besty,              /**< pointer to store the best y-coordinate */
   int                   ncalls              /**< total number of calls of the packing heuristic */
   )
{
   SCIP_Real rext;
   int i;

   rext = rexts[elements[pos]];

   /* consider all pairs of already packed circles */
   for( i = 0; i < nelements - 1; ++i )
   {
      SCIP_Real alpha, a, b, h, u, v, n1, n2;
      SCIP_Real Ri;
      int j;

      if( !ispacked[i] )
         continue;

      Ri = rexts[elements[i]];

      for( j = i + 1; j < nelements; ++j )
      {
         SCIP_Real Rj;
         SCIP_Real dist;

         if( !ispacked[j] )
            continue;

         Rj = rexts[elements[j]];
         dist = SQRT(SQR(xs[i] - xs[j]) + SQR(ys[i] - ys[j]));

         /* circles are too far away */
         if( SCIPisGE(scip, dist, Ri + Rj + 2.0 * rext) )
            continue;

         a = Ri + rext;
         b = Rj + rext;
         assert(dist != 0.0);
         alpha = (dist*dist - b*b + a*a) / (2.0*dist);
         h = SQRT(a*a - alpha*alpha);
         u = xs[i] + (alpha / dist) * (xs[j] - xs[i]);
         v = ys[i] + (alpha / dist) * (ys[j] - ys[i]);
         n1 = h * ((ys[j] - ys[i]) / dist);
         n2 = h * ((xs[i] - xs[j]) / dist);
         assert(n1*n1 + n2*n2 > 0.0);

         updateBestCandidate(scip, xs, ys, rexts, rext, rbound, width, height, rmax, patterntype, ispacked, elements,
            nelements, bestx, besty, u + n1, v + n2, ncalls);
         updateBestCandidate(scip, xs, ys, rexts, rext, rbound, width, height, rmax, patterntype, ispacked, elements,
            nelements, bestx, besty, u - n1, v - n2, ncalls);
      }
   }
}

/** array to compute the score of each element */
static
void computeScores(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int*                  elements,           /**< type of each element */
   int                   nelements,          /**< total number of elements */
   SCIP_Real*            scores,             /**< array to store the score of each element */
   int                   iter                /**< iteration round */
   )
{
   SCIP_Real* rexts;
   int i;

   rexts = SCIPprobdataGetRexts(probdata);
   assert(rexts != NULL);

   for( i = 0; i < nelements; ++i )
   {
      SCIP_Real rext = rexts[elements[i]];
      /* use largest elements first */
      if( iter == 0 )
         scores[i] = rext;

      /* use smallest elements first */
      else if( iter == 1 )
         scores[i] = -rext;

      /* use [0,1] * radius */
      else if( iter <= 10 )
         scores[i] = SCIPrandomGetReal(probdata->randnumgen, 0.0, 1.0) * rext;

      /* use [-1,0] * radius */
      else if( iter <= 20 )
         scores[i] = SCIPrandomGetReal(probdata->randnumgen, -1.0, 0.0) * rext;

      /* use a random order */
      else
         scores[i] = SCIPrandomGetReal(probdata->randnumgen, 0.0, 1.0);
   }
}

/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigRingpacking)
{
   SCIPdebugMessage("free original problem data\n");
   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransRingpacking)
{
   SCIPdebugMessage("free transformed problem data\n");
   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransRingpacking)
{  /*lint --e{715}*/
   /* create transformed problem data */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->patternconss, sourcedata->cpatterns, sourcedata->cvars,
         sourcedata->ncpatterns, sourcedata->rpatterns, sourcedata->rvars, sourcedata->nrpatterns,
         sourcedata->demands, sourcedata->rints, sourcedata->rexts, sourcedata->ntypes,
         sourcedata->width, sourcedata->height) );

   /* transform pattern constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->ntypes, (*targetdata)->patternconss,
         (*targetdata)->patternconss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->ncpatterns, (*targetdata)->cvars, (*targetdata)->cvars) );
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nrpatterns, (*targetdata)->rvars, (*targetdata)->rvars) );

   /* copy statistics to transformed problem data */
   (*targetdata)->ncppatternsunknownbeg = sourcedata->ncppatternsunknownbeg;
   (*targetdata)->enumtime = sourcedata->enumtime;
   (*targetdata)->dualbound = sourcedata->dualbound;
   (*targetdata)->isdualinvalid = sourcedata->isdualinvalid;

   return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   int*                  demands,            /**< array containing the demands */
   SCIP_Real*            rints,              /**< internal radii of each ring */
   SCIP_Real*            rexts,              /**< external radii of each ring (assumed to be sorted) */
   int                   ntypes,             /**< number of different types */
   SCIP_Real             width,              /**< width of each rectangle */
   SCIP_Real             height              /**< height of each rectangle */
   )
{
   SCIP_PROBDATA* probdata;

#ifndef NDEBUG
   {
      int t;

      for( t = 0; t < ntypes -1; ++t )
         assert(rexts[t] >= rexts[t+1]);
   }
#endif

   /* create SCIP problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   /* create and set problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, NULL, NULL, NULL, 0, NULL, NULL, 0, demands, rints, rexts, ntypes, width,
         height) );
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigRingpacking) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransRingpacking) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransRingpacking) );

   /* activate pricer */
   SCIP_CALL( SCIPpricerRpaActivate(scip) );

   /* add table output */
   assert(SCIPfindTable(scip, TABLE_NAME_RPA) == NULL);
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_RPA, TABLE_DESC_RPA, TRUE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputRpa,
         NULL, TABLE_POSITION_RPA, TABLE_EARLIEST_STAGE_RPA) );

   return SCIP_OKAY;
}

/** enumerates circular patterns and creates restricted master problem */
SCIP_RETCODE SCIPprobdataSetupProblem(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* collect parameters for verification */
   SCIP_CALL( SCIPgetRealParam(scip, "ringpacking/verification/nlptilimsoft", &probdata->nlptilimsoft) );
   SCIP_CALL( SCIPgetRealParam(scip, "ringpacking/verification/heurtilimsoft", &probdata->heurtilimsoft) );
   SCIP_CALL( SCIPgetLongintParam(scip, "ringpacking/verification/nlpnodelimsoft", &probdata->nlpnodelimsoft) );
   SCIP_CALL( SCIPgetIntParam(scip, "ringpacking/verification/heuriterlimsoft", &probdata->heuriterlimsoft) );
   SCIP_CALL( SCIPgetRealParam(scip, "ringpacking/verification/totaltilimsoft", &probdata->totaltilimsoft) );

   SCIP_CALL( setupProblem(scip, probdata) );

   return SCIP_OKAY;
}

/** enumerate all non-dominated circular patterns */
SCIP_RETCODE SCIPprobdataEnumeratePatterns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Real             nlptilim,           /**< time limit for each NLP verification */
   SCIP_Real             heurtilim,          /**< time limit for each call of the heuristics */
   SCIP_Real             totaltilim,         /**< total time limit for enumeration */
   SCIP_Longint          nlpnodelim,         /**< node limit for each NLP verification */
   int                   heuriterlim         /**< iteration limit for each call of the heuristics */
   )
{
   SCIP_PATTERN* pattern;
   int* ms;
   int* nselected;
   SCIP_Real timeleft;
   int ntypes;
   int t;

   assert(probdata != NULL);
   ntypes = SCIPprobdataGetNTypes(probdata);
   assert(ntypes > 0);

   /* create data that is used for the whole enumeration algorithm */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &pattern, 0) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ms, ntypes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nselected, ntypes) );
   BMSclearMemoryArray(nselected, ntypes);
   BMSclearMemoryArray(ms, ntypes);

   /* compute time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timeleft) );
   timeleft = MAX(0.0, MIN(timeleft - SCIPgetTotalTime(scip), totaltilim)); /*lint !e666*/

   /* find all circlular patterns of each type separately */
   for( t = 0; t < ntypes; ++t )
   {
      int k;

      for( k = t+1; k < ntypes; ++k )
         ms[k] = maxCircles(scip, probdata, t, k);

      SCIPpatternSetType(pattern, t);
      SCIP_CALL( enumeratePatterns(scip, probdata, pattern, ms, nselected, nlptilim, heurtilim, nlpnodelim,
         heuriterlim, &timeleft) );
   }

   /* release memory */
   SCIPfreeBufferArray(scip, &nselected);
   SCIPfreeBufferArray(scip, &ms);
   SCIPpatternRelease(scip, &pattern);

   /* filter circular patterns */
   SCIP_CALL( filterPatterns(scip, probdata) );

   return SCIP_OKAY;
}

/** returns number of different types */
int SCIPprobdataGetNTypes(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->ntypes;
}

/** returns all external radii */
SCIP_Real* SCIPprobdataGetRexts(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->rexts;
}

/** returns all internal radii */
SCIP_Real* SCIPprobdataGetRints(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->rints;
}

/** returns all demands */
int* SCIPprobdataGetDemands(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->demands;
}

/** returns the width of each rectangle */
SCIP_Real SCIPprobdataGetWidth(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->width;
}


/** returns the height of each rectangle */
SCIP_Real SCIPprobdataGetHeight(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->height;
}

/** returns all information about circular patterns */
void SCIPprobdataGetCInfos(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN***       cpatterns,          /**< pointer to store the circular patterns (might be NULL) */
   SCIP_VAR***           cvars,              /**< pointer to store the variables corresponding circular patterns (might be NULL) */
   int*                  ncpatterns          /**< pointer to store the number of circular patterns (might be NULL) */
   )
{
   assert(probdata != NULL);

   if( cpatterns != NULL )
      *cpatterns = probdata->cpatterns;
   if( cvars != NULL )
      *cvars= probdata->cvars;
   if( ncpatterns != NULL )
      *ncpatterns = probdata->ncpatterns;
}

/** returns all information about rectangular patterns */
void SCIPprobdataGetRInfos(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN***       rpatterns,          /**< pointer to store the rectangular patterns (might be NULL) */
   SCIP_VAR***           rvars,              /**< pointer to store the variables corresponding rectangular patterns (might be NULL) */
   int*                  nrpatterns          /**< pointer to store the number of rectangular patterns (might be NULL) */
   )
{
   assert(probdata != NULL);

   if( rpatterns != NULL )
      *rpatterns = probdata->rpatterns;
   if( rvars != NULL )
      *rvars= probdata->rvars;
   if( nrpatterns != NULL )
      *nrpatterns = probdata->nrpatterns;
}

/** returns array of set pattern constraints */
SCIP_CONS** SCIPprobdataGetPatternConss(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->patternconss;
}

/** adds given variable to the problem data */
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern,            /**< pattern */
   SCIP_VAR*             var                 /**< variables to add */
   )
{
   SCIP_PATTERN* copy;

   assert(probdata != NULL);
   assert(pattern != NULL);
   assert(SCIPpatternGetPackableStatus(pattern) != SCIP_PACKABLE_NO);

   /* copy pattern */
   SCIP_CALL( SCIPpatternCopy(scip, pattern, &copy) );
   SCIPcheckPattern(scip, probdata, copy);

   if( SCIPpatternGetPatternType(pattern) == SCIP_PATTERNTYPE_CIRCULAR )
   {
      SCIP_CALL( ensureSize(scip, probdata, SCIP_PATTERNTYPE_CIRCULAR, probdata->ncpatterns + 1) );
      probdata->cpatterns[probdata->ncpatterns] = copy;
      probdata->cvars[probdata->ncpatterns] = var;
      ++(probdata->ncpatterns);
   }
   else
   {
      SCIP_CALL( ensureSize(scip, probdata, SCIP_PATTERNTYPE_RECTANGULAR, probdata->nrpatterns + 1) );
      probdata->rpatterns[probdata->nrpatterns] = copy;
      probdata->rvars[probdata->nrpatterns] = var;
      ++(probdata->nrpatterns);
   }

   /* capture variable and pattern */
   if( var != NULL )
   {
      SCIP_CALL( SCIPcaptureVar(scip, var) );
   }

   return SCIP_OKAY;
}

/** updates the dual bound */
void SCIPprobdataUpdateDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Real             dualbound           /**< new dual bound */
   )
{
   assert(probdata != NULL);

   if( !probdata->isdualinvalid && SCIPisFeasLT(scip, probdata->dualbound, dualbound) )
   {
      SCIPinfoMessage(scip, NULL, "+++++++++++++ update dual bound to %g\n", dualbound);
      probdata->dualbound = dualbound;
   }
}

/** marks that further reported dual bounds are not valid */
void SCIPprobdataInvalidateDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   if( !probdata->isdualinvalid )
   {
      SCIPinfoMessage(scip, NULL, "+++++++++++++ invalidate dual bound\n");
      probdata->isdualinvalid = TRUE;
   }
}

/** returns whether dual bound is marked to be invalid */
SCIP_Bool SCIPprobdataIsDualboundInvalid(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->isdualinvalid;
}

/** Tries to pack a list of elements into a specified boundary circle by using a simple left-first bottom-second
 *  heuristic. Returns the number of elements that could be stored and indicated which ones these are in the buffer
 *  parameter ispacked. This auxiliary method can be used both to find such a packing or to verify a certain pattern.
 */
void SCIPpackCirclesGreedy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            rexts,              /**< outer radii of elements (in original order of probdata) */
   SCIP_Real*            xs,                 /**< buffer to store the resulting x-coordinates */
   SCIP_Real*            ys,                 /**< buffer to store the resulting y-coordinates */
   SCIP_Real             rbounding,          /**< inner radius of bounding circle (ignored for rectangular patterns) */
   SCIP_Real             width,              /**< width of the rectangle */
   SCIP_Real             height,             /**< height of the rectangle */
   SCIP_Bool*            ispacked,           /**< buffer to store which elements could be packed */
   int*                  elements,           /**< the order of the elements in the pattern */
   int                   nelements,          /**< number of elements in the pattern */
   SCIP_PATTERNTYPE      patterntype,        /**< the pattern type (rectangular or circular) */
   int*                  npacked,            /**< pointer to store the number of packed elements */
   int                   ncalls              /**< total number of calls of the packing heuristic */
   )
{
   SCIP_Real rmax;
   SCIP_Bool added;
   int i;

   assert(rexts != NULL);
   assert(xs != NULL);
   assert(ys != NULL);
   assert(ispacked != NULL);
   assert(elements != NULL);
   assert(nelements > 0);
   assert(npacked != NULL);

   /* no element packed so far */
   BMSclearMemoryArray(ispacked, nelements);

   /* place first element at left-most position */
   if( patterntype == SCIP_PATTERNTYPE_CIRCULAR )
   {
      assert(rexts[elements[0]] <= rbounding);
      xs[0] = rexts[elements[0]] - rbounding;
      ys[0] = 0.0;
   }
   else
   {
      assert(2.0 * rexts[elements[0]] <= width);
      assert(2.0 * rexts[elements[0]] <= height);
      xs[0] = rexts[elements[0]];
      ys[0] = rexts[elements[0]];
   }

   /* initialize results */
   (*npacked) = 1;
   ispacked[0] = TRUE;
   added = TRUE;

   /* find max radius */
   rmax = rexts[elements[0]];
   for( i = 1; i < nelements; ++i )
   {
      if( rexts[elements[i]] > rmax )
         rmax = rexts[elements[i]];
   }

   /* iterate over all elements and try to pack them */
   while( added )
   {
      added = FALSE;

      for( i = 1; i < nelements; ++i )
      {
         SCIP_Real bestx = SCIP_INVALID;
         SCIP_Real besty = SCIP_INVALID;

         /* skip packed elements */
         if( ispacked[i] )
            continue;

         /* use trivial candidates */
         computePosTrivial(scip, elements, nelements, rexts, xs, ys, i, ispacked, rmax, rbounding, width, height,
            patterntype, &bestx, &besty, ncalls);

         /* consider circles intersection a previous circle and the boundary ring */
         if( patterntype == SCIP_PATTERNTYPE_CIRCULAR )
            computePosRingCircle(scip, elements, nelements, rexts, xs, ys, i, ispacked, rmax, rbounding, &bestx,
               &besty, ncalls);
         else
            computePosRectangleCircle(scip, elements, nelements, rexts, xs, ys, i, ispacked, rmax, width, height, &bestx,
               &besty, ncalls);

         /* consider circles that have been packed already */
         computePosCircleCircle(scip, elements, nelements, rexts, xs, ys, i, ispacked, rmax, rbounding, width, height,
            patterntype, &bestx, &besty, ncalls);

         /* pack circle if a possible position has been found */
         if( bestx != SCIP_INVALID && besty != SCIP_INVALID ) /*lint !e777*/
         {
            assert(!ispacked[i]);
            ispacked[i] = TRUE;
            xs[i] = bestx;
            ys[i] = besty;
            ++(*npacked);
            added = TRUE;
         }
      }
   }

   return;
}

/** verifies a circular pattern heuristically */
SCIP_RETCODE SCIPverifyCircularPatternHeuristic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern,            /**< pattern */
   SCIP_Real             timelim,            /**< time limit */
   int                   iterlim             /**< iteration limit */
   )
{
   SCIP_Real* rexts;
   SCIP_Real* rints;
   SCIP_Real* scores;
   SCIP_Real* xs;
   SCIP_Real* ys;
   SCIP_Bool* ispacked;
   int* elements;
   int* pos;
   SCIP_Real timestart;
   int nelements;
   int niters;
   int type;
   int i;

   assert(probdata != NULL);
   assert(pattern != NULL);
   assert(iterlim > 0);
   assert(SCIPpatternGetPatternType(pattern) == SCIP_PATTERNTYPE_CIRCULAR);
   assert(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);
   assert(SCIPpatternGetCircleType(pattern) < SCIPprobdataGetNTypes(probdata));

   /* check whether there is any time left */
   if( timelim <= 0.0 )
      return SCIP_OKAY;

   rexts = SCIPprobdataGetRexts(probdata);
   rints = SCIPprobdataGetRints(probdata);
   nelements = SCIPpatternGetNElemens(pattern);
   type = SCIPpatternGetCircleType(pattern);
   assert(type >= 0 && type < SCIPprobdataGetNTypes(probdata));

   /* pattern is empty -> set status to packable */
   if( SCIPpatternGetNElemens(pattern) == 0 )
   {
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);
      SCIPcheckPattern(scip, probdata, pattern);
      return SCIP_OKAY;
   }

   /* pattern contains only one element -> compare radii */
   if( SCIPpatternGetNElemens(pattern) == 1 )
   {
      int elemtype;

      elemtype = SCIPpatternGetElementType(pattern, 0);
      assert(elemtype >= 0 && elemtype < SCIPprobdataGetNTypes(probdata));

      /* check whether element fits into the circular pattern */
      if( SCIPisGE(scip, rints[type], rexts[elemtype]) )
      {
         SCIPpatternSetElementPos(pattern, 0, rexts[elemtype]-rints[type], 0.0);
         SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);
      }
      else
         SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_NO);

      SCIPcheckPattern(scip, probdata, pattern);
      return SCIP_OKAY;
   }

   timestart = SCIPgetTotalTime(scip);
   niters = 0;

   /* store elements in a separate array; remember positions of elements in the pattern */
   SCIP_CALL( SCIPallocBufferArray(scip, &pos, nelements) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, nelements) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ispacked, nelements) );
   SCIP_CALL( SCIPallocBufferArray(scip, &elements, nelements) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xs, nelements) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ys, nelements) );
   for( i = 0; i < nelements; ++i )
   {
      elements[i] = SCIPpatternGetElementType(pattern, i);
      ispacked[i] = FALSE;
      pos[i] = i;
   }

   /* main loop for calling heuristic verification */
   while( SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN
      && niters < iterlim
      && SCIPgetTotalTime(scip) - timestart <= timelim )
   {
      int npacked;

      /* compute scores depending on iteration counter */
      computeScores(scip, probdata, elements, nelements, scores, niters);

      /* sort elements in non-increasing order */
      SCIPsortDownRealIntInt(scores, elements, pos, nelements);

      /* call heuristic */
      SCIPpackCirclesGreedy(scip, rexts, xs, ys, rints[type], SCIPprobdataGetWidth(probdata),
         SCIPprobdataGetHeight(probdata), ispacked, elements, nelements, SCIP_PATTERNTYPE_CIRCULAR, &npacked, niters);

      /* check whether all elements could have been packed */
      if( npacked == nelements )
      {
         for( i = 0; i < nelements; ++i )
         {
            assert(elements[i] == SCIPpatternGetElementType(pattern, pos[i]));
            SCIPpatternSetElementPos(pattern, pos[i], xs[i], ys[i]);
         }
         SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);

         SCIPdebugMsg(scip, "heuristic verified pattern after %d iterations\n", niters + 1);
      }

      ++niters;
   }

   SCIPcheckPattern(scip, probdata, pattern);

   /* free memory */
   SCIPfreeBufferArray(scip, &ys);
   SCIPfreeBufferArray(scip, &xs);
   SCIPfreeBufferArray(scip, &elements);
   SCIPfreeBufferArray(scip, &ispacked);
   SCIPfreeBufferArray(scip, &scores);
   SCIPfreeBufferArray(scip, &pos);

   return SCIP_OKAY;
}

/** verifies a circular pattern via a verification NLP */
SCIP_RETCODE SCIPverifyCircularPatternNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern,            /**< pattern */
   SCIP_Real             timelim,            /**< time limit */
   SCIP_Longint          nodelim             /**< node limit */
   )
{
   SCIP* subscip;
   SCIP_CONS* cons;
   SCIP_VAR** xvars;
   SCIP_VAR** yvars;
   SCIP_VAR* quadvars1[6];
   SCIP_VAR* quadvars2[6];
   SCIP_Real quadcoefs[6];
   SCIP_Real* rexts;
   SCIP_Real* rints;
   char name[SCIP_MAXSTRLEN];
   int nelems;
   int type;
   int k;

   assert(probdata != NULL);
   assert(pattern != NULL);
   assert(SCIPpatternGetPatternType(pattern) == SCIP_PATTERNTYPE_CIRCULAR);
   assert(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* check whether there is any time left */
   if( timelim <= 0.0 )
      return SCIP_OKAY;

   rexts = SCIPprobdataGetRexts(probdata);
   rints = SCIPprobdataGetRints(probdata);
   type = SCIPpatternGetCircleType(pattern);
   nelems = SCIPpatternGetNElemens(pattern);

   /* set up the sub-SCIP */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPcreateProbBasic(subscip, "verify") );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* allocate memory for (x,y) variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &xvars, nelems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &yvars, nelems) );

   /* set feasibility emphasis settings */
   SCIP_CALL( SCIPsetEmphasis(subscip, SCIP_PARAMEMPHASIS_FEASIBILITY, TRUE) );

   /* set working limit */
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/solutions", 1) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelim) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelim) );

#ifndef SCIP_DEBUG
   SCIPsetMessagehdlrQuiet(subscip, TRUE);
#endif

   /* create (x,y) variables */
   for( k = 0; k < nelems; ++k )
   {
      int elemtype;

      elemtype = SCIPpatternGetElementType(pattern, k);
      assert(elemtype >= 0 && elemtype < SCIPprobdataGetNTypes(probdata));

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", k);
      SCIP_CALL( SCIPcreateVarBasic(subscip, &xvars[k], name, rexts[elemtype] - rints[type],
            rints[type] - rexts[elemtype], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(subscip, xvars[k]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d", k);
      SCIP_CALL( SCIPcreateVarBasic(subscip, &yvars[k], name, rexts[elemtype] - rints[type],
            rints[type] - rexts[elemtype], 1.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(subscip, yvars[k]) );
   }

   /* create non-overlapping constraints */
   for( k = 0; k < nelems; ++k )
   {
      int elemtype1;
      int l;

      elemtype1 = SCIPpatternGetElementType(pattern, k);
      assert(elemtype1 >= 0 && elemtype1 < SCIPprobdataGetNTypes(probdata));

      for( l = k + 1; l < nelems; ++l )
      {
         int elemtype2;

         elemtype2 = SCIPpatternGetElementType(pattern, l);
         assert(elemtype2 >= 0 && elemtype2 < SCIPprobdataGetNTypes(probdata));

         quadvars1[0] = xvars[k]; quadvars2[0] = xvars[k]; quadcoefs[0] =  1.0;
         quadvars1[1] = xvars[k]; quadvars2[1] = xvars[l]; quadcoefs[1] = -2.0;
         quadvars1[2] = xvars[l]; quadvars2[2] = xvars[l]; quadcoefs[2] =  1.0;
         quadvars1[3] = yvars[k]; quadvars2[3] = yvars[k]; quadcoefs[3] =  1.0;
         quadvars1[4] = yvars[k]; quadvars2[4] = yvars[l]; quadcoefs[4] = -2.0;
         quadvars1[5] = yvars[l]; quadvars2[5] = yvars[l]; quadcoefs[5] =  1.0;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "over_%d_%d", k, l);
         SCIP_CALL( SCIPcreateConsBasicQuadratic(subscip, &cons, name, 0, NULL, NULL, 6, quadvars1, quadvars2,
               quadcoefs, SQR(rexts[elemtype1] + rexts[elemtype2]), SCIPinfinity(subscip)) );

         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      }
   }

   /* create non-overlapping constraints with outer ring */
   for( k = 0; k < nelems; ++k )
   {
      int elemtype;

      elemtype = SCIPpatternGetElementType(pattern, k);
      assert(elemtype >= 0 && elemtype < SCIPprobdataGetNTypes(probdata));

      quadvars1[0] = xvars[k]; quadvars2[0] = xvars[k]; quadcoefs[0] = 1.0;
      quadvars1[1] = yvars[k]; quadvars2[1] = yvars[k]; quadcoefs[1] = 1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "bound_%d", k);
      SCIP_CALL( SCIPcreateConsBasicQuadratic(subscip, &cons, name, 0, NULL, NULL, 2, quadvars1, quadvars2, quadcoefs,
            0.0, SQR(rints[type] - rexts[elemtype])) );

      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
   }

   /* sort circles in x direction if they have the same type */
   for( k = 0; k < nelems - 1; ++k )
   {
      int elemtype1;
      int l;

      elemtype1 = SCIPpatternGetElementType(pattern, k);
      assert(elemtype1 >= 0 && elemtype1 < SCIPprobdataGetNTypes(probdata));

      for( l = k + 1; l < nelems; ++l )
      {
         int elemtype2;

         elemtype2 = SCIPpatternGetElementType(pattern, k+1);
         assert(elemtype2 >= 0 && elemtype2 < SCIPprobdataGetNTypes(probdata));

         if( elemtype1 != elemtype2 )
            continue;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sortx_%d_%d", k, l);
         SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &cons, name, 0, NULL, NULL, -SCIPinfinity(subscip), 0.0) );
         SCIP_CALL( SCIPaddCoefLinear(subscip, cons, xvars[k], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(subscip, cons, xvars[l], -1.0) );

         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      }
   }

   /* solve verification NLP */
   SCIPdebugMsg(scip, "--------------------- SOLVE VERIFICATION NLP -------------------\n");
   SCIP_CALL( SCIPsolve(subscip) );
   SCIPdebugMsg(scip, "----------------------------------------------------------------\n");

   SCIPdebugMsg(scip, "result of verification NLP: nsols=%d solstat=%d\n",
      SCIPgetNSols(subscip), SCIPgetStatus(subscip));

   /* check whether a solution could be found or whether the problem is proven to be infeasible */
   if( SCIPgetNSols(subscip) > 0 )
   {
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);

      for( k = 0; k < nelems; ++k )
      {
         SCIP_Real solx = SCIPgetSolVal(subscip, SCIPgetBestSol(subscip), xvars[k]);
         SCIP_Real soly = SCIPgetSolVal(subscip, SCIPgetBestSol(subscip), yvars[k]);

         SCIPpatternSetElementPos(pattern, k, solx, soly);
      }

      SCIPcheckPattern(scip, probdata, pattern);
   }
   else if( SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE )
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_NO);

   /* free all variables */
   for( k = 0; k < nelems; ++k )
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &yvars[k]) );
      SCIP_CALL( SCIPreleaseVar(subscip, &xvars[k]) );
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &yvars);
   SCIPfreeBufferArray(scip, &xvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/** check a pattern for consistency */
void SCIPcheckPattern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern             /**< pattern */
   )
{  /*lint --e{715}*/
#ifndef NDEBUG
   SCIP_Real* rexts;
   SCIP_Real* rints;
   SCIP_Real width;
   SCIP_Real height;
   int i;

   assert(probdata != NULL);
   assert(pattern != NULL);

   rexts = SCIPprobdataGetRexts(probdata);
   rints = SCIPprobdataGetRints(probdata);
   width = SCIPprobdataGetWidth(probdata);
   height = SCIPprobdataGetHeight(probdata);

   /* check types */
   for( i = 0; i < SCIPpatternGetNElemens(pattern); ++i )
   {
      int type = SCIPpatternGetElementType(pattern, i);

      assert(type >= 0);
      assert(type < SCIPprobdataGetNTypes(probdata));
   }

   /* check positions iff packable */
   if( SCIPpatternGetPackableStatus(pattern) != SCIP_PACKABLE_YES )
      return;

   for( i = 0; i < SCIPpatternGetNElemens(pattern); ++i )
   {
      SCIP_Real xi = SCIPpatternGetElementPosX(pattern, i);
      SCIP_Real yi = SCIPpatternGetElementPosY(pattern, i);
      int typei = SCIPpatternGetElementType(pattern, i);
      int j;

      /* check distance between circles */
      for( j = i + 1; j < SCIPpatternGetNElemens(pattern); ++j )
      {
         SCIP_Real xj = SCIPpatternGetElementPosX(pattern, j);
         SCIP_Real yj = SCIPpatternGetElementPosY(pattern, j);
         int typej = SCIPpatternGetElementType(pattern, j);

         assert(SCIPisFeasGE(scip, SQRT(SQR(xi - xj) + SQR(yi - yj)), rexts[typei] + rexts[typej]));
      }

      /* check distance to boundary */
      if( SCIPpatternGetPatternType(pattern) == SCIP_PATTERNTYPE_CIRCULAR )
      {
         SCIP_Real distance = SQRT(SQR(xi) + SQR(yi));
         int patterntype = SCIPpatternGetCircleType(pattern);

         assert(patterntype >= 0);
         assert(patterntype < SCIPprobdataGetNTypes(probdata));
         assert(SCIPisFeasLE(scip, distance, rints[patterntype] - rexts[typei]));
      }
      else
      {
         assert(SCIPisFeasGE(scip, xi, rexts[typei]));
         assert(SCIPisFeasLE(scip, xi, width - rexts[typei]));
         assert(SCIPisFeasGE(scip, yi, rexts[typei]));
         assert(SCIPisFeasLE(scip, yi, height - rexts[typei]));
      }
   }
#endif
}

/**@} */
