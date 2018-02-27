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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_ringpacking.c
 * @brief  Problem data for ringpacking problem
 * @author Benjamin Mueller
 *
 * This file handles the main problem data used in that project. For more details see \ref PROBLEMDATA page.
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

      type = SCIPpatternGetType(pattern);
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
   SCIP_PATTERN*         q                   /**< pattern */
   )
{
   SCIP_Bool res1;
   SCIP_Bool res2;
   int np;
   int nq;
   int ip;
   int iq;

   assert(p != NULL);
   assert(q != NULL);

   /* pattern can only dominate if they are of the same type */
   if( SCIPpatternGetType(p) != SCIPpatternGetType(q) )
      return 0;

   np = SCIPpatternGetNElemens(p);
   nq = SCIPpatternGetNElemens(q);
   res1 = np >= nq;
   res2 = nq >= np;
   ip = 0;
   iq = 0;

   while( ip < np && iq < nq && (res1 || res2) )
   {
      int tp = SCIPpatternGetElementType(p, ip);
      int tq = SCIPpatternGetElementType(q, iq);

      if( tp == tq )
      {
         ++ip;
         ++iq;
      }
      else if( tp < tq )
      {
         res2 = FALSE;
         ++ip;
      }
      else
      {
         res1 = FALSE;
         ++iq;
      }
   }

   assert(!res1 || !res2);

   if( res1 && (SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_YES || SCIPpatternGetPackableStatus(q) == SCIP_PACKABLE_UNKNOWN) )
      return -1;
   else if( res2 && (SCIPpatternGetPackableStatus(q) == SCIP_PACKABLE_YES || SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_UNKNOWN) )
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
   int ncpatterns;
   int i;

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

         res = isPatternDominating(p, q);

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
   int                   heuriterlim         /**< iteration limit for each call of the heuristics */
   )
{
   SCIP_Real* rexts;
   SCIP_Real* _rints;
   SCIP_Real totaltimelim;
   SCIP_Real maxvolume;
   SCIP_Real volume;
   int ntypes;
   int type;
   int lasttype;

   assert(ms != NULL);
   assert(pattern != NULL);

   type = SCIPpatternGetType(pattern);
   assert(type >= 0 && type < SCIPprobdataGetNTypes(probdata));

   /* get problem data */
   rexts = SCIPprobdataGetRexts(probdata);
   _rints = SCIPprobdataGetRints(probdata);
   ntypes = SCIPprobdataGetNTypes(probdata);
   lasttype = ntypes -1;
   volume = 0.0;
   maxvolume = SQR(_rints[SCIPpatternGetType(pattern)]) * M_PI; /*lint !e666*/

   /* check whether there is enough time left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &totaltimelim) );

   /* main loop */
   while( !SCIPisStopped(scip) )
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
         timelim = MIN(heurtilim, totaltimelim - SCIPgetSolvingTime(scip)); /*lint !e666*/

         /* verify pattern */
         SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, timelim, heuriterlim) );

         /*
          * try to verify with NLP
          */
         if( SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN )
         {
            /* compute time limit */
            timelim = MIN(nlptilim, totaltimelim - SCIPgetSolvingTime(scip)); /*lint !e666*/

            /* verify pattern */
            SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, timelim, nlpnodelim) );
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

   /* filter circular pattern */
   SCIP_CALL( filterPatterns(scip, probdata) );

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
   SCIP_Real* _rints;
   int* demands;
   SCIP_Real heurtilim;
   SCIP_Real nlptilim;
   SCIP_Real dualbound;
   SCIP_Real minrext;
   SCIP_Real volume;
   SCIP_Longint nlpnodelim;
   int heuriterlim;
   int ntypes;
   int p;
   int t;

   assert(probdata != NULL);
   assert(SCIPprobdataGetNTypes(probdata) > 0);

   /* set objective sense; tell SCIP that the objective will be always integral */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   SCIP_CALL( SCIPsetObjIntegral(scip) );

   /* collect working limits */
   SCIP_CALL( SCIPgetRealParam(scip, "ringpacking/verification/nlptilimsoft", &nlptilim) );
   SCIP_CALL( SCIPgetLongintParam(scip, "ringpacking/verification/nlpnodelimsoft", &nlpnodelim) );
   SCIP_CALL( SCIPgetRealParam(scip, "ringpacking/verification/heurtilimsoft", &heurtilim) );
   SCIP_CALL( SCIPgetIntParam(scip, "ringpacking/verification/heuriterlimsoft", &heuriterlim) );

   /* get problem data */
   ntypes = SCIPprobdataGetNTypes(probdata);
   rexts = SCIPprobdataGetRexts(probdata);
   _rints = SCIPprobdataGetRints(probdata);
   demands = SCIPprobdataGetDemands(probdata);

   /* compute all non-dominated circular patterns */
   SCIP_CALL( SCIPprobdataEnumeratePatterns(scip, probdata, nlptilim, heurtilim, nlpnodelim, heuriterlim) );
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
         if( SCIPpatternGetType(pattern) == t )
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

      type = SCIPpatternGetType(pattern);
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
      if( SCIPisFeasLT(scip, _rints[t], minrext) )
         vol = M_PI * SQR(rexts[t]);
      else
         vol = M_PI * (SQR(rexts[t]) - SQR(_rints[t]));

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
   SCIP_Real dualbound;
   int* demands;
   int ntypes;
   int nrings;
   int t;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   ntypes = SCIPprobdataGetNTypes(probdata);
   demands = SCIPprobdataGetDemands(probdata);
   nrings = 0;

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
      nrings += demands[t];

   SCIPinfoMessage(scip, file, "Ringpacking        : %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
      "dual", "ntypes", "nrings", "width", "height", "CP", "CP_unk", "CP_unk_end" ,"CP_no", "RP", "CP_time");

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
   SCIP_PATTERNTYPE      patterntype,        /**< pattern type */
   SCIP_Bool*            ispacked,           /**< array indicating which elements are already packed */
   int*                  elements,           /**< the order of the elements in the pattern */
   int                   nelements,          /**< the total number of elements */
   SCIP_Real*            bestx,              /**< buffer to update best x-coordinate */
   SCIP_Real*            besty,              /**< buffer to update best y-coordinate */
   SCIP_Real             x,                  /**< x-coordinate of a candidate point */
   SCIP_Real             y                   /**< y-coordinate of a candidate point */
   )
{

   int i;

   /* candidate is not valid -> skip */
   if( x == SCIP_INVALID || y == SCIP_INVALID )
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

   /* check whether the candidate is better than the best known candidate */
   if( *bestx == SCIP_INVALID || *besty == SCIP_INVALID
      || SCIPisLT(scip, x, *bestx) || (SCIPisEQ(scip, x, *bestx) && SCIPisLT(scip, y, *besty)) )
   {
      *bestx = x;
      *besty = y;
   }
}

/** Auxiliary function to compute the left-most, bottom-second intersection point of two circles.
 *  Optionally, a bounding circle in which the point has to lie can be specified through parameter rbound
 *  If rbound is infinity, all points are valid. Results will be infinity if computation failed.
 */
static
void computeIntersectionCircles(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             x1,                 /**< x-coordinate of first circle */
   SCIP_Real             y1,                 /**< y-coordinate of first circle */
   SCIP_Real             r1,                 /**< radius of first circle */
   SCIP_Real             x2,                 /**< x-coordinate of second circle */
   SCIP_Real             y2,                 /**< y-coordinate of second circle */
   SCIP_Real             r2,                 /**< radius of second circle */
   SCIP_Real             rbound,             /**< radius of bounding circle */
   SCIP_Real*            xres,               /**< buffer for x-coordinate of intersection point */
   SCIP_Real*            yres                /**< buffer for y-coordinate of intersection point */
)
{
   SCIP_Real distsqr;
   SCIP_Real r1sqr;
   SCIP_Real r2sqr;
   SCIP_Real sqrtterm;

   assert(xres != NULL);
   assert(yres != NULL);
   assert(SCIPisGE(scip, r1, 0.0));
   assert(SCIPisGE(scip, r2, 0.0));
   assert(SCIPisGE(scip, rbound, 0.0));

   distsqr = SQR(x2 - x1) + SQR(y2 - y1);
   assert(distsqr != 0.0);

   /* if the distance is too large, the circles don't intersect */
   if( SCIPisGT(scip, distsqr, SQR(r1 + r2)) )
   {
      *xres = SCIPinfinity(scip);
      *yres = SCIPinfinity(scip);
      return;
   }

   r1sqr = SQR(r1);
   r2sqr = SQR(r2);
   sqrtterm = 0.5 * SQRT(2 * (r1sqr + r2sqr) / distsqr - SQR(r1sqr - r2sqr) / SQR(distsqr) - 1);

   *xres = 0.5 * (x1 + x2) + 0.5 * (r1sqr - r2sqr) * (x2 - x1) / distsqr;
   *yres = 0.5 * (y1 + y2) + 0.5 * (r1sqr - r2sqr) * (y2 - y1) / distsqr;

   /* choose intersection according to left-first, bottom-second order */
   if( (!SCIPisEQ(scip, y1, y2) && SCIPisGT(scip, sqrtterm * (y2 - y1), 0.0))
       || (SCIPisEQ(scip, y1, y2) && SCIPisGT(scip, sqrtterm * (x2 - x1), 0.0)) )
   {
      sqrtterm *= -1.0;
   }

   *xres += sqrtterm * (y2 - y1);
   *yres += sqrtterm * (x2 - x1);

   /* if point does not lie in the bounding circle, try the other one */

   if( !SCIPisInfinity(scip, rbound) && SCIPisGT(scip, SQR(*xres) + SQR(*yres), SQR(rbound)) )
   {
      *xres -= 2.0 * sqrtterm * (y2 - y1);
      *yres -= 2.0 * sqrtterm * (x2 - x1);

      if( SCIPisGT(scip, SQR(*xres) + SQR(*yres), SQR(rbound)) )
      {
         *xres = SCIPinfinity(scip);
         *yres = SCIPinfinity(scip);
      }
   }
}

/** auxiliary function for computing the intersection point */
static
void computeIntersectionRingCircle(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  elements,           /**< types of elements that have been packed */
   int                   nelements,          /**< the total number of elements */
   SCIP_Real*            rexts,              /**< external radii */
   SCIP_Real*            xs,                 /**< x-coordinate of circle */
   SCIP_Real*            ys,                 /**< y-coordinate of circle */
   int                   pos,                /**< position of element in the elements array */
   SCIP_Bool*            ispacked,           /**< array indicating whether an element has been packed already */
   SCIP_Real             rbound,             /**< radius of bounding circle */
   int                   npacked,            /**< number of elements packed */
   SCIP_Real*            bestx,              /**< pointer to store the best x-coordinate */
   SCIP_Real*            besty               /**< pointer to store the best y-coordinate */
   )
{
   int i;

   /* consider already packed patterns */
   for( i = 0; i < nelements; ++i )
   {
      SCIP_Real tmpx;
      SCIP_Real tmpy;
      SCIP_Real alpha, a, b, c, h, u, v, n1, n2;

      /* only consider packed elements */
      if( !ispacked[i] )
         continue;

      c = SQRT(xs[i]*xs[i] + ys[i]*ys[i]);

      /* inner ring is too far away from boundary or both rings can not fit */
      if( !SCIPisGE(scip, c + rexts[elements[i]] + 2.0*rexts[elements[pos]], rbound) || SCIPisGT(scip, rexts[elements[pos]] + rexts[elements[i]], rbound) )
         continue;

      a = rexts[elements[pos]] + rexts[elements[i]];
      b = rbound - rexts[elements[pos]];

      /* if a ring is in the center than there are infinitely many solutions; take an arbitrary point */
      if( SCIPisZero(scip, c) )
      {
         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], rbound, -1.0, -1.0, SCIP_PATTERNTYPE_CIRCULAR,
            ispacked, elements, nelements, bestx, besty, -rbound + rexts[elements[pos]], 0.0);
         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], rbound, -1.0, -1.0, SCIP_PATTERNTYPE_CIRCULAR,
            ispacked, elements, nelements, bestx, besty, +rbound - rexts[elements[pos]], 0.0);
      }
      else
      {
         assert(c != 0.0);
         alpha = (c*c - b*b + a*a) / (2*c);
         h = SQRT(a*a - alpha*alpha);
         u = (c - alpha) * xs[i] / c;
         v = (c - alpha) * ys[i] / c;
         n1 = h * (v / SQRT(v*v + u*u));
         n2 = h * (-u / SQRT(v*v + u*u));

         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], rbound, -1.0, -1.0, SCIP_PATTERNTYPE_CIRCULAR,
            ispacked, elements, nelements, bestx, besty, u + n1, v + n2);
         updateBestCandidate(scip, xs, ys, rexts, rexts[elements[pos]], rbound, -1.0, -1.0, SCIP_PATTERNTYPE_CIRCULAR,
            ispacked, elements, nelements, bestx, besty, u - n1, v - n2);
      }
   }
}

static
void computeIntersectionRectangleCircle(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             xcircle,            /**< x-coordinate of circle midpoint */
   SCIP_Real             ycircle,            /**< y-coordinate of circle midpoint */
   SCIP_Real             rcircle,            /**< radius of circle */
   SCIP_Real             xrect,              /**< x-coordinate of lower left corner of rectangle */
   SCIP_Real             yrect,              /**< y-coordinate of lower left corner of rectangle */
   SCIP_Real             wrect,              /**< width of rectangle */
   SCIP_Real             hrect,              /**< height of rectangle */
   SCIP_Real*            xres,               /**< buffer for x-coordinate of intersection point */
   SCIP_Real*            yres                /**< buffer for y-coordinate of intersection point */
   )
{
   assert(xres != NULL);
   assert(yres != NULL);
   assert(SCIPisGE(scip, rcircle, 0.0));
   assert(SCIPisGE(scip, xrect, 0.0));
   assert(SCIPisGE(scip, yrect, 0.0));

   /* TODO: implement */
}

/** Tries to pack a list of elements into a specified boundary circle by using a simple left-first bottom-second
 *  heuristic. Returns the number of elements that could be stored and indicated which ones these are in the buffer
 *  parameter ispacked. This auxiliary method can be used both to find such a packing or to verify a certain pattern.
 */
static
int packCirclesHeuristically(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            rexts,              /**< outer radii of elements (in original order of probdata) */
   SCIP_Real*            xs,                 /**< buffer to store the resulting x-coordinates */
   SCIP_Real*            ys,                 /**< buffer to store the resulting y-coordinates */
   SCIP_Real             rbounding,          /**< inner radius of bounding circle (ignored for rectangular patterns) */
   SCIP_Real             wbounding,          /**< width of bounding rectangular (ignored for circular patterns) */
   SCIP_Real             hbounding,          /**< height of bounding rectangular (ignored for circular patterns) */
   SCIP_Bool*            ispacked,           /**< buffer to store which elements could be packed */
   int*                  elements,           /**< the order of the elements in the pattern */
   int                   nelements,          /**< number of elements in the pattern */
   SCIP_PATTERNTYPE      patterntype         /**< the pattern type (rectangular or circular) */
   )
{
   int npacked;
   int i;
   int j;
   int k;

   assert(rexts != NULL);
   assert(xs != NULL);
   assert(ys != NULL);
   assert(ispacked != NULL);
   assert(elements != NULL);
   assert(nelements > 0);
   assert(patterntype == SCIP_PATTERNTYPE_CIRCULAR || (wbounding > 0.0 && hbounding > 0.0));
   assert(patterntype == SCIP_PATTERNTYPE_RECTANGULAR || rbounding >= 0.0);

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
      assert(2.0 * rexts[elements[0]] <= wbounding);
      assert(2.0 * rexts[elements[0]] <= hbounding);
      xs[0] = rexts[elements[0]];
      ys[0] = rexts[elements[0]];
   }

   /* initialize results */
   npacked = 1;
   ispacked[0] = TRUE;

   /* iterate over all elements and try to pack them */
   for( i = 1; i < nelements; ++i )
   {
      SCIP_Real xbest = SCIP_INVALID;
      SCIP_Real ybest = SCIP_INVALID;
      SCIP_Real rcurr = rexts[elements[i]];

      /* consider circles intersection a previous circle and the boundary ring */
      computeIntersectionRingCircle(scip, elements, nelements, rexts, xs, ys, i, ispacked, rbounding, npacked, &xbest, &ybest);

      if( xbest != SCIP_INVALID && ybest != SCIP_INVALID )
      {
         assert(!ispacked[i]);
         ispacked[i] = TRUE;
         xs[i] = xbest;
         ys[i] = ybest;
         ++npacked;
      }

//      /* for all combinations of two packed circles (and boundary) try to place next element as close as possible */
//      for( j = 0; j < i; ++j )
//      {
//         for( k = j+1; k < i; ++k )
//         {
//            /* compute position implied by circles elem1 and elem2 */
//            computeIntersectionCircles(scip, xvalues[j], yvalues[j], rexts[elements[j]], xvalues[k], yvalues[k],
//               rexts[elements[k]], rbounding, &xtmp, &ytmp);
//
//            /* check whether new position has higher priority and is valid w.r.t already packed circles */
//            if( (SCIPisLT(scip, xtmp, xbest) || (SCIPisEQ(scip, xtmp, xbest) && SCIPisLT(scip, ytmp, ybest)))
//                && isValidPosition(scip, xvalues, yvalues, rexts, xtmp, ytmp, rcurr, ispacked, elements, nelements) )
//            {
//               xbest = xtmp;
//               ybest = ytmp;
//            }
//         }
//
//         /* compute position implied by circle elements[j] and boundary */
//         if( patterntype == SCIP_PATTERNTYPE_CIRCULAR )
//         {
//            computeIntersectionCircles(scip, xvalues[j], yvalues[j], rexts[elements[j]], 0.0, 0.0, rbounding - rcurr,
//               SCIPinfinity(scip), &xtmp, &ytmp);
//         }
//         else
//         {
//            /* TODO: call computeIntersectionRectangle */
//         }
//
//         /* check whether new position has higher priority and is valid w.r.t already packed circles */
//         if( (SCIPisLT(scip, xtmp, xbest) || (SCIPisEQ(scip, xtmp, xbest) && SCIPisLT(scip, ytmp, ybest)))
//             && isValidPosition(scip, xvalues, yvalues, rexts, xtmp, ytmp, rcurr, ispacked, elements, nelements) )
//         {
//            xbest = xtmp;
//            ybest = ytmp;
//         }
//      }

//      /* if a valid position could be found, place the element there and update results */
//      if( !SCIPisInfinity(scip, xbest) )
//      {
//         xs[i] = xbest;
//         ys[i] = ybest;
//         ispacked[i] = TRUE;
//         ++npacked;
//      }
   }

   return npacked;
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

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsolRingpacking)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsolRingpacking)
{  /*lint --e{715}*/
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
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolRingpacking) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolRingpacking) );

   /* activate pricer */
   SCIP_CALL( SCIPpricerRingpackingActivate(scip) );

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

   SCIP_CALL( setupProblem(scip, probdata) );

   return SCIP_OKAY;
}

/** enumerate all non-dominated circular patterns */
SCIP_RETCODE SCIPprobdataEnumeratePatterns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_Real             nlptilim,           /**< time limit for each NLP verification */
   SCIP_Real             heurtilim,          /**< time limit for each call of the heuristics */
   SCIP_Longint          nlpnodelim,         /**< node limit for each NLP verification */
   int                   heuriterlim         /**< iteration limit for each call of the heuristics */
   )
{
   SCIP_PATTERN* pattern;
   int* ms;
   int* nselected;
   int ntypes;
   int t;

   assert(probdata != NULL);
   ntypes = SCIPprobdataGetNTypes(probdata);
   assert(ntypes > 0);

   probdata->enumtime = -SCIPgetTotalTime(scip);

   /* create data that is used for the whole enumeration algorithm */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &pattern, 0) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ms, ntypes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nselected, ntypes) );
   BMSclearMemoryArray(nselected, ntypes);
   BMSclearMemoryArray(ms, ntypes);

   /* find all circlular patterns of each type separately */
   for( t = 0; t < ntypes; ++t )
   {
      int k;

      for( k = t+1; k < ntypes; ++k )
         ms[k] = maxCircles(scip, probdata, t, k);

      SCIPpatternSetType(pattern, t);
      SCIP_CALL( enumeratePatterns(scip, probdata, pattern, ms, nselected,
            nlptilim, heurtilim, nlpnodelim, heuriterlim) );
   }

   /* release memory */
   SCIPfreeBufferArray(scip, &nselected);
   SCIPfreeBufferArray(scip, &ms);
   SCIPpatternRelease(scip, &pattern);

   /* update statistics */
   probdata->enumtime += SCIPgetTotalTime(scip);
   probdata->ncppatternsunknownbeg = getNCPatterns(scip, probdata, SCIP_PACKABLE_UNKNOWN);

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
   assert(SCIPpatternGetType(pattern) < SCIPprobdataGetNTypes(probdata));

   /* check whether there is any time left */
   if( timelim <= 0.0 )
      return SCIP_OKAY;

   rexts = SCIPprobdataGetRexts(probdata);
   rints = SCIPprobdataGetRints(probdata);
   nelements = SCIPpatternGetNElemens(pattern);
   type = SCIPpatternGetType(pattern);
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

      /* compute scores TODO use callback functions */
      for( i = 0; i < nelements; ++i )
      {
         int elemtype = SCIPpatternGetElementType(pattern, pos[i]);
         scores[ pos[i] ] = 0.0;
      }

      /* sort elements */
      SCIPsortRealIntInt(scores, elements, pos, nelements);

      /* call heuristic */
      npacked = packCirclesHeuristically(scip, rexts, xs, ys, rints[type], SCIPprobdataGetWidth(probdata),
         SCIPprobdataGetHeight(probdata), ispacked, elements, nelements, SCIP_PATTERNTYPE_CIRCULAR);

      /* check whether all elements could have been packed */
      if( npacked == nelements )
      {
         for( i = 0; i < nelements; ++i )
         {
            assert(elements[i] == SCIPpatternGetElementType(pattern, pos[i]));
            SCIPpatternSetElementPos(pattern, pos[i], xs[i], ys[i]);
         }
         SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);

         SCIPdebugMsg(scip, "heuristic verified pattern after %d iterations\n", niters);
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
   type = SCIPpatternGetType(pattern);
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

/** check whether a pattern for consistency */
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
         int patterntype = SCIPpatternGetType(pattern);

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
