/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_rinpacking.c
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
   for( i = 0; i < SCIPprobdataGetNTypes(*probdata); ++i )
   {
      assert((*probdata)->patternconss[i] != NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->patternconss[i]));
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->patternconss, (*probdata)->ntypes);

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->rexts, (*probdata)->ntypes);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->rints, (*probdata)->ntypes);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->demands, (*probdata)->ntypes);

   /* release circular patterns */
   for( i = 0; i < (*probdata)->ncpatterns; ++i )
   {
      assert((*probdata)->cpatterns[i] != NULL);
      assert((*probdata)->cvars[i] != NULL);

      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->cvars[i]) );
      SCIPpatternRelease(scip, &(*probdata)->cpatterns[i]);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->cpatterns, (*probdata)->cpatternsize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->cvars, (*probdata)->cpatternsize);
   (*probdata)->cpatternsize = 0;
   (*probdata)->ncpatterns = 0;

   /* release rectangular patterns */
   for( i = 0; i < (*probdata)->nrpatterns; ++i )
   {
      assert((*probdata)->rpatterns[i] != NULL);
      assert((*probdata)->rvars[i] != NULL);

      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->rvars[i]) );
      SCIPpatternRelease(scip, &(*probdata)->rpatterns[i]);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->rpatterns, (*probdata)->rpatternsize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->rvars, (*probdata)->rpatternsize);
   (*probdata)->rpatternsize = 0;
   (*probdata)->nrpatterns = 0;

   SCIPfreeBlockMemory(scip, probdata);
   SCIP_CALL( SCIPsetProbData(scip, NULL) );

   return SCIP_OKAY;
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

/** computes all non-dominated circular patterns and stores them into the problem data */
static
SCIP_RETCODE computeCircularPatterns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int* demands;
   int ntypes;
   int t;

   assert(probdata != NULL);

   ntypes = SCIPprobdataGetNTypes(probdata);
   assert(ntypes > 0);

   demands = SCIPprobdataGetDemands(probdata);
   assert(demands != NULL);

   /* create an empty circular pattern for each type */
   for( t = 0; t < ntypes; ++t )
   {
      SCIP_PATTERN* pattern;
      SCIP_VAR* var;

      SCIP_CALL( SCIPpatternCreateCircular(scip, &pattern, ntypes, t) );
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_UNKNOWN);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "b%d", t);
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, (SCIP_Real)demands[t], 0.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, var) );

      /* add variable and pattern to the problem data */
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, pattern, var, SCIP_PATTERNTYPE_CIRCULAR) );

      SCIP_CALL( SCIPreleaseVar(scip, &var) );
      SCIPpatternRelease(scip, &pattern);
   }

   /* create an empty circular pattern for each type */
   for( t = 0; t < ntypes; ++t )
   {
      SCIP_PATTERN* pattern;
      SCIP_VAR* var;

      SCIP_CALL( SCIPpatternCreateCircular(scip, &pattern, ntypes, t) );
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", t);
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, (SCIP_Real)demands[t], 0.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, var) );

      /* add variable and pattern to the problem data */
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, pattern, var, SCIP_PATTERNTYPE_CIRCULAR) );

      SCIP_CALL( SCIPreleaseVar(scip, &var) );
      SCIPpatternRelease(scip, &pattern);
   }

   /* create an empty circular pattern for each type */
   for( t = 0; t < ntypes; ++t )
   {
      SCIP_PATTERN* pattern;
      SCIP_VAR* var;

      SCIP_CALL( SCIPpatternCreateCircular(scip, &pattern, ntypes, t) );
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_UNKNOWN);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "d%d", t);
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, (SCIP_Real)demands[t], 0.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, var) );

      /* add variable and pattern to the problem data */
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, pattern, var, SCIP_PATTERNTYPE_CIRCULAR) );

      SCIP_CALL( SCIPreleaseVar(scip, &var) );
      SCIPpatternRelease(scip, &pattern);
   }

   /* TODO compute all circular patterns */

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
   int* demands;
   int ntypes;
   int p;
   int t;

   assert(probdata != NULL);

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* tell SCIP that the objective will be always integral */
   SCIP_CALL( SCIPsetObjIntegral(scip) );

   ntypes = SCIPprobdataGetNTypes(probdata);
   demands = SCIPprobdataGetDemands(probdata);

   /* compute all non-dominated circular patterns */
   SCIP_CALL( computeCircularPatterns(scip, probdata) );

   /* create initial rectangular patterns */
   for( t = 0; t < ntypes; ++t )
   {
      SCIP_VAR* var;
      SCIP_PATTERN* pattern;

      /* create a pattern containing a single circle of type t */
      SCIP_CALL( SCIPpatternCreateRectangular(scip, &pattern, ntypes) );
      SCIPpatternAddElement(pattern, t);
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);

      /* create variable */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "r_%d", t);
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, var) );

      /* add pattern and variable to problem data */
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, pattern, var, SCIP_PATTERNTYPE_RECTANGULAR) );

      /* release pattern and variable */
      SCIPpatternRelease(scip, &pattern);
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

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
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->patternconss[t], name, 0, NULL, NULL, 0.0, SCIPinfinity(scip) ) );

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
         int nelems = SCIPpatternGetNElemens(pattern, t);

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
         int nelems = SCIPpatternGetNElemens(pattern, t);

         if( nelems > 0 )
         {
            /* + P_t z_P */
            SCIP_CALL( SCIPaddCoefLinear(scip, probdata->patternconss[t], var, (SCIP_Real)nelems) );
         }
      }
   }

   return SCIP_OKAY;
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
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->ntypes, (*targetdata)->patternconss, (*targetdata)->patternconss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->ncpatterns, (*targetdata)->cvars, (*targetdata)->cvars) );
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nrpatterns, (*targetdata)->rvars, (*targetdata)->rvars) );

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
   SCIP_Real*            rexts,              /**< external radii of each ring */
   int                   ntypes,             /**< number of different types */
   SCIP_Real             width,              /**< width of each rectangle */
   SCIP_Real             height              /**< height of each rectangle */
   )
{
   SCIP_PROBDATA* probdata;

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

   /* setup master problem*/
   SCIP_CALL( setupProblem(scip, probdata) );

   /* activate pricer */
   SCIP_CALL( SCIPpricerRingpackingActivate(scip) );

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

/** adds given variable to the problem data
 *
 * @note this function captures the variable and pattern
 */
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_PATTERN*         pattern,            /**< pattern */
   SCIP_VAR*             var,                /**< variables to add */
   SCIP_PATTERNTYPE      patterntype         /**< pattern type */
   )
{
   assert(probdata != NULL);
   assert(pattern != NULL);
   assert(var != NULL);
   assert(SCIPpatternGetPatternType(pattern) == patterntype);

   if( patterntype == SCIP_PATTERNTYPE_CIRCULAR )
   {
      SCIP_CALL( ensureSize(scip, probdata, patterntype, probdata->ncpatterns + 1) );
      probdata->cpatterns[probdata->ncpatterns] = pattern;
      probdata->cvars[probdata->ncpatterns] = var;
      ++(probdata->ncpatterns);
   }
   else
   {
      SCIP_CALL( ensureSize(scip, probdata, patterntype, probdata->nrpatterns + 1) );
      probdata->rpatterns[probdata->nrpatterns] = pattern;
      probdata->rvars[probdata->nrpatterns] = var;
      ++(probdata->nrpatterns);
   }

   /* capture variable and pattern */
   SCIP_CALL( SCIPcaptureVar(scip, var) );
   SCIPpatternCapture(pattern);

   return SCIP_OKAY;
}

/**@} */
