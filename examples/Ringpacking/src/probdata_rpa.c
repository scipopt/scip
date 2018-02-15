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

   probdata->enumtime = -SCIPgetTotalTime(scip);

   /* TODO compute all circular patterns (use soft working limits for the verification) */

   /* create an empty circular pattern for each type
    *
    * TODO remove this
    */
   for( t = 0; t < ntypes; ++t )
   {
      SCIP_PATTERN* pattern;
      SCIP_VAR* var;

      SCIP_CALL( SCIPpatternCreateCircular(scip, &pattern, t) );
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_UNKNOWN);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", t);
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, (SCIP_Real)demands[t], 0.0, SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip, var) );

      /* add variable and pattern to the problem data */
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, pattern, var, SCIP_PATTERNTYPE_CIRCULAR) );

      SCIP_CALL( SCIPreleaseVar(scip, &var) );
      SCIPpatternRelease(scip, &pattern);
   }

   /* update statistics */
   probdata->enumtime += SCIPgetTotalTime(scip);
   probdata->ncppatternsunknownbeg = getNCPatterns(scip, probdata, SCIP_PACKABLE_UNKNOWN);

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
      SCIP_CALL( SCIPpatternCreateRectangular(scip, &pattern) );
      SCIP_CALL( SCIPpatternAddElement(pattern, t) );
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

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputRpa)
{ /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   int* demands;
   int ntypes;
   int nrings;
   int t;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   ntypes = SCIPprobdataGetNTypes(probdata);
   demands = SCIPprobdataGetDemands(probdata);
   nrings = 0;

   /* count the number of rings */
   for( t = 0; t < ntypes; ++t )
      nrings += demands[t];

   SCIPinfoMessage(scip, file, "Ringpacking        : %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
      "ntypes", "nrings", "width", "height", "CP", "CP_unk", "CP_unk_end" ,"CP_no", "RP", "CP_time");

   SCIPinfoMessage(scip, file, "  %-17s:", SCIPgetProbName(scip));
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

   /* copy statistics to transformed problem data */
   (*targetdata)->ncppatternsunknownbeg = sourcedata->ncppatternsunknownbeg;
   (*targetdata)->enumtime = sourcedata->enumtime;

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

   /* add table output */
   assert(SCIPfindTable(scip, TABLE_NAME_RPA) == NULL);
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_RPA, TABLE_DESC_RPA, TRUE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputRpa,
         NULL, TABLE_POSITION_RPA, TABLE_EARLIEST_STAGE_RPA) );

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
   assert(SCIPpatternGetPackableStatus(pattern) != SCIP_PACKABLE_NO);

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
   SCIP_Real timestart;
   int niters;

   assert(probdata != NULL);
   assert(pattern != NULL);
   assert(iterlim > 0);
   assert(SCIPpatternGetPatternType(pattern) == SCIP_PATTERNTYPE_CIRCULAR);
   assert(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   rexts = SCIPprobdataGetRexts(probdata);
   rints = SCIPprobdataGetRints(probdata);

   /* pattern is empty -> set status to packable */
   if( SCIPpatternGetNElemens(pattern) == 0 )
   {
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);
      return SCIP_OKAY;
   }

   /* pattern contains only one element -> compare radii */
   if( SCIPpatternGetNElemens(pattern) == 1 )
   {
      int elemtype;
      int type;

      elemtype = SCIPpatternGetElementType(pattern, 0);
      assert(elemtype >= 0 && elemtype < SCIPprobdataGetNTypes(probdata));

      type = SCIPpatternGetType(pattern);
      assert(type >= 0 && type < SCIPprobdataGetNTypes(probdata));

      /* check whether element fits into the circular pattern */
      if( SCIPisGE(scip, rints[type], rexts[elemtype]) )
         SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);
      else
         SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_NO);

      return SCIP_OKAY;
   }

   timestart = SCIPgetTotalTime(scip);
   niters = 0;

   /* main loop for calling heuristic verification */
   while( SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN
      && niters < iterlim
      && SCIPgetTotalTime(scip) - timestart >= timelim )
   {
      /* TODO */

      ++niters;
   }

   return SCIP_OKAY;
}

/** verifies a circular pattern via a verification NLP */
extern
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
      assert(elemtype >= 0 && elemtype < nelems);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", k);
      SCIP_CALL( SCIPcreateVarBasic(subscip, &xvars[k], name, rexts[elemtype] - rints[type], rints[type] - rexts[elemtype], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(subscip, xvars[k]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d", k);
      SCIP_CALL( SCIPcreateVarBasic(subscip, &yvars[k], name, rexts[elemtype] - rints[type], rints[type] - rexts[elemtype], 1.0, SCIP_VARTYPE_CONTINUOUS) );
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

   SCIPdebugMsg(scip, "result of verification NLP: nsols=%d solstat=%d\n", SCIPgetNSols(subscip), SCIPgetStatus(subscip));

   /* check whether
    *
    *   -(at least) one solution could be found or
    *   - the problem is proven to be infeasible
    */
   if( SCIPgetNSols(subscip) > 0 )
      SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_YES);
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

/**@} */
