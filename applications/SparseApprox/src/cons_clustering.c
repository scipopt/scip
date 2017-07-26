/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_clustering.c
 * @brief  constraint handler for clustering constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cons_clustering.h"
#include "probdata_spa.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "clustering"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_ENFOPRIORITY       -2000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY      -2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          1  /**< frequency for using all instead of only the useful constraints in separation,
 *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY        1000000 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */


/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for clustering constraints */
struct SCIP_ConsData
{
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
};


/*
 * Local methods
 */

static
SCIP_RETCODE sepaClustering(
   SCIP*                 scip,
   SCIP_CONSHDLR*        conshdlr,
   SCIP_SOL*             sol,
   SCIP_RESULT*          result
)
{
   SCIP_VAR**** edgevars;
   SCIP_VAR*** binvars;
   SCIP_ROW* cut;
   SCIP_VAR* var;
   int nbins;
   int ncluster;
   int i;
   int j;
   int c1;
   int c2;
   int l;
   char cutname[SCIP_MAXSTRLEN];
   int** sign;
   SCIP_Real violation;

   edgevars = SCIPspaGetEdgevars(scip);
   binvars = SCIPspaGetBinvars(scip);
   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);

   assert(nbins > 0);
   assert(ncluster > 0 && ncluster < nbins);
   assert(NULL != edgevars);

   *result = SCIP_DIDNOTFIND;
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &sign, 3) );
   for( i = 0; i < 3; ++i )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &sign[i], 3) );
      for( j = 0; j < 3; ++j )
      {
         sign[i][j] = 1;
      }
      sign[i][i] = -1;
   }

   for( i = 0; i < nbins && *result != SCIP_CUTOFF; ++i )
   {
      for( j = 0; j < i && *result != SCIP_CUTOFF; ++j )
      {
         if( NULL == edgevars[i][j] )
            continue;
         for( c1 = 0; c1 < ncluster && *result != SCIP_CUTOFF; ++c1 )
         {
            for( l = 0; l < 1 && *result != SCIP_CUTOFF; ++l )
            {
               violation = sign[l][0] * SCIPgetSolVal(scip, sol, edgevars[i][j][0]) + sign[l][1] * SCIPgetSolVal(scip, sol, binvars[i][c1]) + sign[l][2] * SCIPgetSolVal(scip, sol, binvars[j][c1]) - 1;
               if( SCIPisPositive(scip, violation) )
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d_var%d", i, j, c1, l);
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, cutname, -SCIPinfinity(scip), 1, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][j][0], sign[l][0]) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, binvars[i][c1], sign[l][1]) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, binvars[j][c1], sign[l][2]) );

                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  // add cut
                  if( SCIPisCutEfficacious(scip, sol, cut) )
                  {
                     SCIP_Bool infeasible;
                     SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE, &infeasible) );
                     if ( infeasible )
                        *result = SCIP_CUTOFF;
                     else
                        *result = SCIP_SEPARATED;
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }
               for( c2 = 0; c2 < ncluster && *result != SCIP_CUTOFF; ++c2 )
               {
                  if( c2 == c1 + 1 || ( c2 == 0 && c1 == ncluster -1) )
                     var = edgevars[i][j][1];
                  else if( c2 == c1 - 1 || ( c1 == 0 && c2 == ncluster -1) )
                     var = edgevars[j][i][1];
                  else
                     continue;

                  violation = sign[l][0] * SCIPgetSolVal(scip, sol, var) + sign[l][1] * SCIPgetSolVal(scip, sol, binvars[i][c1]) + sign[l][2] * SCIPgetSolVal(scip, sol, binvars[j][c2]) - 1;

                  if( SCIPisPositive(scip, violation) )
                  {
                     (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "bins_%d_%d_incluster_%d_%d_var%d", i, j, c1, c2, l);
                     SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, cutname, -SCIPinfinity(scip), 1, FALSE, FALSE, TRUE) );
                     SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cut, var, sign[l][0]) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cut, binvars[i][c1], sign[l][1]) );
                     SCIP_CALL( SCIPaddVarToRow(scip, cut, binvars[j][c2], sign[l][2]) );

                     SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                     // add cut
                     if( SCIPisCutEfficacious(scip, sol, cut) )
                     {
                        SCIP_Bool infeasible;
                        SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE, &infeasible) );
                        if ( infeasible )
                           *result = SCIP_CUTOFF;
                        else
                           *result = SCIP_SEPARATED;
                     }
                     SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                  }
               }
            }
         }
      }
   }

   for( i = 0; i < 3; ++i )
   {
      SCIPfreeMemoryArray(scip, &sign[i]);
   }
   SCIPfreeMemoryArray(scip, &sign);
   return SCIP_OKAY;
}


static
SCIP_Bool findViolation(
   SCIP*                 scip,
   SCIP_SOL*             sol
)
{
   SCIP_VAR**** edgevars;
   SCIP_VAR*** binvars;
   SCIP_VAR* var;
   int nbins;
   int ncluster;
   int i;
   int j;
   int c1;
   int c2;
   int l;
   int sign[3][3];
   SCIP_Real violation;

   edgevars = SCIPspaGetEdgevars(scip);
   binvars = SCIPspaGetBinvars(scip);
   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);

   assert(nbins > 0);
   assert(ncluster > 0 && ncluster < nbins);
   assert(NULL != edgevars);

   for( i = 0; i < 3; ++i )
   {
      for( j = 0; j < 3; ++j )
      {
         sign[i][j] = 1;
      }
      sign[i][i] = -1;
   }

   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( NULL == edgevars[i][j] )
            continue;
         for( c1 = 0; c1 < ncluster; ++c1 )
         {
            for( l = 0; l < 1; ++l )
            {
               violation = sign[l][0] * SCIPgetSolVal(scip, sol, edgevars[i][j][0]) + sign[l][1] * SCIPgetSolVal(scip, sol, binvars[i][c1]) + sign[l][2] * SCIPgetSolVal(scip, sol, binvars[j][c1]) - 1;
               if( SCIPisPositive(scip, violation) )
                  return TRUE;
               for( c2 = 0; c2 < ncluster; ++c2 )
               {
                  if( c2 == c1 + 1 || ( c2 == 0 && c1 == ncluster -1) )
                     var = edgevars[i][j][1];
                  else if( c2 == c1 - 1 || ( c1 == 0 && c2 == ncluster -1) )
                     var = edgevars[j][i][1];
                  else
                     continue;

                  violation = sign[l][0] * SCIPgetSolVal(scip, sol, var) + sign[l][1] * SCIPgetSolVal(scip, sol, binvars[i][c1]) + sign[l][2] * SCIPgetSolVal(scip, sol, binvars[j][c2]) - 1;

                  if( SCIPisPositive(scip, violation) )
                     return TRUE;
               }
            }
         }
      }
   }

   return FALSE;
}

/* put your local methods here, and declare them static */


/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyClustering NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_CONSFREE(consFreeClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeClustering NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitClustering NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitClustering NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreClustering NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreClustering NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolClustering NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolClustering NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteClustering NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransClustering NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpClustering NULL
#endif


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpClustering)
{
   SCIP_CALL( sepaClustering(scip, conshdlr, NULL, result) );
   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */

static
SCIP_DECL_CONSSEPASOL(consSepasolClustering)
{
   SCIP_CALL( sepaClustering(scip, conshdlr, sol, result) );
   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpClustering)
{
   SCIP_Bool found;

   *result = SCIP_FEASIBLE;
   found = findViolation(scip, NULL);

   if( found )
      *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsClustering)
{
   SCIP_Bool found;

   *result = SCIP_FEASIBLE;
   found = findViolation(scip, NULL);

   if( found )
      *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckClustering)
{
   SCIP_Bool found;

   *result = SCIP_FEASIBLE;
   found = findViolation(scip, sol);

   if( found )
      *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropClustering NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolClustering NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropClustering NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockClustering)
{
   int nbins;
   int ncluster;
   SCIP_VAR**** edgevars;
   SCIP_VAR*** binvars;
   int i;
   int j;
   int c;

   nbins = SCIPspaGetNrBins(scip);
   ncluster = SCIPspaGetNrCluster(scip);
   edgevars = SCIPspaGetEdgevars(scip);
   binvars = SCIPspaGetBinvars(scip);
   for( i = 0; i < nbins; ++i )
   {
      for( c = 0; c < ncluster; ++c )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, binvars[i][c], nlockspos + nlocksneg, nlockspos + nlocksneg) );
      }
   }
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         if( NULL != edgevars[i][j] )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, edgevars[i][j][0], nlockspos + nlocksneg, nlockspos + nlocksneg) );
            SCIP_CALL( SCIPaddVarLocks(scip, edgevars[i][j][1], nlockspos + nlocksneg, nlockspos + nlocksneg) );
            SCIP_CALL( SCIPaddVarLocks(scip, edgevars[j][i][1], nlockspos + nlocksneg, nlockspos + nlocksneg) );
         }
      }
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveClustering NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveClustering NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableClustering NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableClustering NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsClustering NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintClustering NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyClustering NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseClustering NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsClustering NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsClustering NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsClustering)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of clustering constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsClustering NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for clustering constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrClustering(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create clustering constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   conshdlr = NULL;

   /* include constraint handler */
#if 0
   /* use SCIPincludeConshdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
      CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
      CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
      CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
      CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
      conshdlrCopyClustering,
      consFreeClustering, consInitClustering, consExitClustering,
      consInitpreClustering, consExitpreClustering, consInitsolClustering, consExitsolClustering,
      consDeleteClustering, consTransClustering, consInitlpClustering,
      consSepalpClustering, consSepasolClustering, consEnfolpClustering, consEnfopsClustering, consCheckClustering,
      consPropClustering, consPresolClustering, consRespropClustering, consLockClustering,
      consActiveClustering, consDeactiveClustering,
      consEnableClustering, consDisableClustering, consDelvarsClustering,
      consPrintClustering, consCopyClustering, consParseClustering,
      consGetVarsClustering, consGetNVarsClustering, consGetDiveBdChgsClustering, conshdlrdata) );
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
      CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
      consEnfolpClustering, consEnfopsClustering, consCheckClustering, consLockClustering,
      conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveClustering) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyClustering, consCopyClustering) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveClustering) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteClustering) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsClustering) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableClustering) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableClustering) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitClustering) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreClustering) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolClustering) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeClustering) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsClustering) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsClustering) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsClustering) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitClustering) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreClustering) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolClustering) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpClustering) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseClustering) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolClustering, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintClustering) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropClustering, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
      CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropClustering) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpClustering, consSepasolClustering, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransClustering) );
#endif

#ifdef LINCONSUPGD_PRIORITY
   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdClustering, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }
#endif

   /* add clustering constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a clustering constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsClustering(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
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
    *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
    *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
    *   if it may be moved to a more global node?
    *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
)
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsClustering() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the clustering constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("clustering constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   /* TODO: create and store constraint specific data here */

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
      local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}
