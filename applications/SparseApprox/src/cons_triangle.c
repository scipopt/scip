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

/**@file   cons_triangle.c
 * @brief  constraint handler for triangle constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cons_triangle.h"
#include "probdata_spa.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "triangle"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_ENFOPRIORITY       -1000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY      -1000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          -1  /**< frequency for using all instead of only the useful constraints in separation,
 *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY        1000000 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            10 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */


/*
 * Data structures
 */

/*
 * Local methods
 */

static
SCIP_RETCODE sepaTriangle(
   SCIP*                 scip,
   SCIP_CONSHDLR*        conshdlr,
   SCIP_SOL*             sol,
   SCIP_RESULT*          result
)
{
   SCIP_VAR**** edgevars;
   SCIP_ROW* cut;
   int nbins;
   int i;
   int j;
   int k;
   int l;
   char cutname[SCIP_MAXSTRLEN];
   int** sign;
   int maxcuts;
   int ncuts;
   SCIP_Real violation;

   edgevars = SCIPspaGetEdgevars(scip);
   nbins = SCIPspaGetNrBins(scip);

   assert(nbins > 0);
   assert(NULL != edgevars);

   *result = SCIP_DIDNOTFIND;
   maxcuts = 500;
   ncuts = 0;
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

   for( i = 0; i < nbins && ncuts < maxcuts; ++i )
   {
      for( j = 0; j < nbins && ncuts < maxcuts; ++j )
      {
         for( k = 0; k < j && ncuts < maxcuts; ++k )
         {
            if( NULL == edgevars[i][j] || NULL == edgevars[i][k] || NULL == edgevars[j][k] )
               continue;
            if( (i != j && i != k && j > k) )
            {

               violation = SCIPgetSolVal(scip, sol, edgevars[i][j][1]) + SCIPgetSolVal(scip, sol, edgevars[i][k][1]) - SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][k][0], -1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][j][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][k][1], 1.0) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, sol, cut) )
                  {
                     SCIP_Bool infeasible;
                     ncuts++;
                     SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE, &infeasible) );
                     if ( infeasible )
                        *result = SCIP_CUTOFF;
                     else
                        *result = SCIP_SEPARATED;
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }

               violation = SCIPgetSolVal(scip, sol, edgevars[j][i][1]) + SCIPgetSolVal(scip, sol, edgevars[k][i][1]) - SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &(cut), conshdlr, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][k][0], -1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][i][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[k][i][1], 1.0) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, sol, cut) )
                  {
                     SCIP_Bool infeasible;
                     ncuts++;
                     SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE, &infeasible) );
                     if ( infeasible )
                        *result = SCIP_CUTOFF;
                     else
                        *result = SCIP_SEPARATED;
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }

               violation = SCIPgetSolVal(scip, sol, edgevars[i][j][1]) - SCIPgetSolVal(scip, sol, edgevars[i][k][1]) + SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &(cut), conshdlr, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][k][0], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][j][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][k][1], -1.0) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, sol, cut) )
                  {
                     SCIP_Bool infeasible;
                     ncuts++;
                     SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE, &infeasible) );
                     if ( infeasible )
                        *result = SCIP_CUTOFF;
                     else
                        *result = SCIP_SEPARATED;
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }

               violation = SCIPgetSolVal(scip, sol, edgevars[j][i][1]) - SCIPgetSolVal(scip, sol, edgevars[k][i][1]) + SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &(cut), conshdlr, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][k][0], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][i][1], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[k][i][1], -1.0) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, sol, cut) )
                  {
                     SCIP_Bool infeasible;
                     ncuts++;
                     SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE, &infeasible) );
                     if ( infeasible )
                        *result = SCIP_CUTOFF;
                     else
                        *result = SCIP_SEPARATED;
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }

               violation = - SCIPgetSolVal(scip, sol, edgevars[i][j][1]) + SCIPgetSolVal(scip, sol, edgevars[i][k][1]) + SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &(cut), conshdlr, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][k][0], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][j][1], -1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][k][1], 1.0) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, sol, cut) )
                  {
                     SCIP_Bool infeasible;
                     ncuts++;
                     SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE, &infeasible) );
                     if ( infeasible )
                        *result = SCIP_CUTOFF;
                     else
                        *result = SCIP_SEPARATED;
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }

               violation = -SCIPgetSolVal(scip, sol, edgevars[j][i][1]) + SCIPgetSolVal(scip, sol, edgevars[k][i][1]) + SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &(cut), conshdlr, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][k][0], 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][i][1], -1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[k][i][1], 1.0) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, sol, cut) )
                  {
                     SCIP_Bool infeasible;
                     ncuts++;
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

   /* separate edges by the valid inequality y_ij0 + y_ik0 - y_jk0 <= 1 */
   for( i = 0; i < nbins && ncuts < maxcuts; ++i )
   {
      for( j = 0; j < i && ncuts < maxcuts; ++j )
      {
         for( k = 0; k < j && ncuts < maxcuts; ++k )
         {
            if( NULL == edgevars[i][j] || NULL == edgevars[j][k] || NULL == edgevars[i][k] )
               continue;
            for( l = 0; l < 3; ++l )
            {
               violation = sign[l][0] * SCIPgetSolVal(scip, sol, edgevars[i][j][0]) + sign[l][1] * SCIPgetSolVal(scip, sol, edgevars[i][k][0]) + sign[l][2] * SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation) )
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "edgecut_%d_%d_%d", i+1, j+1, k+1 );
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &(cut), conshdlr, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][j][0], sign[l][0]) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[i][k][0], sign[l][1]) );
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, edgevars[j][k][0], sign[l][2]) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, sol, cut) )
                  {
                     SCIP_Bool infeasible;
                     ncuts++;
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
   int nbins;
   int i;
   int j;
   int k;
   int l;
   int sign[3][3];
   SCIP_Real violation;

   edgevars = SCIPspaGetEdgevars(scip);
   nbins = SCIPspaGetNrBins(scip);

   assert(nbins > 0);
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
      for( j = 0; j < nbins; ++j )
      {
         for( k = 0; k < j; ++k )
         {
            if( NULL == edgevars[i][j] || NULL == edgevars[i][k] || NULL == edgevars[j][k] )
               continue;
            if( (i != j && i != k && j > k) )
            {

               violation = SCIPgetSolVal(scip, sol, edgevars[i][j][1]) + SCIPgetSolVal(scip, sol, edgevars[i][k][1]) - SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
                  return TRUE;

               violation = SCIPgetSolVal(scip, sol, edgevars[j][i][1]) + SCIPgetSolVal(scip, sol, edgevars[k][i][1]) - SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
                  return TRUE;

               violation = SCIPgetSolVal(scip, sol, edgevars[i][j][1]) - SCIPgetSolVal(scip, sol, edgevars[i][k][1]) + SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
                  return TRUE;

               violation = SCIPgetSolVal(scip, sol, edgevars[j][i][1]) - SCIPgetSolVal(scip, sol, edgevars[k][i][1]) + SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
                  return TRUE;

               violation = - SCIPgetSolVal(scip, sol, edgevars[i][j][1]) + SCIPgetSolVal(scip, sol, edgevars[i][k][1]) + SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
                  return TRUE;

               violation = -SCIPgetSolVal(scip, sol, edgevars[j][i][1]) + SCIPgetSolVal(scip, sol, edgevars[k][i][1]) + SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation))
                  return TRUE;
            }
         }
      }
   }

   /* separate edges by the valid inequality y_ij0 + y_ik0 - y_jk0 <= 1 */
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         for( k = 0; k < j; ++k )
         {
            if( NULL == edgevars[i][j] || NULL == edgevars[j][k] || NULL == edgevars[i][k] )
               continue;
            for( l = 0; l < 3; ++l )
            {
               violation = sign[l][0] * SCIPgetSolVal(scip, sol, edgevars[i][j][0]) + sign[l][1] * SCIPgetSolVal(scip, sol, edgevars[i][k][0]) + sign[l][2] * SCIPgetSolVal(scip, sol, edgevars[j][k][0]) - 1;
               if( SCIPisPositive(scip, violation) )
                  return TRUE;
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
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyTriangle NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_CONSFREE(consFreeTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeTriangle NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitTriangle NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitTriangle NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreTriangle NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreTriangle NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolTriangle NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolTriangle NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteTriangle NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransTriangle NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpTriangle NULL
#endif


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpTriangle)
{
   SCIP_CALL( sepaTriangle(scip, conshdlr, NULL, result) );
   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */

static
SCIP_DECL_CONSSEPASOL(consSepasolTriangle)
{
   SCIP_CALL( sepaTriangle(scip, conshdlr, sol, result) );
   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpTriangle)
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
SCIP_DECL_CONSENFOPS(consEnfopsTriangle)
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
SCIP_DECL_CONSCHECK(consCheckTriangle)
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
SCIP_DECL_CONSPROP(consPropTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropTriangle NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolTriangle NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropTriangle NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockTriangle)
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
SCIP_DECL_CONSACTIVE(consActiveTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveTriangle NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveTriangle NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableTriangle NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableTriangle NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsTriangle NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintTriangle NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyTriangle NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseTriangle NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsTriangle NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsTriangle NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsTriangle)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of triangle constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsTriangle NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for triangle constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrTriangle(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create triangle constraint handler data */
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
      conshdlrCopyTriangle,
      consFreeTriangle, consInitTriangle, consExitTriangle,
      consInitpreTriangle, consExitpreTriangle, consInitsolTriangle, consExitsolTriangle,
      consDeleteTriangle, consTransTriangle, consInitlpTriangle,
      consSepalpTriangle, consSepasolTriangle, consEnfolpTriangle, consEnfopsTriangle, consCheckTriangle,
      consPropTriangle, consPresolTriangle, consRespropTriangle, consLockTriangle,
      consActiveTriangle, consDeactiveTriangle,
      consEnableTriangle, consDisableTriangle, consDelvarsTriangle,
      consPrintTriangle, consCopyTriangle, consParseTriangle,
      consGetVarsTriangle, consGetNVarsTriangle, consGetDiveBdChgsTriangle, conshdlrdata) );
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
      CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
      consEnfolpTriangle, consEnfopsTriangle, consCheckTriangle, consLockTriangle,
      conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveTriangle) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyTriangle, consCopyTriangle) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveTriangle) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteTriangle) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsTriangle) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableTriangle) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableTriangle) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitTriangle) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreTriangle) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolTriangle) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeTriangle) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsTriangle) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsTriangle) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsTriangle) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitTriangle) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreTriangle) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolTriangle) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpTriangle) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseTriangle) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolTriangle, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintTriangle) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropTriangle, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
      CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropTriangle) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpTriangle, consSepasolTriangle, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransTriangle) );
#endif

#ifdef LINCONSUPGD_PRIORITY
   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdTriangle, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }
#endif

   /* add triangle constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a triangle constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsTriangle(
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsTriangle() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the triangle constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("triangle constraint handler not found\n");
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
