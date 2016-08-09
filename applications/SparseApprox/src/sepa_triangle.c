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

/**@file   sepa_Triangle.c
 * @brief  Triangle separator
 * @author Leon Eifler
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "sepa_triangle.h"
#include "probdata_spa.h"
#include "scip/cons_linear.h"

#define SEPA_NAME              "triangle"
#define SEPA_DESC              "separator for graph partitioning in sparse approximation project"
#define SEPA_PRIORITY               -5000
#define SEPA_FREQ                   0
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */

static
SCIP_DECL_SEPACOPY(sepaCopyTriangle)
{   /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaTriangle(scip) );

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpTriangle)
{
   SCIP_VAR***** edgevars;
   int i;               /* the bin in the first partition */
   int j;
   int h;
   int nbins;
   int ncuts = 0;
   SCIP_ROW* cut;
   SCIP_Bool infeasible;
   char cutname[SCIP_MAXSTRLEN];

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   nbins = SCIPspaGetNrBins(scip);
   assert(nbins > 0);

   /* get the variables from scip */

   edgevars = SCIPspaGetEdgevars(scip);

   /* check and add violated triangle inequalities */
   for( i = 0; i < nbins; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         for( h = 0; h < j; ++h )
         {
            SCIP_Real triangleval1 = 0;
            SCIP_Real triangleval2= 0;
            SCIP_Real triangleval3 = 0;
            {
               if( NULL == edgevars[i][j][0][0] || NULL == edgevars[j][h][0][0] || NULL == edgevars[i][h][0][0] )
                  break;
               triangleval1 = SCIPvarGetLPSol(edgevars[i][j][0][0]) + SCIPvarGetLPSol(edgevars[j][h][0][0]) - SCIPvarGetLPSol(edgevars[i][h][0][0]);
               triangleval2 = SCIPvarGetLPSol(edgevars[i][j][0][0]) - SCIPvarGetLPSol(edgevars[j][h][0][0]) + SCIPvarGetLPSol(edgevars[i][h][0][0]);
               triangleval3 = - SCIPvarGetLPSol(edgevars[i][j][0][0]) + SCIPvarGetLPSol(edgevars[j][h][0][0]) + SCIPvarGetLPSol(edgevars[i][h][0][0]);

               if( SCIPisGT(scip, triangleval1, 1.0) )
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "triangle1_%d_%d_%d", i+1,j+1,h+1 );
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIPaddVarToRow(scip, cut, edgevars[i][j][0][0], 1.0);
                  SCIPaddVarToRow(scip, cut, edgevars[j][h][0][0], 1.0);
                  SCIPaddVarToRow(scip, cut, edgevars[i][h][0][0], -1.0);

                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, NULL, cut) )
                  {
                     SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
                     if( !infeasible )
                     {
                        SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                        *result = SCIP_SEPARATED;
                     }
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                  ncuts++;
               }
               if( SCIPisGT(scip, triangleval2, 1.0) )
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "triangle2_%d_%d_%d", i+1,j+1,h+1 );
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIPaddVarToRow(scip, cut, edgevars[i][j][0][0], 1.0);
                  SCIPaddVarToRow(scip, cut, edgevars[j][h][0][0], -1.0);
                  SCIPaddVarToRow(scip, cut, edgevars[i][h][0][0], +1.0);

                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, NULL, cut) )
                  {
                     SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
                     if( !infeasible )
                     {
                        SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                        *result = SCIP_SEPARATED;
                     }
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                  ncuts++;
               }
               if( SCIPisGT(scip, triangleval3, 1.0) )
               {
                  (void)SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "triangle3_%d_%d_%d", i+1,j+1,h+1 );
                  SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
                  SCIPaddVarToRow(scip, cut, edgevars[i][j][0][0], -1.0);
                  SCIPaddVarToRow(scip, cut, edgevars[j][h][0][0], +1.0);
                  SCIPaddVarToRow(scip, cut, edgevars[i][h][0][0], +1.0);

                  SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
                  if( SCIPisCutEfficacious(scip, NULL, cut) )
                  {
                     SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
                     if( !infeasible )
                     {
                        SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                        *result = SCIP_SEPARATED;
                     }
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                  ncuts++;
               }
            }
         }
      }
   }
   return SCIP_OKAY;
}

/** creates the Triangle separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaTriangle(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_SEPA* sepa;


   /* include separator */

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP, SEPA_DELAY,
      sepaExeclpTriangle, NULL,
      NULL) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyTriangle) );

   return SCIP_OKAY;
}
