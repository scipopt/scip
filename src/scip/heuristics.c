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

/**@file   heuristics.c
 * @brief  methods commonly used by primal heuristics
 * @author Gregor Hendel
 */

#include "scip/heuristics.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"

/** creates the rows of the subproblem */
static
SCIP_RETCODE createRows(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_HASHMAP*         varmap              /**< a hashmap to store the mapping of source variables to the corresponding
                                               *   target variables */
   )
{
   SCIP_ROW** rows;                          /* original scip rows                       */
   SCIP_CONS* cons;                          /* new constraint                           */
   SCIP_VAR** consvars;                      /* new constraint's variables               */
   SCIP_COL** cols;                          /* original row's columns                   */

   SCIP_Real constant;                       /* constant added to the row                */
   SCIP_Real lhs;                            /* left hand side of the row                */
   SCIP_Real rhs;                            /* left right side of the row               */
   SCIP_Real* vals;                          /* variables' coefficient values of the row */

   int nrows;
   int nnonz;
   int i;
   int j;

   /* get the rows and their number */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* copy all rows to linear constraints */
   for( i = 0; i < nrows; i++ )
   {
      /* ignore rows that are only locally valid */
      if( SCIProwIsLocal(rows[i]) )
         continue;

      /* get the row's data */
      constant = SCIProwGetConstant(rows[i]);
      lhs = SCIProwGetLhs(rows[i]) - constant;
      rhs = SCIProwGetRhs(rows[i]) - constant;
      vals = SCIProwGetVals(rows[i]);
      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);

      assert(lhs <= rhs);

      /* allocate memory array to be filled with the corresponding subproblem variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ )
         consvars[j] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, (SCIPcolGetVar(cols[j])));

      /* create a new linear constraint and add it to the subproblem */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}


/** get a sub-SCIP copy of the transformed problem */
SCIP_RETCODE SCIPheuristicsInstanceCopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP used by the heuristic */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables */
   const char*           suffix,             /**< suffix for the problem name */
   SCIP_VAR**            fixedvars,          /**< source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Real*            fixedvals,          /**< array of fixing values for target SCIP variables, or NULL */
   int                   nfixedvars,         /**< number of source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Bool             uselprows,          /**< should the linear relaxation of the problem defined by LP rows be copied? */
   SCIP_Bool             copycuts,           /**< should cuts be copied (only if uselprows == FALSE) */
   SCIP_Bool*            success             /**< was the copying successful? */
   )
{
   assert(sourcescip != NULL);
   assert(suffix != NULL);
   assert(subscip != NULL);
   assert(varmap != NULL);

   if( uselprows )
   {
      char probname[SCIP_MAXSTRLEN];

      /* copy all plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

      /* get name of the original problem and add the string "_crossoversub" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_%s", SCIPgetProbName(sourcescip), suffix);

      /* create the subproblem */
      SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* copy all variables */
      SCIP_CALL( SCIPcopyVars(sourcescip, subscip, varmap, NULL,  fixedvars, fixedvals, nfixedvars, TRUE) );

      /* create linear constraints from LP rows of the source problem */
      SCIP_CALL( createRows(sourcescip, subscip, varmap) );
   }
   else
   {
      SCIP_CALL( SCIPcopy(sourcescip, subscip, varmap, NULL, suffix, fixedvars, fixedvals, nfixedvars, TRUE, FALSE, TRUE, success) );

      if( copycuts )
      {
         /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(sourcescip, subscip, varmap, NULL, TRUE, NULL) );
      }
   }

   return SCIP_OKAY;
}
