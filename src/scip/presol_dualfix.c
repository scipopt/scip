/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_dualfix.c
 * @brief  fixing roundable variables to best bound
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "presol_dualfix.h"


#define PRESOL_NAME              "dualfix"
#define PRESOL_DESC              "roundable variables dual fixing"
#define PRESOL_PRIORITY          100000



/*
 * Callback methods
 */

static
DECL_PRESOLEXEC(presolExecDualfix)
{
   VAR** vars;
   Real bound;
   Bool infeasible;
   int nvars;
   int v;

   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);
   assert(result != NULL);
   
   /* get active problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* look for fixable variables
    * loop backwards, since a variable fixing can change the actual and the subsequent slots in the vars array
    */
   for( v = nvars-1; v >= 0; --v )
   {
      /* if it is always possible to round variable in direction of objective value,
       * fix it to its proper bound
       */
      if( SCIPvarMayRoundDown(vars[v]) && !SCIPisNegative(scip, SCIPvarGetObj(vars[v])) )
      {
         bound = SCIPvarGetLbGlobal(vars[v]);
         debugMessage("variable <%s> with objective %g fixed to lower bound %g\n",
            SCIPvarGetName(vars[v]), SCIPvarGetObj(vars[v]), bound);
      }
      else if( SCIPvarMayRoundUp(vars[v]) && !SCIPisPositive(scip, SCIPvarGetObj(vars[v])) )
      {
         bound = SCIPvarGetUbGlobal(vars[v]);
         debugMessage("variable <%s> with objective %g fixed to upper bound %g\n",
            SCIPvarGetName(vars[v]), SCIPvarGetObj(vars[v]), bound);
      }
      else
         continue;
      
      /* apply the fixing */
      if( SCIPisInfinity(scip, ABS(bound)) )
      {
         debugMessage(" -> unbounded fixing\n");
         *result = SCIP_UNBOUNDED;
         return SCIP_OKAY;
      }
      CHECK_OKAY( SCIPfixVar(scip, vars[v], bound, &infeasible) );
      if( infeasible )
      {
         debugMessage(" -> infeasible fixing\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      (*nfixedvars)++;
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}





/*
 * presolver specific interface methods
 */

/** creates the most infeasible LP braching rule and includes it in SCIP */
RETCODE SCIPincludePresolDualfix(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY,
                  NULL, NULL, NULL, presolExecDualfix,
                  NULL) );

   return SCIP_OKAY;
}
