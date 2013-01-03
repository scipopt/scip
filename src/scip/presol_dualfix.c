/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
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

#include "scip/presol_dualfix.h"


#define PRESOL_NAME            "dualfix"
#define PRESOL_DESC            "roundable variables dual fixing"
#define PRESOL_PRIORITY         +8000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY              FALSE /**< should presolver be delayed, if other presolvers found reductions? */


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDualfix)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolDualfix(scip) );
 
   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDualfix)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real bound;
   SCIP_Real roundbound;
   SCIP_Real obj;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;
   int nvars;
   int v;

   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* get active problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* look for fixable variables
    * loop backwards, since a variable fixing can change the current and the subsequent slots in the vars array
    */
   for( v = nvars - 1; v >= 0; --v )
   {
      /* don't perform dual presolving operations on deleted variables */
      if( SCIPvarIsDeleted(vars[v]) )
         continue;

      obj = SCIPvarGetObj(vars[v]);

      /* if the objective coefficient of the variable is 0 and it may be rounded both
       * up and down, then fix it to the closest feasible value to 0 */
      if( SCIPisZero(scip, obj) && SCIPvarMayRoundDown(vars[v]) && SCIPvarMayRoundUp(vars[v]) )
      {
         bound = SCIPvarGetLbGlobal(vars[v]);
         if( SCIPisLT(scip, bound, 0.0) )
         {
            if( SCIPisLE(scip, 0.0, SCIPvarGetUbGlobal(vars[v])) )
               bound = 0.0;
            else
            {
               /* try to take an integer value, only for polishing */
               roundbound = SCIPfloor(scip, SCIPvarGetUbGlobal(vars[v]));
               
               if( roundbound < bound )
                  bound = SCIPvarGetUbGlobal(vars[v]);
               else
                  bound = roundbound;
            }
         }
         else
         {
            /* try to take an integer value, only for polishing */
            roundbound = SCIPceil(scip, bound);

            if( roundbound < SCIPvarGetUbGlobal(vars[v]) )
               bound = roundbound;
         }
         SCIPdebugMessage("variable <%s> with objective 0 fixed to %g\n",
            SCIPvarGetName(vars[v]), bound);
      }
      else
      {
         /* if it is always possible to round variable in direction of objective value,
          * fix it to its proper bound
          */
         if( SCIPvarMayRoundDown(vars[v]) && !SCIPisNegative(scip, obj) )
         {
            bound = SCIPvarGetLbGlobal(vars[v]);
            if( SCIPisZero(scip, obj) && SCIPvarGetNLocksUp(vars[v]) == 1 && SCIPisInfinity(scip, -bound) )
            {
               /* variable can be set to -infinity, and it is only contained in one constraint:
                * we hope that the corresponding constraint handler is clever enough to set/aggregate the variable
                * to something more useful than -infinity and do nothing here
                */
               continue;
            }
            SCIPdebugMessage("variable <%s> with objective %g and %d uplocks fixed to lower bound %g\n",
               SCIPvarGetName(vars[v]), SCIPvarGetObj(vars[v]), SCIPvarGetNLocksUp(vars[v]), bound);
         }
         else if( SCIPvarMayRoundUp(vars[v]) && !SCIPisPositive(scip, obj) )
         {
            bound = SCIPvarGetUbGlobal(vars[v]);
            if( SCIPisZero(scip, obj) && SCIPvarGetNLocksDown(vars[v]) == 1 && SCIPisInfinity(scip, bound) )
            {
               /* variable can be set to +infinity, and it is only contained in one constraint:
                * we hope that the corresponding constraint handler is clever enough to set/aggregate the variable
                * to something more useful than +infinity and do nothing here
                */
               continue;
            }
            SCIPdebugMessage("variable <%s> with objective %g and %d downlocks fixed to upper bound %g\n",
               SCIPvarGetName(vars[v]), SCIPvarGetObj(vars[v]), SCIPvarGetNLocksDown(vars[v]), bound);
         }
         else
            continue;
      }

      /* apply the fixing */
      if( SCIPisInfinity(scip, REALABS(bound)) && !SCIPisZero(scip, obj) )
      {
         SCIPdebugMessage(" -> unbounded fixing\n");
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
            "problem infeasible or unbounded: variable <%s> with objective %.15g can be made infinitely %s\n",
            SCIPvarGetName(vars[v]), SCIPvarGetObj(vars[v]), bound < 0.0 ? "small" : "large");
         *result = SCIP_UNBOUNDED;
         return SCIP_OKAY;
      }
      SCIP_CALL( SCIPfixVar(scip, vars[v], bound, &infeasible, &fixed) );
      if( infeasible )
      {
         SCIPdebugMessage(" -> infeasible fixing\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      assert(fixed);
      (*nfixedvars)++;
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the dual fixing presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDualfix(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presolptr;

   /* create dualfix presolver data */
   presoldata = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presolptr, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_DELAY,
         presolExecDualfix,
         presoldata) );

   assert(presolptr != NULL);

   SCIP_CALL( SCIPsetPresolCopy(scip, presolptr, presolCopyDualfix) );

   return SCIP_OKAY;
}
