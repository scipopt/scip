/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: presol_dualfix.c,v 1.22 2005/07/15 17:20:12 bzfpfend Exp $"

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
#define PRESOL_PRIORITY         +100000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY              FALSE /**< should presolver be delayed, if other presolvers found reductions? */




/*
 * Callback methods of presolver
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
#define presolFreeDualfix NULL


/** initialization method of presolver (called after problem was transformed) */
#define presolInitDualfix NULL


/** deinitialization method of presolver (called before transformed problem is freed) */
#define presolExitDualfix NULL


/** presolving initialization method of presolver (called when presolving is about to begin) */
#define presolInitpreDualfix NULL


/** presolving deinitialization method of presolver (called after presolving has been finished) */
#define presolExitpreDualfix NULL


/** execution method of presolver */
static
DECL_PRESOLEXEC(presolExecDualfix)
{  /*lint --e{715}*/
   VAR** vars;
   Real bound;
   Real obj;
   Bool infeasible;
   Bool fixed;
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
   for( v = nvars-1; v >= 0; --v )
   {
      /* don't perform dual presolving operations on deleted variables */
      if( SCIPvarIsDeleted(vars[v]) )
         continue;

      obj = SCIPvarGetObj(vars[v]);

      /* if it is always possible to round variable in direction of objective value,
       * fix it to its proper bound
       */
      if( SCIPvarMayRoundDown(vars[v]) && !SCIPisNegative(scip, obj) )
      {
         bound = SCIPvarGetLbGlobal(vars[v]);
         debugMessage("variable <%s> with objective %g fixed to lower bound %g\n",
            SCIPvarGetName(vars[v]), SCIPvarGetObj(vars[v]), bound);
      }
      else if( SCIPvarMayRoundUp(vars[v]) && !SCIPisPositive(scip, obj) )
      {
         bound = SCIPvarGetUbGlobal(vars[v]);
         debugMessage("variable <%s> with objective %g fixed to upper bound %g\n",
            SCIPvarGetName(vars[v]), SCIPvarGetObj(vars[v]), bound);
      }
      else
         continue;
      
      /* apply the fixing */
      if( SCIPisInfinity(scip, REALABS(bound)) && !SCIPisZero(scip, obj) )
      {
         debugMessage(" -> unbounded fixing\n");
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
            "problem infeasible or unbounded: variable <%s> with objective %g can be made infinitely %s\n",
            SCIPvarGetName(vars[v]), SCIPvarGetObj(vars[v]), bound < 0.0 ? "small" : "large");
         *result = SCIP_UNBOUNDED;
         return SCIP_OKAY;
      }
      CHECK_OKAY( SCIPfixVar(scip, vars[v], bound, &infeasible, &fixed) );
      if( infeasible )
      {
         debugMessage(" -> infeasible fixing\n");
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
RETCODE SCIPincludePresolDualfix(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   PRESOLDATA* presoldata;

   /* create dualfix presolver data */
   presoldata = NULL;

   /* include presolver */
   CHECK_OKAY( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_DELAY,
         presolFreeDualfix, presolInitDualfix, presolExitDualfix, 
         presolInitpreDualfix, presolExitpreDualfix, presolExecDualfix,
         presoldata) );

   return SCIP_OKAY;
}
