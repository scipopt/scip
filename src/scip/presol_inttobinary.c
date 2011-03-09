/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_inttobinary.c
 * @ingroup PRESOLVERS
 * @brief  presolver that converts integer variables with domain [a,a+1] to binaries
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/presol_inttobinary.h"
#include "scip/pub_misc.h"


#define PRESOL_NAME            "inttobinary"
#define PRESOL_DESC            "converts integer variables with domain [a,a+1] to binaries"
#define PRESOL_PRIORITY        +7000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY              FALSE /**< should presolver be delayed, if other presolvers found reductions? */




/*
 * Callback methods of presolver
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
#define presolFreeInttobinary NULL


/** initialization method of presolver (called after problem was transformed) */
#define presolInitInttobinary NULL


/** deinitialization method of presolver (called before transformed problem is freed) */
#define presolExitInttobinary NULL


/** presolving initialization method of presolver (called when presolving is about to begin) */
#define presolInitpreInttobinary NULL


/** presolving deinitialization method of presolver (called after presolving has been finished) */
#define presolExitpreInttobinary NULL


/** presolving execution method */
static
SCIP_DECL_PRESOLEXEC(presolExecInttobinary)
{  /*lint --e{715}*/
   SCIP_VAR** scipvars;
   SCIP_VAR** vars;
   int nbinvars;
   int nintvars;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get the problem variables */
   scipvars = SCIPgetVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   nintvars = SCIPgetNIntVars(scip);
   if( nintvars == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* copy the integer variables into an own array, since adding binary variables affects the left-most slots in the
    * array and thereby interferes with our search loop
    */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, &scipvars[nbinvars], nintvars) );

   /* scan the integer variables for possible conversion into binaries;
    * we have to collect the variables first in an own 
    */
   for( v = 0; v < nintvars; ++v )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER);

      /* get variable's bounds */
      lb = SCIPvarGetLbGlobal(vars[v]);
      ub = SCIPvarGetUbGlobal(vars[v]);

      /* check if bounds are exactly one apart */
      if( SCIPisEQ(scip, lb, ub - 1.0) )
      {
         SCIP_VAR* binvar;
         char binvarname[SCIP_MAXSTRLEN];
         SCIP_Bool infeasible;
         SCIP_Bool redundant;
         SCIP_Bool aggregated;

         SCIPdebugMessage("converting <%s>[%g,%g] into binary variable\n", SCIPvarGetName(vars[v]), lb, ub);

         /* create binary variable */
         (void) SCIPsnprintf(binvarname, SCIP_MAXSTRLEN, "%s_bin", SCIPvarGetName(vars[v]));
         SCIP_CALL( SCIPcreateVar(scip, &binvar, binvarname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
               SCIPvarIsInitial(vars[v]), SCIPvarIsRemovable(vars[v]), NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, binvar) );

         /* aggregate integer and binary variable */
         SCIP_CALL( SCIPaggregateVars(scip, vars[v], binvar, 1.0, -1.0, lb, &infeasible, &redundant, &aggregated) );
         assert(!infeasible);
         assert(redundant);
         assert(aggregated);
         
         /* release binary variable */
         SCIP_CALL( SCIPreleaseVar(scip, &binvar) );

         (*nchgvartypes)++;
         *result = SCIP_SUCCESS;
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}





/*
 * presolver specific interface methods
 */

/** creates the inttobinary presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolInttobinary(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;

   /* create inttobinary presolver data */
   presoldata = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_DELAY,
         presolFreeInttobinary, presolInitInttobinary, presolExitInttobinary, 
         presolInitpreInttobinary, presolExitpreInttobinary, presolExecInttobinary,
         presoldata) );

   return SCIP_OKAY;
}
