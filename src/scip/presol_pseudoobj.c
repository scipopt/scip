/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: presol_pseudoobj.c,v 1.42 2011/01/02 11:10:44 bzfheinz Exp $"

/**@file   presol_pseudoobj.c
 * @ingroup PRESOLVERS
 * @brief  pseudoobj presolver
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/presol_pseudoobj.h"
#include "scip/prop_pseudoobj.h"


#define PRESOL_NAME            "pseudoobj"
#define PRESOL_DESC            "pseudoobj presolver: round fractional bounds on integers, fix variables with equal bounds"
#define PRESOL_PRIORITY        +6000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY              FALSE /**< should presolver be delayed, if other presolvers found reductions? */


/** presolver data */
struct SCIP_PresolData
{
   SCIP_Real             cutoffbound;        /**< last cutoff bound used */
   SCIP_Real             pseudoobjval;       /**< last pseudo objective used */
};


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyPseudoobj)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolPseudoobj(scip) );
 
   return SCIP_OKAY;
}


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreePseudoobj)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   
   SCIPfreeMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
#define presolInitPseudoobj NULL


/** deinitialization method of presolver (called before transformed problem is freed) */
#define presolExitPseudoobj NULL


/** presolving initialization method of presolver (called when presolving is about to begin) */
#define presolInitprePseudoobj NULL


/** presolving deinitialization method of presolver (called after presolving has been finished) */
#define presolExitprePseudoobj NULL


/** presolving execution method */
static
SCIP_DECL_PRESOLEXEC(presolExecPseudoobj)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   SCIP_VAR** vars;
   SCIP_Real cutoffbound;
   SCIP_Real pseudoobjval;
   int oldnchgbds;
   int nvars;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get current pseudo objective value and cutoff bound */
   pseudoobjval = SCIPgetPseudoObjval(scip);
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;
   cutoffbound = SCIPgetCutoffbound(scip);
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;
   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   
   if( cutoffbound < presoldata->cutoffbound || pseudoobjval > presoldata->pseudoobjval )
   {
      *result = SCIP_DIDNOTFIND;
      oldnchgbds = *nchgbds;

      /* get the problem variables */
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);
      
      /*  scan the variables for pseudoobj bound reductions
       * (loop backwards, since a variable fixing can change the current and the subsequent slots in the vars array)
       */
      for( v = nvars-1; v >= 0; --v )
      {
         if( SCIPisZero(scip, SCIPvarGetObj(vars[v])) )
            continue;
         
         SCIP_CALL( SCIPpropagateCutoffboundVar(scip, NULL, vars[v], cutoffbound, pseudoobjval, nchgbds) );
      }
      
      if( *nchgbds > oldnchgbds )
         *result = SCIP_SUCCESS;

      /* stroe the old values */
      presoldata->cutoffbound = cutoffbound;
      presoldata->pseudoobjval = pseudoobjval;
   }
   
   return SCIP_OKAY;
}





/*
 * presolver specific interface methods
 */

/** creates the pseudoobj presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolPseudoobj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;

   /* create pseudoobj presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );

   presoldata->cutoffbound = SCIPinfinity(scip);
   presoldata->pseudoobjval = -SCIPinfinity(scip);
   
   /* include presolver */
   SCIP_CALL( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_DELAY,
         presolCopyPseudoobj,
         presolFreePseudoobj, presolInitPseudoobj, presolExitPseudoobj, 
         presolInitprePseudoobj, presolExitprePseudoobj, presolExecPseudoobj,
         presoldata) );

   return SCIP_OKAY;
}
