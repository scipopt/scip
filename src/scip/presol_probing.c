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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: presol_probing.c,v 1.5 2005/02/04 12:51:35 bzfpfend Exp $"

/**@file   presol_probing.c
 * @brief  probing presolver
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "presol_probing.h"


#define PRESOL_NAME            "probing"
#define PRESOL_DESC            "probing presolver on binary variables"
#define PRESOL_PRIORITY        -100000
#define PRESOL_MAXROUNDS       -1




/*
 * Default parameter settings
 */

#define DEFAULT_PROPROUNDS          -1  /**< maximal number of propagation rounds in probing subproblems */




/*
 * Data structures
 */

/** presolver data */
struct PresolData
{
   int              proprounds;         /**< maximal number of propagation rounds in probing subproblems */
};




/*
 * Local methods
 */

/** applies and evaluates probing of a single variable in the given direction */
static
RETCODE applyProbing(
   SCIP*            scip,               /**< SCIP data structure */
   PRESOLDATA*      presoldata,         /**< presolver data */
   VAR**            vars,               /**< problem variables */
   int              nvars,              /**< number of problem variables */
   int              probingpos,         /**< variable number to apply probing on */
   Bool             probingdir,         /**< value to fix probing variable to */
   Real*            lbs,                /**< array to store propagated lower bounds */
   Real*            ubs,                /**< array to store propagated upper bounds */
   Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   )
{
   assert(presoldata != NULL);
   assert(lbs != NULL);
   assert(ubs != NULL);
   assert(cutoff != NULL);
   assert(0 <= probingpos && probingpos < nvars);
   assert(SCIPvarGetType(vars[probingpos]) == SCIP_VARTYPE_BINARY);
   assert(SCIPvarGetLbLocal(vars[probingpos]) < 0.5);
   assert(SCIPvarGetUbLocal(vars[probingpos]) > 0.5);

   /*debugMessage("applying probing on variable <%s> == %d\n", SCIPvarGetName(vars[probingpos]), probingdir);*/

   /* start probing mode */
   CHECK_OKAY( SCIPstartProbing(scip) );

   /* fix variable */
   if( probingdir == FALSE )
   {
      CHECK_OKAY( SCIPchgVarUbProbing(scip, vars[probingpos], 0.0) );
   }
   else
   {
      CHECK_OKAY( SCIPchgVarLbProbing(scip, vars[probingpos], 1.0) );
   }

   /* apply propagation */
   CHECK_OKAY( SCIPpropagateProbing(scip, presoldata->proprounds, cutoff) );
   
   /* evaluate propagation */
   if( !(*cutoff) )
   {
      int i;

      for( i = 0; i < nvars; ++i )
      {
         lbs[i] = SCIPvarGetLbLocal(vars[i]);
         ubs[i] = SCIPvarGetUbLocal(vars[i]);
#if 0
#ifdef DEBUG
         if( !SCIPisEQ(scip, lbs[i], SCIPvarGetLbGlobal(vars[i]))
            || !SCIPisEQ(scip, ubs[i], SCIPvarGetUbGlobal(vars[i])) )
         {
            debugMessage("  -> <%s>[%g,%g]\n", SCIPvarGetName(vars[i]), lbs[i], ubs[i]);
         }
#endif
#endif
      }
   }
   else
   {
      /*debugMessage("  -> probing determined a cutoff\n");*/
   }

   /* exit probing mode */
   CHECK_OKAY( SCIPendProbing(scip) );

   return SCIP_OKAY;
}




/*
 * Callback methods of presolver
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
DECL_PRESOLFREE(presolFreeProbing)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   SCIPfreeMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
#define presolInitProbing NULL


/** deinitialization method of presolver (called before transformed problem is freed) */
#define presolExitProbing NULL


/** presolving initialization method of presolver (called when presolving is about to begin) */
#define presolInitpreProbing NULL


/** presolving deinitialization method of presolver (called after presolving has been finished) */
#define presolExitpreProbing NULL


/** execution method of presolver */
static
DECL_PRESOLEXEC(presolExecProbing)
{  /*lint --e{715}*/
   PRESOLDATA* presoldata;
   VAR** probvars;
   VAR** vars;
   Real* zerolbs;
   Real* zeroubs;
   Real* onelbs;
   Real* oneubs;
   int nvars;
   int nbinvars;
   int oldnfixedvars;
   int oldnaggrvars;
   int i;
   Bool cutoff;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get variable data */
   CHECK_OKAY( SCIPgetVarsData(scip, &probvars, &nvars, &nbinvars, NULL, NULL, NULL) );
   if( nbinvars == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;

   /* get presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* copy the vars array, because SCIPfixVar() renders the vars of SCIPgetVarsData() array invalid */
   CHECK_OKAY( SCIPduplicateBufferArray(scip, &vars, probvars, nvars) );

   /* get temporary memory for storing probing results */
   CHECK_OKAY( SCIPallocBufferArray(scip, &zerolbs, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &zeroubs, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &onelbs, nvars) );
   CHECK_OKAY( SCIPallocBufferArray(scip, &oneubs, nvars) );

   /* for each binary variable, probe fixing the variable to zero and one */
   cutoff = FALSE;
   for( i = 0; i < nbinvars && !cutoff && !SCIPpressedCtrlC(scip); ++i )
   {
      Bool localcutoff;
      int j;

#if 0
      if( (i+1) % 1000 == 0 )
      {
         SCIPmessage(scip, SCIP_VERBLEVEL_FULL, "probing: %d/%d (%.1f%%)\n", 
            i+1, nbinvars, 100.0*(Real)(i+1)/(Real)nbinvars);
      }
#endif

      /* ignore variables, that were fixed or aggregated in prior probings */
      if( !SCIPvarIsActive(vars[i]) )
         continue;

      /* apply probing for fixing the variable to zero */
      CHECK_OKAY( applyProbing(scip, presoldata, vars, nvars, i, FALSE, zerolbs, zeroubs, &localcutoff) );
      if( localcutoff )
      {
         Bool fixed;
         
         /* the variable can be fixed to TRUE */
         debugMessage("fixing variable <%s> to 1.0\n", SCIPvarGetName(vars[i]));
         CHECK_OKAY( SCIPfixVar(scip, vars[i], 1.0, &cutoff, &fixed) );
         assert(fixed);
         (*nfixedvars)++;
         if( cutoff )
            break;
         continue; /* don't try upwards direction, because the variable is already fixed */
      }

      /* apply probing for fixing the variable to one */
      CHECK_OKAY( applyProbing(scip, presoldata, vars, nvars, i, TRUE, onelbs, oneubs, &localcutoff) );
      if( localcutoff )
      {
         Bool fixed;
         
         /* the variable can be fixed to FALSE */
         debugMessage("fixing variable <%s> to 0.0\n", SCIPvarGetName(vars[i]));
         CHECK_OKAY( SCIPfixVar(scip, vars[i], 0.0, &cutoff, &fixed) );
         assert(fixed);
         (*nfixedvars)++;
         if( cutoff )
            break;
         continue; /* don't analyze probing deductions, because the variable is already fixed */
      }

      /* analyze probing deductions */
      for( j = 0; j < nbinvars && !cutoff; ++j )
      {
         Bool fixed;
         Bool redundant;
         Bool aggregated;

         if( j == i )
            continue;

         if( zeroubs[j] < 0.5 && oneubs[j] < 0.5 )
         {
            /* in both probings, variable j is deduced to 0: fix variable to 0 */
            debugMessage("fixing variable <%s> to 0\n", SCIPvarGetName(vars[j]));
            CHECK_OKAY( SCIPfixVar(scip, vars[j], 0.0, &cutoff, &fixed) );
            if( fixed )
               (*nfixedvars)++;
         }
         else if( zerolbs[j] > 0.5 && onelbs[j] > 0.5 )
         {
            /* in both probings, variable j is deduced to 1: fix variable to 1 */
            debugMessage("fixing variable <%s> to 1\n", SCIPvarGetName(vars[j]));
            CHECK_OKAY( SCIPfixVar(scip, vars[j], 1.0, &cutoff, &fixed) );
            if( fixed )
               (*nfixedvars)++;
         }
         else if( zeroubs[j] < 0.5 && onelbs[j] > 0.5 )
         {
            /* variable j is always deduced to the same value as probing variable i:
             * both variables can be aggregated with x_i - x_j == 0
             */
            debugMessage("aggregating variables <%s> == <%s>\n", SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));
            CHECK_OKAY( SCIPaggregateVars(scip, vars[i], vars[j], 1.0, -1.0, 0.0, &cutoff, &redundant, &aggregated) );
            if( aggregated )
               (*naggrvars)++;
         }
         else if( zerolbs[j] > 0.5 && oneubs[j] < 0.5 )
         {
            /* variable j is always deduced to the opposite value of probing variable i:
             * both variables can be aggregated with x_i + x_j == 1
             */
            debugMessage("aggregating variables <%s> == 1 - <%s>\n", SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));
            CHECK_OKAY( SCIPaggregateVars(scip, vars[i], vars[j], 1.0, 1.0, 1.0, &cutoff, &redundant, &aggregated) );
            if( aggregated )
               (*naggrvars)++;
         }
         else if( zeroubs[j] < 0.5 )
         {
            /* insert implication: x_i == 0  =>  x_j == 0 */
            debugMessage("found implication <%s> == 0  =>  <%s> == 0\n", SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));
            CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], FALSE, vars[j], SCIP_BOUNDTYPE_UPPER, 0.0, &cutoff) );
         }
         else if( zerolbs[j] > 0.5 )
         {
            /* insert implication: x_i == 0  =>  x_j == 1 */
            debugMessage("found implication <%s> == 0  =>  <%s> == 1\n", SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));
            CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], FALSE, vars[j], SCIP_BOUNDTYPE_LOWER, 1.0, &cutoff) );
         }
         else if( oneubs[j] < 0.5 )
         {
            /* insert implication: x_i == 1  =>  x_j == 0 */
            debugMessage("found implication <%s> == 1  =>  <%s> == 0\n", SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));
            CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], TRUE, vars[j], SCIP_BOUNDTYPE_UPPER, 0.0, &cutoff) );
         }
         else if( onelbs[j] > 0.5 )
         {
            /* insert implication: x_i == 1  =>  x_j == 1 */
            debugMessage("found implication <%s> == 1  =>  <%s> == 1\n", SCIPvarGetName(vars[i]), SCIPvarGetName(vars[j]));
            CHECK_OKAY( SCIPaddVarImplication(scip, vars[i], TRUE, vars[j], SCIP_BOUNDTYPE_LOWER, 1.0, &cutoff) );
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &oneubs);
   SCIPfreeBufferArray(scip, &onelbs);
   SCIPfreeBufferArray(scip, &zeroubs);
   SCIPfreeBufferArray(scip, &zerolbs);
   SCIPfreeBufferArray(scip, &vars);

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}





/*
 * presolver specific interface methods
 */

/** creates the probing presolver and includes it in SCIP */
RETCODE SCIPincludePresolProbing(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   PRESOLDATA* presoldata;

   /* create probing presolver data */
   CHECK_OKAY( SCIPallocMemory(scip, &presoldata) );

   /* include presolver */
   CHECK_OKAY( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         presolFreeProbing, presolInitProbing, presolExitProbing, 
         presolInitpreProbing, presolExitpreProbing, presolExecProbing,
         presoldata) );

   /* add probing presolver parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
         "presolving/probing/proprounds", 
         "maximal number of propagation rounds in probing subproblems (-1: no limit, 0: auto)",
         &presoldata->proprounds, DEFAULT_PROPROUNDS, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
