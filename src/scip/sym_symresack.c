/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sym_symresack.c
 * @brief  symmetry handler for symresack constraints
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_symresack.h"
#include "scip/pub_sym.h"
#include "scip/scip_sym.h"
#include "scip/sym_symresack.h"

/* symmetry handler properties */
#define SYM_NAME            "symresack"
#define SYM_DESC            "symmetry handler for symresack constraint"
#define SYM_PRIORITY          -1000000           /**< propagator priority */
#define SYM_FREQ                     1           /**< propagator frequency */


/** symmetry component data */
struct SCIP_SymCompData
{
   SCIP_CONS**           conss;              /**< constraints added by the symmetry handler */
   int                   nconss;             /**< number of constraints added by the symmetry handler */
};


/** addition method for symmetry method handler plugins (tries to add symmetry handling method for given symmetries) */
static
SCIP_DECL_SYMHDLRTRYADD(symhdlrTryaddSymresack)
{  /*lint --e{715}*/
   char name[SCIP_MAXSTRLEN];
   int p;

   if( symtype != SYM_SYMTYPE_PERM )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   *success = TRUE;

   SCIP_CALL( SCIPallocBlockMemory(scip, symcompdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*symcompdata)->conss, nperms) );
   (*symcompdata)->nconss = nperms;

   for( p = 0; p < nperms; ++p )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_compnent_%d_%d", id, p);

      SCIP_CALL( SCIPcreateSymbreakCons(scip, &(*symcompdata)->conss[p], "cons", perms[p],
            permvars, npermvars, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, (*symcompdata)->conss[p]) );
      ++(*naddedconss);
      /* do not release constraints here, this will be done later */
   }

   return SCIP_OKAY;
}

/** destructor of symmetry handler to free symmetry handler data (called when SCIP is exiting) */
static
SCIP_DECL_SYMHDLRFREE(symhdlrFreeSymresack)
{  /*lint --e{715}*/
   SCIP_SYMHDLRDATA* symhdlrdata;

   symhdlrdata = SCIPsymhdlrGetData(symhdlr);
   assert(symhdlrdata == NULL);

   return SCIP_OKAY;
}

/** deinitialization method of symmetry handler (called before transformed problem is freed) */
static
SCIP_DECL_SYMHDLREXIT(symhdlrExitSymresack)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   for( s = 0; s < nsymcomps; ++s )
   {
      symdata = symcompdata[s];
      assert(symdata->conss == NULL || symdata->nconss > 0);

      for( c = 0; c < symdata->nconss; ++c )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &symdata->conss[c]) );
      }

      SCIPfreeBlockMemoryArrayNull(scip, &symdata->conss, symdata->nconss);
      SCIPfreeBlockMemory(scip, &symdata);
   }

   return SCIP_OKAY;
}

/** presolving method of symmetry handler */
static
SCIP_DECL_SYMHDLRPRESOL(symhdlrPresolSymreack)
{  /*lint --e{715}*/
   SCIP_SYMCOMPDATA* symdata;
   int s;
   int c;

   *result = SCIP_DIDNOTFIND;

   for( s = 0; s < nsymcomps; ++s )
   {
      symdata = symcompdata[s];

      if( symdata->nconss == 0 )
         return SCIP_OKAY;

      for( c = 0; c < symdata->nconss; ++c )
      {
         SCIP_CALL( SCIPpresolCons(scip, symdata->conss[c], nrounds, presoltiming, nnewfixedvars, nnewaggrvars,
               nnewchgvartypes, nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs,
               nnewchgsides, nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes, ndelconss, naddconss,
               nupgdconss, nchgcoefs, nchgsides, result) );

         /* exit if cutoff or unboundedness has been detected */
         if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
         {
            SCIPdebugMsg(scip, "Presolving constraint <%s> detected cutoff or unboundedness.\n",
               SCIPconsGetName(symdata->conss[c]));
            return SCIP_OKAY;
         }
      }

   }

   return SCIP_OKAY;
}


/** include symmetry handler for symresack constraints */
SCIP_RETCODE SCIPincludeSymhdlrSymresack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYMHDLRDATA* symhdlrdata = NULL;

   assert(scip != NULL);

   SCIP_CALL( SCIPincludeSymhdlrBasic(scip, SYM_NAME, SYM_DESC,
         1, 1, 1, 1, -1, -1, FALSE, FALSE, -1, SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_FAST,
         symhdlrTryaddSymresack, NULL, symhdlrFreeSymresack, NULL, symhdlrExitSymresack,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, symhdlrPresolSymreack, symhdlrdata) );

   return SCIP_OKAY;
}
