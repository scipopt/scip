/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_symbreak.cpp
 * @brief  presolver for adding symmetry breaking constraints
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Christopher Hojny
 *
 * This presolver adds the following symmetry breaking constraints:
 *
 * - minimal cover inequalities for symresacks via a constraint handler
 *
 * @note It is important to control the order of presolvers within SCIP in order to avoid contraditions. First, one needs
 * to take care of presolvers that have an effect on symmetry, e.g., presol_domcol. If symmetry is computed on the
 * original formulation, we perform this presolver at the very beginning. Otherwise, we try to compute symmetry as late
 * as possible and then add constraints based on this information.
 *
 * @note Currently, we only allocate memory for pointers to symresack constraints for group generators. If further
 * constraints are considered, we have to reallocate memory.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include <scip/cons_linear.h>
#include <scip/cons_setppc.h>
#include <scip/cons_symresack.h>
#include <scip/prop_probing.h>
#include <scip/presol_symbreak.h>
#include <scip/presol_symmetry.h>
#include <symmetry/type_symmetry.h>

/* presolver properties */
#define PRESOL_NAME            "symbreak"
#define PRESOL_DESC            "presolver for adding symmetry breaking constraints"
#define PRESOL_PRIORITY         -10000000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               -1    /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY                 TRUE    /**< should presolver be delayed, if other presolvers found reductions? */
#define PRESOL_TIMING   SCIP_PRESOLTIMING_EXHAUSTIVE   /**< timing for presolving */

/* default parameters */
#define DEFAULT_CONSSADDLP           TRUE    /**< Should the symmetry breaking constraints be added to the LP? */
#define DEFAULT_ADDSYMRESACKS        TRUE    /**< Add inequalities for symresacks for each generator? */

/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   int**                 perms;              /**< list of permutations */
   int                   nperms;             /**< number of permutations in perms */
   int                   npermvars;          /**< number of variables affected by permutations */
   SCIP_VAR**            permvars;           /**< array of variables on which permutations act */
   int                   nvars;              /**< number of variables for symmetry computations */
   SCIP_VAR**            vars;               /**< variables used for symmetry computations */
   SCIP_Bool             addedconss;         /**< whether we already added symmetry breaking constraints */
   SCIP_Bool             computedsymmetry;   /**< whether symmetry has been computed already */
   SCIP_Bool             conssaddlp;         /**< Should the symmetry breaking constraints be added to the LP? */
   SCIP_Bool             addsymresacks;      /**< Add symresack constraints for each generator? */
   SCIP_PRESOL*          symmetrypresol;     /**< pointer to symmetry presolver */
   SCIP_Bool             enabled;            /**< run presolver? */
   SCIP_Bool             early;              /**< run presolver as early as possible if symmetry has been detected in initpre() */
   SCIP_CONS**           genconss;           /**< list of generated constraints */
   int                   ngenconss;          /**< number of generated constraints */
   int                   nsymresacks;        /**< number of symresack constraints */
};


/*
 * Local methods
 */


/** add symresack constraints */
static
SCIP_RETCODE addSymresackConss(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOL*          presol              /**< symmetry breaking presolver */
   )
{
   SCIP_PRESOLDATA* presoldata;
   int nperms;
   int** perms;
   int nsymresackcons = 0;
   SCIP_VAR** permvars;
   int npermvars;
   SCIP_Bool conssaddlp;
   int p;

   assert( scip != 0 );
   assert( presol != 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   perms = presoldata->perms;
   nperms = presoldata->nperms;
   permvars = presoldata->permvars;
   npermvars = presoldata->npermvars;
   conssaddlp = presoldata->conssaddlp;

   assert( nperms <= 0 || (nperms > 0 && perms != NULL) );
   assert( permvars != NULL );
   assert( npermvars > 0 );

   /* loop through perms and add symresack constraints */
   for (p = 0; p < nperms; ++p)
   {
      SCIP_CONS* cons;
      SCIP_Bool success = FALSE;

      SCIP_CALL( SCIPcreateConsSymresack(scip, &cons, "symresack", (unsigned int*) perms[p], permvars, npermvars, &success,
            conssaddlp, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      /* add the constraint only if the constraint is not trivial */
      if ( success )
      {
         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* do not release constraint here - will be done later */
         presoldata->genconss[presoldata->ngenconss] = cons;
         ++presoldata->ngenconss;
         ++presoldata->nsymresacks;
         ++nsymresackcons;
         SCIPdebugMsg(scip, "Added symresack constraint: %u.\n", nsymresackcons);
      }
      else
      {
         /* otherwise the constraint was not generated */
         assert( cons == NULL );
      }
   }

   return SCIP_OKAY;
}


/** analyze generators and add symmetry breaking constraints */
static
SCIP_RETCODE addSymmetryBreakingConstraints(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_Bool decision;

   assert( scip != 0 );
   assert( presol != 0 );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   /* exit if no or only trivial symmetry group is available */
   if ( presoldata->nperms < 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetBoolParam(scip, "presolving/" PRESOL_NAME"/addsymresacks", &decision) );

   if ( decision )
   {
      SCIP_CALL( addSymresackConss(scip, presol) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeSymbreak)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* free pointers to symmetry group and binary variables */
   if ( presoldata->ngenconss > 0 )
      SCIPfreeBlockMemoryArrayNull(scip, &presoldata->genconss, presoldata->ngenconss);

   SCIPfreeMemory(scip, &presoldata);

   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitSymbreak)
{
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );

   SCIPdebugMsg(scip, "Initpre method of symmetry breaking presolver ...\n");

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   /* check whether we should run */
#if 0
   SCIP_CALL( SCIPgetBoolParam(scip, "usesymbreak", &presoldata->enabled) );
#else
   presoldata->enabled = TRUE;
#endif

   if ( presoldata->enabled )
   {
      /* allow all problem specifications, since we handle them in the code above */
      SYMsetSpecRequirement(presoldata->symmetrypresol, SYM_SPEC_BINARY);
      SYMsetSpecRequirement(presoldata->symmetrypresol, SYM_SPEC_INTEGER);
      SYMsetSpecRequirement(presoldata->symmetrypresol, SYM_SPEC_REAL);
   }

   return SCIP_OKAY;
}


/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitSymbreak)
{
   SCIP_PRESOLDATA* presoldata;
   int i;
   SCIP_CONS** genconss;
   int ngenconss;

   assert( scip != NULL );
   assert( presol != NULL );

   SCIPdebugMessage("Exit method of symmetry breaking presolver ...\n");

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   ngenconss = presoldata->ngenconss;

   if ( ngenconss > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Number of symresacks handled:\t\t\t\t\t%u\n", presoldata->nsymresacks);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Total number of constraints added:\t\t\t\t%u\n\n", presoldata->ngenconss);
   }

   /* release constraints */
   genconss = presoldata->genconss;
   for (i = 0; i < ngenconss; ++i)
   {
      SCIP_CONS* cons = genconss[i];
      assert( cons != NULL );

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}


/** presolving initialization method of presolver (called when presolving is about to begin) */
static
SCIP_DECL_PRESOLINITPRE(presolInitpreSymbreak)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   assert( scip != NULL );
   assert( presol != NULL );

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* check whether we have to run the presolver at the beginning of presolving */
   presoldata->early = ! SYMcomputeSymmetryPresolved(presoldata->symmetrypresol);

   if ( presoldata->early && presoldata->enabled )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Executing presolver <%s> early, since symmetries are computed early.\n\n", SCIPpresolGetName(presol));

      SCIP_CALL( SCIPsetIntParam(scip, "presolving/" PRESOL_NAME"/priority", 90000000) );
   }

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecSymbreak)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   int i;
   int noldfixedvars = *nfixedvars;
   int noldaggrvars = *naggrvars;
   int noldbdchgs = *nchgbds;
   int noldaddconss = *naddconss;


   assert( scip != NULL );
   assert( presol != NULL );
   assert( result != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING );

   *result = SCIP_DIDNOTRUN;

   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != 0 );

   /* possibly skip presolver */
   if ( ! presoldata->enabled )
      return SCIP_OKAY;

   /* skip presolving if we are not at the end */
   if ( ! presoldata->early && ! SCIPisPresolveFinished(scip) )
      return SCIP_OKAY;

   /* possibly stop */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Presolving method of symmetry breaking presolver ...\n");

   /* get symmetry information, if not already computed */
   if ( ! presoldata->computedsymmetry )
   {
      SCIPdebugMsg(scip, "Symmetry breaking presolver: computing symmetry ...\n");
      assert( presoldata->nperms < 0 );

      /* get symmetries */
      SCIP_CALL( SCIPgetSymmetryGenerators(scip, presoldata->symmetrypresol, &(presoldata->npermvars),
            &(presoldata->permvars), &(presoldata->nperms), &(presoldata->perms)) );

      presoldata->computedsymmetry = TRUE;

      if ( presoldata->nperms < 0 )
      {
         SCIPdebugMessage("Symmetry breaking presolver: no symmetry has been found, turning presolver off.\n");
         presoldata->enabled = FALSE;
         return SCIP_OKAY;
      }
      else if ( presoldata->nperms > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata->genconss), presoldata->nperms) );
      }
   }

   /* at this point, the symmetry group should be computed and nontrivial */
   assert( presoldata->nperms > 0 );
   *result = SCIP_DIDNOTFIND;

   /* possibly stop */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* if not already done, add symmetry breaking constraints */
   if ( ! presoldata->addedconss )
   {
      /* add symmetry breaking constraints */
      int noldngenconns = presoldata->ngenconss;

      SCIP_CALL( addSymmetryBreakingConstraints(scip, presol) );

      presoldata->addedconss = TRUE;
      *naddconss += presoldata->ngenconss - noldngenconns;
      SCIPdebugMsg(scip, "Added symmetry breaking constraints: %u.\n", presoldata->ngenconss - noldngenconns);

      /* if constraints have been added, loop through generated constraints and presolve each */
      for (i = 0; i < presoldata->ngenconss; ++i)
      {
         SCIP_CALL( SCIPpresolCons(scip, presoldata->genconss[i], nrounds, SCIP_PRESOLTIMING_ALWAYS, nnewfixedvars, nnewaggrvars, nnewchgvartypes,
               nnewchgbds, nnewholes, nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides, nfixedvars, naggrvars,
               nchgvartypes, nchgbds, naddholes, ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

         /* exit if cutoff has been detected */
         if ( *result == SCIP_CUTOFF || *result == SCIP_UNBOUNDED )
         {
            SCIPdebugMessage("Presolving constraint <%s> detected cutoff or unboundedness.\n", SCIPconsGetName(presoldata->genconss[i]));
            return SCIP_OKAY;
         }
      }
      SCIPdebugMsg(scip, "Presolved %u constraints generated by symbreak.\n", presoldata->ngenconss);

      /* possibly stop */
      if ( SCIPisStopped(scip) )
         return SCIP_OKAY;
   }

   /* determine success */
   if ( noldaddconss + noldbdchgs + noldaggrvars + noldfixedvars > 0 )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the symbreak presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolSymbreak(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata = 0;
   SCIP_PRESOL* presol = 0;

   /* create presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );

   /* we must call the constructor explictly, because memory was m'alloced and not new'ed */
   presoldata->vars = NULL;
   presoldata->nvars = 0;
   presoldata->addedconss = FALSE;
   presoldata->computedsymmetry = FALSE;
   presoldata->enabled = TRUE;
   presoldata->early = FALSE;
   presoldata->nsymresacks = 0;
   presoldata->ngenconss = 0;
   presoldata->genconss = NULL;
   presoldata->nperms = -1;

   /* determine cons_symmetries constraint handler (preuse presol) */
   presol = SCIPfindPresol(scip, "symmetry");
   if ( presol == 0 )
   {
      SCIPerrorMessage("Could not find symmetry presolver.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   presoldata->symmetrypresol = presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecSymbreak, presoldata) );

   /* set additional callbacks */
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeSymbreak) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitSymbreak) );
   SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitSymbreak) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreSymbreak) );

   /* add symmetry breaking presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME"/conssaddlp",
         "Should the symmetry breaking constraints be added to the LP?",
         &presoldata->conssaddlp, TRUE, DEFAULT_CONSSADDLP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME"/addsymresacks",
         "Add inequalities for symresacks for each generator?",
         &presoldata->addsymresacks, TRUE, DEFAULT_ADDSYMRESACKS, NULL, NULL) );

   return SCIP_OKAY;
}
