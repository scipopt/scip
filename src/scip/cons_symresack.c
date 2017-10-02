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

/**@file   cons_symresack.c
 * @brief  constraint handler for symresack constraints
 * @author Christopher Hojny
 *
 * @TODO: Currently, the copy methods of the constraint handler are deactivated, as otherwise, we would reduce the
 * effect of heuristics. The trade-off for this is that we cannot transfer the symmetry information to the sub-scips
 * of the components presolver.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/cons_setppc.h"
#include "scip/cons_orbisack.h"
#include "scip/cons_symresack.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "symresack"
#define CONSHDLR_DESC          "symmetry breaking constraint handler relying on symresacks"
#define CONSHDLR_SEPAPRIORITY    +40100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1005200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1005200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             5 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             5 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_EXHAUSTIVE

#define DEFAULT_ENFORCING         FALSE /**< whether we use enforcing methods of constraint handler */
#define DEFAULT_CHECK             FALSE /**< whether we use check methods of constraint handler */
#define DEFAULT_UPGRADE            TRUE /**< whether we allow upgrading to orbisack constraints */
#define DEFAULT_PPSYMRESACK       FALSE /**< whether we allow upgrading to packing/partitioning symresacks */

/* macros for getting bounds of pseudo solutions in propagation */
#define ISFIXED0(x)   (SCIPvarGetUbLocal(x) < 0.5 ? TRUE : FALSE)
#define ISFIXED1(x)   (SCIPvarGetLbLocal(x) > 0.5 ? TRUE : FALSE)


/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             symresackEnforcing; /**< whether we use enforcing methods of constraint hanlder */
   SCIP_Bool             symresackCheck;     /**< whether we use check methods of constraint hanlder */
   SCIP_Bool             symresackUpgrade;   /**< whether we allow upgrading symresack constraints to orbisack constraints */
   SCIP_Bool             checkPPsymresack;   /**< whether we allow upgrading to packing/partitioning symresacks */
};


/** constraint data for symresack constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables */
   unsigned int          nVars;              /**< number of variables */
   SCIP_Real*            vals;               /**< LP-solution for the variables */
   unsigned int*         perm;               /**< permutation associated to the symresack */
   unsigned int*         invPerm;            /**< inverse permutation */
   SCIP_Bool             ppUpgrade;          /**< whether constraint is upgraded to packing/partitioning symresack */

   /* data for upgraded symresack constraints */
   unsigned int          nCycles;            /**< number of cycles in permutation */
   unsigned int**        cycleDecomposition; /**< cycle decomposition */
};


/*
 * Local methods
 */

/** frees a symresack constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to symresack constraint data */
   )
{
   unsigned int nVars;
   unsigned int i;

   assert( consdata != NULL );
   assert( *consdata != NULL );

   nVars = (*consdata)->nVars;

   if ( (*consdata)->ppUpgrade )
   {
      for (i = 0; i < (*consdata)->nCycles; ++i)
      {
         SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cycleDecomposition[i]), nVars + 1);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cycleDecomposition), (*consdata)->nCycles);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vals), nVars);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->invPerm), nVars);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->perm), nVars);

   for (i = 0; i < nVars; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars), nVars);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** check whether constraint can be upgraded to packing/partitioning symresack */
static
SCIP_RETCODE packingUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   unsigned int*         perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables affected by permutation */
   unsigned int          n,                  /**< length of permutation */
   SCIP_Bool*            upgrade             /**< whether upgrade was succesful */
   )
{
   SCIP_Bool* covered;
   unsigned int nCycles = 0;
   unsigned int i;
   unsigned int j;
   SCIP_Bool descent;
   unsigned int** cycleDecomposition;
   unsigned int curCycle;
   unsigned int maxCycleLength;
   unsigned int cycleLength;
   SCIP_CONSHDLR* setppcconshdlr;
   SCIP_CONS** setppcconss;
   int nsetppcconss;
   int* indicesInCycle;
   SCIP_VAR* var;
   int c;
   SCIP_Bool terminated = FALSE;

   assert( scip != NULL );
   assert( perm != NULL );
   assert( vars != NULL );
   assert( n > 0 );
   assert( upgrade != NULL );

   *upgrade = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &covered, n) );

   for (i = 0; i < n; ++i)
      covered[i] = FALSE;

   /* check wether permutation is monotone */
   for (i = 0; i < n; ++i)
   {
      /* skip checked indices */
      if ( covered[i] )
         continue;

      ++nCycles;
      j = i;
      descent = FALSE;

      do
      {
         covered[j] = TRUE;

         if ( perm[j] < j )
         {
            if ( ! descent )
               descent = TRUE;
            else
               break;
         }

         j = perm[j];
      }
      while ( j != i );

      /* if cycle is not monotone */
      if ( j != i )
      {
         SCIPfreeBufferArray(scip, &covered);

         return SCIP_OKAY;
      }
   }
   assert( nCycles <= n / 2 );

   /* each cycle is monotone; check for packing/partitioning type */
   for (i = 0; i < n; ++i)
      covered[i] = FALSE;

   /* compute cycle decomposition: row i stores in entry 0 the length of the cycle,
    * the remaining entries are the coordinates in the cycle */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cycleDecomposition, nCycles) );
   for (i = 0; i < nCycles; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cycleDecomposition[i], n + 1) );
   }

   curCycle = 0;
   maxCycleLength = 0;
   for (i = 0; i < n; ++i)
   {
      /* skip checked indices */
      if ( covered[i] )
         continue;

      j = i;
      cycleLength = 0;
      do
      {
         covered[j] = TRUE;
         cycleDecomposition[curCycle][++cycleLength] = j;
         j = perm[j];
      }
      while ( j != i );

      cycleDecomposition[curCycle][0] = cycleLength;
      ++curCycle;

      if ( maxCycleLength < cycleLength )
         maxCycleLength = cycleLength;
   }

   /* permutation can be upgraded -> check whether the symresack is of packing/partitioning type */
   setppcconshdlr = SCIPfindConshdlr(scip, "setppc");
   assert( setppcconshdlr != 0 );
   setppcconss = SCIPconshdlrGetConss(setppcconshdlr);
   nsetppcconss = SCIPconshdlrGetNConss(setppcconshdlr);

   /* Check whether each cycle of the symresack is contained in a set packing/partitioning constraint.
    * To this end, we have to guarantee that all affected variables are not negated since permutations
    * are given w.r.t. original variables. */
   *upgrade = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &indicesInCycle, maxCycleLength) );

   for (i = 0; i < nCycles && *upgrade && (! terminated); ++i)
   {
      /* get indices of variables in current cycle */
      for (j = 0; j < cycleDecomposition[i][0]; ++ j)
      {
         var = vars[cycleDecomposition[i][j + 1]];

         if ( SCIPvarIsNegated(var) )
         {
            terminated = TRUE;
            break;
         }

         indicesInCycle[j] = SCIPvarGetProbindex(var);
      }

      cycleLength = cycleDecomposition[i][0];

      /* iterate over constraints */
      for (c = 0; c < nsetppcconss; ++c)
      {
         int nsetppcvars;
         SCIP_VAR** setppcvars;
         int varidx;
         unsigned int nFound = 0;
         unsigned int k;

         /* check type */
         if ( SCIPgetTypeSetppc(scip, setppcconss[c]) != SCIP_SETPPCTYPE_PARTITIONING &&
            SCIPgetTypeSetppc(scip, setppcconss[c]) != SCIP_SETPPCTYPE_PACKING )
            continue;

         /* get set packing/partitioning variables */
         nsetppcvars = SCIPgetNVarsSetppc(scip, setppcconss[c]);
         assert( nsetppcvars > 0 );

         setppcvars = SCIPgetVarsSetppc(scip, setppcconss[c]);
         assert( setppcvars != NULL );

         /* check whether all variables of the cycle are contained in setppc constraint */
         for (j = 0; (int) j < nsetppcvars && nFound < cycleLength; ++j)
         {
            var = setppcvars[j];

            if ( SCIPvarIsNegated(var) )
               continue;

            varidx = SCIPvarGetProbindex(var);

            for (k = 0; k < cycleLength; ++k)
            {
               if ( varidx == indicesInCycle[k] )
               {
                  ++nFound;
                  break;
               }
            }
         }

         if ( nFound == cycleLength )
            break;
      }

      /* row is not contained in a set packing/partitioning constraint */
      if ( c >= nsetppcconss )
         *upgrade = FALSE;
   }

   if ( *upgrade )
   {
      (*consdata)->nCycles = nCycles;
      (*consdata)->cycleDecomposition = cycleDecomposition;

      SCIPfreeBufferArray(scip, &indicesInCycle);
      SCIPfreeBufferArray(scip, &covered);

      return SCIP_OKAY;
   }

   SCIPfreeBufferArray(scip, &indicesInCycle);
   for (i = 0; i < nCycles; ++i)
   {
      SCIPfreeBlockMemoryArray(scip, &cycleDecomposition[i], n + 1);
   }
   SCIPfreeBlockMemoryArray(scip, &cycleDecomposition, nCycles);
   SCIPfreeBufferArray(scip, &covered);

   return SCIP_OKAY;
}


/** creates symresack constraint data
 *
 *  If the input data contain non-binary variables of fixed
 *  points, we delete these variables in a preprocessing step.
 */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_VAR*const*       inputVars,          /**< input variables of the constraint handler */
   unsigned int          inputNVars,         /**< input number of variables of the constraint handler*/
   unsigned int*         inputPerm           /**< input permutation of the constraint handler */
   )
{
   unsigned int i;
   unsigned int j = 0;
   int* indexCorrection;
   int nAffectedVariables;
   SCIP_VAR** vars;
   unsigned int* perm;
   unsigned int* invPerm;
   SCIP_Bool upgrade;

   assert( consdata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   /* count the number of binary variables which are affected by the permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &indexCorrection, inputNVars) );
   indexCorrection[0] = -1;
   for (i = 0; i < inputNVars; ++i)
   {
      if ( inputPerm[i] != i && SCIPvarIsBinary(inputVars[i]) )
      {
         if ( i == 0 )
            indexCorrection[i] = 0;
         else
            indexCorrection[i] = indexCorrection[i - 1] + 1;
      }
      else
      {
         if ( i > 0 )
            indexCorrection[i] = indexCorrection[i - 1];
      }
   }
   nAffectedVariables = indexCorrection[inputNVars - 1] + 1;

   (*consdata)->nVars = nAffectedVariables;

   /* Stop if we detect that the permutation fixes each binary point. */
   if ( nAffectedVariables == 0 )
   {
      SCIPfreeBufferArrayNull(scip, &indexCorrection);
      return SCIP_OKAY;
   }

   /* remove fixed points from permutation representation */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, nAffectedVariables) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &perm, nAffectedVariables) );
   for (i = 0; i < inputNVars; ++i)
   {
      if ( i == 0 )
      {
         if ( indexCorrection[i] > -1 )
         {
            vars[j] = inputVars[i];
            perm[j++] = indexCorrection[inputPerm[i]];
         }
      }
      else
      {
         if ( indexCorrection[i] > indexCorrection[i - 1] )
         {
            vars[j] = inputVars[i];
            perm[j++] = indexCorrection[inputPerm[i]];
         }
      }
   }
   (*consdata)->vars = vars;
   (*consdata)->perm = perm;

   for (i = 0; (int) i < nAffectedVariables; ++i)
   {
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[i]) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &invPerm, nAffectedVariables) );
   for (i = 0; (int) i < nAffectedVariables; ++i)
      invPerm[perm[i]] = i;
   (*consdata)->invPerm = invPerm;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vals, nAffectedVariables) );

   SCIPfreeBufferArrayNull(scip, &indexCorrection);

   /* check whether an upgrade to packing/partitioning symresacks is possible */
   upgrade = FALSE;
   SCIP_CALL( SCIPgetBoolParam(scip, "cons/symresack/ppsymresack", &upgrade) );

   if ( upgrade )
   {
      upgrade = FALSE;
      SCIP_CALL( packingUpgrade(scip, consdata, perm, vars, nAffectedVariables, &upgrade) );
   }

   (*consdata)->ppUpgrade = upgrade;

   return SCIP_OKAY;
}


/** generate initial LP cut
 *
 *  We generate the ordering inequality for the pair \f$(1, \gamma^{-1}(1))\f$, i.e.,
 *  the inequality \f$-x_{1} + x_{\gamma^{-1}(1)} \leq 0\f$. This inequality is valid,
 *  because we guaranteed in a preprocessing stept that all variables are binary.
 *
 *  Furthermore, we add facet inequalities of packing/partitioning symresacks if
 *  we deal with packing/partitioning symresacks.
 */
static
SCIP_RETCODE initLP(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   SCIP_CONSDATA* consdata;
   unsigned int nVars;
   SCIP_VAR** vars;
   SCIP_ROW* row;
   unsigned int i, j, k;
   unsigned int nCycles;
   unsigned int** cycleDecomposition;
   SCIP_VAR** varsInCons;
   SCIP_Real* coeffs;
   unsigned int nVarsInCons;
   unsigned int nVarsInCycle;
   unsigned int firstElemInCycle;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != 0 );

   nVars = consdata->nVars;
   vars = consdata->vars;

   /* avoid stupid problems */
   if ( nVars <= 1 )
      return SCIP_OKAY;

   /* there are no fixed points */
   assert( consdata->invPerm[0] != 0 );

   /* add ordering inequality */
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "symresack_init", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, vars[0], -1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, vars[consdata->invPerm[0]], 1.0) );

   SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, infeasible) );

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   /* check whether we have a packing/partioning symresack */
   if ( consdata->ppUpgrade )
   {
      nCycles = consdata->nCycles;
      cycleDecomposition = consdata->cycleDecomposition;

      SCIP_CALL( SCIPallocBufferArray(scip, &varsInCons, nVars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &coeffs, nVars) );

      coeffs[0] = 1.0;

      /* add packing/partioning symresack constraints */
      for (i = 0; i < nCycles; ++i)
      {
         assert( cycleDecomposition[i][0] > 0 );

         nVarsInCycle = cycleDecomposition[i][0];
         varsInCons[0] = vars[cycleDecomposition[i][nVarsInCycle]];
         firstElemInCycle = cycleDecomposition[i][1];

         assert( firstElemInCycle == consdata->perm[cycleDecomposition[i][nVarsInCycle]] );

         nVarsInCons = 1;

         /* add variables of other cycles to the constraint */
         for (j = 0; j < i; ++j)
         {
            nVarsInCycle = cycleDecomposition[j][0];
            for (k = 1; k <= nVarsInCycle; ++k)
            {
               if ( cycleDecomposition[j][k] < firstElemInCycle )
               {
                  varsInCons[nVarsInCons] = vars[cycleDecomposition[j][k]];
                  coeffs[nVarsInCons++] = -1.0;
               }
               else
                  continue;
            }
         }

         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "ppSymresack", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPaddVarsToRow(scip, row, nVarsInCons, varsInCons, coeffs) );

         SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, infeasible) );
      }

      SCIPfreeBufferArray(scip, &coeffs);
      SCIPfreeBufferArray(scip, &varsInCons);
   }

   return SCIP_OKAY;
}


/** propagation */
static
SCIP_RETCODE propVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_Bool*            infeasible,         /**< whether it was detected that the node is infeasible */
   SCIP_Bool*            found,              /**< whether a new propagation could be found */
   unsigned int*         nGen                /**< number of generated bound strengthenings */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   unsigned int nVars;
   unsigned int* invPerm;
   SCIP_Bool tightened;
   unsigned int i;
   SCIP_VAR* var;
   SCIP_VAR* var2;
   unsigned int r;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nGen != NULL );
   assert( infeasible != NULL );
   assert( found != NULL );

   SCIPdebugMessage("Propagating variables of constraint <%s>.\n", SCIPconsGetName(cons));

   *nGen = 0;
   *infeasible = FALSE;
   *found = FALSE;

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nVars != 0 );
   assert( consdata->invPerm != NULL );

   vars = consdata->vars;
   nVars = consdata->nVars;
   invPerm = consdata->invPerm;

   /* avoid trivial problems */
   if ( nVars < 2 )
      return SCIP_OKAY;

   /* loop through all variables */
   for (i = 0; i < nVars; ++i)
   {
      /* there are no fixed points */
      assert( invPerm[i] != i );

      /* get variables of first and second column */
      var = vars[i];
      var2 = vars[invPerm[i]];
      assert( var != NULL );
      assert( var2 != NULL );

      /* if first part of variable pair fixed to 0 and second part is fixed to 1 */
      if ( ISFIXED0(var) && ISFIXED1(var2) )
      {
         SCIPdebugMessage("Check variable pair (%u,%u).\n", i, invPerm[i]);

         SCIPdebugMessage(" -> node infeasible (pair was fixed to (0,1) but there was no pair of type (1,0) before).\n");

         /* perform conflict analysis */
         if ( SCIPisConflictAnalysisApplicable(scip) )
         {
#if ( SCIP_VERSION > 321 || (SCIP_VERSION == 321 && SCIP_SUBVERSION >= 2 ))
            SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
#else
            SCIP_CALL( SCIPinitConflictAnalysis(scip) );
#endif

            for (r = 0; r <= i; ++r)
            {
               /* there are no fixed points */
               assert( invPerm[r] != r );

               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[r]) );
               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[invPerm[r]]) );
            }

            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

            *infeasible = TRUE;
            break;
         }
      }
      /* if first part of the variable pair is fixed to 0 and the second part is free --> fix second part to 0 */
      else if ( ISFIXED0(var) && ( ! ISFIXED0(var2) ) )
      {
         assert( SCIPvarGetUbLocal(var) < 0.5 );
         assert( SCIPvarGetLbLocal(var2) < 0.5 );
         assert( SCIPvarGetUbLocal(var2) > 0.5 );

         SCIPdebugMessage("Check variable pair (%u,%u).\n", i, invPerm[i]);

         assert( SCIPvarGetLbLocal(var2) < 0.5 );
         SCIP_CALL( SCIPinferVarUbCons(scip, var2, 0.0, cons, i, FALSE, infeasible, &tightened) ); /*lint !e713*/
         assert( ! *infeasible );

         *found = *found || tightened;
         if ( tightened )
            ++(*nGen);
      }
      /* if second part of the variable pair is fixed to 1 and the first part is free --> fix first part to 1 */
      else if ( ( ! ISFIXED1(var) ) && ISFIXED1(var2) )
      {
         assert( SCIPvarGetLbLocal(var) < 0.5 );
         assert( SCIPvarGetUbLocal(var) > 0.5 );
         assert( SCIPvarGetLbLocal(var2) > 0.5 );

         SCIPdebugMessage("Check variable pair (%u,%u).\n", i, invPerm[i]);

         assert( SCIPvarGetUbLocal(var) > 0.5 );
         SCIP_CALL( SCIPinferVarLbCons(scip, var, 1.0, cons, i + nVars, FALSE, infeasible, &tightened) ); /*lint !e713*/
         assert( ! *infeasible );

         *found = *found || tightened;
         if ( tightened )
            ++(*nGen);
      }
      /* if solution is lexicographically maximal */
      else if ( ISFIXED1(var) && ISFIXED0(var2) )
      {
         assert( SCIPvarGetLbLocal(var) > 0.5 );
         assert( SCIPvarGetUbLocal(var2) < 0.5 );

         SCIPdebugMessage("Check variable pair (%u,%u).\n", i, invPerm[i]);
         SCIPdebugMessage(" -> node is feasible (pair was fixed to (1,0) and every earlier pair is constant).\n");

         break;
      }
      /* cannot apply propagation */
      else
         break;
   }

   return SCIP_OKAY;
}


/** Get the current lp solution from SCIP solution @p sol.
 *
 *  We store the solution in the data of the constraint.
 */
static
SCIP_RETCODE getValues(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_SOL*             sol,                /**< solution (may be 0) */
   const SCIP_CONSDATA*  consdata            /**< constraint data */
   )
{
   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nVars > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );

   SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nVars, consdata->vars, consdata->vals) );

   return SCIP_OKAY;
}


/** add symresack cover inequality */
static
SCIP_RETCODE addSymresackInequality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   const SCIP_CONSDATA*  consdata,           /**< constraint data */
   int*                  coeffs,             /**< coefficient vector of inequality to be added */
   SCIP_Real             rhs,                /**< right-hand side of inequality to be added */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   SCIP_ROW* row;
   unsigned int i;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nVars > 0 );
   assert( consdata->vars != NULL );
   assert( coeffs != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "symresack", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
   for (i = 0; i < (consdata->nVars); ++i)
   {
      if ( coeffs[i] == 1 )
         SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[i], 1.0) );
      else if ( coeffs[i] == -1 )
         SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[i], -1.0) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );
   SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, infeasible) ); /* note that the solution is not used in SCIPaddCut() */
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


/** separate symresack cover inequalities
 *
 *  We currently do NOT enter cuts into the pool.
 */
static
SCIP_RETCODE separateSymresackCovers(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   const SCIP_CONSDATA*  consdata,           /**< constraint data */
   unsigned int*         nGen,               /**< number of separated covers */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   unsigned int nVars;
   SCIP_Real* vals;
   unsigned int* perm;
   unsigned int* invPerm;
   SCIP_Real constObjective;
   SCIP_Real* sepaObjective;
   unsigned int i;
   int* tmpSolu;
   SCIP_Real tmpSoluObj;
   int* maxSolu;
   SCIP_Real maxSoluObj;
   unsigned int c;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nVars > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );
   assert( consdata->perm != NULL );
   assert( consdata->invPerm != NULL );
   assert( infeasible != NULL );
   assert( nGen != NULL );

   *infeasible = FALSE;
   *nGen = 0;

   nVars = consdata->nVars;
   vals = consdata->vals;
   perm = consdata->perm;
   invPerm = consdata->invPerm;

   /* initialize objective */
   SCIP_CALL( SCIPallocBufferArray(scip, &sepaObjective, nVars) );
   constObjective = 1.0; /* constant part of separation objective */

   for (i = 0; i < nVars; ++i)
   {
      if ( i < perm[i] )
      {
         sepaObjective[i] = vals[i];
         constObjective -= vals[i];
      }
      else
         sepaObjective[i] = vals[i] - 1.0;
   }

   /* allocate memory for temporary and global solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpSolu, nVars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxSolu, nVars) );
   tmpSoluObj = 0.0;
   maxSoluObj = 0.0;

   /* start separation procedure by iterating over critical rows */
   for (c = 0; c < nVars; ++c)
   {
      /* there are no fixed points */
      assert( perm[c] != c );

      /* initialize temporary solution */
      for (i = 0; i < nVars; ++i)
         tmpSolu[i] = 2;
      tmpSoluObj = 0.0;

      /* perform fixings implied by the critical row */
      tmpSolu[c] = 0;
      assert( invPerm[c] < nVars );

      tmpSolu[invPerm[c]] = 1;
      tmpSoluObj += sepaObjective[invPerm[c]];

      /* perform 1-fixings */
      i = invPerm[c];
      while ( i < c )
      {
         i = invPerm[i];
         tmpSolu[i] = 1;
         tmpSoluObj += sepaObjective[i];
      }

      /* row c cannot be critical */
      if ( i == c )
         continue;

      assert( tmpSolu[c] == 0 );

      /* perform 0-fixing */
      i = perm[c];
      while ( i < c )
      {
         tmpSolu[i] = 0;
         i = perm[i];
      }

      /* iterate over rows above the critical row */
      for (i = 0; i < c; ++i)
      {
         unsigned int j;
         SCIP_Real objImpact = 0.0;

         /* skip already fixed entries */
         if ( tmpSolu[i] != 2 )
            continue;

         /* Check effect of fixing entry i to 1 and apply all implied fixing to other entries.
          *
          * Observe: Experiments indicate that entries are more often fixed to 1 than to 0.
          * For this reason, we apply the 1-fixings directly. If it turs out that the 1-fixings
          * have a negative impact on the objective, we undo these fixings afterwards and apply
          * 0-fixings instead. */
         /* check fixings in invPerm direction */
         j = i;
         do
         {
            assert( tmpSolu[j] == 2 );
            tmpSolu[j] = 1;
            objImpact += sepaObjective[j];
            j = invPerm[j];
         }
         while ( j < c && j != i );

         /* if we do not detect a cycle */
         if ( j != i )
         {
            /* fix entry j since this is not done in the above do-while loop */
            assert( tmpSolu[j] == 2 );
            tmpSolu[j] = 1;
            objImpact += sepaObjective[j];

            /* check fixings in perm direction */
            j = perm[i];
            while ( j < c )
            {
               assert( j != i );
               assert( tmpSolu[j] == 2 );
               tmpSolu[j] = 1;
               objImpact += sepaObjective[j];
               j = perm[j];
            }

            assert( j != c );
         }

         /* if fixing entry i has a positive impact -> keep above fixings of entries to 1 */
         /* otherwise -> reset entries to 0 */
         if ( SCIPisEfficacious(scip, objImpact) )
            tmpSoluObj += objImpact;
         else
         {
            j = i;
            do
            {
               assert( tmpSolu[j] == 1 );
               tmpSolu[j] = 0;
               j = invPerm[j];
            }
            while ( j < c && j != i );

            /* if we do not detect a cycle */
            if ( j != i )
            {
               /* fix entry j since this is not done in the above do-while loop */
               assert( tmpSolu[j] == 1 );
               tmpSolu[j] = 0;

               /* check fixings in perm direction */
               j = perm[i];
               while ( j < c )
               {
                  assert( j != i );
                  assert( tmpSolu[j] == 1 );
                  tmpSolu[j] = 0;
                  j = perm[j];
               }

               assert( j != c );
            }
         }
      }

      /* iterate over unfixed entries below the critical row */
      for (i = c + 1; i < nVars; ++i)
      {
         /* skip already fixed entries */
         if ( tmpSolu[i] != 2 )
            continue;

         if ( SCIPisEfficacious(scip, sepaObjective[i]) )
         {
            assert( tmpSolu[i] == 2 );
            tmpSolu[i] = 1;
            tmpSoluObj += sepaObjective[i];
         }
         else
         {
            assert( tmpSolu[i] == 2 );
            tmpSolu[i] = 0;
         }
      }

      /* check whether we have found a better solution which has positive separation objective*/
      if ( SCIPisEfficacious(scip, tmpSoluObj + constObjective - maxSoluObj) )
      {
         assert( SCIPisEfficacious(scip, tmpSoluObj + constObjective) );
         for (i = 0; i < nVars; ++i)
            maxSolu[i] = tmpSolu[i];
         maxSoluObj = tmpSoluObj + constObjective;
      }
   }

   /* Check whether the separation objective is positive, i.e., a violated cover was found. */
   if ( SCIPisEfficacious(scip, maxSoluObj) )
   {
      SCIP_Real rhs = -1.0;
      SCIP_Real lhs = 0.0;

      for (i = 0; i < nVars; ++i)
      {
         if ( i < perm[i] )
         {
            maxSolu[i] = maxSolu[i] - 1;
            lhs += vals[i] * maxSolu[i];
         }
         else
         {
            lhs += vals[i] * maxSolu[i];
            rhs += maxSolu[i];
         }
      }

      assert( SCIPisGT(scip, lhs, rhs) );

      /* add cover inequality */
      SCIP_CALL( addSymresackInequality(scip, cons, consdata, maxSolu, rhs, infeasible) );

      if ( ! *infeasible )
         ++(*nGen);
   }

   SCIPfreeBufferArrayNull(scip, &maxSolu);
   SCIPfreeBufferArrayNull(scip, &tmpSolu);
   SCIPfreeBufferArrayNull(scip, &sepaObjective);

   return SCIP_OKAY;
}


/** Upgrade symresack constraints to orbisacks */
static
SCIP_RETCODE orbisackUpgrade(
   SCIP*                 scip,               /**< SCIP pointer */
   unsigned int*         perm,               /**< permutation */
   unsigned int          nVars,              /**< size of perm array */
   SCIP_VAR**            inputVars,          /**< permuted variables array */
   SCIP_Bool*            success,            /**< whether constraint was upgraded */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   unsigned int i;
   unsigned int nrows = 0;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   unsigned int maxVars;

   assert( scip != NULL );
   assert( perm != NULL );
   assert( nVars > 0 );
   assert( inputVars != NULL );
   assert( success != NULL );

   *success = TRUE;

   maxVars = nVars / 2;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars1, maxVars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars2, maxVars) );

   /* check whether permutation is a composition of 2-cycles */
   for (i = 0; i < nVars; ++i)
   {
      if ( perm[perm[i]] != i || ! SCIPvarIsBinary(inputVars[i]) )
      {
         *success = FALSE;
         break;
      }

      if ( perm[i] > i )
      {
         vars1[nrows] = inputVars[i];
         vars2[nrows++] = inputVars[perm[i]];

         assert( nrows <= maxVars );
      }
   }

   /* if permutation can be upgraded to an orbisack */
   if ( *success )
   {
      SCIP_CALL( SCIPcreateConsOrbisack(scip, cons, "orbisack", vars1, vars2, nrows, FALSE,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   SCIPfreeBufferArray(scip, &vars2);
   SCIPfreeBufferArray(scip, &vars1);

   return SCIP_OKAY;
}


/*--------------------------------------------------------------------------------------------
 *--------------------------------- SCIP functions -------------------------------------------
 *--------------------------------------------------------------------------------------------*/

#if 0
/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySymresack)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( valid != NULL );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSymresack(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}
#endif


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSymresack)
{
   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( consdata != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** frees constraint handler */
static
SCIP_DECL_CONSFREE(consFreeSymresack)
{   /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSymresack)
{
   SCIP_CONSDATA* sourcedata;
   unsigned int nVars;
   SCIP_CONSDATA* consdata = NULL;
   unsigned int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   SCIPdebugMessage("Transforming constraint.\n");

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);
   assert( sourcedata->nVars != 0 );
   assert( sourcedata->vars != NULL );
   assert( sourcedata->vals != NULL );
   assert( sourcedata->perm != NULL );
   assert( sourcedata->invPerm != NULL );
   if ( sourcedata->ppUpgrade )
   {
      assert( sourcedata->nCycles != 0 );
      assert( sourcedata->cycleDecomposition != NULL );
      for (i = 0; i < sourcedata->nCycles; ++i)
      {
         assert( sourcedata->cycleDecomposition[i] != NULL );
         assert( sourcedata->cycleDecomposition[i][0] != 0 );
      }
   }

   /* create transformed constraint data (copy data where necessary) */
   nVars = sourcedata->nVars;

   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nVars = nVars;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, nVars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vals, nVars) );
   SCIP_CALL( SCIPgetTransformedVars(scip, nVars, sourcedata->vars, consdata->vars) );
   for (i = 0; i < nVars; ++i)
   {
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
   }

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->perm, sourcedata->perm, nVars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->invPerm, sourcedata->invPerm, nVars) );

   consdata->ppUpgrade = sourcedata->ppUpgrade;

   if ( sourcedata->ppUpgrade )
   {
      consdata->nCycles = sourcedata->nCycles;
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->cycleDecomposition, sourcedata->cycleDecomposition, sourcedata->nCycles) );
      for (i = 0; i < sourcedata->nCycles; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->cycleDecomposition[i], sourcedata->cycleDecomposition[i], nVars + 1) );
      }
   }

   /* create transformed constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpSymresack)
{
   int c;
#if ( SCIP_VERSION < 321 || ( SCIP_VERSION == 321 && SCIP_SUBVERSION == 0 ) )
   SCIP_Bool infeasible = FALSE;
#else
   assert( infeasible != NULL );
   *infeasible = FALSE;
#endif

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != 0 );

      SCIPdebugMessage("Generating initial symresack cut for constraint <%s> ...\n", SCIPconsGetName(conss[c]));

#if ( SCIP_VERSION < 321 || ( SCIP_VERSION == 321 && SCIP_SUBVERSION == 0 ) )
      SCIP_CALL( initLP(scip, conss[c], &infeasible) );
      if ( infeasible )
         break;
#else
      SCIP_CALL( initLP(scip, conss[c], infeasible) );
      if ( *infeasible )
         break;
#endif
   }
   SCIPdebugMessage("Generated initial symresack cuts.\n");

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpSymresack)
{
   SCIP_CONSDATA* consdata;
   int c;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Separation method for symresack constraints\n");

   *result = SCIP_DIDNOTRUN;

   /* if solution is not integer */
   if ( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      unsigned int nGen = 0;

      /* get data of constraint */
      assert( conss[c] != 0 );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      SCIP_CALL( getValues(scip, NULL, consdata) );

      SCIPdebugMessage("Separating symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, &nGen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if ( nGen > 0 )
         *result = SCIP_SEPARATED;

      if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solution */
static
SCIP_DECL_CONSSEPASOL(consSepasolSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Separation method for symresack constraints\n");

   *result = SCIP_DIDNOTRUN;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      unsigned int nGen = 0;

      /* get data of constraint */
      assert( conss[c] != 0 );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      SCIP_CALL( getValues(scip, sol, consdata) );

      SCIPdebugMessage("Separating symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, &nGen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if ( nGen > 0 )
         *result = SCIP_SEPARATED;

      if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions.
 *
 *  To check feasibility, we separate cover inequalities.
 *
 *  @pre It is assumed that the solution is integral (this can be ensured by appropriate priorities).
 */
static
SCIP_DECL_CONSENFOLP(consEnfolpSymresack)
{
   SCIP_CONSDATA* consdata;
   int c;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Enforcing method for symresack constraints (lp solutions) ...\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      unsigned int nGen = 0;

      /* get data of constraint */
      assert( conss[c] != 0 );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      SCIP_CALL( getValues(scip, NULL, consdata) );

      SCIPdebugMessage("Enforcing symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, &nGen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* SCIPdebugMessage("Generated symresack inequalities for <%s>: %u\n", SCIPconsGetName(conss[c]), nGen); */

      if ( nGen > 0 )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSymresack)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Enforcing method for symresack constraints (pseudo solutions) ...\n");

   *result = SCIP_FEASIBLE;

   if ( objinfeasible || solinfeasible )
      return SCIP_OKAY;

   if ( solinfeasible ) /*lint !e774*/
      return SCIP_OKAY;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool terminated = FALSE;
      SCIP_CONSDATA* consdata;
      unsigned int nVars;
      SCIP_VAR** vars;
      unsigned int* invPerm;
      unsigned int i;
      unsigned int* solu;
      SCIP_Real val;

      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL);
      assert( consdata->nVars > 0 );
      assert( consdata->vars != NULL );
      assert( consdata->invPerm != NULL );
      nVars = consdata->nVars;
      vars = consdata->vars;
      invPerm = consdata->invPerm;

      /* determine solution */
      SCIP_CALL( SCIPallocBufferArray(scip, &solu, (int) nVars) );

      for (i = 0; i < nVars; ++i)
      {
         /* there are no fixed points */
         assert( invPerm[i] != i );

         /* get value of variables */
         val = SCIPgetSolVal(scip, NULL, vars[i]);
         assert( SCIPisFeasIntegral(scip, val) );

         /* if variable is fixed to 1 -> solu[i] = 1, else = 0 */
         if ( val > 0.5 )
            solu[i] = 1;
         else
            solu[i] = 0;
      }

      /* check whether solution is lexicographically not smaller than its permutation */
      for (i = 0; i < nVars; ++i)
      {
         /* there are no fixed points */
         assert( invPerm[i] != i );

         /* if pair (i,invPerm[i]) is constant */
         if ( solu[i] == solu[invPerm[i]] )
            continue;

         /* if first non-constant pair is (1,0): feasible */
         if ( solu[i] == 1 )
            break;
         else /* infeasible */
         {
            SCIPdebugMessage("Solution is infeasible.\n");
            *result = SCIP_INFEASIBLE;
            terminated = TRUE;
            break;
         }
      }

      /* free buffers */
      SCIPfreeBufferArrayNull(scip, &solu);

      if ( terminated )
         break;
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckSymresack)
{   /*lint --e{715}*/
   int c;
   SCIP_CONSDATA* consdata;
   SCIP_Bool terminated = FALSE;
   unsigned int nVars;
   SCIP_VAR** vars;
   unsigned int* invPerm;
   unsigned int i;
   SCIP_Real solVal;
   int val1;
   int val2;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   /* loop through constraints */
   for (c = 0; c < nconss && ! terminated; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL);
      assert( consdata->nVars > 0 );
      assert( consdata->vars != NULL );
      assert( consdata->invPerm != NULL );

      SCIPdebugMessage("Check method for symresack constraint <%s> (%u rows) ...\n", SCIPconsGetName(conss[c]), consdata->nVars);

      nVars = consdata->nVars;
      vars = consdata->vars;
      invPerm = consdata->invPerm;

      /* detect first non-constant pair of variables */
      for (i = 0; i < nVars; ++i)
      {
         /* there are no fixed points */
         assert( invPerm[i] != i );

         /* get value of variable i and its inverse */
         solVal = SCIPgetSolVal(scip, sol, vars[i]);
         assert( SCIPisFeasIntegral(scip, solVal) );
         if ( solVal > 0.5 )
            val1 = 1;
         else
            val1 = 0;

         solVal = SCIPgetSolVal(scip, sol, vars[invPerm[i]]);
         assert( SCIPisFeasIntegral(scip, solVal) );
         if ( solVal > 0.5 )
            val2 = 1;
         else
            val2 = 0;

         /* if we detected a constant pair */
         if ( val1 == val2 )
            continue;
         /* pair is (1,0) --> lexicographically maximal */
         else if ( val1 > val2 )
            break;

         /* pair is (0,1) --> solution is infeasible */
         assert( val2 > val1 );
         SCIPdebugMessage("Solution is infeasible.\n");
         *result = SCIP_INFEASIBLE;
         terminated = TRUE;

         if ( printreason )
            SCIPinfoMessage(scip, NULL, "First non-constant pair (%u, %u) of variables has pattern (0,1).\n", i, invPerm[i]);

         break;
      }
   }

   if ( ! terminated )
      SCIPdebugMessage("Solution is feasible.\n");

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSymresack)
{
   int c;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMessage("Propagation method of symresack constraint handler.\n");

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_Bool found = FALSE;
      unsigned int nGen = 0;

      assert( conss[c] != 0 );

      SCIP_CALL( propVariables(scip, conss[c], &infeasible, &found, &nGen) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if ( found )
      {
         *result = SCIP_REDUCEDDOM;
          return SCIP_OKAY;
      }

      *result = SCIP_DIDNOTFIND;
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSymresack)
{
   int c;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Presolving method of symresack constraint handler. Propagating symresack inequalities.\n");
   *result = SCIP_DIDNOTRUN;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_Bool found = FALSE;
      unsigned int nGen = 0;

      assert( conss[c] != 0 );
      SCIP_CALL( propVariables(scip, conss[c], &infeasible, &found, &nGen) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if ( nGen > 0 )
      {
         *nfixedvars += (int) nGen;
         *result = SCIP_SUCCESS;
      }

      if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;
   }

   return SCIP_OKAY;
}


/** Propagation resolution for conflict analysis */
static
SCIP_DECL_CONSRESPROP(consRespropSymresack)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   unsigned int nVars;
   unsigned int* invPerm;
   unsigned int i;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIPdebugMessage("Propagation resolution method of symresack constraint handler.\n");

   *result = SCIP_DIDNOTFIND;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nVars > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->invPerm != NULL );

   vars = consdata->vars;
   nVars = consdata->nVars;
   invPerm = consdata->invPerm;

   assert( 0 <= inferinfo && inferinfo < (int) (2 * nVars - 1) );

   /* if first part of variable pair was fixed to 0 */
   if ( inferinfo < (int) nVars )
   {
      assert( vars[invPerm[inferinfo]] == infervar );
      assert( SCIPvarGetUbAtIndex(vars[invPerm[inferinfo]], bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(vars[invPerm[inferinfo]], bdchgidx, TRUE) < 0.5 );

      if ( SCIPvarGetUbAtIndex(vars[invPerm[inferinfo]], bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(vars[invPerm[inferinfo]], bdchgidx, TRUE) < 0.5 )
      {
         SCIPdebugMessage(" -> reason for setting x[%u] = 0 was fixing x[%u] to 0 and each pair of binary variables before (%u,%u) which are not fixed points is constant.\n",
            invPerm[inferinfo], inferinfo, inferinfo, invPerm[inferinfo]);

         SCIP_CALL( SCIPaddConflictUb(scip, vars[inferinfo], bdchgidx) );
         for (i = 0; (int) i < inferinfo; ++i)
         {
            /* there are no fixed points */
            assert( invPerm[i] != i );

            SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[invPerm[i]], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[invPerm[i]], bdchgidx) );
         }
      }
   }
   /* if second part of variable pair was fixed to 1 */
   else
   {
      int inferinfo2;

      inferinfo2 = inferinfo - nVars;
      assert( vars[inferinfo2] == infervar );
      assert( SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, FALSE) < 0.5 && SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, TRUE) > 0.5 );

      if ( SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, FALSE) < 0.5 && SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, TRUE) > 0.5 )
      {
         SCIPdebugMessage(" -> reason for setting x[%u] = 1 was fixing x[%u] to 0 and each pair of binary variables before (%u,%u) which are not fixed points is constant.\n",
            inferinfo2, invPerm[inferinfo2], inferinfo2, invPerm[inferinfo2]);

         SCIP_CALL( SCIPaddConflictLb(scip, vars[invPerm[inferinfo2]], bdchgidx) );
         for (i = 0; (int) i < inferinfo2; ++i)
         {
            /* there are no fixed points */
            assert( invPerm[i] != i );

            SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[invPerm[i]], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[invPerm[i]], bdchgidx) );
         }
      }
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** lock variables
 *
 *  We assume we have only one global (void) constraint and lock all binary variables
 *  which do not correspond to fixed points of the permutation.
 *
 * - Symresack constraints may get violated if the variables with a negative coefficient
 *   in the Friedman inequality are rounded down, we therefor call
 *   SCIPaddVarLocks(..., nlockspos, nlocksneg).
 * - Symresack constraints may get violated if the variables with a positive coefficient
 *   in the Friedman inequality are rounded up, we therefor call
 *   SCIPaddVarLocks(..., nlocksneg, nlockspo ).
 */
static
SCIP_DECL_CONSLOCK(consLockSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   unsigned int nVars;
   unsigned int i;
   SCIP_VAR** vars;
   unsigned int* perm;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMessage("Locking method for symresack constraint handler.\n");

   /* get data of original constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nVars > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->perm != NULL );

   nVars = consdata->nVars;
   vars = consdata->vars;
   perm = consdata->perm;

   for (i = 0; i < nVars; ++i)
   {
      /* there are no fixed points */
      assert( perm[i] != i );

      if ( perm[i] > i )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, vars[i], nlockspos, nlocksneg) );
      }
      else
      {
         assert( perm[i] != i );
         SCIP_CALL( SCIPaddVarLocks(scip, vars[i], nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler
 *
 *  The constraint handler should output a representation of the constraint into the given text file.
 */
static
SCIP_DECL_CONSPRINT(consPrintSymresack)
{
   /*lint --e{715}*/
   return SCIP_OKAY;
}


#if 0
/** constraint copying method of constraint handler
 *
 *  The constraint handler can provide a copy method which copy a constraint from one SCIP data structure into another
 *  SCIP data structure.
 */
static
SCIP_DECL_CONSCOPY(consCopySymresack)
{
   SCIP_CONSDATA* sourcedata;
   unsigned int nVars;
   SCIP_VAR** sourceVars;
   SCIP_VAR** vars = NULL;
   unsigned int i;
   SCIP_Bool success;

   /*lint --e{715}*/
   assert( scip != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( sourcescip != NULL );
   assert( sourcecons != NULL );
   assert( varmap != NULL );
   assert( valid != NULL );

   *valid = TRUE;

   SCIPdebugMessage("Copying method for symresack constraint handler.\n");

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->vars != NULL );
   assert( sourcedata->nVars  > 0 );

   /* copy constraint data */
   nVars = sourcedata->nVars;
   sourceVars = sourcedata->vars;

   /* separately allocate space to account for unsuccessful copying */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, (int) nVars) );

   /* get copies */
   for (i = 0; i < nVars; ++i)
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourceVars[i], &vars[i], varmap, consmap, global, valid) );
      assert( *valid );
      assert( vars[i] != NULL );
   }

   /* only create target constraint if all variables could be copied */
   if ( *valid )
   {
      /* create copied constraint */
      if ( name == 0 )
         name = SCIPconsGetName(sourcecons);

      SCIP_CALL( SCIPcreateConsSymresack(scip, cons, name, sourcedata->perm, vars, nVars, &success,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic,
            removable, stickingatnode) );
   }

   SCIPfreeBufferArrayNull(scip, &vars);

   return SCIP_OKAY;
}
#endif


/** creates the handler for symresack constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSymresack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSymresack, consEnfopsSymresack, consCheckSymresack, consLockSymresack,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
#if 0
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySymresack, consCopySymresack) );
#endif
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSymresack) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSymresack) );
#if ( SCIP_VERSION >= 320 )
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSymresack, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
#else
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSymresack, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
#endif
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSymresack) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSymresack, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSymresack) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSymresack, consSepasolSymresack, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSymresack) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSymresack) );

   /* get enforcing setting */
   SCIP_CALL( SCIPaddBoolParam(scip, "cons/symresack/enforcing", "Enforce symresack constraints?", &conshdlrdata->symresackEnforcing, TRUE,
         DEFAULT_ENFORCING, NULL, NULL) );

   /* get check setting */
   SCIP_CALL( SCIPaddBoolParam(scip, "cons/symresack/check", "Check symresack constraints?", &conshdlrdata->symresackCheck, TRUE,
         DEFAULT_CHECK, NULL, NULL) );

   /* whether we allow upgrading to orbisack constraints*/
   SCIP_CALL( SCIPaddBoolParam(scip, "cons/symresack/upgrade", "Upgrade symresack constraints to orbisack constraints?", &conshdlrdata->symresackUpgrade, TRUE,
         DEFAULT_UPGRADE, NULL, NULL) );

   /* whether we allow upgrading to packing/partioning symresack constraints*/
   SCIP_CALL( SCIPaddBoolParam(scip, "cons/symresack/ppsymresack", "Upgrade symresack constraints to packing/partioning symresacks?", &conshdlrdata->checkPPsymresack, TRUE,
         DEFAULT_PPSYMRESACK, NULL, NULL) );

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates and captures a symresack constraint
 *
 *  In a presolving step, we check whether the permutation acts only on binary points. Otherwise, we eliminate
 *  the non-binary variables from the permutation. If the permutation is the identity (after variable elimination)
 *  the boolean success is set to false.
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   unsigned int*         perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   unsigned int          nVars,              /**< number of variables in problem */
   SCIP_Bool*            success,            /**< whether permutation is acting only on binary points */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_Bool upgrade;

   assert( cons != NULL );
   assert( success != NULL );
   *success = FALSE;

   /* find the symresack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("Symresack constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   assert( nVars > 0 );

   /* get enforce and check settings */
   SCIP_CALL( SCIPgetBoolParam(scip, "cons/symresack/enforcing", &enforce) );
   SCIP_CALL( SCIPgetBoolParam(scip, "cons/symresack/check", &check) );

   /* check whether constraint can be upgraded to an orbisack constraint */
   SCIP_CALL( SCIPgetBoolParam(scip, "cons/symresack/upgrade", &upgrade) );

   if ( upgrade )
   {
      SCIP_CALL( orbisackUpgrade(scip, perm, nVars, vars, success, cons, initial, separate, enforce, check, propagate,
            local, modifiable, dynamic, removable, stickingatnode) );

      if ( *success )
      {
         SCIPdebugMessage("Upgraded symresack constraint to orbisack constraint.\n");

         return SCIP_OKAY;
      }
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, vars, nVars, perm) );

   /* delete constraint if there are no restrictions, i.e., the input permutation does not act on binary variables */
   if ( consdata->nVars == 0 )
   {
      SCIPfreeBlockMemory(scip, &consdata);
      *cons = NULL;
   }
   else
   {
      /* create constraint */
      SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate && (! consdata->ppUpgrade), enforce, check, propagate,
            local, modifiable, dynamic, removable, stickingatnode) );

      *success = TRUE;
   }

   return SCIP_OKAY;
}


/** creates and captures a symresack constraint
 *  in its most basic variant, i.e., with all constraint flags set to their default values
 *
 *  In a presolving step, we check whether the permutation acts only on binary points. Otherwise, the boolean
 *  success is set to false.
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsBasicSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   unsigned int*         perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   unsigned int          nVars,              /**< number of variables in problem */
   SCIP_Bool*            success             /**< whether permutation is acting only on binary points */
   )
{
   SCIP_CALL( SCIPcreateConsSymresack(scip, cons, name, perm, vars, nVars, success,
         TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
