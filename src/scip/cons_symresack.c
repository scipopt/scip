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
 * The type of constraints of this constraint handler is described in cons_symresack.h.
 *
 * The details of the method implemented here are described in the following papers:
 *
 * Fundamental Domains for Integer Programs with Symmetries@n
 * Eric J. Friedman,@n
 * Combinatorial Optimization, volume 4616 of LNCS, 146-153 (2007)
 *
 * This paper describes an inequality to handle symmetries of a single permutation. This
 * so-called FD-inequality is the basic for the propagation routine of our implementation.
 *
 * Polytopes Associated with Symmetry Handling@n
 * Christopher Hojny and Marc E. Pfetsch,@n
 * (2017), preprint available at http://www.optimization-online.org/DB_HTML/2017/01/5835.html
 *
 * This paper describes an almost linear time separation routine for so-called cove
 * inequalities of symresacks. In our implementation, however, we use a separation routine with
 * quadratic worst case running time.
 *
 * Packing, Partitioning, and Covering Symresacks@n
 * Christopher Hojny,@n
 * (2017), preprint available at http://www.optimization-online.org/DB_HTML/2017/05/5990.html
 *
 * This paper introduces linearly many inequalities with ternary coefficients that suffice to
 * characterize the binary points contained in a packing and partitioning symresack completely.
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
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_EXHAUSTIVE

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
   SCIP_Bool             symresackupgrade;   /**< whether we allow upgrading symresack constraints to orbisack constraints */
   SCIP_Bool             checkppsymresack;   /**< whether we allow upgrading to packing/partitioning symresacks */
};


/** constraint data for symresack constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables */
   int                   nvars;              /**< number of variables */
   SCIP_Real*            vals;               /**< LP-solution for the variables */
   int*                  perm;               /**< permutation associated to the symresack */
   int*                  invperm;            /**< inverse permutation */
   SCIP_Bool             ppupgrade;          /**< whether constraint is upgraded to packing/partitioning symresack */

   /* data for upgraded symresack constraints */
   int                   ncycles;            /**< number of cycles in permutation */
   int**                 cycledecomposition; /**< cycle decomposition */
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
   int nvars;
   int i;

   assert( consdata != NULL );
   assert( *consdata != NULL );

   nvars = (*consdata)->nvars;

   if ( nvars == 0 )
   {
      SCIPfreeBlockMemory(scip, consdata);

      return SCIP_OKAY;
   }

   if ( (*consdata)->ppupgrade )
   {
      for (i = 0; i < (*consdata)->ncycles; ++i)
      {
         SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cycledecomposition[i]), nvars + 1);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cycledecomposition), (*consdata)->ncycles);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vals), nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->invperm), nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->perm), nvars);

   for (i = 0; i < nvars; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars), nvars);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** check whether constraint can be upgraded to packing/partitioning symresack */
static
SCIP_RETCODE packingUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables affected by permutation */
   int                   n,                  /**< length of permutation */
   SCIP_Bool*            upgrade             /**< pointer to store whether upgrade was successful */
   )
{
   SCIP_Bool* covered;
   SCIP_Bool descent;
   SCIP_CONSHDLR* setppcconshdlr;
   SCIP_CONS** setppcconss;
   SCIP_VAR* var;
   SCIP_Bool terminated = FALSE;
   int** cycledecomposition;
   int* indicesincycle;
   int nsetppcconss;
   int curcycle;
   int maxcyclelength;
   int ncycles = 0;
   int c;
   int i;
   int j;

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

      ++ncycles;
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
   assert( ncycles <= n / 2 );

   /* each cycle is monotone; check for packing/partitioning type */
   for (i = 0; i < n; ++i)
      covered[i] = FALSE;

   /* compute cycle decomposition: row i stores in entry 0 the length of the cycle,
    * the remaining entries are the coordinates in the cycle */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cycledecomposition, ncycles) );
   for (i = 0; i < ncycles; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cycledecomposition[i], n + 1) );
   }

   curcycle = 0;
   maxcyclelength = 0;
   for (i = 0; i < n; ++i)
   {
      int cyclelength = 0;

      /* skip checked indices */
      if ( covered[i] )
         continue;

      j = i;
      do
      {
         covered[j] = TRUE;
         cycledecomposition[curcycle][++cyclelength] = j;
         j = perm[j];
      }
      while ( j != i );

      cycledecomposition[curcycle][0] = cyclelength;
      ++curcycle;

      if ( maxcyclelength < cyclelength )
         maxcyclelength = cyclelength;
   }

   /* permutation can be upgraded -> check whether the symresack is of packing/partitioning type */
   setppcconshdlr = SCIPfindConshdlr(scip, "setppc");
   if ( setppcconshdlr == NULL )
   {
      SCIPerrorMessage("Setppc constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   setppcconss = SCIPconshdlrGetConss(setppcconshdlr);
   nsetppcconss = SCIPconshdlrGetNConss(setppcconshdlr);

   /* Check whether each cycle of the symresack is contained in a set packing/partitioning constraint.
    * To this end, we have to guarantee that all affected variables are not negated since permutations
    * are given w.r.t. original variables. */
   *upgrade = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &indicesincycle, maxcyclelength) );

   for (i = 0; i < ncycles && *upgrade && ! terminated; ++i)
   {
      int cyclelength;

      /* get indices of variables in current cycle */
      for (j = 0; j < cycledecomposition[i][0]; ++ j)
      {
         var = vars[cycledecomposition[i][j + 1]];

         if ( SCIPvarIsNegated(var) )
         {
            terminated = TRUE;
            break;
         }

         indicesincycle[j] = SCIPvarGetProbindex(var);
      }

      cyclelength = cycledecomposition[i][0];

      /* iterate over constraints */
      for (c = 0; c < nsetppcconss; ++c)
      {
         int nsetppcvars;
         SCIP_VAR** setppcvars;
         int varidx;
         int nfound = 0;
         int k;

         /* check type */
         if ( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_COVERING )
            continue;
         assert( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PARTITIONING || SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PACKING );

         /* get set packing/partitioning variables */
         nsetppcvars = SCIPgetNVarsSetppc(scip, setppcconss[c]);
         assert( nsetppcvars > 0 );

         setppcvars = SCIPgetVarsSetppc(scip, setppcconss[c]);
         assert( setppcvars != NULL );

         /* check whether all variables of the cycle are contained in setppc constraint */
         for (j = 0; j < nsetppcvars && nfound < cyclelength; ++j)
         {
            var = setppcvars[j];

            if ( SCIPvarIsNegated(var) )
               continue;

            varidx = SCIPvarGetProbindex(var);

            for (k = 0; k < cyclelength; ++k)
            {
               if ( varidx == indicesincycle[k] )
               {
                  ++nfound;
                  break;
               }
            }
         }

         if ( nfound == cyclelength )
            break;
      }

      /* row is not contained in a set packing/partitioning constraint */
      if ( c >= nsetppcconss )
         *upgrade = FALSE;
   }

   if ( *upgrade )
   {
      (*consdata)->ncycles = ncycles;
      (*consdata)->cycledecomposition = cycledecomposition;

      SCIPfreeBufferArray(scip, &indicesincycle);
      SCIPfreeBufferArray(scip, &covered);
   }
   else
   {
      SCIPfreeBufferArray(scip, &indicesincycle);
      for (i = 0; i < ncycles; ++i)
      {
         SCIPfreeBlockMemoryArray(scip, &cycledecomposition[i], n + 1);
      }
      SCIPfreeBlockMemoryArray(scip, &cycledecomposition, ncycles);
      SCIPfreeBufferArray(scip, &covered);
   }

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
   SCIP_VAR*const*       inputvars,          /**< input variables of the constraint handler */
   int                   inputnvars,         /**< input number of variables of the constraint handler*/
   int*                  inputperm           /**< input permutation of the constraint handler */
   )
{
   SCIP_VAR** vars;
   SCIP_Bool upgrade;
   int* indexcorrection;
   int* invperm;
   int* perm;
   int naffectedvariables;
   int i;
   int j = 0;

   assert( consdata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   /* count the number of binary variables which are affected by the permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &indexcorrection, inputnvars) );
   indexcorrection[0] = -1;
   for (i = 0; i < inputnvars; ++i)
   {
      if ( inputperm[i] != i && SCIPvarIsBinary(inputvars[i]) )
      {
         if ( i == 0 )
            indexcorrection[i] = 0;
         else
            indexcorrection[i] = indexcorrection[i - 1] + 1;
      }
      else
      {
         if ( i > 0 )
            indexcorrection[i] = indexcorrection[i - 1];
      }
   }
   naffectedvariables = indexcorrection[inputnvars - 1] + 1;

   (*consdata)->nvars = naffectedvariables;

   /* Stop if we detect that the permutation fixes each binary point. */
   if ( naffectedvariables == 0 )
   {
      SCIPfreeBufferArrayNull(scip, &indexcorrection);
      return SCIP_OKAY;
   }

   /* remove fixed points from permutation representation */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, naffectedvariables) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &perm, naffectedvariables) );
   for (i = 0; i < inputnvars; ++i)
   {
      if ( i == 0 )
      {
         if ( indexcorrection[i] > -1 )
         {
            vars[j] = inputvars[i];
            perm[j++] = indexcorrection[inputperm[i]];
         }
      }
      else
      {
         if ( indexcorrection[i] > indexcorrection[i - 1] )
         {
            vars[j] = inputvars[i];
            perm[j++] = indexcorrection[inputperm[i]];
         }
      }
   }
   (*consdata)->vars = vars;
   (*consdata)->perm = perm;

   for (i = 0; i < naffectedvariables; ++i)
   {
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[i]) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &invperm, naffectedvariables) );
   for (i = 0; i < naffectedvariables; ++i)
      invperm[perm[i]] = i;
   (*consdata)->invperm = invperm;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*consdata)->vals, naffectedvariables) );

   SCIPfreeBufferArrayNull(scip, &indexcorrection);

   /* check whether an upgrade to packing/partitioning symresacks is possible */
   upgrade = FALSE;
   SCIP_CALL( SCIPgetBoolParam(scip, "cons/symresack/ppsymresack", &upgrade) );

   if ( upgrade )
   {
      upgrade = FALSE;
      SCIP_CALL( packingUpgrade(scip, consdata, perm, vars, naffectedvariables, &upgrade) );
   }

   (*consdata)->ppupgrade = upgrade;

   return SCIP_OKAY;
}


/** generate initial LP cut
 *
 *  We generate the ordering inequality for the pair \f$(1, \gamma^{-1}(1))\f$, i.e.,
 *  the inequality \f$-x_{1} + x_{\gamma^{-1}(1)} \leq 0\f$. This inequality is valid,
 *  because we guaranteed in a preprocessing step that all variables are binary.
 *
 *  Furthermore, we add facet inequalities of packing/partitioning symresacks if
 *  we deal with packing/partitioning symresacks.
 */
static
SCIP_RETCODE initLP(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_ROW* row;
   SCIP_VAR** varsincons;
   SCIP_Real* coeffs;
   int** cycledecomposition;
   int nvarsincons;
   int nvarsincycle;
   int firstelemincycle;
   int nvars;
   int ncycles;
   int i;
   int j;
   int k;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != 0 );

   nvars = consdata->nvars;
   vars = consdata->vars;

   /* avoid stupid problems */
   if ( nvars <= 1 )
      return SCIP_OKAY;

   /* there are no fixed points */
   assert( consdata->invperm[0] != 0 );

   /* add ordering inequality */
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "symresack_init", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, vars[0], -1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, vars[consdata->invperm[0]], 1.0) );

   SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, infeasible) );

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   /* check whether we have a packing/partioning symresack */
   if ( consdata->ppupgrade && ! *infeasible )
   {
      ncycles = consdata->ncycles;
      cycledecomposition = consdata->cycledecomposition;

      SCIP_CALL( SCIPallocBufferArray(scip, &varsincons, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &coeffs, nvars) );

      coeffs[0] = 1.0;

      /* add packing/partitioning symresack constraints */
      for (i = 0; i < ncycles; ++i)
      {
         assert( cycledecomposition[i][0] > 0 );

         nvarsincycle = cycledecomposition[i][0];
         varsincons[0] = vars[cycledecomposition[i][nvarsincycle]];
         firstelemincycle = cycledecomposition[i][1];

         assert( firstelemincycle == consdata->perm[cycledecomposition[i][nvarsincycle]] );

         nvarsincons = 1;

         /* add variables of other cycles to the constraint */
         for (j = 0; j < i; ++j)
         {
            nvarsincycle = cycledecomposition[j][0];
            for (k = 1; k <= nvarsincycle; ++k)
            {
               if ( cycledecomposition[j][k] < firstelemincycle )
               {
                  varsincons[nvarsincons] = vars[cycledecomposition[j][k]];
                  coeffs[nvarsincons++] = -1.0;
               }
               else
                  continue;
            }
         }

         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "ppSymresack", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPaddVarsToRow(scip, row, nvarsincons, varsincons, coeffs) );

         SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, infeasible) );
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }

      SCIPfreeBufferArray(scip, &coeffs);
      SCIPfreeBufferArray(scip, &varsincons);
   }

   return SCIP_OKAY;
}


/** propagation */
static
SCIP_RETCODE propVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_Bool*            infeasible,         /**< pointer to store whether it was detected that the node is infeasible */
   int*                  ngen                /**< pointer to store number of generated bound strengthenings */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool tightened;
   SCIP_VAR** vars;
   SCIP_VAR* var2;
   SCIP_VAR* var;
   int* invperm;
   int nvars;
   int r;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   SCIPdebugMessage("Propagating variables of constraint <%s>.\n", SCIPconsGetName(cons));

   *ngen = 0;
   *infeasible = FALSE;

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nvars != 0 );
   assert( consdata->invperm != NULL );

   vars = consdata->vars;
   nvars = consdata->nvars;
   invperm = consdata->invperm;

   /* avoid trivial problems */
   if ( nvars < 2 )
      return SCIP_OKAY;

   /* loop through all variables */
   for (i = 0; i < nvars; ++i)
   {
      /* there are no fixed points */
      assert( invperm[i] != i );

      /* get variables of first and second column */
      var = vars[i];
      var2 = vars[invperm[i]];
      assert( var != NULL );
      assert( var2 != NULL );

      /* if first part of variable pair fixed to 0 and second part is fixed to 1 */
      if ( ISFIXED0(var) && ISFIXED1(var2) )
      {
         SCIPdebugMessage("Check variable pair (%d,%d).\n", i, invperm[i]);

         SCIPdebugMessage(" -> node infeasible (pair was fixed to (0,1) but there was no pair of type (1,0) before).\n");

         /* perform conflict analysis */
         if ( SCIPisConflictAnalysisApplicable(scip) )
         {
            SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

            for (r = 0; r <= i; ++r)
            {
               /* there are no fixed points */
               assert( invperm[r] != r );

               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[r]) );
               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[invperm[r]]) );
            }

            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }

         *infeasible = TRUE;
         break;
      }
      /* if first part of the variable pair is fixed to 0 and the second part is free --> fix second part to 0 */
      else if ( ISFIXED0(var) && ( ! ISFIXED0(var2) ) )
      {
         assert( SCIPvarGetUbLocal(var) < 0.5 );
         assert( SCIPvarGetLbLocal(var2) < 0.5 );
         assert( SCIPvarGetUbLocal(var2) > 0.5 );

         SCIPdebugMessage("Check variable pair (%d,%d).\n", i, invperm[i]);

         assert( SCIPvarGetLbLocal(var2) < 0.5 );
         SCIP_CALL( SCIPinferVarUbCons(scip, var2, 0.0, cons, i, FALSE, infeasible, &tightened) ); /*lint !e713*/
         assert( ! *infeasible );

         if ( tightened )
            ++(*ngen);
      }
      /* if second part of the variable pair is fixed to 1 and the first part is free --> fix first part to 1 */
      else if ( ( ! ISFIXED1(var) ) && ISFIXED1(var2) )
      {
         assert( SCIPvarGetLbLocal(var) < 0.5 );
         assert( SCIPvarGetUbLocal(var) > 0.5 );
         assert( SCIPvarGetLbLocal(var2) > 0.5 );

         SCIPdebugMessage("Check variable pair (%d,%d).\n", i, invperm[i]);

         assert( SCIPvarGetUbLocal(var) > 0.5 );
         SCIP_CALL( SCIPinferVarLbCons(scip, var, 1.0, cons, i + nvars, FALSE, infeasible, &tightened) ); /*lint !e713*/
         assert( ! *infeasible );

         if ( tightened )
            ++(*ngen);
      }
      /* if solution is lexicographically maximal */
      else if ( ISFIXED1(var) && ISFIXED0(var2) )
      {
         assert( SCIPvarGetLbLocal(var) > 0.5 );
         assert( SCIPvarGetUbLocal(var2) < 0.5 );

         SCIPdebugMessage("Check variable pair (%d,%d).\n", i, invperm[i]);
         SCIPdebugMessage(" -> node is feasible (pair was fixed to (1,0) and every earlier pair is constant).\n");

         break;
      }
      /* cannot apply propagation */
      else
         break;
   }

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
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_ROW* row;
   int i;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nvars > 0 );
   assert( consdata->vars != NULL );
   assert( coeffs != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "symresack", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for (i = 0; i < (consdata->nvars); ++i)
   {
      if ( coeffs[i] == 1 )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[i], 1.0) );
      }
      else if ( coeffs[i] == -1 )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[i], -1.0) );
      }
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
   int*                  ngen,               /**< pointer to store the number of separated covers */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_Real* vals;
   SCIP_Real constobjective;
   SCIP_Real* sepaobjective;
   SCIP_Real tmpsoluobj = 0.0;
   SCIP_Real maxsoluobj = 0.0;
   int* tmpsolu;
   int* maxsolu;
   int* invperm;
   int* perm;
   int nvars;
   int c;
   int i;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->nvars > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->vals != NULL );
   assert( consdata->perm != NULL );
   assert( consdata->invperm != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   *infeasible = FALSE;
   *ngen = 0;

   nvars = consdata->nvars;
   vals = consdata->vals;
   perm = consdata->perm;
   invperm = consdata->invperm;

   /* initialize objective */
   SCIP_CALL( SCIPallocBufferArray(scip, &sepaobjective, nvars) );

   constobjective = 1.0; /* constant part of separation objective */
   for (i = 0; i < nvars; ++i)
   {
      if ( i < perm[i] )
      {
         sepaobjective[i] = vals[i];
         constobjective -= vals[i];
      }
      else
         sepaobjective[i] = vals[i] - 1.0;
   }

   /* allocate memory for temporary and global solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpsolu, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxsolu, nvars) );

   /* start separation procedure by iterating over critical rows */
   for (c = 0; c < nvars; ++c)
   {
      /* there are no fixed points */
      assert( perm[c] != c );

      /* initialize temporary solution */
      for (i = 0; i < nvars; ++i)
         tmpsolu[i] = 2;
      tmpsoluobj = 0.0;

      /* perform fixings implied by the critical row */
      tmpsolu[c] = 0;
      assert( invperm[c] < nvars );

      tmpsolu[invperm[c]] = 1;
      tmpsoluobj += sepaobjective[invperm[c]];

      /* perform 1-fixings */
      i = invperm[c];
      while ( i < c )
      {
         i = invperm[i];
         tmpsolu[i] = 1;
         tmpsoluobj += sepaobjective[i];
      }

      /* row c cannot be critical */
      if ( i == c )
         continue;

      assert( tmpsolu[c] == 0 );

      /* perform 0-fixing */
      i = perm[c];
      while ( i < c )
      {
         tmpsolu[i] = 0;
         i = perm[i];
      }

      /* iterate over rows above the critical row */
      for (i = 0; i < c; ++i)
      {
         SCIP_Real objimpact = 0.0;
         int j;

         /* skip already fixed entries */
         if ( tmpsolu[i] != 2 )
            continue;

         /* Check effect of fixing entry i to 1 and apply all implied fixing to other entries.
          *
          * Observe: Experiments indicate that entries are more often fixed to 1 than to 0.
          * For this reason, we apply the 1-fixings directly. If it turns out that the 1-fixings
          * have a negative impact on the objective, we undo these fixings afterwards and apply
          * 0-fixings instead. */

         /* check fixings in invperm direction */
         j = i;
         do
         {
            assert( tmpsolu[j] == 2 );
            tmpsolu[j] = 1;
            objimpact += sepaobjective[j];
            j = invperm[j];
         }
         while ( j < c && j != i );

         /* if we do not detect a cycle */
         if ( j != i )
         {
            /* fix entry j since this is not done in the above do-while loop */
            assert( tmpsolu[j] == 2 );
            tmpsolu[j] = 1;
            objimpact += sepaobjective[j];

            /* check fixings in perm direction */
            j = perm[i];
            while ( j < c )
            {
               assert( j != i );
               assert( tmpsolu[j] == 2 );
               tmpsolu[j] = 1;
               objimpact += sepaobjective[j];
               j = perm[j];
            }

            assert( j != c );
         }

         /* if fixing entry i has a positive impact -> keep above fixings of entries to 1 */
         /* otherwise -> reset entries to 0 */
         if ( SCIPisEfficacious(scip, objimpact) )
            tmpsoluobj += objimpact;
         else
         {
            j = i;
            do
            {
               assert( tmpsolu[j] == 1 );
               tmpsolu[j] = 0;
               j = invperm[j];
            }
            while ( j < c && j != i );

            /* if we do not detect a cycle */
            if ( j != i )
            {
               /* fix entry j since this is not done in the above do-while loop */
               assert( tmpsolu[j] == 1 );
               tmpsolu[j] = 0;

               /* check fixings in perm direction */
               j = perm[i];
               while ( j < c )
               {
                  assert( j != i );
                  assert( tmpsolu[j] == 1 );
                  tmpsolu[j] = 0;
                  j = perm[j];
               }

               assert( j != c );
            }
         }
      }

      /* iterate over unfixed entries below the critical row */
      for (i = c + 1; i < nvars; ++i)
      {
         /* skip already fixed entries */
         if ( tmpsolu[i] != 2 )
            continue;

         if ( SCIPisEfficacious(scip, sepaobjective[i]) )
         {
            assert( tmpsolu[i] == 2 );
            tmpsolu[i] = 1;
            tmpsoluobj += sepaobjective[i];
         }
         else
         {
            assert( tmpsolu[i] == 2 );
            tmpsolu[i] = 0;
         }
      }

      /* check whether we have found a better solution which has positive separation objective*/
      if ( SCIPisEfficacious(scip, tmpsoluobj + constobjective - maxsoluobj) )
      {
         assert( SCIPisEfficacious(scip, tmpsoluobj + constobjective) );
         for (i = 0; i < nvars; ++i)
            maxsolu[i] = tmpsolu[i];
         maxsoluobj = tmpsoluobj + constobjective;
      }
   }

   /* Check whether the separation objective is positive, i.e., a violated cover was found. */
   if ( SCIPisEfficacious(scip, maxsoluobj) )
   {
      SCIP_Real rhs = -1.0;
      SCIP_Real lhs = 0.0;

      for (i = 0; i < nvars; ++i)
      {
         if ( i < perm[i] )
         {
            maxsolu[i] = maxsolu[i] - 1;
            lhs += vals[i] * maxsolu[i];
         }
         else
         {
            lhs += vals[i] * maxsolu[i];
            rhs += maxsolu[i];
         }
      }

      assert( SCIPisGT(scip, lhs, rhs) );

      /* add cover inequality */
      SCIP_CALL( addSymresackInequality(scip, cons, consdata, maxsolu, rhs, infeasible) );

      if ( ! *infeasible )
         ++(*ngen);
   }

   SCIPfreeBufferArrayNull(scip, &maxsolu);
   SCIPfreeBufferArrayNull(scip, &tmpsolu);
   SCIPfreeBufferArrayNull(scip, &sepaobjective);

   return SCIP_OKAY;
}


/** Upgrade symresack constraints to orbisacks */
static
SCIP_RETCODE orbisackUpgrade(
   SCIP*                 scip,               /**< SCIP pointer */
   int*                  perm,               /**< permutation */
   int                   nvars,              /**< size of perm array */
   SCIP_VAR**            inputvars,          /**< permuted variables array */
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
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   int maxvars;
   int nrows = 0;
   int i;

   assert( scip != NULL );
   assert( perm != NULL );
   assert( nvars > 0 );
   assert( inputvars != NULL );
   assert( success != NULL );

   *success = TRUE;

   maxvars = nvars / 2;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars1, maxvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars2, maxvars) );

   /* check whether permutation is a composition of 2-cycles */
   for (i = 0; i < nvars; ++i)
   {
      /* ignore non-binary variables */
      if ( ! SCIPvarIsBinary(inputvars[i]) )
         continue;

      if ( perm[perm[i]] != i )
      {
         *success = FALSE;
         break;
      }

      if ( perm[i] > i )
      {
         vars1[nrows] = inputvars[i];
         vars2[nrows++] = inputvars[perm[i]];

         assert( nrows <= maxvars );
      }
   }

   /* if permutation can be upgraded to an orbisack */
   if ( *success )
   {
      SCIP_CALL( SCIPcreateConsOrbisack(scip, cons, "orbisack", vars1, vars2, nrows, FALSE, FALSE,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   SCIPfreeBufferArray(scip, &vars2);
   SCIPfreeBufferArray(scip, &vars1);

   return SCIP_OKAY;
}


/*--------------------------------------------------------------------------------------------
 *--------------------------------- SCIP functions -------------------------------------------
 *--------------------------------------------------------------------------------------------*/

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSymresack)
{   /*lint --e{715}*/
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
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

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
   SCIP_CONSDATA* consdata = NULL;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   SCIPdebugMessage("Transforming constraint.\n");

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);
   assert( sourcedata->nvars != 0 );
   assert( sourcedata->vars != NULL );
   assert( sourcedata->vals != NULL );
   assert( sourcedata->perm != NULL );
   assert( sourcedata->invperm != NULL );
   if ( sourcedata->ppupgrade )
   {
      assert( sourcedata->ncycles != 0 );
      assert( sourcedata->cycledecomposition != NULL );
      for (i = 0; i < sourcedata->ncycles; ++i)
      {
         assert( sourcedata->cycledecomposition[i] != NULL );
         assert( sourcedata->cycledecomposition[i][0] != 0 );
      }
   }

   /* create transformed constraint data (copy data where necessary) */
   nvars = sourcedata->nvars;

   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nvars = nvars;

   if ( nvars > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vals, nvars) );
      SCIP_CALL( SCIPgetTransformedVars(scip, nvars, sourcedata->vars, consdata->vars) );
      for (i = 0; i < nvars; ++i)
      {
         SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
      }

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->perm, sourcedata->perm, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->invperm, sourcedata->invperm, nvars) );

      consdata->ppupgrade = sourcedata->ppupgrade;

      if ( sourcedata->ppupgrade )
      {
         consdata->ncycles = sourcedata->ncycles;
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->cycledecomposition, sourcedata->cycledecomposition, sourcedata->ncycles) );
         for (i = 0; i < sourcedata->ncycles; ++i)
         {
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->cycledecomposition[i], sourcedata->cycledecomposition[i], nvars + 1) ); /*lint !e866*/
         }
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

   assert( infeasible != NULL );
   *infeasible = FALSE;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != 0 );

      SCIPdebugMessage("Generating initial symresack cut for constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( initLP(scip, conss[c], infeasible) );
      if ( *infeasible )
         break;
   }
   SCIPdebugMessage("Generated initial symresack cuts.\n");

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

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
      int ngen = 0;

      /* get data of constraint */
      assert( conss[c] != 0 );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nvars, consdata->vars, consdata->vals) );

      SCIPdebugMessage("Separating symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, &ngen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if ( ngen > 0 )
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
      int ngen = 0;

      /* get data of constraint */
      assert( conss[c] != 0 );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, consdata->vals) );

      SCIPdebugMessage("Separating symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, &ngen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if ( ngen > 0 )
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
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

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
      int ngen = 0;

      /* get data of constraint */
      assert( conss[c] != 0 );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nvars, consdata->vars, consdata->vals) );

      SCIPdebugMessage("Enforcing symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, &ngen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* SCIPdebugMessage("Generated symresack inequalities for <%s>: %d\n", SCIPconsGetName(conss[c]), ngen); */

      if ( ngen > 0 )
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
      SCIP_VAR** vars;
      SCIP_Real val;
      int* invperm;
      int* solu;
      int nvars;
      int i;

      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL);
      assert( consdata->nvars > 0 );
      assert( consdata->vars != NULL );
      assert( consdata->invperm != NULL );

      nvars = consdata->nvars;
      vars = consdata->vars;
      invperm = consdata->invperm;

      /* determine solution */
      SCIP_CALL( SCIPallocBufferArray(scip, &solu, nvars) );

      for (i = 0; i < nvars; ++i)
      {
         /* there are no fixed points */
         assert( invperm[i] != i );

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
      for (i = 0; i < nvars; ++i)
      {
         /* there are no fixed points */
         assert( invperm[i] != i );

         /* if pair (i,invperm[i]) is constant */
         if ( solu[i] == solu[invperm[i]] )
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


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSymresack)
{   /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Enforcing method for symresack constraints (relaxation solutions) ...\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      int ngen = 0;

      /* get data of constraint */
      assert( conss[c] != 0 );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, consdata->vals) );

      SCIPdebugMessage("Enforcing symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, &ngen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if ( ngen > 0 )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckSymresack)
{   /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool terminated = FALSE;
   SCIP_VAR** vars;
   SCIP_Real solVal;
   int* invperm;
   int nvars;
   int val1;
   int val2;
   int c;
   int i;

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
      assert( consdata->nvars > 0 );
      assert( consdata->vars != NULL );
      assert( consdata->invperm != NULL );

      SCIPdebugMessage("Check method for symresack constraint <%s> (%d rows) ...\n", SCIPconsGetName(conss[c]), consdata->nvars);

      nvars = consdata->nvars;
      vars = consdata->vars;
      invperm = consdata->invperm;

      /* detect first non-constant pair of variables */
      for (i = 0; i < nvars; ++i)
      {
         /* there are no fixed points */
         assert( invperm[i] != i );

         /* get value of variable i and its inverse */
         solVal = SCIPgetSolVal(scip, sol, vars[i]);
         assert( SCIPisFeasIntegral(scip, solVal) );
         if ( solVal > 0.5 )
            val1 = 1;
         else
            val1 = 0;

         solVal = SCIPgetSolVal(scip, sol, vars[invperm[i]]);
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
            SCIPinfoMessage(scip, NULL, "First non-constant pair (%d, %d) of variables has pattern (0,1).\n", i, invperm[i]);

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
{  /*lint --e{715}*/
   int c;
   int ngen = 0;

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
      int localngen = 0;

      assert( conss[c] != 0 );

      SCIP_CALL( propVariables(scip, conss[c], &infeasible, &localngen) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      ngen += localngen;

      *result = SCIP_DIDNOTFIND;
   }

   if ( ngen > 0 )
   {
      *result = SCIP_REDUCEDDOM;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSymresack)
{  /*lint --e{715}*/
   int c;
   int ngen = 0;
   int oldndelconss;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   oldndelconss = *ndelconss;

   SCIPdebugMessage("Presolving method of symresack constraint handler. Propagating symresack inequalities.\n");
   *result = SCIP_DIDNOTRUN;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_CONSDATA* consdata;
      int localngen = 0;

      assert( conss[c] != 0 );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* avoid trivial problems */
      if ( consdata->nvars == 0 )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         (*ndelconss)++;
      }
      else
      {
         SCIP_CALL( propVariables(scip, conss[c], &infeasible, &localngen) );
      }

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if ( localngen > 0 )
      {
         *nfixedvars += localngen;
         ngen += localngen;
      }

      *result = SCIP_DIDNOTFIND;
   }

   if ( *ndelconss > oldndelconss ||  ngen > 0 )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** Propagation resolution for conflict analysis */
static
SCIP_DECL_CONSRESPROP(consRespropSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* invperm;
   int nvars;
   int i;

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
   assert( consdata->nvars > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->invperm != NULL );

   vars = consdata->vars;
   nvars = consdata->nvars;
   invperm = consdata->invperm;

   assert( 0 <= inferinfo && inferinfo < (2 * nvars - 1) );

   /* if first part of variable pair was fixed to 0 */
   if ( inferinfo < nvars )
   {
      assert( vars[invperm[inferinfo]] == infervar );
      assert( SCIPvarGetUbAtIndex(vars[invperm[inferinfo]], bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(vars[invperm[inferinfo]], bdchgidx, TRUE) < 0.5 );

      if ( SCIPvarGetUbAtIndex(vars[invperm[inferinfo]], bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(vars[invperm[inferinfo]], bdchgidx, TRUE) < 0.5 )
      {
         SCIPdebugMessage(" -> reason for setting x[%d] = 0 was fixing x[%d] to 0 and each pair of binary variables before (%d,%d) which are not fixed points is constant.\n",
            invperm[inferinfo], inferinfo, inferinfo, invperm[inferinfo]);

         SCIP_CALL( SCIPaddConflictUb(scip, vars[inferinfo], bdchgidx) );

         for (i = 0; i < inferinfo; ++i)
         {
            /* there are no fixed points */
            assert( invperm[i] != i );

            SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[invperm[i]], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[invperm[i]], bdchgidx) );
         }
      }
   }
   /* if second part of variable pair was fixed to 1 */
   else
   {
      int inferinfo2;

      inferinfo2 = inferinfo - nvars;
      assert( vars[inferinfo2] == infervar );
      assert( SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, FALSE) < 0.5 && SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, TRUE) > 0.5 );

      if ( SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, FALSE) < 0.5 && SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, TRUE) > 0.5 )
      {
         SCIPdebugMessage(" -> reason for setting x[%d] = 1 was fixing x[%d] to 1 and each pair of binary variables before (%d,%d) which are not fixed points is constant.\n",
            inferinfo2, invperm[inferinfo2], inferinfo2, invperm[inferinfo2]);

         SCIP_CALL( SCIPaddConflictLb(scip, vars[invperm[inferinfo2]], bdchgidx) );

         for (i = 0; i < inferinfo2; ++i)
         {
            /* there are no fixed points */
            assert( invperm[i] != i );

            SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[invperm[i]], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[invperm[i]], bdchgidx) );
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
 *   in the FD inequality are rounded down, we therefor call
 *   SCIPaddVarLocks(..., nlockspos, nlocksneg).
 * - Symresack constraints may get violated if the variables with a positive coefficient
 *   in the FD inequality are rounded up, we therefor call
 *   SCIPaddVarLocks(..., nlocksneg, nlockspo ).
 */
static
SCIP_DECL_CONSLOCK(consLockSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* perm;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMessage("Locking method for symresack constraint handler.\n");

   /* get data of original constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nvars > 0 );
   assert( consdata->vars != NULL );
   assert( consdata->perm != NULL );

   nvars = consdata->nvars;
   vars = consdata->vars;
   perm = consdata->perm;

   for (i = 0; i < nvars; ++i)
   {
      /* due to clean-up in consdataCreate, there are no fixed points */
      assert( perm[i] != i );

      if ( perm[i] > i )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, vars[i], nlockspos, nlocksneg) );
      }
      else
      {
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
{  /*lint --e{715}*/

   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int* perm;
   SCIP_Bool* covered;
   int i;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   assert( consdata->nvars > 0 );
   assert( consdata->perm != NULL );

   vars = consdata->vars;
   nvars = consdata->nvars;
   perm = consdata->perm;

   SCIPdebugMsg(scip, "Printing method for symresack constraint handler\n");

   SCIP_CALL( SCIPallocBufferArray(scip, &covered, nvars) );
   for (i = 0; i < nvars; ++i)
      covered[i] = FALSE;

   if ( consdata->ppupgrade )
      SCIPinfoMessage(scip, file, "ppSymresack(");
   else
      SCIPinfoMessage(scip, file, "symresack(");

   for (i = 0; i < nvars; ++i)
   {
      if ( covered[i] )
         continue;

      /* print cycle of perm containing i */
      SCIPinfoMessage(scip, file, "[%s", SCIPvarGetName(vars[i]));
      covered[i] = TRUE;
      j = perm[i];
      while ( j != i )
      {
         SCIPinfoMessage(scip, file, ",%s", SCIPvarGetName(vars[j]));
         covered[j] = TRUE;
         j = perm[j];
      }
      SCIPinfoMessage(scip, file, "]");
   }
   SCIPinfoMessage(scip, file, ")\n");

   SCIPfreeBufferArray(scip, &covered);

   return SCIP_OKAY;
}


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
   assert( conshdlr != NULL );

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSymresack) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSymresack) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSymresack) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSymresack, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSymresack) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSymresack, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSymresack) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSymresack, consSepasolSymresack, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSymresack) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSymresack) );

   /* whether we allow upgrading to orbisack constraints*/
   SCIP_CALL( SCIPaddBoolParam(scip, "cons/" CONSHDLR_NAME "/upgrade",
         "Upgrade symresack constraints to orbisack constraints?",
         &conshdlrdata->symresackupgrade, TRUE, DEFAULT_UPGRADE, NULL, NULL) );

   /* whether we allow upgrading to packing/partioning symresack constraints*/
   SCIP_CALL( SCIPaddBoolParam(scip, "cons/" CONSHDLR_NAME "/ppsymresack",
         "Upgrade symresack constraints to packing/partioning symresacks?",
         &conshdlrdata->checkppsymresack, TRUE, DEFAULT_PPSYMRESACK, NULL, NULL) );

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates and captures a symresack constraint
 *
 *  In a presolving step, we check whether the permutation acts only on binary points. Otherwise, we eliminate
 *  the non-binary variables from the permutation.
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in problem */
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
   SCIP_Bool success;

   assert( cons != NULL );
   assert( nvars > 0 );

   success = FALSE;

   /* find the symresack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("Symresack constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* check whether constraint can be upgraded to an orbisack constraint */
   SCIP_CALL( SCIPgetBoolParam(scip, "cons/symresack/upgrade", &upgrade) );

   if ( upgrade )
   {
      SCIP_CALL( orbisackUpgrade(scip, perm, nvars, vars, &success, cons, initial, separate, enforce, check, propagate,
            local, modifiable, dynamic, removable, stickingatnode) );

      if ( success )
      {
         SCIPdebugMessage("Upgraded symresack constraint to orbisack constraint.\n");

         return SCIP_OKAY;
      }
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, vars, nvars, perm) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate && (! consdata->ppupgrade), enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}


/** creates and captures a symresack constraint
 *  in its most basic variant, i.e., with all constraint flags set to their default values
 *
 *  In a presolving step, we remove all fixed points and cycles that act on non-binary variables of the permutation
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsBasicSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars               /**< number of variables in problem */
   )
{
   SCIP_CALL( SCIPcreateConsSymresack(scip, cons, name, perm, vars, nvars,
         TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
