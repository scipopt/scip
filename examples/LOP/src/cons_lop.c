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

/* uncomment for debug output: */
/* #define SCIP_DEBUG */

/**@file   cons_lop.c
 * @brief  constraint handler for linear ordering constraints
 * @author Marc Pfetsch
 *
 * We handle the following system of linear constraints:
 * - \f$ x_{ij} + x_{ji} = 1 \f$ for \f$i < j\f$                               (symmetry equations - added initially)
 * - \f$ x_{ij} + x_{jk} + x_{ki} \leq 2 \f$ for \f$i < j, i < k, j \neq k\f$  (triangle inequalities - separated)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "cons_lop.h"

#include <assert.h>
#include <string.h>


/* constraint handler properties */
#define CONSHDLR_NAME          "lop"
#define CONSHDLR_DESC          "linear ordering constraint handler"
#define CONSHDLR_SEPAPRIORITY       100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY      -100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY     -100 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP


/** constraint data for linear ordering constraints */
struct SCIP_ConsData
{
   int                   n;                  /**< number of elements */
   SCIP_VAR***           vars;               /**< variables */
};


/** separate symmetry equations and triangle inequalities */
static
SCIP_RETCODE LOPseparate(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   n,                  /**< number of elements */
   SCIP_VAR***           vars,               /**< n x n matrix of variables */
   SCIP_SOL*             sol,                /**< solution to be separated */
   int*                  ngen,               /**< output: pointer to store number of added rows */
   SCIP_Bool*            cutoff              /**< output: pointer to store whether we detected a cutoff */
   )
{
   char s[SCIP_MAXSTRLEN];
   int i;
   int j;
   int k;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( ngen != NULL );
   assert( cutoff != NULL );

   /* Consider all (i,j,k) with i < j, i < k, j != k; since the inequalities are symmetric under cyclic shifts, we can
    * assume i to be the smallest index. */
   *cutoff = FALSE;
   for (i = 0; i < n && ! (*cutoff); ++i)
   {
      for (j = i+1; j < n && ! (*cutoff); ++j)
      {
	 SCIP_Real valIJ;

	 valIJ = SCIPgetSolVal(scip, sol, vars[i][j]);

	 /* if symmetry equations are violated - should not be the case, if they are added in the beginning */
	 if ( ! SCIPisFeasEQ(scip, valIJ + SCIPgetSolVal(scip, sol, vars[j][i]), 1.0) )
	 {
	    SCIP_ROW *row;

	    (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sym#%d#%d", i, j);

	    SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, s, 1.0, 1.0, FALSE, FALSE, TRUE) );
	    SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
	    SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], 1.0) );
	    SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j][i], 1.0) );
	    SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
	    SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
	    SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
	    SCIP_CALL( SCIPreleaseRow(scip, &row));
	    ++(*ngen);

            if ( *cutoff )
               break;
	 }

	 /* check triangle inequalities */
	 for (k = i+1; k < n; ++k)
	 {
	    SCIP_Real sum;

	    if ( k == j )
	       continue;

	    sum = valIJ + SCIPgetSolVal(scip, sol, vars[j][k]) + SCIPgetSolVal(scip, sol, vars[k][i]);

	    /* if sum - 2.0 > 0, i.e., the cut is violated */
	    if ( SCIPisEfficacious(scip, sum - 2.0) )
	    {
	       SCIP_ROW *row;

	       (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "triangle#%d#%d#%d", i, j, k);

	       SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, s, -SCIPinfinity(scip), 2.0, FALSE, FALSE, TRUE) );
	       SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], 1.0) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j][k], 1.0) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[k][i], 1.0) );
	       SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
	       SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
	       SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
	       SCIP_CALL( SCIPreleaseRow(scip, &row));
	       ++(*ngen);

               if ( *cutoff )
                  break;
	    }
	 }
      }
   }

   return SCIP_OKAY;
}


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyLOP)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( valid != NULL );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrLOP(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteLOP)
{  /*lint --e{715}*/
   int i;
   int n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( consdata != NULL);
   assert( *consdata != NULL);
   assert( (*consdata)->vars != NULL );

   SCIPdebugMsg(scip, "deleting linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

   n = (*consdata)->n;
   for (i = 0; i < n; ++i)
      SCIPfreeBlockMemoryArray(scip, &((*consdata)->vars[i]), n); /*lint !e866*/
   SCIPfreeBlockMemoryArray(scip, &((*consdata)->vars), n);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed)
 *
 *  We output the final linear ordering.
 */
static
SCIP_DECL_CONSEXIT(consExitLOP)
{  /*lint --e{715}*/
   SCIP_SOL* sol;
   int c;
   int i;
   int j;
   int n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMsg(scip, "exiting linear ordering constraint handler <%s>.\n", SCIPconshdlrGetName(conshdlr));

   /* avoid output for subscips */
   if ( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   /* get best solution */
   sol = SCIPgetBestSol(scip);
   if ( sol == NULL )
      return SCIP_OKAY;

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR*** vars;
      int* outdeg;
      int* indices;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      SCIPdebugMsg(scip, "solution for for linear ordering constraint <%s>.\n", SCIPconsGetName(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );
      n = consdata->n;
      vars = consdata->vars;

      SCIP_CALL( SCIPallocBufferArray(scip, &outdeg, n) );
      SCIP_CALL( SCIPallocBufferArray(scip, &indices, n) );

      /* compute out-degree */
      for (i = 0; i < n; ++i)
      {
	 int deg = 0;
	 for (j = 0; j < n; ++j)
	 {
	    SCIP_Real val;

	    if (j == i)
	       continue;

	    val = SCIPgetSolVal(scip, sol, vars[i][j]);
	    assert( SCIPisFeasIntegral(scip, val) );
	    if ( val < 0.5 )
	       ++deg;
	 }
	 outdeg[i] = deg;
	 indices[i] = i;
      }

      /* sort such that degrees are non-decreasing */
      SCIPsortIntInt(outdeg, indices, n);

      /* output */
      SCIPinfoMessage(scip, NULL, "\nFinal order of linear ordering constraint <%s>:\n", SCIPconsGetName(conss[c]));
      for (i = 0; i < n; ++i)
	 SCIPinfoMessage(scip, NULL, "%d ", indices[i]);
      SCIPinfoMessage(scip, NULL, "\n");

      SCIPfreeBufferArray(scip, &indices);
      SCIPfreeBufferArray(scip, &outdeg);
   }

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransLOP)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_CONSDATA* sourcedata;
   int i;
   int j;
   int n;
   char s[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   SCIPdebugMsg(scip, "transforming linear ordering constraint <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   n = sourcedata->n;
   consdata->n = n;

   /* transform variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, n) );
   for (i = 0; i < n; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->vars[i]), n) ); /*lint !e866*/
      for (j = 0; j < n; ++j)
      {
	 if (j != i)
	 {
	    assert( sourcedata->vars[i][j] != NULL );
	    SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->vars[i][j], &(consdata->vars[i][j])) );
	 }
      }
   }

   /* create constraint */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));

   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpLOP)
{  /*lint --e{715}*/
   char s[SCIP_MAXSTRLEN];
   int c;
   int ngen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR*** vars;
      int i;
      int j;
      int n;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      SCIPdebugMsg(scip, "adding initial rows for linear ordering constraint <%s>.\n", SCIPconsGetName(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );
      n = consdata->n;
      vars = consdata->vars;

      /* add symmetry equation */
      for (i = 0; i < n; ++i)
      {
	 for (j = i+1; j < n; ++j)
	 {
	    SCIP_ROW* row;

	    (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sym#%d#%d", i, j);
	    SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, s, 1.0, 1.0, FALSE, FALSE, FALSE) );
	    SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
	    SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], 1.0) );
	    SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j][i], 1.0) );
	    SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
	    SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
	    SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
	    SCIP_CALL( SCIPreleaseRow(scip, &row));
	    ++ngen;

            /* cannot handle infeasible case here - just exit */
            if ( *infeasible )
               return SCIP_OKAY;
	 }
      }
   }
   SCIPdebugMsg(scip, "added %d equations.\n", ngen);

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpLOP)
{  /*lint --e{715}*/
   int ngen = 0;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( conss != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      SCIP_Bool cutoff;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "separating LP solution for linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( LOPseparate(scip, conshdlr, consdata->n, consdata->vars, NULL, &ngen, &cutoff) );
      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }
   if ( ngen > 0 )
      *result = SCIP_SEPARATED;
   SCIPdebugMsg(scip, "separated %d cuts.\n", ngen);

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolLOP)
{  /*lint --e{715}*/
   int ngen = 0;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( conss != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      SCIP_Bool cutoff;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "separating solution for linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( LOPseparate(scip, conshdlr, consdata->n, consdata->vars, sol, &ngen, &cutoff) );
      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }
   if ( ngen > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpLOP)
{  /*lint --e{715}*/
   char s[SCIP_MAXSTRLEN];
   int ngen = 0;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( conss != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      SCIP_VAR*** vars;
      int i;
      int j;
      int k;
      int n;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "enforcing lp solution for linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      n = consdata->n;
      vars = consdata->vars;
      assert( vars != NULL );

      for (i = 0; i < n; ++i)
      {
	 for (j = i + 1; j < n; ++j)
	 {
	    SCIP_Real valIJ;

	    valIJ = SCIPgetSolVal(scip, NULL, vars[i][j]);

	    /* if symmetry equations are violated - should not be the case, if they are added in the beginning */
	    if ( ! SCIPisFeasEQ(scip, 1.0 - valIJ, SCIPgetSolVal(scip, NULL, vars[j][i])) )
	    {
	       SCIP_ROW *row;
               SCIP_Bool infeasible;

	       (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sym#%d#%d", i, j);

	       SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, s, 1.0, 1.0, FALSE, FALSE, TRUE) );
	       SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], 1.0) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j][i], 1.0) );
	       SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
	       SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
	       SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
	       SCIP_CALL( SCIPreleaseRow(scip, &row));
	       ++ngen;

               if ( infeasible )
               {
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
	    }

	    /* enforce triangle inequalities */
	    for (k = i + 1; k < n; ++k)
	    {
	       SCIP_Real sum;

	       if ( k == j )
		  continue;

	       sum = valIJ + SCIPgetSolVal(scip, NULL, vars[j][k]) + SCIPgetSolVal(scip, NULL, vars[k][i]);

	       /* if sum > 2.0, i.e., the cut is violated */
	       if ( SCIPisFeasGT(scip, sum, 2.0) ) /* this is the only difference to the separation call */
	       {
		  SCIP_ROW *row;
                  SCIP_Bool infeasible;

		  (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "triangle#%d#%d#%d", i, j, k);

		  SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, s, -SCIPinfinity(scip), 2.0, FALSE, FALSE, TRUE) );
		  SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
		  SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], 1.0) );
		  SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j][k], 1.0) );
		  SCIP_CALL( SCIPaddVarToRow(scip, row, vars[k][i], 1.0) );
		  SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
		  SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
		  SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
		  SCIP_CALL( SCIPreleaseRow(scip, &row));
		  ++ngen;

                  if ( infeasible )
                  {
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
	       }
	    }
	 }
      }

      if ( ngen > 0 )
      {
	 *result = SCIP_SEPARATED;
	 return SCIP_OKAY;
      }

   }
   SCIPdebugMsg(scip, "all linear ordering constraints are feasible.\n");
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsLOP)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( conss != NULL );
   assert( result != NULL );

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      SCIP_VAR*** vars;
      int i;
      int j;
      int k;
      int n;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "enforcing pseudo solution for linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );
      vars = consdata->vars;
      n = consdata->n;

      /* check triangle inequalities */
      for (i = 0; i < n; ++i)
      {
	 for (j = i + 1; j < n; ++j)
	 {
	    SCIP_Bool oneIJ;
            SCIP_Bool oneJI;

	    /* the priorities should ensure that the solution is integral */
	    assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, vars[i][j])) );
	    assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, vars[j][i])) );

            oneIJ = SCIPgetSolVal(scip, NULL, vars[i][j]) > 0.5 ? TRUE : FALSE;
            oneJI = SCIPgetSolVal(scip, NULL, vars[j][i]) > 0.5 ? TRUE : FALSE;

	    if ( oneIJ == oneJI )
	    {
	       SCIPdebugMsg(scip, "constraint <%s> infeasible (violated equation).\n", SCIPconsGetName(cons));
	       *result = SCIP_INFEASIBLE;
	       return SCIP_OKAY;
	    }

	    for (k = i + 1; k < n; ++k)
	    {
	       SCIP_Bool oneJK;
               SCIP_Bool oneKI;

	       if ( k == j )
		  continue;

	       assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, vars[j][k])) );
	       assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, vars[k][i])) );

	       oneJK = SCIPgetSolVal(scip, NULL, vars[j][k]) > 0.5 ? TRUE : FALSE;
	       oneKI = SCIPgetSolVal(scip, NULL, vars[k][i]) > 0.5 ? TRUE : FALSE;

	       /* if triangle inequality is violated */
	       if ( oneIJ && oneJK && oneKI )
	       {
		  SCIPdebugMsg(scip, "constraint <%s> infeasible (violated triangle ineq.).\n", SCIPconsGetName(cons));
		  *result = SCIP_INFEASIBLE;
		  return SCIP_OKAY;
	       }
	    }
	 }
      }
   }
   SCIPdebugMsg(scip, "all linear ordering constraints are feasible.\n");
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckLOP)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( conss != NULL );
   assert( result != NULL );

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      SCIP_VAR*** vars;
      int i;
      int j;
      int k;
      int n;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "checking linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );
      vars = consdata->vars;
      n = consdata->n;

      /* check triangle inequalities and symmetry equations */
      for (i = 0; i < n; ++i)
      {
	 for (j = i + 1; j < n; ++j)
	 {
	    SCIP_Bool oneIJ;
            SCIP_Bool oneJI;

	    /* the priorities should ensure that the solution is integral */
	    assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[i][j])) );
	    assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[j][i])) );

	    oneIJ = SCIPgetSolVal(scip, sol, vars[i][j]) > 0.5 ? TRUE : FALSE;
            oneJI = SCIPgetSolVal(scip, sol, vars[j][i]) > 0.5 ? TRUE : FALSE;

	    /* check symmetry equations */
	    if ( oneIJ == oneJI )
	    {
	       SCIPdebugMsg(scip, "constraint <%s> infeasible (violated equation).\n", SCIPconsGetName(cons));
	       *result = SCIP_INFEASIBLE;
               if( printreason )
               {
                  SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
                  SCIPinfoMessage(scip, NULL, "violation: symmetry equation violated <%s> = %.15g and <%s> = %.15g\n",
                     SCIPvarGetName(vars[i][j]), SCIPgetSolVal(scip, sol, vars[i][j]),
                     SCIPvarGetName(vars[j][i]), SCIPgetSolVal(scip, sol, vars[j][i]));
               }
	       return SCIP_OKAY;
	    }

	    for (k = i + 1; k < n; ++k)
	    {
	       SCIP_Bool oneJK;
               SCIP_Bool oneKI;

	       if ( k == j )
		  continue;

	       assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[j][k])) );
	       assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[k][i])) );

	       oneJK = SCIPisGT(scip, SCIPgetSolVal(scip, sol, vars[j][k]), 0.5);
	       oneKI = SCIPisGT(scip, SCIPgetSolVal(scip, sol, vars[k][i]), 0.5);

	       /* if triangle inequality is violated */
	       if ( oneIJ && oneJK && oneKI )
	       {
		  SCIPdebugMsg(scip, "constraint <%s> infeasible (violated triangle ineq.).\n", SCIPconsGetName(cons));
		  *result = SCIP_INFEASIBLE;
                  if( printreason )
                  {
                     SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
                     SCIPinfoMessage(scip, NULL,
                        "violation: triangle inequality violated <%s> = %.15g, <%s> = %.15g, <%s> = %.15g\n",
                        SCIPvarGetName(vars[i][j]), SCIPgetSolVal(scip, sol, vars[i][j]),
                        SCIPvarGetName(vars[j][k]), SCIPgetSolVal(scip, sol, vars[j][k]),
                        SCIPvarGetName(vars[k][i]), SCIPgetSolVal(scip, sol, vars[k][i]));
                  }
		  return SCIP_OKAY;
	       }
	    }
	 }
      }
   }
   SCIPdebugMsg(scip, "all linear ordering constraints are feasible.\n");
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropLOP)
{  /*lint --e{715}*/
   int c;
   int ngen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( conss != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      SCIP_VAR*** vars;
      int i;
      int j;
      int k;
      int n;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "propagating linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

      *result = SCIP_DIDNOTFIND;

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );

      vars = consdata->vars;
      n = consdata->n;

      /* check triangle inequalities */
      for (i = 0; i < n; ++i)
      {
	 for (j = i + 1; j < n; ++j)
	 {
            SCIP_Bool infeasible;
            SCIP_Bool tightened;

            /* for consistency make sure that the complementarity constraints are satisfied */

            /* if x[i][j] == 1 then x[j][i] = 0 */
	    if ( SCIPvarGetLbLocal(vars[i][j]) > 0.5 )
	    {
	       SCIP_CALL( SCIPinferBinvarCons(scip, vars[j][i], FALSE, cons, i*n + j, &infeasible, &tightened) );
	       if ( infeasible )
	       {
		  SCIPdebugMsg(scip, " -> node infeasible.\n");
                  SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[j][i]) );
                  SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
		  *result = SCIP_CUTOFF;
		  return SCIP_OKAY;
	       }
	       if ( tightened )
		  ++ngen;
	    }

	    /* if x[i][j] == 0 then x[j][i] = 1 */
	    if ( SCIPvarGetUbLocal(vars[i][j]) < 0.5 )
	    {
	       SCIP_CALL( SCIPinferBinvarCons(scip, vars[j][i], TRUE, cons, i*n + j, &infeasible, &tightened) );
	       if ( infeasible )
	       {
		  SCIPdebugMsg(scip, " -> node infeasible.\n");
                  SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars[j][i]) );
                  SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
		  *result = SCIP_CUTOFF;
		  return SCIP_OKAY;
	       }
	       if ( tightened )
		  ++ngen;
	    }

            /* check whether triangle inequality allows to fix variables */
	    for (k = i + 1; k < n; ++k)
	    {
	       if ( k == j )
		  continue;

	       if ( SCIPvarGetLbLocal(vars[i][j]) > 0.5 )
	       {
                  if ( SCIPvarGetLbLocal(vars[j][k]) > 0.5 )
                  {
                     /* if x[i][j] == 1 and x[j][k] == 1 then x[k][i] = 0 */
                     SCIP_CALL( SCIPinferBinvarCons(scip, vars[k][i], FALSE, cons, n*n + i*n*n + j*n + k, &infeasible, &tightened) );
                     if ( infeasible )
                     {
                        SCIPdebugMsg(scip, " -> node infeasible.\n");
                        SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
                        SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) );
                        SCIP_CALL( SCIPaddConflictBinvar(scip, vars[j][k]) );
                        SCIP_CALL( SCIPaddConflictBinvar(scip, vars[k][i]) );
                        SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
                        *result = SCIP_CUTOFF;
                        return SCIP_OKAY;
                     }
                     if ( tightened )
                        ++ngen;
                  }

                  if ( SCIPvarGetLbLocal(vars[k][i]) > 0.5 )
                  {
                     /* if x[k][i] == 1 and x[i][j] = 1 then x[j][k] = 0 */
                     SCIP_CALL( SCIPinferBinvarCons(scip, vars[j][k], FALSE, cons, n*n + i*n*n + j*n + k, &infeasible, &tightened) );
                     if ( infeasible )
                     {
                        SCIPdebugMsg(scip, " -> node infeasible.\n");
                        SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
                        SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) );
                        SCIP_CALL( SCIPaddConflictBinvar(scip, vars[j][k]) );
                        SCIP_CALL( SCIPaddConflictBinvar(scip, vars[k][i]) );
                        SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
                        *result = SCIP_CUTOFF;
                        return SCIP_OKAY;
                     }
                     if ( tightened )
                        ++ngen;
                  }
               }

               /* if x[j][k] == 1 and x[k][i] == 1 then x[i][j] = 0 */
	       if ( SCIPvarGetLbLocal(vars[j][k]) > 0.5 && SCIPvarGetLbLocal(vars[k][i]) > 0.5 )
	       {
		  SCIP_CALL( SCIPinferBinvarCons(scip, vars[i][j], FALSE, cons, n*n + i*n*n + j*n + k, &infeasible, &tightened) );
		  if ( infeasible )
		  {
		     SCIPdebugMsg(scip, " -> node infeasible.\n");
                     SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );
                     SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) );
                     SCIP_CALL( SCIPaddConflictBinvar(scip, vars[j][k]) );
                     SCIP_CALL( SCIPaddConflictBinvar(scip, vars[k][i]) );
                     SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
		     *result = SCIP_CUTOFF;
		     return SCIP_OKAY;
		  }
		  if ( tightened )
		     ++ngen;
	       }

	       /* all other implications occur with other indices i, j, k */
	    }
	 }
      }
   }
   SCIPdebugMsg(scip, "propagated %d domains.\n", ngen);
   if ( ngen > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}

/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropLOP)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int n;
   int nsqrd;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Propagation resolution of constraint <%s>.\n", SCIPconsGetName(cons));
   *result = SCIP_DIDNOTFIND;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->vars != NULL );

   n = consdata->n;
   nsqrd = n * n;
   vars = consdata->vars;

   assert( 0 <= inferinfo && inferinfo < n*n + n*n*n );

   /* if the conflict came from an equation */
   if ( inferinfo < nsqrd )
   {
      int index1;
      int index2;

      index1 = inferinfo/n;
      index2 = inferinfo % n;
      assert( 0 <= index1 && index1 < n );
      assert( 0 <= index2 && index2 < n );
      assert( vars[index2][index1] == infervar );

      /* if the variable was fixed to 0 */
      if ( SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, FALSE) > 0.5 && SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5 )
      {
	 SCIPdebugMsg(scip, " -> reason for x[%d][%d] == 0 was x[%d][%d] = 1.\n", index2, index1, index1, index2);
	 /* the reason was that x[i][j] was fixed to 1 */
	 SCIP_CALL( SCIPaddConflictLb(scip, vars[index1][index2], bdchgidx) );
	 *result = SCIP_SUCCESS;
	 return SCIP_OKAY;
      }

      /* if the variable was fixed to 1 */
      if ( SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, FALSE) < 0.5 && SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE) > 0.5 )
      {
	 SCIPdebugMsg(scip, " -> reason for x[%d][%d] == 1 was x[%d][%d] = 0.\n", index2, index1, index1, index2);
	 /* the reason was that x[i][j] was fixed to 0 */
	 SCIP_CALL( SCIPaddConflictUb(scip, vars[index1][index2], bdchgidx) );
	 *result = SCIP_SUCCESS;
	 return SCIP_OKAY;
      }
   }
   else
   {
      /* otherwise the conflict came from a triangle inequality */
      int index1;
      int index2;
      int index3;

      index1 = (inferinfo - nsqrd) / nsqrd;
      index2 = (inferinfo - nsqrd - index1 * nsqrd) / n;
      index3 = (inferinfo - nsqrd) % n;

      assert( 0 <= index1 && index1 < n );
      assert( 0 <= index2 && index2 < n );
      assert( 0 <= index3 && index3 < n );
      assert( index1 < index2 );
      assert( index1 < index3 );
      assert( index2 != index3 );

      if ( vars[index3][index1] == infervar )
      {
         /* the variable should have been fixed to 0 */
         assert( SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, FALSE) > 0.5 && SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5 );

         /* the reason was that x[index1][index2] and x[index2][index3] were fixed to 1 */
         SCIPdebugMsg(scip, " -> reason for x[%d][%d] == 0 was x[%d][%d] = x[%d][%d] = 1.\n", index3, index1, index1, index2, index2, index3);
         SCIP_CALL( SCIPaddConflictLb(scip, vars[index1][index2], bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, vars[index2][index3], bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else if ( vars[index2][index3] == infervar )
      {
         /* the variable should have been fixed to 0 */
         assert( SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, FALSE) > 0.5 && SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5 );

         /* the reason was that x[index1][index2] and x[index3][index1] were fixed to 1 */
         SCIPdebugMsg(scip, " -> reason for x[%d][%d] == 0 was x[%d][%d] = x[%d][%d] = 1.\n", index2, index3, index1, index2, index3, index1);
         SCIP_CALL( SCIPaddConflictLb(scip, vars[index1][index2], bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, vars[index3][index1], bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else if ( vars[index1][index2] == infervar )
      {
         /* the variable should have been fixed to 0 */
         assert( SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, FALSE) > 0.5 && SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) < 0.5 );

         /* the reason was that x[index2][index3] and x[index3][index1] were fixed to 1 */
         SCIPdebugMsg(scip, " -> reason for x[%d][%d] == 0 was x[%d][%d] = x[%d][%d] = 1.\n", index1, index2, index2, index3, index3, index1);
         SCIP_CALL( SCIPaddConflictLb(scip, vars[index2][index3], bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, vars[index3][index1], bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         /* should not happen */
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockLOP)
{  /*lint --e{715}*/
   int i;
   int j;
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMsg(scip, "Locking linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->vars != NULL );
   n = consdata->n;
   vars = consdata->vars;

   for (i = 0; i < n; ++i)
   {
      for (j = 0; j < n; ++j)
      {
	 if ( i != j )
	 {
	    /* the constraint may be violated in any way */
	    SCIP_CALL( SCIPaddVarLocksType(scip, vars[i][j], SCIP_LOCKTYPE_MODEL, nlockspos + nlocksneg, nlockspos + nlocksneg) );
	 }
      }
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintLOP)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int i;
   int j;
   int n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   n = consdata->n;
   vars = consdata->vars;

   SCIPinfoMessage(scip, file, "LOP[");
   for (i = 0; i < n; ++i)
   {
      if ( i > 0 )
	 SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "(");
      for (j = 0; j < n; ++j)
      {
	 if ( j != i )
	 {
	    if ( j > 0 && (i > 0 || j > 1) )
	       SCIPinfoMessage(scip, file, ",");
	    SCIPinfoMessage(scip, file, "%s", SCIPvarGetName(vars[i][j]));
	 }
      }
      SCIPinfoMessage(scip, file, ")");
   }
   SCIPinfoMessage(scip, file, "]\n");

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyLOP)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_VAR*** sourcevars;
   SCIP_VAR*** vars;
   int i;
   int j;
   int n;

   assert( scip != NULL );
   assert( sourceconshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(sourceconshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( sourcescip != NULL );
   assert( sourcecons != NULL );
   assert( varmap != NULL );
   assert( valid != NULL );

   *valid = TRUE;

   SCIPdebugMsg(scip, "Copying method for linear ordering constraint handler.\n");

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   n = sourcedata->n;
   sourcevars = sourcedata->vars;
   assert( sourcevars != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, n) );
   BMSclearMemoryArray(vars, n);

   for (i = 0; i < n; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(vars[i]), n) ); /*lint !e866*/

      for (j = 0; j < n && *valid; ++j)
      {
         if ( i != j )
         {
            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[i][j], &vars[i][j], varmap, consmap, global, valid) );
            assert( !(*valid) || vars[i][j] != NULL );
         }
      }
   }

   if ( *valid )
   {
      /* create copied constraint */
      if ( name == 0 )
         name = SCIPconsGetName(sourcecons);

      SCIP_CALL( SCIPcreateConsLOP(scip, cons, name, n, vars,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   /* free memory in reverse order */
   for (i = n-1; i >= 0; --i)
      SCIPfreeBufferArrayNull(scip, &vars[i]);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** creates the handler for linear ordering constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrLOP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpLOP, consEnfopsLOP, consCheckLOP, consLockLOP, NULL) );
   assert( conshdlr != NULL );

   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteLOP) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitLOP) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyLOP, consCopyLOP) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransLOP) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpLOP) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpLOP, consSepasolLOP,
         CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropLOP, CONSHDLR_PROPFREQ,
         CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropLOP) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintLOP) );

   return SCIP_OKAY;
}

/** creates and captures a linear ordering constraint */
SCIP_RETCODE SCIPcreateConsLOP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   n,                  /**< number of elements */
   SCIP_VAR***           vars,               /**< n x n matrix of binary variables */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int i;
   int j;

   /* find the linear ordering constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if (conshdlr == NULL)
   {
      SCIPerrorMessage("linear ordering constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->n = n;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, n) );
   for (i = 0; i < n; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->vars[i]), n) ); /*lint !e866*/
      for (j = 0; j < n; ++j)
      {
	 if ( j != i )
	 {
	    assert( vars[i][j] != NULL );
	    consdata->vars[i][j] = vars[i][j];
	 }
      }
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}
