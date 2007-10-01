/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_LO.c,v 1.1 2007/10/01 13:41:53 bzfpfets Exp $"
//#define SCIP_DEBUG

/**@file   cons_LO.c
 * @brief  example constraint handler for linear ordering constraints
 * @author Marc Pfetsch
 *
 * We handle the following inequality system:
 * x[i][j] + x[j][i] == 1  (added initially)
 * x[i][j] + x[j][k] + x[k][i] <= 2 (triangle inequalities)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "cons_LO.h"

#include <assert.h>
#include <string.h>
#include <scip/cons_setppc.h>


/* constraint handler properties */
#define CONSHDLR_NAME          "LO"
#define CONSHDLR_DESC          "linear ordering constraint handler"
#define CONSHDLR_SEPAPRIORITY       100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY      -100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY     -100 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            10 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */


/** constraint data for LO constraints */
struct SCIP_ConsData
{
   int n;             /**< number of elements */
   SCIP_VAR*** Vars;  /**< variables */
};



/** separate triangle inequalities */
static
SCIP_RETCODE LOseparateTriangle(
		SCIP* scip,         /**< SCIP pointer */
		int n,              /**< number of elements */
		SCIP_VAR*** Vars,   /**< n x n matrix of variables */
		SCIP_SOL* sol,      /**< solution to be separated */
		int* nGen           /**< output: number of added rows */
		)
{
   int i, j, k;

   assert( scip != NULL );
   assert( Vars != NULL );
   assert( nGen != NULL );

   for (i = 0; i < n; ++i)
   {
      for (j = 0; j < n; ++j)
      {
	 SCIP_Real valIJ = 0.0;
	 if (j == i)
	    continue;

	 valIJ = SCIPgetSolVal(scip, sol, Vars[i][j]);
	 assert( SCIPisFeasEQ(scip, valIJ + SCIPgetSolVal(scip, sol, Vars[j][i]), 1.0) );
	 for (k = 0; k < n; ++k)
	 {
	    SCIP_Real sum = 0.0;
	    if (k == i || k == j)
	       continue;

	    sum = valIJ + SCIPgetSolVal(scip, sol, Vars[j][k]) + SCIPgetSolVal(scip, sol, Vars[k][i]);

	    /* if sum - 2.0 > 0, i.e., the cut is violated */
	    if ( SCIPisEfficacious(scip, sum - 2.0) )
	    {
	       SCIP_ROW *row;
	       char s[SCIP_MAXSTRLEN];

	       sprintf(s, "triangle#%d#%d#%d", i, j, k);

	       SCIP_CALL( SCIPcreateEmptyRow(scip, &row, s, -SCIPinfinity(scip), 2.0, FALSE, FALSE, TRUE) );
	       SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, Vars[i][j], 1.0) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, Vars[j][k], 1.0) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, Vars[k][i], 1.0) );
	       SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
	       SCIPdebug( SCIProwPrint(row, NULL) );
#endif
	       SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );
	       SCIP_CALL( SCIPreleaseRow(scip, &row));
	       ++(*nGen);
	    }
	 }
      }
   }

   return SCIP_OKAY;
}


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#define consFreeLO NULL

/** initialization method of constraint handler (called after problem was transformed) */
#define consInitLO NULL

/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitLO NULL

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreLO NULL

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreLO NULL

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolLO NULL

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#define consExitsolLO NULL

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteLO)
{
   int i, n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( consdata != NULL);
   assert( *consdata != NULL);
   assert( (*consdata)->Vars != NULL );

   SCIPdebugMessage("deleting LO constraint <%s>.\n", SCIPconsGetName(cons));

   n = (*consdata)->n;
   for (i = 0; i < n; ++i)
      SCIPfreeBlockMemoryArray(scip, &((*consdata)->Vars[i]), n);
   SCIPfreeBlockMemoryArray(scip, &((*consdata)->Vars), n);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransLO)
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSDATA* sourcedata;
   int i, j, n;
   char s[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   SCIPdebugMessage("transforming LO constraint <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   n = sourcedata->n;
   consdata->n = n;

   /* transfrom variables */
   SCIPallocBlockMemoryArray(scip, &consdata->Vars, n);
   for (i = 0; i < n; ++i)
   {
      SCIPallocBlockMemoryArray(scip, &(consdata->Vars[i]), n);
      for (j = 0; j < n; ++j)
      {
	 if (j == i)
	    continue;
	 assert( sourcedata->Vars[i][j] != NULL );
	 SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->Vars[i][j], &(consdata->Vars[i][j])) );
      }
   }

   /* create constraint */
   sprintf(s, "t_%s", SCIPconsGetName(sourcecons));

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
SCIP_DECL_CONSINITLP(consInitlpLO)
{
   int c;
   int nGen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      int i, j, n;
      SCIP_VAR*** Vars;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      SCIPdebugMessage("adding initial rows for LO constraint <%s>.\n", SCIPconsGetName(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      assert( consdata->Vars != NULL );
      n = consdata->n;
      Vars = consdata->Vars;

      /* add symmetry equation */
      for (i = 0; i < n; ++i)
      {
	 for (j = i+1; j < n; ++j)
	 {
	    char s[SCIP_MAXSTRLEN];
	    SCIP_ROW* row;

	    sprintf(s, "sym#%d#%d", i, j);
	    SCIP_CALL( SCIPcreateEmptyRow(scip, &row, s, 1.0, 1.0, FALSE, FALSE, FALSE) );
	    SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
	    SCIP_CALL( SCIPaddVarToRow(scip, row, Vars[i][j], 1.0) );
	    SCIP_CALL( SCIPaddVarToRow(scip, row, Vars[j][i], 1.0) );
	    SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
	    SCIPdebug( SCIProwPrint(row, NULL) );
#endif
	    SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE) );
	    SCIP_CALL( SCIPreleaseRow(scip, &row));
	    ++nGen;
	 }
      }
   }
   SCIPdebugMessage("added %d equations.\n", nGen);

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpLO)
{
   int c;
   int nGen = 0;

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

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("separating LP solution for LO constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( LOseparateTriangle(scip, consdata->n, consdata->Vars, NULL, &nGen) );
   }
   if (nGen > 0)
      *result = SCIP_SEPARATED;
   SCIPdebugMessage("separated %d cuts.\n", nGen);

   return SCIP_OKAY;
}



/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolLO)
{
   int c;
   int nGen = 0;

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

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("separating solution for LO constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( LOseparateTriangle(scip, consdata->n, consdata->Vars, NULL, &nGen) );
   }
   if (nGen > 0)
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpLO)
{
   int c;
   int nGen = 0;

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

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("enforcing lp solution for LO constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      /* we separate triangle inequalities */
      SCIP_CALL( LOseparateTriangle(scip, consdata->n, consdata->Vars, NULL, &nGen) );
      if (nGen > 0)
      {
	 *result = SCIP_SEPARATED;
	 return SCIP_OKAY;
      }
   }
   SCIPdebugMessage("all LO constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsLO)
{
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckLO)
{
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
      SCIP_VAR*** Vars;
      int i, j, k, n;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("checking LO constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->Vars != NULL );
      Vars = consdata->Vars;
      n = consdata->n;

      /* check triangle inequalities */
      for (i = 0; i < n; ++i)
      {
	 for (j = 0; j < n; ++j)
	 {
	    SCIP_Bool valIJ;
	    if (j == i)
	       continue;

	    /* the priorities should ensure that the solution is integral */
	    assert( SCIPisIntegral(scip, SCIPgetSolVal(scip, sol, Vars[i][j])) );
	    assert( SCIPisIntegral(scip, SCIPgetSolVal(scip, sol, Vars[j][i])) );
	    valIJ = SCIPisGT(scip, SCIPgetSolVal(scip, sol, Vars[i][j]), 0.5);

	    if ( valIJ && SCIPisGT(scip, SCIPgetSolVal(scip, sol, Vars[j][i]), 0.5) )
	    {
	       SCIPdebugMessage("constraint <%s> infeasible (violated equation).\n", SCIPconsGetName(cons));
	       abort();
	       *result = SCIP_INFEASIBLE;
	       return SCIP_OKAY;
	    }

	    for (k = 0; k < n; ++k)
	    {
	       SCIP_Bool valJK, valKI;
	       if (k == i || k == j)
		  continue;

	       assert( SCIPisIntegral(scip, SCIPgetSolVal(scip, sol, Vars[j][k])) );
	       assert( SCIPisIntegral(scip, SCIPgetSolVal(scip, sol, Vars[k][i])) );
	       valJK = SCIPisGT(scip, SCIPgetSolVal(scip, sol, Vars[j][k]), 0.5);
	       valKI = SCIPisGT(scip, SCIPgetSolVal(scip, sol, Vars[k][i]), 0.5);

	       /* if triangle inequality is violated */
	       if ( valIJ && valJK && valKI )
	       {
		  SCIPdebugMessage("constraint <%s> infeasible (violated triangle ineq.).\n", SCIPconsGetName(cons));
		  *result = SCIP_INFEASIBLE;
		  return SCIP_OKAY;
	       }
	    }
	 }
      }
   }
   SCIPdebugMessage("all LO constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropLO)
{
   int c;
   int nGen = 0;

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
      SCIP_VAR*** Vars;
      int i, j, k, n;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("propagating LO constraint <%s>.\n", SCIPconsGetName(cons));

      *result = SCIP_DIDNOTFIND;
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->Vars != NULL );
      Vars = consdata->Vars;
      n = consdata->n;

      /* check triangle inequalities */
      for (i = 0; i < n; ++i)
      {
	 for (j = 0; j < n; ++j)
	 {
	    if (j == i)
	       continue;

	    /* if x[i][j] == 1 then x[j][i] = 0 */
	    if ( (SCIPvarGetLbLocal(Vars[i][j]) > 0.5) )
	    {
	       SCIP_Bool infeasible, tightened;
	       SCIP_CALL( SCIPinferBinvarCons(scip, Vars[j][i], FALSE, cons, i*n + j, &infeasible, &tightened) );
	       if ( infeasible )
	       {
		  SCIPdebugMessage(" -> node infeasible.\n");
		  *result = SCIP_CUTOFF;
		  return SCIP_OKAY;
	       }
	       if ( tightened )
		  ++nGen;
	    }

	    /* if x[i][j] == 0 then x[j][i] = 1 */
	    if ( (SCIPvarGetUbLocal(Vars[i][j]) < 0.5) )
	    {
	       SCIP_Bool infeasible, tightened;
	       SCIP_CALL( SCIPinferBinvarCons(scip, Vars[j][i], TRUE, cons, i*n + j, &infeasible, &tightened) );
	       if ( infeasible )
	       {
		  SCIPdebugMessage(" -> node infeasible.\n");
		  *result = SCIP_CUTOFF;
		  return SCIP_OKAY;
	       }
	       if ( tightened )
		  ++nGen;
	    }

	    for (k = 0; k < n; ++k)
	    {
	       if (k == i || k == j)
		  continue;

	       /* if x[i][j] == 1 and x[j][k] == 1 then x[k][i] = 0 */
	       if ( (SCIPvarGetLbLocal(Vars[i][j]) > 0.5) && (SCIPvarGetLbLocal(Vars[j][k]) > 0.5))
	       {
		  SCIP_Bool infeasible, tightened;
		  SCIP_CALL( SCIPinferBinvarCons(scip, Vars[k][i], FALSE, cons, n*n + i*n*n + j*n + k, &infeasible, &tightened) );
		  if ( infeasible )
		  {
		     SCIPdebugMessage(" -> node infeasible.\n");
		     *result = SCIP_CUTOFF;
		     return SCIP_OKAY;
		  }
		  if ( tightened )
		     ++nGen;
	       }

	       /* all other implications occur with other indices i, j, k */
	    }
	 }
      }
   }
   if (nGen > 0)
      *result = SCIP_REDUCEDDOM;
   SCIPdebugMessage("propagated %d domains.\n", nGen);

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
#define consPresolLO NULL


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropLO)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** Vars;
   int n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIPdebugMessage("Propagation resolution of constraint <%s>.\n", SCIPconsGetName(cons));
   *result = SCIP_DIDNOTFIND;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->Vars != NULL );

   n = consdata->n;
   Vars = consdata->Vars;

   assert( 0 <= inferinfo && inferinfo < n*n + n*n*n );

   /* if the conflict came from an equation */
   if ( inferinfo < (n*n) )
   {
      int index1, index2;
      index1 = inferinfo/n;
      index2 = inferinfo % n;
      assert( 0 <= index1 && index1 < n );
      assert( 0 <= index2 && index2 < n );
      assert( Vars[index2][index1] == infervar );

      /* if the variable was fixed to 0 */
      if ( SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < 0.5 )
      {
	 SCIPdebugMessage(" -> reason for x[%d][%d] == 0 was x[%d][%d] = 1.\n", index2, index1, index1, index2);
	 /* the reason was that x[i][j] was fixed to 1 */
	 SCIP_CALL( SCIPaddConflictLb(scip, Vars[index1][index2], bdchgidx) );
	 *result = SCIP_SUCCESS;
	 return SCIP_OKAY;
      }

      /* if the variable was fixed to 1 */
      if ( SCIPvarGetLbAtIndex(infervar, bdchgidx, FALSE) < 0.5 && SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) > 0.5 )
      {
	 SCIPdebugMessage(" -> reason for x[%d][%d] == 1 was x[%d][%d] = 0.\n", index2, index1, index1, index2);
	 /* the reason was that x[i][j] was fixed to 0 */
	 SCIP_CALL( SCIPaddConflictUb(scip, Vars[index1][index2], bdchgidx) );
	 *result = SCIP_SUCCESS;
	 return SCIP_OKAY;
      }
   }
   else
   {
      /* otherwise the conflict came from a triangle inequality */
      int index1, index2, index3;
      index1 = (inferinfo - n*n)/(n*n);
      index2 = (inferinfo - n*n - index1 * n*n)/n;
      index3 = (inferinfo - n*n) % n;

      assert( 0 <= index1 && index1 < n );
      assert( 0 <= index2 && index2 < n );
      assert( 0 <= index3 && index3 < n );
      assert( index1 != index2 && index2 != index3 && index1 != index3 );
      assert( Vars[index3][index1] == infervar );

      /* the variable should have been fixed to 0 */
      assert( SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < 0.5 );

      /* the reason was that x[index1][index2] and x[index2][index3] were fixed to 1 */
      SCIPdebugMessage(" -> reason for x[%d][%d] == 0 was x[%d][%d] = x[%d][%d] = 0.\n", index3, index1, index1, index2, index2, index3);
      SCIP_CALL( SCIPaddConflictLb(scip, Vars[index1][index2], bdchgidx) );
      SCIP_CALL( SCIPaddConflictLb(scip, Vars[index2][index3], bdchgidx) );
      *result = SCIP_SUCCESS;
   }
   
   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockLO)
{
   int i, j;
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** Vars;
   int n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMessage("Locking LO constraint <%s>.\n", SCIPconsGetName(cons));

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->Vars != NULL );
   n = consdata->n;
   Vars = consdata->Vars;

   for (i = 0; i < n; ++i)
   {
      for (j = 0; j < n; ++j)
      {
	 if (i == j)
	    continue;

	 /* the constaint may be violated in any way */
	 SCIPaddVarLocks(scip, Vars[i][j], nlockspos + nlocksneg, nlockspos + nlocksneg);
      }
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveLO NULL

/** constraint deactivation notification method of constraint handler */
#define consDeactiveLO NULL

/** constraint enabling notification method of constraint handler */
#define consEnableLO NULL

/** constraint disabling notification method of constraint handler */
#define consDisableLO NULL

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintLO)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** Vars;
   int i, j, n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->Vars != NULL );
   n = consdata->n;
   Vars = consdata->Vars;

   SCIPinfoMessage(scip, file, "LO[");
   for (i = 0; i < n; ++i)
   {
      if ( i > 0 )
	 SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "(");
      for (j = 0; j < n; ++j)
      {
	 if (j == i)
	    continue;
	 if ( j > 0 && (i > 0 || j > 1) )
	    SCIPinfoMessage(scip, file, ",");
	 SCIPinfoMessage(scip, file, "%s", SCIPvarGetName(Vars[i][j]));
      }
      SCIPinfoMessage(scip, file, ")");
   }
   SCIPinfoMessage(scip, file, "]\n");

   return SCIP_OKAY;
}





/** creates the handler for LO constraints and includes it in SCIP */
SCIP_RETCODE LOincludeConshdlr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeLO, consInitLO, consExitLO,
         consInitpreLO, consExitpreLO, consInitsolLO, consExitsolLO,
         consDeleteLO, consTransLO, consInitlpLO,
         consSepalpLO, consSepasolLO, consEnfolpLO, consEnfopsLO, consCheckLO,
         consPropLO, consPresolLO, consRespropLO, consLockLO,
         consActiveLO, consDeactiveLO,
         consEnableLO, consDisableLO,
         consPrintLO,
         NULL) );

   return SCIP_OKAY;
}


/** creates and captures a LO constraint */
SCIP_RETCODE LOcreateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   n,                  /**< number of elements */
   SCIP_VAR***           Vars,               /**< n x n matrix of binary variables */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local               /**< is constraint only valid locally? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int i, j;

   /* find the LO constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if (conshdlr == NULL)
   {
      SCIPerrorMessage("LO constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->n = n;
   SCIPallocBlockMemoryArray(scip, &consdata->Vars, n);
   for (i = 0; i < n; ++i)
   {
      SCIPallocBlockMemoryArray(scip, &(consdata->Vars[i]), n);
      for (j = 0; j < n; ++j)
      {
	 if (j == i)
	    continue;
	 assert( Vars[i][j] != NULL );
	 consdata->Vars[i][j] = Vars[i][j];
      }
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
			     local, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
