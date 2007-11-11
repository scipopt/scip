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
#pragma ident "@(#) $Id: cons_indicator.c,v 1.1 2007/11/11 17:29:58 bzfpfets Exp $"

/**@file   cons_indicator.c
 * @brief  constraint handler for indicator constraints
 * @author Marc Pfetsch
 *
 * See also the comments in the .h file.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"
#include <string.h>


/* constraint handler properties */
#define CONSHDLR_NAME          "indicator"
#define CONSHDLR_DESC          "indicator constraint handler"
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY      -10 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* event handler properties */
#define EVENTHDLR_NAME         "indicator"
#define EVENTHDLR_DESC         "bound change event handler for indicator constraints"





/** constraint data for indicator constraints */
struct SCIP_ConsData
{
   SCIP_VAR*   binvar;             /**< binary variable for indicator constraint */
   SCIP_VAR*   slackvar;           /**< slack variable of inequality of indicator constraint */
   SCIP_CONS*  lincons;            /**< linear constraint corresponding to indicator constraint */
   int         nFixedNonzero;      /**< number of variables fixed to be nonzero */
};

/** indicator constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR* eventhdlr;     /**< event handler for bound change events */
   SCIP_LPI*   altLP;             /**< alternative LP for cut separation */
   int         nCols;             /**< number of columns in alternative LP = number of indicator constraints */
   int         nvars;             /**< number of variables in all linear constraints */
   SCIP_VAR**  vars;              /**< variables appearing in linear constraints */
   SCIP_HASHMAP* varHash;         /**< hash map from variable to column index in alternative LP */
};


/* ---------------------------- separation methods ----------------------*/

/** hash key retrieval function for variables */
static
SCIP_DECL_HASHGETKEY(hashGetKeyVar)
{
   return elem;
}

/** returns TRUE iff the indices of both variables are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqVar)
{
   if ( key1 == key2 )
      return TRUE;
   return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVar)
{
   return SCIPvarGetIndex((SCIP_VAR*) key);
}

/** coount variables appearing in linear constraints corresponding to the indicator constraints */
static
SCIP_RETCODE countLinearVars(
      SCIP* scip,               /**< SCIP pointer */
      int nconss,               /**< number of constraints */
      SCIP_CONS** conss,        /**< indicator constraints */
      int* nvars                /**< number of variables */
)
{
   int c;
   SCIP_HASHTABLE* varHash;

   assert( scip != NULL );
   assert( conss != NULL );
   assert( nvars != NULL );

   /* init */
   SCIP_CALL( SCIPhashtableCreate(&varHash, SCIPblkmem(scip), 10 * SCIPgetNVars(scip), hashGetKeyVar, hashKeyEqVar, hashKeyValVar) );

   *nvars = 0;
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* lincons;
      SCIP_VAR** linvars;
      SCIP_VAR* slackvar;
      int nlinvars;
      int v;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      lincons = consdata->lincons;
      assert( lincons != NULL );

      slackvar = consdata->slackvar;
      linvars = SCIPgetVarsLinear(scip, lincons);
      nlinvars = SCIPgetNVarsLinear(scip, lincons);

      for (v = 0; v < nlinvars; ++v)
      {
	 SCIP_VAR* var;
	 var = linvars[v];
	 if ( var != slackvar && ! SCIPhashtableExists(varHash, (void*) var) )
	 {
	    ++(*nvars);
	    SCIP_CALL( SCIPhashtableInsert(varHash, (void*) var) );
	 }
      }
   }
   SCIPhashtableFree(&varHash);

   return SCIP_OKAY;
}


/** collect all variables appearing in linear constraints corresponding to the indicator constraints */
static
SCIP_RETCODE collectAllLinearVars(
      SCIP* scip,               /**< SCIP pointer */
      int nconss,               /**< number of constraints */
      SCIP_CONS** conss,        /**< indicator constraints */
      int* nvars,               /**< number of variables */
      SCIP_VAR** vars,          /**< variables */
      SCIP_HASHMAP* varHash,    /**< hash map storing the variable */
      int* nnonz                /**< total number of nonzeros */
)
{
   int c;

   assert( scip != NULL );
   assert( conss != NULL );
   assert( vars != NULL );
   assert( nvars != NULL );
   assert( nnonz != NULL );

   /* init */
   *nnonz = 0;
   *nvars = 0;
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* lincons;
      SCIP_VAR** linvars;
      SCIP_VAR* slackvar;
      int nlinvars;
      int v;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      lincons = consdata->lincons;
      assert( lincons != NULL );

      slackvar = consdata->slackvar;
      linvars = SCIPgetVarsLinear(scip, lincons);
      nlinvars = SCIPgetNVarsLinear(scip, lincons);
      *nnonz += nlinvars;

      for (v = 0; v < nlinvars; ++v)
      {
	 SCIP_VAR* var;
	 var = linvars[v];
	 if ( var != slackvar && ! SCIPhashmapExists(varHash, var) )
	 {
	    vars[*nvars] = var;
	    SCIP_CALL( SCIPhashmapInsert(varHash, var, (void*) *nvars) );
	    ++(*nvars);
	 }
      }
   }

   return SCIP_OKAY;
}


/** generate the alternative LP corresponding to all currently present indicator constraints */
static
SCIP_RETCODE generateAlternativeLP(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_CONSHDLR* conshdlr,  /**< constraint handler */
      int nconss,               /**< number of constraints */
      SCIP_CONS** conss         /**< indicator constraints */
      )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nvars = 0;
   int nnonz = 0;
   SCIP_Real sumLastRow = 0.0;
   int lastnnonz = 0;
   int i, c, cnt;

   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* obj;
   int* matbeg;
   int* matind;
   SCIP_Real* matval;
   SCIP_Real* lhs;
   SCIP_Real* rhs;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );

   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   assert( conshdlrdata->altLP == NULL );
   assert( conshdlrdata->vars == NULL );
   assert( conshdlrdata->varHash == NULL );

   /* create hash map of variables */
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->varHash, SCIPblkmem(scip), 10 * SCIPgetNVars(scip)) );

   /* init array of linear variables */
   SCIP_CALL( countLinearVars(scip, nconss, conss, &conshdlrdata->nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrdata->vars), conshdlrdata->nvars) );

   /* get linear variables */
   SCIP_CALL( collectAllLinearVars(scip, nconss, conss, &nvars, conshdlrdata->vars, conshdlrdata->varHash, &nnonz) );
   assert( nvars == conshdlrdata->nvars );

   /* create alternative LP */
   SCIP_CALL( SCIPlpiCreate(&conshdlrdata->altLP, "altLP", SCIP_OBJSEN_MINIMIZE) );

   /* reserve space */
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matbeg, conshdlrdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matind, nnonz + conshdlrdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matval, nnonz + conshdlrdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, conshdlrdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, conshdlrdata->nvars) );
   conshdlrdata->nCols = nconss;

   /* set up rows of alternative LP */
   cnt = 0;
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* lincons;
      SCIP_VAR** linvars;
      SCIP_Real* linvals;
      SCIP_VAR* slackvar;
      SCIP_Real val;
      SCIP_Real sign = 1.0;
      int nlinvars;
      int v;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      lincons = consdata->lincons;
      assert( lincons != NULL );

      slackvar = consdata->slackvar;
      linvars = SCIPgetVarsLinear(scip, lincons);
      linvals = SCIPgetValsLinear(scip, lincons);
      nlinvars = SCIPgetNVarsLinear(scip, lincons);

      val = SCIPgetRhsLinear(scip, lincons);
      if ( val == SCIPinfinity(scip) )
      {
	 val = SCIPgetLhsLinear(scip, lincons);
	 assert( val > -SCIPinfinity(scip) );
	 sign = -1.0;
      }

      obj[c] = 1.0;
      lb[c] = 0.0;
      ub[c] = SCIPlpiInfinity(conshdlrdata->altLP);
      matbeg[c] = cnt;

      for (v = 0; v < nlinvars; ++v)
      {
	 SCIP_VAR* var;
	 var = linvars[v];

	 if ( var != slackvar )
	 {
	    assert( SCIPhashmapExists(conshdlrdata->varHash, var) );
	    matind[cnt] = (int) SCIPhashmapGetImage(conshdlrdata->varHash, var);
	    matval[cnt] = sign * linvals[v];
	    ++cnt;
	 }
      }

      /* last row: */
      if (! SCIPisFeasZero(scip, val) )
      {
	 sumLastRow += fabs(val);
	 ++lastnnonz;
	 matind[cnt] = conshdlrdata->nvars;
	 matval[cnt] = sign * val;
	 ++cnt;
      }
   }
   assert( cnt <= nnonz + conshdlrdata->nvars );

   /* set up rows */
   for (i = 0; i < conshdlrdata->nvars; ++i)
   {
      lhs[i] = 0.0;
      rhs[i] = 0.0;
   }
   /* last row: */
   lhs[conshdlrdata->nvars] = -sumLastRow/( (double) lastnnonz);
   rhs[conshdlrdata->nvars] = lhs[conshdlrdata->nvars];

   SCIP_CALL( SCIPlpiLoadColLP(conshdlrdata->altLP, SCIP_OBJSEN_MINIMIZE, nconss, obj, lb, ub, NULL, conshdlrdata->nvars+1,
			       lhs, rhs, NULL, cnt, matbeg, matind, matval) );

   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArray(scip, &matval);
   SCIPfreeBufferArray(scip, &matind);
   SCIPfreeBufferArray(scip, &matbeg);
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);

   SCIP_CALL( SCIPlpiWriteLP(conshdlrdata->altLP, "alt.lp") );

   /*
   CHECK_PARAM_OK( SCIPlpiSetIntpar(altLP_, SCIP_LPPAR_FROMSCRATCH, altFromScratch()) );
   CHECK_PARAM_OK( SCIPlpiSetIntpar(altLP_, SCIP_LPPAR_PRESOLVING, altPreproc()) );
   CHECK_PARAM_OK( SCIPlpiSetIntpar(altLP_, SCIP_LPPAR_SCALING, altScaling()) );
   CHECK_PARAM_OK( SCIPlpiSetIntpar(altLP_, SCIP_LPPAR_LPINFO, altLPInfo()) );
   */

   return SCIP_OKAY;
}








/* ---------------------------- constraint handler local methods ----------------------*/

/** propagate indicator constraint */
static
SCIP_RETCODE propIndicator(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_CONS* cons,          /**< constraint */
      SCIP_CONSDATA* consdata,  /**< constraint data */
      SCIP_Bool* cutoff,        /**< whether a cutoff happend */
      int* nGen                 /**< number of domain changes */
      )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( cutoff != NULL );
   assert( nGen != NULL );

   *cutoff = FALSE;

   /* if both slackvar and binvar are fixed to be nonzero */
   if ( consdata->nFixedNonzero > 1 )
   {
      SCIPdebugMessage("the node is infeasible, both the slackvariable and the binary variable are fixed to be nonzero.\n");
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if exactly one of the variables is fixed to be nonzero */
   if ( consdata->nFixedNonzero == 1 )
   {
      SCIP_Bool infeasible, tightened;

      /* if binvar is fixed to be nonzero */
      if ( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
      {
	 SCIPdebugMessage("binary variable <%s> is fixed to be nonzero, fixing slack variable <%s> to 0.\n",
			  SCIPvarGetName(consdata->binvar), SCIPvarGetName(consdata->slackvar));

	 /* fix slack variable to 0 */
	 assert( SCIPvarGetStatus(consdata->slackvar) != SCIP_VARSTATUS_MULTAGGR );
	 assert( SCIPvarGetStatus(consdata->slackvar) != SCIP_VARSTATUS_AGGREGATED );

	 SCIP_CALL( SCIPinferVarUbCons(scip, consdata->slackvar, 0.0, cons, 0, &infeasible, &tightened) );
	 assert( ! infeasible );
	 if ( tightened )
	    ++(*nGen);
      }

      /* if slackvar is fixed to be nonzero */
      if ( SCIPisPositive(scip, SCIPvarGetLbLocal(consdata->slackvar)) )
      {
	 SCIPdebugMessage("slack variable <%s> is fixed to be nonzero, fixing binary variable <%s> to 0.\n",
			  SCIPvarGetName(consdata->slackvar), SCIPvarGetName(consdata->binvar));

	 /* fix binary variable to 0 */
	 SCIP_CALL( SCIPinferVarUbCons(scip, consdata->binvar, 0.0, cons, 1, &infeasible, &tightened) );
	 assert( ! infeasible );
	 if ( tightened )
	    ++(*nGen);
      }

      /* reset constraint age counter */
      if ( *nGen > 0 )
	 SCIP_CALL( SCIPresetConsAge(scip, cons) );

      /* delete constraint locally */
      assert( !SCIPconsIsModifiable(cons) );
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
   }

   return SCIP_OKAY;
}




/** enforcement method
 *
 *  We check whether the current solution is feasible, i.e., if binvar = 1
 *  implies that slackvar = 0. If not, we branch as follows:
 *
 *  In one branch we fix binvar = 1 and slackvar = 0. In the other branch
 *  we fix binvar = 0 and leave slackvar unchanged.
 */
static
SCIP_RETCODE enforceIndicator(
	 SCIP* scip,               /**< SCIP pointer */
	 SCIP_CONS* cons,          /**< constraint */
	 SCIP_CONSDATA* consdata,  /**< constraint data */
	 SCIP_RESULT* result       /**< result */
	 )
{
   SCIP_Bool cutoff;
   SCIP_NODE* node1;
   SCIP_NODE* node2;
   SCIP_VAR* slackvar;
   SCIP_VAR* binvar;
   int cnt = 0;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( result != NULL );

   /* first perform propagation (it might happen that standard propagation is turned off) */
   SCIP_CALL( propIndicator(scip, cons, consdata, &cutoff, &cnt) );
   SCIPdebugMessage("propagation in enforcing (cutoff: %d, domain reductions: %d).\n", cutoff, cnt);
   if ( cutoff )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }
   if ( cnt > 0 )
   {
      *result = SCIP_REDUCEDDOM;
      return SCIP_OKAY;
   }

   /* if constraint is infeasible */
   binvar = consdata->binvar;
   slackvar = consdata->slackvar;
   if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, binvar)) &&
	! SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, slackvar)) )
   {
      /* binary variable is not fixed - otherwise we would not be infeasible */
      assert( SCIPvarGetLbLocal(binvar) < 0.5 && SCIPvarGetUbLocal(binvar) > 0.5 );

      /* create branches */
      SCIPdebugMessage("Creating two branches.\n");

      /* node1: binvar = 1, slackvar = 0 */
      SCIP_CALL( SCIPcreateChild(scip, &node1, 0.0, SCIPcalcChildEstimate(scip, binvar, 1.0) ) );

      if ( ! SCIPisFeasEQ(scip, SCIPvarGetLbLocal(binvar), 1.0) )
	 SCIP_CALL( SCIPchgVarLbNode(scip, node1, binvar, 1.0) );

      if ( ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(slackvar)) )
	 SCIP_CALL( SCIPchgVarUbNode(scip, node1, slackvar, 0.0) );

      /* node2: binvar = 0, no restriction on slackvar */
      SCIP_CALL( SCIPcreateChild(scip, &node2, 0.0, SCIPcalcChildEstimate(scip, binvar, 0.0) ) );

      if ( ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(binvar)) )
	 SCIP_CALL( SCIPchgVarUbNode(scip, node2, binvar, 0.0) );

      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}





/* ---------------------------- constraint handler callback methods ----------------------*/

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeIndicator)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->vars != NULL )
   {
      assert( conshdlrdata->varHash != NULL );
      assert( conshdlrdata->altLP != NULL );

      SCIPfreeMemoryArray(scip, &conshdlrdata->vars);
      SCIPhashmapFree(&conshdlrdata->varHash);
      SCIP_CALL( SCIPlpiFree(&conshdlrdata->altLP) );
   }

   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}



/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolIndicator)
{
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* generat alternative LP */
   SCIP_CALL( generateAlternativeLP(scip, conshdlr, nconss, conss) );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* lincons;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      SCIPdebugMessage("Initializing indicator constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* if not happend already, get transformed linear constraint */
      lincons = consdata->lincons;
      SCIP_CALL( SCIPgetTransformedCons(scip, lincons, &consdata->lincons) );

      /* init separation routines */
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolIndicator)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->vars != NULL )
   {
      assert( conshdlrdata->varHash != NULL );
      assert( conshdlrdata->altLP != NULL );

      SCIPfreeMemoryArray(scip, &conshdlrdata->vars);
      SCIPhashmapFree(&conshdlrdata->varHash);
      SCIP_CALL( SCIPlpiFree(&conshdlrdata->altLP) );
   }

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      SCIPdebugMessage("Exiting indicator constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* free separation routines */
   }
   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteIndicator)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMessage("Deleting indicator constraint <%s>.\n", SCIPconsGetName(cons) );

   /* drop events on transfromed variables */
   if ( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get constraint handler data */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );
      assert( conshdlrdata->eventhdlr != NULL );

      SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->binvar, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
				  (SCIP_EVENTDATA*)*consdata, -1) );

      SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->slackvar, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
				  (SCIP_EVENTDATA*)*consdata, -1) );
   }
   else
   {
      /* release linear constraint and slack variable only for nontransformed constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &(*consdata)->lincons) );
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->slackvar) );
   }

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransIndicator)
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   char s[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   assert( conshdlrdata->eventhdlr != NULL );

   SCIPdebugMessage("Transforming indicator constraint: <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->lincons != NULL );
   assert( sourcedata->binvar != NULL );
   assert( sourcedata->slackvar != NULL );

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->binvar, &(consdata->binvar)) );
   SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->slackvar, &(consdata->slackvar)) );
   assert( consdata->binvar != NULL );
   assert( consdata->slackvar != NULL );

   /* if binary variable is fixed to be nonzero */
   consdata->nFixedNonzero = 0;
   if ( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
      ++(consdata->nFixedNonzero);

   /* if slack variable is fixed to be nonzero */
   if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(consdata->slackvar) > 0.5 ) )
      ++(consdata->nFixedNonzero);

   /* try to obtain transformed linear constraint - this may fail because it might be transformed later */
   SCIP_CALL( SCIPgetTransformedCons(scip, sourcedata->lincons, &consdata->lincons) );
   /* if failed store original constraint */
   if ( consdata->lincons == NULL )
      consdata->lincons = sourcedata->lincons;

   /* create transformed constraint with the same flags */
   snprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
			     SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
			     SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
			     SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
			     SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
			     SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* catch bound change events on variables */
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->binvar, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
				(SCIP_EVENTDATA*)consdata, NULL) );
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->slackvar, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
				(SCIP_EVENTDATA*)consdata, NULL) );

#ifdef SCIP_DEBUG
   if ( consdata->nFixedNonzero > 0 )
   {
      SCIPdebugMessage("constraint <%s> has %d variables fixed to be nonzero.\n", SCIPconsGetName(*targetcons),
		       consdata->nFixedNonzero );
   }
#endif

   return SCIP_OKAY;
}




/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolIndicator)
{  /*lint --e{715}*/
   int c;
   int oldnfixedvars = 0;
   int oldndelconss = 0;
   int removedvars = 0;
   SCIP_EVENTHDLR* eventhdlr;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;
   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;

   /* get constraint handler data */
   assert( SCIPconshdlrGetData(conshdlr) != NULL );
   eventhdlr = SCIPconshdlrGetData(conshdlr)->eventhdlr;
   assert( eventhdlr != NULL );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->lincons != NULL );
      assert( consdata->binvar != NULL );
      assert( consdata->slackvar != NULL );
      assert( ! SCIPconsIsModifiable(cons) );

      SCIPdebugMessage("Presolving indicator constraint <%s>.\n", SCIPconsGetName(cons) );

      *result = SCIP_DIDNOTFIND;

      /* only run if sucess if possible */
      if ( nrounds == 0 || nnewfixedvars > 0 || nnewaggrvars > 0 || *nfixedvars > oldnfixedvars )
      {
	 SCIP_Bool infeasible, fixed;

	 /* if the binary variable if fixed to nonzero */
	 if ( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
	 {
	    /* if slack variable is fixed to nonzero, we are infeasible */
	    if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(consdata->slackvar)) )
	    {
	       SCIPdebugMessage("The problem is infeasible: binary and slack variable are fixed to be nonzero.\n");
	       *result = SCIP_CUTOFF;
	       return SCIP_OKAY;
	    }

	    /* otherwise fix slack variable to 0 */
	    SCIPdebugMessage("Fix slack variable to 0 and delete constraint.\n");
	    SCIP_CALL( SCIPfixVar(scip, consdata->slackvar, 0.0, &infeasible, &fixed) );
	    assert( ! infeasible );
	    if ( fixed )
	       ++(*nfixedvars);

	    /* delete constraint */
	    assert( ! SCIPconsIsModifiable(cons) );
	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);
	    *result = SCIP_SUCCESS;
	    continue;
	 }

	 /* if the slack variable if fixed to nonzero */
	 if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(consdata->slackvar)) )
	 {
	    /* if binary variable is fixed to nonzero, we are infeasible */
	    if ( SCIPvarGetLbLocal(consdata->slackvar) > 0.5 )
	    {
	       SCIPdebugMessage("The problem is infeasible: binary and slack variable are fixed to be nonzero.\n");
	       *result = SCIP_CUTOFF;
	       return SCIP_OKAY;
	    }

	    /* otherwise fix binary variable to 0 */
	    SCIPdebugMessage("Fix binary variable to 0 and delete constraint.\n");
	    SCIP_CALL( SCIPfixVar(scip, consdata->binvar, 0.0, &infeasible, &fixed) );
	    assert( ! infeasible );
	    if ( fixed )
	       ++(*nfixedvars);

	    /* delete constraint */
	    assert( ! SCIPconsIsModifiable(cons) );
	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);
	    *result = SCIP_SUCCESS;
	    continue;
	 }
      }
   }

   SCIPdebugMessage("presolving fixed %d variables, removed %d variables, and deleted %d constraints.\n",
		    *nfixedvars - oldnfixedvars, removedvars, *ndelconss - oldndelconss);

   return SCIP_OKAY;
}




/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpIndicator)
{
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      SCIPdebugMessage("Checking for initial rows for indicator constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* check separation routines */
   }

   return SCIP_OKAY;
}



/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpIndicator)
{  /*lint --e{715}*/
   int c;
   int nGen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      *result = SCIP_DIDNOTFIND;
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMessage("Separating inequalities for indicator constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* start separation */
   }
   SCIPdebugMessage("Separated %d indicator constraints.\n", nGen);
   if ( nGen > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolIndicator)
{  /*lint --e{715}*/
   int c;
   int nGen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      *result = SCIP_DIDNOTFIND;
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMessage("Separating solution for indicator constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* start separation */
   }
   SCIPdebugMessage("Separated %d indicator constraints.\n", nGen);
   if ( nGen > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpIndicator)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMessage("Enforcing indicator constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      SCIP_CALL( enforceIndicator(scip, conss[c], consdata, result) );

      if ( *result != SCIP_FEASIBLE )
	 return SCIP_OKAY;
   }
   SCIPdebugMessage("All indicator constraints are feasible.\n");

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsIndicator)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMessage("Enforcing indicator constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      SCIP_CALL( enforceIndicator(scip, conss[c], consdata, result) );

      if ( *result != SCIP_FEASIBLE )
	 return SCIP_OKAY;
   }
   SCIPdebugMessage("All indicator constraints are feasible.\n");

   return SCIP_OKAY;
}




/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckIndicator)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      *result = SCIP_DIDNOTFIND;
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMessage("Checking indicator constraint <%s>.\n", SCIPconsGetName(conss[c]) );
      assert( consdata->binvar != NULL );
      assert( consdata->slackvar != NULL );

      if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->binvar)) &&
	   ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->slackvar)) )
      {
	 SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
	 *result = SCIP_INFEASIBLE;
	 return SCIP_OKAY;
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropIndicator)
{  /*lint --e{715}*/
   int c;
   int nGen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );
   *result = SCIP_DIDNOTRUN;

   assert( SCIPisTransformed(scip) );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      SCIP_Bool cutoff;

      *result = SCIP_DIDNOTFIND;
      assert( conss[c] != NULL );
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      SCIPdebugMessage("Propagating indicator constraint <%s>.\n", SCIPconsGetName(cons) );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( propIndicator(scip, cons, consdata, &cutoff, &nGen) );
      if ( cutoff )
      {
	 *result = SCIP_CUTOFF;
	 return SCIP_OKAY;
      }
   }
   SCIPdebugMessage("Propagated %d domains.\n", nGen);
   if ( nGen > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}




/** propagation conflict resolving method of constraint handler
 *
 *  We check which bound changes were the reason for infeasibility. We
 *  use that @a inferinfo is 0 if the binary variable has bounds that
 *  fix it to be nonzero (these bounds are the reason). Likewise
 *  @inferinfo is 1 if the slack variable * has bounds that fix it to
 *  be nonzero.
 */
static
SCIP_DECL_CONSRESPROP(consRespropIndicator)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;
   SCIPdebugMessage("Propagation resolution method of indicator constraint <%s>.\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( inferinfo == 0 || inferinfo == 1 );

   /* if the binary variable was the reason */
   if ( inferinfo == 0 )
   {
      assert( SCIPvarGetLbAtIndex(consdata->binvar, bdchgidx, FALSE) > 0.5 );
      assert( infervar != consdata->binvar );

      SCIP_CALL( SCIPaddConflictLb(scip, consdata->binvar, bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   /* if the slack variable was the reason */
   if ( inferinfo == 1 )
   {
      assert( SCIPisFeasPositive(scip, SCIPvarGetLbAtIndex(consdata->slackvar, bdchgidx, FALSE)) );
      assert( infervar != consdata->slackvar );

      SCIP_CALL( SCIPaddConflictLb(scip, consdata->slackvar, bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}




/** variable rounding lock method of constraint handler
 *
 *  The up-rounding of the binary and slack variable may violate the
 *  constraint.
 */
static
SCIP_DECL_CONSLOCK(consLockIndicator)
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->binvar != NULL );
   assert( consdata->slackvar != NULL );

   SCIPdebugMessage("Locking constraint <%s>.\n", SCIPconsGetName(cons));

   SCIP_CALL( SCIPaddVarLocks(scip, consdata->binvar, nlocksneg, nlockspos) );
   SCIP_CALL( SCIPaddVarLocks(scip, consdata->slackvar, nlocksneg, nlockspos) );

   return SCIP_OKAY;
}



/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintIndicator)
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->binvar != NULL );
   assert( consdata->slackvar != NULL );
   assert( consdata->lincons != NULL );

   SCIPinfoMessage(scip, file, "[%s] <%s>: Indicator(", CONSHDLR_NAME, SCIPconsGetName(cons));
   SCIPinfoMessage(scip, file, "%s = 1", SCIPvarGetName(consdata->binvar));
   SCIPinfoMessage(scip, file, " -> %s = 0)\n", SCIPvarGetName(consdata->slackvar));

   return SCIP_OKAY;
}






/** constraint activation notification method of constraint handler */
#define consActiveIndicator NULL

/** constraint deactivation notification method of constraint handler */
#define consDeactiveIndicator NULL

/** constraint enabling notification method of constraint handler */
#define consEnableIndicator NULL

/** constraint disabling notification method of constraint handler */
#define consDisableIndicator NULL

/** initialization method of constraint handler (called after problem was transformed) */
#define consInitIndicator NULL

/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitIndicator NULL

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreIndicator NULL

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreIndicator NULL





/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * We update the number of variables fixed to be nonzero
 */
static
SCIP_DECL_EVENTEXEC(eventExecIndicator)
{
   SCIP_EVENTTYPE eventtype;
   SCIP_CONSDATA* consdata;
   SCIP_Real oldbound, newbound;

   assert( eventhdlr != NULL );
   assert( eventdata != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );
   assert( event != NULL );

   consdata = (SCIP_CONSDATA*)eventdata;
   assert( consdata != NULL );
   assert( 0 <= consdata->nFixedNonzero && consdata->nFixedNonzero <= 2 );

   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);

   eventtype = SCIPeventGetType(event);
   switch ( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisFeasPositive(scip, oldbound) && SCIPisFeasPositive(scip, newbound) )
	 ++(consdata->nFixedNonzero);
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisFeasNegative(scip, oldbound) && SCIPisFeasNegative(scip, newbound) )
	 ++(consdata->nFixedNonzero);
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisFeasPositive(scip, oldbound) && ! SCIPisFeasPositive(scip, newbound) )
	 --(consdata->nFixedNonzero);
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisFeasNegative(scip, oldbound) && ! SCIPisFeasNegative(scip, newbound) )
	 --(consdata->nFixedNonzero);
      break;
   default:
      SCIPerrorMessage("invalid event type.\n");
      return SCIP_INVALIDDATA;
   }
   assert( 0 <= consdata->nFixedNonzero && consdata->nFixedNonzero <= 2 );

   SCIPdebugMessage("changed bound of variable <%s> from %f to %f (nFixedNonzero: %d).\n",
		    SCIPvarGetName(SCIPeventGetVar(event)),
		    oldbound, newbound, consdata->nFixedNonzero);

   return SCIP_OKAY;
}




/* ---------------- Constraint specific interface methods ---------------- */

/** creates the handler for indicator constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrIndicator(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
			 NULL, NULL, NULL, NULL, NULL, NULL, eventExecIndicator, NULL) );

   /* create constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   /* get event handler for bound change events */
   conshdlrdata->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if ( conshdlrdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for indicator constraints not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   conshdlrdata->altLP = NULL;
   conshdlrdata->nCols = 0;
   conshdlrdata->nvars = 0;
   conshdlrdata->vars = 0;
   conshdlrdata->varHash = 0;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
			  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
			  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
			  CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
			  consFreeIndicator, consInitIndicator, consExitIndicator,
			  consInitpreIndicator, consExitpreIndicator, consInitsolIndicator, consExitsolIndicator,
			  consDeleteIndicator, consTransIndicator, consInitlpIndicator, consSepalpIndicator,
			  consSepasolIndicator, consEnfolpIndicator, consEnfopsIndicator, consCheckIndicator,
			  consPropIndicator, consPresolIndicator, consRespropIndicator, consLockIndicator,
			  consActiveIndicator, consDeactiveIndicator, consEnableIndicator, consDisableIndicator,
			  consPrintIndicator, conshdlrdata) );

   return SCIP_OKAY;
}



/** creates and captures a indicator constraint
 *
 *  We set the constraint to not be modifable. If the weights are non
 *  NULL, the variables are ordered according to these weights (in
 *  ascending order).
 */
SCIP_RETCODE SCIPcreateConsIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< binary indicator variable */
   int                   nvars,              /**< number of variables in the inequality */
   SCIP_VAR**            vars,               /**< array with variables of inequality */
   SCIP_Real*            vals,               /**< values of variables in inequality */
   SCIP_Real             rhs,                /**< rhs of the inequality */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? Usually set to TRUE. */
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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? Usually set to FALSE. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* lincons;
   SCIP_VAR* slackvar;
   SCIP_Bool modifiable = FALSE;
   char s[SCIP_MAXSTRLEN];

   /* find the indicator constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("<%s> constraint handler not found\n", CONSHDLR_NAME);
      return SCIP_PLUGINNOTFOUND;
   }

   if ( SCIPvarGetType(binvar) != SCIP_VARTYPE_BINARY )
   {
      SCIPerrorMessage("indicator variable is not binary.\n");
      return SCIP_ERROR;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->binvar = binvar;
   consdata->nFixedNonzero = 0;

   /* create slack variable */
   snprintf(s, SCIP_MAXSTRLEN, "indslack_%s", name);
   SCIP_CALL( SCIPcreateVar(scip, &slackvar, s, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE,
			    NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, slackvar) );
   consdata->slackvar = slackvar;

   /* create linear constraint */
   snprintf(s, SCIP_MAXSTRLEN, "indlin_%s", name);

   /* the constraint is inital, enforced, separated, and checked */
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, s, nvars, vars, vals, -SCIPinfinity(scip), rhs,
				   TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   /* add slack variable */
   SCIP_CALL( SCIPaddCoefLinear(scip, lincons, slackvar, -1.0) );

   SCIP_CALL( SCIPaddCons(scip, lincons) );
   consdata->lincons = lincons;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
			     local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}
