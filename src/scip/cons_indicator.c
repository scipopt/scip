/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_indicator.c,v 1.19 2008/06/22 18:49:26 bzfpfets Exp $"
//#define SCIP_DEBUG
//#define SCIP_OUTPUT

/**@file   cons_indicator.c
 * @brief  constraint handler for indicator constraints
 * @author Marc Pfetsch
 *
 * An indicator constraint is given by a binary variable \f$z\f$ and an inequality \f$ax \leq
 * b\f$. It states that if \f$z = 1\f$ then \f$ax \leq b\f$ holds.
 *
 * This constraint is handled by adding a slack variable \f$s\f$: \f$ax - s \leq b\f$ with \f$s \leq
 * 0\f$. The constraint is enforced by fixing \f$s\f$ to 0 if \f$z = 1\f$.
 *
 * @note: The constraint only implements an implication not an equivalence, i.e., it does not ensure
 * that \f$z = 1\f$ if \f$ax \leq b\f$ or equivalently if \f$s = 0\f$.
 *
 * This constraint is quivalent to a linear constraint \f$ax - s \leq b\f$ and an SOS1 constraint on
 * \f$z\f$ and \f$s\f$ (at most one should be nonzero). In this context we can, however, separate
 * more inequalities.
 *
 * The name indicator apparently comes from ILOG CPLEX.
 *
 * @section SEPARATION Separation Methods
 *
 * We now explain the handling of indicator constraints in more detail.  The indicator constraint
 * handler adds an inequality for each indicator constraint. We assume that this system (with added
 * slack variables) is \f$ Ax - s \leq b \f$, where \f$ x \f$ are the original variables and \f$ s
 * \f$ are the slack variables added by the indicator constraint. Variables \f$ y \f$ are the binary
 * variables corresponding to the indicator constraints.
 *
 * @note In the implementation, we assume that bounds on the original variables \f$x\f$ cannot be
 * influenced by the indicator constraint. If it should be possible to relax these constraints as
 * well, then these constraints have to be added as indicator constraints.
 *
 * We separate inequalities by two methods:
 * - The first uses the so-called alternative polyhedron.
 * - The second uses SCIP's conflict analysis.
 *
 * @section ALTERNATIVEPOLYHEDRON Separation via the Alternative Polyhedron
 *
 * We now describe the separation method of the first method in more detail.
 *
 * Consider the LP-relaxation of the current subproblem:
 * \f[
 * \begin{array}{ll}
 *  min & c^T x + d^T y\\
 *      & A x - s \leq b \\
 *      & D x + C z \leq f \\
 *      & l \leq x \leq u \\
 *      & u \leq z \leq v \\
 *      & 0 \leq s
 * \end{array}
 * \f]
 * As above \f$Ax - s \leq b\f$ contains all inequalities corresponding to indicator constraints,
 * while the system \f$Dx + Cy \leq f\f$ contains all other inequalities (which are ignored in the
 * following). Additional bounds for each variable can be given, in particular, variables can be
 * fixed.
 *
 * To generate cuts, we construct the so-called @a alternative @a polyhedron:
 * \f[
 * \begin{array}{ll}
 *      P = \{ y : & w^T A - r + t = 0,\\
 *                 & w^T b - l^T r + u^T t = -1,\\
 *                 & w, r, t \geq 0 \}.
 * \end{array}
 * \f]
 * Here, \f$r\f$ and \f$t\f$ correspond to the lower and upper bounds on \f$x\f$, respectively.
 *
 * It turns out that the vertices of \f$P\f$ correspond to minimal infeasible subsystems of \f$A x
 * \leq b\f$. If \f$I\f$ is the index set of such a system, it follows that not all \f$s_i\f$ for
 * \f$i \in I\f$ can be 1. In other words, the following cut is valid:
 * \f[
 *      \sum_{i \in I} y_i \leq |I| - 1.
 * \f]
 *
 * @subsection DETAIL Separation heuristic
 *
 * We separate the above inequalities by a heuristic described in
 *
 * Branch-And-Cut for the Maximum Feasible Subsystem Problem,@n
 * Marc Pfetsch, SIAM Journal on Optimization 19, No.1, 21-38 (2008)
 *
 * The first step in the separation heuristic is to apply the transformation \f$y' = 1 - y\f$, which
 * transforms the above inequality into the set covering constraint
 * \f[
 *      \sum_{i \in I} y_i' \geq 1,
 * \f]
 * that is, it is a set covering constraint on the negated variables.
 *
 * The basic idea is to use the current solution to the LP relaxation and use it as the objective,
 * when optimizing of the alternative polyhedron. Since any vertex corresponds to such an
 * inequality, we can check whether it is violated. To enlarge the chance that we find a @em
 * violated inequality, we perform a fixing procedure, in which the variable corresponding to an
 * arbitrary element of the last IIS \f$I\f$ is fixed to zero, i.e., cannot be used in the next
 * IISs. This is repeated until the corresponding alternative polyhedron is infeasible, i.e., we
 * have obtained an IIS-cover. For more details see the paper above.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
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


/* maximal value for coupling inequalities */
#define maxCouplingValue  1e3


/** constraint data for indicator constraints */
struct SCIP_ConsData
{
   SCIP_VAR*   binvar;              /**< binary variable for indicator constraint */
   SCIP_VAR*   slackvar;            /**< slack variable of inequality of indicator constraint */
   SCIP_CONS*  lincons;             /**< linear constraint corresponding to indicator constraint */
   int         nFixedNonzero;       /**< number of variables among binvar and slackvar fixed to be nonzero */
   int         colIndex;            /**< column index in alternative LP */
};

/** indicator constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR* eventhdlr;       /**< event handler for bound change events */
   SCIP_Bool   removable;           /**< whether the separated cuts should be removable */
   SCIP_LPI*   altLP;               /**< alternative LP for cut separation */
   int         nvars;               /**< total number of linear constraint variables = number of rows in alt LP - 1 */
   int         nLbBounds;           /**< number of lower bounds of original variables */
   int         nUbBounds;           /**< number of upper bounds of original variables */
   SCIP_HASHMAP* varHash;           /**< hash map from variable to row index in alternative LP */
   SCIP_HASHMAP* lbHash;            /**< hash map from variable to index of lower bound column in alternative LP */
   SCIP_HASHMAP* ubHash;            /**< hash map from variable to index of upper bound column in alternative LP */
   int roundingRounds;              /**< number of rounds in separation */
   SCIP_Real roundingMinThreshold;  /**< minimal value for rounding in separation */
   SCIP_Real roundingMaxThreshold;  /**< maximal value for rounding in separation */
   SCIP_Real roundingOffset;        /**< offset for rounding in separation */
   SCIP_Bool branchIndicators;      /**< Branch on indicator constraints in enforcing? */
   SCIP_Bool genLogicor;            /**< Generate logicor constraints instead of cuts? */
   SCIP_Bool sepaAlternativeLP;     /**< Separate using the alternative LP? */
   SCIP_Bool addCoupling;           /**< whether the coupling inequalities should be added */
};


/* ------------------------ operations on the alternative LP -------------------*/

/** initialize alternative LP */
static
SCIP_RETCODE initAlternativeLP(
      SCIP* scip,                 /**< SCIP pointer */
      SCIP_CONSHDLR* conshdlr     /**< constraint handler */
      )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RETCODE retcode;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert( scip != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   assert( conshdlrdata->altLP == NULL );
   assert( conshdlrdata->varHash == NULL );
   assert( conshdlrdata->lbHash == NULL );
   assert( conshdlrdata->ubHash == NULL );

   SCIPdebugMessage("Initializing alternative LP ...\n");

   /* create hash map of variables */
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->varHash, SCIPblkmem(scip), 10 * SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->lbHash, SCIPblkmem(scip), 10 * SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->ubHash, SCIPblkmem(scip), 10 * SCIPgetNVars(scip)) );

   /* create alternative LP */
   SCIP_CALL( SCIPlpiCreate(&conshdlrdata->altLP, "altLP", SCIP_OBJSEN_MINIMIZE) );

   /* add first row */
   lhs = -1.0;
   rhs = -1.0;
   SCIP_CALL( SCIPlpiAddRows(conshdlrdata->altLP, 1, &lhs, &rhs, NULL, 0, NULL, NULL, NULL) );

   /* set parameters */
   retcode = SCIPlpiSetIntpar(conshdlrdata->altLP, SCIP_LPPAR_FROMSCRATCH, FALSE);
   if ( retcode != SCIP_OKAY && retcode != SCIP_PARAMETERUNKNOWN )
   {
      SCIPerrorMessage("Error <%d> in function call\n", retcode);
      return retcode;
   }

   retcode = SCIPlpiSetIntpar(conshdlrdata->altLP, SCIP_LPPAR_PRESOLVING, TRUE);
   if ( retcode != SCIP_OKAY && retcode != SCIP_PARAMETERUNKNOWN )
   {
      SCIPerrorMessage("Error <%d> in function call\n", retcode);
      return retcode;
   }

   retcode = SCIPlpiSetIntpar(conshdlrdata->altLP, SCIP_LPPAR_SCALING, TRUE);
   if ( retcode != SCIP_OKAY && retcode != SCIP_PARAMETERUNKNOWN )
   {
      SCIPerrorMessage("Error <%d> in function call\n", retcode);
      return retcode;
   }

   /* set constraint handler data */
   SCIPconshdlrSetData(conshdlr, conshdlrdata);

   /* uncomment the following for debugging */
   /* SCIP_CALL( SCIPlpiSetIntpar(conshdlrdata->altLP, SCIP_LPPAR_LPINFO, TRUE) ); */

   return SCIP_OKAY;
}


/** Check whether the bounds are set correctly (for debugging) */
#ifndef NDEBUG
static
SCIP_RETCODE checkLPBoundsClean(
        SCIP* scip,         /**< SCIP pointer */
	SCIP_LPI* lp,       /**< lp for which bounds should be checked */
	int nconss,         /**< number of constraints */
	SCIP_CONS** conss   /**< constraints */
	)
{
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Bool* covered;
   int nCols;
   int j;

   assert( scip != NULL );
   assert( lp != NULL );

   SCIP_CALL( SCIPlpiGetNCols(lp, &nCols) );

   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nCols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nCols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &covered, nCols) );

   for (j = 0; j < nCols; ++j)
      covered[j] = FALSE;

   /* check columns used by contraints */
   SCIP_CALL( SCIPlpiGetBounds(lp, 0, nCols-1, lb, ub) );
   for (j = 0; j < nconss; ++j)
   {
      SCIP_CONSDATA* consdata;
      int ind;

      assert( conss[j] != NULL );
      consdata = SCIPconsGetData(conss[j]);
      assert( consdata != NULL );
      ind = consdata->colIndex;
      assert( 0 <= ind && ind < nCols );
      covered[ind] = TRUE;
      if ( ! SCIPisZero(scip, lb[ind]) || ! SCIPlpiIsInfinity(lp, ub[ind]) )
	 SCIPABORT();
   }

   /* check other columns */
   for (j = 0; j < nCols; ++j)
   {
      if (! covered[j] )
      {
	 /* some columns can be fixed to 0, since they correspond to disabled constraints */
	 if ( ! SCIPisZero(scip, lb[j]) || (! SCIPlpiIsInfinity(lp, ub[j]) && ! SCIPisZero(scip, ub[j])) )
	    SCIPABORT();
      }
   }

   SCIPfreeBufferArray(scip, &covered);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &ub);

   return SCIP_OKAY;
}
#endif


/** Set the alternative system objective function
 *
 *  We assume that the objective function coefficients of the
 *  variables other than the binary indicators are always 0 and
 *  hence do not have to be changed.
 *
 *  We already use the tranformation \f$y' = 1 - y\f$.
 */
static
SCIP_RETCODE setAltLPObj(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_LPI* lp,             /**< alternative LP */
      SCIP_SOL* sol,            /**< solution to be dealt with */
      int nconss,               /**< number of constraints */
      SCIP_CONS** conss         /**< indicator constraints */
      )
{
   int j;
   SCIP_Real* obj;
   int* indices;

   assert( scip != NULL );
   assert( lp != NULL );
   assert( conss != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indices, nconss) );

   for (j = 0; j < nconss; ++j)
   {
      SCIP_CONSDATA* consdata;

      assert( conss[j] != NULL );
      consdata = SCIPconsGetData(conss[j]);
      assert( consdata != NULL );

      obj[j] = 1.0 - SCIPgetSolVal(scip, sol, consdata->binvar);
      assert( consdata->colIndex >= 0 );
      indices[j] = consdata->colIndex;
   }

   SCIP_CALL( SCIPlpiChgObj(lp, nconss, indices, obj) );

   SCIPfreeBufferArray(scip, &indices);
   SCIPfreeBufferArray(scip, &obj);

   return SCIP_OKAY;
}


/** Fix variable given by @a S to 0 */
static
SCIP_RETCODE fixAltLPVariables(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_LPI* lp,             /**< alternative LP */
      int nconss,               /**< number of constraints */
      SCIP_CONS** conss,        /**< indicator constraints */
      SCIP_Bool* S              /**< bitset of variables */
      )
{
   int j;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int* indices;
   int cnt = 0;

   assert( scip != NULL );
   assert( lp != NULL );
   assert( conss != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indices, nconss) );

   /* collect bounds to be changed */
   for (j = 0; j < nconss; ++j)
   {
      SCIP_CONSDATA* consdata;

      assert( conss[j] != NULL );
      consdata = SCIPconsGetData(conss[j]);
      assert( consdata != NULL );
      assert( consdata->colIndex >= 0 );

      if ( S[j] )
      {
	 indices[cnt] = consdata->colIndex;
	 lb[cnt] = 0.0;
	 ub[cnt] = 0.0;
	 ++cnt;
      }
   }
   /* change bounds */
   SCIP_CALL( SCIPlpiChgBounds(lp, cnt, indices, lb, ub) );

   SCIPfreeBufferArray(scip, &indices);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);

   return SCIP_OKAY;
}


/** Fix variable @p ind to 0 */
static
SCIP_RETCODE fixAltLPVariable(
      SCIP_LPI* lp,             /**< alternative LP */
      int ind                   /**< variable that should be fixed to 0 */
      )
{
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 0.0;

   /* change bounds */
   SCIP_CALL( SCIPlpiChgBounds(lp, 1, &ind, &lb, &ub) );

   return SCIP_OKAY;
}


/** unfix variable @p ind to 0 */
static
SCIP_RETCODE unfixAltLPVariable(
      SCIP_LPI* lp,             /**< alternative LP */
      int ind                   /**< variable that should be fixed to 0 */
      )
{
   SCIP_Real lb = 0.0;
   SCIP_Real ub = SCIPlpiInfinity(lp);

   /* change bounds */
   SCIP_CALL( SCIPlpiChgBounds(lp, 1, &ind, &lb, &ub) );

   return SCIP_OKAY;
}


/** unfix variable given by @a S to 0 */
static
SCIP_RETCODE unfixAltLPVariables(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_LPI* lp,             /**< alternative LP */
      int nconss,               /**< number of constraints */
      SCIP_CONS** conss,        /**< indicator constraints */
      SCIP_Bool* S              /**< bitset of variables */
      )
{
   int j;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int* indices;
   int cnt = 0;

   assert( scip != NULL );
   assert( lp != NULL );
   assert( conss != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indices, nconss) );

   /* collect bounds to be changed */
   for (j = 0; j < nconss; ++j)
   {
      if ( S[j] )
      {
	 SCIP_CONSDATA* consdata;

	 assert( conss[j] != NULL );
	 consdata = SCIPconsGetData(conss[j]);
	 assert( consdata != NULL );
	 assert( consdata->colIndex >= 0 );

	 indices[cnt] = consdata->colIndex;
	 lb[cnt] = 0.0;
	 ub[cnt] = SCIPlpiInfinity(lp);
	 ++cnt;
      }
   }
   /* change bounds */
   SCIP_CALL( SCIPlpiChgBounds(lp, cnt, indices, lb, ub) );

   SCIPfreeBufferArray(scip, &indices);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);

   return SCIP_OKAY;
}



/** add column corresponding to constraint to alternative LP */
static
SCIP_RETCODE addAltLPConstraint(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_CONSHDLR* conshdlr,  /**< constraint handler */
      SCIP_CONS* cons           /**< new indicator constraint */
      )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* lincons;
   SCIP_VAR** linvars;
   SCIP_Real* linvals;
   SCIP_VAR* slackvar;
   int nlinvars;
   SCIP_Real val;
   SCIP_Real sign = 1.0;
   int v;
   int* matbeg;
   int* matind;
   SCIP_Real* matval;
   int nNewRows = 0;
   SCIP_VAR** newVars;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int nCols;
   int cnt = 0;
   int nNewCols = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMessage("Adding column to alternative LP ...\n");

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->altLP == NULL )
      SCIP_CALL( initAlternativeLP(scip, conshdlr) );
   assert( conshdlrdata->varHash != NULL );
   assert( conshdlrdata->lbHash != NULL );
   assert( conshdlrdata->ubHash != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   lincons = consdata->lincons;
   assert( lincons != NULL );

   slackvar = consdata->slackvar;
   linvars = SCIPgetVarsLinear(scip, lincons);
   linvals = SCIPgetValsLinear(scip, lincons);
   nlinvars = SCIPgetNVarsLinear(scip, lincons);

#ifndef NDEBUG
   {
      int nRows;
      SCIP_CALL( SCIPlpiGetNRows(conshdlrdata->altLP, &nRows) );
      assert( nRows == conshdlrdata->nvars+1 );
   }
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &matbeg, nlinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matind, 2*nlinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matval, 2*nlinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newVars, nlinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, 2*nlinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, 2*nlinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, 2*nlinvars) );

   /* store index of column in constraint */
   SCIP_CALL( SCIPlpiGetNCols(conshdlrdata->altLP, &nCols) );
   assert( consdata->colIndex < 0 );
   consdata->colIndex = nCols;

   /* adapt rhs of linear constraint */
   val = SCIPgetRhsLinear(scip, lincons);
   if ( SCIPisInfinity(scip, val) )
   {
      val = SCIPgetLhsLinear(scip, lincons);
      assert( val > -SCIPinfinity(scip) );
      sign = -1.0;
   }

   /* handle first row */
   if (! SCIPisFeasZero(scip, val) )
   {
      matind[cnt] = 0;
      matval[cnt] = sign * val;
      ++cnt;
   }

   /* set up column (recognize new original variables) */
   for (v = 0; v < nlinvars; ++v)
   {
      SCIP_VAR* var;
      var = linvars[v];
      assert( var != NULL );

      if ( var != slackvar )
      {
	 /* if variable is new */
	 if ( ! SCIPhashmapExists(conshdlrdata->varHash, var) )
	 {
	    /* add variable in map and array and remember to add a new row */
	    SCIP_CALL( SCIPhashmapInsert(conshdlrdata->varHash, var, (void*) (size_t) conshdlrdata->nvars) );
	    assert( conshdlrdata->nvars == (int) (size_t) SCIPhashmapGetImage(conshdlrdata->varHash, var) );
	    SCIPdebugMessage("inserted variable <%s> into hashmap (%d)\n", SCIPvarGetName(var), conshdlrdata->nvars);

	    matind[cnt] = conshdlrdata->nvars + 1;

	    /* store new variables */
	    newVars[nNewRows] = var;
	    ++(conshdlrdata->nvars);
	    ++nNewRows;
	 }
	 else
	    matind[cnt] = (int) (size_t) SCIPhashmapGetImage(conshdlrdata->varHash, var) + 1;

	 assert( SCIPhashmapExists(conshdlrdata->varHash, var) );
	 matval[cnt] = sign * linvals[v];
	 ++cnt;
      }
   }

   /* if we added new rows */
   if ( nNewRows > 0 )
   {
      SCIP_Real* lhs;
      SCIP_Real* rhs;
      int i;

      SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nNewRows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nNewRows) );
      for (i = 0; i < nNewRows; ++i)
      {
	 lhs[i] = 0.0;
	 rhs[i] = 0.0;
      }
      /* add new rows */
      SCIP_CALL( SCIPlpiAddRows(conshdlrdata->altLP, nNewRows, lhs, rhs, NULL, 0, NULL, NULL, NULL) );

      SCIPfreeBufferArray(scip, &lhs);
      SCIPfreeBufferArray(scip, &rhs);
   }

   /* now add column */
   obj[0] = 1.0;
   lb[0] = 0.0;
   ub[0] = SCIPlpiInfinity(conshdlrdata->altLP);
   matbeg[0] = 0;

   SCIP_CALL( SCIPlpiAddCols(conshdlrdata->altLP, 1, obj, lb, ub, NULL, cnt, matbeg, matind, matval) );

   /* add columns corresponding to bounds of original variables */
   cnt = 0;
   nNewCols = 0;
   for (v = 0; v < nNewRows; ++v)
   {
      SCIP_VAR* var = newVars[v];

      /* if the lower bound is finite */
      val  = SCIPvarGetLbLocal(var);
      if ( ! SCIPisInfinity(scip, -val) )
      {
	 matbeg[nNewCols] = cnt;
	 if ( ! SCIPisZero(scip, val) )
	 {
	    matind[cnt] = 0;
	    matval[cnt] = -val;
	    ++cnt;
	 }
	 assert( SCIPhashmapExists(conshdlrdata->varHash, var) );
	 matind[cnt] = (int) (size_t) SCIPhashmapGetImage(conshdlrdata->varHash, var) + 1;
	 matval[cnt] = -1.0;
	 ++cnt;
	 obj[nNewCols] = 0.0;
	 lb[nNewCols] = 0.0;
	 ub[nNewCols] = SCIPlpiInfinity(conshdlrdata->altLP);
	 ++conshdlrdata->nLbBounds;
	 SCIP_CALL( SCIPhashmapInsert(conshdlrdata->lbHash, var, (void*) (size_t) (nCols + nNewCols)) );
	 assert( SCIPhashmapExists(conshdlrdata->lbHash, var) );
	 SCIPdebugMessage("added column corr. to lower bound (%f) of variable <%s> to alternative polyhedron\n", val, SCIPvarGetName(var));
	 ++nNewCols;
      }

      /* if the upper bound is finite */
      val = SCIPvarGetUbLocal(var);
      if ( ! SCIPisInfinity(scip, val) )
      {
	 matbeg[nNewCols] = cnt;
	 if ( ! SCIPisZero(scip, val) )
	 {
	    matind[cnt] = 0;
	    matval[cnt] = val;
	    ++cnt;
	 }
	 assert( SCIPhashmapExists(conshdlrdata->varHash, var) );
	 matind[cnt] = (int) (size_t) SCIPhashmapGetImage(conshdlrdata->varHash, var) + 1;
	 matval[cnt] = 1.0;
	 ++cnt;
	 obj[nNewCols] = 0.0;
	 lb[nNewCols] = 0.0;
	 ub[nNewCols] = SCIPlpiInfinity(conshdlrdata->altLP);
	 ++conshdlrdata->nUbBounds;
	 SCIP_CALL( SCIPhashmapInsert(conshdlrdata->ubHash, var, (void*) (size_t) (nCols + nNewCols)) );
	 assert( SCIPhashmapExists(conshdlrdata->ubHash, var) );
	 SCIPdebugMessage("added column corr. to upper bound (%f) of variable <%s> to alternative polyhedron\n", val, SCIPvarGetName(var));
	 ++nNewCols;
      }
   }

   /* add columns if necessary */
   if ( nNewCols > 0 )
   {
      SCIP_CALL( SCIPlpiAddCols(conshdlrdata->altLP, nNewCols, obj, lb, ub, NULL, cnt, matbeg, matind, matval) );
   }

#ifndef NDEBUG
   SCIP_CALL( SCIPlpiGetNCols(conshdlrdata->altLP, &cnt) );
   assert( cnt == nCols + nNewCols + 1 );
#endif

   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &newVars);
   SCIPfreeBufferArray(scip, &matind);
   SCIPfreeBufferArray(scip, &matval);
   SCIPfreeBufferArray(scip, &matbeg);

   SCIP_CALL( SCIPlpiWriteLP(conshdlrdata->altLP, "alt.lp") );

   return SCIP_OKAY;
}




/** delete column corresponding to constraint in alternative LP
 *
 *  We currently just fix the corresponding variable to 0.
 */
static
SCIP_RETCODE deleteAltLPConstraint(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_CONSHDLR* conshdlr,  /**< constraint handler */
      SCIP_CONS* cons           /**< indicator constraint */
      )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->altLP != NULL )
   {
      SCIP_CONSDATA* consdata;

      SCIPdebugMessage("Deleting column from alternative LP ...\n");

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      SCIP_CALL( fixAltLPVariable(conshdlrdata->altLP, consdata->colIndex) );
      consdata->colIndex = -1;
   }

   return SCIP_OKAY;
}









/** Check whether the given LP is infeasible
 *
 * If @a primal is false we assume that the problem is <em>dual feasible</em>, e.g.,
 * the problem was only changed by fixing bounds!
 *
 * This is the workhorse for all methods that have to solve the alternative LP.
 * We try in several ways to recover from possible stability problems.
 *
 * @pre It is assumed that all parameters for the alternative LP are set.
 */
static
SCIP_RETCODE checkAltLPInfeasible(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_LPI* lp,             /**< LP */
      SCIP_Bool primal,         /**< whether we are using the primal or dual simplex */
      SCIP_Bool* infeasible,    /**< output: whether the LP is infeasible */
      SCIP_Bool* error          /**< output: whether an error occured */
      )
{
   assert( scip != NULL );
   assert( lp != NULL );
   assert( infeasible != NULL );
   assert( error != NULL );

   *error = FALSE;

   /* solve LP */
   if ( primal )
      SCIP_CALL( SCIPlpiSolvePrimal(lp) );  /* use primal simplex */
   else
      SCIP_CALL( SCIPlpiSolveDual(lp) );    /* use dual simplex */

   /* resolve if LP is not stable */
   if ( ! SCIPlpiIsStable(lp) )
   {
      SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, TRUE) );
      SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_PRESOLVING, FALSE) );
      SCIPwarningMessage("Numerical problems, retrying ...\n");

      /* re-solve LP */
      if ( primal )
	 SCIP_CALL( SCIPlpiSolvePrimal(lp) );  /* use primal simplex */
      else
	 SCIP_CALL( SCIPlpiSolveDual(lp) );    /* use dual simplex */

      /* reset parameters */
      SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, FALSE) );
      SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_PRESOLVING, TRUE) );
   }

   /* check whether we are in the paradoxical situtation that
    * - the primal is not infeasible
    * - the primal is not unbounded
    * - the LP is not optimal
    * - we have a primal ray
    *
    * If we ran the dual simplex algorithm, then we run again with the primal simplex
    */
   if ( ! SCIPlpiIsPrimalInfeasible(lp) && ! SCIPlpiIsPrimalUnbounded(lp) &&
	! SCIPlpiIsOptimal(lp) && SCIPlpiExistsPrimalRay(lp) && !primal )
   {
      SCIPwarningMessage("The dual simplex produced a primal ray. Retrying with primal ...\n");
      /* the following settings might be changed: */
      SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, TRUE) );
      SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_PRESOLVING, TRUE) );
      SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_SCALING, TRUE) );

      SCIP_CALL( SCIPlpiSolvePrimal(lp) );   /* use primal simplex */

      /* reset parameters */
      SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, FALSE) );
      SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_PRESOLVING, TRUE) );
      SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_SCALING, TRUE) );
   }

   /* examine LP solution status */
   if ( SCIPlpiIsPrimalInfeasible(lp) )     /* the LP is provably infeasible */
   {
      assert( ! SCIPlpiIsPrimalUnbounded(lp) );   /* can't be unbounded or optimal */
      assert( ! SCIPlpiIsOptimal(lp) );           /* if it is infeasible! */
      *infeasible = TRUE;                         /* LP is infeasible */
      return SCIP_OKAY;
   }
   else
   {
      /* By assumption the dual is feasible if the dual simplex is run, therefore
       * the status has to be primal unbounded or optimal. */
      if ( ! SCIPlpiIsPrimalUnbounded(lp) && ! SCIPlpiIsOptimal(lp) )
      {
	 /* We have a status different from unbounded or optimal. This should not be the case ... */
	 if (primal)
	 {
	    SCIPerrorMessage("Primal simplex returned with unknown status: %d\n", SCIPlpiGetInternalStatus(lp));
	 }
	 else
	 {
	    SCIPerrorMessage("Dual simplex returned with unknown status: %d\n", SCIPlpiGetInternalStatus(lp));
	 }
	 /* SCIP_CALL( SCIPlpiWriteLP(lp, "debug.lp") ); */
	 *error = TRUE;
	 return SCIP_OKAY;
      }
   }

   /* at this point we have a feasible solution */
   *infeasible = FALSE;
   return SCIP_OKAY;
}



/** Tries to extend a given set of variables to a cover.
 *
 * At each step we include a variable which covers a new IIS. Ties are
 * broken according to the number of IISs a variable is contained in.
 * The corresponding IIS inequalities are added to the LP if this not
 * already happend.
 *
 * @pre It is assumed that all parameters for the alternative LP are
 * set and that the variables corresponding to @a S are
 * fixed. Furthermore @c xVal_ should contain the current LP solution.
 */
static
SCIP_RETCODE extendToCover(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_LPI* lp,             /**< LP */
      SCIP_SOL* sol,            /**< solution to be separated */
      SCIP_Bool removable,      /**< whether cuts should be removable */
      SCIP_Bool genLogicor,     /**< should logicor constraints be generated? */
      int nconss,               /**< number of constraints */
      SCIP_CONS** conss,        /**< indicator constraints */
      SCIP_Bool* S,             /**< bitset of variables */
      int* size,                /**< size of S */
      SCIP_Real* value,         /**< objective value of S */
      int* nGen                 /**< number of generated cuts */
      )
{
   int step = 0;
   SCIP_Real* primsol = NULL;
   int nCols = 0;

   assert( scip != NULL );
   assert( lp != NULL );
   assert( conss != NULL );
   assert( S != NULL );
   assert( size != NULL );
   assert( value != NULL );
   assert( nGen != NULL );

   SCIP_CALL( SCIPlpiGetNCols(lp, &nCols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &primsol, nCols) );
   assert( nconss <= nCols );

   do
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_Bool error = FALSE;
      SCIP_Real sum = 0.0;
      SCIP_Real sizeIIS = 0;
      int candidate = -1;
      int candIndex = -1;
      SCIP_Real candObj = -1.0;
      int j;

      if ( step == 0 )
      {
	 /* the first LP is solved without warm start, after that we use a warmstart. */
	 SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, TRUE) );
	 SCIP_CALL( checkAltLPInfeasible(scip, lp, TRUE, &infeasible, &error) );
	 SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_FROMSCRATCH, FALSE) );
      }
      else
	 SCIP_CALL( checkAltLPInfeasible(scip, lp, FALSE, &infeasible, &error) );

      if ( error )
	 break;

      if ( infeasible )
	 break;

      /* get solution of alternative LP */
      SCIP_CALL( SCIPlpiGetSol(lp, NULL, primsol, NULL, NULL, NULL) );

      /* get value of cut and find candidate for variable to add */
      for (j = 0; j < nconss; ++j)
      {
	 SCIP_CONSDATA* consdata;
	 int ind;

	 consdata = SCIPconsGetData(conss[j]);
	 assert( consdata != NULL );
	 ind = consdata->colIndex;
	 assert( 0 <= ind && ind < nCols );

	 /* check support of the solution, i.e., the corresponding IIS */
	 if ( ! SCIPisFeasZero(scip, primsol[ind]) )
	 {
	    assert( ! S[j] );
	    ++sizeIIS;
	    sum += SCIPgetSolVal(scip, sol, consdata->binvar);
	    /* take first element */
	    if ( candidate < 0 )
	    {
	       candidate = j;
	       candIndex = ind;
	       candObj = SCIPvarGetObj(consdata->binvar);
	    }
	 }
      }
      assert( candidate >= 0 );
      assert( ! S[candidate] );

      /* update new set S */
      SCIPdebugMessage("   size: %4d  add %4d with value %f\n", *size, candidate, candObj);
      S[candidate] = TRUE;
      ++(*size);
      *value += candObj;

      /* fix chosen variable to 0 */
      SCIP_CALL( fixAltLPVariable(lp, candIndex) );

      /* if cut is violated, i.e., sum - sizeIIS + 1 > 0 */
      if ( SCIPisEfficacious(scip, sum - (SCIP_Real) sizeIIS + 1.0) )
      {
	 if ( genLogicor )
	 {
	    SCIP_CONS* cons;
	    SCIP_VAR** vars;
	    int cnt =0;

	    SCIP_CALL( SCIPallocBufferArray(scip, &vars, nconss) );

	    /* collect variables corresponding to support to cut */
	    for (j = 0; j < nconss; ++j)
	    {
	       int ind;
	       SCIP_CONSDATA* consdata;
	       consdata = SCIPconsGetData(conss[j]);
	       ind = consdata->colIndex;
	       assert( 0 <= ind && ind < nCols );
	       assert( consdata->binvar != NULL );

	       /* check support of the solution, i.e., the corresponding IIS */
	       if ( ! SCIPisFeasZero(scip, primsol[ind]) )
	       {
		  SCIP_VAR* var;
		  SCIP_CALL( SCIPgetNegatedVar(scip, consdata->binvar, &var) );
		  vars[cnt++] = var;
	       }
	    }
	    assert( cnt == sizeIIS );

	    SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, "iis", cnt, vars, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, removable, FALSE) );

#ifdef SCIP_OUTPUT
	    SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
#endif

	    SCIP_CALL( SCIPaddCons(scip, cons) );
	    SCIP_CALL( SCIPreleaseCons(scip, &cons) );

	    SCIPfreeBufferArray(scip, &vars);
	    ++(*nGen);
	 }
	 else
	 {
	    SCIP_ROW* row;

	    /* create row */
	    SCIP_CALL( SCIPcreateEmptyRow(scip, &row, "iis", -SCIPinfinity(scip), (SCIP_Real) sizeIIS - 1.0, FALSE, FALSE, removable) );
	    SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

	    /* add variables corresponding to support to cut */
	    for (j = 0; j < nconss; ++j)
	    {
	       int ind;
	       SCIP_CONSDATA* consdata;
	       consdata = SCIPconsGetData(conss[j]);
	       ind = consdata->colIndex;
	       assert( 0 <= ind && ind < nCols );
	       assert( consdata->binvar != NULL );

	       /* check support of the solution, i.e., the corresponding IIS */
	       if ( ! SCIPisFeasZero(scip, primsol[ind]) )
		  SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->binvar, 1.0) );
	    }
	    SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_OUTPUT
	    SCIProwPrint(row, NULL);
#endif
	    SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );

            /* cut should be violated: */
            assert( SCIPisFeasNegative(scip, SCIPgetRowLPFeasibility(scip, row)) );

	    /* add cuts to pool if we are in the root -> cuts are globally valid */
	    if ( SCIPgetDepth(scip) == 0 )
	       SCIP_CALL( SCIPaddPoolCut(scip, row) );
	    SCIP_CALL( SCIPreleaseRow(scip, &row));
	    ++(*nGen);
	 }
      }
      ++step;
   }
   while (step < nconss);

   SCIPfreeBufferArray(scip, &primsol);

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

	 /* if variable is not multi-aggregated */
	 if ( SCIPvarGetStatus(consdata->slackvar) != SCIP_VARSTATUS_MULTAGGR )
	 {
	    SCIP_CALL( SCIPinferVarUbCons(scip, consdata->slackvar, 0.0, cons, 0, &infeasible, &tightened) );
	    assert( ! infeasible );
	    if ( tightened )
	       ++(*nGen);
	 }
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
SCIP_RETCODE enforceIndicators(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_CONSHDLR* conshdlr,  /**< constraint handler */
      int nconss,               /**< number of constraints */
      SCIP_CONS** conss,        /**< indicator constraints */
      SCIP_RESULT* result       /**< result */
      )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_NODE* node1;
   SCIP_NODE* node2;
   SCIP_VAR* slackvar;
   SCIP_VAR* binvar;
   SCIP_CONS* branchCons = NULL;
   SCIP_Real maxSlack = -1.0;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   SCIPdebugMessage("Enforcing indicator constraints <%s>.\n", SCIPconshdlrGetName(conshdlr) );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool cutoff;
      SCIP_Real valSlack = 0.0;
      int cnt = 0;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* first perform propagation (it might happen that standard propagation is turned off) */
      SCIP_CALL( propIndicator(scip, conss[c], consdata, &cutoff, &cnt) );
      SCIPdebugMessage("propagation in enforcing <%s> (cutoff: %d, domain reductions: %d).\n", SCIPconsGetName(conss[c]), cutoff, cnt);
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

      /* check whether constraint is infeasible */
      binvar = consdata->binvar;
      valSlack = SCIPgetSolVal(scip, NULL, consdata->slackvar);
      assert( ! SCIPisFeasNegative(scip, valSlack) );
      if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, binvar)) &&
	   ! SCIPisFeasZero(scip, valSlack) )
      {
	 /* binary variable is not fixed - otherwise we would not be infeasible */
	 assert( SCIPvarGetLbLocal(binvar) < 0.5 && SCIPvarGetUbLocal(binvar) > 0.5 );

	 if ( valSlack > maxSlack )
	 {
	    maxSlack = valSlack;
	    branchCons = conss[c];
	 }
      }
   }

   /* if all constraints are feasible */
   if ( branchCons == NULL )
   {
      SCIPdebugMessage("All indicator constraints are feasible.\n");
      return SCIP_OKAY;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* skip branching if required */
   if ( ! conshdlrdata->branchIndicators )
   {
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   /* otherwise create branches */
   SCIPdebugMessage("Branching on constraint <%s> (slack value: %f).\n", SCIPconsGetName(branchCons), maxSlack);
   consdata = SCIPconsGetData(branchCons);
   assert( consdata != NULL );
   binvar = consdata->binvar;
   slackvar = consdata->slackvar;

   /* node1: binvar = 1, slackvar = 0 */
   SCIP_CALL( SCIPcreateChild(scip, &node1, 0.0, SCIPcalcChildEstimate(scip, binvar, 1.0) ) );

   if ( SCIPvarGetLbLocal(binvar) < 0.5 )
      SCIP_CALL( SCIPchgVarLbNode(scip, node1, binvar, 1.0) );

   /* if slack-variable is multi-aggregated */
   if ( SCIPvarGetStatus(slackvar) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_CONS* cons;
      SCIP_Real val = 1.0;

      if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(slackvar)) || ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(slackvar)) )
      {
	 SCIPdebugMessage("creating constraint to force multi-aggregated variable <%s> to 0.\n", SCIPvarGetName(slackvar));
	 /* we have to insert a local constraint var = 0 */
	 SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &slackvar, &val, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
	       TRUE, FALSE, FALSE, FALSE, FALSE) );
	 SCIP_CALL( SCIPaddConsNode(scip, node1, cons, NULL) );
	 SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }
   else
   {
      assert( SCIPvarGetStatus(slackvar) != SCIP_VARSTATUS_MULTAGGR );
      if ( ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(slackvar)) )
	 SCIP_CALL( SCIPchgVarUbNode(scip, node1, slackvar, 0.0) );
   }

   /* node2: binvar = 0, no restriction on slackvar */
   SCIP_CALL( SCIPcreateChild(scip, &node2, 0.0, SCIPcalcChildEstimate(scip, binvar, 0.0) ) );

   if ( SCIPvarGetUbLocal(binvar) > 0.5 )
      SCIP_CALL( SCIPchgVarUbNode(scip, node2, binvar, 0.0) );

   SCIP_CALL( SCIPresetConsAge(scip, branchCons) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}





/** separate IIS-cuts via rounding */
static
SCIP_RETCODE separateIISRounding(
      SCIP* scip,               /**< SCIP pointer */
      SCIP_CONSHDLR* conshdlr,  /**< constraint handler */
      SCIP_SOL* sol,            /**< solution to be separated */
      int nconss,               /**< number of constraints */
      SCIP_CONS** conss,        /**< indicator constraints */
      int* nGen                 /**< number of domain changes */
      )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_LPI* lp;
   int rounds = 0;
   SCIP_Real threshold = 0.0;
   SCIP_Bool* S;
   int nGenOld;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( nGen != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   lp = conshdlrdata->altLP;
   assert( lp != NULL );

   nGenOld = *nGen;
   SCIPdebugMessage("Separating IIS-cuts by rounding ...\n");

#ifndef NDEBUG
   SCIP_CALL( checkLPBoundsClean(scip, lp, nconss, conss) );
#endif

   /* ???????????????????????????????????????????????? */
   /* TODO: change coefficients of bounds in alternative LP */

   /* set obj. func. to current solution */
   SCIP_CALL( setAltLPObj(scip, lp, sol, nconss, conss) );

   SCIP_CALL( SCIPallocBufferArray(scip, &S, nconss) );

   /* loop through the possible thresholds */
   for (threshold = conshdlrdata->roundingMaxThreshold; rounds < conshdlrdata->roundingRounds && threshold >= conshdlrdata->roundingMinThreshold;
	threshold -= conshdlrdata->roundingOffset)
   {
      int size = 0;
      SCIP_Real value = 0.0;
      int nCuts = 0;
      int j;

      SCIPdebugMessage("Threshold: %f\n", threshold);

      /* choose variables that have a value < current threshold value */
      for (j = 0; j < nconss; ++j)
      {
	 SCIP_CONSDATA* consdata;

	 assert( conss[j] != NULL );
	 consdata = SCIPconsGetData(conss[j]);
	 assert( consdata != NULL );

	 if ( SCIPisFeasLT(scip, SCIPgetVarSol(scip, consdata->binvar), threshold) )
	 {
	    S[j] = TRUE;
	    value += SCIPvarGetObj(consdata->binvar);
	    ++size;
	 }
	 else
	    S[j] = FALSE;
      }

      if (size == nconss)
      {
	 SCIPdebugMessage("All variables in the set. Continue ...\n");
	 continue;
      }

      /* fix the variables in S */
      SCIP_CALL( fixAltLPVariables(scip, lp, nconss, conss, S) );

      /* extend set S to a cover and generate cuts */
      SCIP_CALL( extendToCover(scip, lp, sol, conshdlrdata->removable, conshdlrdata->genLogicor, nconss, conss, S, &size, &value, &nCuts) );

      if ( nCuts > 0 )
      {
	 *nGen += nCuts;
	 ++rounds;
      }

      SCIPdebugMessage("Produced cover of size %d with value %f\n", size, value);

      /* todo: check whether the cover is a feasible solution */

      /* reset bounds */
      SCIP_CALL( unfixAltLPVariables(scip, lp, nconss, conss, S) );
   }
   SCIPdebugMessage("Generated %d IISs.\n", *nGen - nGenOld);

#ifndef NDEBUG
   SCIP_CALL( checkLPBoundsClean(scip, lp, nconss, conss) );
#endif

   SCIPfreeBufferArray(scip, &S);

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
   assert( conshdlrdata->altLP == NULL );
   assert( conshdlrdata->varHash == NULL );
   assert( conshdlrdata->lbHash == NULL );
   assert( conshdlrdata->ubHash == NULL );
   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}



/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolIndicator)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      assert( SCIPconsIsTransformed(conss[c]) );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      SCIPdebugMessage("Initializing indicator constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* if not happend already, get transformed linear constraint */
      if ( ! SCIPconsIsTransformed(consdata->lincons) )
      {
	 SCIP_CALL( SCIPgetTransformedCons(scip, consdata->lincons, &consdata->lincons) );
	 assert( consdata->lincons != NULL );
      }

      /* add constraint to alternative LP if not already done */
      if ( conshdlrdata->sepaAlternativeLP && consdata->colIndex < 0 )
	 SCIP_CALL( addAltLPConstraint(scip, conshdlr, conss[c]) );
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolIndicator)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->altLP != NULL )
   {
      assert( conshdlrdata->varHash != NULL );
      assert( conshdlrdata->lbHash != NULL );
      assert( conshdlrdata->ubHash != NULL );
      assert( conshdlrdata->altLP != NULL );

#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "\nStatistics for var hash:\n");
      SCIPhashmapPrintStatistics(conshdlrdata->varHash);
      SCIPinfoMessage(scip, NULL, "\nStatistics for lower bound hash:\n");
      SCIPhashmapPrintStatistics(conshdlrdata->lbHash);
      SCIPinfoMessage(scip, NULL, "\nStatistics for upper bound hash:\n");
      SCIPhashmapPrintStatistics(conshdlrdata->ubHash);
#endif

      SCIPhashmapFree(&conshdlrdata->varHash);
      SCIPhashmapFree(&conshdlrdata->lbHash);
      SCIPhashmapFree(&conshdlrdata->ubHash);

      SCIP_CALL( SCIPlpiFree(&conshdlrdata->altLP) );
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

      if ( conshdlrdata->sepaAlternativeLP )
	 SCIP_CALL( deleteAltLPConstraint(scip, conshdlr, cons) );
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
   consdata->colIndex = -1;
   consdata->lincons = sourcedata->lincons;

   /* if binary variable is fixed to be nonzero */
   consdata->nFixedNonzero = 0;
   if ( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
      ++(consdata->nFixedNonzero);

   /* if slack variable is fixed to be nonzero */
   if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(consdata->slackvar)) )
      ++(consdata->nFixedNonzero);

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

   /* add corresponding column to alternative LP if the constraint is new */
   if ( conshdlrdata->sepaAlternativeLP && SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE )
   {
      SCIP_CALL( addAltLPConstraint(scip, conshdlr, *targetcons) );
   }

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

      /* get check for transformed linear constraint */
      if ( ! SCIPconsIsTransformed(consdata->lincons) )
      {
	 SCIP_CALL( SCIPgetTransformedCons(scip, consdata->lincons, &consdata->lincons) );
	 assert( consdata->lincons != NULL );
      }

      SCIPdebugMessage("Presolving indicator constraint <%s>.\n", SCIPconsGetName(cons) );

      *result = SCIP_DIDNOTFIND;

      /* only run if sucess is possible */
      if ( nrounds == 0 || nnewfixedvars > 0 || nnewaggrvars > 0 || *nfixedvars > oldnfixedvars )
      {
	 SCIP_Bool infeasible, fixed;

	 /* if the binary variable if fixed to nonzero */
	 if ( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
	 {
	    SCIPdebugMessage("Presolving: Binary variable fixed to 1.\n");

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

	 /* if the binary variable if fixed to zero */
	 if ( SCIPvarGetUbLocal(consdata->binvar) < 0.5 )
	 {
	    SCIPdebugMessage("Presolving: Binary variable fixed to 0.\n");

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
	    SCIPdebugMessage("Presolving: Slack variable fixed to nonzero.\n");

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

	 /* if the slack variable if fixed to zero */
	 if ( SCIPisFeasZero(scip, SCIPvarGetUbLocal(consdata->slackvar)) )
	 {
	    SCIPdebugMessage("Presolving: Slack variable fixed to zero.\n");

	    /* otherwise fix binary variable to 0/1 */
	    if ( SCIPvarGetObj(consdata->binvar) >= 0 )
	    {
	       SCIPdebugMessage("Fix binary variable to 0 and delete constraint.\n");
	       SCIP_CALL( SCIPfixVar(scip, consdata->binvar, 0.0, &infeasible, &fixed) );
	    }
	    else
	    {
	       SCIPdebugMessage("Fix binary variable to 1 and delete constraint.\n");
	       SCIP_CALL( SCIPfixVar(scip, consdata->binvar, 1.0, &infeasible, &fixed) );
	    }
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




/** LP initialization method of constraint handler 
 *
 *  For an indicator constraint with binary variable \f$y\f$ and slack variable \f$s\f$ the coupling
 *  inequality \f$s \le M (1-y)\f$ (equivalently: \f$s + M y \le M\f$) is inserted, where \f$M\f$ is
 *  an upper bound on the value of \f$s\f$. If \f$M\f$ is too large the inequality is not inserted.
 */
static
SCIP_DECL_CONSINITLP(consInitlpIndicator)
{
   int c;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* add coupling constraint if required */
      if ( conshdlrdata->addCoupling )
      {
	 SCIP_Real ub;
	 SCIP_ROW* row;
	 char name[50];

	 SCIPdebugMessage("Checking for initial rows for indicator constraint <%s>.\n", SCIPconsGetName(conss[c]) );

	 /* get upper bound for slack variable in linear constraint */
	 ub = SCIPvarGetUbGlobal(consdata->slackvar);
	 assert( ! SCIPisNegative(scip, ub) );

	 /* insert corresponding row if helpful and coefficient is not too large */
	 if ( ub < maxCouplingValue )
	 {
#ifndef NDEBUG
	    sprintf(name, "couple%d", c);
#else
	    name[0] = '\0';
#endif
	    SCIP_CALL( SCIPcreateEmptyRow(scip, &row, name, -SCIPinfinity(scip), ub, FALSE, FALSE, FALSE) );
	    SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

	    SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->slackvar, 1.0) );
	    SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->binvar, ub) );
	    SCIP_CALL( SCIPflushRowExtensions(scip, row) );

	    SCIPdebugMessage("Insert coupling inequality for indicator constraint <%s> (coeff: %f).\n", SCIPconsGetName(conss[c]), ub);
#ifdef SCIP_OUTPUT
	    SCIProwPrint(row, NULL);
#endif
	    SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE) );

	    SCIP_CALL( SCIPaddPoolCut(scip, row) );
	    SCIP_CALL( SCIPreleaseRow(scip, &row));
	 }
      }
   }

   return SCIP_OKAY;
}



/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpIndicator)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   if ( nconss > 0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int nGen = 0;

      SCIPdebugMessage("Separating inequalities for indicator constraints.\n");

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );

      if ( conshdlrdata->sepaAlternativeLP )
      {
	 *result = SCIP_DIDNOTFIND;

	 /* start separation */
	 SCIP_CALL( separateIISRounding(scip, conshdlr, NULL, nconss, conss, &nGen) );
	 SCIPdebugMessage("Separated %d cuts from indicator constraints.\n", nGen);

	 if ( nGen > 0 )
	 {
	    if ( conshdlrdata->genLogicor )
	       *result = SCIP_CONSADDED;
	    else
	       *result = SCIP_SEPARATED;
	 }
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolIndicator)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   if ( nconss > 0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int nGen = 0;

      SCIPdebugMessage("Separating inequalities for indicator constraints.\n");

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );

      if ( conshdlrdata->sepaAlternativeLP )
      {
	 *result = SCIP_DIDNOTFIND;
	 /* start separation */
	 SCIP_CALL( separateIISRounding(scip, conshdlr, sol, nconss, conss, &nGen) );
	 SCIPdebugMessage("Separated %d cuts from indicator constraints.\n", nGen);

	 if ( nGen > 0 )
	 {
	    if ( conshdlrdata->genLogicor )
	       *result = SCIP_CONSADDED;
	    else
	       *result = SCIP_SEPARATED;
	 }
      }
   }

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpIndicator)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( enforceIndicators(scip, conshdlr, nconss, conss, result) );

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsIndicator)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( enforceIndicators(scip, conshdlr, nconss, conss, result) );

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

         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, "violation:  <%s> = %g and <%s> = %.15g\n",
               SCIPvarGetName(consdata->binvar), SCIPgetSolVal(scip, sol, consdata->binvar),
               SCIPvarGetName(consdata->slackvar), SCIPgetSolVal(scip, sol, consdata->slackvar));
         }

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
 *  @a inferinfo is 1 if the slack variable * has bounds that fix it to
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

   SCIPinfoMessage(scip, file, "<%s> = 1", SCIPvarGetName(consdata->binvar));
   SCIPinfoMessage(scip, file, " -> <%s> = 0)\n", SCIPvarGetName(consdata->slackvar));

   return SCIP_OKAY;
}


/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableIndicator)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMessage("Enabling constraint <%s>.\n", SCIPconsGetName(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->altLP != NULL )
   {
      SCIP_CONSDATA* consdata;
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( conshdlrdata->sepaAlternativeLP );

      SCIP_CALL( unfixAltLPVariable(conshdlrdata->altLP, consdata->colIndex) );
   }

   return SCIP_OKAY;
}


/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableIndicator)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMessage("Disabling constraint <%s>.\n", SCIPconsGetName(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->altLP != NULL )
   {
      SCIP_CONSDATA* consdata;
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( conshdlrdata->sepaAlternativeLP );

      SCIP_CALL( fixAltLPVariable(conshdlrdata->altLP, consdata->colIndex) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveIndicator NULL

/** constraint deactivation notification method of constraint handler */
#define consDeactiveIndicator NULL

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
      SCIPdebugMessage("changed lower bound of variable <%s> from %f to %f (nFixedNonzero: %d).\n",
	 SCIPvarGetName(SCIPeventGetVar(event)), oldbound, newbound, consdata->nFixedNonzero);
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisFeasNegative(scip, oldbound) && SCIPisFeasNegative(scip, newbound) )
	 ++(consdata->nFixedNonzero);
      SCIPdebugMessage("changed upper bound of variable <%s> from %f to %f (nFixedNonzero: %d).\n",
	 SCIPvarGetName(SCIPeventGetVar(event)), oldbound, newbound, consdata->nFixedNonzero);
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisFeasPositive(scip, oldbound) && ! SCIPisFeasPositive(scip, newbound) )
	 --(consdata->nFixedNonzero);
      SCIPdebugMessage("changed lower bound of variable <%s> from %f to %f (nFixedNonzero: %d).\n",
	 SCIPvarGetName(SCIPeventGetVar(event)), oldbound, newbound, consdata->nFixedNonzero);
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisFeasNegative(scip, oldbound) && ! SCIPisFeasNegative(scip, newbound) )
	 --(consdata->nFixedNonzero);
      SCIPdebugMessage("changed upper bound of variable <%s> from %f to %f (nFixedNonzero: %d).\n",
	 SCIPvarGetName(SCIPeventGetVar(event)), oldbound, newbound, consdata->nFixedNonzero);
      break;
   default:
      SCIPerrorMessage("invalid event type.\n");
      return SCIP_INVALIDDATA;
   }
   assert( 0 <= consdata->nFixedNonzero && consdata->nFixedNonzero <= 2 );

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
   conshdlrdata->removable = TRUE;
   conshdlrdata->altLP = NULL;
   conshdlrdata->nvars = 0;
   conshdlrdata->varHash = NULL;
   conshdlrdata->lbHash = NULL;
   conshdlrdata->ubHash = NULL;
   conshdlrdata->nLbBounds = 0;
   conshdlrdata->nUbBounds = 0;
   conshdlrdata->roundingMinThreshold =	0.1;
   conshdlrdata->roundingMaxThreshold =	0.6;
   conshdlrdata->roundingRounds = 1;
   conshdlrdata->roundingOffset = 0.1;
   conshdlrdata->branchIndicators = TRUE;
   conshdlrdata->genLogicor = TRUE;
   conshdlrdata->sepaAlternativeLP = TRUE;

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

   /* add indicator constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/indicator/branchIndicators",
         "Branch on indicator constraints in enforcing?",
         &conshdlrdata->branchIndicators, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/indicator/genLogicor",
         "Generate logicor constraints instead of cuts?",
         &conshdlrdata->genLogicor, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/indicator/sepaAlternativeLP",
         "Separate using the alternative LP?",
         &conshdlrdata->sepaAlternativeLP, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/indicator/addCoupling",
         "add initial coupling inequalities",
         &conshdlrdata->addCoupling, TRUE, FALSE, NULL, NULL) );

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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
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
   consdata->colIndex = 0;

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



/** gets binary variable corresponding to indicator constraint */
SCIP_VAR* SCIPgetBinaryVarIndicator(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->binvar;
}


/** gets slack variable corresponding to indicator constraint */
SCIP_VAR* SCIPgetSlackVarIndicator(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->slackvar;
}


/** checks whether indicator constraint is violated w.r.t. sol */
SCIP_Bool SCIPisViolatedIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   )
{
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return(
      SCIPisFeasPositive(scip, SCIPgetSolVal(scip, sol, consdata->slackvar)) &&
      SCIPisFeasPositive(scip, SCIPgetSolVal(scip, sol, consdata->binvar)) );
}
