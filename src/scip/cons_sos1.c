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
#pragma ident "@(#) $Id: cons_sos1.c,v 1.1 2007/10/12 20:47:47 bzfpfets Exp $"
//#define SCIP_DEBUG

/**@file   cons_sos1.c
 * @brief  constraint handler for SOS type 1 constraints
 * @author Marc Pfetsch
 *
 * A specially ordered set of type 1 (SOS1) is a set of variables such
 * that at most one variable is nonzero. The special case of two
 * variables arises, for instance, from equilibrium or complementary
 * conditions like x * y = 0.
 *
 * This implementation of this constraint handler is based on classical ideas, see e.g.@n
 *  "Special Facilities in General Mathematical Programming System for
 *  Non-Convex Problems Using Ordered Sets of Variables"@n
 *  E. Beale and J. Tomlin, Proc. 5th IFORS Conference, 447-454 (1970)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_sos1.h"
#include "scip/cons_linear.h"
#include <string.h>


/* constraint handler properties */
#define CONSHDLR_NAME          "SOS1"
#define CONSHDLR_DESC          "SOS1 constraint handler"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* event handler properties */
#define EVENTHDLR_NAME         "SOS1"
#define EVENTHDLR_DESC         "bound change event handler for SOS1 constraints"





/** constraint data for SOS1 constraints */
struct SCIP_ConsData
{
   int nVars;              /**< number of variables in the constraint */
   int maxVars;            /**< maximal number of variables */
   SCIP_VAR** Vars;        /**< variables in constraint */
   SCIP_ROW* row;          /**< row corresponding to upper and lower bound inequalities */
   int nFixedNonzero;      /**< number of variables fixed to be nonzero */
};

/** SOS1 constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR* eventhdlr;     /**< event handler for bound change events */
};






/** fix variable in given node to 0 */
static
SCIP_RETCODE fixVariableZeroNode(
	     SCIP* scip,            /**< SCIP pointer */
	     SCIP_VAR* var,         /**< variable to be fixed to 0*/
	     SCIP_NODE* node,       /**< node */
	     SCIP_Bool* infeasible  /**< if fixing is infeasible */
	     )
{
   /* if variable cannot be nonzero */
   *infeasible = FALSE;
   if ( SCIPisPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisNegative(scip, SCIPvarGetUbLocal(var)) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* if variable is multi aggregated */
   if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_CONS* cons;
      SCIP_Real val = 1.0;

      if ( ! SCIPisZero(scip, SCIPvarGetLbLocal(var)) || ! SCIPisZero(scip, SCIPvarGetUbLocal(var)) )
      {
	 SCIPdebugMessage("creating constraint to force multi aggregated variable <%s> to 0.\n", SCIPvarGetName(var));
	 /* we have to insert a local constraint var = 0 */
	 SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &var, &val, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
					 TRUE, FALSE, FALSE, FALSE, FALSE) );
	 SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
	 SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }
   else
   {
      if ( ! SCIPisZero(scip, SCIPvarGetLbLocal(var)) )
	 SCIP_CALL( SCIPchgVarLbNode(scip, node, var, 0.0) );
      if ( ! SCIPisZero(scip, SCIPvarGetUbLocal(var)) )
	 SCIP_CALL( SCIPchgVarUbNode(scip, node, var, 0.0) );
   }

   return SCIP_OKAY;
}


/** fix variable in local node to 0 */
static
SCIP_RETCODE inferVariableZero(
	     SCIP* scip,            /**< SCIP pointer */
	     SCIP_VAR* var,         /**< variable to be fixed to 0*/
	     SCIP_CONS* cons,       /**< constraint */
	     int inferinfo,         /**< info for reverse prop. */
	     SCIP_Bool* infeasible, /**< if fixing is infeasible */
	     SCIP_Bool* tightened   /**< if fixing was performed */
	     )
{
   *infeasible = FALSE;
   *tightened = FALSE;

   /* if variable cannot be nonzero */
   if ( SCIPisPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisNegative(scip, SCIPvarGetUbLocal(var)) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* do nothing if variable is multi aggregated */
   if ( SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_Bool tighten;

      /* fix lower bound */
      SCIP_CALL( SCIPinferVarLbCons(scip, var, 0.0, cons, inferinfo, infeasible, &tighten) );
      *tightened = *tightened || tighten;

      /* fix upper bound */
      SCIP_CALL( SCIPinferVarUbCons(scip, var, 0.0, cons, inferinfo, infeasible, &tighten) );
      *tightened = *tightened || tighten;
   }

   return SCIP_OKAY;
}




/** fix all variables except given one to 0 */
static
SCIP_RETCODE fixOtherZero(
	     SCIP* scip,           /**< SCIP pointer */
	     SCIP_CONS* cons,      /**< constraint */
	     int nVars,            /**< number of variables */
	     SCIP_VAR** Vars,      /**< variables in constraint */
	     int ind,              /**< index of variable which is nonzero */
	     SCIP_RESULT* result   /**< result */
	     )
{
   SCIP_Bool tightened = FALSE;
   int j;
   assert( SCIPisPositive(scip, SCIPvarGetLbLocal(Vars[ind])) || SCIPisNegative(scip, SCIPvarGetUbLocal(Vars[ind])) );

   for (j = 0; j < nVars; ++j)
   {
      if ( j != ind )
      {
	 SCIP_Bool infeasible, tighten;

	 /* fix lower bound */
	 SCIP_CALL( SCIPinferVarLbCons(scip, Vars[j], 0.0, cons, 0, &infeasible, &tighten) );
	 if ( infeasible )
	 {
	    *result = SCIP_CUTOFF;
	    return SCIP_OKAY;
	 }
	 tightened = tightened || tighten;

	 /* fix upper bound */
	 SCIP_CALL( SCIPinferVarUbCons(scip, Vars[j], 0.0, cons, 0, &infeasible, &tighten) );
	 if ( infeasible )
	 {
	    *result = SCIP_CUTOFF;
	    return SCIP_OKAY;
	 }
	 tightened = tightened || tighten;
      }
   }
   if ( tightened )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}



/** enforcement method
 *
 *  We check whether the current solution is feasible, i.e., contains
 *  at most one nonzero variable. If not, we branch along the lines
 *  indicated by Beale and Tomlin:
 *
 *  We first compute \f$W = \sum_{j=1}^n |x_i|\f$ and \f$w =
 *  \sum_{j=1}^n j\, |x_i|\f$. Then we search for the index k that
 *  satisfies
 *  \f[
 *        k \leq \frac{w}{W} < k+1.
 *  \f]
 *  The branches are then
 *  \f[
 *        x_1 = 0, \ldots, x_{k-1} = 0 \qquad \mbox{and}\qquad
 *        x_{k+1} = 0, \ldots, x_n = 0.
 *  \f]
 */
static
SCIP_RETCODE enforceSOS1(
			 SCIP* scip,           /**< SCIP pointer */
			 SCIP_CONS* cons,      /**< constraint */
			 int nVars,            /**< number of variables */
			 SCIP_VAR** Vars,      /**< variables in constraint */
			 SCIP_RESULT* result   /**< result */
			 )
{
   SCIP_Bool infeasible, tightened;
   SCIP_NODE* node1;
   SCIP_NODE* node2;
   int j;
   SCIP_Real weight1 = 0.0;
   SCIP_Real weight2 = 0.0;
   SCIP_Real w = 0.0;
   int ind;
   int cnt = 0;

   assert( scip != NULL );
   assert( Vars != NULL );
   assert( result != NULL );

   if ( nVars == 2 )
   {
      /* if the bounds of the first variable are such that it cannot be 0 */
      if ( SCIPisPositive(scip, SCIPvarGetLbLocal(Vars[0])) || SCIPisNegative(scip, SCIPvarGetUbLocal(Vars[0])) )
      {
	 SCIPdebugMessage("Variable <%s> is fixed to be nonzero.\n", SCIPvarGetName(Vars[0]));
	 SCIP_CALL( inferVariableZero(scip, Vars[1], cons, 0, &infeasible, &tightened) );
	 if ( infeasible )
	 {
	    *result = SCIP_CUTOFF;
	    return SCIP_OKAY;
	 }
	 if ( tightened )
	 {
	    *result = SCIP_REDUCEDDOM;
	    return SCIP_OKAY;
	 }
      }

      /* if the bounds of the second variable are such that it cannot be 0 */
      if ( SCIPisPositive(scip, SCIPvarGetLbLocal(Vars[1])) || SCIPisNegative(scip, SCIPvarGetUbLocal(Vars[1])) )
      {
	 SCIPdebugMessage("Variable <%s> is fixed to be nonzero.\n", SCIPvarGetName(Vars[1]));
	 SCIP_CALL( inferVariableZero(scip, Vars[0], cons, 1, &infeasible, &tightened) );
	 if ( infeasible )
	 {
	    *result = SCIP_CUTOFF;
	    return SCIP_OKAY;
	 }
	 if ( tightened )
	 {
	    *result = SCIP_REDUCEDDOM;
	    return SCIP_OKAY;
	 }
      }

      /* if constraint is infeasible */
      if ( ! SCIPisZero(scip, SCIPgetSolVal(scip, NULL, Vars[0])) && ! SCIPisZero(scip, SCIPgetSolVal(scip, NULL, Vars[1])) )
      {
	 /* create branches */
	 SCIPdebugMessage("Creating two branches.\n");

	 SCIP_CALL( SCIPcreateChild(scip, &node1, 0.0, SCIPgetSolTransObj(scip, NULL) ) );
	 SCIP_CALL( fixVariableZeroNode(scip, Vars[0], node1, &infeasible) );
	 assert( ! infeasible );

	 SCIP_CALL( SCIPcreateChild(scip, &node2, 0.0, SCIPgetSolTransObj(scip, NULL) ) );
	 SCIP_CALL( fixVariableZeroNode(scip, Vars[1], node2, &infeasible) );
	 assert( ! infeasible );

	 *result = SCIP_BRANCHED;
      }
   }
   else
   {
      /* compute weight */
      for (j = 0; j < nVars; ++j)
      {
	 SCIP_Real val = fabs(SCIPgetSolVal(scip, NULL, Vars[j]));
	 weight1 += val * (SCIP_Real) j;
	 weight2 += val;
	 if ( ! SCIPisFeasZero(scip, val) )
	    ++cnt;

	 /* if the bounds of the variable are such that it cannot be 0 */
	 if ( SCIPisPositive(scip, SCIPvarGetLbLocal(Vars[j])) || SCIPisNegative(scip, SCIPvarGetUbLocal(Vars[j])) )
	 {
	    SCIP_CALL( fixOtherZero(scip, cons, nVars, Vars, j, result) );
	    if (*result == SCIP_REDUCEDDOM || *result == SCIP_CUTOFF )
	       return SCIP_OKAY;
	 }
      }

      /* if at most one variable is nonzero, the constraint is feasible -> return */
      if ( cnt <= 1 )
	 return SCIP_OKAY;

      assert( cnt >= 2 );
      assert( !SCIPisFeasZero(scip, weight2) );
      w = weight1/weight2;

      ind = SCIPfloor(scip, w);
      if ( ind >= nVars )
	 ind = nVars;
      assert( 0 <= ind && ind < nVars );

      /* create branches */
      SCIPdebugMessage("Branching on variable %d.\n", ind);

      /* branch on variable ind: either all variables before ind or all variables after ind are zero */
      SCIP_CALL( SCIPcreateChild(scip, &node1, 0.0, SCIPgetSolTransObj(scip, NULL) ) );
      for (j = 0; j < ind; ++j)
      {
	 SCIP_CALL( fixVariableZeroNode(scip, Vars[j], node1, &infeasible) );
	 assert( ! infeasible );
      }

      SCIP_CALL( SCIPcreateChild(scip, &node2, 0.0, SCIPgetSolTransObj(scip, NULL) ) );
      for (j = ind+1; j < nVars; ++j)
      {
	 SCIP_CALL( fixVariableZeroNode(scip, Vars[j], node2, &infeasible) );
	 assert( ! infeasible );
      }

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}



/** Generate row
 *
 *  We generate the rows corresponding to the following simple valid inequalities:
 *  \f[
 *         x_1 + \ldots + x_n \leq \max\{u_1, \ldots, u_n\}\qquad\mbox{and}\qquad
 *         x_1 + \ldots + x_n \geq \min\{l_1, \ldots, l_n\},
 *  \f]
 *  where \f$l_1, \ldots, l_n\f$ and \f$u_1, \ldots, u_n\f$ are the
 *  lower and upper bounds of the variables \f$x_1, \ldots, x_n\f$.
 *  Of course, these inequalities are only added if the upper and
 *  lower bounds are all finite and at least one lower bound is less
 *  than 0 (lower bounds > 0 are detected in presolving).
 */
static
SCIP_RETCODE generateRowsSOS1(
		  SCIP* scip,               /**< SCIP pointer */
		  SCIP_CONSDATA* consdata   /**< constraint data */
		  )
{
   SCIP_VAR** Vars;
   SCIP_Real minLb = SCIPinfinity(scip);
   SCIP_Real maxUb = -SCIPinfinity(scip);
   SCIP_ROW* row;
   int j, nVars;

   assert( scip != NULL );
   assert( consdata != NULL );

   nVars = consdata->nVars;
   Vars = consdata->Vars;
   assert( Vars != NULL );

   /* find minimum and maximum lower and upper bounds */
   for (j = 0; j < nVars; ++j)
   {
      minLb = MIN( SCIPvarGetLbLocal(Vars[j]), minLb );
      maxUb = MAX( SCIPvarGetUbLocal(Vars[j]), maxUb );
   }
   /* ignore trivial inequality if all lower bounds are 0 */
   if ( ! SCIPisNegative(scip, minLb) )
      minLb = -SCIPinfinity(scip);

   /* create upper and lower bound inequality */
   if ( maxUb < SCIPinfinity(scip) || minLb > -SCIPinfinity(scip) )
   {
      SCIP_CALL( SCIPcreateEmptyRow(scip, &row, "sosbnd", minLb, maxUb, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
      for (j = 0; j < nVars; ++j)
	 SCIP_CALL( SCIPaddVarToRow(scip, row, Vars[j], 1.0) );
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
      consdata->row = row;
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
   }

   return SCIP_OKAY;
}





/** propagate changed variables */
static
SCIP_RETCODE propSOS1(SCIP* scip, SCIP_CONS* cons, SCIP_CONSDATA* consdata, SCIP_RESULT* result, int* nGen)
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( result != NULL );
   assert( nGen != NULL );

   /* if more than one variable is fixed to be nonzero */
   if ( consdata->nFixedNonzero > 1 )
   {
      SCIPdebugMessage("The node is infeasible, more than 1 variable is fixed to be nonzero.\n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* if exactly one variable is fixed to be nonzero */
   if ( consdata->nFixedNonzero == 1 )
   {
      int j, nVars, firstFixedNonzero = -1;
      SCIP_VAR** Vars;

      nVars = consdata->nVars;
      Vars = consdata->Vars;
      assert( Vars != NULL );

      /* search nonzero variable - it is needed for propinfo */
      for (j = 0; j < nVars; ++j)
      {
	 if ( SCIPisPositive(scip, SCIPvarGetLbLocal(Vars[j])) || SCIPisNegative(scip, SCIPvarGetUbLocal(Vars[j])) )
	 {
	    firstFixedNonzero = j;
	    break;
	 }
      }
      assert( firstFixedNonzero >= 0 );

      SCIPdebugMessage("variable <%s> is nonzero, fixing other variables to 0.\n", SCIPvarGetName(Vars[firstFixedNonzero]));

      /* fix variables before firstFixedNonzero to 0 */
      for (j = 0; j < firstFixedNonzero; ++j)
      {
	 SCIP_Bool infeasible, tightened;

	 /* fix variable */
	 SCIP_CALL( inferVariableZero(scip, Vars[j], cons, firstFixedNonzero, &infeasible, &tightened) );
	 assert( ! infeasible );

	 if ( tightened )
	    ++(*nGen);
      }

      /* fix variables after firstFixedNonzero to 0 */
      for (j = firstFixedNonzero+1; j < nVars; ++j)
      {
	 SCIP_Bool infeasible, tightened;

	 /* fix variable */
	 SCIP_CALL( inferVariableZero(scip, Vars[j], cons, firstFixedNonzero, &infeasible, &tightened) );
	 /* the node can be infeasible, since we did not check variables after firstFixedNonzero */
	 if ( infeasible )
	 {
	    SCIPdebugMessage("The node is infeasible.\n");
	    *result = SCIP_CUTOFF;
	    return SCIP_OKAY;
	 }

	 if ( tightened )
	    ++(*nGen);
      }
   }
   return SCIP_OKAY;
}





/* ---------------------------- constraint handler callback methods ----------------------*/

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeSOS1)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}



/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolSOS1)
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

      SCIPdebugMessage("Initializing SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* generate rows */
      generateRowsSOS1(scip, consdata);
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolSOS1)
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

      SCIPdebugMessage("Exiting SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* free row */
      if ( consdata->row != NULL )
	 SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
   }
   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSOS1)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( (*consdata)->nVars > 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMessage("Deleting SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

   /* drop events on transfromed variables */
   if ( (*consdata)->nVars > 0 && SCIPvarIsTransformed((*consdata)->Vars[0]) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get constraint handler data */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );
      assert( conshdlrdata->eventhdlr != NULL );

      int j;
      for (j = 0; j < (*consdata)->nVars; ++j)
	 SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->Vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)*consdata, -1) );
   }


   SCIPfreeBlockMemoryArray(scip, &(*consdata)->Vars, (*consdata)->maxVars);
   assert( (*consdata)->row == NULL );
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSOS1)
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   char s[SCIP_MAXSTRLEN];
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   assert( conshdlrdata->eventhdlr != NULL );

   SCIPdebugMessage("Transforming SOS1 constraint: <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->nVars > 0 );

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nVars = sourcedata->nVars;
   consdata->maxVars = sourcedata->maxVars;
   consdata->row = NULL;
   consdata->nFixedNonzero = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->Vars, consdata->nVars) );
   for (j = 0; j < sourcedata->nVars; ++j)
   {
      assert( sourcedata->Vars[j] != 0 );
      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->Vars[j], &(consdata->Vars[j])) );
      /* if variable is fixed to be nonzero */
      if ( SCIPisPositive(scip, SCIPvarGetLbLocal(consdata->Vars[j])) || SCIPisNegative(scip, SCIPvarGetUbLocal(consdata->Vars[j])) )
	 ++(consdata->nFixedNonzero);
   }

   /* create transformed constraint with the same flags */
   snprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
			     SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
			     SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
			     SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
			     SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
			     SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* catch bound change events on variable */
   for (j = 0; j < consdata->nVars; ++j)
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->Vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );
   }

#ifdef SCIP_DEBUG
   if ( consdata->nFixedNonzero > 0 )
   {
      SCIPdebugMessage("constraint <%s> has %d variables fixed to be nonzero.\n", SCIPconsGetName(*targetcons), consdata->nFixedNonzero );
   }
#endif

   return SCIP_OKAY;
}




/* remove lock on variable */
static
SCIP_RETCODE unlockVariableSOS1(SCIP* scip, SCIP_VAR* var, SCIP_CONS* cons)
{
   if ( SCIPisFeasZero(scip, SCIPvarGetLbGlobal(var)) )
   {
      /* rounding up == bad, rounding down = o.k. */
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, FALSE, TRUE) );
   }
   else
   {
      if ( ! SCIPisPositive(scip, SCIPvarGetUbGlobal(var)) )
      {
	 /* rounding up == o.k., rounding down = bad */
	 SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, FALSE) );
      }
      else
      {
	 /* any rounding == bad */
	 SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, TRUE) );
      }
   }
   return SCIP_OKAY;
}



/** presolving method of constraint handler
 *
 *  We perform the following presolving steps.
 *
 *  - If the bounds of some variable force it to be nonzero, we can
 *    fix all other variables to zero and remove the SOS1 constraints
 *    that contain it.
 *  - If a variable is fixed to zero, we can remove the constraint.
 */
static
SCIP_DECL_CONSPRESOL(consPresolSOS1)
{
   int c;
   int oldnfixedvars = 0;
   int oldndelconss = 0;
   int removedvars = 0;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;
   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );
   assert( conshdlrdata->eventhdlr != NULL );

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
      assert( ! SCIPconsIsModifiable(cons) );

      SCIPdebugMessage("Presolving SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

      *result = SCIP_DIDNOTFIND;

      /* only run if sucess if possible */
      if ( nrounds == 0 || nnewfixedvars > 0 || nnewaggrvars > 0 || *nfixedvars > oldnfixedvars )
      {
	 SCIP_VAR** Vars;
	 int nFixedNonzero = 0;
	 int lastFixedNonzero = -1;
	 int nVars;
	 int j = 0;

	 nVars = consdata->nVars;
	 Vars = consdata->Vars;

	 /* check for variables fixed to 0 and bounds that fix a variable to be nonzero */
	 while (j < nVars)
	 {
	    SCIP_Real lb,ub;

	    lb = SCIPvarGetLbGlobal(Vars[j]);
	    ub = SCIPvarGetUbGlobal(Vars[j]);

	    /* if the variable if fixed to nonzero */
	    if ( SCIPisPositive(scip, lb) || SCIPisNegative(scip, ub) )
	    {
	       ++nFixedNonzero;
	       lastFixedNonzero = j;
	    }

	    /* if the variable is fixed to 0 */
	    if ( SCIPisFeasZero(scip, lb) && SCIPisFeasZero(scip, ub) )
	    {
	       int l;
	       /* remove lock of variable */
	       SCIP_CALL( unlockVariableSOS1(scip, Vars[j], cons) );

	       /* drop events on variable */
	       SCIP_CALL( SCIPdropVarEvent(scip, Vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

	       /* delete variable - need to copy since order is important */
	       ++removedvars;
	       for (l = j; l < nVars-1; ++l)
		  Vars[l] = Vars[l+1];
	       --nVars;
	    }
	    else
	       ++j;
	 }
	 consdata->nVars = nVars;

	 /* if the number of variables is less than 2 */
	 if ( nVars < 2 )
	 {
	    SCIPdebugMessage("Deleting constraint with < 2 variables.\n");
	    /* delete constraint */
	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);
	    *result = SCIP_SUCCESS;
	    continue;
	 }

	 /* if more than one variable are fixed to be nonzero, we are infeasible */
	 if ( nFixedNonzero > 1 )
	 {
	    SCIPdebugMessage("The problem is infeasible: more than one variable have bounds that keep it from being 0.\n");
	    assert( lastFixedNonzero >= 0 );
	    *result = SCIP_CUTOFF;
	    return SCIP_OKAY;
	 }

	 /* if there is exactly one fixed nonzero variable */
	 if ( nFixedNonzero == 1 )
	 {
	    assert( lastFixedNonzero >= 0 );
	    /* fix all other variables to zero */
	    for (j = 0; j < nVars; ++j)
	    {
	       if ( j != lastFixedNonzero )
	       {
		  SCIP_Bool infeasible, fixed;

		  SCIP_CALL( SCIPfixVar(scip, Vars[j], 0.0, &infeasible, &fixed) );
		  assert( ! infeasible );
		  if ( fixed )
		     ++(*nfixedvars);
	       }
	    }
	    /* delete original constraint */
	    SCIP_CALL( SCIPdelCons(scip, cons) );
	    ++(*ndelconss);
	    *result = SCIP_SUCCESS;
	 }
      }
   }

   SCIPdebugMessage("presolving fixed %d variables, removed %d variables, and deleted %d constraints.\n",
		    *nfixedvars - oldnfixedvars, removedvars, *ndelconss - oldndelconss);

   return SCIP_OKAY;
}




/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpSOS1)
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

      SCIPdebugMessage("Checking for initial rows for SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* put corresponding rows into LP if they are usefull */
      if ( consdata->row != NULL )
      {
	 SCIP_CALL( SCIPaddCut(scip, NULL, consdata->row, TRUE) );
#ifdef SCIP_DEBUG
	 SCIP_CALL( SCIPprintRow(scip, consdata->row, NULL) );
#endif
      }
   }

   return SCIP_OKAY;
}



/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpSOS1)
{
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
      SCIP_ROW* row;

      *result = SCIP_DIDNOTFIND;
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMessage("Separating inequalities for SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* put corresponding rows into LP if they are usefull */
      row = consdata->row;
      if ( row != NULL && ! SCIProwIsInLP(row) && SCIPisCutEfficacious(scip, NULL, row) )
      {
	 SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE) );
#ifdef SCIP_DEBUG
	 SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
	 ++nGen;
      }
   }
   SCIPdebugMessage("Separated %d SOS1 constraints.\n", nGen);
   if ( nGen > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolSOS1)
{
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
      SCIP_ROW* row;

      *result = SCIP_DIDNOTFIND;
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMessage("Separating solution for SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* put corresponding row into LP if it is usefull */
      row = consdata->row;
      if ( row != NULL && SCIPisCutEfficacious(scip, sol, row) )
      {
	 SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );
#ifdef SCIP_DEBUG
	 SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
	 ++nGen;
      }
   }
   SCIPdebugMessage("Separated %d SOS1 constraints.\n", nGen);
   if ( nGen > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpSOS1)
{
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
      SCIPdebugMessage("Enforcing SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      SCIP_CALL( enforceSOS1(scip, conss[c], consdata->nVars, consdata->Vars, result) );

      if ( *result == SCIP_BRANCHED )
	 return SCIP_OKAY;
   }
   SCIPdebugMessage("All SOS1 constraints are feasible.\n");

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSOS1)
{
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
      SCIPdebugMessage("Enforcing SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      SCIP_CALL( enforceSOS1(scip, conss[c], consdata->nVars, consdata->Vars, result) );

      if ( *result == SCIP_BRANCHED )
	 return SCIP_OKAY;
   }
   SCIPdebugMessage("All SOS1 constraints are feasible.\n");

   return SCIP_OKAY;
}




/** feasibility check method of constraint handler for integral solutions
 *
 *  We simply check whether at most one variable is nonzero in the given solution.
 */
static
SCIP_DECL_CONSCHECK(consCheckSOS1)
{
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
      int j;
      int cnt = 0;

      *result = SCIP_DIDNOTFIND;
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMessage("Checking SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* check all variables */
      for (j = 0; j < consdata->nVars; ++j)
      {
	 /* if variable is nonzero */
	 if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->Vars[j])) )
	 {
	    /* if more than one variable is nonzero */
	    if ( cnt > 0 )
	    {
	       *result = SCIP_INFEASIBLE;
	       return SCIP_OKAY;
	    }
	    ++cnt;
	 }
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler
 *
 *  If the bounds for one variable fix it to be nonzero, the other variables can be fixed to 0.
 */
static
SCIP_DECL_CONSPROP(consPropSOS1)
{
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

      *result = SCIP_DIDNOTFIND;
      assert( conss[c] != NULL );
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      SCIPdebugMessage("Propagating SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( propSOS1(scip, cons, consdata, result, &nGen) );
   }
   SCIPdebugMessage("Propagated %d domains.\n", nGen);
   if ( nGen > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}




/** propagation conflict resolving method of constraint handler
 *
 *  We check which bound changes were the reason for infeasibility. We
 *  use that @a inferinfo stores the index of the variable that has
 *  bounds that fix it to be nonzero (these bounds are the reason). */
static
SCIP_DECL_CONSRESPROP(consRespropSOS1)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;
   SCIPdebugMessage("Propagation resolution method of SOS1 constraint <%s>.\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( 0 <= inferinfo && inferinfo < consdata->nVars );
   var = consdata->Vars[inferinfo];
   assert( var != infervar );

   /* check if lower bound was changed */
   if ( SCIPisLT(scip, SCIPvarGetLbAtIndex(var, bdchgidx, FALSE), SCIPvarGetLbAtIndex(var, bdchgidx, TRUE)) )
   {
      SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   /* check if upper bound was changed */
   if ( SCIPisLT(scip, SCIPvarGetUbAtIndex(var, bdchgidx, FALSE), SCIPvarGetUbAtIndex(var, bdchgidx, TRUE)) )
   {
      SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}




/** variable rounding lock method of constraint handler
 *
 *  Let lb and ub be the lower and upper bounds of a
 *  variable. Preprocessing makes sure that lb <= 0 <= ub.
 *
 *  - If lb == 0 then rounding up the variable may violate the
 *    constraint, but rounding down does not (as this will yield 0).
 *  - If ub == 0 then rounding down the variable may violate the
 *    constraint, but rounding up does not (as this will yield 0).
 *  - In all other cases rounding the variable in any direction
 *    may violated the constraint.
 */
static
SCIP_DECL_CONSLOCK(consLockSOS1)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** Vars;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIPdebugMessage("Locking constraint <%s>.\n", SCIPconsGetName(cons));

   Vars = consdata->Vars;
   assert( Vars != NULL );

   for (j = 0; j < consdata->nVars; ++j)
   {
      if ( SCIPisFeasZero(scip, SCIPvarGetLbGlobal(Vars[j])) )
      {
	 /* rounding up == bad, rounding down = o.k. */
	 SCIP_CALL( SCIPaddVarLocks(scip, consdata->Vars[j], nlocksneg, nlockspos) );
      }
      else
      {
	 if ( ! SCIPisPositive(scip, SCIPvarGetUbGlobal(Vars[j])) )
	 {
	    /* rounding up == o.k., rounding down = bad */
	    SCIP_CALL( SCIPaddVarLocks(scip, consdata->Vars[j], nlockspos, nlocksneg) );
	 }
	 else
	 {
	    /* any rounding == bad */
	    SCIP_CALL( SCIPaddVarLocks(scip, consdata->Vars[j], nlockspos + nlocksneg, nlockspos + nlocksneg) );
	 }
      }
   }

   return SCIP_OKAY;
}



/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintSOS1)
{
   SCIP_CONSDATA* consdata;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIPinfoMessage(scip, file, "SOS1[");
   for (j = 0; j < consdata->nVars; ++j)
   {
      if ( j > 0 )
	 SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "%s", SCIPvarGetName(consdata->Vars[j]));
   }
   SCIPinfoMessage(scip, file, "]\n");

   return SCIP_OKAY;
}






/** constraint activation notification method of constraint handler */
#define consActiveSOS1 NULL

/** constraint deactivation notification method of constraint handler */
#define consDeactiveSOS1 NULL

/** constraint enabling notification method of constraint handler */
#define consEnableSOS1 NULL

/** constraint disabling notification method of constraint handler */
#define consDisableSOS1 NULL

/** initialization method of constraint handler (called after problem was transformed) */
#define consInitSOS1 NULL

/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitSOS1 NULL

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreSOS1 NULL

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreSOS1 NULL





/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * We update the number of variables fixed to be nonzero
 */
static
SCIP_DECL_EVENTEXEC(eventExecSOS1)
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
   assert( 0 <= consdata->nFixedNonzero && consdata->nFixedNonzero <= consdata->nVars );

   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);

   eventtype = SCIPeventGetType(event);
   switch( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisPositive(scip, oldbound) && SCIPisPositive(scip, newbound) )
	 ++(consdata->nFixedNonzero);
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisNegative(scip, oldbound) && SCIPisNegative(scip, newbound) )
	 ++(consdata->nFixedNonzero);
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisPositive(scip, oldbound) && ! SCIPisPositive(scip, newbound) )
	 --(consdata->nFixedNonzero);
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisNegative(scip, oldbound) && ! SCIPisNegative(scip, newbound) )
	 --(consdata->nFixedNonzero);
      break;
   default:
      SCIPerrorMessage("invalid event type.\n");
      return SCIP_INVALIDDATA;
   }
   assert( 0 <= consdata->nFixedNonzero && consdata->nFixedNonzero <= consdata->nVars );

   SCIPdebugMessage("changed bound of variable <%s> from %f to %f (nFixedNonzero: %d).\n", SCIPvarGetName(SCIPeventGetVar(event)),
		    oldbound, newbound, consdata->nFixedNonzero);

   return SCIP_OKAY;
}



/*
 * constraint specific interface methods
 */

/** creates the handler for SOS1 constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSOS1(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
			 NULL, NULL, NULL, NULL, NULL, NULL, eventExecSOS1, NULL) );

   /* create constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* get event handler for bound change events */
   conshdlrdata->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if ( conshdlrdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for SOS1 constraints not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
			  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
			  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
			  CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
			  consFreeSOS1, consInitSOS1, consExitSOS1,
			  consInitpreSOS1, consExitpreSOS1, consInitsolSOS1, consExitsolSOS1,
			  consDeleteSOS1, consTransSOS1, consInitlpSOS1,
			  consSepalpSOS1, consSepasolSOS1, consEnfolpSOS1, consEnfopsSOS1, consCheckSOS1,
			  consPropSOS1, consPresolSOS1, consRespropSOS1, consLockSOS1,
			  consActiveSOS1, consDeactiveSOS1,
			  consEnableSOS1, consDisableSOS1,
			  consPrintSOS1,
			  conshdlrdata) );

   return SCIP_OKAY;
}


/** creates and captures a SOS1 constraint
 *
 *  We set the constraint to not be modifable and non be sticking at a node
 */
SCIP_RETCODE SCIPcreateConsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to TRUE. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   if ( nvars < 2 )
   {
      SCIPerrorMessage("Need at least 2 variables.\n");
      return SCIP_ERROR;
   }

   /* find the SOS1 constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("<%s> constraint handler not found\n", CONSHDLR_NAME);
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIPallocBlockMemory(scip, &consdata);
   consdata->nVars = nvars;
   consdata->maxVars = nvars;
   consdata->row = NULL;
   consdata->nFixedNonzero = -1;
   SCIPduplicateBlockMemoryArray(scip, &consdata->Vars, vars, nvars);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
			     local, FALSE, dynamic, removable, FALSE) );

   return SCIP_OKAY;
}
