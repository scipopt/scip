/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_sos1.c
 * @brief  constraint handler for SOS type 1 constraints
 * @author Marc Pfetsch
 *
 * A specially ordered set of type 1 (SOS1) is a sequence of variables such that at most one
 * variable is nonzero. The special case of two variables arises, for instance, from equilibrium or
 * complementary conditions like \f$x \cdot y = 0\f$. Note that it is in principle allowed that a
 * variables appears twice, but it then can be fixed to 0.
 *
 * This implementation of this constraint handler is based on classical ideas, see e.g.@n
 *  "Special Facilities in General Mathematical Programming System for
 *  Non-Convex Problems Using Ordered Sets of Variables"@n
 *  E. Beale and J. Tomlin, Proc. 5th IFORS Conference, 447-454 (1970)
 *
 *
 * The order of the variables is determined as follows:
 *
 * - If the constraint is created with SCIPcreateConsSOS1() and weights are given, the weights
 *   determine the order (decreasing weights). Additional variables can be added with
 *   SCIPaddVarSOS1(), which adds a variable with given weight.
 *
 * - If an empty constraint is created and then variables are added with SCIPaddVarSOS1(), weights
 *   are needed and stored.
 *
 * - All other calls ignore the weights, i.e., if a nonempty constraint is created or variables are
 *   added with SCIPappendVarSOS1().
 *
 * The validity of the constraint is enforced by the classical SOS branching. Depending on the
 * parameters there are two ways to choose the branching constraint. Either the constraint with the
 * most number of nonzeros is chosen or the constraint with the largest nonzero-variable
 * weight. The later version allows the user to specify an order for the branching importance of the
 * constraints. Constraint branching can also be turned off.
 *
 * @todo Possibly allow to generate local cuts via strengthened local cuts (would affect lhs/rhs of rows)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_sos1.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/pub_misc.h"
#include <string.h>
#include <ctype.h>


/* constraint handler properties */
#define CONSHDLR_NAME          "SOS1"
#define CONSHDLR_DESC          "SOS1 constraint handler"
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

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP

/* event handler properties */
#define EVENTHDLR_NAME         "SOS1"
#define EVENTHDLR_DESC         "bound change event handler for SOS1 constraints"


/** constraint data for SOS1 constraints */
struct SCIP_ConsData
{
   int                   nvars;              /**< number of variables in the constraint */
   int                   maxvars;            /**< maximal number of variables (= size of storage) */
   int                   nfixednonzeros;     /**< number of variables fixed to be nonzero */
   SCIP_VAR**            vars;               /**< variables in constraint */
   SCIP_ROW*             row;                /**< row corresponding to upper and lower bound inequalities, or NULL if not yet created */
   SCIP_Real*            weights;            /**< weights determining the order (ascending), or NULL if not used */
};

/** SOS1 constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             branchsos;          /**< Branch on SOS condition in enforcing? */
   SCIP_Bool             branchnonzeros;     /**< Branch on SOS cons. with most number of nonzeros? */
   SCIP_Bool             branchweight;       /**< Branch on SOS cons. with highest nonzero-variable weight for branching - needs branchnonzeros to be false */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
};


/** fix variable in given node to 0 or add constraint if variable is multi-aggregated */
static
SCIP_RETCODE fixVariableZeroNode(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable to be fixed to 0*/
   SCIP_NODE*            node,               /**< node */
   SCIP_Bool*            infeasible          /**< if fixing is infeasible */
   )
{
   /* if variable cannot be nonzero */
   *infeasible = FALSE;
   if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* if variable is multi-aggregated */
   if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_CONS* cons;
      SCIP_Real val;

      val = 1.0;

      if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
      {
         SCIPdebugMessage("creating constraint to force multi-aggregated variable <%s> to 0.\n", SCIPvarGetName(var));
         /* we have to insert a local constraint var = 0 */
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &var, &val, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
               TRUE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }
   else
   {
      if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) )
         SCIP_CALL( SCIPchgVarLbNode(scip, node, var, 0.0) );
      if ( ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
         SCIP_CALL( SCIPchgVarUbNode(scip, node, var, 0.0) );
   }

   return SCIP_OKAY;
}


/** fix variable in local node to 0, and return whether the operation was feasible
 *
 *  @note We do not add a linear constraint if the variable is multi-aggregated as in
 *  fixVariableZeroNode(), since this would be too time consuming.
 */
static
SCIP_RETCODE inferVariableZero(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable to be fixed to 0*/
   SCIP_CONS*            cons,               /**< constraint */
   int                   inferinfo,          /**< info for reverse prop. */
   SCIP_Bool*            infeasible,         /**< if fixing is infeasible */
   SCIP_Bool*            tightened,          /**< if fixing was performed */
   SCIP_Bool*            success             /**< whether fixing was successful, i.e., variable is not multi-aggregated */
   )
{
   *infeasible = FALSE;
   *tightened = FALSE;
   *success = FALSE;

   /* if variable cannot be nonzero */
   if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* directly fix variable if it is not multi-aggregated */
   if ( SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_Bool tighten;

      /* fix lower bound */
      SCIP_CALL( SCIPinferVarLbCons(scip, var, 0.0, cons, inferinfo, FALSE, infeasible, &tighten) );
      *tightened = *tightened || tighten;

      /* fix upper bound */
      SCIP_CALL( SCIPinferVarUbCons(scip, var, 0.0, cons, inferinfo, FALSE, infeasible, &tighten) );
      *tightened = *tightened || tighten;

      *success = TRUE;
   }

   return SCIP_OKAY;
}


/** add lock on variable */
static
SCIP_RETCODE lockVariableSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( var != NULL );

   /* rounding down == bad if lb < 0, rounding up == bad if ub > 0 */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)), SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var))) );

   return SCIP_OKAY;
}


/* remove lock on variable */
static
SCIP_RETCODE unlockVariableSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( var != NULL );

   /* rounding down == bad if lb < 0, rounding up == bad if ub > 0 */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)), SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var))) );

   return SCIP_OKAY;
}


/** ensures that the vars and weights array can store at least num entries */
static
SCIP_RETCODE consdataEnsurevarsSizeSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   num,                /**< minimum number of entries to store */
   SCIP_Bool             reserveWeights      /**< whether the weights array is handled */
   )
{
   assert( consdata != NULL );
   assert( consdata->nvars <= consdata->maxvars );

   if ( num > consdata->maxvars )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->maxvars, newsize) );
      if ( reserveWeights )
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->weights, consdata->maxvars, newsize) );
      consdata->maxvars = newsize;
   }
   assert( num <= consdata->maxvars );

   return SCIP_OKAY;
}


/** handle new variable */
static
SCIP_RETCODE handleNewVariableSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Bool             transformed         /**< whether original variable was transformed */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( var != NULL );

   /* if we are in transformed problem, catch the variable's events */
   if ( transformed )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );
      assert( conshdlrdata->eventhdlr != NULL );

      /* catch bound change events of variable */
      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
            (SCIP_EVENTDATA*)consdata, NULL) );

      /* if the variable if fixed to nonzero */
      assert( consdata->nfixednonzeros >= 0 );
      if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
         ++consdata->nfixednonzeros;
   }

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockVariableSOS1(scip, cons, var) );

   /* add the new coefficient to the LP row, if necessary */
   if ( consdata->row != NULL )
   {
      /* this is currently dead code, since the constraint is not modifiable */
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, 1.0) );

      /* update lhs and rhs if necessary */
      if ( SCIPisFeasGT(scip, SCIPvarGetUbLocal(var), SCIProwGetRhs(consdata->row)) )
         SCIP_CALL( SCIPchgRowRhs(scip, consdata->row, SCIPvarGetUbLocal(var) ) );
      if ( SCIPisFeasLT(scip, SCIPvarGetLbLocal(var), SCIProwGetLhs(consdata->row)) )
         SCIP_CALL( SCIPchgRowLhs(scip, consdata->row, SCIPvarGetLbLocal(var) ) );
   }

   return SCIP_OKAY;
}


/** adds a variable to an SOS1 constraint, at position given by weight - ascending order */
static
SCIP_RETCODE addVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_Real             weight              /**< weight to determine position */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;
   int pos;
   int j;

   assert( var != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( consdata->weights == NULL && consdata->maxvars > 0 )
   {
      SCIPerrorMessage("cannot add variable to SOS1 constraint <%s> that does not contain weights.\n", SCIPconsGetName(cons));
      return SCIP_INVALIDCALL;
   }

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if ( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert( var != NULL );
   assert( transformed == SCIPvarIsTransformed(var) );

   SCIP_CALL( consdataEnsurevarsSizeSOS1(scip, consdata, consdata->nvars + 1, TRUE) );
   assert( consdata->weights != NULL );
   assert( consdata->maxvars >= consdata->nvars+1 );

   /* find variable position */
   for (pos = 0; pos < consdata->nvars; ++pos)
   {
      if ( consdata->weights[pos] > weight )
         break;
   }
   assert( 0 <= pos && pos <= consdata->nvars );

   /* move other variables, if necessary */
   for (j = consdata->nvars; j > pos; --j)
   {
      consdata->vars[j] = consdata->vars[j-1];
      consdata->weights[j] = consdata->weights[j-1];
   }

   /* insert variable */
   consdata->vars[pos] = var;
   consdata->weights[pos] = weight;
   ++consdata->nvars;

   /* handle the new variable */
   SCIP_CALL( handleNewVariableSOS1(scip, cons, consdata, var, transformed) );

   return SCIP_OKAY;
}


/** appends a variable to an SOS1 constraint */
static
SCIP_RETCODE appendVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert( var != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if ( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert( var != NULL );
   assert( transformed == SCIPvarIsTransformed(var) );

   SCIP_CALL( consdataEnsurevarsSizeSOS1(scip, consdata, consdata->nvars + 1, FALSE) );

   /* insert variable */
   consdata->vars[consdata->nvars] = var;
   assert( consdata->weights != NULL || consdata->nvars > 0 );
   if ( consdata->weights != NULL && consdata->nvars > 0 )
      consdata->weights[consdata->nvars] = consdata->weights[consdata->nvars-1] + 1.0;
   ++consdata->nvars;

   /* handle the new variable */
   SCIP_CALL( handleNewVariableSOS1(scip, cons, consdata, var, transformed) );

   return SCIP_OKAY;
}


/** deletes a variable of an SOS1 constraint */
static
SCIP_RETCODE deleteVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< corresponding event handler */
   int                   pos                 /**< position of variable in array */
   )
{
   int j;

   assert( 0 <= pos && pos < consdata->nvars );

   /* remove lock of variable */
   SCIP_CALL( unlockVariableSOS1(scip, cons, consdata->vars[pos]) );

   /* drop events on variable */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );

   /* delete variable - need to copy since order is important */
   for (j = pos; j < consdata->nvars-1; ++j)
   {
      consdata->vars[j] = consdata->vars[j+1];
      if ( consdata->weights != NULL )
         consdata->weights[j] = consdata->weights[j+1];
   }
   --consdata->nvars;

   return SCIP_OKAY;
}


/** perform one presolving round
 *
 *  We perform the following presolving steps.
 *
 *  - If the bounds of some variable force it to be nonzero, we can
 *    fix all other variables to zero and remove the SOS1 constraints
 *    that contain it.
 *  - If a variable is fixed to zero, we can remove the variable.
 *  - If a variable appears twice, it can be fixed to 0.
 *  - We substitute appregated variables.
 */
static
SCIP_RETCODE presolRoundSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_Bool*            cutoff,             /**< whether a cutoff happened */
   SCIP_Bool*            success,            /**< whether we performed a successful reduction */
   int*                  ndelconss,          /**< number of deleted constraints */
   int*                  nupgdconss,         /**< number of upgraded constraints */
   int*                  nfixedvars,         /**< number of fixed variables */
   int*                  nremovedvars        /**< number of variables removed */
   )
{
   SCIP_VAR** vars;
   SCIP_Bool allvarsbinary;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;
   int nfixednonzeros;
   int lastFixedNonzero;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( eventhdlr != NULL );
   assert( cutoff != NULL );
   assert( success != NULL );
   assert( ndelconss != NULL );
   assert( nfixedvars != NULL );
   assert( nremovedvars != NULL );

   *cutoff = FALSE;
   *success = FALSE;

   SCIPdebugMessage("Presolving SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

   j = 0;
   nfixednonzeros = 0;
   lastFixedNonzero = -1;
   allvarsbinary = TRUE;
   vars = consdata->vars;

   /* check for variables fixed to 0 and bounds that fix a variable to be nonzero */
   while ( j < consdata->nvars )
   {
      int l;
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real scalar;
      SCIP_Real constant;

      scalar = 1.0;
      constant = 0.0;

      /* check for aggregation: if the constant is zero the variable is zero iff the aggregated
       * variable is 0 */
      var = vars[j];
      SCIP_CALL( SCIPgetProbvarSum(scip, &var, &scalar, &constant) );

      /* if constant is zero and we get a different variable, substitute variable */
      if ( SCIPisZero(scip, constant) && ! SCIPisZero(scip, scalar) && var != vars[j] )
      {
         SCIPdebugMessage("substituted variable <%s> by <%s>.\n", SCIPvarGetName(vars[j]), SCIPvarGetName(var));
         SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)consdata, -1) );
         SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

         /* change the rounding locks */
         SCIP_CALL( unlockVariableSOS1(scip, cons, consdata->vars[j]) );
         SCIP_CALL( lockVariableSOS1(scip, cons, var) );

         vars[j] = var;
      }

      /* check whether the variable appears again later */
      for (l = j+1; l < consdata->nvars; ++l)
      {
         /* if variable appeared before, we can fix it to 0 and remove it */
         if ( vars[j] == vars[l] )
         {
            SCIPdebugMessage("variable <%s> appears twice in constraint, fixing it to 0.\n", SCIPvarGetName(vars[j]));
            SCIP_CALL( SCIPfixVar(scip, vars[j], 0.0, &infeasible, &fixed) );

            if ( infeasible )
            {
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if ( fixed )
               ++(*nfixedvars);
         }
      }

      /* get bounds */
      lb = SCIPvarGetLbLocal(vars[j]);
      ub = SCIPvarGetUbLocal(vars[j]);

      /* if the variable if fixed to nonzero */
      if ( SCIPisFeasPositive(scip, lb) || SCIPisFeasNegative(scip, ub) )
      {
         ++nfixednonzeros;
         lastFixedNonzero = j;
      }

      /* if the variable is fixed to 0 */
      if ( SCIPisFeasZero(scip, lb) && SCIPisFeasZero(scip, ub) )
      {
         SCIPdebugMessage("deleting variable <%s> fixed to 0.\n", SCIPvarGetName(vars[j]));
         SCIP_CALL( deleteVarSOS1(scip, cons, consdata, eventhdlr, j) );
         ++(*nremovedvars);
      }
      else
      {
         /* check whether all variables are binary */
         if ( ! SCIPvarIsBinary(vars[j]) )
            allvarsbinary = FALSE;

         ++j;
      }
   }

   /* if the number of variables is less than 2 */
   if ( consdata->nvars < 2 )
   {
      SCIPdebugMessage("Deleting SOS1 constraint <%s> with < 2 variables.\n", SCIPconsGetName(cons));

      /* delete constraint */
      assert( ! SCIPconsIsModifiable(cons) );
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
      *success = TRUE;
      return SCIP_OKAY;
   }

   /* if more than one variable are fixed to be nonzero, we are infeasible */
   if ( nfixednonzeros > 1 )
   {
      SCIPdebugMessage("The problem is infeasible: more than one variable has bounds that keep it from being 0.\n");
      assert( lastFixedNonzero >= 0 );
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if there is exactly one fixed nonzero variable */
   if ( nfixednonzeros == 1 )
   {
      assert( lastFixedNonzero >= 0 );

      /* fix all other variables to zero */
      for (j = 0; j < consdata->nvars; ++j)
      {
         if ( j != lastFixedNonzero )
         {
            SCIP_CALL( SCIPfixVar(scip, vars[j], 0.0, &infeasible, &fixed) );
            assert( ! infeasible );
            if ( fixed )
               ++(*nfixedvars);
         }
      }

      SCIPdebugMessage("Deleting redundant SOS1 constraint <%s> with one variable.\n", SCIPconsGetName(cons));

      /* delete original constraint */
      assert( ! SCIPconsIsModifiable(cons) );
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
      *success = TRUE;
   }
   /* note: there is no need to update consdata->nfixednonzeros, since the constraint is deleted as soon nfixednonzeros > 0. */
   else
   {
      /* if all variables are binary create a set packing constraint */
      if ( allvarsbinary )
      {
         SCIP_CONS* setpackcons;

         /* create, add, and release the logicor constraint */
         SCIP_CALL( SCIPcreateConsSetpack(scip, &setpackcons, SCIPconsGetName(cons), consdata->nvars, consdata->vars,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), 
               SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
         SCIP_CALL( SCIPaddCons(scip, setpackcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &setpackcons) );

         SCIPdebugMessage("Upgrading SOS1 constraint <%s> to set packing constraint.\n", SCIPconsGetName(cons));

         /* remove the SOS1 constraint globally */
         assert( ! SCIPconsIsModifiable(cons) );
         SCIP_CALL( SCIPdelCons(scip, cons) );
         ++(*nupgdconss);
         *success = TRUE;
      }
   }

   return SCIP_OKAY;
}


/** propagate variables */
static
SCIP_RETCODE propSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool*            cutoff,             /**< whether a cutoff happened */
   int*                  nGen                /**< number of domain changes */
   )
{
   assert( scip != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( cutoff != NULL );
   assert( nGen != NULL );

   *cutoff = FALSE;

   /* if more than one variable is fixed to be nonzero */
   if ( consdata->nfixednonzeros > 1 )
   {
      SCIPdebugMessage("the node is infeasible, more than 1 variable is fixed to be nonzero.\n");
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if exactly one variable is fixed to be nonzero */
   if ( consdata->nfixednonzeros == 1 )
   {
      SCIP_VAR** vars;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;
      SCIP_Bool success;
      SCIP_Bool allVarFixed;
      int firstFixedNonzero;
      int nvars;
      int j;

      firstFixedNonzero = -1;
      nvars = consdata->nvars;
      vars = consdata->vars;
      assert( vars != NULL );

      /* search nonzero variable - is needed for propinfo */
      for (j = 0; j < nvars; ++j)
      {
         if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(vars[j])) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(vars[j])) )
         {
            firstFixedNonzero = j;
            break;
         }
      }
      assert( firstFixedNonzero >= 0 );

      SCIPdebugMessage("variable <%s> is fixed nonzero, fixing other variables to 0.\n", SCIPvarGetName(vars[firstFixedNonzero]));

      /* fix variables before firstFixedNonzero to 0 */
      allVarFixed = TRUE;
      for (j = 0; j < firstFixedNonzero; ++j)
      {
         /* fix variable */
         SCIP_CALL( inferVariableZero(scip, vars[j], cons, firstFixedNonzero, &infeasible, &tightened, &success) );
         assert( ! infeasible );
         allVarFixed = allVarFixed && success;
         if ( tightened )
            ++(*nGen);
      }

      /* fix variables after firstFixedNonzero to 0 */
      for (j = firstFixedNonzero+1; j < nvars; ++j)
      {
         /* fix variable */
         SCIP_CALL( inferVariableZero(scip, vars[j], cons, firstFixedNonzero, &infeasible, &tightened, &success) );
         assert( ! infeasible ); /* there should be no variables after firstFixedNonzero that are fixed to be nonzero */
         allVarFixed = allVarFixed && success;
         if ( tightened )
            ++(*nGen);
      }

      /* reset constraint age counter */
      if ( *nGen > 0 )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }

      /* delete constraint locally */
      if ( allVarFixed )
      {
         assert( !SCIPconsIsModifiable(cons) );
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      }
   }

   return SCIP_OKAY;
}


/** enforcement method
 *
 *  We check whether the current solution is feasible, i.e., contains at most one nonzero
 *  variable. If not, we branch along the lines indicated by Beale and Tomlin:
 *
 *  We first compute \f$W = \sum_{j=1}^n |x_i|\f$ and \f$w = \sum_{j=1}^n j\, |x_i|\f$. Then we
 *  search for the index \f$k\f$ that satisfies
 *  \f[
 *        k \leq \frac{w}{W} < k+1.
 *  \f]
 *  The branches are then
 *  \f[
 *        x_1 = 0, \ldots, x_k = 0 \qquad \mbox{and}\qquad x_{k+1} = 0, \ldots, x_n = 0.
 *  \f]
 *
 *  If the constraint contains two variables, the branching of course simplifies.
 *
 *  Depending on the parameters (@c branchnonzeros, @c branchweight) there are three ways to choose
 *  the branching constraint.
 *
 *  <TABLE>
 *  <TR><TD>@c branchnonzeros</TD><TD>@c branchweight</TD><TD>constraint chosen</TD></TR>
 *  <TR><TD>@c true          </TD><TD> ?             </TD><TD>most number of nonzeros</TD></TR>
 *  <TR><TD>@c false         </TD><TD> @c true       </TD><TD>maximal weight corresponding to nonzero variable</TD></TR>
 *  <TR><TD>@c false         </TD><TD> @c true       </TD><TD>largest sum of variable values</TD></TR>
 *  </TABLE>
 *
 *  @c branchnonzeros = @c false, @c branchweight = @c true allows the user to specify an order for
 *  the branching importance of the constraints (setting the weights accordingly).
 *
 *  Constraint branching can also be turned off using parameter @c branchsos.
 */
static
SCIP_RETCODE enforceSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           conss,              /**< indicator constraints */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_NODE* node1;
   SCIP_NODE* node2;
   SCIP_CONS* branchCons;
   SCIP_Real maxWeight;
   SCIP_VAR** vars;
   int nvars;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   maxWeight = -SCIP_REAL_MAX;
   branchCons = NULL;

   SCIPdebugMessage("Enforcing SOS1 constraints <%s>.\n", SCIPconshdlrGetName(conshdlr) );
   *result = SCIP_FEASIBLE;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* check each constraint */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONS* cons;
      SCIP_Bool cutoff;
      SCIP_Real weight;
      int nGen;
      int cnt;
      int j;

      cons = conss[c];
      assert( cons != NULL );
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      nGen = 0;
      cnt = 0;
      nvars = consdata->nvars;
      vars = consdata->vars;

      /* do nothing if there are not enough variables - this is usually eliminated by preprocessing */
      if ( nvars < 2 )
         continue;

      /* first perform propagation (it might happen that standard propagation is turned off) */
      SCIP_CALL( propSOS1(scip, cons, consdata, &cutoff, &nGen) );
      SCIPdebugMessage("propagating <%s> in enforcing (cutoff: %u, domain reductions: %d).\n", SCIPconsGetName(cons), cutoff, nGen);
      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if ( nGen > 0 )
      {
         *result = SCIP_REDUCEDDOM;
         return SCIP_OKAY;
      }
      assert( nGen == 0 );

      /* check constraint */
      weight = 0.0;
      for (j = 0; j < nvars; ++j)
      {
         SCIP_Real val = REALABS(SCIPgetSolVal(scip, NULL, vars[j]));

         if ( ! SCIPisFeasZero(scip, val) )
         {
            if ( conshdlrdata->branchnonzeros )
               weight += 1.0;
            else
            {
               if ( conshdlrdata->branchweight )
               {
                  /* choose maximum nonzero-variable weight */
                  if ( consdata->weights[j] > weight )
                     weight = consdata->weights[j];
               }
               else
                  weight += val;
            }
            ++cnt;
         }
      }
      /* if constraint is violated */
      if ( cnt > 1 && weight > maxWeight )
      {
         maxWeight = weight;
         branchCons = cons;
      }
   }

   /* if all constraints are feasible */
   if ( branchCons == NULL )
   {
      SCIPdebugMessage("All SOS1 constraints are feasible.\n");
      return SCIP_OKAY;
   }

   /* if we should leave branching decision to branching rules */
   if ( ! conshdlrdata->branchsos )
   {
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   /* otherwise create branches */
   SCIPdebugMessage("Branching on constraint <%s> (weight: %f).\n", SCIPconsGetName(branchCons), maxWeight);
   consdata = SCIPconsGetData(branchCons);
   assert( consdata != NULL );
   nvars = consdata->nvars;
   vars = consdata->vars;

   if ( nvars == 2 )
   {
      SCIP_Bool infeasible;

      /* constraint is infeasible: */
      assert( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, vars[0])) && ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, vars[1])) );

      /* create branches */
      SCIPdebugMessage("Creating two branches.\n");

      SCIP_CALL( SCIPcreateChild(scip, &node1, SCIPcalcNodeselPriority(scip, vars[0], SCIP_BRANCHDIR_DOWNWARDS, 0.0), SCIPcalcChildEstimate(scip, vars[0], 0.0) ) );
      SCIP_CALL( fixVariableZeroNode(scip, vars[0], node1, &infeasible) );
      assert( ! infeasible );

      SCIP_CALL( SCIPcreateChild(scip, &node2, SCIPcalcNodeselPriority(scip, vars[1], SCIP_BRANCHDIR_DOWNWARDS, 0.0), SCIPcalcChildEstimate(scip, vars[1], 0.0) ) );
      SCIP_CALL( fixVariableZeroNode(scip, vars[1], node2, &infeasible) );
      assert( ! infeasible );
   }
   else
   {
      SCIP_Bool infeasible;
      SCIP_Real weight1;
      SCIP_Real weight2;
      SCIP_Real nodeselest;
      SCIP_Real objest;
      SCIP_Real w;
      int j;
      int ind;
      int cnt;

      cnt = 0;

      weight1 = 0.0;
      weight2 = 0.0;

      /* compute weight */
      for (j = 0; j < nvars; ++j)
      {
         SCIP_Real val = REALABS(SCIPgetSolVal(scip, NULL, vars[j]));
         weight1 += val * (SCIP_Real) j;
         weight2 += val;

         if ( ! SCIPisFeasZero(scip, val) )
            ++cnt;
      }

      assert( cnt >= 2 );
      assert( !SCIPisFeasZero(scip, weight2) );
      w = weight1/weight2;  /*lint !e795*/

      ind = (int) SCIPfloor(scip, w);
      assert( 0 <= ind && ind < nvars-1 );

      /* branch on variable ind: either all variables up to ind or all variables after ind are zero */
      SCIPdebugMessage("Branching on variable <%s>.\n", SCIPvarGetName(vars[ind]));

      /* calculate node selection and objective estimate for node 1 */
      nodeselest = 0.0;
      objest = 0.0;
      for (j = 0; j <= ind; ++j)
      {
         nodeselest += SCIPcalcNodeselPriority(scip, vars[j], SCIP_BRANCHDIR_DOWNWARDS, 0.0);
         objest += SCIPcalcChildEstimate(scip, vars[j], 0.0);
      }
      /* take the average of the individual estimates */
      objest = objest/((SCIP_Real) ind + 1.0);

      /* create node 1 */
      SCIP_CALL( SCIPcreateChild(scip, &node1, nodeselest, objest) );
      for (j = 0; j <= ind; ++j)
      {
         SCIP_CALL( fixVariableZeroNode(scip, vars[j], node1, &infeasible) );
         assert( ! infeasible );
      }

      /* calculate node selection and objective estimate for node 1 */
      nodeselest = 0.0;
      objest = 0.0;
      for (j = ind+1; j < nvars; ++j)
      {
         nodeselest += SCIPcalcNodeselPriority(scip, vars[j], SCIP_BRANCHDIR_DOWNWARDS, 0.0);
         objest += SCIPcalcChildEstimate(scip, vars[j], 0.0);
      }
      /* take the average of the individual estimates */
      objest = objest/((SCIP_Real) (nvars - ind - 1));

      /* create node 2 */
      SCIP_CALL( SCIPcreateChild(scip, &node2, nodeselest, objest) );
      for (j = ind+1; j < nvars; ++j)
      {
         SCIP_CALL( fixVariableZeroNode(scip, vars[j], node2, &infeasible) );
         assert( ! infeasible );
      }
   }
   SCIP_CALL( SCIPresetConsAge(scip, branchCons) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** Generate row
 *
 *  We generate the row corresponding to the following simple valid inequalities:
 *  \f[
 *         x_1 + \ldots + x_n \leq \max\{u_1, \ldots, u_n\}\qquad\mbox{and}\qquad
 *         x_1 + \ldots + x_n \geq \min\{l_1, \ldots, l_n\},
 *  \f]
 *  where \f$l_1, \ldots, l_n\f$ and \f$u_1, \ldots, u_n\f$ are the
 *  lower and upper bounds of the variables \f$x_1, \ldots, x_n\f$.
 *  Of course, these inequalities are only added if the upper and
 *  lower bounds are all finite and at least one lower bound is < 0 or
 *  one upper bound is > 0 (lower bounds > 0 and upper bounds < 0 are
 *  usually detected in presolving).
 */
static
SCIP_RETCODE generateRowSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   SCIP_VAR** vars;
   SCIP_Real minLb;
   SCIP_Real maxUb;
   SCIP_ROW* row;
   int nvars;
   int j;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->row == NULL );

   minLb = SCIPinfinity(scip);
   maxUb = -SCIPinfinity(scip);
   nvars = consdata->nvars;
   vars = consdata->vars;
   assert( vars != NULL );

   /* find minimum and maximum lower and upper bounds */
   for (j = 0; j < nvars; ++j)
   {
      if ( SCIPvarGetLbLocal(vars[j]) < minLb )
         minLb = SCIPvarGetLbLocal(vars[j]);
      if ( SCIPvarGetUbLocal(vars[j]) > maxUb )
         maxUb = SCIPvarGetUbLocal(vars[j]);
   }

   /* ignore trivial inequality if all lower bounds are 0 */
   if ( SCIPisFeasZero(scip, minLb) )
      minLb = -SCIPinfinity(scip);

   /* ignore trivial inequality if all upper bounds are 0 */
   if ( SCIPisFeasZero(scip, maxUb) )
      maxUb = SCIPinfinity(scip);

   /* create upper and lower bound inequality if one of the bounds is finite */
   if ( ! SCIPisInfinity(scip, REALABS(minLb)) || ! SCIPisInfinity(scip, REALABS(maxUb)) )
   {
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "sosbnd", minLb, maxUb, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, row, nvars, vars, 1.0) );
      consdata->row = row;
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
   }

   return SCIP_OKAY;
}

/* ---------------------------- constraint handler callback methods ----------------------*/

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSOS1(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


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

   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolSOS1)
{  /*lint --e{715}*/
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
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
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
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMessage("Deleting SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

   /* drop events on transformed variables */
   if ( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int j;

      /* get constraint handler data */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );
      assert( conshdlrdata->eventhdlr != NULL );

      for (j = 0; j < (*consdata)->nvars; ++j)
      {
         SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
               (SCIP_EVENTDATA*)*consdata, -1) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->maxvars);
   if ( (*consdata)->weights != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->weights, (*consdata)->maxvars);
   }

   /* free row */
   if ( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }
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
   assert( sourcedata->nvars > 0 );
   assert( sourcedata->nvars <= sourcedata->maxvars );

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nvars = sourcedata->nvars;
   consdata->maxvars = sourcedata->nvars;
   consdata->row = NULL;
   consdata->nfixednonzeros = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, consdata->nvars) );
   /* if weights were used */
   if ( sourcedata->weights != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->weights, sourcedata->weights, consdata->nvars) );
   }
   else
      consdata->weights = NULL;

   for (j = 0; j < sourcedata->nvars; ++j)
   {
      assert( sourcedata->vars[j] != 0 );
      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->vars[j], &(consdata->vars[j])) );

      /* if variable is fixed to be nonzero */
      if ( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(consdata->vars[j])) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(consdata->vars[j])) )
         ++(consdata->nfixednonzeros);
   }

   /* create transformed constraint with the same flags */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* catch bound change events on variable */
   for (j = 0; j < consdata->nvars; ++j)
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[j], SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr,
            (SCIP_EVENTDATA*)consdata, NULL) );
   }

#ifdef SCIP_DEBUG
   if ( consdata->nfixednonzeros > 0 )
   {
      SCIPdebugMessage("constraint <%s> has %d variables fixed to be nonzero.\n", SCIPconsGetName(*targetcons),
         consdata->nfixednonzeros );
   }
#endif

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSOS1)
{  /*lint --e{715}*/
   int oldnfixedvars;
   int oldndelconss;
   int oldnupgdconss;
   int nremovedvars;
   SCIP_EVENTHDLR* eventhdlr;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMessage("Presolving SOS1 constraints.\n");

   *result = SCIP_DIDNOTRUN;
   oldnfixedvars = *nfixedvars;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;
   nremovedvars = 0;

   /* only run if success if possible */
   if( nrounds == 0 || nnewfixedvars > 0 || nnewaggrvars > 0 )
   {
      /* get constraint handler data */
      assert( SCIPconshdlrGetData(conshdlr) != NULL );
      eventhdlr = SCIPconshdlrGetData(conshdlr)->eventhdlr;
      assert( eventhdlr != NULL );

      *result = SCIP_DIDNOTFIND;

      /* check each constraint */
      for (c = 0; c < nconss; ++c)
      {
         SCIP_CONSDATA* consdata;
         SCIP_CONS* cons;
         SCIP_Bool cutoff;
         SCIP_Bool success;

         assert( conss != NULL );
         assert( conss[c] != NULL );
         cons = conss[c];
         consdata = SCIPconsGetData(cons);

         assert( consdata != NULL );
         assert( consdata->nvars >= 0 );
         assert( consdata->nvars <= consdata->maxvars );
         assert( ! SCIPconsIsModifiable(cons) );

         /* perform one presolving round */
         SCIP_CALL( presolRoundSOS1(scip, cons, consdata, eventhdlr, &cutoff, &success, ndelconss, nupgdconss, nfixedvars, &nremovedvars) );

         if ( cutoff )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         if ( success )
            *result = SCIP_SUCCESS;
      }
   }
   (*nchgcoefs) += nremovedvars;

   SCIPdebugMessage("presolving fixed %d variables, removed %d variables, deleted %d constraints, and upgraded %d constraints.\n",
      *nfixedvars - oldnfixedvars, nremovedvars, *ndelconss - oldndelconss, *nupgdconss - oldnupgdconss);

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
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

      /* possibly generate row if not yet done */
      if ( consdata->row == NULL )
      {
         SCIP_CALL( generateRowSOS1(scip, conshdlr, consdata) );
      }

      /* put corresponding rows into LP */
      if ( consdata->row != NULL && ! SCIProwIsInLP(consdata->row) )
      {
         assert( ! SCIPisInfinity(scip, REALABS(SCIProwGetLhs(consdata->row))) || ! SCIPisInfinity(scip, REALABS(SCIProwGetRhs(consdata->row))) );

         SCIP_CALL( SCIPaddCut(scip, NULL, consdata->row, FALSE) );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, consdata->row, NULL) ) );
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpSOS1)
{  /*lint --e{715}*/
   int nGen = 0;
   int c;

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

      /* put corresponding rows into LP if they are useful */
      row = consdata->row;

      /* possibly generate row if not yet done */
      if ( row == NULL )
      {
         SCIP_CALL( generateRowSOS1(scip, conshdlr, consdata) );
      }

      /* possibly add row to LP if it is useful */
      if ( row != NULL && ! SCIProwIsInLP(row) && SCIPisCutEfficacious(scip, NULL, row) )
      {
         assert( ! SCIPisInfinity(scip, REALABS(SCIProwGetLhs(consdata->row))) || ! SCIPisInfinity(scip, REALABS(SCIProwGetRhs(consdata->row))) );

         SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE) );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
         SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
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
{  /*lint --e{715}*/
   int nGen = 0;
   int c;

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

      /* put corresponding row into LP if it is useful */
      row = consdata->row;

      /* possibly generate row if not yet done */
      if ( row == NULL )
      {
         SCIP_CALL( generateRowSOS1(scip, conshdlr, consdata) );
      }

      /* possibly add row to LP if it is useful */
      if ( row != NULL && ! SCIProwIsInLP(row) && SCIPisCutEfficacious(scip, NULL, row) )
      {
         assert( ! SCIPisInfinity(scip, REALABS(SCIProwGetLhs(consdata->row))) || ! SCIPisInfinity(scip, REALABS(SCIProwGetRhs(consdata->row))) );

         SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
         SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
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
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( enforceSOS1(scip, conshdlr, nconss, conss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSOS1)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( enforceSOS1(scip, conshdlr, nconss, conss, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions
 *
 *  We simply check whether at most one variable is nonzero in the given solution.
 */
static
SCIP_DECL_CONSCHECK(consCheckSOS1)
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
      int j;
      int cnt;

      cnt = 0;
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      SCIPdebugMessage("Checking SOS1 constraint <%s>.\n", SCIPconsGetName(conss[c]));

      /* check all variables */
      for (j = 0; j < consdata->nvars; ++j)
      {
         /* if variable is nonzero */
         if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[j])) )
         {
            ++cnt;

            /* if more than one variable is nonzero */
            if ( cnt > 1 )
            {
               SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
               *result = SCIP_INFEASIBLE;

               if ( printreason )
               {
                  int l;

                  SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
                  SCIPinfoMessage(scip, NULL, ";\nviolation: ");

                  for (l = 0; l < consdata->nvars; ++l)
                  {
                     /* if variable is nonzero */
                     if ( ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[l])) )
                     {
                        SCIPinfoMessage(scip, NULL, "<%s> = %.15g ",
                           SCIPvarGetName(consdata->vars[l]), SCIPgetSolVal(scip, sol, consdata->vars[l]));
                     }
                  }
                  SCIPinfoMessage(scip, NULL, "\n");
               }
               return SCIP_OKAY;
            }
         }
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSOS1)
{  /*lint --e{715}*/
   int nGen = 0;
   int c;

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
      SCIPdebugMessage("Propagating SOS1 constraint <%s>.\n", SCIPconsGetName(cons) );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( propSOS1(scip, cons, consdata, &cutoff, &nGen) );
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
 *  use that @a inferinfo stores the index of the variable that has
 *  bounds that fix it to be nonzero (these bounds are the reason). */
static
SCIP_DECL_CONSRESPROP(consRespropSOS1)
{  /*lint --e{715}*/
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
   assert( 0 <= inferinfo && inferinfo < consdata->nvars );
   var = consdata->vars[inferinfo];
   assert( var != infervar );

   /* check if lower bound of var was the reason */
   if ( SCIPisFeasPositive(scip, SCIPvarGetLbAtIndex(var, bdchgidx, FALSE)) )
   {
      SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   /* check if upper bound of var was the reason */
   if ( SCIPisFeasNegative(scip, SCIPvarGetUbAtIndex(var, bdchgidx, FALSE)) )
   {
      SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler
 *
 *  Let lb and ub be the lower and upper bounds of a
 *  variable. Preprocessing usually makes sure that lb <= 0 <= ub.
 *
 *  - If lb < 0 then rounding down may violate the constraint.
 *  - If ub > 0 then rounding up may violated the constraint.
 *  - If lb > 0 or ub < 0 then the constraint is infeasible and we do
 *    not have to deal with it here.
 *  - If lb == 0 then rounding down does not violate the constraint.
 *  - If ub == 0 then rounding up does not violate the constraint.
 */
static
SCIP_DECL_CONSLOCK(consLockSOS1)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIPdebugMessage("Locking constraint <%s>.\n", SCIPconsGetName(cons));

   vars = consdata->vars;
   nvars = consdata->nvars;
   assert( vars != NULL );

   for (j = 0; j < nvars; ++j)
   {
      SCIP_VAR* var;
      var = vars[j];

      /* if lower bound is negative, rounding down may violate constraint */
      if ( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)) )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, var, nlockspos, nlocksneg) );
      }

      /* additionally: if upper bound is positive, rounding up may violate constraint */
      if ( SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var)) )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, var, nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintSOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   for (j = 0; j < consdata->nvars; ++j)
   {
      if ( j > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->vars[j], FALSE) );
      if ( consdata->weights == NULL )
         SCIPinfoMessage(scip, file, " (%d)", j+1);
      else
         SCIPinfoMessage(scip, file, " (%3.2f)", consdata->weights[j]);
   }

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourceconsdata;
   SCIP_VAR** sourcevars;
   SCIP_VAR** targetvars;
   SCIP_Real* sourceweights;
   SCIP_Real* targetweights;
   const char* consname;
   int nvars;
   int v;

   assert( scip != NULL );
   assert( sourcescip != NULL );
   assert( sourcecons != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0 );

   *valid = TRUE;

   if ( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   SCIPdebugMessage("Copying SOS1 constraint <%s> ...\n", consname);

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert( sourceconsdata != NULL );

   /* get variables and weights of the source constraint */
   nvars = sourceconsdata->nvars;

   if ( nvars == 0 )
      return SCIP_OKAY;

   sourcevars = sourceconsdata->vars;
   assert( sourcevars != NULL );
   sourceweights = sourceconsdata->weights;
   assert( sourceweights != NULL );

   /* duplicate variable array */
   SCIP_CALL( SCIPallocBufferArray(sourcescip, &targetvars, nvars) );
   SCIP_CALL( SCIPduplicateBufferArray(sourcescip, &targetweights, sourceweights, nvars) );

   /* get copied variables in target SCIP */
   for( v = 0; v < nvars && *valid; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &(targetvars[v]), varmap, consmap, global, valid) );
   }

    /* only create the target constraint, if all variables could be copied */
   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsSOS1(scip, cons, consname, nvars, targetvars, targetweights,
            initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );
   }

   /* free buffer array */
   SCIPfreeBufferArray(sourcescip, &targetweights);
   SCIPfreeBufferArray(sourcescip, &targetvars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseSOS1)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_Real weight;
   const char* s;
   char* t;

   *success = TRUE;
   s = str;

   /* create empty SOS1 constraint */
   SCIP_CALL( SCIPcreateConsSOS1(scip, cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );

   /* loop through string */
   do
   {
      /* parse variable name */
      SCIP_CALL( SCIPparseVarName(scip, s, &var, &t) );
      s = t;

      /* skip until beginning of weight */
      while ( *s != '\0' && *s != '(' )
         ++s;

      if ( *s == '\0' )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: expected weight at input: %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      /* skip '(' */
      ++s;

      /* find weight */
      weight = strtod(s, &t);
      if ( t == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error during parsing of the weight: %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      s = t;

      /* skip white space, ',', and ')' */
      while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' || *s == ')' ) )
         ++s;

      /* add variable */
      SCIP_CALL( SCIPaddVarSOS1(scip, *cons, var, weight) );
   }
   while ( *s != '\0' );

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSOS1)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}


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
   SCIP_Real oldbound;
   SCIP_Real newbound;

   assert( eventhdlr != NULL );
   assert( eventdata != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );
   assert( event != NULL );

   consdata = (SCIP_CONSDATA*)eventdata;
   assert( consdata != NULL );
   assert( 0 <= consdata->nfixednonzeros && consdata->nfixednonzeros <= consdata->nvars );

   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);

   eventtype = SCIPeventGetType(event);
   switch ( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisFeasPositive(scip, oldbound) && SCIPisFeasPositive(scip, newbound) )
         ++(consdata->nfixednonzeros);
      break;
   case SCIP_EVENTTYPE_UBTIGHTENED:
      /* if variable is now fixed to be nonzero */
      if ( ! SCIPisFeasNegative(scip, oldbound) && SCIPisFeasNegative(scip, newbound) )
         ++(consdata->nfixednonzeros);
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisFeasPositive(scip, oldbound) && ! SCIPisFeasPositive(scip, newbound) )
         --(consdata->nfixednonzeros);
      break;
   case SCIP_EVENTTYPE_UBRELAXED:
      /* if variable is not fixed to be nonzero anymore */
      if ( SCIPisFeasNegative(scip, oldbound) && ! SCIPisFeasNegative(scip, newbound) )
         --(consdata->nfixednonzeros);
      break;
   default:
      SCIPerrorMessage("invalid event type.\n");
      return SCIP_INVALIDDATA;
   }
   assert( 0 <= consdata->nfixednonzeros && consdata->nfixednonzeros <= consdata->nvars );

   SCIPdebugMessage("changed bound of variable <%s> from %f to %f (nfixednonzeros: %d).\n", SCIPvarGetName(SCIPeventGetVar(event)),
                    oldbound, newbound, consdata->nfixednonzeros);

   return SCIP_OKAY;
}


/* ---------------- Constraint specific interface methods ---------------- */

/** creates the handler for SOS1 constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSOS1(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->branchsos = TRUE;
   conshdlrdata->eventhdlr = NULL;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecSOS1, NULL) );
   if ( conshdlrdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for SOS1 constraints not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSOS1, consEnfopsSOS1, consCheckSOS1, consLockSOS1, conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySOS1, consCopySOS1) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSOS1) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolSOS1) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSOS1) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSOS1) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSOS1) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSOS1) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSOS1) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSOS1, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSOS1) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSOS1, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSOS1) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSOS1, consSepasolSOS1, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSOS1) );

   /* add SOS1 constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SOS1/branchsos",
         "Use SOS1 branching in enforcing (otherwise leave decision to branching rules)?",
         &conshdlrdata->branchsos, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SOS1/branchnonzeros",
         "Branch on SOS constraint with most number of nonzeros?",
         &conshdlrdata->branchnonzeros, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SOS1/branchweight",
         "Branch on SOS cons. with highest nonzero-variable weight for branching (needs branchnonzeros = false)?",
         &conshdlrdata->branchweight, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}


/** creates and captures a SOS1 constraint
 *
 *  We set the constraint to not be modifable. If the weights are non
 *  NULL, the variables are ordered according to these weights (in
 *  ascending order).
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            weights,            /**< weights determining the variable order, or NULL if natural order should be used */
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
   SCIP_Bool modifiable;

   modifiable = FALSE;

   /* find the SOS1 constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("<%s> constraint handler not found\n", CONSHDLR_NAME);
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->vars = NULL;
   consdata->nvars = nvars;
   consdata->maxvars = nvars;
   consdata->row = NULL;
   consdata->nfixednonzeros = -1;
   consdata->weights = NULL;

   if ( nvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->vars, vars, nvars) );

      /* check weights */
      if ( weights != NULL )
      {
         /* store weights */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->weights, weights, nvars) );

         /* sort variables - ascending order */
         SCIPsortRealPtr(consdata->weights, (void**)consdata->vars, nvars);
      }
   }
   else
   {
      assert( weights == NULL );
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}


/** creates and captures a SOS1 constraint with all constraint flags set to their default values.
 *
 *  @warning Do NOT set the constraint to be modifiable manually, because this might lead
 *  to wrong results as the variable array will not be resorted
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            weights             /**< weights determining the variable order, or NULL if natural order should be used */
   )
{
   SCIP_CALL( SCIPcreateConsSOS1( scip, cons, name, nvars, vars, weights, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}


/** adds variable to SOS1 constraint, the position is determined by the given weight */
SCIP_RETCODE SCIPaddVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_Real             weight              /**< weight determining position of variable */
   )
{
   assert( scip != NULL );
   assert( var != NULL );
   assert( cons != NULL );

   SCIPdebugMessage("adding variable <%s> to constraint <%s> with weight %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), weight);

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( addVarSOS1(scip, cons, var, weight) );

   return SCIP_OKAY;
}


/** appends variable to SOS1 constraint */
SCIP_RETCODE SCIPappendVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{
   assert( scip != NULL );
   assert( var != NULL );
   assert( cons != NULL );

   SCIPdebugMessage("appending variable <%s> to constraint <%s>\n", SCIPvarGetName(var), SCIPconsGetName(cons));

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( appendVarSOS1(scip, cons, var) );

   return SCIP_OKAY;
}


/** gets number of variables in SOS1 constraint */
int SCIPgetNVarsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->nvars;
}


/** gets array of variables in SOS1 constraint */
SCIP_VAR** SCIPgetVarsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->vars;
}


/** gets array of weights in SOS1 constraint (or NULL if not existent) */
SCIP_Real* SCIPgetWeightsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an SOS1 constraint.\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->weights;
}
