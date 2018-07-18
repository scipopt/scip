/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   benders.c
 * @brief  methods for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/pricestore.h"
#include "scip/scip.h"
#include "scip/benders.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/cons_linear.h"

#include "scip/struct_benders.h"
#include "scip/struct_benderscut.h"

#include "scip/benderscut.h"

/* Defaults for parameters */
#define SCIP_DEFAULT_TRANSFERCUTS          TRUE  /** should Benders' cuts generated in LNS heuristics be transferred to the main SCIP instance? */
#define SCIP_DEFAULT_CUTSASCONSS           TRUE  /** should the transferred cuts be added as constraints? */
#define SCIP_DEFAULT_LNSCHECK              TRUE  /** should the Benders' decomposition be used in LNS heuristics */
#define SCIP_DEFAULT_LNSMAXDEPTH             -1  /** maximum depth at which the LNS check is performed */
#define SCIP_DEFAULT_SUBPROBFRAC            1.0  /** fraction of subproblems that are solved in each iteration */
#define SCIP_DEFAULT_UPDATEAUXVARBOUND     TRUE  /** should the auxiliary variable lower bound be updated by solving the subproblem */

#define BENDERS_MAXPSEUDOSOLS                 5  /** the maximum number of pseudo solutions checked before suggesting
                                                     merge candidates */

#define AUXILIARYVAR_NAME     "##bendersauxiliaryvar" /** the name for the Benders' auxiliary variables in the master problem */

/* event handler properties */
#define NODEFOCUS_EVENTHDLR_NAME         "bendersnodefocus"
#define NODEFOCUS_EVENTHDLR_DESC         "node focus event handler for Benders' decomposition"

#define MIPNODEFOCUS_EVENTHDLR_NAME      "bendersmipsolvenodefocus"
#define MIPNODEFOCUS_EVENTHDLR_DESC      "node focus event handler for the MIP solve method for Benders' decomposition"

#define UPPERBOUND_EVENTHDLR_NAME        "bendersupperbound"
#define UPPERBOUND_EVENTHDLR_DESC        "found solution event handler to terminate subproblem solve for a given upper bound"

#define NODESOLVED_EVENTHDLR_NAME        "bendersnodesolved"
#define NODESOLVED_EVENTHDLR_DESC        "node solved event handler for the Benders' integer cuts"


/** event handler data */
struct SCIP_EventhdlrData
{
   int                   filterpos;          /**< the event filter entry */
   int                   numruns;            /**< the number of times that the problem has been solved */
   SCIP_Real             upperbound;         /**< an upper bound for the problem */
   SCIP_Bool             solvecip;           /**< is the event called from a MIP subproblem solve*/
};


/* ---------------- Local methods for event handlers ---------------- */

/** initialises the members of the eventhandler data */
static
SCIP_RETCODE initEventhandlerData(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< the event handler data */
   )
{
   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->filterpos = -1;
   eventhdlrdata->numruns = 0;
   eventhdlrdata->upperbound = -SCIPinfinity(scip);
   eventhdlrdata->solvecip = FALSE;

   return SCIP_OKAY;
}

/** initsol method for the event handlers */
static
SCIP_RETCODE initsolEventhandler(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< the event handlers data structure */
   SCIP_EVENTTYPE        eventtype           /**< event type mask to select events to catch */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   SCIP_CALL(SCIPcatchEvent(scip, eventtype, eventhdlr, NULL, &eventhdlrdata->filterpos));

   return SCIP_OKAY;
}

/** the exit sol method for the event handlers */
static
SCIP_RETCODE exitsolEventhandler(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< the event handlers data structure */
   SCIP_EVENTTYPE        eventtype           /**< event type mask to select events to catch */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   if( eventhdlrdata->filterpos >= 0 )
   {
      SCIP_CALL(SCIPdropEvent(scip, eventtype, eventhdlr, NULL, eventhdlrdata->filterpos));
      eventhdlrdata->filterpos = -1;
   }

   return SCIP_OKAY;
}

/** the exit method for the event handlers */
static
SCIP_RETCODE exitEventhandler(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr           /**< the event handlers data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   /* reinitialise the event handler data */
   SCIP_CALL( initEventhandlerData(scip, eventhdlrdata) );

   return SCIP_OKAY;
}

/** free method for the event handler */
static
SCIP_RETCODE freeEventhandler(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr           /**< the event handlers data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &eventhdlrdata);

   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}



/* ---------------- Callback methods of node focus event handler ---------------- */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBendersNodefocus)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   /* sending an interrupt solve signal to return the control back to the Benders' decomposition plugin.
    * This will ensure the SCIP stage is SCIP_STAGE_SOLVING, allowing the use of probing mode. */
   SCIP_CALL( SCIPinterruptSolve(scip) );

   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, eventhdlrdata->filterpos));
   eventhdlrdata->filterpos = -1;

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersNodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( initsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_NODEFOCUSED) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolBendersNodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_NODEFOCUSED) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitBendersNodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTFREE(eventFreeBendersNodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( freeEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}


/* ---------------- Callback methods of MIP solve node focus event handler ---------------- */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBendersMipnodefocus)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   /* interrupting the solve so that the control is returned back to the Benders' core. */
   if( eventhdlrdata->numruns == 0 && !eventhdlrdata->solvecip )
   {
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, eventhdlrdata->filterpos));
   eventhdlrdata->filterpos = -1;

   eventhdlrdata->numruns++;

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersMipnodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( initsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_NODEFOCUSED) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolBendersMipnodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_NODEFOCUSED) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitBendersMipnodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTFREE(eventFreeBendersMipnodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( freeEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}

/* ---------------- Callback methods of solution found event handler ---------------- */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBendersUpperbound)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_SOL* bestsol;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   bestsol = SCIPgetBestSol(scip);

   if( SCIPisLT(scip, SCIPgetSolOrigObj(scip, bestsol), eventhdlrdata->upperbound) )
   {
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersUpperbound)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   SCIP_CALL( initsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_BESTSOLFOUND) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolBendersUpperbound)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_BESTSOLFOUND) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitBendersUpperbound)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTFREE(eventFreeBendersUpperbound)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   SCIP_CALL( freeEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}

/** updates the upper bound in the event handler data */
static
SCIP_RETCODE updateEventhdlrUpperbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             upperbound          /**< the upper bound value */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   eventhdlr = SCIPfindEventhdlr(SCIPbendersSubproblem(benders, probnumber), UPPERBOUND_EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->upperbound = upperbound;

   return SCIP_OKAY;
}

/* ---------------- Callback methods of the node solved event handler ---------------- */

/** Updates the cut constant of the Benders' cuts data.
 *  This function solves the master problem with only the auxiliary variables in the objective function.
 */
static
SCIP_RETCODE updateSubproblemLowerbound(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int nsubproblems;
   int i;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;

   assert(masterprob != NULL);
   assert(benders != NULL);

   /* don't run in probing or in repropagation */
   if( SCIPinProbing(masterprob) || SCIPinRepropagation(masterprob) || SCIPinDive(masterprob) )
      return SCIP_OKAY;

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   SCIP_CALL( SCIPstartProbing(masterprob) );

   /* change the master problem variables to 0 */
   nvars = SCIPgetNVars(masterprob);
   vars = SCIPgetVars(masterprob);

   /* setting the objective function coefficient to 0 for all variables */
   for( i = 0; i < nvars; i++ )
   {
      if( SCIPvarGetStatus(vars[i]) == SCIP_VARSTATUS_COLUMN )
      {
         SCIP_CALL( SCIPchgVarObjProbing(masterprob, vars[i], 0.0) );
      }
   }

   /* solving an LP for all subproblems to find the lower bound */
   for( i = 0; i < nsubproblems; i++)
   {
      SCIP_VAR* auxiliaryvar;

      auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, i);

      if( SCIPvarGetStatus(auxiliaryvar) != SCIP_VARSTATUS_COLUMN )
         continue;

      SCIP_CALL( SCIPchgVarObjProbing(masterprob, auxiliaryvar, 1.0) );

      /* solving the probing LP to get a lower bound on the auxiliary variables */
      SCIP_CALL( SCIPsolveProbingLP(masterprob, -1, &lperror, &cutoff) );

      if( !SCIPisInfinity(masterprob, -SCIPgetSolTransObj(masterprob, NULL)) )
         SCIPbendersUpdateSubproblemLowerbound(benders, i, SCIPgetSolTransObj(masterprob, NULL));

      SCIPdebugMsg(masterprob, "Cut constant for subproblem %d: %g\n", i,
         SCIPbendersGetSubproblemLowerbound(benders, i));

      SCIP_CALL( SCIPchgVarObjProbing(masterprob, auxiliaryvar, 0.0) );
   }

   SCIP_CALL( SCIPendProbing(masterprob) );

   return SCIP_OKAY;
}

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBendersNodesolved)
{  /*lint --e{715}*/
   SCIP_BENDERS* benders;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODESOLVED_EVENTHDLR_NAME) == 0);

   benders = (SCIP_BENDERS*)SCIPeventhdlrGetData(eventhdlr);   /*lint !e826*/

   if( SCIPbendersGetNSubproblems(benders) > 0
      && SCIPbendersGetNSubproblems(benders) > SCIPbendersGetNConvexSubproblems(benders) )
   {
      SCIP_CALL( updateSubproblemLowerbound(scip, benders) );
   }

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersNodesolved)
{
   SCIP_BENDERS* benders;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODESOLVED_EVENTHDLR_NAME) == 0);

   /* getting the Benders' decomposition data structure */
   benders = (SCIP_BENDERS*)SCIPeventhdlrGetData(eventhdlr);   /*lint !e826*/

   /* The event is only caught if there is an active Benders' decomposition */
   if( SCIPbendersIsActive(benders) && !SCIPbendersOnlyCheckConvexRelax(benders) )
   {
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, NULL) );
   }

   return SCIP_OKAY;
}



/* Local methods */

/** A workaround for GCG. This is a temp vardata that is set for the auxiliary variables */
struct SCIP_VarData
{
   int                   vartype;             /**< the variable type. In GCG this indicates whether the variable is a
                                               *   master problem or subproblem variable. */
};


/** adds the auxiliary variables to the Benders' decomposition master problem */
static
SCIP_RETCODE addAuxiliaryVariablesToMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition structure */
   )
{
   SCIP_BENDERS* topbenders;        /* the highest priority Benders' decomposition */
   SCIP_VAR* auxiliaryvar;
   SCIP_VARDATA* vardata;
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */
   SCIP_Bool shareauxvars;
   int i;

   /* this is a workaround for GCG. GCG expects that the variable has vardata when added. So a dummy vardata is created */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = -1;

   /* getting the highest priority Benders' decomposition */
   topbenders = SCIPgetBenders(scip)[0];

   /* if the current Benders is the highest priority Benders, then we need to create the auxiliary variables.
    * Otherwise, if the shareauxvars flag is set, then the auxiliary variables from the highest priority Benders' are
    * stored with this Benders. */
   shareauxvars = FALSE;
   if( topbenders != benders && SCIPbendersShareAuxVars(benders) )
      shareauxvars = TRUE;

   for( i = 0; i < SCIPbendersGetNSubproblems(benders); i++ )
   {
      /* if the auxiliary variables are shared, then a pointer to the variable is retrieved from topbenders,
       * otherwise the auxiliaryvariable is created. */
      if( shareauxvars )
      {
         auxiliaryvar = SCIPbendersGetAuxiliaryVar(topbenders, i);

         SCIP_CALL( SCIPcaptureVar(scip, auxiliaryvar) );
      }
      else
      {
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s_%d_%s", AUXILIARYVAR_NAME, i, SCIPbendersGetName(benders) );
         SCIP_CALL( SCIPcreateVarBasic(scip, &auxiliaryvar, varname, -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
               SCIP_VARTYPE_CONTINUOUS) );

         SCIPvarSetData(auxiliaryvar, vardata);

         SCIP_CALL( SCIPaddVar(scip, auxiliaryvar) );
      }

      benders->auxiliaryvars[i] = auxiliaryvar;
   }

   SCIPfreeBlockMemory(scip, &vardata);

   return SCIP_OKAY;
}

/** assigns the copied auxiliary variables in the target SCIP to the target Benders' decomposition data */
static
SCIP_RETCODE assignAuxiliaryVariables(
   SCIP*                 scip,               /**< SCIP data structure, the target scip */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_BENDERS* topbenders;        /* the highest priority Benders' decomposition */
   SCIP_VAR* targetvar;
   SCIP_VARDATA* vardata;
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */
   SCIP_Bool shareauxvars;
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   /* this is a workaround for GCG. GCG expects that the variable has vardata when added. So a dummy vardata is created */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = -1;

   /* getting the highest priority Benders' decomposition */
   topbenders = SCIPgetBenders(scip)[0];

   /* if the auxiliary variable are shared, then the variable name will have a suffix of the highest priority Benders'
    * name. So the shareauxvars flag indicates how to search for the auxiliary variables */
   shareauxvars = FALSE;
   if( topbenders != benders && SCIPbendersShareAuxVars(benders) )
      shareauxvars = TRUE;

   for( i = 0; i < SCIPbendersGetNSubproblems(benders); i++ )
   {
      if( shareauxvars )
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s_%d_%s", AUXILIARYVAR_NAME, i, SCIPbendersGetName(topbenders));
      else
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s_%d_%s", AUXILIARYVAR_NAME, i, SCIPbendersGetName(benders));

      /* finding the variable in the copied problem that has the same name as the auxiliary variable */
      targetvar = SCIPfindVar(scip, varname);
      assert(targetvar != NULL);

      SCIPvarSetData(targetvar, vardata);

      benders->auxiliaryvars[i] = SCIPvarGetTransVar(targetvar);

      SCIP_CALL( SCIPcaptureVar(scip, benders->auxiliaryvars[i]) );
   }

   SCIPfreeBlockMemory(scip, &vardata);

   return SCIP_OKAY;
}

/** sets the subproblem objective value array to -infinity */
static
void resetSubproblemObjectiveValue(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   SCIP* subproblem;
   int nsubproblems;
   int i;

   assert(benders != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   for( i = 0; i < nsubproblems; i++ )
   {
      subproblem = SCIPbendersSubproblem(benders, i);
      SCIPbendersSetSubproblemObjval(benders, i, SCIPinfinity(subproblem));
   }
}

/** compares two Benders' decompositions w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPbendersComp)
{  /*lint --e{715}*/
   return ((SCIP_BENDERS*)elem2)->priority - ((SCIP_BENDERS*)elem1)->priority;
}

/** comparison method for sorting Benders' decompositions w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPbendersCompName)
{
   return strcmp(SCIPbendersGetName((SCIP_BENDERS*)elem1), SCIPbendersGetName((SCIP_BENDERS*)elem2));
}

/** method to call, when the priority of a Benders' decomposition was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdBendersPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetBendersPriority() to mark the Benders' decompositions as unsorted */
   SCIPsetBendersPriority(scip, (SCIP_BENDERS*)paramdata, SCIPparamGetInt(param)); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a variable mapping between the master problem variables of the source scip and the sub scip */
static
SCIP_RETCODE createMasterVarMapping(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition of the target SCIP instance */
   SCIP_SET*             sourceset,          /**< global SCIP settings from the source SCIP */
   SCIP_HASHMAP*         varmap              /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables; must not be NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* targetvar;
   int nvars;
   int i;

   assert(benders != NULL);
   assert(sourceset != NULL);
   assert(benders->iscopy);
   assert(benders->mastervarsmap == NULL);

   /* getting the master problem variable data */
   vars = SCIPgetVars(sourceset->scip);
   nvars = SCIPgetNVars(sourceset->scip);

   /* creating the hashmap for the mapping between the master variable of the target and source scip */
   SCIP_CALL( SCIPhashmapCreate(&benders->mastervarsmap, SCIPblkmem(sourceset->scip), nvars) );

   for( i = 0; i < nvars; i++ )
   {
      /* getting the variable pointer for the target SCIP variables. The variable mapping returns the target SCIP
       * varibale for a given source SCIP variable. */
      targetvar = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);
      if( targetvar != NULL )
      {
         SCIP_CALL( SCIPhashmapInsert(benders->mastervarsmap, targetvar, vars[i]) );
         SCIP_CALL( SCIPcaptureVar(sourceset->scip, vars[i]) );
      }
   }

   return SCIP_OKAY;
}

/** copies the given Benders' decomposition to a new SCIP */
SCIP_RETCODE SCIPbendersCopyInclude(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             sourceset,          /**< SCIP_SET of SCIP to copy from */
   SCIP_SET*             targetset,          /**< SCIP_SET of SCIP to copy to */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables; must not be NULL */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   )
{
   SCIP_BENDERS* targetbenders;  /* the copy of the Benders' decomposition struct in the target set */
   int i;

   assert(benders != NULL);
   assert(targetset != NULL);
   assert(varmap != NULL);
   assert(valid != NULL);
   assert(targetset->scip != NULL);

   (*valid) = FALSE;

   if( benders->benderscopy != NULL && targetset->benders_copybenders && SCIPbendersIsActive(benders) )
   {
      SCIPsetDebugMsg(targetset, "including Benders' decomposition %s in subscip %p\n", SCIPbendersGetName(benders), (void*)targetset->scip);
      SCIP_CALL( benders->benderscopy(targetset->scip, benders) );

      /* copying the Benders' cuts */
      targetbenders = SCIPsetFindBenders(targetset, SCIPbendersGetName(benders));

      /* storing the pointer to the source scip instance */
      targetbenders->sourcescip = sourceset->scip;

      /* the flag is set to indicate that the Benders' decomposition is a copy */
      targetbenders->iscopy = TRUE;

      /* calling the copy method for the Benders' cuts */
      SCIPbendersSortBenderscuts(benders);
      for( i = 0; i < benders->nbenderscuts; i++ )
      {
         SCIP_CALL( SCIPbenderscutCopyInclude(targetbenders, benders->benderscuts[i], targetset) );
      }

      /* When the Benders' decomposition is copied then a variable mapping between the master problem variables is
       * required. This variable mapping is used to transfer the cuts generated in the target SCIP to the source SCIP.
       * The variable map is stored in the target Benders' decomposition. This will be freed when the sub-SCIP is freed.
       */
      SCIP_CALL( createMasterVarMapping(targetbenders, sourceset, varmap) );
   }

   /* if the Benders' decomposition is active, then copy is not valid. */
   (*valid) = !SCIPbendersIsActive(benders);

   return SCIP_OKAY;
}

/** creates a Benders' decomposition structure
 *
 *  To use the Benders' decomposition for solving a problem, it first has to be activated with a call to SCIPactivateBenders().
 */
SCIP_RETCODE SCIPbendersCreate(
   SCIP_BENDERS**        benders,            /**< pointer to Benders' decomposition data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   SCIP_Bool             cutlp,              /**< should Benders' cuts be generated for LP solutions */
   SCIP_Bool             cutpseudo,          /**< should Benders' cuts be generated for pseudo solutions */
   SCIP_Bool             cutrelax,           /**< should Benders' cuts be generated for relaxation solutions */
   SCIP_Bool             shareauxvars,       /**< should this Benders' use the highest priority Benders aux vars */
   SCIP_DECL_BENDERSCOPY ((*benderscopy)),   /**< copy method of Benders' decomposition or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree)),   /**< destructor of Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit)),   /**< initialize Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit)),   /**< deinitialize Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre)),/**< presolving initialization method for Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre)),/**< presolving deinitialization method for Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol)),/**< solving process initialization method of Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol)),/**< solving process deinitialization method of Benders' decomposition */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub)),/**< creates a Benders' decomposition subproblem */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve)),/**< called prior to the subproblem solving loop */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex)),/**< the solving method for convex Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub)),/**< the solving method for the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve)),/**< called after the subproblems are solved. */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub)),/**< the freeing method for the Benders' decomposition subproblems */
   SCIP_BENDERSDATA*     bendersdata         /**< Benders' decomposition data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(benders != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   /* Checking whether the benderssolvesub and the bendersfreesub are both implemented or both are not implemented */
   if( (benderssolvesubconvex == NULL && benderssolvesub == NULL && bendersfreesub != NULL)
      || ((benderssolvesubconvex != NULL || benderssolvesub != NULL) && bendersfreesub == NULL) )
   {
      SCIPerrorMessage("Benders' decomposition <%s> requires that if bendersFreesub%s is implemented, then at least "
         "one of bendersSolvesubconvex%s or bendersSolvesub%s are implemented.\n", name, name, name, name);
      return SCIP_INVALIDCALL;
   }

   SCIP_ALLOC( BMSallocMemory(benders) );
   BMSclearMemory(*benders);
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->desc, desc, strlen(desc)+1) );
   (*benders)->priority = priority;
   (*benders)->cutlp = cutlp;
   (*benders)->cutpseudo = cutpseudo;
   (*benders)->cutrelax = cutrelax;
   (*benders)->shareauxvars = shareauxvars;
   (*benders)->benderscopy = benderscopy;
   (*benders)->bendersfree = bendersfree;
   (*benders)->bendersinit = bendersinit;
   (*benders)->bendersexit = bendersexit;
   (*benders)->bendersinitpre = bendersinitpre;
   (*benders)->bendersexitpre = bendersexitpre;
   (*benders)->bendersinitsol = bendersinitsol;
   (*benders)->bendersexitsol = bendersexitsol;
   (*benders)->bendersgetvar = bendersgetvar;
   (*benders)->benderscreatesub = benderscreatesub;
   (*benders)->benderspresubsolve = benderspresubsolve;
   (*benders)->benderssolvesubconvex = benderssolvesubconvex;
   (*benders)->benderssolvesub = benderssolvesub;
   (*benders)->benderspostsolve = benderspostsolve;
   (*benders)->bendersfreesub = bendersfreesub;
   (*benders)->bendersdata = bendersdata;
   SCIP_CALL( SCIPclockCreate(&(*benders)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*benders)->bendersclock, SCIP_CLOCKTYPE_DEFAULT) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of Benders' decomposition <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*benders)->priority, FALSE, priority, INT_MIN/4, INT_MAX/4,
         paramChgdBendersPriority, (SCIP_PARAMDATA*)(*benders)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutlp", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' cuts be generated for LP solutions?", &(*benders)->cutlp, FALSE, cutlp, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutpseudo", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' cuts be generated for pseudo solutions?", &(*benders)->cutpseudo, FALSE, cutpseudo, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutrelax", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' cuts be generated for relaxation solutions?", &(*benders)->cutrelax, FALSE, cutrelax, NULL, NULL) ); /*lint !e740*/

   /* These parameters are left for the user to decide in a settings file. This departs from the usual SCIP convention
    * where the settings available at the creation of the plugin can be set in the function call.
    */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/transfercuts", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' cuts from LNS heuristics be transferred to the main SCIP instance?", &(*benders)->transfercuts,
         FALSE, SCIP_DEFAULT_TRANSFERCUTS, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/lnscheck", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' decomposition be used in LNS heurisics?", &(*benders)->lnscheck, FALSE, SCIP_DEFAULT_LNSCHECK,
         NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/lnsmaxdepth", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "maximum depth at which the LNS check is performed (-1: no limit)", &(*benders)->lnsmaxdepth, TRUE,
         SCIP_DEFAULT_LNSMAXDEPTH, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutsasconss", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should the transferred cuts be added as constraints?", &(*benders)->cutsasconss, FALSE,
         SCIP_DEFAULT_CUTSASCONSS, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/subprobfrac", name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname,
         "fraction of subproblems that are solved in each iteration", &(*benders)->subprobfrac, FALSE,
         SCIP_DEFAULT_SUBPROBFRAC, 0.0, 1.0, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/updateauxvarbound", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should the auxiliary variable bound be updated by solving the subproblem?", &(*benders)->updateauxvarbound,
         FALSE, SCIP_DEFAULT_UPDATEAUXVARBOUND, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}


/** releases the variables that have been captured in the hashmap */
static
SCIP_RETCODE releaseVarMappingHashmapVars(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   int nentries;
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   assert(benders->mastervarsmap != NULL);

   nentries = SCIPhashmapGetNEntries(benders->mastervarsmap);

   for( i = 0; i < nentries; ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      entry = SCIPhashmapGetEntry(benders->mastervarsmap, i);

      if( entry != NULL )
      {
         SCIP_VAR* var;
         var = (SCIP_VAR*) SCIPhashmapEntryGetImage(entry);

         SCIP_CALL( SCIPreleaseVar(scip, &var) );
      }
   }

   return SCIP_OKAY;
}


/** calls destructor and frees memory of Benders' decomposition */
SCIP_RETCODE SCIPbendersFree(
   SCIP_BENDERS**        benders,            /**< pointer to Benders' decomposition data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(*benders != NULL);
   assert(!(*benders)->initialized);
   assert(set != NULL);

   /* call destructor of Benders' decomposition */
   if( (*benders)->bendersfree != NULL )
   {
      SCIP_CALL( (*benders)->bendersfree(set->scip, *benders) );
   }

   /* if the Benders' decomposition is a copy, then the variable map between the source and the target SCIP needs to be
    * freed.
    */
   if( (*benders)->iscopy )
   {
      SCIP_CALL( releaseVarMappingHashmapVars((*benders)->sourcescip, (*benders)) );
      SCIPhashmapFree(&(*benders)->mastervarsmap);
   }

   /* freeing the Benders' cuts */
   for( i = 0; i < (*benders)->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutFree(&((*benders)->benderscuts[i]), set) );
   }
   BMSfreeMemoryArrayNull(&(*benders)->benderscuts);

   SCIPclockFree(&(*benders)->bendersclock);
   SCIPclockFree(&(*benders)->setuptime);
   BMSfreeMemoryArray(&(*benders)->name);
   BMSfreeMemoryArray(&(*benders)->desc);
   BMSfreeMemory(benders);

   return SCIP_OKAY;
}

/** initialises a MIP subproblem by putting the problem into SCIP_STAGE_SOLVING. This is achieved by calling SCIPsolve
 *  and then interrupting the solve in a node focus event handler.
 *  The LP subproblem is also initialised using this method; however, a different event handler is added. This event
 *  handler will put the LP subproblem into probing mode.
 *  The MIP solving function is called to initialise the subproblem because this function calls SCIPsolve with the
 *  appropriate parameter settings for Benders' decomposition.
 */
static
SCIP_RETCODE initialiseSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            success             /**< was the initialisation process successful */
   )
{
   SCIP* subproblem;
   SCIP_Bool infeasible;
   SCIP_Bool cutoff;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));
   assert(success != NULL);

   (*success) = FALSE;

   subproblem = SCIPbendersSubproblem(benders, probnumber);
   assert(subproblem != NULL);

   /* Getting the problem into the right SCIP stage for solving */
   SCIP_CALL( SCIPbendersSolveSubproblemCIP(set->scip, benders, probnumber, &infeasible, SCIP_BENDERSENFOTYPE_LP,
         FALSE) );

   /* Constructing the LP that can be solved in later iterations */
   if( SCIPgetStatus(subproblem) != SCIP_STATUS_BESTSOLLIMIT && SCIPgetStatus(subproblem) != SCIP_STATUS_TIMELIMIT
      && SCIPgetStatus(subproblem) != SCIP_STATUS_MEMLIMIT )
   {
      assert(SCIPgetStage(subproblem) == SCIP_STAGE_SOLVING);

      SCIP_CALL( SCIPconstructLP(subproblem, &cutoff) );
      (*success) = TRUE;
   }


   return SCIP_OKAY;
}


/** initialises an LP subproblem by putting the problem into probing mode. The probing mode is invoked in a node focus
 *  event handler. This event handler is added just prior to calling the initialise subproblem function.
 */
static
SCIP_RETCODE initialiseLPSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_Bool success;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   subproblem = SCIPbendersSubproblem(benders, probnumber);
   assert(subproblem != NULL);

   /* include event handler into SCIP */
   SCIP_CALL( SCIPallocBlockMemory(subproblem, &eventhdlrdata) );

   SCIP_CALL( initEventhandlerData(subproblem, eventhdlrdata) );

   SCIP_CALL( SCIPincludeEventhdlrBasic(subproblem, &eventhdlr, NODEFOCUS_EVENTHDLR_NAME, NODEFOCUS_EVENTHDLR_DESC,
         eventExecBendersNodefocus, eventhdlrdata) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(subproblem, eventhdlr, eventInitsolBendersNodefocus) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(subproblem, eventhdlr, eventExitsolBendersNodefocus) );
   SCIP_CALL( SCIPsetEventhdlrExit(subproblem, eventhdlr, eventExitBendersNodefocus) );
   SCIP_CALL( SCIPsetEventhdlrFree(subproblem, eventhdlr, eventFreeBendersNodefocus) );
   assert(eventhdlr != NULL);

   /* calling an initial solve to put the problem into probing mode */
   SCIP_CALL( initialiseSubproblem(benders, set, probnumber, &success) );

   return SCIP_OKAY;
}

/** creates the subproblems and registers it with the Benders' decomposition struct */
static
SCIP_RETCODE createSubproblems(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP* subproblem;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_VAR* mastervar;
   SCIP_VAR** vars;
   int nvars;
   int nbinvars;
   int nintvars;
   int nimplintvars;
   int nsubproblems;
   int i;
   int j;

   assert(benders != NULL);
   assert(set != NULL);

   /* if the subproblems have already been created, then they will not be created again. This is the case if the
    * transformed problem has been freed and then retransformed. The subproblems should only be created when the problem
    * is first transformed. */
   if( benders->subprobscreated )
      return SCIP_OKAY;

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* creating all subproblems */
   for( i = 0; i < nsubproblems; i++ )
   {
      /* calling the create subproblem call back method */
      SCIP_CALL( benders->benderscreatesub(set->scip, benders, i) );

      subproblem = SCIPbendersSubproblem(benders, i);

      assert(subproblem != NULL);

      /* setting global limits for the subproblems. This overwrites the limits set by the user */
      SCIP_CALL( SCIPsetIntParam(subproblem, "limits/maxorigsol", 0) );

      /* getting the number of integer and binary variables to determine the problem type */
      SCIP_CALL( SCIPgetVarsData(subproblem, &vars, &nvars, &nbinvars, &nintvars, &nimplintvars, NULL) );

      /* The objective function coefficients of the master problem are set to zero. This is necessary for the Benders'
       * decomposition algorithm, since the cut methods and the objective function check assumes that the objective
       * coefficients of the master problem variables are zero.
       *
       * This only occurs if the Benders' decomposition is not a copy. It is assumed that the correct objective
       * coefficients are given during the first subproblem creation.
       */
      if( !benders->iscopy )
      {
         SCIP_Bool objchanged = FALSE;

         assert(SCIPgetStage(subproblem) == SCIP_STAGE_PROBLEM);
         for( j = 0; j < nvars; j++ )
         {
            /* retrieving the master problem variable */
            SCIP_CALL( SCIPbendersGetVar(benders, set, vars[j], &mastervar, -1) );

            /* if mastervar is not NULL, then the subproblem variable has a corresponding master problem variable */
            if( mastervar != NULL && !SCIPisZero(subproblem, SCIPvarGetObj(vars[j])) )
            {
               SCIPverbMessage(subproblem, SCIP_VERBLEVEL_FULL, NULL, "Benders' decomposition: Changing the objective "
                  "coefficient of copy of master problem variable <%s> in subproblem %d to zero.\n",
                  SCIPvarGetName(mastervar), i);
               /* changing the subproblem variable objective coefficient to zero */
               SCIP_CALL( SCIPchgVarObj(subproblem, vars[j], 0.0) );

               objchanged = TRUE;
            }
         }

         if( objchanged )
         {
            SCIPverbMessage(subproblem, SCIP_VERBLEVEL_HIGH, NULL, "Benders' decomposition: Objective coefficients of "
               "copied of master problem variables has been changed to zero.\n");
         }
      }

      /* if there are no binary and integer variables, then the subproblem is an LP.
       * In this case, the SCIP instance is put into probing mode via the use of an event handler. */
      if( nbinvars == 0 && nintvars == 0 && nimplintvars == 0 )
      {
         SCIPbendersSetSubproblemIsConvex(benders, i, TRUE);

         /* if the user has not implemented a solve subproblem callback, then the subproblem solves are performed
          * internally. To be more efficient the subproblem is put into probing mode. */
         if( benders->benderssolvesubconvex == NULL && benders->benderssolvesub == NULL
            && SCIPgetStage(subproblem) <= SCIP_STAGE_PROBLEM )
         {
            SCIP_CALL( initialiseLPSubproblem(benders, set, i) );
         }
      }
      else
      {
         SCIP_EVENTHDLRDATA* eventhdlrdata_mipnodefocus;
         SCIP_EVENTHDLRDATA* eventhdlrdata_upperbound;

         SCIPbendersSetSubproblemIsConvex(benders, i, FALSE);

         /* because the subproblems could be reused in the copy, the event handler is not created again.
          * NOTE: This currently works with the benders_default implementation. It may not be very general. */
         if( benders->benderssolvesubconvex == NULL && benders->benderssolvesub == NULL && !benders->iscopy )
         {
            SCIP_CALL( SCIPallocBlockMemory(subproblem, &eventhdlrdata_mipnodefocus) );
            SCIP_CALL( SCIPallocBlockMemory(subproblem, &eventhdlrdata_upperbound) );

            SCIP_CALL( initEventhandlerData(subproblem, eventhdlrdata_mipnodefocus) );
            SCIP_CALL( initEventhandlerData(subproblem, eventhdlrdata_upperbound) );

            /* include the first LP solved event handler into the subproblem */
            SCIP_CALL( SCIPincludeEventhdlrBasic(subproblem, &eventhdlr, MIPNODEFOCUS_EVENTHDLR_NAME,
                  MIPNODEFOCUS_EVENTHDLR_DESC, eventExecBendersMipnodefocus, eventhdlrdata_mipnodefocus) );
            SCIP_CALL( SCIPsetEventhdlrInitsol(subproblem, eventhdlr, eventInitsolBendersMipnodefocus) );
            SCIP_CALL( SCIPsetEventhdlrExitsol(subproblem, eventhdlr, eventExitsolBendersMipnodefocus) );
            SCIP_CALL( SCIPsetEventhdlrExit(subproblem, eventhdlr, eventExitBendersMipnodefocus) );
            SCIP_CALL( SCIPsetEventhdlrFree(subproblem, eventhdlr, eventFreeBendersMipnodefocus) );
            assert(eventhdlr != NULL);

            /* include the upper bound interrupt event handler into the subproblem */
            SCIP_CALL( SCIPincludeEventhdlrBasic(subproblem, &eventhdlr, UPPERBOUND_EVENTHDLR_NAME,
                  UPPERBOUND_EVENTHDLR_DESC, eventExecBendersUpperbound, eventhdlrdata_upperbound) );
            SCIP_CALL( SCIPsetEventhdlrInitsol(subproblem, eventhdlr, eventInitsolBendersUpperbound) );
            SCIP_CALL( SCIPsetEventhdlrExitsol(subproblem, eventhdlr, eventExitsolBendersUpperbound) );
            SCIP_CALL( SCIPsetEventhdlrExit(subproblem, eventhdlr, eventExitBendersUpperbound) );
            SCIP_CALL( SCIPsetEventhdlrFree(subproblem, eventhdlr, eventFreeBendersUpperbound) );
            assert(eventhdlr != NULL);
         }
      }
   }

   benders->subprobscreated = TRUE;

   return SCIP_OKAY;
}


/** initializes Benders' decomposition */
SCIP_RETCODE SCIPbendersInit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   if( benders->initialized )
   {
      SCIPerrorMessage("Benders' decomposition <%s> already initialized\n", benders->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(benders->setuptime);
      SCIPclockReset(benders->bendersclock);

      benders->ncalls = 0;
      benders->ncutsfound = 0;
      benders->ntransferred = 0;
   }

   /* start timing */
   SCIPclockStart(benders->setuptime, set);

   if( benders->bendersinit != NULL )
   {
      SCIP_CALL( benders->bendersinit(set->scip, benders) );
   }

   benders->initialized = TRUE;

   /* creates the subproblems and sets up the probing mode for LP subproblems. This function calls the benderscreatesub
    * callback. */
   SCIP_CALL( createSubproblems(benders, set) );

   /* initialising the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutInit(benders->benderscuts[i], set) );
   }

   /* stop timing */
   SCIPclockStop(benders->setuptime, set);

   return SCIP_OKAY;
}


/** Transfers Benders' cuts that were generated while solving a sub-SCIP to the original SCIP instance. This involves
 *  creating a constraint/cut that is equivalent to the generated cut in the sub-SCIP. This new constraint/cut is then
 *  added to the original SCIP instance.
 */
static
SCIP_RETCODE createAndAddTransferredCut(
   SCIP*                 sourcescip,         /**< the source SCIP from when the Benders' decomposition was copied */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure of the sub SCIP */
   SCIP_VAR**            vars,               /**< the variables from the source constraint */
   SCIP_Real*            vals,               /**< the coefficients of the variables in the source constriant */
   SCIP_Real             lhs,                /**< the LHS of the source constraint */
   SCIP_Real             rhs,                /**< the RHS of the source constraint */
   int                   nvars               /**< the number of variables in the source constraint */
   )
{
   SCIP_BENDERS* sourcebenders;     /* the Benders' decomposition of the source SCIP */
   SCIP_CONSHDLR* consbenders;      /* a helper variable for the Benders' decomposition constraint handler */
   SCIP_CONS* transfercons;         /* the constraint that is generated to transfer the constraints/cuts */
   SCIP_ROW* transfercut;           /* the cut that is generated to transfer the constraints/cuts */
   SCIP_VAR* sourcevar;             /* the source variable that will be added to the transferred cut */
   SCIP_VAR* origvar;
   SCIP_Real scalar;
   SCIP_Real constant;
   char cutname[SCIP_MAXSTRLEN];    /* the name of the transferred cut */
   int i;
   SCIP_Bool fail;

   assert(sourcescip != NULL);
   assert(benders != NULL);
   assert(vars != NULL);
   assert(vals != NULL);

   /* retrieving the source Benders' decomposition structure */
   sourcebenders = SCIPfindBenders(sourcescip, SCIPbendersGetName(benders));

   /* retrieving the Benders' decomposition constraint handler */
   consbenders = SCIPfindConshdlr(sourcescip, "benders");

   /* setting the name of the transferred cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "transferredcut_%d",
      SCIPbendersGetNTransferredCuts(sourcebenders) );

   /* TODO: It could be more efficient to pass an updated vars array with the vals array to the
    * SCIPcreateConsBasicLinear/SCIPcreateEmptyRowCons. This should be implemented to improve the performance of the
    * Large Neighbourhood Benders Search.
    */

   /* creating an empty row/constraint for the transferred cut */
   if( sourcebenders->cutsasconss )
   {
      SCIP_CALL( SCIPcreateConsBasicLinear(sourcescip, &transfercons, cutname, 0, NULL, NULL, lhs, rhs) );
      SCIP_CALL( SCIPsetConsRemovable(sourcescip, transfercons, TRUE) );
   }
   else
   {
      SCIP_CALL( SCIPcreateEmptyRowCons(sourcescip, &transfercut, consbenders, cutname, lhs, rhs, FALSE,
            FALSE, TRUE) );
   }

   fail = FALSE;
   for( i = 0; i < nvars; i++ )
   {
      /* getting the original variable for the transformed variable */
      origvar = vars[i];
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      /* getting the source var from the hash map */
      sourcevar = (SCIP_VAR*) SCIPhashmapGetImage(benders->mastervarsmap, origvar);

      /* if the source variable is not found, then the mapping in incomplete. So the constraint can not be
       * transferred. */
      if( sourcevar == NULL )
      {
         fail = TRUE;
         break;
      }

      if( sourcebenders->cutsasconss )
      {
         SCIP_CALL( SCIPaddCoefLinear(sourcescip, transfercons, sourcevar, vals[i]) );    /*lint !e644*/
      }
      else
      {
         SCIP_CALL( SCIPaddVarToRow(sourcescip, transfercut, sourcevar, vals[i]) );       /*lint !e644*/
      }
   }

   /* if all of the source variables were found to generate the cut */
   if( !fail )
   {
      if( sourcebenders->cutsasconss )
      {
         SCIP_CALL( SCIPaddCons(sourcescip, transfercons) );
      }
      else
      {
         SCIP_CALL( SCIPaddPoolCut(sourcescip, transfercut) );
      }

      sourcebenders->ntransferred++;
   }

   /* release the row/constraint */
   if( sourcebenders->cutsasconss )
   {
      SCIP_CALL( SCIPreleaseCons(sourcescip, &transfercons) );
   }
   else
   {
      SCIP_CALL( SCIPreleaseRow(sourcescip, &transfercut) );
   }

   return SCIP_OKAY;
}


/** transfers the cuts generated in a subscip to the source scip */
static
SCIP_RETCODE transferBendersCuts(
   SCIP*                 sourcescip,         /**< the source SCIP from when the Benders' decomposition was copied */
   SCIP*                 subscip,            /**< the sub SCIP where the Benders' cuts were generated */
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure of the sub SCIP */
   )
{
   SCIP_BENDERS* sourcebenders;     /* the Benders' decomposition of the source SCIP */
   SCIP_BENDERSCUT* benderscut;     /* a helper variable for the Benders' cut plugin */
   SCIP_CONS** addedcons;           /* the constraints added by the Benders' cut */
   SCIP_ROW** addedcuts;            /* the cuts added by the Benders' cut */
   SCIP_VAR** vars;                 /* the variables of the added constraint/row */
   SCIP_Real* vals;                 /* the values of the added constraint/row */
   SCIP_Real lhs;                   /* the LHS of the added constraint/row */
   SCIP_Real rhs;                   /* the RHS of the added constraint/row */
   int naddedcons;
   int naddedcuts;
   int nvars;
   int i;
   int j;
   int k;

   assert(subscip != NULL);
   assert(benders != NULL);

   /* retrieving the source Benders' decomposition structure */
   sourcebenders = SCIPfindBenders(sourcescip, SCIPbendersGetName(benders));

   /* exit if the cuts should not be transferred from the sub SCIP to the source SCIP. */
   if( !sourcebenders->transfercuts )
      return SCIP_OKAY;

   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      benderscut = benders->benderscuts[i];

      /* retreiving the Benders' cuts constraints */
      SCIP_CALL( SCIPbenderscutGetAddedConss(benderscut, &addedcons, &naddedcons) );

      /* looping over all added constraints to construct the cut for the source scip */
      for( j = 0; j < naddedcons; j++ )
      {
         SCIP_CONSHDLR* conshdlr;
         const char * conshdlrname;

         conshdlr = SCIPconsGetHdlr(addedcons[j]);
         assert(conshdlr != NULL);
         conshdlrname = SCIPconshdlrGetName(conshdlr);

         /* it is only possible to transfer linear constraints. If the Benders' cut has been added as another
          * constraint, then this will not be transferred to the source SCIP */
         if( strcmp(conshdlrname, "linear") == 0 )
         {
            /* collecting the variable information from the constraint */
            nvars = SCIPgetNVarsLinear(subscip, addedcons[j]);
            vars = SCIPgetVarsLinear(subscip, addedcons[j]);
            vals = SCIPgetValsLinear(subscip, addedcons[j]);

            /* collecting the bounds from the constraint */
            lhs = SCIPgetLhsLinear(subscip, addedcons[j]);
            rhs = SCIPgetRhsLinear(subscip, addedcons[j]);

            if( nvars > 0 )
            {
               /* create and add the cut to be transferred from the sub SCIP to the source SCIP */
               SCIP_CALL( createAndAddTransferredCut(sourcescip, benders, vars, vals, lhs, rhs, nvars) );
            }
         }
      }

      /* retreiving the Benders' cuts added cuts */
      SCIP_CALL( SCIPbenderscutGetAddedCuts(benderscut, &addedcuts, &naddedcuts) );

      /* looping over all added constraints to costruct the cut for the source scip */
      for( j = 0; j < naddedcuts; j++ )
      {
         SCIP_COL** cols;
         int ncols;

         cols = SCIProwGetCols(addedcuts[j]);
         ncols = SCIProwGetNNonz(addedcuts[j]);

         /* get all variables of the row */
         SCIP_CALL( SCIPallocBufferArray(subscip, &vars, ncols) );
         for( k = 0; k < ncols; ++k )
            vars[k] = SCIPcolGetVar(cols[k]);

         /* collecting the variable information from the constraint */
         vals = SCIProwGetVals(addedcuts[j]);

         /* collecting the bounds from the constraint */
         lhs = SCIProwGetLhs(addedcuts[j]) - SCIProwGetConstant(addedcuts[j]);
         rhs = SCIProwGetRhs(addedcuts[j]) - SCIProwGetConstant(addedcuts[j]);

         /* create and add the cut to be transferred from the sub SCIP to the source SCIP */
         SCIP_CALL( createAndAddTransferredCut(sourcescip, benders, vars, vals, lhs, rhs, ncols) );

         SCIPfreeBufferArray(subscip, &vars);
      }
   }

   return SCIP_OKAY;
}


/** calls exit method of Benders' decomposition */
SCIP_RETCODE SCIPbendersExit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int nsubproblems;
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   if( !benders->initialized )
   {
      SCIPerrorMessage("Benders' decomposition <%s> not initialized\n", benders->name);
      return SCIP_INVALIDCALL;
   }

   /* start timing */
   SCIPclockStart(benders->setuptime, set);

   if( benders->bendersexit != NULL )
   {
      SCIP_CALL( benders->bendersexit(set->scip, benders) );
   }

   /* if the Benders' decomposition is a copy, then the generated cuts will be transferred to the source scip */
   if( benders->iscopy )
   {
      SCIP_CALL( transferBendersCuts(benders->sourcescip, set->scip, benders) );
   }

   /* releasing all of the auxiliary variables */
   nsubproblems = SCIPbendersGetNSubproblems(benders);
   for( i = 0; i < nsubproblems; i++ )
   {
      /* it is possible that the master problem is not solved. As such, the auxiliary variables will not be created. So
       * we don't need to release the variables
       */
      if( benders->auxiliaryvars[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(set->scip, &benders->auxiliaryvars[i]) );
      }
   }

   /* calling the exit method for the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutExit(benders->benderscuts[i], set) );
   }

   benders->initialized = FALSE;

   /* stop timing */
   SCIPclockStop(benders->setuptime, set);

   return SCIP_OKAY;
}

/** Checks whether a subproblem is independent. */
static
SCIP_RETCODE checkSubproblemIndependence(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int nsubproblems;
   int i;
   int j;

   assert(scip != NULL);
   assert(benders != NULL);

   /* retrieving the master problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* looping over all subproblems to check whether there exists at least one master problem variable */
   for( i = 0; i < nsubproblems; i++ )
   {
      SCIP_Bool independent = FALSE;

      /* if there are user defined solving or freeing functions, then it is not possible to declare the independence of
       * the subproblems.
       */
      if( benders->benderssolvesubconvex == NULL && benders->benderssolvesub == NULL
         && benders->bendersfreesub == NULL )
      {
         independent = TRUE;

         for( j = 0; j < nvars; j++ )
         {
            SCIP_VAR* subprobvar;

            /* getting the subproblem problem variable corresponding to the master problem variable */
            SCIP_CALL( SCIPgetBendersSubproblemVar(scip, benders, vars[j], &subprobvar, i) );

            /* if the subporblem variable is not NULL, then the subproblem depends on the master problem */
            if( subprobvar != NULL )
            {
               independent = FALSE;
               break;
            }
         }

         /* setting the independent flag */
         SCIPbendersSetSubproblemIsIndependent(benders, i, independent);
      }
   }

   return SCIP_OKAY;
}

/** informs the Benders' decomposition that the presolving process is being started */
SCIP_RETCODE SCIPbendersInitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(stat != NULL);

   if( !benders->iscopy )
   {
      /* check the subproblem independence. This check is only performed if the user has not implemented a solve
       * subproblem function.
       */
      if( benders->benderssolvesubconvex == NULL && benders->benderssolvesub == NULL )
        SCIP_CALL( checkSubproblemIndependence(set->scip, benders) );

      /* adding the auxiliary variables to the master problem */
      SCIP_CALL( addAuxiliaryVariablesToMaster(set->scip, benders) );
   }
   else
   {
      /* the copied auxiliary variables must be assigned to the target Benders' decomposition */
      SCIP_CALL( assignAuxiliaryVariables(set->scip, benders) );
   }

   /* call presolving initialization method of Benders' decomposition */
   if( benders->bendersinitpre != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersinitpre(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   return SCIP_OKAY;
}


/** informs the Benders' decomposition that the presolving process has completed */
SCIP_RETCODE SCIPbendersExitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(stat != NULL);

   /* call presolving  deinitialization method of Benders' decomposition */
   if( benders->bendersexitpre != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersexitpre(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs Benders' decomposition that the branch and bound process is being started */
SCIP_RETCODE SCIPbendersInitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   /* call solving process initialization method of Benders' decomposition */
   if( benders->bendersinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersinitsol(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   /* calling the initsol method for the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutInitsol(benders->benderscuts[i], set) );
   }

   return SCIP_OKAY;
}

/** informs Benders' decomposition that the branch and bound process data is being freed */
SCIP_RETCODE SCIPbendersExitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int nsubproblems;
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);
   /* freeing all subproblems that are independent, this is because they have not bee freed during the subproblem
    * solving loop.
    */
   for( i = 0; i < nsubproblems; i++ )
   {
      if( SCIPbendersSubproblemIsIndependent(benders, i) )
      {
         /* disabling the independence of the subproblem so that it can be freed */
         SCIPbendersSetSubproblemIsIndependent(benders, i, FALSE);

         /* freeing the independent subproblem */
         SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, i) );
      }
   }

   /* call solving process deinitialization method of Benders' decomposition */
   if( benders->bendersexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersexitsol(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   /* sorting the Benders' decomposition cuts in order of priority. Only a single cut is generated for each subproblem
    * per solving iteration. This is particularly important in the case of the optimality and feasibility cuts. Since
    * these work on two different solutions to the subproblem, it is not necessary to generate both cuts. So, once the
    * feasibility cut is generated, then no other cuts will be generated.
    */
   SCIPbendersSortBenderscuts(benders);

   /* calling the exitsol method for the Benders' cuts */
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutExitsol(benders->benderscuts[i], set) );
   }

   return SCIP_OKAY;
}

/** activates Benders' decomposition such that it is called in LP solving loop */
SCIP_RETCODE SCIPbendersActivate(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nsubproblems        /**< the number subproblems used in this decomposition */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   int i;

   assert(benders != NULL);
   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_INIT || set->stage == SCIP_STAGE_PROBLEM);

   if( !benders->active )
   {
      benders->active = TRUE;
      set->nactivebenders++;
      set->benderssorted = FALSE;

      benders->nsubproblems = nsubproblems;
      benders->nactivesubprobs = nsubproblems;

      /* allocating memory for the subproblems arrays */
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subproblems, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->auxiliaryvars, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobobjval, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->bestsubprobobjval, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subproblowerbound, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobisconvex, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobsetup, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->indepsubprob, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobenabled, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->mastervarscont, benders->nsubproblems) );

      for( i = 0; i < benders->nsubproblems; i++ )
      {
         benders->subproblems[i] = NULL;
         benders->auxiliaryvars[i] = NULL;
         benders->subprobobjval[i] = SCIPsetInfinity(set);
         benders->bestsubprobobjval[i] = SCIPsetInfinity(set);
         benders->subproblowerbound[i] = -SCIPsetInfinity(set);
         benders->subprobisconvex[i] = FALSE;
         benders->subprobsetup[i] = FALSE;
         benders->indepsubprob[i] = FALSE;
         benders->subprobenabled[i] = TRUE;
         benders->mastervarscont[i] = FALSE;
      }

      /* adding an eventhandler for updating the lower bound when the root node is solved. */
      eventhdlrdata = (SCIP_EVENTHDLRDATA*)benders;

      /* include event handler into SCIP */
      SCIP_CALL( SCIPincludeEventhdlrBasic(set->scip, &eventhdlr, NODESOLVED_EVENTHDLR_NAME, NODESOLVED_EVENTHDLR_DESC,
            eventExecBendersNodesolved, eventhdlrdata) );
      SCIP_CALL( SCIPsetEventhdlrInitsol(set->scip, eventhdlr, eventInitsolBendersNodesolved) );
      assert(eventhdlr != NULL);
   }

   return SCIP_OKAY;
}

/** deactivates Benders' decomposition such that it is no longer called in LP solving loop */
void SCIPbendersDeactivate(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_INIT || set->stage == SCIP_STAGE_PROBLEM);

   if( benders->active )
   {
#ifndef NDEBUG
      int nsubproblems;
      int i;

      nsubproblems = SCIPbendersGetNSubproblems(benders);

      /* checking whether the auxiliary variables and subproblems are all NULL */
      for( i = 0; i < nsubproblems; i++ )
         assert(benders->auxiliaryvars[i] == NULL);
#endif

      benders->active = FALSE;
      set->nactivebenders--;
      set->benderssorted = FALSE;

      /* freeing the memory allocated during the activation of the Benders' decomposition */
      BMSfreeMemoryArray(&benders->mastervarscont);
      BMSfreeMemoryArray(&benders->subprobenabled);
      BMSfreeMemoryArray(&benders->indepsubprob);
      BMSfreeMemoryArray(&benders->subprobsetup);
      BMSfreeMemoryArray(&benders->subprobisconvex);
      BMSfreeMemoryArray(&benders->subproblowerbound);
      BMSfreeMemoryArray(&benders->bestsubprobobjval);
      BMSfreeMemoryArray(&benders->subprobobjval);
      BMSfreeMemoryArray(&benders->auxiliaryvars);
      BMSfreeMemoryArray(&benders->subproblems);
   }
}

/** returns whether the given Benders' decomposition is in use in the current problem */
SCIP_Bool SCIPbendersIsActive(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   assert(benders != NULL);

   return benders->active;
}

/** updates the lower bound for all auxiliary variables. This is called if the first LP enforced is unbounded. */
static
SCIP_RETCODE updateAuxiliaryVarLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RESULT*          result              /**< the result from updating the auxiliary variable lower bound */
   )
{
   int nsubproblems;
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   (*result) = SCIP_DIDNOTRUN;

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   for( i = 0; i < nsubproblems; i++ )
   {
      SCIP_VAR* auxiliaryvar;
      SCIP_Real lowerbound;
      SCIP_Bool infeasible;

      infeasible = FALSE;

      /* computing the lower bound of the subproblem by solving it without any variable fixings */
      SCIP_CALL( SCIPbendersComputeSubproblemLowerbound(benders, set, i, &lowerbound, &infeasible) );

      /* if the subproblem is infeasible, then the original problem is infeasible */
      if( infeasible )
      {
         (*result) = SCIP_INFEASIBLE;
         break;
      }

      /* retrieving the auxiliary variable */
      auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, i);

      /* only update the lower bound if it is greater than the current lower bound */
      if( SCIPsetIsGT(set, lowerbound, SCIPvarGetLbGlobal(auxiliaryvar)) )
      {
         SCIPsetDebugMsg(set, "Tightened lower bound of <%s> to %g\n", SCIPvarGetName(auxiliaryvar), lowerbound);
         /* updating the lower bound of the auxiliary variable */
         SCIP_CALL( SCIPchgVarLb(set->scip, auxiliaryvar, lowerbound) );
         (*result) = SCIP_REDUCEDDOM;
      }

      /* stores the lower bound for the subproblem */
      SCIPbendersUpdateSubproblemLowerbound(benders, i, lowerbound);
   }

   return SCIP_OKAY;
}

/** Returns whether only the convex relaxations will be checked in this solve loop
 *  when Benders' is used in the LNS heuristics, only the convex relaxations of the master/subproblems are checked,
 *  i.e. no integer cuts are generated. In this case, then Benders' decomposition is performed under the assumption
 *  that all subproblems are convex relaxations.
 */
SCIP_Bool SCIPbendersOnlyCheckConvexRelax(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   return benders->iscopy && benders->lnscheck;
}

/** returns the number of subproblems that will be checked in this iteration */
static
int numSubproblemsToCheck(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERSENFOTYPE  type                /**< the type of solution being enforced */
   )
{
   if( benders->ncalls == 0 || type == SCIP_BENDERSENFOTYPE_CHECK || SCIPbendersOnlyCheckConvexRelax(benders) )
      return SCIPbendersGetNSubproblems(benders);
   else
      return (int) SCIPsetCeil(set, (SCIP_Real) SCIPbendersGetNSubproblems(benders)*benders->subprobfrac);
}

/** returns whether the solving of the given subproblem needs to be executed */
static
SCIP_Bool subproblemIsActive(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem index */
   )
{
   return (!SCIPbendersSubproblemIsIndependent(benders, probnumber)
      && SCIPbendersSubproblemIsEnabled(benders, probnumber));
}

/** Solves each of the Benders' decomposition subproblems for the given solution. All, or a fraction, of subproblems are
 *  solved before the Benders' decomposition cuts are generated.
 *  Since a convex relaxation of the subproblem could be solved to generate cuts, a parameter nverified is used to
 *  identified the number of subproblems that have been solved in their "original" form. For example, if the subproblem
 *  is a MIP, then if the LP is solved to generate cuts, this does not constitute a verification. The verification is
 *  only performed when the MIP is solved.
 */
static
SCIP_RETCODE solveBendersSubproblems(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_BENDERSSOLVELOOP solveloop,          /**< the current solve loop */
   SCIP_Bool             checkint,           /**< are the subproblems called during a check/enforce of integer sols? */
   int*                  nchecked,           /**< the number of subproblems checked in this solve loop, they may not be solved */
   int*                  nverified,          /**< the number of subproblems verified in the current loop */
   SCIP_Bool**           subprobsolved,      /**< an array indicating the subproblems that were solved in this loop. */
   SCIP_BENDERSSUBSTATUS** substatus,        /**< array to store the status of the subsystem */
   SCIP_Bool*            infeasible,         /**< is the master problem infeasible with respect to the Benders' cuts? */
   SCIP_Bool*            optimal,            /**< is the current solution optimal? */
   SCIP_Bool*            stopped             /**< was the solving process stopped? */
   )
{
   SCIP_Bool onlyconvexcheck;
   int nsubproblems;
   int numtocheck;
   int numnotopt;
   int subproblemcount;
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   (*stopped) = FALSE;

   /* getting the number of subproblems in the Benders' decompsition */
   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* in the case of an LNS check, only the convex relaxations of the subproblems will be solved. This is a performance
    * feature, since solving the convex relaxation is typically much faster than solving the corresponding CIP. While
    * the CIP is not solved during the LNS check, the solutions are still of higher quality than when Benders' is not
    * employed.
    */
   onlyconvexcheck = SCIPbendersOnlyCheckConvexRelax(benders);

   /* it is possible to only solve a subset of subproblems. This is given by a parameter. */
   numtocheck = numSubproblemsToCheck(benders, set, type);

   SCIPsetDebugMsg(set, "Performing the subproblem solving process. Number of subproblems to check %d\n", numtocheck);

   SCIPsetDebugMsg(set, "Benders' decomposition - solve loop %d\n", solveloop);
   numnotopt = 0;
   subproblemcount = 0;

   if( type == SCIP_BENDERSENFOTYPE_CHECK && sol == NULL )
   {
      /* TODO: Check whether this is absolutely necessary. I think that this if statment can be removed. */
      (*infeasible) = TRUE;
   }
   else
   {
      /* solving each of the subproblems for Benders' decomposition */
      /* TODO: ensure that the each of the subproblems solve and update the parameters with the correct return values
       */
      i = benders->firstchecked;
      /*for( i = 0; i < nsubproblems; i++ )*/
      while( subproblemcount < nsubproblems && numnotopt < numtocheck && !(*stopped) )
      {
         SCIP_Bool subinfeas = FALSE;
         SCIP_Bool convexsub = SCIPbendersSubproblemIsConvex(benders, i);
         SCIP_Bool solvesub = TRUE;
         SCIP_Bool solved;

         /* the subproblem is initially flagged as not solved for this solving loop */
         (*subprobsolved)[i] = FALSE;

         /* setting the subsystem status to UNKNOWN at the start of each solve loop */
         (*substatus)[i] = SCIP_BENDERSSUBSTATUS_UNKNOWN;

         /* for the second solving loop, if the problem is an LP, it is not solved again. If the problem is a MIP,
          * then the subproblem objective function value is set to infinity. However, if the subproblem is proven
          * infeasible from the LP, then the IP loop is not performed.
          * If the solve loop is SCIP_BENDERSSOLVELOOP_USERCIP, then nothing is done. It is assumed that the user will
          * correctly update the objective function within the user-defined solving function.
          */
         if( solveloop == SCIP_BENDERSSOLVELOOP_CIP )
         {
            if( convexsub || (*substatus)[i] == SCIP_BENDERSSUBSTATUS_INFEAS )
               solvesub = FALSE;
            else
               SCIPbendersSetSubproblemObjval(benders, i, SCIPinfinity(SCIPbendersSubproblem(benders, i)));
         }

         /* if the subproblem is independent, then it does not need to be solved. In this case, the nverified flag will
          * increase by one. When the subproblem is not independent, then it needs to be checked.
          */
         if( !subproblemIsActive(benders, i) )
         {
            /* NOTE: There is no need to update the optimal flag. This is because optimal is always TRUE until a
             * non-optimal subproblem is found.
             */
            /* if the auxiliary variable value is infinity, then the subproblem has not been solved yet. Currently the
             * subproblem statue is unknown. */
            if( SCIPsetIsInfinity(set, SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i))
               || SCIPsetIsInfinity(set, -SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i))
               || SCIPsetIsInfinity(set, -SCIPbendersGetSubproblemLowerbound(benders, i)) )
            {
               SCIPbendersSetSubproblemObjval(benders, i, SCIPinfinity(SCIPbendersSubproblem(benders, i)));
               (*substatus)[i] = SCIP_BENDERSSUBSTATUS_UNKNOWN;
               (*optimal) = FALSE;

               SCIPsetDebugMsg(set, "Benders' decomposition: subproblem %d is not active, but has not been solved."
                 " setting status to UNKNOWN\n", i);
            }
            else
            {
               SCIP_Real soltol;

               SCIP_CALL( SCIPsetGetRealParam(set, "benders/solutiontol", &soltol) );

               if( SCIPrelDiff(SCIPbendersGetSubproblemLowerbound(benders, i),
                     SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i)) < soltol )
               {
                  SCIPbendersSetSubproblemObjval(benders, i, SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i));
                  (*substatus)[i] = SCIP_BENDERSSUBSTATUS_OPTIMAL;
               }
               else
               {
                  SCIPbendersSetSubproblemObjval(benders, i, SCIPbendersGetSubproblemLowerbound(benders, i));
                  (*substatus)[i] = SCIP_BENDERSSUBSTATUS_AUXVIOL;
               }

               SCIPsetDebugMsg(set, "Benders' decomposition: subproblem %d is not active, setting status to OPTIMAL\n",
                  i);

            }

            (*subprobsolved)[i] = TRUE;

            /* the nverified counter is only increased in the convex solving loop */
            if( solveloop == SCIP_BENDERSSOLVELOOP_CONVEX || solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX )
               (*nverified)++;
         }
         else if( solvesub )
         {
            SCIP_CALL( SCIPbendersExecSubproblemSolve(benders, set, sol, i, solveloop, FALSE, &solved, &subinfeas, type) );

#ifdef SCIP_DEBUG
            if( type == SCIP_BENDERSENFOTYPE_LP )
            {
               SCIPsetDebugMsg(set, "LP: Subproblem %d (%f < %f)\n", i, SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i),
                  SCIPbendersGetSubproblemObjval(benders, i));
            }
#endif
            (*subprobsolved)[i] = solved;

            (*infeasible) = (*infeasible) || subinfeas;
            if( subinfeas )
               (*substatus)[i] = SCIP_BENDERSSUBSTATUS_INFEAS;

            /* if the subproblems are solved to check integer feasibility, then the optimality check must be performed.
             * This will only be performed if checkint is TRUE and the subproblem was solved. The subproblem may not be
             * solved if the user has defined a solving function
             */
            if( checkint && (*subprobsolved)[i] )
            {
               /* if the subproblem is feasible, then it is necessary to update the value of the auxiliary variable to the
                * objective function value of the subproblem.
                */
               if( !subinfeas )
               {
                  SCIP_Bool subproboptimal = FALSE;

                  SCIP_CALL( SCIPbendersCheckSubproblemOptimality(benders, set, sol, i, &subproboptimal) );

                  if( subproboptimal )
                     (*substatus)[i] = SCIP_BENDERSSUBSTATUS_OPTIMAL;
                  else
                     (*substatus)[i] = SCIP_BENDERSSUBSTATUS_AUXVIOL;

                  /* It is only possible to determine the optimality of a solution within a given subproblem in four
                   * different cases:
                   * i) solveloop == SCIP_BENDERSSOLVELOOP_CONVEX or USERCONVEX and the subproblem is convex.
                   * ii) solveloop == SCIP_BENDERSOLVELOOP_CONVEX  and only the convex relaxations will be checked.
                   * iii) solveloop == SCIP_BENDERSSOLVELOOP_USERCIP and the subproblem was solved, since the user has
                   * defined a solve function, it is expected that the solving is correctly executed.
                   * iv) solveloop == SCIP_BENDERSSOLVELOOP_CIP and the MIP for the subproblem has been solved.
                   */
                  if( convexsub || onlyconvexcheck
                     || solveloop == SCIP_BENDERSSOLVELOOP_CIP
                     || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP )
                     (*optimal) = (*optimal) && subproboptimal;

#ifdef SCIP_DEBUG
                  if( convexsub || solveloop >= SCIP_BENDERSSOLVELOOP_CIP )
                  {
                     if( subproboptimal )
                     {
                        SCIPsetDebugMsg(set, "Subproblem %d is Optimal (%f >= %f)\n", i,
                           SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i), SCIPbendersGetSubproblemObjval(benders, i));
                     }
                     else
                     {
                        SCIPsetDebugMsg(set, "Subproblem %d is NOT Optimal (%f < %f)\n", i,
                           SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i), SCIPbendersGetSubproblemObjval(benders, i));
                     }
                  }
#endif

                  /* the nverified variable is only incremented when the original form of the subproblem has been solved.
                   * What is meant by "original" is that the LP relaxation of CIPs are solved to generate valid cuts. So
                   * if the subproblem is defined as a CIP, then it is only classified as checked if the CIP is solved.
                   * There are three cases where the "original" form is solved are:
                   * i) solveloop == SCIP_BENDERSSOLVELOOP_CONVEX or USERCONVEX and the subproblem is an LP
                   *    - the original form has been solved.
                   * ii) solveloop == SCIP_BENDERSSOLVELOOP_CIP or USERCIP and the CIP for the subproblem has been
                   *    solved.
                   * iii) or, only a convex check is performed.
                   */
                  if( ((solveloop == SCIP_BENDERSSOLVELOOP_CONVEX || solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX)
                        && convexsub)
                     || ((solveloop == SCIP_BENDERSSOLVELOOP_CIP || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP)
                        && !convexsub)
                     || onlyconvexcheck )
                     (*nverified)++;

                  if( !subproboptimal )
                  {
                     numnotopt++;
                     assert(numnotopt <= nsubproblems);
                  }
               }
               else
               {
                  numnotopt++;
                  assert(numnotopt <= nsubproblems);
               }
            }
         }

         subproblemcount++;
         i++;
         if( i >= nsubproblems )
            i = 0;
         benders->lastchecked = i;

         /* checking whether the limits have been exceeded in the master problem */
         (*stopped) = SCIPisStopped(set->scip);
      }
   }

   (*nchecked) = subproblemcount;

   return SCIP_OKAY;
}

/** Calls the Benders' decompsition cuts for the given solve loop. There are four cases:
 *  i) solveloop == SCIP_BENDERSSOLVELOOP_CONVEX - only the LP Benders' cuts are called
 *  ii) solveloop == SCIP_BENDERSSOLVELOOP_CIP - only the CIP Benders' cuts are called
 *  iii) solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX - only the LP Benders' cuts are called
 *  iv) solveloop == SCIP_BENDERSSOLVELOOP_USERCIP - only the CIP Benders' cuts are called
 *
 *  The priority of the results are: SCIP_CONSADDED (SCIP_SEPARATED), SCIP_DIDNOTFIND, SCIP_FEASIBLE, SCIP_DIDNOTRUN. In
 *  this function, there are four levels of results that need to be assessed. These are:
 *  i) The result from the individual cut for the subproblem
 *  ii) The overall result for the subproblem from all cuts
 *  iii) the overall result for the solve loop from all cuts
 *  iv) the over all result from all solve loops.
 *  In each level, the priority of results must be adhered to.
 */
static
SCIP_RETCODE generateBendersCuts(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RESULT*          result,             /**< result of the pricing process */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_BENDERSSOLVELOOP solveloop,          /**< the current solve loop */
   SCIP_Bool             checkint,           /**< are the subproblems called during a check/enforce of integer sols? */
   int                   nchecked,           /**< the number of subproblems checked in this solve loop, they may not be solved */
   SCIP_Bool*            subprobsolved,      /**< an array indicating the subproblems that were solved in this loop. */
   SCIP_BENDERSSUBSTATUS* substatus,         /**< array to store the status of the subsystem */
   int**                 mergecands,         /**< the subproblems that are merge candidates */
   int*                  npriomergecands,    /**< the number of priority merge candidates. */
   int*                  nmergecands,        /**< the number of merge candidates. */
   int*                  nsolveloops         /**< the number of solve loops, is updated w.r.t added cuts */
   )
{
   SCIP_BENDERSCUT** benderscuts;
   SCIP_RESULT solveloopresult;
   int nbenderscuts;
   int nsubproblems;
   int subproblemcount;
   SCIP_Longint addedcuts = 0;
   int i;
   int j;
   SCIP_Bool onlyconvexcheck;

   assert(benders != NULL);
   assert(set != NULL);

   /* getting the Benders' decomposition cuts */
   benderscuts = SCIPbendersGetBenderscuts(benders);
   nbenderscuts = SCIPbendersGetNBenderscuts(benders);

   solveloopresult = SCIP_DIDNOTRUN;

   /* getting the number of subproblems in the Benders' decomposition */
   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* in the case of an LNS check, only the convex relaxations of the subproblems will be solved. This is a performance
    * feature, since solving the convex relaxation is typically much faster than solving the corresponding CIP. While
    * the CIP is not solved during the LNS check, the solutions are still of higher quality than when Benders' is not
    * employed.
    */
   onlyconvexcheck = SCIPbendersOnlyCheckConvexRelax(benders);

   /* It is only possible to add cuts to the problem if it has not already been solved */
   if( SCIPsetGetStage(set) < SCIP_STAGE_SOLVED && type != SCIP_BENDERSENFOTYPE_CHECK )
   {
      /* This is done in two loops. The first is by subproblem and the second is by cut type. */
      i = benders->firstchecked;
      subproblemcount = 0;
      while( subproblemcount < nchecked )
      {
         SCIP_RESULT subprobresult;
         SCIP_Bool convexsub = SCIPbendersSubproblemIsConvex(benders, i);

         /* cuts can only be generated if the subproblem is not independent and if it has been solved. The subproblem
          * solved flag is important for the user-defined subproblem solving methods
          */
         if( subproblemIsActive(benders, i) && subprobsolved[i] )
         {
            subprobresult = SCIP_DIDNOTRUN;
            for( j = 0; j < nbenderscuts; j++ )
            {
               SCIP_RESULT cutresult;
               SCIP_Longint prevaddedcuts;

               assert(benderscuts[j] != NULL);

               prevaddedcuts = SCIPbenderscutGetNFound(benderscuts[j]);
               cutresult = SCIP_DIDNOTRUN;

               /* the result is updated only if a Benders' cut is generated or one was not found. However, if a cut has
                * been found in a previous iteration, then the result is returned as SCIP_CONSADDED or SCIP_SEPARATED.
                * This result is permitted because if a constraint was added, the solution that caused the error in the cut
                * generation will be cutoff from the master problem.
                */
               if( (SCIPbenderscutIsLPCut(benderscuts[j]) && (solveloop == SCIP_BENDERSSOLVELOOP_CONVEX
                        || solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX))
                  || (!SCIPbenderscutIsLPCut(benderscuts[j]) && ((solveloop == SCIP_BENDERSSOLVELOOP_CIP && !convexsub)
                        || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP)) )
                  SCIP_CALL( SCIPbenderscutExec(benderscuts[j], set, benders, sol, i, type, &cutresult) );

               addedcuts += (SCIPbenderscutGetNFound(benderscuts[j]) - prevaddedcuts);

               /* the result is updated only if a Benders' cut is generated */
               if( cutresult == SCIP_CONSADDED || cutresult == SCIP_SEPARATED )
               {
                  subprobresult = cutresult;

                  benders->ncutsfound++;

                  /* at most a single cut is generated for each subproblem */
                  break;
               }
               else
               {
                  /* checking from lowest priority result */
                  if( subprobresult == SCIP_DIDNOTRUN )
                     subprobresult = cutresult;
                  else if( subprobresult == SCIP_FEASIBLE && cutresult == SCIP_DIDNOTFIND )
                     subprobresult = cutresult;
                  /* if the subprobresult is SCIP_DIDNOTFIND, then it can't be updated. */
               }
            }

            /* the highest priority for the results is CONSADDED and SEPARATED. The solveloopresult will always be
             * updated if the subprobresult is either of these.
             */
            if( subprobresult == SCIP_CONSADDED || subprobresult == SCIP_SEPARATED )
            {
               solveloopresult = subprobresult;
            }
            else if( subprobresult == SCIP_FEASIBLE )
            {
               /* updating the solve loop result based upon the priority */
               if( solveloopresult == SCIP_DIDNOTRUN )
                  solveloopresult = subprobresult;
            }
            else if( subprobresult == SCIP_DIDNOTFIND )
            {
               /* updating the solve loop result based upon the priority */
               if( solveloopresult == SCIP_DIDNOTRUN || solveloopresult == SCIP_FEASIBLE )
                  solveloopresult = subprobresult;

               /* since a cut was not found, then merging could be useful to avoid this in subsequent iterations. The
                * candidate is labelled as a non-priority merge candidate
                */
               if( substatus[i] != SCIP_BENDERSSUBSTATUS_OPTIMAL )
               {
                  (*mergecands)[(*nmergecands)] = i;
                  (*nmergecands)++;
               }
            }
            else if( subprobresult == SCIP_DIDNOTRUN )
            {
               /* if the subproblem is infeasible and no cut generation methods were run, then the infeasibility will
                * never be resolved. As such, the subproblem will be merged into the master problem. If the subproblem
                * was not infeasible, then it is added as a possible merge candidate
                */
               if( substatus[i] == SCIP_BENDERSSUBSTATUS_INFEAS )
               {
                  (*mergecands)[(*nmergecands)] = (*mergecands)[(*npriomergecands)];
                  (*mergecands)[(*npriomergecands)] = i;
                  (*npriomergecands)++;
                  (*nmergecands)++;
               }
               else if( substatus[i] != SCIP_BENDERSSUBSTATUS_OPTIMAL )
               {
                  (*mergecands)[(*nmergecands)] = i;
                  (*nmergecands)++;
               }
            }
         }

         subproblemcount++;
         i++;
         if( i >= nsubproblems )
            i = 0;
      }
   }

   /* updating the overall result based upon the priorities */
   if( solveloopresult == SCIP_CONSADDED || solveloopresult == SCIP_SEPARATED )
   {
      (*result) = solveloopresult;
   }
   else if( solveloopresult == SCIP_FEASIBLE )
   {
      /* updating the solve loop result based upon the priority */
      if( (*result) == SCIP_DIDNOTRUN )
         (*result) = solveloopresult;
   }
   else if( solveloopresult == SCIP_DIDNOTFIND )
   {
      /* updating the solve loop result based upon the priority */
      if( (*result) == SCIP_DIDNOTRUN || (*result) == SCIP_FEASIBLE )
         (*result) = solveloopresult;
   }

   /* if no cuts were added, then the number of solve loops is increased */
   if( addedcuts == 0 && SCIPbendersGetNConvexSubproblems(benders) < SCIPbendersGetNSubproblems(benders)
      && checkint && !onlyconvexcheck )
      (*nsolveloops) = 2;

   return SCIP_OKAY;
}

/** Solves the subproblem using the current master problem solution.
 *
 *  The checkint flag indicates whether integer feasibility can be assumed. If it is not assumed, i.e. checkint ==
 *  FALSE, then only the convex relaxations of the subproblems are solved. If integer feasibility is assumed, i.e.
 *  checkint == TRUE, then the convex relaxations and the full CIP are solved to generate Benders' cuts and check
 *  solution feasibility.
 *
 *  TODO: consider allowing the possibility to pass solution information back from the subproblems instead of the scip
 *  instance. This would allow the use of different solvers for the subproblems, more importantly allowing the use of an
 *  LP solver for LP subproblems.
 */
SCIP_RETCODE SCIPbendersExec(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RESULT*          result,             /**< result of the pricing process */
   SCIP_Bool*            infeasible,         /**< is the master problem infeasible with respect to the Benders' cuts? */
   SCIP_Bool*            auxviol,            /**< set to TRUE only if the solution is feasible but the aux vars are violated */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint            /**< should the integer solution be checked by the subproblems */
   )
{
   int nsubproblems;
   int subproblemcount;
   int nchecked;
   int nsolveloops;
   int nverified;
   int* mergecands;
   int npriomergecands;
   int nmergecands;
   SCIP_Bool* subprobsolved;
   SCIP_BENDERSSUBSTATUS* substatus;
   SCIP_Bool optimal;
   SCIP_Bool allverified;
   SCIP_Bool success;
   SCIP_Bool stopped;
   int i;
   int l;

   success = TRUE;
   stopped = FALSE;

   SCIPsetDebugMsg(set, "Starting Benders' decomposition subproblem solving. type %d checkint %d\n", type, checkint);

   /* start timing */
   SCIPclockStart(benders->bendersclock, set);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   (*auxviol) = FALSE;
   (*infeasible) = FALSE;

   /* It is assumed that the problem is optimal, until a subproblem is found not to be optimal. However, not all
    * subproblems could be checked in each iteration. As such, it is not possible to state that the problem is optimal
    * if not all subproblems are checked. Situations where this may occur is when a subproblem is a MIP and only the LP
    * is solved. Also, in a distributed computation, then it may be advantageous to only solve some subproblems before
    * resolving the master problem. As such, for a problem to be optimal, then (optimal && allverified) == TRUE
    */
   optimal = TRUE;
   nverified = 0;

   assert(benders != NULL);
   assert(result != NULL);
   assert(infeasible != NULL);
   assert(auxviol != NULL);

   /* if the Benders' decomposition is called from a sub-scip, it is assumed that this is an LNS heuristic. As such, the
    * check is not performed and the solution is assumed to be feasible
    */
   if( benders->iscopy
      && (!benders->lnscheck
         || (benders->lnsmaxdepth > -1 && SCIPgetDepth(benders->sourcescip) > benders->lnsmaxdepth)) )
   {
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* it is not necessary to check all primal solutions by solving the Benders' decomposition subproblems.
    * Only the improving solutions are checked to improve efficiency of the algorithm.
    * If the solution is non-improving, the result FEASIBLE is returned. While this may be incorrect w.r.t to the
    * Benders' subproblems, this solution will never be the optimal solution. A non-improving solution may be used
    * within LNS primal heuristics. If this occurs, the improving solution, if found, will be checked by the solving
    * the Benders' decomposition subproblems.
    * TODO: Add a parameter to control this behaviour.
    */
   if( checkint && SCIPsetIsFeasLE(set, SCIPgetPrimalbound(set->scip), SCIPgetSolOrigObj(set->scip, sol)) )
   {
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* if the enforcement type is SCIP_BENDERSENFOTYPE_LP and the LP is currently unbounded. This could mean that there
    * is no lower bound on the auxiliary variables. In this case, we try to update the lower bound for the auxiliary
    * variables.
    */
   if( type == SCIP_BENDERSENFOTYPE_LP && SCIPgetLPSolstat(set->scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY
      && benders->updateauxvarbound )
   {
      SCIP_CALL( updateAuxiliaryVarLowerbound(benders, set, result) );

      /* the auxiliary variable bound will only be updated once. */
      benders->updateauxvarbound = FALSE;
   }

   /* setting the first subproblem to check in this round of subproblem checks */
   benders->firstchecked = benders->lastchecked;

   /* sets the stored objective function values of the subproblems to infinity */
   resetSubproblemObjectiveValue(benders);

   *result = SCIP_DIDNOTRUN;

   if( benders->benderspresubsolve != NULL )
   {
      SCIP_Bool skipsolve;

      skipsolve = FALSE;
      SCIP_CALL( benders->benderspresubsolve(set->scip, benders, sol, type, checkint, &skipsolve, result) );

      /* evaluate result */
      if( (*result) != SCIP_DIDNOTRUN
         && (*result) != SCIP_FEASIBLE
         && (*result) != SCIP_INFEASIBLE
         && (*result) != SCIP_CONSADDED
         && (*result) != SCIP_SEPARATED )
      {
         SCIPerrorMessage("the user-defined pre subproblem solving method for the Benders' decomposition <%s> returned "
            "invalid result <%d>\n", benders->name, *result);
         return SCIP_INVALIDRESULT;
      }

      /* if the solve must be skipped, then the solving loop is exited and the user defined result is returned */
      if( skipsolve )
      {
         SCIPsetDebugMsg(set, "skipping the subproblem solving for Benders' decomposition <%s>. "
            "returning result <%d>\n", benders->name, *result);
         return SCIP_OKAY;
      }
   }

   /* allocating memory for the infeasible subproblem array */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(set->scip, &subprobsolved, nsubproblems) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(set->scip, &substatus, nsubproblems) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(set->scip, &mergecands, nsubproblems) );
   npriomergecands = 0;
   nmergecands = 0;

   /* by default the number of solve loops is 1. This is the case if all subproblems are LP or the user has defined a
    * benderssolvesub callback. If there is a subproblem that is not an LP, then 2 solve loops are performed. The first
    * loop is the LP solving loop, the second solves the subproblem to integer optimality.
    */
   nsolveloops = 1;

   for( l = 0; l < nsolveloops; l++ )
   {
      SCIP_BENDERSSOLVELOOP solveloop;    /* identifies what problem type is solve in this solve loop */

      /* if either benderssolvesubconvex or benderssolvesub are implemented, then the user callbacks are invoked */
      if( benders->benderssolvesubconvex != NULL || benders->benderssolvesub != NULL )
      {
         if( l == 0 )
            solveloop = SCIP_BENDERSSOLVELOOP_USERCONVEX;
         else
            solveloop = SCIP_BENDERSSOLVELOOP_USERCIP;
      }
      else
         solveloop = (SCIP_BENDERSSOLVELOOP) l;

      /* solving the subproblems for this round of enforcement/checking. */
      SCIP_CALL( solveBendersSubproblems(benders, set, sol, type, solveloop, checkint, &nchecked, &nverified,
            &subprobsolved, &substatus, infeasible, &optimal, &stopped) );

      /* if the solving has been stopped, then the subproblem solving and cut generation must terminate */
      if( stopped )
         goto TERMINATE;

      /* Generating cuts for the subproblems. Cuts are only generated when the solution is from primal heuristics,
       * relaxations or the LP
       */
      if( type != SCIP_BENDERSENFOTYPE_PSEUDO )
      {
         SCIP_CALL( generateBendersCuts(benders, set, sol, result, type, solveloop, checkint, nchecked,
            subprobsolved, substatus, &mergecands, &npriomergecands, &nmergecands, &nsolveloops) );
      }
      else
      {
         /* The first solving loop solves the convex subproblems and the convex relaxations of the CIP subproblems. The
          * second solving loop solves the CIP subproblems. The second solving loop is only called if the integer
          * feasibility is being checked and if the convex subproblems and convex relaxations are not infeasible.
          */
         if( !(*infeasible) && checkint && !SCIPbendersOnlyCheckConvexRelax(benders)
            && SCIPbendersGetNConvexSubproblems(benders) < SCIPbendersGetNSubproblems(benders))
            nsolveloops = 2;
      }
   }

   allverified = (nverified == nsubproblems);

   SCIPsetDebugMsg(set, "End Benders' decomposition subproblem solve. result %d infeasible %d auxviol %d nverified %d\n",
      *result, *infeasible, *auxviol, nverified);

#ifdef SCIP_DEBUG
   if( (*result) == SCIP_CONSADDED )
   {
      SCIPsetDebugMsg(set, "Benders' decomposition: Cut added\n");
   }
#endif

   /* if the number of checked pseudo solutions exceeds a set limit, then all subproblems are passed as merge
    * candidates. Currently, merging subproblems into the master problem is the only method for resolving numerical
    * troubles.
    *
    * We are only interested in the pseudo solutions that have been checked completely for integrality. This is
    * identified by checkint == TRUE. This means that the Benders' decomposition constraint is one of the last
    * constraint handlers that must resolve the infeasibility. If the Benders' decomposition framework can't resolve the
    * infeasibility, then this will result in an error.
    */
   if( type == SCIP_BENDERSENFOTYPE_PSEUDO && checkint )
   {
      benders->npseudosols++;

      if( benders->npseudosols > BENDERS_MAXPSEUDOSOLS )
      {
         /* if a priority merge candidate already exists, then no other merge candidates need to be added.*/
         if( npriomergecands == 0 )
         {
            /* all subproblems are added to the merge candidate list. The first active subproblem is added as a
             * priority merge candidate
             */
            nmergecands = 0;
            npriomergecands = 1;
            for( i = 0; i < nsubproblems; i++ )
            {
               /* only active subproblems are added to the merge candidate list */
               if( subproblemIsActive(benders, i) )
               {
                  mergecands[nmergecands] = i;
                  nmergecands++;
               }
            }

            SCIPverbMessage(set->scip, SCIP_VERBLEVEL_HIGH, NULL, "   The number of checked pseudo solutions exceeds the "
              "limit of %d. All active subproblems are merge candidates, with subproblem %d a priority candidate.\n",
              BENDERS_MAXPSEUDOSOLS, mergecands[0]);
         }
      }
   }
   else
      benders->npseudosols = 0;

   /* if the result is SCIP_DIDNOTFIND, then there was a error in generating cuts in all subproblems that are not
    * optimal. This result does not cutoff any solution, so the Benders' decomposition algorithm will fail.
    * TODO: Work out a way to ensure Benders' decomposition does not terminate due to a SCIP_DIDNOTFIND result.
    */
   if( (*result) == SCIP_DIDNOTFIND )
   {
      if( type == SCIP_BENDERSENFOTYPE_PSEUDO )
         (*result) = SCIP_SOLVELP;
      else
         (*result) = SCIP_INFEASIBLE;

      SCIPerrorMessage("An error was found when generating cuts for non-optimal subproblems of Benders' "
         "decomposition <%s>. Consider merging the infeasible subproblems into the master problem.\n", SCIPbendersGetName(benders));

      /* since no other cuts are generated, then this error will result in a crash. It is possible to avoid the error,
       * by merging the affected subproblem into the master problem.
       */
      success = FALSE;

      goto POSTSOLVE;
   }

   if( type == SCIP_BENDERSENFOTYPE_PSEUDO )
   {
      if( (*infeasible) || !allverified )
         (*result) = SCIP_SOLVELP;
      else
      {
         (*result) = SCIP_FEASIBLE;

         /* if the subproblems are not infeasible, but they are also not optimal. This means that there is a violation
          * in the auxiliary variable values. In this case, a feasible result is returned with the auxviol flag set to
          * TRUE.
          */
         (*auxviol) = !optimal;
      }
   }
   else if( checkint && (type == SCIP_BENDERSENFOTYPE_CHECK || (*result) != SCIP_CONSADDED) )
   {
      /* if the subproblems are being solved as part of conscheck, then the results flag must be returned after the solving
       * has completed.
       */
      if( (*infeasible) || !allverified )
         (*result) = SCIP_INFEASIBLE;
      else
      {
         (*result) = SCIP_FEASIBLE;

         /* if the subproblems are not infeasible, but they are also not optimal. This means that there is a violation
          * in the auxiliary variable values. In this case, a feasible result is returned with the auxviol flag set to
          * TRUE.
          */
         (*auxviol) = !optimal;
      }
   }

POSTSOLVE:
   /* calling the post-solve call back for the Benders' decomposition algorithm. This allows the user to work directly
    * with the solved subproblems and the master problem */
   if( benders->benderspostsolve != NULL )
   {
      SCIP_Bool merged;

      merged = FALSE;

      SCIP_CALL( benders->benderspostsolve(set->scip, benders, sol, type, mergecands, npriomergecands, nmergecands,
            checkint, (*infeasible), &merged) );

      if( merged )
      {
         (*result) = SCIP_CONSADDED;

         /* since subproblems have been merged, then constraints have been added. This could resolve the unresolved
          * infeasibility, so the error has been corrected.
          */
         success = TRUE;
      }
      else if( !success )
      {
         SCIPerrorMessage("An error occurred during Benders' decomposition cut generations and no merging had been "
            "performed. It is not possible to continue solving the problem by Benders' decomposition\n");
      }
   }

TERMINATE:
   /* freeing the subproblems after the cuts are generated */
   i = benders->firstchecked;
   subproblemcount = 0;

   /* if the solving process has stopped, then all subproblems need to be freed */
   if( stopped )
      nchecked = nsubproblems;

   while( subproblemcount < nchecked )
   {
      SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, i) );

      subproblemcount++;
      i++;
      if( i >= nsubproblems )
         i = 0;
   }

#ifndef NDEBUG
   for( i = 0; i < nsubproblems; i++ )
      assert(SCIPgetStage(SCIPbendersSubproblem(benders, i)) < SCIP_STAGE_TRANSFORMED
         || !SCIPinProbing(SCIPbendersSubproblem(benders, i)) || !subproblemIsActive(benders, i));
#endif

   /* increment the number of calls to the Benders' decomposition subproblem solve */
   benders->ncalls++;

   SCIPsetDebugMsg(set, "End Benders' decomposition execution method. result %d infeasible %d auxviol %d\n", *result,
      *infeasible, *auxviol);

   /* end timing */
   SCIPclockStop(benders->bendersclock, set);

   /* freeing memory */
   SCIPfreeBlockMemoryArray(set->scip, &mergecands, nsubproblems);
   SCIPfreeBlockMemoryArray(set->scip, &substatus, nsubproblems);
   SCIPfreeBlockMemoryArray(set->scip, &subprobsolved, nsubproblems);

   if( !success )
      return SCIP_ERROR;
   else
      return SCIP_OKAY;
}

/** solves the user-defined subproblem solving function */
static
SCIP_RETCODE executeUserDefinedSolvesub(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the subproblem number */
   SCIP_BENDERSSOLVELOOP solveloop,          /**< the solve loop iteration. The first iter is for LP, the second for IP */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_Real*            objective,          /**< the objective function value of the subproblem */
   SCIP_RESULT*          result              /**< the result from solving the subproblem */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);
   assert(benders->benderssolvesubconvex != NULL || benders->benderssolvesub != NULL);

   assert(solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP);

   (*objective) = -SCIPsetInfinity(set);

   /* calls the user defined subproblem solving method. Only the convex relaxations are solved during the Large
    * Neighbourhood Benders' Search. */
   if( solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX )
   {
      if( benders->benderssolvesubconvex != NULL )
      {
         SCIP_CALL( benders->benderssolvesubconvex(set->scip, benders, sol, probnumber,
               SCIPbendersOnlyCheckConvexRelax(benders), objective, result) );
      }
      else
         (*result) = SCIP_DIDNOTRUN;
   }
   else if( solveloop == SCIP_BENDERSSOLVELOOP_USERCIP )
   {
      if( benders->benderssolvesub != NULL )
      {
         SCIP_CALL( benders->benderssolvesub(set->scip, benders, sol, probnumber, objective, result) );
      }
      else
         (*result) = SCIP_DIDNOTRUN;
   }

   /* evaluate result */
   if( (*result) != SCIP_DIDNOTRUN
      && (*result) != SCIP_FEASIBLE
      && (*result) != SCIP_INFEASIBLE
      && (*result) != SCIP_UNBOUNDED )
   {
      SCIPerrorMessage("the user-defined solving method for the Benders' decomposition <%s> returned invalid result <%d>\n",
         benders->name, *result);
      return SCIP_INVALIDRESULT;
   }

   if( (*result) == SCIP_INFEASIBLE )
      (*infeasible) = TRUE;

   if( (*result) == SCIP_FEASIBLE
      && (SCIPsetIsInfinity(set, -(*objective)) || SCIPsetIsInfinity(set, (*objective))) )
   {
      SCIPerrorMessage("the user-defined solving method for the Benders' decomposition <%s> returned objective value %g\n",
         benders->name, (*objective));
      return SCIP_ERROR;
   }

   /* if the result is SCIP_DIDNOTFIND, then an error is returned and SCIP will terminate. */
   if( (*result) == SCIP_DIDNOTFIND )
      return SCIP_ERROR;
   else
      return SCIP_OKAY;
}

/** executes the subproblem solving process */
SCIP_RETCODE SCIPbendersExecSubproblemSolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the subproblem number */
   SCIP_BENDERSSOLVELOOP solveloop,          /**< the solve loop iteration. The first iter is for LP, the second for IP */
   SCIP_Bool             enhancement,        /**< is the solve performed as part of and enhancement? */
   SCIP_Bool*            solved,             /**< flag to indicate whether the subproblem was solved */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type                /**< the enforcement type calling this function */
   )
{
   SCIP* subproblem;
   SCIP_SOL* bestsol;
   SCIP_RESULT result;
   SCIP_Real objective;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   SCIPsetDebugMsg(set, "Benders' decomposition: solving subproblem %d\n", probnumber);

   result = SCIP_DIDNOTRUN;
   objective = SCIPsetInfinity(set);

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* initially setting the solved flag to FALSE */
   (*solved) = FALSE;

   /* if the subproblem solve callback is implemented, then that is used instead of the default setup */
   if( solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP )
   {
      /* calls the user defined subproblem solving method. Only the convex relaxations are solved during the Large
       * Neighbourhood Benders' Search. */
      SCIP_CALL( executeUserDefinedSolvesub(benders, set, sol, probnumber, solveloop, infeasible, &objective, &result) );

      /* if the result is DIDNOTRUN, then the subproblem was not solved */
      (*solved) = (result != SCIP_DIDNOTRUN);
   }
   else
   {
      /* setting up the subproblem */
      if( solveloop == SCIP_BENDERSSOLVELOOP_CONVEX )
      {
         SCIP_CALL( SCIPbendersSetupSubproblem(benders, set, sol, probnumber) );

         /* if the limits of the master problem were hit during the setup process, then the subproblem will not have
          * been setup. In this case, the solving function must be exited.
          */
         if( !SCIPbendersSubproblemIsSetup(benders, probnumber) )
         {
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPsetInfinity(set));
            (*solved) = FALSE;
            return SCIP_OKAY;
         }
      }
      else
      {
         SCIP_CALL( updateEventhdlrUpperbound(benders, probnumber, SCIPbendersGetAuxiliaryVarVal(benders, set, sol, probnumber)) );
      }

      /* solving the subproblem
       * the LP of the subproblem is solved in the first solveloop.
       * In the second solve loop, the MIP problem is solved */
      if( solveloop == SCIP_BENDERSSOLVELOOP_CONVEX || SCIPbendersSubproblemIsConvex(benders, probnumber) )
      {
         SCIP_CALL( SCIPbendersSolveSubproblemLP(set->scip, benders, probnumber, infeasible) );

         /* if the LP was solved without error, then the subproblem is labelled as solved */
         if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_OPTIMAL
            || SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE )
            (*solved) = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPbendersSolveSubproblemCIP(set->scip, benders, probnumber, infeasible, type, FALSE) );

         /* if the generic subproblem solving methods are used, then the CIP subproblems are always solved. */
         (*solved) = TRUE;
      }
   }

   bestsol = SCIPgetBestSol(subproblem);

   if( !enhancement )
   {
      /* The following handles the cases when the subproblem is OPTIMAL, INFEASIBLE and UNBOUNDED.
       * If a subproblem is unbounded, then the auxiliary variables are set to -infinity and the unbounded flag is
       * returned as TRUE. No cut will be generated, but the result will be set to SCIP_FEASIBLE.
       */
      if( solveloop == SCIP_BENDERSSOLVELOOP_CONVEX )
      {
         if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_OPTIMAL )
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPgetSolOrigObj(subproblem, NULL));
         else if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE )
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPsetInfinity(set));
         else if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
         {
            SCIPerrorMessage("The LP of Benders' decomposition subproblem %d is unbounded. This should not happen.\n",
               probnumber);
            SCIPABORT();
         }
         else if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_ERROR
            || SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_NOTSOLVED
            || SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_TIMELIMIT )
         {
            SCIPverbMessage(set->scip, SCIP_VERBLEVEL_FULL, NULL, "   Benders' decomposition: Error solving LP "
               "relaxation of subproblem %d. No cut will be generated for this subproblem.\n", probnumber);
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPsetInfinity(set));
         }
         else
         {
            SCIPerrorMessage("Invalid status returned from solving the LP of Benders' decomposition subproblem %d. LP status: %d\n",
               probnumber, SCIPgetLPSolstat(subproblem));
            SCIPABORT();
         }
      }
      else if( solveloop == SCIP_BENDERSSOLVELOOP_CIP )
      {
         /* TODO: Consider whether other solutions status should be handled */
         if( SCIPgetStatus(subproblem) == SCIP_STATUS_OPTIMAL )
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPgetSolOrigObj(subproblem, bestsol));
         else if( SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE )
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPsetInfinity(set));
         else if( SCIPgetStatus(subproblem) == SCIP_STATUS_USERINTERRUPT || SCIPgetStatus(subproblem) == SCIP_STATUS_BESTSOLLIMIT )
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPgetSolOrigObj(subproblem, bestsol));
         else if( SCIPgetStatus(subproblem) == SCIP_STATUS_MEMLIMIT
            || SCIPgetStatus(subproblem) == SCIP_STATUS_TIMELIMIT )
         {
            SCIPverbMessage(set->scip, SCIP_VERBLEVEL_FULL, NULL, "   Benders' decomposition: Error solving CIP "
               "of subproblem %d. No cut will be generated for this subproblem.\n", probnumber);
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPsetInfinity(set));
         }
         else if( SCIPgetStatus(subproblem) == SCIP_STATUS_UNBOUNDED )
         {
            SCIPerrorMessage("The Benders' decomposition subproblem %d is unbounded. This should not happen.\n",
               probnumber);
            SCIPABORT();
         }
         else
         {
            SCIPerrorMessage("Invalid status returned from solving Benders' decomposition subproblem %d. Solution status: %d\n",
               probnumber, SCIPgetStatus(subproblem));
            SCIPABORT();
         }
      }
      else
      {
         assert(solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP);
         if( result == SCIP_FEASIBLE )
            SCIPbendersSetSubproblemObjval(benders, probnumber, objective);
         else if( result == SCIP_INFEASIBLE )
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPsetInfinity(set));
         else if( result == SCIP_UNBOUNDED )
         {
            SCIPerrorMessage("The Benders' decomposition subproblem %d is unbounded. This should not happen.\n",
               probnumber);
            SCIPABORT();
         }
         else if( result != SCIP_DIDNOTRUN )
         {
            SCIPerrorMessage("Invalid result <%d> from user-defined subproblem solving method. This should not happen.\n",
               result);
         }
      }
   }

   return SCIP_OKAY;
}

/** sets up the subproblem using the solution to the master problem  */
SCIP_RETCODE SCIPbendersSetupSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_VAR** vars;
   SCIP_VAR* mastervar;
   SCIP_Real solval;
   int nvars;
   int i;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* changing all of the master problem variable to continuous. */
   SCIP_CALL( SCIPbendersChgMastervarsToCont(benders, set, probnumber) );

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* if the Benders' decomposition subproblem is an LP, then probing mode must be started.
    * If the subproblem is a MIP, the problem must be initialised, put into SCIP_STAGE_SOLVING to be able to change the
    * variable bounds. The probing mode is entered once the variable bounds are set.
    * In the MIP case, the transformed problem is freed after each subproblem solve round. */
   if( SCIPbendersSubproblemIsConvex(benders, probnumber) )
   {
      SCIP_CALL( SCIPstartProbing(subproblem) );
   }
   else
   {
      SCIP_Bool success;

      SCIP_CALL( initialiseSubproblem(benders, set, probnumber, &success) );

      if( !success )
      {
         /* set the flag to indicate that the subproblems have been set up */
         SCIPbendersSetSubproblemIsSetup(benders, probnumber, FALSE);

         return SCIP_OKAY;
      }
   }

   vars = SCIPgetVars(subproblem);
   nvars = SCIPgetNVars(subproblem);

   /* looping over all variables in the subproblem to find those corresponding to the master problem variables. */
   /* TODO: It should be possible to store the pointers to the master variables to speed up the subproblem setup */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPbendersGetVar(benders, set, vars[i], &mastervar, -1) );

      if( mastervar != NULL )
      {
         /* It is possible due to numerics that the solution value exceeds the upper or lower bounds. When this
          * happens, it causes an error in the LP solver as a result of inconsistent bounds. So the following statements
          * are used to ensure that the bounds are not exceeded when applying the fixings for the Benders'
          * decomposition subproblems
          */
         solval = SCIPgetSolVal(set->scip, sol, mastervar);
         if( solval > SCIPvarGetUbLocal(vars[i]) )
            solval = SCIPvarGetUbLocal(vars[i]);
         else if( solval < SCIPvarGetLbLocal(vars[i]) )
            solval = SCIPvarGetLbLocal(vars[i]);

         /* fixing the variable in the subproblem */
         if( !SCIPisEQ(subproblem, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) )
         {
            if( SCIPisGT(subproblem, solval, SCIPvarGetLbLocal(vars[i])) )
            {
               SCIP_CALL( SCIPchgVarLb(subproblem, vars[i], solval) );
            }
            if( SCIPisLT(subproblem, solval, SCIPvarGetUbLocal(vars[i])) )
            {
               SCIP_CALL( SCIPchgVarUb(subproblem, vars[i], solval) );
            }
         }

         assert(SCIPisEQ(subproblem, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])));
      }
   }

   /* if the subproblem is a MIP, the probing mode is entered after setting up the subproblem */
   if( !SCIPbendersSubproblemIsConvex(benders, probnumber) )
   {
      SCIP_CALL( SCIPstartProbing(subproblem) );
   }

   /* set the flag to indicate that the subproblems have been set up */
   SCIPbendersSetSubproblemIsSetup(benders, probnumber, TRUE);

   return SCIP_OKAY;
}

/** Solve a Benders' decomposition subproblems. This will either call the user defined method or the generic solving
 *  methods. If the generic method is called, then the subproblem must be set up before calling this method. */
SCIP_RETCODE SCIPbendersSolveSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_Bool             solvecip,           /**< directly solve the CIP subproblem */
   SCIP_Real*            objective           /**< the objective function value of the subproblem, can be NULL */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* the subproblem must be set up before this function is called. */
   if( !SCIPbendersSubproblemIsSetup(benders, probnumber) && !SCIPbendersSubproblemIsIndependent(benders, probnumber) )
   {
      SCIPerrorMessage("Benders' decomposition subproblem %d must be set up before calling SCIPbendersSolveSubproblem(). Call SCIPsetupSubproblem() first.\n", probnumber);
      return SCIP_ERROR;
   }

   /* if the subproblem solve callback is implemented, then that is used instead of the default setup */
   if( benders->benderssolvesubconvex != NULL ||  benders->benderssolvesub != NULL)
   {
      SCIP_BENDERSSOLVELOOP solveloop;
      SCIP_RESULT result;
      SCIP_Real subobj;

      if( solvecip )
         solveloop = SCIP_BENDERSSOLVELOOP_USERCIP;
      else
         solveloop = SCIP_BENDERSSOLVELOOP_USERCONVEX;

      SCIP_CALL( executeUserDefinedSolvesub(benders, set, sol, probnumber, solveloop, infeasible, &subobj, &result) );

      if( objective != NULL )
         (*objective) = subobj;
   }
   else
   {
      SCIP* subproblem;

      subproblem = SCIPbendersSubproblem(benders, probnumber);

      /* solving the subproblem */
      if( solvecip && !SCIPbendersSubproblemIsConvex(benders, probnumber) )
      {
         SCIP_CALL( SCIPbendersSolveSubproblemCIP(set->scip, benders, probnumber, infeasible, type, solvecip) );

         if( objective != NULL )
            (*objective) = SCIPgetSolOrigObj(subproblem, SCIPgetBestSol(subproblem));
      }
      else
      {
         SCIP_Bool success;

         /* if the subproblem is an LP, then it should have been initialised and in SCIP_STAGE_SOLVING.
          * in this case, the subproblem only needs to be put into probing mode. */
         if( SCIPbendersSubproblemIsConvex(benders, probnumber) )
         {
            /* if the subproblem is not in probing mode, then it must be put into that mode for the LP solve. */
            if( !SCIPinProbing(subproblem) )
            {
               SCIP_CALL( SCIPstartProbing(subproblem) );
            }

            success = TRUE;
         }
         else
         {
            SCIP_CALL( initialiseSubproblem(benders, set, probnumber, &success) );
         }

         /* if setting up the subproblem was successful */
         if( success )
         {
            SCIP_CALL( SCIPbendersSolveSubproblemLP(set->scip, benders, probnumber, infeasible) );

            if( objective != NULL )
               (*objective) = SCIPgetSolOrigObj(subproblem, NULL);
         }
         else
         {
            if( objective != NULL )
               (*objective) = SCIPinfinity(subproblem);
         }
      }
   }

   return SCIP_OKAY;
}

/** copies the time and memory limit from the master problem to the subproblem */
static
SCIP_RETCODE copyMemoryAndTimeLimits(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP*                 subproblem          /**< the Benders' decomposition subproblem */
   )
{
   SCIP_Real mastertimelimit;
   SCIP_Real subtimelimit;
   SCIP_Real maxsubtimelimit;
   SCIP_Real mastermemorylimit;
   SCIP_Real submemorylimit;
   SCIP_Real maxsubmemorylimit;

   assert(scip != NULL);

   /* setting the time limit for the Benders' decomposition subproblems. It is set to 102% of the remaining time. */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &mastertimelimit) );
   maxsubtimelimit = SCIPparamGetRealMax(SCIPgetParam(subproblem, "limits/time"));
   subtimelimit = (mastertimelimit - SCIPgetSolvingTime(scip)) * 1.02;
   subtimelimit = MIN(subtimelimit, maxsubtimelimit);
   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/time", MAX(0.0, subtimelimit)) );

   /* setting the memory limit for the Benders' decomposition subproblems. */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &mastermemorylimit) );
   maxsubmemorylimit = SCIPparamGetRealMax(SCIPgetParam(subproblem, "limits/memory"));
   submemorylimit = mastermemorylimit - (SCIPgetMemUsed(scip) + SCIPgetMemExternEstim(scip))/1048576.0;
   submemorylimit = MIN(submemorylimit, maxsubmemorylimit);
   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/memory", MAX(0.0, submemorylimit)) );

   return SCIP_OKAY;
}

/** stores the original parameters from the subproblem */
static
SCIP_RETCODE storeOrigSubproblemParams(
   SCIP*                 subproblem,         /**< the SCIP data structure */
   SCIP_SUBPROBPARAMS*   origparams          /**< the original subproblem parameters */
   )
{
   assert(subproblem != NULL);
   assert(origparams != NULL);

   SCIP_CALL( SCIPgetRealParam(subproblem, "limits/memory", &origparams->limits_memory) );
   SCIP_CALL( SCIPgetRealParam(subproblem, "limits/time", &origparams->limits_time) );
   SCIP_CALL( SCIPgetBoolParam(subproblem, "conflict/enable", &origparams->conflict_enable) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/disablecutoff", &origparams->lp_disablecutoff) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/scaling", &origparams->lp_scaling) );
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/initalgorithm", &origparams->lp_initalg) );
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/resolvealgorithm", &origparams->lp_resolvealg) );
   SCIP_CALL( SCIPgetBoolParam(subproblem, "lp/alwaysgetduals", &origparams->lp_alwaysgetduals) );
   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/scaleobj", &origparams->misc_scaleobj) );
   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/catchctrlc", &origparams->misc_catchctrlc) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxrounds", &origparams->prop_maxrounds) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxroundsroot", &origparams->prop_maxroundsroot) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "constraints/linear/propfreq", &origparams->cons_linear_propfreq) );

   return SCIP_OKAY;
}

/** sets the parameters for the subproblem */
static
SCIP_RETCODE setSubproblemParams(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP*                 subproblem          /**< the subproblem SCIP instance */
   )
{
   assert(scip != NULL);
   assert(subproblem != NULL);

   /* copying memory and time limits */
   SCIP_CALL( copyMemoryAndTimeLimits(scip, subproblem) );

   /* Do we have to disable presolving? If yes, we have to store all presolving parameters. */
   SCIP_CALL( SCIPsetPresolving(subproblem, SCIP_PARAMSETTING_OFF, TRUE) );

   /* Disabling heuristics so that the problem is not trivially solved */
   SCIP_CALL( SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_OFF, TRUE) );

   /* store parameters that are changed for the generation of the subproblem cuts */
   SCIP_CALL( SCIPsetParam(subproblem, "conflict/enable", FALSE) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/disablecutoff", 1) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/scaling", 0) );

   SCIP_CALL( SCIPsetCharParam(subproblem, "lp/initalgorithm", 'd') );
   SCIP_CALL( SCIPsetCharParam(subproblem, "lp/resolvealgorithm", 'd') );

   SCIP_CALL( SCIPsetBoolParam(subproblem, "lp/alwaysgetduals", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/scaleobj", FALSE) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/catchctrlc", FALSE) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxroundsroot", 0) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "constraints/linear/propfreq", -1) );

   return SCIP_OKAY;
}

/** resets the original parameters from the subproblem */
static
SCIP_RETCODE resetOrigSubproblemParams(
   SCIP*                 subproblem,               /**< the SCIP data structure */
   SCIP_SUBPROBPARAMS*   origparams          /**< the original subproblem parameters */
   )
{
   assert(subproblem != NULL);
   assert(origparams != NULL);

   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/memory", origparams->limits_memory) );
   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/time", origparams->limits_time) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "conflict/enable", origparams->conflict_enable) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/disablecutoff", origparams->lp_disablecutoff) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/scaling", origparams->lp_scaling) );
   SCIP_CALL( SCIPsetCharParam(subproblem, "lp/initalgorithm", origparams->lp_initalg) );
   SCIP_CALL( SCIPsetCharParam(subproblem, "lp/resolvealgorithm", origparams->lp_resolvealg) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "lp/alwaysgetduals", origparams->lp_alwaysgetduals) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/scaleobj", origparams->misc_scaleobj) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/catchctrlc", origparams->misc_catchctrlc) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxrounds", origparams->prop_maxrounds) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxroundsroot", origparams->prop_maxroundsroot) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "constraints/linear/propfreq", origparams->cons_linear_propfreq) );

   return SCIP_OKAY;
}

/** solves the LP of the Benders' decomposition subproblem
 *
 *  This requires that the subproblem is in probing mode.
 */
SCIP_RETCODE SCIPbendersSolveSubproblemLP(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible          /**< a flag to indicate whether all subproblems are feasible */
   )
{
   SCIP* subproblem;
   SCIP_SUBPROBPARAMS* origparams;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;

   assert(benders != NULL);
   assert(infeasible != NULL);
   assert(SCIPbendersSubproblemIsSetup(benders, probnumber));

   (*infeasible) = FALSE;

   /* TODO: This should be solved just as an LP, so as a MIP. There is too much overhead with the MIP.
    * Need to change status check for checking the LP. */
   subproblem = SCIPbendersSubproblem(benders, probnumber);

   assert(SCIPisLPConstructed(subproblem));
   assert(SCIPinProbing(subproblem));

   /* allocating memory for the parameter storage */
   SCIP_CALL( SCIPallocBlockMemory(subproblem, &origparams) );

   /* store the original parameters of the subproblem */
   SCIP_CALL( storeOrigSubproblemParams(subproblem, origparams) );

   /* setting the subproblem parameters */
   SCIP_CALL( setSubproblemParams(scip, subproblem) );

   SCIP_CALL( SCIPsolveProbingLP(subproblem, -1, &lperror, &cutoff) );

   if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE )
      (*infeasible) = TRUE;
   else if( SCIPgetLPSolstat(subproblem) != SCIP_LPSOLSTAT_OPTIMAL
      && SCIPgetLPSolstat(subproblem) != SCIP_LPSOLSTAT_UNBOUNDEDRAY
      && SCIPgetLPSolstat(subproblem) != SCIP_LPSOLSTAT_NOTSOLVED
      && SCIPgetLPSolstat(subproblem) != SCIP_LPSOLSTAT_ERROR
      && SCIPgetLPSolstat(subproblem) != SCIP_LPSOLSTAT_TIMELIMIT )
   {
      SCIPerrorMessage("Invalid status: %d. Solving the LP relaxation of Benders' decomposition subproblem %d.\n",
         SCIPgetLPSolstat(subproblem), probnumber);
      SCIPABORT();
   }

   /* resetting the subproblem parameters */
   SCIP_CALL( resetOrigSubproblemParams(subproblem, origparams) );

   /* freeing the parameter storage */
   SCIPfreeBlockMemory(subproblem, &origparams);

   return SCIP_OKAY;
}

/** solves the Benders' decomposition subproblem */
SCIP_RETCODE SCIPbendersSolveSubproblemCIP(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_Bool             solvecip            /**< directly solve the CIP subproblem */
   )
{
   SCIP* subproblem;
   SCIP_SUBPROBPARAMS* origparams;

   assert(benders != NULL);
   assert(infeasible != NULL);

   (*infeasible) = FALSE;

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* allocating memory for the parameter storage */
   SCIP_CALL( SCIPallocBlockMemory(subproblem, &origparams) );

   /* store the original parameters of the subproblem */
   SCIP_CALL( storeOrigSubproblemParams(subproblem, origparams) );

   /* If the solve has been stopped for the subproblem, then we need to restart it to complete the solve. The subproblem
    * is stopped when it is a MIP so that LP cuts and IP cuts can be generated. */
   if( SCIPgetStage(subproblem) == SCIP_STAGE_SOLVING )
   {
      /* the subproblem should be in probing mode. Otherwise, the event handler did not work correctly */
      assert( SCIPinProbing(subproblem) );

      /* the probing mode needs to be stopped so that the MIP can be solved */
      SCIP_CALL( SCIPendProbing(subproblem) );

      /* the problem was interrupted in the event handler, so SCIP needs to be informed that the problem is to be restarted */
      SCIP_CALL( SCIPrestartSolve(subproblem) );

      /* if the solve type is for CHECK, then the FEASIBILITY emphasis setting is used. */
      if( type == SCIP_BENDERSENFOTYPE_CHECK )
      {
         SCIP_CALL( SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_FAST, TRUE) );

         /* the number of solution improvements is limited to try and prove feasibility quickly */
         /* NOTE: This should be a parameter */
         /* SCIP_CALL( SCIPsetIntParam(subproblem, "limits/bestsol", 5) ); */
      }
   }
   else if( solvecip )
   {
      /* if the MIP will be solved directly, then the probing mode needs to be skipped.
       * This is achieved by setting the solvecip flag in the event handler data to TRUE
       */
      SCIP_EVENTHDLR* eventhdlr;
      SCIP_EVENTHDLRDATA* eventhdlrdata;

      eventhdlr = SCIPfindEventhdlr(subproblem, MIPNODEFOCUS_EVENTHDLR_NAME);
      eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

      eventhdlrdata->solvecip = TRUE;
   }
   else
   {
      /* if the problem is not in probing mode, then we need to solve the LP. That requires all methods that will
       * modify the structure of the problem need to be deactivated */

      /* setting the subproblem parameters */
      SCIP_CALL( setSubproblemParams(scip, subproblem) );

#ifdef SCIP_MOREDEBUG
      SCIP_CALL( SCIPsetBoolParam(subproblem, "display/lpinfo", TRUE) );
#endif
   }

#ifdef SCIP_MOREDEBUG
      SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_FULL) );
#endif

   SCIP_CALL( SCIPsolve(subproblem) );

   if( SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE )
      (*infeasible) = TRUE;
   else if( SCIPgetStatus(subproblem) != SCIP_STATUS_OPTIMAL && SCIPgetStatus(subproblem) != SCIP_STATUS_UNBOUNDED
      && SCIPgetStatus(subproblem) != SCIP_STATUS_USERINTERRUPT && SCIPgetStatus(subproblem) != SCIP_STATUS_BESTSOLLIMIT
      && SCIPgetStatus(subproblem) != SCIP_STATUS_TIMELIMIT && SCIPgetStatus(subproblem) != SCIP_STATUS_MEMLIMIT )
   {
      SCIPerrorMessage("Invalid status: %d. Solving the CIP of Benders' decomposition subproblem %d.\n",
         SCIPgetStatus(subproblem), probnumber);
      SCIPABORT();
   }

   /* resetting the subproblem parameters */
   SCIP_CALL( resetOrigSubproblemParams(subproblem, origparams) );

   /* freeing the parameter storage */
   SCIPfreeBlockMemory(subproblem, &origparams);

   return SCIP_OKAY;
}

/** frees the subproblems */
SCIP_RETCODE SCIPbendersFreeSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(benders->bendersfreesub != NULL
      || (benders->bendersfreesub == NULL && benders->benderssolvesubconvex == NULL && benders->benderssolvesub == NULL));
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   if( benders->bendersfreesub != NULL )
   {
      SCIP_CALL( benders->bendersfreesub(set->scip, benders, probnumber) );
   }
   else
   {
      /* the subproblem is only freed if it is not independent */
      if( subproblemIsActive(benders, probnumber) )
      {
         SCIP* subproblem = SCIPbendersSubproblem(benders, probnumber);

         if( SCIPbendersSubproblemIsConvex(benders, probnumber) )
         {
            /* ending probing mode to reset the current node. The probing mode will be restarted at the next solve */
            if( SCIPinProbing(subproblem) )
            {
               SCIP_CALL( SCIPendProbing(subproblem) );
            }
         }
         else
         {
            /* if the subproblems were solved as part of an enforcement stage, then they will still be in probing mode. The
             * probing mode must first be finished and then the problem can be freed */
            if( SCIPgetStage(subproblem) >= SCIP_STAGE_TRANSFORMED && SCIPinProbing(subproblem) )
            {
               SCIP_CALL( SCIPendProbing(subproblem) );
            }

            SCIP_CALL( SCIPfreeTransform(subproblem) );
         }
      }
   }

   /* setting the setup flag for the subproblem to FALSE */
   SCIPbendersSetSubproblemIsSetup(benders, probnumber, FALSE);
   return SCIP_OKAY;
}

/** compares the subproblem objective value with the auxiliary variable value for optimality */
SCIP_RETCODE SCIPbendersCheckSubproblemOptimality(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            optimal             /**< flag to indicate whether the current subproblem is optimal for the master */
   )
{
   SCIP_Real auxiliaryvarval;
   SCIP_Real soltol;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   (*optimal) = FALSE;

   auxiliaryvarval = SCIPbendersGetAuxiliaryVarVal(benders, set, sol, probnumber);

   SCIP_CALL( SCIPsetGetRealParam(set, "benders/solutiontol", &soltol) );

   SCIPsetDebugMsg(set, "Subproblem %d - Auxiliary Variable: %g Subproblem Objective: %g Reldiff: %g Soltol: %g\n",
      probnumber, auxiliaryvarval, SCIPbendersGetSubproblemObjval(benders, probnumber),
      SCIPrelDiff(SCIPbendersGetSubproblemObjval(benders, probnumber), auxiliaryvarval), soltol);

   if( SCIPrelDiff(SCIPbendersGetSubproblemObjval(benders, probnumber), auxiliaryvarval) < soltol )
      (*optimal) = TRUE;

   return SCIP_OKAY;
}

/** returns the value of the auxiliary variable value in a master problem solution */
SCIP_Real SCIPbendersGetAuxiliaryVarVal(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP_VAR* auxiliaryvar;

   assert(benders != NULL);
   assert(set != NULL);

   auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, probnumber);
   assert(auxiliaryvar != NULL);

   return SCIPgetSolVal(set->scip, sol, auxiliaryvar);
}

/** Solves an independent subproblem to identify its lower bound. The lower bound is then used to update the bound on
 *  the auxiliary variable.
 */
SCIP_RETCODE SCIPbendersComputeSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber,         /**< the subproblem to be evaluated */
   SCIP_Real*            lowerbound,         /**< the lowerbound for the subproblem */
   SCIP_Bool*            infeasible          /**< was the subproblem found to be infeasible? */
   )
{
   SCIP* subproblem;
   SCIP_Real memorylimit;
   SCIP_Real timelimit;
   SCIP_Longint totalnodes;
   int disablecutoff;
   int verblevel;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;

   assert(benders != NULL);
   assert(set != NULL);

   /* getting the subproblem to evaluate */
   subproblem = SCIPbendersSubproblem(benders, probnumber);

   (*lowerbound) = -SCIPinfinity(subproblem);
   (*infeasible) = FALSE;

   SCIPverbMessage(set->scip, SCIP_VERBLEVEL_FULL, NULL, "Benders' decomposition: Computing a lower bound for"
      " subproblem %d\n", probnumber);

   SCIP_CALL( SCIPgetIntParam(subproblem, "display/verblevel", &verblevel) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
#ifdef SCIP_MOREDEBUG
   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_HIGH) );
#endif

   /* copying memory and time limits */
   SCIP_CALL( SCIPgetRealParam(subproblem, "limits/time", &timelimit) );
   SCIP_CALL( SCIPgetRealParam(subproblem, "limits/memory", &memorylimit) );
   SCIP_CALL( copyMemoryAndTimeLimits(set->scip, subproblem) );

   /* if the subproblem is independent, then the default SCIP settings are used. Otherwise, only the root node is solved
    * to compute a lower bound on the subproblem
    */
   SCIP_CALL( SCIPgetLongintParam(subproblem, "limits/totalnodes", &totalnodes) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/disablecutoff", &disablecutoff) );
   if( !SCIPbendersSubproblemIsIndependent(benders, probnumber) )
   {
      SCIP_CALL( SCIPsetLongintParam(subproblem, "limits/totalnodes", 1LL) );
      SCIP_CALL( SCIPsetIntParam(subproblem, "lp/disablecutoff", 1) );
   }

   /* if the subproblem not independent and is convex, then the probing LP is solved. Otherwise, the MIP is solved */
   if( SCIPbendersSubproblemIsConvex(benders, probnumber) )
   {
      assert(SCIPisLPConstructed(subproblem));

      SCIP_CALL( SCIPstartProbing(subproblem) );
      SCIP_CALL( SCIPsolveProbingLP(subproblem, -1, &lperror, &cutoff) );

      if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE )
         (*infeasible) = TRUE;
   }
   else
   {
      SCIP_EVENTHDLRDATA* eventhdlrdata;

      /* if the subproblem is not convex, then event handlers have been added to interrupt the solve. These must be
       * disabled
       */
      eventhdlrdata = SCIPeventhdlrGetData(SCIPfindEventhdlr(subproblem, MIPNODEFOCUS_EVENTHDLR_NAME));
      eventhdlrdata->solvecip = TRUE;

      SCIP_CALL( SCIPsolve(subproblem) );

      if( SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE )
         (*infeasible) = TRUE;
   }

   /* getting the lower bound value */
   if( !(*infeasible) )
      (*lowerbound) = SCIPgetDualbound(subproblem);
   else
      (*lowerbound) = -SCIPinfinity(subproblem);

   if( !SCIPbendersSubproblemIsIndependent(benders, probnumber) )
   {
      SCIP_CALL( SCIPsetLongintParam(subproblem, "limits/totalnodes", totalnodes) );
      SCIP_CALL( SCIPsetIntParam(subproblem, "lp/disablecutoff", disablecutoff) );
   }
   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", verblevel) );
   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/time", timelimit) );

   /* the subproblem must be freed so that it is reset for the subsequent Benders' decomposition solves. If the
    * subproblems are independent, they are not freed. SCIPfreeBendersSubproblem must still be called, but in this
    * function the independent subproblems are not freed. However, they will still be freed at the end of the
    * solving process for the master problem.
    */
   SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, probnumber) );

   return SCIP_OKAY;
}

/** Merges a subproblem into the master problem. This process just adds a copy of the subproblem variables and
 *  constraints to the master problem, but keeps the subproblem stored in the Benders' decomposition data structure. The reason for
 *  keeping the subproblem available is for when it is queried for solutions after the problem is solved.
 *
 *  Once the subproblem is merged into the master problem, then the subproblem is flagged as disabled. This means that
 *  it will not be solved in the subsequent subproblem solving loops.
 *
 *  The associated auxiliary variables are kept in the master problem. The objective function of the merged subproblem
 *  is added as an underestimator constraint.
 */
SCIP_RETCODE SCIPbendersMergeSubproblemIntoMaster(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of subproblem variables corresponding
                                              *   to the newly created master variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of subproblem constraints to the
                                                  corresponding newly created constraints, or NULL */
   int                   probnumber          /**< the number of the subproblem that will be merged into the master problem*/
   )
{
   SCIP* subproblem;
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_VAR** vars;
   SCIP_VAR* auxiliaryvar;
   SCIP_CONS** conss;
   SCIP_CONS* objcons;
   int nvars;
   int nconss;
   int i;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;
   char varname[SCIP_MAXSTRLEN];
   char consname[SCIP_MAXSTRLEN];
   const char* origvarname;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   SCIPverbMessage(set->scip, SCIP_VERBLEVEL_HIGH, NULL, "   Benders' decomposition: Infeasibility of subproblem %d can't "
      "be resolved. Subproblem %d is being merged into the master problem.\n", probnumber, probnumber);

   /* freeing the subproblem because it will be flagged as independent. Since the subproblem is flagged as independent,
    * it will no longer be solved or freed within the solving loop.
    */
   SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, probnumber) );

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(set->scip), SCIPgetNVars(subproblem)) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(set->scip), SCIPgetNConss(subproblem)) );
   }
   else
      localconsmap = consmap;

   /* retrieving the subproblem variable to build a subproblem mapping */
   vars = SCIPgetVars(subproblem);
   nvars = SCIPgetNVars(subproblem);

   /* creating the objective function constraint that will be added to the master problem */
   /* setting the name of the transferred cut */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "objectivecons_%d", probnumber );
   SCIP_CALL( SCIPcreateConsBasicLinear(set->scip, &objcons, consname, 0, NULL, NULL, -SCIPsetInfinity(set), 0.0) );
   SCIP_CALL( SCIPsetConsRemovable(set->scip, objcons, TRUE) );

   for( i = 0; i < nvars; i++ )
   {
      SCIP_VAR* mastervar = NULL;
      SCIP_Bool releasevar = FALSE;

      SCIP_CALL( SCIPgetBendersMasterVar(set->scip, benders, vars[i], &mastervar) );

      /* if the master problem variable is not NULL, then there is a corresponding variable in the master problem for
       * the given subproblem variable. In this case, the variable is added to the hashmap.
       */
      if( mastervar == NULL )
      {
         SCIP_VAR* origvar;
         SCIP_Real scalar;
         SCIP_Real constant;

         /* This is following the same process as in createVariableMappings. The original variable is used to map
          * between the subproblem and the master problem
          */
         origvar = vars[i];
         scalar = 1.0;
         constant = 0.0;
         SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

         /* retrieving the var name */
         origvarname = SCIPvarGetName(origvar);
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s", origvarname);

         /* creating and adding the variable to the Benders' decomposition master problem */
         SCIP_CALL( SCIPcreateVarBasic(set->scip, &mastervar, varname, SCIPvarGetLbOriginal(origvar),
            SCIPvarGetUbOriginal(origvar), 0.0, SCIPvarGetType(origvar)) );

         /* adding the variable to the master problem */
         SCIP_CALL( SCIPaddVar(set->scip, mastervar) );

         /* adds the variable to the objective function constraint */
         SCIP_CALL( SCIPaddCoefLinear(set->scip, objcons, mastervar, SCIPvarGetObj(origvar)) );

         /* the variable must be released */
         releasevar = TRUE;
      }

      /* creating the mapping betwen the subproblem var and the master var for the constraint copying */
      SCIP_CALL( SCIPhashmapInsert(localvarmap, vars[i], mastervar) );

      /* releasing the variable */
      if( releasevar )
      {
         SCIP_CALL( SCIPreleaseVar(set->scip, &mastervar) );
      }
   }

   /* getting the constraints from the subproblem that will be added to the master problem */
   conss = SCIPgetConss(subproblem);
   nconss = SCIPgetNConss(subproblem);

   /* getting a copy of all constraints and adding it to the master problem */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CONS* targetcons;
      SCIP_Bool initial;
      SCIP_Bool valid;

      /* NOTE: adding all subproblem constraints appears to cause an error when resolving the LP, which results in the
       * current incumbent being reported as optimal. To avoid this, only half of the subproblem constraints are added
       * the master problem. The remaining half are marked as lazy and are separated as required.
       */
      initial = (i < nconss/2);

      SCIP_CALL( SCIPgetConsCopy(subproblem, set->scip, conss[i], &targetcons, SCIPconsGetHdlr(conss[i]),
         localvarmap, localconsmap, NULL, initial, SCIPconsIsSeparated(conss[i]),
         SCIPconsIsEnforced(conss[i]), SCIPconsIsChecked(conss[i]), SCIPconsIsPropagated(conss[i]), FALSE,
         SCIPconsIsModifiable(conss[i]), SCIPconsIsDynamic(conss[i]), SCIPconsIsRemovable(conss[i]),
         FALSE, TRUE, &valid) );
      assert(SCIPconsIsInitial(conss[i]));
      assert(valid);

      SCIP_CALL( SCIPaddCons(set->scip, targetcons) );

      SCIP_CALL( SCIPreleaseCons(set->scip, &targetcons) );
   }

   /* freeing the hashmaps */
   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   /* adding the auxiliary variable to the objective constraint */
   auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, probnumber);
   SCIP_CALL( SCIPaddCoefLinear(set->scip, objcons, auxiliaryvar, -1.0) );

   /* adding the objective function constraint to the master problem */
   SCIP_CALL( SCIPaddCons(set->scip, objcons) );

   SCIP_CALL( SCIPreleaseCons(set->scip, &objcons) );

   /* the merged subproblem is no longer solved. This is indicated by setting the subproblem as disabled. The
    * subproblem still exists, but it is not solved in the solving loop.
    */
   SCIPbendersSetSubproblemEnabled(benders, probnumber, FALSE);

   return SCIP_OKAY;
}

/** Returns the corresponding master or subproblem variable for the given variable.
 *  This provides a call back for the variable mapping between the master and subproblems. */
SCIP_RETCODE SCIPbendersGetVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< the variable for which the corresponding variable is desired */
   SCIP_VAR**            mappedvar,          /**< the variable that is mapped to var */
   int                   probnumber          /**< the problem number for the desired variable, -1 for the master problem */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(mappedvar != NULL);
   assert(benders->bendersgetvar != NULL);

   (*mappedvar) = NULL;

   /* if the variable name matches the auxiliary variable, then the master variable is returned as NULL */
   if( strstr(SCIPvarGetName(var), AUXILIARYVAR_NAME) != NULL )
      return SCIP_OKAY;

   SCIP_CALL( benders->bendersgetvar(set->scip, benders, var, mappedvar, probnumber) );

   return SCIP_OKAY;
}

/** gets user data of Benders' decomposition */
SCIP_BENDERSDATA* SCIPbendersGetData(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->bendersdata;
}

/** sets user data of Benders' decomposition; user has to free old data in advance! */
void SCIPbendersSetData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSDATA*     bendersdata         /**< new Benders' decomposition user data */
   )
{
   assert(benders != NULL);

   benders->bendersdata = bendersdata;
}

/** sets copy callback of Benders' decomposition */
void SCIPbendersSetCopy(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSCOPY ((*benderscopy))    /**< copy callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->benderscopy = benderscopy;
}

/** sets destructor callback of Benders' decomposition */
void SCIPbendersSetFree(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSFREE ((*bendersfree))    /**< destructor of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersfree = bendersfree;
}

/** sets initialization callback of Benders' decomposition */
void SCIPbendersSetInit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINIT((*bendersinit))     /**< initialize the Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersinit = bendersinit;
}

/** sets deinitialization callback of Benders' decomposition */
void SCIPbendersSetExit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXIT((*bendersexit))     /**< deinitialize the Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersexit = bendersexit;
}

/** sets presolving initialization callback of Benders' decomposition */
void SCIPbendersSetInitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre))/**< initialize presolving for Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersinitpre = bendersinitpre;
}

/** sets presolving deinitialization callback of Benders' decomposition */
void SCIPbendersSetExitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre))/**< deinitialize presolving for Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersexitpre = bendersexitpre;
}

/** sets solving process initialization callback of Benders' decomposition */
void SCIPbendersSetInitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol))/**< solving process initialization callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersinitsol = bendersinitsol;
}

/** sets solving process deinitialization callback of Benders' decomposition */
void SCIPbendersSetExitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol))/**< solving process deinitialization callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersexitsol = bendersexitsol;
}

/** sets the pre subproblem solve callback of Benders' decomposition */
void SCIPbendersSetPresubsolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve))/**< called prior to the subproblem solving loop */
   )
{
   assert(benders != NULL);

   benders->benderspresubsolve = benderspresubsolve;
}

/** sets convex solve callback of Benders' decomposition */
void SCIPbendersSetSolvesubconvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex))/**< solving method for the convex Benders' decomposition subproblem */
   )
{
   assert(benders != NULL);

   benders->benderssolvesubconvex = benderssolvesubconvex;
}

/** sets solve callback of Benders' decomposition */
void SCIPbendersSetSolvesub(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub))/**< solving method for a Benders' decomposition subproblem */
   )
{
   assert(benders != NULL);

   benders->benderssolvesub = benderssolvesub;
}

/** sets post-solve callback of Benders' decomposition */
void SCIPbendersSetPostsolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve))/**< solving process deinitialization callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->benderspostsolve = benderspostsolve;
}

/** sets free subproblem callback of Benders' decomposition */
void SCIPbendersSetFreesub(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub))/**< the freeing callback for the subproblem */
   )
{
   assert(benders != NULL);

   benders->bendersfreesub = bendersfreesub;
}

/** gets name of Benders' decomposition */
const char* SCIPbendersGetName(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->name;
}

/** gets description of Benders' decomposition */
const char* SCIPbendersGetDesc(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->desc;
}

/** gets priority of Benders' decomposition */
int SCIPbendersGetPriority(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->priority;
}

/** sets priority of Benders' decomposition */
void SCIPbendersSetPriority(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the Benders' decomposition */
   )
{
   assert(benders != NULL);
   assert(set != NULL);

   benders->priority = priority;
   set->benderssorted = FALSE;
}

/** gets the number of subproblems for the Benders' decomposition */
int SCIPbendersGetNSubproblems(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition data structure */
   )
{
   assert(benders != NULL);

   return benders->nsubproblems;
}

/** returns the SCIP instance for a given subproblem */
SCIP* SCIPbendersSubproblem(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   return benders->subproblems[probnumber];
}

/** gets the number of times, the Benders' decomposition was called and tried to find a variable with negative reduced costs */
int SCIPbendersGetNCalls(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->ncalls;
}

/** gets the number of optimality cuts found by the collection of Benders' decomposition subproblems */
int SCIPbendersGetNCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->ncutsfound;
}

/** gets time in seconds used in this Benders' decomposition for setting up for next stages */
SCIP_Real SCIPbendersGetSetupTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return SCIPclockGetTime(benders->setuptime);
}

/** gets time in seconds used in this Benders' decomposition */
SCIP_Real SCIPbendersGetTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return SCIPclockGetTime(benders->bendersclock);
}

/** enables or disables all clocks of the Benders' decomposition, depending on the value of the flag */
void SCIPbendersEnableOrDisableClocks(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the Benders' decomposition be enabled? */
   )
{
   assert(benders != NULL);

   SCIPclockEnableOrDisable(benders->setuptime, enable);
   SCIPclockEnableOrDisable(benders->bendersclock, enable);
}

/** is Benders' decomposition initialized? */
SCIP_Bool SCIPbendersIsInitialized(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->initialized;
}

/** Are Benders' cuts generated from the LP solutions? */
SCIP_Bool SCIPbendersCutLP(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutlp;
}

/** Are Benders' cuts generated from the pseudo solutions? */
SCIP_Bool SCIPbendersCutPseudo(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutpseudo;
}

/** Are Benders' cuts generated from the relaxation solutions? */
SCIP_Bool SCIPbendersCutRelaxation(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutrelax;
}

/** should this Benders' use the auxiliary variables from the highest priority Benders' */
SCIP_Bool SCIPbendersShareAuxVars(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->shareauxvars;
}

/** adds a subproblem to the Benders' decomposition data */
SCIP_RETCODE SCIPbendersAddSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP*                 subproblem          /**< subproblem to be added to the data storage */
   )
{
   assert(benders != NULL);
   assert(subproblem != NULL);
   assert(benders->subproblems != NULL);
   assert(benders->naddedsubprobs + 1 <= benders->nsubproblems);

   benders->subproblems[benders->naddedsubprobs] = subproblem;

   benders->naddedsubprobs++;

   return SCIP_OKAY;
}

/** removes the subproblems from the Benders' decomposition data */
void SCIPbendersRemoveSubproblems(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);
   assert(benders->subproblems != NULL);

   BMSclearMemoryArray(&benders->subproblems, benders->naddedsubprobs);
   benders->naddedsubprobs = 0;
}

/** returns the auxiliary variable for the given subproblem */
SCIP_VAR* SCIPbendersGetAuxiliaryVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->auxiliaryvars[probnumber];
}

/** returns all auxiliary variables */
SCIP_VAR** SCIPbendersGetAuxiliaryVars(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->auxiliaryvars;
}

/** stores the objective function value of the subproblem for use in cut generation */
void SCIPbendersSetSubproblemObjval(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             objval              /**< the objective function value for the subproblem */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* updating the best objval */
   if( objval < benders->bestsubprobobjval[probnumber] )
      benders->bestsubprobobjval[probnumber] = objval;

   benders->subprobobjval[probnumber] = objval;
}

/** returns the objective function value of the subproblem for use in cut generation */
SCIP_Real SCIPbendersGetSubproblemObjval(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobobjval[probnumber];
}

/** sets the flag indicating whether a subproblem is convex
 *
 *  It is possible that this can change during the solving process. One example is when the three-phase method is
 *  employed, where the first phase solves the convex relaxation of both the master and subproblems, the second phase
 *  reintroduces the integrality constraints to the master problem and the third phase then reintroduces integrality
 *  constraints to the subproblems.
 */
void SCIPbendersSetSubproblemIsConvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isconvex            /**< flag to indicate whether the subproblem is convex */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   if( isconvex && !benders->subprobisconvex[probnumber] )
      benders->nconvexsubprobs++;
   else if( !isconvex && benders->subprobisconvex[probnumber] )
      benders->nconvexsubprobs--;

   benders->subprobisconvex[probnumber] = isconvex;

   assert(benders->nconvexsubprobs >= 0 && benders->nconvexsubprobs <= benders->nsubproblems);
}

/** returns whether the subproblem is convex
 *
 *  This means that the dual solution can be used to generate cuts.
 */
SCIP_Bool SCIPbendersSubproblemIsConvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobisconvex[probnumber];
}

/** returns the number of subproblems that are convex */
int SCIPbendersGetNConvexSubproblems(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nconvexsubprobs;
}

/** changes all of the master problem variables in the given subproblem to continuous. */
SCIP_RETCODE SCIPbendersChgMastervarsToCont(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_VAR** vars;
   int nbinvars;
   int nintvars;
   int nimplvars;
   int chgvarscount;
   int origintvars;
   int i;
   SCIP_Bool infeasible;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   subproblem = SCIPbendersSubproblem(benders, probnumber);
   assert(subproblem != NULL);

   /* only set the master problem variable to continuous if they have not already been changed. */
   if( !SCIPbendersGetMastervarsCont(benders, probnumber) )
   {
      SCIP_VAR* mastervar;

      /* retrieving the variable data */
      SCIP_CALL( SCIPgetVarsData(subproblem, &vars, NULL, &nbinvars, &nintvars, &nimplvars, NULL) );

      origintvars = nbinvars + nintvars + nimplvars;

      chgvarscount = 0;

      /* looping over all integer variables to change the master variables to continuous */
      i = 0;
      while( i < nbinvars + nintvars + nimplvars )
      {
         SCIP_CALL( SCIPbendersGetVar(benders, set, vars[i], &mastervar, -1) );

         if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS && mastervar != NULL )
         {
            /* changing the type of the subproblem variable corresponding to mastervar to CONTINUOUS */
            SCIP_CALL( SCIPchgVarType(subproblem, vars[i], SCIP_VARTYPE_CONTINUOUS, &infeasible) );

            assert(!infeasible);

            chgvarscount++;
            SCIP_CALL( SCIPgetVarsData(subproblem, NULL, NULL, &nbinvars, &nintvars, &nimplvars, NULL) );
         }
         else
            i++;
      }

      /* if all of the integer variables have been changed to continuous, then the subproblem must now be an LP. In this
       * case, the subproblem is initialised and then put into probing mode */
      if( chgvarscount > 0 && chgvarscount == origintvars )
      {
         SCIP_CALL( initialiseLPSubproblem(benders, set, probnumber) );
         SCIPbendersSetSubproblemIsConvex(benders, probnumber, TRUE);
      }

      SCIP_CALL( SCIPbendersSetMastervarsCont(benders, probnumber, TRUE) );
   }

   return SCIP_OKAY;
}

/** sets the subproblem setup flag */
void SCIPbendersSetSubproblemIsSetup(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             issetup             /**< flag to indicate whether the subproblem has been setup */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   benders->subprobsetup[probnumber] = issetup;
}

/** returns the subproblem setup flag */
SCIP_Bool SCIPbendersSubproblemIsSetup(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobsetup[probnumber];
}

/** sets the independent subproblem flag */
void SCIPbendersSetSubproblemIsIndependent(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isindep             /**< flag to indicate whether the subproblem is independent */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* if the user has defined solving or freeing functions, then it is not possible to declare a subproblem as
    * independent. This is because declaring a subproblem as independent changes the solving loop, so it would change
    * the expected behaviour of the user defined plugin. If a user calls this function, then an error will be returned.
    */
   if( benders->benderssolvesubconvex != NULL || benders->benderssolvesub != NULL || benders->bendersfreesub != NULL )
   {
      SCIPerrorMessage("The user has defined either bendersSolvesubconvex%d, bendersSolvesub%d or bendersFreesub%s. "
         "Thus, it is not possible to declare the independence of a subproblem.\n", benders->name, benders->name,
         benders->name);
      SCIPABORT();
   }
   else
   {
      SCIP_Bool activesubprob;

      /* if the active status of the subproblem changes, then we must update the activesubprobs counter */
      activesubprob = subproblemIsActive(benders, probnumber);

      benders->indepsubprob[probnumber] = isindep;

      /* updating the activesubprobs counter */
      if( activesubprob && !subproblemIsActive(benders, probnumber) )
         benders->nactivesubprobs--;
      else if( !activesubprob && subproblemIsActive(benders, probnumber) )
         benders->nactivesubprobs++;

      assert(benders->nactivesubprobs >= 0 && benders->nactivesubprobs <= SCIPbendersGetNSubproblems(benders));
   }
}

/** returns whether the subproblem is independent */
SCIP_Bool SCIPbendersSubproblemIsIndependent(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->indepsubprob[probnumber];
}

/** Sets whether the subproblem is enabled or disabled. A subproblem is disabled if it has been merged into the master
 *  problem.
 */
void SCIPbendersSetSubproblemEnabled(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             enabled             /**< flag to indicate whether the subproblem is enabled */
   )
{
   SCIP_Bool activesubprob;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* if the active status of the subproblem changes, then we must update the activesubprobs counter */
   activesubprob = subproblemIsActive(benders, probnumber);

   benders->subprobenabled[probnumber] = enabled;

   /* updating the activesubprobs counter */
   if( activesubprob && !subproblemIsActive(benders, probnumber) )
      benders->nactivesubprobs--;
   else if( !activesubprob && subproblemIsActive(benders, probnumber) )
      benders->nactivesubprobs++;

   assert(benders->nactivesubprobs >= 0 && benders->nactivesubprobs <= SCIPbendersGetNSubproblems(benders));
}

/** returns whether the subproblem is enabled, i.e. the subproblem is still solved in the solving loop. */
SCIP_Bool SCIPbendersSubproblemIsEnabled(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobenabled[probnumber];
}

/** sets a flag to indicate whether the master variables are all set to continuous */
SCIP_RETCODE SCIPbendersSetMastervarsCont(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             arecont             /**< flag to indicate whether the master problem variables are continuous */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* if the master variables were all continuous and now are not, then the subproblem must exit probing mode and be
    * changed to non-LP subproblem */
   if( benders->mastervarscont[probnumber] && !arecont )
   {
      if( SCIPinProbing(SCIPbendersSubproblem(benders, probnumber)) )
      {
         SCIP_CALL( SCIPendProbing(SCIPbendersSubproblem(benders, probnumber)) );
      }

      SCIPbendersSetSubproblemIsConvex(benders, probnumber, FALSE);
   }

   benders->mastervarscont[probnumber] = arecont;

   return SCIP_OKAY;
}

/** returns whether the master variables are all set to continuous */
SCIP_Bool SCIPbendersGetMastervarsCont(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->mastervarscont[probnumber];
}

/** returns the number of cuts that have been transferred from sub SCIPs to the master SCIP */
int SCIPbendersGetNTransferredCuts(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition data structure */
   )
{
   assert(benders != NULL);

   return benders->ntransferred;
}

/** updates the lower bound for the subproblem. If the lower bound is not greater than the previously stored lowerbound,
 *  then no update occurs.
 */
void SCIPbendersUpdateSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             lowerbound          /**< the lower bound */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   if( lowerbound > benders->subproblowerbound[probnumber] )
      benders->subproblowerbound[probnumber] = lowerbound;
   else
   {
      SCIPdebugMessage("The lowerbound %g for subproblem %d is less than the currently stored lower bound %g\n",
         lowerbound, probnumber, benders->subproblowerbound[probnumber]);
   }
}

/** returns the stored lower bound for the given subproblem */
SCIP_Real SCIPbendersGetSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subproblowerbound[probnumber];
}

/** sets the sorted flags in the Benders' decomposition */
void SCIPbendersSetBenderscutsSorted(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition structure */
   SCIP_Bool             sorted              /**< the value to set the sorted flag to */
   )
{
   assert(benders != NULL);

   benders->benderscutssorted = sorted;
   benders->benderscutsnamessorted = sorted;
}

/** inserts a Benders' cut into the Benders' cuts list */
SCIP_RETCODE SCIPbendersIncludeBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERSCUT*      benderscut          /**< Benders' cut */
   )
{
   assert(benders != NULL);
   assert(benderscut != NULL);

   if( benders->nbenderscuts >= benders->benderscutssize )
   {
      benders->benderscutssize = SCIPsetCalcMemGrowSize(set, benders->nbenderscuts+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&benders->benderscuts, benders->benderscutssize) );
   }
   assert(benders->nbenderscuts < benders->benderscutssize);

   benders->benderscuts[benders->nbenderscuts] = benderscut;
   benders->nbenderscuts++;
   benders->benderscutssorted = FALSE;

   return SCIP_OKAY;
}

/** returns the Benders' cut of the given name, or NULL if not existing */
SCIP_BENDERSCUT* SCIPfindBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   const char*           name                /**< name of Benderscut' decomposition */
   )
{
   int i;

   assert(benders != NULL);
   assert(name != NULL);

   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      if( strcmp(SCIPbenderscutGetName(benders->benderscuts[i]), name) == 0 )
         return benders->benderscuts[i];
   }

   return NULL;
}

/** returns the array of currently available Benders' cuts; active Benders' decomposition are in the first slots of
 * the array
 */
SCIP_BENDERSCUT** SCIPbendersGetBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   if( !benders->benderscutssorted )
   {
      SCIPsortPtr((void**)benders->benderscuts, SCIPbenderscutComp, benders->nbenderscuts);
      benders->benderscutssorted = TRUE;
      benders->benderscutsnamessorted = FALSE;
   }

   return benders->benderscuts;
}

/** returns the number of currently available Benders' cuts */
int SCIPbendersGetNBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nbenderscuts;
}

/** sets the priority of a Benders' decomposition */
SCIP_RETCODE SCIPbendersSetBenderscutPriority(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' cut */
   int                   priority            /**< new priority of the Benders' decomposition */
   )
{
   assert(benders != NULL);
   assert(benderscut != NULL);

   benderscut->priority = priority;
   benders->benderscutssorted = FALSE;

   return SCIP_OKAY;
}

/** sorts Benders' decomposition cuts by priorities */
void SCIPbendersSortBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   if( !benders->benderscutssorted )
   {
      SCIPsortPtr((void**)benders->benderscuts, SCIPbenderscutComp, benders->nbenderscuts);
      benders->benderscutssorted = TRUE;
      benders->benderscutsnamessorted = FALSE;
   }
}

/** sorts Benders' decomposition cuts by name */
void SCIPbendersSortBenderscutsName(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   if( !benders->benderscutsnamessorted )
   {
      SCIPsortPtr((void**)benders->benderscuts, SCIPbenderscutCompName, benders->nbenderscuts);
      benders->benderscutssorted = FALSE;
      benders->benderscutsnamessorted = TRUE;
   }
}
