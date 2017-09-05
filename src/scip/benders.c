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

/**@file   benders.c
 * @brief  methods for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
//#define SCIP_MOREDEBUG
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
#include "scip/misc_benders.h"

/* Defaults for parameters */
#define SCIP_DEFAULT_USEMW               FALSE  /** Should Magnanti-Wong cut strengthening be used? */
#define SCIP_DEFAULT_COMPUTERELINT       FALSE  /** Should the relative interior point be computed? */
#define SCIP_DEFAULT_UPDATEFACTOR          0.5  /** The factor to update the relative interior point? */
#define SCIP_DEFAULT_TRANSFERCUTS         TRUE  /** Should Benders' cuts generated in LNS heuristics be transferred to the main SCIP instance? */
#define SCIP_DEFAULT_CUTSASCONS           TRUE  /** Should the transferred cuts be added as constraints? */
//#define SCIP_DEFAULT_CUTSASCONS          FALSE  /** Should the transferred cuts be added as constraints? */
#define SCIP_DEFAULT_MIPCHECKFREQ            5  /** the number of iterations that the MIP is checked, -1 for always. */
#define SCIP_DEFAULT_LNSCHECK              TRUE /** should the Benders' decomposition be used in LNS heuristics */

#define AUXILIARYVAR_NAME     "##bendersauxiliaryvar"
#define MW_AUXILIARYVAR_NAME  "##MWauxiliaryvar##"

/* event handler properties */
#define NODEFOCUS_EVENTHDLR_NAME         "bendersnodefocus"
#define NODEFOCUS_EVENTHDLR_DESC         "node focus event handler for Benders' decomposition"
#define MIPNODEFOCUS_EVENTHDLR_NAME      "bendersmipsolvenodefocus"
#define MIPNODEFOCUS_EVENTHDLR_DESC      "node focus event handler the MIP solve method for Benders' decomposition"
#define UPPERBOUND_EVENTHDLR_NAME        "bendersupperbound"
#define UPPERBOUND_EVENTHDLR_DESC        "found solution event handler to terminate subproblem solve for a given upper bound"

#define HEUR_RENS "rens"
#define HEUR_RINS "rins"


struct SCIP_EventhdlrData
{
   int                   numruns;            /**< the number of times that the problem has been solved */
   SCIP_Real             upperbound;         /**< an upper bound for the problem */
};

/* ---------------- Callback methods of node focus event handler ---------------- */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBendersNodefocus)
{  /*lint --e{715}*/

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   /* sending an interrupt solve signal to return the control back to the Benders' decomposition plugin.
    * This will ensure the SCIP stage is SCIP_STAGE_SOLVING, allowing the use of probing mode. */
   SCIP_CALL( SCIPinterruptSolve(scip) );

   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1));

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersNodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL));

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolBendersNodefocus)
{
   assert(scip != NULL);

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

#if 0 /* not sure whether this is actually needed */
   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1));
#endif

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
   if( eventhdlrdata->numruns == 0 )
      SCIP_CALL( SCIPinterruptSolve(scip) );

   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1));

   eventhdlrdata->numruns++;

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitBendersMipnodefocus)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->numruns = 0;

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersMipnodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL));

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolBendersMipnodefocus)
{
   assert(scip != NULL);

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   //SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1));

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTFREE(eventFreeBendersMipnodefocus)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &eventhdlrdata);

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
      SCIP_CALL( SCIPinterruptSolve(scip) );

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitBendersUpperbound)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->upperbound = -SCIPinfinity(scip);

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersUpperbound)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL));

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolBendersUpperbound)
{
   assert(scip != NULL);

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, -1));

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTFREE(eventFreeBendersUpperbound)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &eventhdlrdata);

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
   SCIPeventhdlrSetData(eventhdlr, eventhdlrdata);

   return SCIP_OKAY;
}

/* Local methods */

/* A workaround for GCG. This is a temp vardata that is set for the auxiliary variables */
struct SCIP_VarData
{
   int vartype;
};


#if 0
/** fix master problem variables to core point */
static
SCIP_RETCODE fixVarsToCorePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< the event handler data */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* mastervar;
   SCIP_Real solval;
   int nvars;
   int i;
   SCIP_Bool fixed;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   for( i = 0; i < nvars; i++ )
   {
      mastervar = SCIPgetBendersMasterVar(eventhdlrdata->masterprob, eventhdlrdata->benders, vars[i]);

      if( mastervar != NULL )
      {
         solval = SCIPgetSolVal(eventhdlrdata->masterprob, eventhdlrdata->relintsol, mastervar);

         /* fixing the variable in the subproblem */
         SCIP_CALL( SCIPfixVar(scip, vars[i], solval, &infeasible, &fixed) );
      }
   }

   return SCIP_OKAY;
}
#endif

/** perform the Magnanti-Wong technique */
/*  TODO: If the problem contains SetPPC or logic or constraints, then it is not possible to apply M-W.  */
static
SCIP_RETCODE performMagnantiWongTechnique(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnum             /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_CONS** conss;
   SCIP_CONS** origconss;
   SCIP_Real* vals;
   SCIP_Real objval;
   int nconss;
   int norigconss;
   int i;
   SCIP_Bool infeasible;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnum >= 0 && probnum < benders->nsubproblems);

   /* a relative interior point is required for the Magnanti-Wong technique */
   if( benders->relintsol == NULL )
      return SCIP_OKAY;

   subproblem = SCIPbendersSubproblem(benders, probnum);

   if( SCIPgetStatus(subproblem) != SCIP_STATUS_OPTIMAL )
      return SCIP_OKAY;

   assert(subproblem != NULL);

   /* retreiving the coefficients of the auxiliary variable */

   /* getting all constraints from the subproblem */
   conss = SCIPgetConss(subproblem);
   nconss = SCIPgetNConss(subproblem);

   SCIP_CALL( SCIPallocBlockMemoryArray(subproblem, &vals, nconss) );

   /* looping over all constraints to find the currently active rows. This will give the LHS/RHS that is necessary for
    * the Magnanti-Wong auxiliary variable. */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_Real dualsol;
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      int nconsvars;
      int j;

      vals[i] = 0;

      dualsol = BDconsGetDualsol(subproblem, conss[i]);

      if( SCIPisPositive(subproblem, dualsol) )
         vals[i] = BDconsGetLhs(subproblem, conss[i]);
      else
         vals[i] = BDconsGetRhs(subproblem, conss[i]);

      /* retreiving the bound contribution of the fixed variables */
      nconsvars = BDconsGetNVars(subproblem, conss[i]);
      SCIP_CALL( SCIPallocBufferArray(subproblem, &consvars, nconsvars) );
      SCIP_CALL( SCIPallocBufferArray(subproblem, &consvals, nconsvars) );
      SCIP_CALL( BDconsGetVars(subproblem, conss[i], consvars, nconsvars) );
      SCIP_CALL( BDconsGetVals(subproblem, conss[i], consvals, nconsvars) );

      /* loop over all variables with non-zero coefficient */
      for( j = 0; j < nconsvars; j++ )
      {
         SCIP_VAR* consvar;
         SCIP_Real consval;

         consvar = consvars[j];
         consval = consvals[j];

         /* TODO: Do we need the problem variable? */
         consvar = SCIPvarGetProbvar(consvars[j]);

         if( SCIPisEQ(subproblem, SCIPvarGetLbLocal(consvar), SCIPvarGetUbLocal(consvar)) )
            vals[i] -= consval * SCIPvarGetLbLocal(consvar);
      }

      SCIPfreeBufferArray(subproblem, &consvars);
      SCIPfreeBufferArray(subproblem, &consvals);
   }

   /* Getting the objective function value for the current solution */
   objval = SCIPgetSolTransObj(subproblem, SCIPgetBestSol(subproblem));

   /* freeing the solution data from the previous run. */
   SCIP_CALL( SCIPendProbing(subproblem) );
   SCIP_CALL( SCIPstartProbing(subproblem) );

   /* updating the coefficients of the auxiliary variable */

   /* getting all constraints from the subproblem */
   origconss = SCIPgetOrigConss(subproblem);
   norigconss = SCIPgetNOrigConss(subproblem);
   assert(norigconss == nconss);

   /* looping over all constraints to add the auxiliary variable. */
   for( i = 0; i < norigconss; i++ )
      SCIP_CALL( BDconsChgCoef(subproblem, origconss[i], benders->mwauxiliaryvars[probnum], vals[i]) );

   /* changing the objective coefficient of the auxiliary variable */
   SCIP_CALL( SCIPchgVarObj(subproblem, benders->mwauxiliaryvars[probnum], objval) );

   /* changing the upper and lower bound of the auxiliary variable */
   SCIP_CALL( SCIPchgVarLb(subproblem, benders->mwauxiliaryvars[probnum], -SCIPinfinity(subproblem)) );
   SCIP_CALL( SCIPchgVarUb(subproblem, benders->mwauxiliaryvars[probnum], SCIPinfinity(subproblem)) );

   /* solving the subproblem with the updated solution and additional variable */
   SCIP_CALL( SCIPbendersExecSubproblemSolve(benders, set, benders->relintsol, probnum, 0, TRUE, &infeasible, LP) );
   assert(!infeasible);

   /* freeing the values array memory */
   SCIPfreeBlockMemoryArray(subproblem, &vals, nconss);

   return SCIP_OKAY;
}



/** clean up after the Magnanti-Wong technique
 *  TODO: Need to look at if keeping the auxiliary variable is affecting the cuts generated. */
static
SCIP_RETCODE cleanupMagnantiWong(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnum             /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnum >= 0 && probnum < benders->nsubproblems);

   subproblem = SCIPbendersSubproblem(benders, probnum);

   assert(subproblem != NULL);

   /* a relative interior point is required for the Magnanti-Wong technique */
   if( !benders->usemagnantiwong || SCIPgetStatus(subproblem) != SCIP_STATUS_OPTIMAL )
      return SCIP_OKAY;

   /* freeing the solution data from the previous run. */
   SCIP_CALL( SCIPendProbing(subproblem) );
   SCIP_CALL( SCIPstartProbing(subproblem) );

   /* creating the auxiliary variable */
   SCIP_CALL( SCIPfixVar(subproblem, benders->mwauxiliaryvars[probnum], 0.0, &infeasible, &fixed) );
   assert(!infeasible);
   assert(fixed);

   return SCIP_OKAY;
}


/** Adds the Magnanti-Wong auxiliary variables to the subproblems */
static
SCIP_RETCODE addMagnantiWongAuxiliaryVars(
   SCIP_BENDERS*         benders             /**< benders */
   )
{
   SCIP* subproblem;
   SCIP_VAR* auxiliaryvar;
   int nsubproblems;
   int i;
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */

   assert(benders != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* creating the auxiliary variable */
   for( i = 0; i < nsubproblems; i++ )
   {
      subproblem = SCIPbendersSubproblem(benders, i);
      assert(subproblem != NULL);

      /* if no optimality cuts have been added for this subproblem, then the auxiliary variable will be created and
       * added */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s", MW_AUXILIARYVAR_NAME );
      SCIP_CALL( SCIPcreateVarBasic(subproblem, &auxiliaryvar, varname, 0.0, 0.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(subproblem, auxiliaryvar) );

      benders->mwauxiliaryvars[i] = auxiliaryvar;

      SCIP_CALL( SCIPreleaseVar(subproblem, &auxiliaryvar) );
   }

   return SCIP_OKAY;
}



/** updates the core point for the Magnanti-Wong method */
static
SCIP_RETCODE updateRelativeInteriorPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_VAR** vars;
   SCIP_Real updateval;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   for( i = 0; i < nvars; i++ )
   {
      updateval = benders->updatefactor * SCIPgetSolVal(scip, benders->relintsol, vars[i])
         + (1 - benders->updatefactor) * SCIPgetSolVal(scip, sol, vars[i]);
      SCIP_CALL( SCIPsetSolVal(scip, benders->relintsol, vars[i], updateval) );
   }

   return SCIP_OKAY;
}




/** constructs a core point for the Magnanti-Wong method */
static
SCIP_RETCODE constructRelativeInteriorPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol                 /**< the solution to be constructed */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(*sol == NULL);

   SCIP_CALL( SCIPcreateSol(scip, sol, NULL) );

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   for( i = 0; i < nvars; i++ )
      SCIP_CALL( SCIPsetSolVal(scip, (*sol), vars[i], SCIPgetSolVal(scip, NULL, vars[i])) );
      //SCIP_CALL( SCIPsetSolVal(scip, (*sol), vars[i], 1.0) );

   return SCIP_OKAY;
}




/** compute the relative interior point for the Magnanti-Wong method */
static
SCIP_RETCODE computeRelativeInteriorPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   /* check whether we have to compute a relative interior point */
   if( benders->computerelint )
   {
      /* if relative interior point is not available ... */
      if ( benders->relintsol == NULL )
      {
         //SCIP_Longint nlpiters;
         SCIP_Real timelimit;
         int iterlimit;

         /* prepare time limit */
         SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
         if ( ! SCIPisInfinity(scip, timelimit) )
            timelimit -= SCIPgetSolvingTime(scip);
         /* exit if no time left */
         if ( timelimit <= 0.0 )
            return SCIP_OKAY;

         iterlimit = INT_MAX;
#if 0
         /* determine iteration limit */
         if ( benders->maxlpiterfactor < 0.0 || SCIPisInfinity(scip, benders->maxlpiterfactor) )
            iterlimit = INT_MAX;
         else
         {
            /* determine iteration limit; the number of iterations in the root is only set after its solution, but the
             * total number of LP iterations is always updated. */
            if ( SCIPgetDepth(scip) == 0 )
               nlpiters = SCIPgetNLPIterations(scip);
            else
               nlpiters = SCIPgetNRootLPIterations(scip);
            iterlimit = (int)(sepadata->maxlpiterfactor * nlpiters);
            iterlimit = MAX(iterlimit, SCIP_MIN_LPITERS);
            if ( iterlimit <= 0 )
               return SCIP_OKAY;
         }
#endif

         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Computing relative interior point (time limit: %g, iter limit: %d) ...\n", timelimit, iterlimit);
         SCIP_CALL( SCIPcomputeLPRelIntPoint(scip, TRUE, TRUE, timelimit, iterlimit,
               &benders->relintsol) );
      }
   }
   else
   {
      /* get best solution (NULL if not present) */
      if( benders->relintsol == NULL )
      {
         if( SCIPgetBestSol(scip) == NULL )
         {
            if( sol == NULL )
               SCIP_CALL( constructRelativeInteriorPoint(scip, &benders->relintsol) );
            else
               benders->relintsol = sol;
         }
         else
            benders->relintsol = SCIPgetBestSol(scip);
      }
      else
      {
         SCIP_CALL( updateRelativeInteriorPoint(scip, benders, sol) );
      }
   }

   return SCIP_OKAY;
}

/* prepares the data for Magnanti-Wong cut strengthening technique */
static
SCIP_RETCODE prepareMagnantiWongTechnique(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);

   if( benders->usemagnantiwong )
   {
      /* compute the relative interior point. */
      SCIP_CALL( computeRelativeInteriorPoint(scip, benders, sol) );

      /* setting current solution in the event handler data */
      benders->currentsol = sol;
   }

   return SCIP_OKAY;
}


/** adds the auxiliary variables to the Benders' decomposition master problem */
static
SCIP_RETCODE addAuxiliaryVariablesToMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< benders */
   )
{
   SCIP_VAR* auxiliaryvar;
   SCIP_VARDATA* vardata;
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */
   int i;

   /* this is a workaround for GCG. GCG expects that the variable has vardata when added. So a dummy vardata is created */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = -1;

   for( i = 0; i < SCIPbendersGetNSubproblems(benders); i++ )
   {
      /* if no optimality cuts have been added for this subproblem, then the auxiliary variable will be created and
       * added */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s_%d", AUXILIARYVAR_NAME, i );
      SCIP_CALL( SCIPcreateVarBasic(scip, &auxiliaryvar, varname, -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
            SCIP_VARTYPE_CONTINUOUS) );

      SCIPvarSetData(auxiliaryvar, vardata);

      SCIP_CALL( SCIPaddVar(scip, auxiliaryvar) );

      benders->auxiliaryvars[i] = auxiliaryvar;

      SCIP_CALL( SCIPreleaseVar(scip, &auxiliaryvar) );
   }

   SCIPfreeBlockMemory(scip, &vardata);

   return SCIP_OKAY;
}

/** assigns the copied auxiliary variables in the target scip to the target benders data */
static
SCIP_RETCODE assignAuxiliaryVariables(
   SCIP*                 scip,               /**< SCIP data structure, the target scip */
   SCIP_BENDERS*         benders             /**< benders */
   )
{
   SCIP_VAR* targetvar;
   SCIP_VARDATA* vardata;
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   /* this is a workaround for GCG. GCG expects that the variable has vardata when added. So a dummy vardata is created */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = -1;

   for( i = 0; i < SCIPbendersGetNSubproblems(benders); i++ )
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s_%d", AUXILIARYVAR_NAME, i );

      /* finding the variable in the copied problem that has the same name as the auxiliary variable */
      targetvar = SCIPfindVar(scip, varname);

      SCIPvarSetData(targetvar, vardata);

      benders->auxiliaryvars[i] = targetvar;
   }

   SCIPfreeBlockMemory(scip, &vardata);

   return SCIP_OKAY;
}

/* sets the subproblem objective value array to -infinity */
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
      SCIPbendersSetSubprobObjval(benders, SCIPinfinity(subproblem), i);
   }
}


/** compares two benders w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPbendersComp)
{  /*lint --e{715}*/
   return ((SCIP_BENDERS*)elem2)->priority - ((SCIP_BENDERS*)elem1)->priority;
}

/** comparison method for sorting benders w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPbendersCompName)
{
   return strcmp(SCIPbendersGetName((SCIP_BENDERS*)elem1), SCIPbendersGetName((SCIP_BENDERS*)elem2));
}

/** method to call, when the priority of a benders was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdBendersPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetBendersPriority() to mark the benderss unsorted */
   SCIP_CALL( SCIPsetBendersPriority(scip, (SCIP_BENDERS*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given benders to a new scip */
SCIP_RETCODE SCIPbendersCopyInclude(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_SET*             sourceset,          /**< SCIP_SET of SCIP to copy from */
   SCIP_SET*             set,                /**< SCIP_SET of SCIP to copy to */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   )
{
   SCIP_BENDERS* targetbenders;  /* the copy of the Benders' decomposition struct in the target set */
   int i;

   assert(benders != NULL);
   assert(set != NULL);
   assert(valid != NULL);
   assert(set->scip != NULL);

   if( benders->benderscopy != NULL )
   {
      SCIPsetDebugMsg(set, "including benders %s in subscip %p\n", SCIPbendersGetName(benders), (void*)set->scip);
      SCIP_CALL( benders->benderscopy(set->scip, benders, valid) );

      /* if the copy was valid, then the Benders cuts are copied. */
      if( valid )
      {
         targetbenders = SCIPsetFindBenders(set, SCIPbendersGetName(benders));

         /* storing the pointer to the source scip instance */
         targetbenders->sourcescip = sourceset->scip;

         /* the flag is set to indicate that the Benders' decomposition is a copy */
         targetbenders->iscopy = TRUE;

         /* calling the copy method for the Benders' cuts */
         SCIPbendersSortBenderscuts(benders);
         for( i = 0; i < benders->nbenderscuts; i++ )
         {
            SCIP_CALL( SCIPbenderscutCopyInclude(benders->benderscuts[i], set) );
         }
      }
   }

   return SCIP_OKAY;
}

/** creates a Benders' decomposition structure
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
   SCIP_DECL_BENDERSCOPY ((*benderscopy)),   /**< copy method of benders or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree)),   /**< destructor of Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit)),   /**< initialize Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit)),   /**< deinitialize Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre)),/**< presolving initialization method for Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre)),/**< presolving deinitialization method for Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol)),/**< solving process initialization method of Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol)),/**< solving process deinitialization method of Benders' decomposition */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSEXEC ((*bendersexec)),   /**< the execution method of the Benders' decomposition algorithm */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub)),/**< creates a Benders' decomposition subproblem */
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

   SCIP_ALLOC( BMSallocMemory(benders) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->desc, desc, strlen(desc)+1) );
   (*benders)->priority = priority;
   (*benders)->cutlp = cutlp;
   (*benders)->cutpseudo = cutpseudo;
   (*benders)->cutrelax = cutrelax;
   (*benders)->benderscopy = benderscopy;
   (*benders)->bendersfree = bendersfree;
   (*benders)->bendersinit = bendersinit;
   (*benders)->bendersexit = bendersexit;
   (*benders)->bendersinitpre = bendersinitpre;
   (*benders)->bendersexitpre = bendersexitpre;
   (*benders)->bendersinitsol = bendersinitsol;
   (*benders)->bendersexitsol = bendersexitsol;
   (*benders)->bendersgetvar = bendersgetvar;
   (*benders)->bendersexec = bendersexec;
   (*benders)->benderscreatesub = benderscreatesub;
   (*benders)->benderssolvesub = benderssolvesub;
   (*benders)->benderspostsolve = benderspostsolve;
   (*benders)->bendersfreesub = bendersfreesub;
   (*benders)->bendersdata = bendersdata;
   SCIP_CALL( SCIPclockCreate(&(*benders)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*benders)->bendersclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*benders)->ncalls = 0;
   (*benders)->noptcutsfound = 0;
   (*benders)->nfeascutsfound = 0;
   (*benders)->ntransferred = 0;
   (*benders)->active = FALSE;
   (*benders)->initialized = FALSE;
   (*benders)->cutsascons = SCIP_DEFAULT_CUTSASCONS;
   (*benders)->sourcescip = NULL;
   (*benders)->iscopy = FALSE;
   (*benders)->mastervarsmap = NULL;
   (*benders)->maxlpiterfactor = 1.0;
   (*benders)->updatefactor = SCIP_DEFAULT_UPDATEFACTOR;
   (*benders)->coreptupdated = FALSE;
   (*benders)->relintsol = NULL;
   (*benders)->currentsol = NULL;
   (*benders)->addedsubprobs = 0;
   (*benders)->nlpsubprobs = 0;
   (*benders)->subprobscreated = FALSE;

   (*benders)->benderscuts = NULL;
   (*benders)->nbenderscuts = 0;
   (*benders)->benderscutssize = 0;
   (*benders)->benderscutssorted = FALSE;
   (*benders)->benderscutsnamessorted = FALSE;

   /* setting the subproblem arrays to NULL */
   (*benders)->subproblems = NULL;
   (*benders)->mwauxiliaryvars = NULL;
   (*benders)->auxiliaryvars = NULL;
   (*benders)->subprobobjval = NULL;
   (*benders)->bestsubprobobjval = NULL;
   (*benders)->subprobislp = NULL;
   (*benders)->mastervarscont = NULL;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of benders <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->priority, FALSE, priority, INT_MIN/4, INT_MAX/4,
                  paramChgdBendersPriority, (SCIP_PARAMDATA*)(*benders)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutlp", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "should Benders' cuts be generated for LP solutions?");
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->cutlp, FALSE, cutlp, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutpseudo", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "should Benders' cuts be generated for pseudo solutions?");
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->cutpseudo, FALSE, cutpseudo, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutrelax", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "should Benders' cuts be generated for relaxation solutions?");
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->cutrelax, FALSE, cutrelax, NULL, NULL) ); /*lint !e740*/

   /* These parameters are left for the user to decide in a settings file. This departs from the usual SCIP convention
    * where the settings available at the creation of the plugin can be set in the function call. */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/usemagnantiwong", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "Should the Magnanti-Wong cut strengthening technique be used?");
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->usemagnantiwong, FALSE, SCIP_DEFAULT_USEMW, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/computerelint", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "Should the relative interior point be computed?");
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->computerelint, FALSE, SCIP_DEFAULT_COMPUTERELINT, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/transfercuts", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "Should Benders' cuts from LNS heuristics be transferred to the main SCIP instance?");
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->transfercuts, FALSE, SCIP_DEFAULT_TRANSFERCUTS, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/mipcheckfreq", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "The frequency at which the MIP subproblems are checked, -1 for always");
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->mipcheckfreq, FALSE, SCIP_DEFAULT_MIPCHECKFREQ, -1, INT_MAX, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/lnscheck", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "Should Benders' decomposition be used in LNS heurisics?");
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benders)->lnscheck, FALSE, SCIP_DEFAULT_LNSCHECK, NULL, NULL) ); /*lint !e740*/

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

   /* freeing the Benders' cuts */
   for( i = 0; i < (*benders)->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutFree(&((*benders)->benderscuts[i]), set) );
   }
   BMSfreeMemoryArrayNull(&(*benders)->benderscuts);

   if( (*benders)->relintsol != NULL )
      SCIP_CALL( SCIPfreeSol(set->scip, &(*benders)->relintsol) );

   SCIPclockFree(&(*benders)->bendersclock);
   SCIPclockFree(&(*benders)->setuptime);
   BMSfreeMemoryArray(&(*benders)->name);
   BMSfreeMemoryArray(&(*benders)->desc);
   BMSfreeMemory(benders);

   return SCIP_OKAY;
}


/** initialises an LP subproblem by putting the problem into probing mode */
static
SCIP_RETCODE initialiseLPSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_Bool infeasible;
   SCIP_Bool cutoff;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   subproblem = SCIPbendersSubproblem(benders, probnumber);
   assert(subproblem != NULL);

   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(subproblem, &eventhdlr, NODEFOCUS_EVENTHDLR_NAME, NODEFOCUS_EVENTHDLR_DESC,
         eventExecBendersNodefocus, NULL) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(subproblem, eventhdlr, eventInitsolBendersNodefocus) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(subproblem, eventhdlr, eventExitsolBendersNodefocus) );
   assert(eventhdlr != NULL);

   /* Getting the problem into the right SCIP stage for solving */
   SCIP_CALL( SCIPbendersSolveSubproblemMIP(benders, probnumber, &infeasible, LP, TRUE, FALSE) );

   /* Constructing the LP that can be solved in later iterations */
   SCIP_CALL( SCIPconstructLP(subproblem, &cutoff) );

   assert(SCIPgetStage(subproblem) == SCIP_STAGE_SOLVING);

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
   int nbinvars;
   int nintvars;
   int nimplintvars;
   int nsubproblems;
   int i;

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

      /* getting the number of integer and binary variables to determine the problem type */
      SCIP_CALL( SCIPgetVarsData(subproblem, NULL, NULL, &nbinvars, &nintvars, &nimplintvars, NULL) );

      /* if there are no binary and integer variables, then the subproblem is an LP.
       * In this case, the SCIP instance is put into probing mode via the use of an event handler. */
      if( nbinvars == 0 && nintvars == 0 && nimplintvars == 0 )
      {
         SCIPbendersSetSubprobIsLP(benders, i, TRUE);

         /* if the user has not implemented a solve subproblem callback, then the subproblem solves are performed
          * internally. To be more efficient the subproblem is put into probing mode. */
         if( benders->benderssolvesub == NULL && SCIPgetStage(subproblem) <= SCIP_STAGE_PROBLEM )
            SCIP_CALL( initialiseLPSubproblem(benders, i) );
      }
      else
      {
         SCIP_EVENTHDLRDATA* eventhdlrdata_mipnodefocus;
         SCIP_EVENTHDLRDATA* eventhdlrdata_upperbound;

         /* because the subproblems could be reused in the copy, the event handler is not created again.
          * NOTE: This currently works with the benders_default implementation. It may not be very general. */
         if( benders->benderssolvesub == NULL && !benders->iscopy )
         {
            SCIP_CALL( SCIPallocBlockMemory(subproblem, &eventhdlrdata_mipnodefocus) );
            SCIP_CALL( SCIPallocBlockMemory(subproblem, &eventhdlrdata_upperbound) );
            eventhdlrdata_mipnodefocus->numruns = 0;
            eventhdlrdata_upperbound->upperbound = -SCIPinfinity(subproblem);

            /* include the first LP solved event handler into the subproblem */
            SCIP_CALL( SCIPincludeEventhdlrBasic(subproblem, &eventhdlr, MIPNODEFOCUS_EVENTHDLR_NAME,
                  MIPNODEFOCUS_EVENTHDLR_DESC, eventExecBendersMipnodefocus, eventhdlrdata_mipnodefocus) );
            SCIP_CALL( SCIPsetEventhdlrInit(subproblem, eventhdlr, eventInitBendersMipnodefocus) );
            SCIP_CALL( SCIPsetEventhdlrInitsol(subproblem, eventhdlr, eventInitsolBendersMipnodefocus) );
            SCIP_CALL( SCIPsetEventhdlrExitsol(subproblem, eventhdlr, eventExitsolBendersMipnodefocus) );
            SCIP_CALL( SCIPsetEventhdlrFree(subproblem, eventhdlr, eventFreeBendersMipnodefocus) );
            assert(eventhdlr != NULL);


            /* include the upper bound interrupt event handler into the subproblem */
            SCIP_CALL( SCIPincludeEventhdlrBasic(subproblem, &eventhdlr, UPPERBOUND_EVENTHDLR_NAME,
                  UPPERBOUND_EVENTHDLR_DESC, eventExecBendersUpperbound, eventhdlrdata_upperbound) );
            SCIP_CALL( SCIPsetEventhdlrInit(subproblem, eventhdlr, eventInitBendersUpperbound) );
            SCIP_CALL( SCIPsetEventhdlrInitsol(subproblem, eventhdlr, eventInitsolBendersUpperbound) );
            SCIP_CALL( SCIPsetEventhdlrExitsol(subproblem, eventhdlr, eventExitsolBendersUpperbound) );
            SCIP_CALL( SCIPsetEventhdlrFree(subproblem, eventhdlr, eventFreeBendersUpperbound) );
            assert(eventhdlr != NULL);
         }
      }
   }

   benders->subprobscreated = TRUE;

   return SCIP_OKAY;

}


/** creates a variable mapping between the master problem variables of the source scip and the sub scip */
static
SCIP_RETCODE createMasterVarMapping(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_VAR** vars;  /* the variables from the copied sub SCIP */
   SCIP_VAR* sourcevar;
   SCIP_VAR* origvar;
   SCIP_Real scalar;
   SCIP_Real constant;
   int nvars;
   int i;
   SCIP_RETCODE retcode;

   assert(benders != NULL);
   assert(set != NULL);
   assert(benders->iscopy);
   assert(benders->mastervarsmap == NULL);

   /* getting the master problem variable data */
   vars = SCIPgetVars(set->scip);
   nvars = SCIPgetNVars(set->scip);

   /* creating the hashmap for the mapping between the master variable of the target and source scip */
   SCIP_CALL( SCIPhashmapCreate(&benders->mastervarsmap, SCIPblkmem(set->scip), nvars) );

   for( i = 0; i < nvars; i++ )
   {
      origvar = vars[i];

      /* The variable needs to be transformed back into an original variable. If the variable is already original, then
       * this function just returns the same variable */
      retcode = SCIPvarGetOrigvarSum(&origvar, &scalar, &constant);
      if( retcode != SCIP_OKAY )
         assert(FALSE);

      sourcevar = SCIPfindVar(benders->sourcescip, SCIPvarGetName(origvar));
      if( sourcevar != NULL )
         SCIP_CALL( SCIPhashmapInsert(benders->mastervarsmap, vars[i], sourcevar) );
   }

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

   /* Checking whether the benderssolvesub and the bendersfreesub are both implemented or both are not implemented */
   if( (benders->benderssolvesub == NULL && benders->bendersfreesub != NULL)
      || (benders->benderssolvesub != NULL && benders->bendersfreesub == NULL) )
   {
      SCIPerrorMessage("Benders' decomposition <%s> requires that both bendersSolvesub%s and bendersFreesub%s are \
         implemented or neither\n", benders->name, benders->name, benders->name);
      return SCIP_INVALIDCALL;
   }

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
      benders->noptcutsfound = 0;
      benders->nfeascutsfound = 0;
   }

   /* start timing */
   SCIPclockStart(benders->setuptime, set);

   /* creates the subproblems and sets up the probing mode for LP subproblems. This function calls the benderscreatesub
    * callback. */
   SCIP_CALL( createSubproblems(benders, set) );

   if( benders->bendersinit != NULL )
      SCIP_CALL( benders->bendersinit(set->scip, benders) );

   /* if the Magnanti-Wong technique is used, then auxiliary variable must be added to the subproblem */
   if( benders->usemagnantiwong )
      SCIP_CALL( addMagnantiWongAuxiliaryVars(benders) );

   /* if the Benders' decomposition is a copy, then a variable mapping between the master problem variables is required */
   if( benders->iscopy )
      SCIP_CALL( createMasterVarMapping(benders, set) );


   /* initialising the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutInit(benders->benderscuts[i], set) );
   }

   benders->initialized = TRUE;

   /* stop timing */
   SCIPclockStop(benders->setuptime, set);

   return SCIP_OKAY;
}


/** create and add transferred cut */
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

   /* creating an empty row/constraint for the transferred cut */
   if( sourcebenders->cutsascons )
   {
      SCIP_CALL( SCIPcreateConsBasicLinear(sourcescip, &transfercons, cutname, 0, NULL, NULL, lhs, rhs) );
      SCIPconsSetRemovable(transfercons, TRUE);
   }
   else
      SCIP_CALL( SCIPcreateEmptyRowCons(sourcescip, &transfercut, consbenders, cutname, lhs, rhs, FALSE,
            FALSE, TRUE) );

   fail = FALSE;
   for( i = 0; i < nvars; i++ )
   {
      /* getting the source var from the hash map */
      sourcevar = (SCIP_VAR*) SCIPhashmapGetImage(benders->mastervarsmap, vars[i]);

      /* if the source variable is not found, then the mapping in incomplete. So the constraint can not be
       * transferred. */
      if( sourcevar == NULL )
      {
         fail = TRUE;
         break;
      }

      if( sourcebenders->cutsascons )
         SCIP_CALL( SCIPaddCoefLinear(sourcescip, transfercons, sourcevar, vals[i]) );
      else
         SCIP_CALL( SCIPaddVarToRow(sourcescip, transfercut, sourcevar, vals[i]) );

      /* NOTE: There could be a problem with the auxiliary variables. They may not be copied. */
   }

   /* if all of the source variables were found to generate the cut */
   if( !fail )
   {
      if( sourcebenders->cutsascons )
         SCIP_CALL( SCIPaddCons(sourcescip, transfercons) );
      else
         SCIP_CALL( SCIPaddPoolCut(sourcescip, transfercut) );

      sourcebenders->ntransferred++;
   }

   /* release the row/constraint */
   if( sourcebenders->cutsascons )
      SCIP_CALL( SCIPreleaseCons(sourcescip, &transfercons) );
   else
      SCIP_CALL( SCIPreleaseRow(sourcescip, &transfercut) );

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
      SCIP_CALL( SCIPbenderscutGetCons(benderscut, &addedcons, &naddedcons) );

      /* looping over all added constraints to costruct the cut for the source scip */
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

            /* create and add the cut to be transferred from the sub SCIP to the source SCIP */
            SCIP_CALL( createAndAddTransferredCut(sourcescip, benders, vars, vals, lhs, rhs, nvars) );
         }
      }

      /* retreiving the Benders' cuts added cuts */
      SCIP_CALL( SCIPbenderscutGetCuts(benderscut, &addedcuts, &naddedcuts) );

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
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   if( !benders->initialized )
   {
      SCIPerrorMessage("Benders' decomposition <%s> not initialized\n", benders->name);
      return SCIP_INVALIDCALL;
   }

   if( benders->bendersexit != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersexit(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   /* if the Benders' decomposition is a copy, then
    * - the generated cuts will be transferred to the source scip, and
    * - the hash map must be freed */
   if( benders->iscopy )
   {
      SCIP_CALL( transferBendersCuts(benders->sourcescip, set->scip, benders) );
      SCIPhashmapFree(&benders->mastervarsmap);
   }

   /* calling the exit method for the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutExit(benders->benderscuts[i], set) );
   }

   benders->initialized = FALSE;

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
      /* adding the auxiliary variables to the master problem */
      SCIP_CALL( addAuxiliaryVariablesToMaster(set->scip, benders) );
   }
   else
   {
      /* the copied auxiliary variables must be assigned to the target benders */
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
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of Benders' decomposition */
   if( benders->bendersexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersexitsol(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   /* calling the exitsol method for the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutExitsol(benders->benderscuts[i], set) );
   }

   return SCIP_OKAY;
}

/** activates benders such that it is called in LP solving loop */
SCIP_RETCODE SCIPbendersActivate(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nsubproblems        /**< the number subproblems used in this decomposition */
   )
{
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

      /* allocating memory for the subproblems arrays */
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subproblems, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->mwauxiliaryvars, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->auxiliaryvars, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobobjval, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->bestsubprobobjval, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobislp, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->mastervarscont, benders->nsubproblems) );

      for( i = 0; i < benders->nsubproblems; i++ )
      {
         benders->subproblems[i] = NULL;
         benders->mwauxiliaryvars[i] = NULL;
         benders->auxiliaryvars[i] = NULL;
         benders->subprobobjval[i] = SCIPsetInfinity(set);
         benders->bestsubprobobjval[i] = SCIPsetInfinity(set);
         benders->subprobislp[i] = FALSE;
         benders->mastervarscont[i] = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** deactivates benders such that it is no longer called in LP solving loop */
SCIP_RETCODE SCIPbendersDeactivate(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_INIT || set->stage == SCIP_STAGE_PROBLEM);

   if( benders->active )
   {
      benders->active = FALSE;
      set->nactivebenders--;
      set->benderssorted = FALSE;

      /* freeing the memory allocated during the activation of the Benders' decomposition */
      BMSfreeMemoryArray(&benders->mastervarscont);
      BMSfreeMemoryArray(&benders->subprobislp);
      BMSfreeMemoryArray(&benders->bestsubprobobjval);
      BMSfreeMemoryArray(&benders->subprobobjval);
      BMSfreeMemoryArray(&benders->auxiliaryvars);
      BMSfreeMemoryArray(&benders->mwauxiliaryvars);
      BMSfreeMemoryArray(&benders->subproblems);
   }

   return SCIP_OKAY;
}

/** returns whether the given Benders decomposition is in use in the current problem */
SCIP_Bool SCIPbendersIsActive(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   assert(benders != NULL);

   return benders->active;
}

/** solves the subproblem using the current master problem solution. */
/*  TODO: consider allowing the possibility to pass solution information back from the subproblems instead of the scip
 *  instance. This would allow the use of different solvers for the subproblems, more importantly allowing the use of an
 *  LP solver for LP subproblems. */
SCIP_RETCODE SCIPbendersExec(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RESULT*          result,             /**< result of the pricing process */
   SCIP_Bool*            infeasible,         /**< is the master problem infeasible with respect to the Benders' cuts? */
   SCIP_BENDERSENFOTYPE  type                /**< the enforcement type calling this function */
   )
{
   SCIP_BENDERSCUT** benderscuts;
   int nsubproblems;
   int nbenderscuts;
   int nsolveloops;     /* the number of times the subproblems are solved. An additional loop is required when the subproblem in integer */
   int ncutloops;
   int i;
   int j;
   int k;
   int l;
   SCIP_Bool optimal;
   SCIP_Bool allchecked;      /* flag to indicate whether all subproblems have been checked */
   int nchecked;              /* the number of subproblems that have been checked */
   SCIP_Bool* subisinfeas;
   SCIP_Bool onlylpcheck;     /* should only the LP be checked in the presence of integer subproblems */

   (*infeasible) = FALSE;

   /* It is assumed that the problem is optimal, until a subproblem is found not to be optimal. However, not all
    * subproblems could be checked in each iteration. As such, it is not possible to state that the problem is optimal
    * if not all subproblems are checked. Situations where this may occur is when a subproblem is a MIP and only the LP
    * is solved. Also, in a distributed computation, then it may be adventageous to only solve some subproblems before
    * resolving the master problem. As such, for a problem to be optimal, then (optimal && allchecked) == TRUE */
   optimal = TRUE;
   allchecked = FALSE;
   nchecked = 0;

   assert(benders != NULL);
   assert(result != NULL);
   assert(benders->bendersexec != NULL);

   /* if the Benders' decomposition is called from a sub-scip, it is assumed that this is an LNS heuristic. As such, the
    * check is not performed and the solution is assumed to be feasible */
   if( benders->iscopy && !benders->lnscheck )
   {
      (*result) = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   //if( sol != NULL && SCIPsolGetHeur(sol) != NULL && strcmp(SCIPheurGetName(SCIPsolGetHeur(sol)), HEUR_RINS) == 0 && benders->lnscheck )
   //{
      //(*result) = SCIP_FEASIBLE;
      //return SCIP_OKAY;
   //}

   onlylpcheck = benders->iscopy && benders->lnscheck;

   /* stating that the core point has not been updated for this iteration */
   benders->coreptupdated = FALSE;

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* allocating memory for the infeasible subproblem array */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &subisinfeas, nsubproblems) );

   for( i = 0; i < nsubproblems; i++ )
      subisinfeas[i] = FALSE;

   /* sets the stored objective function values of the subproblems to infinity */
   resetSubproblemObjectiveValue(benders);

   SCIPdebugMessage("Performing the subproblem solving process.\n");

   SCIP_CALL( benders->bendersexec(set->scip, benders) );

   *result = SCIP_DIDNOTRUN;

   /* by default the number of solve loops is 1. This is the case is all subproblems are LP or the user has defined a
    * benderssolvesub callback. If there is a subproblem that is not an LP, then 2 solve loops are performed. The first
    * loop is the LP solving loop, the second solves the subproblem to integer optimality. */
   nsolveloops = 1;
   //if( SCIPbendersGetNLPSubprobs(benders) < SCIPbendersGetNSubproblems(benders) && benders->benderssolvesub == NULL && type == CHECK && benders->ncalls % benders->mipcheckfreq == 0 )
      //nsolveloops = 2;

   for( l = 0; l < nsolveloops; l++ )
   {
      if( type == CHECK && sol == NULL )
      {
         (*infeasible) = TRUE;
      }
      else
      {
         /* solving each of the subproblems for Benders decomposition */
         /* TODO: ensure that the each of the subproblems solve and update the parameters with the correct return values */
         for( i = 0; i < nsubproblems; i++ )
         {
            SCIP_Bool subinfeas = FALSE;
            SCIP_Bool lpsub = SCIPbendersSubprobIsLP(benders, i);

            /* for the second solving loop, if the problem is an LP, it is not solved again. If the problem is a MIP,
             * then the subproblem objective function value is set to infinity. However, if the subproblem is proven
             * infeasible from the LP, then the IP loop is not performed. */
            if( l > 0 )
            {
               if( lpsub || subisinfeas[i] )
                  continue;
               else
                  SCIPbendersSetSubprobObjval(benders, SCIPinfinity(SCIPbendersSubproblem(benders, i)), i);
            }

            SCIP_CALL( SCIPbendersExecSubproblemSolve(benders, set, sol, i, l, FALSE, &subinfeas, type) );

#ifdef SCIP_DEBUG
            if( type == LP )
               SCIPdebugMessage("LP: Subproblem %d (%f < %f)\n", i, SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i),
                  SCIPbendersGetSubprobObjval(benders, i));
#endif
            (*infeasible) = (*infeasible) || subinfeas;
            subisinfeas[i] = subinfeas;

            /* if the subproblems are being solved as part of the conscheck, then we break once an infeasibility is found.
             * The result pointer is set to (*infeasible) and the execution is halted. */
            if( type == CHECK )
            {
               /* if the subproblem is feasible, then it is necessary to update the value of the auxiliary variable to the
                * objective function value of the subproblem. */
               if( !subinfeas )
               {
                  SCIP_Bool suboptimal;

                  SCIP_CALL( SCIPbendersCheckAuxiliaryVar(benders, set, sol, i, &suboptimal) );

                  if( lpsub || benders->benderssolvesub != NULL || l > 0 || onlylpcheck )
                     optimal = optimal && suboptimal;

#ifdef SCIP_DEBUG
                  if( lpsub || l > 0 )
                  {
                     if( suboptimal )
                        SCIPdebugMessage("Subproblem %d is Optimal (%f >= %f)\n", i,
                           SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i), SCIPbendersGetSubprobObjval(benders, i));
                     else
                        SCIPdebugMessage("Subproblem %d is NOT Optimal (%f < %f)\n", i,
                           SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i), SCIPbendersGetSubprobObjval(benders, i));
                  }
#endif

                  /* only increment the checked count if the subproblem is not an LP, or the solve loop is the MIP
                   * solving loop. Hence, the LP are solved once and the MIPs are solved twice */
                  if( lpsub || (l > 0 && !lpsub) || onlylpcheck )
                     nchecked++;
               }
            }
         }
      }


      /* Preparing the data for the Magnanti-Wong technique */
      if( !(*infeasible) )
         prepareMagnantiWongTechnique(set->scip, benders, sol);


      /* the result flag is set to FEASIBLE as default. If a cut is added, then this is changed to CONS_ADDED */
      (*result) = SCIP_FEASIBLE;

      /* Generating cuts for the subproblems. */
      benderscuts = SCIPbendersGetBenderscuts(benders);
      nbenderscuts = SCIPbendersGetNBenderscuts(benders);

      /* if the Magnanti-Wong enhancement technique is used, then an additional cutloop is performed. */
      /* TODO: create a parameter to modify the behaviour of the cutloops. */
      ncutloops = 1;
      if( benders->usemagnantiwong && !(*infeasible) )
         ncutloops = 2;

      /* It is only possible to add cuts to the problem if it has not already been solved */
      if( SCIPsetGetStage(set) < SCIP_STAGE_SOLVED )
      {
         int addedcuts = 0;

         /* This is done in two loops. The first is by subproblem and the second is by cut type. */
         for( i = 0; i < nsubproblems; i++ )
         {
            SCIP_Bool lpsub = SCIPbendersSubprobIsLP(benders, i);

            for( j = 0; j < ncutloops; j++ )
            {
               /* Performing cut enhancement techniques */
               if( j == 1 )
                  SCIP_CALL( performMagnantiWongTechnique(benders, set, i) );

               for( k = 0; k < nbenderscuts; k++ )
               {
                  int prevaddedcuts;

                  assert(benderscuts[k] != NULL);

                  prevaddedcuts = SCIPbenderscutGetNAddedCons(benderscuts[k]) + SCIPbenderscutGetNAddedCuts(benderscuts[k]);

                  /* if the subproblem is an LP, then only LP based cuts are generated. This is also only performed in
                   * the first iteration of the solve loop. */
                  if( (l == 0 && SCIPbenderscutIsLPCut(benderscuts[k]))
                     || (l > 0 && !lpsub && !SCIPbenderscutIsLPCut(benderscuts[k])) )
                     SCIP_CALL( SCIPbenderscutExec(benderscuts[k], set, benders, sol, i, type, result) );

                  addedcuts += (SCIPbenderscutGetNAddedCons(benderscuts[k]) + SCIPbenderscutGetNAddedCuts(benderscuts[k]) - prevaddedcuts);
               }
            }
         }

         /* if no cuts were added, then the number of solve loops is increased */
         if( addedcuts == 0 && SCIPbendersGetNLPSubprobs(benders) < SCIPbendersGetNSubproblems(benders)
            && benders->benderssolvesub == NULL && type == CHECK && !onlylpcheck )
            nsolveloops = 2;
      }
   }

   allchecked = (nchecked == nsubproblems);

   if( type == CHECK )
   {
      /* if the subproblems are being solved as part of conscheck, then the results flag must be returned after the solving
       * has completed. No cut is generated during conscheck. */
      if( (*infeasible) || !allchecked )
         (*result) = SCIP_INFEASIBLE;
      else
      {
         (*result) = SCIP_FEASIBLE;

         /* if the subproblems are not infeasible, but they are also not optimal, then the check return the result of
          * feasible and the flag of infeasible. */
         (*infeasible) = !optimal;
      }
   }
   else if( type == PSEUDO )
   {
      if( (*infeasible) || !(optimal && allchecked) )
         (*result) = SCIP_INFEASIBLE;
      else
         (*result) = SCIP_FEASIBLE;
   }

   /* calling the post-solve call back for the Benders' decomposition algorithm. This allows the user to work directly
    * with the solved subproblems and the master problem */
   if( benders->benderspostsolve != NULL )
      SCIP_CALL( benders->benderspostsolve(set->scip, benders, (*infeasible)) );

   /* freeing the subproblems after the cuts are generated */
   for( i = 0; i < benders->nsubproblems; i++ )
   {
      /* cleaning up after performing Magnanti-Wong */
      SCIP_CALL( cleanupMagnantiWong(benders, set, i) );

      SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, i) );
   }

   /* increment the number of calls to the Benders' decomposition subproblem solve */
   benders->ncalls++;

   /* freeing memory */
   SCIPsetFreeBufferArray(set, &subisinfeas);

   return SCIP_OKAY;
}

/** solves the subproblems. */
SCIP_RETCODE SCIPbendersExecSubproblemSolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnum,            /**< the subproblem number */
   int                   solveloop,          /**< the solve loop iteration. The first iter is for LP, the second for IP */
   SCIP_Bool             enhancement,        /**< is the solve performed as part of and enhancement? */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type                /**< the enforcement type calling this function */
   )
{
   SCIP* subproblem;
   SCIP_SOL* bestsol;

   assert(benders != NULL);
   assert(probnum >= 0 && probnum < benders->nsubproblems);

   SCIPdebugMessage("Benders decomposition: solving subproblem %d\n", probnum);

   /* if the subproblem solve callback is implemented, then that is used instead of the default setup */
   if( benders->benderssolvesub != NULL)
      SCIP_CALL( benders->benderssolvesub(set->scip, benders, sol, probnum, infeasible) );
   else
   {
      /* setting up the subproblem */
      if( solveloop == 0 )
         SCIP_CALL( SCIPbendersSetupSubproblem(benders, set, sol, probnum) );
      else
         SCIP_CALL( updateEventhdlrUpperbound(benders, probnum, SCIPbendersGetAuxiliaryVarVal(benders, set, sol, probnum)) );


      /* solving the subproblem */
      if( SCIPbendersSubprobIsLP(benders, probnum) )
         SCIP_CALL( SCIPbendersSolveSubproblemLP(benders, probnum, infeasible) );
      else
         SCIP_CALL( SCIPbendersSolveSubproblemMIP(benders, probnum, infeasible, type, FALSE, FALSE) );
   }

   subproblem = SCIPbendersSubproblem(benders, probnum);

   bestsol = SCIPgetBestSol(subproblem);

   if( !enhancement )
   {
      if( solveloop == 0 )
      {
         if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_OPTIMAL )
            SCIPbendersSetSubprobObjval(benders, SCIPgetSolOrigObj(subproblem, NULL), probnum);
         else if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE )
            SCIPbendersSetSubprobObjval(benders, SCIPinfinity(set->scip), probnum);
         else
            assert(FALSE);
      }
      else
      {
         assert(solveloop == 1);

         if( SCIPgetStatus(subproblem) == SCIP_STATUS_OPTIMAL )
            SCIPbendersSetSubprobObjval(benders, SCIPgetSolOrigObj(subproblem, bestsol), probnum);
         else if( SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE )
            SCIPbendersSetSubprobObjval(benders, SCIPinfinity(set->scip), probnum);
         else if( SCIPgetStatus(subproblem) == SCIP_STATUS_USERINTERRUPT || SCIPgetStatus(subproblem) == SCIP_STATUS_BESTSOLLIMIT )
            SCIPbendersSetSubprobObjval(benders, SCIPgetSolOrigObj(subproblem, bestsol), probnum);
         else
            assert(FALSE);
      }
   }

   return SCIP_OKAY;
}

/** sets up the subproblem using the solution to the master problem  */
SCIP_RETCODE SCIPbendersSetupSubproblem(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnum             /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_VAR** vars;
   SCIP_VAR* mastervar;
   SCIP_Real solval;
   int nvars;
   int i;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnum >= 0 && probnum < SCIPbendersGetNSubproblems(benders));

   /* changing all of the master problem variable to continuous. */
   SCIP_CALL( SCIPbendersChgMastervarsToCont(benders, set, probnum) );

   subproblem = SCIPbendersSubproblem(benders, probnum);

   /* if the Benders subproblem is an LP, then probing mode must be started */
   if( SCIPbendersSubprobIsLP(benders, probnum) )
      SCIP_CALL( SCIPstartProbing(subproblem) );

   vars = SCIPgetVars(subproblem);
   nvars = SCIPgetNVars(subproblem);

   /* looping over all variables in the subproblem to find those corresponding to the master problem variables. */
   /* TODO: It should be possible to store the pointers to the master variables to speed up the subproblem setup */
   for( i = 0; i < nvars; i++ )
   {
      mastervar = SCIPbendersGetVar(benders, set, vars[i], -1);

      if( mastervar != NULL )
      {
         solval = SCIPgetSolVal(set->scip, sol, mastervar);

         /* fixing the variable in the subproblem */
         SCIP_CALL( SCIPfixVar(subproblem, vars[i], solval, &infeasible, &fixed) );

         assert(SCIPisFeasEQ(subproblem, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])));
         assert(fixed);
         assert(!infeasible);
      }
   }

   return SCIP_OKAY;
}

/** Solve a Benders' decomposition subproblems. This will either call the user defined method or the generic solving
 * methods. If the generic method is called, then the subproblem must be set up before calling this method. */
extern
SCIP_RETCODE SCIPbendersSolveSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_Bool             solvemip            /**< directly solve the MIP subproblem */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* if the subproblem solve callback is implemented, then that is used instead of the default setup */
   if( benders->benderssolvesub != NULL)
      SCIP_CALL( benders->benderssolvesub(set->scip, benders, sol, probnumber, infeasible) );
   else
   {
      /* solving the subproblem */
      if( SCIPbendersSubprobIsLP(benders, probnumber) )
      {
         /* if the subproblem is not in probing mode, then it must be put into that mode for the LP solve. */
         if( !SCIPinProbing(SCIPbendersSubproblem(benders, probnumber)) )
            SCIP_CALL( SCIPstartProbing(SCIPbendersSubproblem(benders, probnumber)) );

         SCIP_CALL( SCIPbendersSolveSubproblemLP(benders, probnumber, infeasible) );
      }
      else
         SCIP_CALL( SCIPbendersSolveSubproblemMIP(benders, probnumber, infeasible, type, FALSE, solvemip) );
   }

   return SCIP_OKAY;
}

/** solves the LP of the Benders' decomposition subproblem. This requires that the subproblem is in probing mode */
SCIP_RETCODE SCIPbendersSolveSubproblemLP(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible          /**< a flag to indicate whether all subproblems are feasible */
   )
{
   SCIP* subproblem;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;

   /* previous parameter settings */
   int prevCutoffParam;
   char prevInitAlgParam;
   char prevResolveAlgParam;
   SCIP_Bool prevDualParam;

   assert(benders != NULL);
   assert(infeasible != NULL);

   (*infeasible) = FALSE;


   /* TODO: This should be solved just as an LP, so as a MIP. There is too much overhead with the MIP.
    * Need to change status check for checking the LP. */
   subproblem = SCIPbendersSubproblem(benders, probnumber);

   assert(SCIPisLPConstructed(subproblem));
   assert(SCIPinProbing(subproblem));

   /* modifying all of the parameters */
   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/disablecutoff", &prevCutoffParam) );
   SCIPsetIntParam(subproblem, "lp/disablecutoff", 1);

   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/initalgorithm", &prevInitAlgParam) );
   SCIPsetCharParam(subproblem, "lp/initalgorithm", 'd');
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/resolvealgorithm", &prevResolveAlgParam) );
   SCIPsetCharParam(subproblem, "lp/resolvealgorithm", 'd');

   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/alwaysgetduals", &prevDualParam) );
   SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", TRUE);

   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
   //SCIP_CALL( SCIPsetBoolParam(subproblem, "display/lpinfo", TRUE) );

   SCIP_CALL( SCIPsolveProbingLP(subproblem, -1, &lperror, &cutoff) );

   assert(!lperror);

   if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE )
      (*infeasible) = TRUE;
   else if( SCIPgetLPSolstat(subproblem) != SCIP_LPSOLSTAT_OPTIMAL )
      assert(FALSE);

   //SCIP_CALL( SCIPprintStatistics(subprob, NULL) );

   SCIPsetIntParam(subproblem, "lp/disablecutoff", prevCutoffParam);
   SCIPsetCharParam(subproblem, "lp/initalgorithm", prevInitAlgParam);
   SCIPsetCharParam(subproblem, "lp/resolvealgorithm", prevResolveAlgParam);
   SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", prevDualParam);

   return SCIP_OKAY;
}

/** solves the Benders' decomposition subproblem. */
SCIP_RETCODE SCIPbendersSolveSubproblemMIP(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_Bool             initialisation,   /**< indicates whether the MIP is solved as part of an initalisation */
   SCIP_Bool             solvemip            /**< directly solve the MIP subproblem */
   )
{
   SCIP* subproblem;

   /* previous parameter settings */
   int prevCutoffParam;
   int prevLPScalingParam;
   int prevPropMaxroundsParam;
   int prevPropMaxroundsRootParam;
   int prevLimitsBestsolParam;
   int prevLimitsMaxOrigSolParam;
   char prevInitAlgParam;
   char prevResolveAlgParam;
   SCIP_Bool prevConfParam;
   SCIP_Bool prevDualParam;
   SCIP_Bool prevScaleObjParam;
   SCIP_Bool mipchecksolve;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;

   assert(benders != NULL);
   assert(infeasible != NULL);

   (*infeasible) = FALSE;
   mipchecksolve = FALSE;


   /* TODO: This should be solved just as an LP, so as a MIP. There is too much overhead with the MIP.
    * Need to change status check for checking the LP. */
   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* retrieving the current parameter settings */
   SCIP_CALL( SCIPgetBoolParam(subproblem, "conflict/enable", &prevConfParam) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/disablecutoff", &prevCutoffParam) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/scaling", &prevLPScalingParam) );
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/initalgorithm", &prevInitAlgParam) );
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/resolvealgorithm", &prevResolveAlgParam) );
   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/alwaysgetduals", &prevDualParam) );
   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/scaleobj", &prevScaleObjParam) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxrounds", &prevPropMaxroundsParam) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxroundsroot", &prevPropMaxroundsRootParam) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "limits/bestsol", &prevLimitsBestsolParam) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "limits/maxorigsol", &prevLimitsMaxOrigSolParam) );


   /* If the solve has been stopped for the subproblem, then we need to restart it to complete the solve. The subproblem
    * is stopped when it is a MIP so that LP cuts and IP cuts can be generated. */
   if( SCIPgetStage(subproblem) == SCIP_STAGE_SOLVING )
   {
      /* the subproblem should be in probing mode. Otherwise, the eventhandler did not work correctly */
      assert( SCIPinProbing(subproblem) );

      /* the probing mode needs to be stopped so that the MIP can be solved */
      SCIP_CALL( SCIPendProbing(subproblem) );

      /* the problem was interrupted in the event handler, so SCIP needs to be informed that the problem is to be restarted */
      SCIP_CALL( SCIPrestartSolve(subproblem) );

      /* if the solve type is for CHECK, then the FEASIBILITY emphasis setting is used. */
      if( type == CHECK )
      {
         //SCIP_CALL( SCIPsetEmphasis(subproblem, SCIP_PARAMEMPHASIS_FEASIBILITY, TRUE) );

         SCIP_CALL( SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_FAST, TRUE) );

         /* the number of solution improvements is limited to try and prove feasibility quickly */
         /* NOTE: This should be a parameter */
         //SCIP_CALL( SCIPsetIntParam(subproblem, "limits/bestsol", 5) );
      }

      mipchecksolve = TRUE;
   }
   else if( solvemip )
   {
      /* if the MIP will be solved directly, then the probing mode needs to be skipped. This is done by dropping the
       * node focus event */
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(subproblem, MIPNODEFOCUS_EVENTHDLR_NAME);

      /* this may throw an error because the event is not being tracked. I am not sure about this yet. */
      SCIP_CALL(SCIPdropEvent(subproblem, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1));
   }
   else
   {
      /* if the problem is not in probing mode, then we need to solve the LP. That requires all methods that will
       * modify the structure of the problem need to be deactivated */

      /* modifying all of the parameters */

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

      SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", TRUE) );
      SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/scaleobj", FALSE) );

      //SCIPinfoMessage(subproblem, NULL, "Pricing problem %d\n", probnumber);
      SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );

      SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxrounds", 0) );
      SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxroundsroot", 0) );

      SCIP_CALL( SCIPsetIntParam(subproblem, "constraints/linear/propfreq", -1) );

      /* disable solution storage ! */
      SCIP_CALL( SCIPsetIntParam(subproblem, "limits/maxorigsol", 0) );

      //if( type != CHECK )
         //SCIP_CALL( SCIPsetIntParam(subproblem, "limits/nodes", 0) );
#ifdef SCIP_MOREDEBUG
      SCIP_CALL( SCIPsetBoolParam(subproblem, "display/lpinfo", TRUE) );
#endif
   }

#ifdef SCIP_MOREDEBUG
      SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_FULL) );
#endif

   SCIP_CALL( SCIPsolve(subproblem) );

   /* if the problem is still in the solving stage, then this indicates that the LP needs to be solved for the Benders'
    * cuts.
    * NOTE: This function needs an additional parameter to indicate whether the MIP should be solved directly without
    * also solving the LP first. */
   if( SCIPgetStage(subproblem) == SCIP_STAGE_SOLVING && !initialisation && !mipchecksolve )
   {
      /* constructing the LP */
      SCIP_CALL( SCIPconstructLP(subproblem, &cutoff) );

      /* When the node is focused, we are able to enter probing mode and solve the first LP. This will store the dual
       * information that is necessary to generate the Benders' cuts. */
      SCIP_CALL( SCIPstartProbing(subproblem) );

      /* solving the probing LP */
      SCIP_CALL( SCIPsolveProbingLP(subproblem, -1, &lperror, &cutoff) );
      assert(!lperror);

      /* the subproblem should be in probing mode. Otherwise, the eventhandler did not work correctly */
      assert( SCIPinProbing(subproblem) );
   }
   else if( SCIPgetStage(subproblem) != SCIP_STAGE_SOLVED )
      assert(initialisation || mipchecksolve || FALSE);

   if( SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE
      || (SCIPinProbing(subproblem) && SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE) )
      (*infeasible) = TRUE;
   else if( SCIPgetStatus(subproblem) != SCIP_STATUS_OPTIMAL && SCIPgetStatus(subproblem) != SCIP_STATUS_USERINTERRUPT
      && SCIPgetStatus(subproblem) != SCIP_STATUS_BESTSOLLIMIT )
      assert(FALSE);

   //SCIP_CALL( SCIPprintStatistics(subprob, NULL) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "display/lpinfo", FALSE) );
   /* resetting the parameter settings to the previous state */
   SCIP_CALL( SCIPsetPresolving(subproblem, SCIP_PARAMSETTING_DEFAULT, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_DEFAULT, TRUE) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "conflict/enable", prevConfParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/disablecutoff", prevCutoffParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/scaling", prevLPScalingParam) );
   SCIP_CALL( SCIPsetCharParam(subproblem, "lp/initalgorithm", prevInitAlgParam) );
   SCIP_CALL( SCIPsetCharParam(subproblem, "lp/resolvealgorithm", prevResolveAlgParam) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", prevDualParam) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/scaleobj", prevScaleObjParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxrounds", prevPropMaxroundsParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxroundsroot", prevPropMaxroundsRootParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "limits/bestsol", prevLimitsBestsolParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "limits/maxorigsol", prevLimitsMaxOrigSolParam) );

   return SCIP_OKAY;
}

/** frees the subproblems. */
SCIP_RETCODE SCIPbendersFreeSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnum             /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(benders->bendersfreesub != NULL || (benders->bendersfreesub == NULL && benders->benderssolvesub == NULL));
   assert(probnum >= 0 && probnum < benders->nsubproblems);

   if( benders->bendersfreesub != NULL )
      SCIP_CALL( benders->bendersfreesub(set->scip, benders, probnum) );
   else
   {
      SCIP* subproblem = SCIPbendersSubproblem(benders, probnum);

      if( SCIPbendersSubprobIsLP(benders, probnum) )
      {
         /* ending probing mode to reset the current node. The probing mode will be restarted at the next solve */
         SCIP_CALL( SCIPendProbing(subproblem) );
      }
      else
      {
         /* if the subproblems were solved as part of an enforcement stage, then they will still be in probing mode. The
          * probing mode must first be finished and then the problem can be freed */
         if( SCIPinProbing(subproblem) )
            SCIP_CALL( SCIPendProbing(subproblem) );

         SCIP_CALL( SCIPfreeTransform(subproblem) );
      }
   }

   return SCIP_OKAY;
}

/** checks the auxiliary variable value for optimality */
SCIP_RETCODE SCIPbendersCheckAuxiliaryVar(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the number of the pricing problem */
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

   SCIPsetDebugMsg(set, "Subproblem %d - Auxiliary Variable: %g Subproblem Objective: %g\n", probnumber, auxiliaryvarval,
      SCIPbendersGetSubprobObjval(benders, probnumber));

   /* if the value of the auxiliary variable in the master problem is greater or equal to the subproblem objective,
    * then a cut is not added by the subproblem.
    */
   if( SCIPsetIsFeasGE(set, auxiliaryvarval + soltol, SCIPbendersGetSubprobObjval(benders, probnumber)) )
      (*optimal) = TRUE;

   return SCIP_OKAY;
}

/** returns the value of the auxiliary variable value in a master problem solution */
SCIP_Real SCIPbendersGetAuxiliaryVarVal(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber          /**< the number of the pricing problem */
   )
{
   SCIP_VAR* auxiliaryvar;

   assert(benders != NULL);
   assert(set != NULL);

   auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, probnumber);
   assert(auxiliaryvar != NULL);

   return SCIPgetSolVal(set->scip, sol, auxiliaryvar);
}

/** returns the corresponding master or subproblem variable for the given variable.
 * This provides a call back for the variable mapping between the master and subproblems */
SCIP_VAR* SCIPbendersGetVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< the variable for which the corresponding variable is desired */
   int                   probnumber          /**< the problem number for the desired variable, -1 for the master problem */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(benders->bendersgetvar != NULL);

   /* if the variable name matches the auxiliary variable, then the master variable is returned as NULL */
   if( strstr(SCIPvarGetName(var), AUXILIARYVAR_NAME) != NULL
      || strstr(SCIPvarGetName(var), MW_AUXILIARYVAR_NAME) != NULL)
      return NULL;

   return benders->bendersgetvar(set->scip, benders, var, probnumber);
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

/** sets copy callback of benders */
void SCIPbendersSetCopy(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSCOPY ((*benderscopy))    /**< copy callback of benders */
   )
{
   assert(benders != NULL);

   benders->benderscopy = benderscopy;
}

/** sets destructor callback of benders */
void SCIPbendersSetFree(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSFREE ((*bendersfree))    /**< destructor of benders */
   )
{
   assert(benders != NULL);

   benders->bendersfree = bendersfree;
}

/** sets initialization callback of benders */
void SCIPbendersSetInit(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSINIT((*bendersinit))     /**< initialize benders */
   )
{
   assert(benders != NULL);

   benders->bendersinit = bendersinit;
}

/** sets deinitialization callback of benders */
void SCIPbendersSetExit(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSEXIT((*bendersexit))     /**< deinitialize benders */
   )
{
   assert(benders != NULL);

   benders->bendersexit = bendersexit;
}

/** sets presolving initialization callback of benders */
void SCIPbendersSetInitpre(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre))     /**< initialize benders */
   )
{
   assert(benders != NULL);

   benders->bendersinitpre = bendersinitpre;
}

/** sets presolving deinitialization callback of benders */
void SCIPbendersSetExitpre(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre))     /**< deinitialize benders */
   )
{
   assert(benders != NULL);

   benders->bendersexitpre = bendersexitpre;
}

/** sets solving process initialization callback of benders */
void SCIPbendersSetInitsol(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol))/**< solving process initialization callback of benders */
   )
{
   assert(benders != NULL);

   benders->bendersinitsol = bendersinitsol;
}

/** sets solving process deinitialization callback of benders */
void SCIPbendersSetExitsol(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol))/**< solving process deinitialization callback of benders */
   )
{
   assert(benders != NULL);

   benders->bendersexitsol = bendersexitsol;
}

/** sets solve callback of benders */
void SCIPbendersSetSolvesub(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub))/**< solving method for a Benders' decomposition subproblem */
   )
{
   assert(benders != NULL);

   benders->benderssolvesub = benderssolvesub;
}
/** sets post-solve callback of benders */
void SCIPbendersSetPostsolve(
   SCIP_BENDERS*         benders,            /**< benders */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve))/**< solving process deinitialization callback of benders */
   )
{
   assert(benders != NULL);

   benders->benderspostsolve = benderspostsolve;
}

/** sets free subproblem callback of benders */
void SCIPbendersSetFreesub(
   SCIP_BENDERS*         benders,            /**< benders */
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

/** gets the number of times, the benders was called and tried to find a variable with negative reduced costs */
int SCIPbendersGetNCalls(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->ncalls;
}

/** gets the number of optimality cuts found by the collection of Benders' decomposition subproblems */
int SCIPbendersGetNOptCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->noptcutsfound;
}

/** gets the number of feasibility cuts found by the collection of Benders' decomposition subproblems */
int SCIPbendersGetNFeasCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nfeascutsfound;
}

/** gets time in seconds used in this benders for setting up for next stages */
SCIP_Real SCIPbendersGetSetupTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return SCIPclockGetTime(benders->setuptime);
}

/** gets time in seconds used in this benders */
SCIP_Real SCIPbendersGetTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return SCIPclockGetTime(benders->bendersclock);
}

/** enables or disables all clocks of \p benders, depending on the value of the flag */
void SCIPbendersEnableOrDisableClocks(
   SCIP_BENDERS*         benders,            /**< the benders for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the benders be enabled? */
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

/** are Benders' cuts generated from the LP solutions? */
SCIP_Bool SCIPbendersCutLP(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutlp;
}

/** are Benders' cuts generated from the pseudo solutions? */
SCIP_Bool SCIPbendersCutPseudo(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutpseudo;
}

/** are Benders' cuts generated from the relaxation solutions? */
SCIP_Bool SCIPbendersCutRelaxation(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutrelax;
}

/** Adds a subproblem to the Benders' decomposition data */
SCIP_RETCODE SCIPbendersAddSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP*                 subproblem          /**< subproblem to be added to the data storage */
   )
{
   assert(benders != NULL);
   assert(subproblem != NULL);
   assert(benders->subproblems != NULL);
   assert(benders->addedsubprobs + 1 <= benders->nsubproblems);

   benders->subproblems[benders->addedsubprobs] = subproblem;

   benders->addedsubprobs++;

   return SCIP_OKAY;
}

/** Removes the subproblems from the Benders' decomposition data */
SCIP_RETCODE SCIPbendersRemoveSubproblems(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);
   assert(benders->subproblems != NULL);

   while( benders->addedsubprobs > 0 )
   {
      benders->addedsubprobs--;
      benders->subproblems[benders->addedsubprobs] = NULL;
   }

   return SCIP_OKAY;
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
void SCIPbendersSetSubprobObjval(
   SCIP_BENDERS*         benders,            /**< variable benders */
   SCIP_Real             objval,             /**< the objective function value for the subproblem */
   int                   probnumber          /**< the subproblem number */
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
SCIP_Real SCIPbendersGetSubprobObjval(
   SCIP_BENDERS*         benders,            /**< variable benders */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobobjval[probnumber];
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

/** returns the array of currently available Benders' cuts; active benders are in the first slots of the array */
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

/** sorts benders cuts by priorities */
void SCIPbendersSortBenderscuts(
   SCIP_BENDERS*         benders             /**< benders */
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

/** sorts benders cuts by name */
void SCIPbendersSortBenderscutsName(
   SCIP_BENDERS*         benders             /**< benders */
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

/** returns whether the Magnanti-Wong method has been applied */
SCIP_Bool SCIPbendersGetUseMagnantiWong(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->usemagnantiwong;
}

/* sets the flag indicating whether a subproblem is an LP. It is possible that this can change during the solving
 * process. One example is when the three-phase method is employed, where the first phase solves the of both the master
 * and subproblems and by the third phase the integer subproblem is solved. */
void SCIPbendersSetSubprobIsLP(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             islp                /**< flag to indicate whether the subproblem is an LP */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   if( islp && !benders->subprobislp[probnumber] )
      benders->nlpsubprobs++;
   else if( !islp && benders->subprobislp[probnumber] )
      benders->nlpsubprobs--;

   benders->subprobislp[probnumber] = islp;

   assert(benders->nlpsubprobs >= 0 && benders->nlpsubprobs <= benders->nsubproblems);
}

/* returns whether the subproblem is an LP. This means that the dual solution can be trusted. */
SCIP_Bool SCIPbendersSubprobIsLP(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobislp[probnumber];
}

/** returns the number of subproblems that are LPs */
int SCIPbendersGetNLPSubprobs(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nlpsubprobs;
}

/** changes all of the master problem variables in the given subproblem to continuous */
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

   /* if the benders is a copy, then the actions of this function should have already been performed. As such, the
    * MastervarsCont flag is being set to TRUE */
   //if( benders->iscopy )
      //SCIP_CALL( SCIPbendersSetMastervarsCont(benders, probnumber, TRUE) );

   /* only set the master problem variable to continuous if they have not already been changed. */
   if( !SCIPbendersGetMastervarsCont(benders, probnumber) )
   {
      /* retreiving the variable data */
      SCIP_CALL( SCIPgetVarsData(subproblem, &vars, NULL, &nbinvars, &nintvars, &nimplvars, NULL) );

      origintvars = nbinvars + nintvars + nimplvars;

      chgvarscount = 0;

      /* looping over all integer variables to change the master variables to continuous */
      i = 0;
      while( i < nbinvars + nintvars + nimplvars )
      {
         if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS
            && SCIPbendersGetVar(benders, set, vars[i], -1) != NULL )
         {
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
         SCIP_CALL( initialiseLPSubproblem(benders, probnumber) );
         SCIPbendersSetSubprobIsLP(benders, probnumber, TRUE);
      }

      SCIP_CALL( SCIPbendersSetMastervarsCont(benders, probnumber, TRUE) );
   }

   return SCIP_OKAY;
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
         SCIP_CALL( SCIPendProbing(SCIPbendersSubproblem(benders, probnumber)) );

      SCIPbendersSetSubprobIsLP(benders, probnumber, FALSE);
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
