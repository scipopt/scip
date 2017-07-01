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

#include "scip/struct_benders.h"
#include "scip/struct_benderscut.h"

#include "scip/benderscut.h"
#include "scip/misc_benders.h"

/* Defaults for parameters */
#define SCIP_DEFAULT_USEMW               FALSE  /** Should Magnanti-Wong cut strengthening be used? */
#define SCIP_DEFAULT_COMPUTERELINT       FALSE  /** Should the relative interior point be computed? */
#define SCIP_DEFAULT_UPDATEFACTOR          0.5  /** The factor to update the relative interior point? */

#define MW_AUXILIARYVAR_NAME "##MWauxiliaryvar##"

/* event handler properties */
#define EVENTHDLR_NAME         "bendersnodefocus"
#define EVENTHDLR_DESC         "node focus event handler for Benders' decomposition"

/* ---------------- Callback methods of event handler ---------------- */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBendersNodefocus)
{  /*lint --e{715}*/

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

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
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL));

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolBendersNodefocus)
{
   assert(scip != NULL);

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1));

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
   SCIP_CALL( SCIPbendersExecSubproblemSolve(benders, set, benders->relintsol, probnum, TRUE, &infeasible) );
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
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "auxiliaryvar_%d", i );
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
   }

   /* if the copy was valid, then the Benders cuts are copied. */
   if( valid )
   {
      targetbenders = SCIPsetFindBenders(set, SCIPbendersGetName(benders));

      /* the flag is set to indicate that the  */
      targetbenders->iscopy = TRUE;

      /* calling the copy method for the Benders' cuts */
      SCIPbendersSortBenderscuts(benders);
      for( i = 0; i < benders->nbenderscuts; i++ )
      {
         SCIP_CALL( SCIPbenderscutCopyInclude(benders->benderscuts[i], set) );
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
   int                   nsubproblems,       /**< the number subproblems used in this decomposition */
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

   int i;

   assert(benders != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   SCIP_ALLOC( BMSallocMemory(benders) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->desc, desc, strlen(desc)+1) );
   (*benders)->priority = priority;
   (*benders)->nsubproblems = nsubproblems;
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
   (*benders)->initialized = FALSE;
   (*benders)->iscopy = FALSE;
   (*benders)->maxlpiterfactor = 1.0;
   (*benders)->updatefactor = SCIP_DEFAULT_UPDATEFACTOR;
   (*benders)->coreptupdated = FALSE;
   (*benders)->relintsol = NULL;
   (*benders)->currentsol = NULL;
   (*benders)->addedsubprobs = 0;

   (*benders)->benderscuts = NULL;
   (*benders)->nbenderscuts = 0;
   (*benders)->benderscutssize = 0;
   (*benders)->benderscutssorted = FALSE;
   (*benders)->benderscutsnamessorted = FALSE;

   /* allocating memory for the subproblems arrays */
   SCIP_ALLOC( BMSallocMemoryArray(&(*benders)->subproblems, (*benders)->nsubproblems) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*benders)->mwauxiliaryvars, (*benders)->nsubproblems) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*benders)->auxiliaryvars, (*benders)->nsubproblems) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*benders)->subprobobjval, (*benders)->nsubproblems) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*benders)->bestsubprobobjval, (*benders)->nsubproblems) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*benders)->subprobislp, (*benders)->nsubproblems) );

   for( i = 0; i < (*benders)->nsubproblems; i++ )
   {
      (*benders)->subproblems[i] = NULL;
      (*benders)->mwauxiliaryvars[i] = NULL;
      (*benders)->auxiliaryvars[i] = NULL;
      (*benders)->subprobobjval[i] = SCIPsetInfinity(set);
      (*benders)->bestsubprobobjval[i] = SCIPsetInfinity(set);
      (*benders)->subprobislp[i] = FALSE;
   }

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

   SCIP_CALL( SCIPfreeSol(set->scip, &(*benders)->relintsol) );

   BMSfreeMemoryArray(&(*benders)->subprobislp);
   BMSfreeMemoryArray(&(*benders)->bestsubprobobjval);
   BMSfreeMemoryArray(&(*benders)->subprobobjval);
   BMSfreeMemoryArray(&(*benders)->auxiliaryvars);
   BMSfreeMemoryArray(&(*benders)->mwauxiliaryvars);
   BMSfreeMemoryArray(&(*benders)->subproblems);

   SCIPclockFree(&(*benders)->bendersclock);
   SCIPclockFree(&(*benders)->setuptime);
   BMSfreeMemoryArray(&(*benders)->name);
   BMSfreeMemoryArray(&(*benders)->desc);
   BMSfreeMemory(benders);

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
   int nbinvars;
   int nintvars;
   int nsubproblems;
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* creating all subproblems */
   for( i = 0; i < nsubproblems; i++ )
   {
      /* calling the create subproblem call back method */
      SCIP_CALL( benders->benderscreatesub(set->scip, benders, i) );

      subproblem = SCIPbendersSubproblem(benders, i);

      assert(subproblem != NULL);

      /* getting the number of integer and binary variables to determine the problem type */
      SCIP_CALL( SCIPgetVarsData(subproblem, NULL, NULL, &nbinvars, &nintvars, NULL, NULL) );

      /* if there are no binary and integer variables, then the subproblem is an LP.
       * In this case, the SCIP instance is put into probing mode via the use of an event handler. */
      if( nbinvars == 0 && nintvars == 0 )
      {
         SCIP_EVENTHDLR* eventhdlr;
         SCIP_Bool infeasible;
         SCIP_Bool cutoff;

         benders->subprobislp[i] = TRUE;

         /* if the user has not implemented a solve subproblem callback, then the subproblem solves are performed
          * internally. To be more efficient the subproblem is put into probing mode. */
         if( benders->benderssolvesub == NULL )
         {

            /* include event handler into SCIP */
            SCIP_CALL( SCIPincludeEventhdlrBasic(subproblem, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
                  eventExecBendersNodefocus, NULL) );
            SCIP_CALL( SCIPsetEventhdlrInitsol(subproblem, eventhdlr, eventInitsolBendersNodefocus) );
            SCIP_CALL( SCIPsetEventhdlrExitsol(subproblem, eventhdlr, eventExitsolBendersNodefocus) );
            assert(eventhdlr != NULL);

            /* Getting the problem into the right SCIP stage for solving */
            SCIP_CALL( SCIPbendersSolveSubproblem(benders, i, &infeasible) );

            /* Constructing the LP that can be solved in later iterations */
            SCIP_CALL( SCIPconstructLP(subproblem, &cutoff) );

            /* starting probing mode */
            SCIP_CALL( SCIPstartProbing(subproblem) );
         }
      }
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

   /* adding the auxiliary variables to the master problem */
   SCIP_CALL( addAuxiliaryVariablesToMaster(set->scip, benders) );

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

/** solves the subproblem using the current master problem solution. */
/*  TODO: consider allowing the possibility to pass solution information back from the subproblems instead of the scip
 *  instance. This would allow the use of different solvers for the subproblems, more importantly allowing the use of an
 *  LP solver for LP subproblems. */
SCIP_RETCODE SCIPbendersExec(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RESULT*          result,             /**< result of the pricing process */
   SCIP_Bool             check               /**< is the execution method called as a check. i.e. no cuts are required */
   )
{
   SCIP_BENDERSCUT** benderscuts;
   int nsubproblems;
   int nbenderscuts;
   int ncutloops;
   int i;
   int j;
   int k;
   SCIP_Bool infeasible = FALSE;

   assert(benders != NULL);
   assert(result != NULL);
   assert(benders->bendersexec != NULL);

   /* stating that the core point has not been updated for this iteration */
   benders->coreptupdated = FALSE;

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* sets the stored objective function values of the subproblems to infinity */
   resetSubproblemObjectiveValue(benders);

   SCIP_CALL( benders->bendersexec(set->scip, benders) );

   *result = SCIP_DIDNOTRUN;

   if( check && sol == NULL )
   {
      infeasible = TRUE;
   }
   else
   {
      /* solving each of the subproblems for Benders decomposition */
      /* TODO: ensure that the each of the subproblems solve and update the parameters with the correct return values */
      for( i = 0; i < nsubproblems; i++ )
      {
         SCIP_Bool subinfeas = FALSE;

         SCIP_CALL( SCIPbendersExecSubproblemSolve(benders, set, sol, i, FALSE, &subinfeas) );

         infeasible = infeasible || subinfeas;

         /* if the subproblems are being solved as part of the conscheck, then we break once an infeasibility is found.
          * The result pointer is set to infeasible and the execution is halted. */
         if( check )
         {
            /* if the subproblem is feasible, then it is necessary to update the value of the auxiliary variable to the
             * objective function value of the subproblem. */
            if( !subinfeas )
            {
               SCIP_Bool optimal;

               SCIP_CALL( SCIPbendersCheckAuxiliaryVar(benders, set, sol, i, &optimal) );

               infeasible = infeasible || !optimal;
            }
         }
      }
   }


   /* Preparing the data for the Magnanti-Wong technique */
   if( !infeasible )
      prepareMagnantiWongTechnique(set->scip, benders, sol);


   /* the result flag is set to FEASIBLE as default. If a cut is added, then this is changed to CONS_ADDED */
   (*result) = SCIP_FEASIBLE;

   /* Generating cuts for the subproblems. */
   benderscuts = SCIPbendersGetBenderscuts(benders);
   nbenderscuts = SCIPbendersGetNBenderscuts(benders);

   /* if the Magnanti-Wong enhancement technique is used, then an additional cutloop is performed. */
   /* TODO: create a parameter to modify the behaviour of the cutloops. */
   ncutloops = 1;
   if( benders->usemagnantiwong && !infeasible )
      ncutloops = 2;

   /* This is done in two loops. The first is by subproblem and the second is by cut type. */
   for( i = 0; i < nsubproblems; i++ )
   {
      for( j = 0; j < ncutloops; j++ )
      {
         /* Performing cut enhancement techniques */
         if( j == 1 )
            SCIP_CALL( performMagnantiWongTechnique(benders, set, i) );

         for( k = 0; k < nbenderscuts; k++ )
         {
            assert(benderscuts[k] != NULL);

            SCIP_CALL( SCIPbenderscutExec(benderscuts[k], set, benders, sol, i, result) );
         }
      }
   }

   if( check )
   {
      /* if the subproblems are being solved as part of conscheck, then the results flag must be returned after the solving
       * has completed. No cut is generated during conscheck. */
      if( infeasible )
         (*result) = SCIP_INFEASIBLE;
      else
         (*result) = SCIP_FEASIBLE;
   }

   /* calling the post-solve call back for the Benders' decomposition algorithm. This allows the user to work directly
    * with the solved subproblems and the master problem */
   if( benders->benderspostsolve != NULL )
      SCIP_CALL( benders->benderspostsolve(set->scip, benders, infeasible) );

   /* freeing the subproblems after the cuts are generated */
   for( i = 0; i < benders->nsubproblems; i++ )
   {
      /* cleaning up after performing Magnanti-Wong */
      SCIP_CALL( cleanupMagnantiWong(benders, set, i) );

      SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, i) );
   }

   /* increment the number of calls to the Benders' decomposition subproblem solve */
   benders->ncalls++;

   return SCIP_OKAY;
}

/** solves the subproblems. */
SCIP_RETCODE SCIPbendersExecSubproblemSolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnum,            /**< the subproblem number */
   SCIP_Bool             enhancement,        /**< is the solve performed as part of and enhancement? */
   SCIP_Bool*            infeasible          /**< returns whether the current subproblem is infeasible */
   )
{
   SCIP* subproblem;
   SCIP_SOL* bestsol;

   assert(benders != NULL);
   assert(probnum >= 0 && probnum < benders->nsubproblems);

   /* if the subproblem solve callback is implemented, then that is used instead of the default setup */
   if( benders->benderssolvesub != NULL)
      SCIP_CALL( benders->benderssolvesub(set->scip, benders, sol, probnum, infeasible) );
   else
   {
      /* setting up the subproblem */
      SCIP_CALL( SCIPbendersSetupSubproblem(benders, set, sol, probnum) );

      /* solving the subproblem */
      if( SCIPbendersSubprobIsLP(benders, probnum) )
         SCIP_CALL( SCIPbendersSolveSubproblemLP(benders, probnum, infeasible) );
      else
         SCIP_CALL( SCIPbendersSolveSubproblem(benders, probnum, infeasible) );
   }

   subproblem = SCIPbendersSubproblem(benders, probnum);

   bestsol = SCIPgetBestSol(subproblem);

   if( !enhancement )
   {
      if( SCIPgetStatus(subproblem) == SCIP_STATUS_OPTIMAL || SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_OPTIMAL )
         SCIPbendersSetSubprobObjval(benders, SCIPgetSolTransObj(subproblem, bestsol), probnum);
      else if( SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE
         || SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE )
         SCIPbendersSetSubprobObjval(benders, SCIPinfinity(set->scip), probnum);
      else
         assert(FALSE);
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

   subproblem = SCIPbendersSubproblem(benders, probnum);

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

         assert(fixed);
         assert(!infeasible);
      }
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
SCIP_RETCODE SCIPbendersSolveSubproblem(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible          /**< a flag to indicate whether all subproblems are feasible */
   )
{
   SCIP* subproblem;

   /* previous parameter settings */
   int prevCutoffParam;
   int prevPropMaxroundsParam;
   int prevPropMaxroundsRootParam;
   char prevInitAlgParam;
   char prevResolveAlgParam;
   SCIP_Bool prevConfParam;
   SCIP_Bool prevDualParam;

   assert(benders != NULL);
   assert(infeasible != NULL);

   (*infeasible) = FALSE;


   /* TODO: This should be solved just as an LP, so as a MIP. There is too much overhead with the MIP.
    * Need to change status check for checking the LP. */
   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* modifying all of the parameters */

   /* Do we have to disable presolving? If yes, we have to store all presolving parameters. */
   SCIPsetPresolving(subproblem, SCIP_PARAMSETTING_OFF, TRUE);

   /* Disabling heuristics so that the problem is not trivially solved */
   SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_OFF, TRUE);

   /* store parameters that are changed for the generation of the subproblem cuts */
   SCIP_CALL( SCIPgetBoolParam(subproblem, "conflict/enable", &prevConfParam) );
   SCIPsetParam(subproblem, "conflict/enable", FALSE);

   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/disablecutoff", &prevCutoffParam) );
   SCIPsetIntParam(subproblem, "lp/disablecutoff", 1);

   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/initalgorithm", &prevInitAlgParam) );
   SCIPsetCharParam(subproblem, "lp/initalgorithm", 'd');
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/resolvealgorithm", &prevResolveAlgParam) );
   SCIPsetCharParam(subproblem, "lp/resolvealgorithm", 'd');

   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/alwaysgetduals", &prevDualParam) );
   SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", TRUE);

   //SCIPinfoMessage(subproblem, NULL, "Pricing problem %d\n", probnumber);
   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
   //SCIP_CALL( SCIPsetBoolParam(subproblem, "display/lpinfo", TRUE) );

   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxrounds", &prevPropMaxroundsParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxrounds", 0) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxroundsroot", &prevPropMaxroundsRootParam) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxroundsroot", 0) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "constraints/linear/propfreq", -1) );

   SCIP_CALL( SCIPsolve(subproblem) );

   if( SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE )
      (*infeasible) = TRUE;
   else if( SCIPgetStatus(subproblem) != SCIP_STATUS_OPTIMAL && SCIPgetStatus(subproblem) != SCIP_STATUS_USERINTERRUPT )
      assert(FALSE);

   //SCIP_CALL( SCIPprintStatistics(subprob, NULL) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
   /* resetting the parameter settings to the previous state */
   SCIPsetPresolving(subproblem, SCIP_PARAMSETTING_DEFAULT, TRUE);
   SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_DEFAULT, TRUE);
   SCIPsetBoolParam(subproblem, "conflict/enable", prevConfParam);
   SCIPsetIntParam(subproblem, "lp/disablecutoff", prevCutoffParam);
   SCIPsetCharParam(subproblem, "lp/initalgorithm", prevInitAlgParam);
   SCIPsetCharParam(subproblem, "lp/resolvealgorithm", prevResolveAlgParam);
   SCIPsetBoolParam(subproblem, "misc/alwaysgetduals", prevDualParam);
   SCIPsetIntParam(subproblem, "propagating/maxrounds", prevPropMaxroundsParam);
   SCIPsetIntParam(subproblem, "propagating/maxroundsroot", prevPropMaxroundsRootParam);

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
         /* ending probing mode to reset the current node */
         SCIP_CALL( SCIPendProbing(subproblem) );

         /* starting probing mode to fix variables for the subproblem */
         SCIP_CALL( SCIPstartProbing(subproblem) );
      }
      else
         SCIP_CALL( SCIPfreeTransform(subproblem) );
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
   SCIP_VAR* auxiliaryvar;
   SCIP_Real auxiliaryvarval;
   SCIP_Real soltol;

   assert(benders != NULL);
   assert(set != NULL);

   (*optimal) = FALSE;

   auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, probnumber);
   assert(auxiliaryvar != NULL);

   auxiliaryvarval = SCIPgetSolVal(set->scip, sol, auxiliaryvar);

   SCIP_CALL( SCIPsetGetRealParam(set, "benders/solutiontol", &soltol) );

   SCIPsetDebugMsg(set, "Auxiliary Variable: %g Subproblem Objective: %g\n", auxiliaryvarval,
      SCIPbendersGetSubprobObjval(benders, probnumber));

   /* if the value of the auxiliary variable in the master problem is greater or equal to the subproblem objective,
    * then a cut is not added by the subproblem.
    */
   if( SCIPsetIsGE(set, auxiliaryvarval + soltol, SCIPbendersGetSubprobObjval(benders, probnumber)) )
      (*optimal) = TRUE;

   return SCIP_OKAY;
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
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */

   assert(benders != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(benders->bendersgetvar != NULL);

   (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "t_%s", MW_AUXILIARYVAR_NAME );
   /* if the variable name matches the auxiliary variable, then the master variable is returned as NULL */
   if( strcmp(SCIPvarGetName(var), varname) == 0 || strcmp(SCIPvarGetName(var), MW_AUXILIARYVAR_NAME) == 0)
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
