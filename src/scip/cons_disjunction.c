/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file cons_disjunction.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for disjunction constraints
 * @author Stefan Heinz
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_disjunction.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "disjunction"
#define CONSHDLR_DESC          "disjunction of constraints (or(cons1, cons2, ..., consn))"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -950000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -900000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in
                                         *   (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP



/*
 * Data structures
 */

/** constraint data for disjunction constraints */
struct SCIP_ConsData
{
   SCIP_CONS**           conss;              /**< constraints in disjunction */
   int                   consssize;          /**< size of conss array */
   int                   nconss;             /**< number of constraints in disjunction */
};




/*
 * Local methods
 */

/** creates disjunction constraint data, captures initial constraints of disjunction */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to constraint data */
   SCIP_CONS**           conss,              /**< initial constraint in disjunction */
   int                   nconss              /**< number of initial constraints in disjunction */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   if( nconss > 0 )
   {
      assert(conss != NULL);
      
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->conss, conss, nconss) );

      (*consdata)->consssize = nconss;
      (*consdata)->nconss = nconss;

      /* we need to capture the constraints to avoid that SCIP deletes them since they are not (yet) added to the
       * problem 
       */ 
      if( SCIPisTransformed(scip) )
      {
         SCIP_CALL( SCIPtransformConss(scip, nconss, (*consdata)->conss, (*consdata)->conss) );
      }
      else
      {
         int c;
      
         for( c = 0; c < nconss; ++c )
         {
            assert(conss[c] != NULL);
            SCIP_CALL( SCIPcaptureCons(scip, conss[c]) );
            
         }
      }
   }
   else
   {
      (*consdata)->conss = NULL;
      (*consdata)->consssize = 0;
      (*consdata)->nconss = 0;
   }

   return SCIP_OKAY;
}

/** frees constraint data and releases all constraints in disjunction */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to constraint data */
   )
{
   int c;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release constraints */
   for( c = 0; c < (*consdata)->nconss; ++c )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*consdata)->conss[c]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->conss, (*consdata)->consssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** adds constraint to disjunction */
static
SCIP_RETCODE consdataAddCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_CONS*            cons                /**< constraint to add to the disjunction */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(cons != NULL);

   /* get memory for additional constraint */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &consdata->conss, &consdata->consssize, consdata->nconss+1) );
   assert(consdata->conss != NULL);
   assert(consdata->nconss < consdata->consssize);

   /* insert constraint in array */
   consdata->conss[consdata->nconss] = cons;
   consdata->nconss++;

   /* capture constraint */
   SCIP_CALL( SCIPcaptureCons(scip, cons) );

   return SCIP_OKAY;
}

/** branches disjuctive constraint */
static
SCIP_RETCODE branchCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< active disjunction constraint */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;
   SCIP_NODE* child;
   SCIP_Real estimate;
   int nconss;
   int i;

   assert(result != NULL);

   /* can branch modifiable constraint */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conss = consdata->conss;
   assert(conss != NULL);

   nconss = consdata->nconss;
   assert(nconss > 0);

   estimate =  SCIPgetLocalTransEstimate(scip);

   /* add all inactive constraints to local subproblem */
   for( i = 0; i < nconss; ++i )
   {
      /* create the branch-and-bound tree child nodes of the current node */
      SCIP_CALL( SCIPcreateChild(scip, &child, 0.0, estimate) );

      /* add constraints to nodes */
      SCIP_CALL( SCIPaddConsNode(scip, child, conss[i], NULL) );

      /* remove disjunction constraint, from child node */
      SCIP_CALL( SCIPdelConsNode(scip, child, cons) );
   }
   
   /* reset constraint age */
   SCIP_CALL( SCIPresetConsAge(scip, cons) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** checks disjunction constraints if at least one is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< active disjunction constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;
   int nconss;
   int i;

   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conss = consdata->conss;
   assert(conss != NULL);

   nconss = consdata->nconss;
   assert(nconss > 0);

   *result = SCIP_INFEASIBLE;
   
   /* check all constraints */
   for( i = 0; i < nconss && *result != SCIP_FEASIBLE; ++i )
   {
      SCIP_CALL( SCIPcheckCons(scip, conss[i], sol, checkintegrality, checklprows, printreason, result) );
      assert(*result == SCIP_FEASIBLE || *result == SCIP_INFEASIBLE);
   }
   
   return SCIP_OKAY;
}

/** propagation method for disjunction constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONS** conss;
   int nconss;
   int c;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conss = consdata->conss;
   assert(conss != NULL);

   nconss = consdata->nconss;
   assert(nconss > 1);
   
   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsActive(conss[c]) )
      {
         (*ndelconss)++;
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         break;
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyDisjunction)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrDisjunction(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#define consFreeDisjunction NULL


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitDisjunction NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitDisjunction NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreDisjunction NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreDisjunction NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolDisjunction NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#define consExitsolDisjunction NULL


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteDisjunction)
{  /*lint --e{715}*/
   SCIP_CALL( consdataFree(scip, consdata) );
   
   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransDisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   /* get constraint data of source constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->conss, sourcedata->nconss) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
#define consInitlpDisjunction NULL


/** separation method of constraint handler for LP solutions */
#define consSepalpDisjunction NULL


/** separation method of constraint handler for arbitrary primal solutions */
#define consSepasolDisjunction NULL


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpDisjunction)
{  /*lint --e{715}*/
   int c;

   *result = SCIP_FEASIBLE;
   
   for( c = 0; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      /* check the disjunction */
      SCIP_CALL( checkCons(scip, conss[c], NULL, FALSE, FALSE, FALSE, result) );

      if( *result == SCIP_INFEASIBLE )
      {
         SCIP_CALL( branchCons(scip, conss[c], result) );
      }
   }
   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsDisjunction)
{  /*lint --e{715}*/
   int c;

   *result = SCIP_FEASIBLE;
   
   for( c = 0; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      /* check the disjunction */
      SCIP_CALL( checkCons(scip, conss[c], NULL, FALSE, FALSE, FALSE, result) );

      if( *result == SCIP_INFEASIBLE )
      {
         SCIP_CALL( branchCons(scip, conss[c], result) );
      }
   }
   
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckDisjunction)
{  /*lint --e{715}*/
   int c;

   *result = SCIP_FEASIBLE;
   
   for( c = 0; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      /* check the disjunction */
      SCIP_CALL( checkCons(scip, conss[c], sol, checkintegrality, checklprows, printreason, result) );
      assert(*result == SCIP_FEASIBLE || *result == SCIP_INFEASIBLE);
   }
   
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropDisjunction)
{  /*lint --e{715}*/
   int ndelconss;
   int c;
   
   ndelconss = 0;
   
   for( c = 0; c < nconss; ++c )
   {
      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, conss[c], &ndelconss) );
   }

   /* adjust result code */
   if( ndelconss > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolDisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int oldndelconss;
   int c;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;
   oldndelconss = *ndelconss;

   /* all disjunction constraints with one constraint can be replaced with that corresponding constraint */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPconsIsModifiable(conss[c]) && consdata->nconss == 1 )
      {
         /* add constaint to the problem */
         if( !SCIPconsIsActive(consdata->conss[0]) )
         {
            SCIP_CALL( SCIPaddCons(scip, consdata->conss[0]) );

            /* release constraint from the disjunction constraint */
            SCIP_CALL( SCIPreleaseCons(scip, &consdata->conss[0]) );
         }
         
         /* remove disjunction constraint */
         SCIP_CALL( SCIPdelCons(scip, conss[0]) );

         *result = SCIP_SUCCESS;
      }

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, conss[c], ndelconss) );
   }

   if( *ndelconss > oldndelconss )
      *result = SCIP_SUCCESS;
   
   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
#define consRespropDisjunction NULL


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockDisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* lock sub constraints */
   for( c = 0; c < consdata->nconss; ++c )
   {
      SCIP_CALL( SCIPaddConsLocks(scip, consdata->conss[c], nlockspos, nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveDisjunction NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveDisjunction NULL


/** constraint enabling notification method of constraint handler */
#define consEnableDisjunction NULL


/** constraint disabling notification method of constraint handler */
#define consDisableDisjunction NULL

/** variable deletion method of constraint handler */
#define consDelVarsDisjunction NULL


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintDisjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPinfoMessage(scip, file, "disjunction(");
      
   for( i = 0; i < consdata->nconss; ++i )
   {
      if( i > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "<%s>", SCIPconsGetName(consdata->conss[i]));
   }
   SCIPinfoMessage(scip, file, ")");
   
   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyConjuction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONS** sourceconss;
   SCIP_CONS** conss;
   int nconss;
   int c;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   nconss = sourcedata->nconss;

   if( nconss == 0 && !SCIPconsIsModifiable(sourcecons) )
   {
      *valid = TRUE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &conss, nconss) );
   sourceconss = sourcedata->conss;

   /* copy each constraint one by one */
   for( c = 0; c < nconss && (*valid); ++c )
   {
      SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourceconss[c], &conss[c], sourceconshdlr, 
            varmap, consmap, SCIPconsGetName(sourceconss[c]),  
            SCIPconsIsInitial(sourceconss[c]), SCIPconsIsSeparated(sourceconss[c]), SCIPconsIsEnforced(sourceconss[c]),
            SCIPconsIsChecked(sourceconss[c]), SCIPconsIsPropagated(sourceconss[c]),
            SCIPconsIsLocal(sourceconss[c]), SCIPconsIsModifiable(sourceconss[c]),
            SCIPconsIsDynamic(sourceconss[c]), SCIPconsIsRemovable(sourceconss[c]), SCIPconsIsStickingAtNode(sourceconss[c]), 
            global, valid) );
      assert(!(*valid) || conss[c] != NULL);
   }

   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsDisjunction(scip, cons, name, nconss, conss, 
            initial, enforce, check, local, modifiable, dynamic) );
   }
   
   SCIPfreeBufferArray(scip, &conss);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
#define consParseConjuction NULL



/*
 * constraint specific interface methods
 */

/** creates the handler for disjunction constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrDisjunction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create disjunction constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyDisjunction,
         consFreeDisjunction, consInitDisjunction, consExitDisjunction, 
         consInitpreDisjunction, consExitpreDisjunction, consInitsolDisjunction, consExitsolDisjunction,
         consDeleteDisjunction, consTransDisjunction, consInitlpDisjunction,
         consSepalpDisjunction, consSepasolDisjunction, consEnfolpDisjunction, consEnfopsDisjunction, 
         consCheckDisjunction, consPropDisjunction, consPresolDisjunction, consRespropDisjunction, consLockDisjunction,
         consActiveDisjunction, consDeactiveDisjunction, 
         consEnableDisjunction, consDisableDisjunction,
         consDelVarsDisjunction, consPrintDisjunction, consCopyConjuction, consParseConjuction,
         conshdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a disjunction constraint */
SCIP_RETCODE SCIPcreateConsDisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nconss,             /**< number of initial constraints in disjunction */
   SCIP_CONS**           conss,              /**< initial constraint in disjunction */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic             /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are seperated as constraints. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the disjunction constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("disjunction constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, conss, nconss) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, FALSE, enforce, check, FALSE,
         local, modifiable, dynamic, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** adds constraint to the disjunction of constraints */
SCIP_RETCODE SCIPaddConsElemDisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< disjunction constraint */
   SCIP_CONS*            addcons             /**< additional constraint in disjunction */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(addcons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a disjunction constraint\n");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( consdataAddCons(scip, consdata, addcons) );

   return SCIP_OKAY;
   
}

