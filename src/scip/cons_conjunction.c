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

/**@file   cons_conjunction.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for conjunction constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_conjunction.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "conjunction"
#define CONSHDLR_DESC          "conjunction of constraints"
#define CONSHDLR_SEPAPRIORITY   +000000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   +900000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -900000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP



/*
 * Data structures
 */

/** constraint data for conjunction constraints */
struct SCIP_ConsData
{
   SCIP_CONS**           conss;              /**< constraints in conjunction */
   int                   consssize;          /**< size of conss array */
   int                   nconss;             /**< number of constraints in conjunction */
};




/*
 * Local methods
 */

/** creates conjunction constraint data, captures initial constraints of conjunction */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to constraint data */
   SCIP_CONS**           conss,              /**< initial constraint in conjunction */
   int                   nconss              /**< number of initial constraints in conjunction */
   )
{
   assert(consdata != NULL);
   
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   if( nconss > 0 )
   {
      int c;

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->conss, conss, nconss) );
      (*consdata)->consssize = nconss;
      (*consdata)->nconss = nconss;
      for( c = 0; c < nconss; ++c )
      {
         if( SCIPconsIsInitial(conss[c]) || SCIPconsIsChecked(conss[c]) )
         {
            SCIPerrorMessage("constraints in a conjunction must not be initial or checked\n");
            return SCIP_INVALIDDATA;
         }
         SCIP_CALL( SCIPcaptureCons(scip, conss[c]) );
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

/** frees constraint data and releases all constraints in conjunction */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to constraint data */
   )
{
   int c;

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

/** adds constraint to conjunction */
static
SCIP_RETCODE consdataAddCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_CONS*            cons                /**< constraint to add to the conjunction */
   )
{
   assert(consdata != NULL);

   if( SCIPconsIsInitial(cons) || SCIPconsIsChecked(cons) )
   {
      SCIPerrorMessage("constraints in a conjunction must not be initial or checked\n");
      return SCIP_INVALIDDATA;
   }

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

/** adds all constraints in conjunction constraints to the problem; disables unmodifiable conjunction constraints */
static
SCIP_RETCODE addAllConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< active conjunction constraints */
   int                   nconss,             /**< number of active conjunction constraints */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSDATA* consdata;
   int c;
   int i;

   assert(result != NULL);

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* add all inactive constraints to local subproblem */
      for( i = 0; i < consdata->nconss; ++i )
      {
         if( !SCIPconsIsActive(consdata->conss[i]) )
         {
            SCIPdebugMessage("adding constraint <%s> from add conjunction <%s>\n",
               SCIPconsGetName(consdata->conss[i]), SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPaddConsLocal(scip, consdata->conss[i], NULL) );
            *result = SCIP_CONSADDED;
         }
      }

      /* disable conjunction constraint, if it is unmodifiable */
      if( !SCIPconsIsModifiable(conss[c]) )
      {
         SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** checks all constraints in conjunction constraints for feasibility */
static
SCIP_RETCODE checkAllConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< active conjunction constraints */
   int                   nconss,             /**< number of active conjunction constraints */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSDATA* consdata;
   int c;
   int i;

   assert(result != NULL);

   for( c = 0; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* check all constraints */
      for( i = 0; i < consdata->nconss && *result == SCIP_FEASIBLE; ++i )
      {
         SCIP_CALL( SCIPcheckCons(scip, consdata->conss[i], sol, checkintegrality, checklprows, printreason, result) );
      }

      /* disable conjunction constraint, if it is unmodifiable */
      if( !SCIPconsIsModifiable(conss[c]) )
      {
         SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
      }
   }
   
   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyConjunction)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrConjunction(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}
#else
#define conshdlrCopyConjunction NULL
#endif


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#define consFreeConjunction NULL


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitConjunction NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitConjunction NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreConjunction NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreConjunction NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolConjunction NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#define consExitsolConjunction NULL


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteConjunction)
{  /*lint --e{715}*/
   SCIP_CALL( consdataFree(scip, consdata) );
   
   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransConjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   int c;

   /* create constraint data for target constraint */
   SCIP_CALL( SCIPallocBlockMemory(scip, &targetdata) );

   /* get constraint data of source constraint */
   sourcedata = SCIPconsGetData(sourcecons);

   if( sourcedata->nconss > 0 )
   {
      targetdata->consssize = sourcedata->nconss;
      targetdata->nconss = sourcedata->nconss;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &targetdata->conss, targetdata->consssize) );
      for( c = 0; c < sourcedata->nconss; ++c )
      {
         SCIP_CALL( SCIPtransformCons(scip, sourcedata->conss[c], &targetdata->conss[c]) );
      }
   }
   else
   {
      targetdata->conss = NULL;
      targetdata->consssize = 0;
      targetdata->nconss = 0;
   }

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
#define consInitlpConjunction NULL


/** separation method of constraint handler for LP solutions */
#define consSepalpConjunction NULL


/** separation method of constraint handler for arbitrary primal solutions */
#define consSepasolConjunction NULL


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpConjunction)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* add all constraints to the current node */
   SCIP_CALL( addAllConss(scip, conss, nconss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsConjunction)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* add all constraints to the current node */
   SCIP_CALL( addAllConss(scip, conss, nconss, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckConjunction)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* check all constraints of the conjunction */
   SCIP_CALL( checkAllConss(scip, conss, nconss, sol, checkintegrality, checklprows, printreason, result) );

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#define consPropConjunction NULL


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolConjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;
   int i;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* all constraints in a conjunction constraint of the global problem can be added directly to the problem and 
    * removed from the conjunction constraint;
    * an unmodifiable conjunction constraint can be deleted
    */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* add all inactive constraints to the global problem */
      for( i = 0; i < consdata->nconss; ++i )
      {
         /* make sure, the constraint is checked for feasibility */
         SCIP_CALL( SCIPsetConsChecked(scip, consdata->conss[i], TRUE) );

         /* add constraint, if it is not active yet */
         if( !SCIPconsIsActive(consdata->conss[i]) )
         {
            SCIPdebugMessage("adding constraint <%s> from add conjunction <%s>\n",
               SCIPconsGetName(consdata->conss[i]), SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPaddCons(scip, consdata->conss[i]) );
            *result = SCIP_SUCCESS;
         }
         
         /* release constraint from the conjunction constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &consdata->conss[i]) );
      }

      /* now, the conjunction constraint is empty, since all constraints are added directly to the problem */
      consdata->nconss = 0;

      /* delete conjunction constraint, if it is unmodifiable */
      if( !SCIPconsIsModifiable(conss[c]) )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
#define consRespropConjunction NULL


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockConjunction)
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
#define consActiveConjunction NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveConjunction NULL


/** constraint enabling notification method of constraint handler */
#define consEnableConjunction NULL


/** constraint disabling notification method of constraint handler */
#define consDisableConjunction NULL

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintConjunction)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPinfoMessage(scip, file, "conjunction(");
      
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
#define consCopyConjuction NULL

/** constraint parsing method of constraint handler */
#define consParseConjuction NULL



/*
 * constraint specific interface methods
 */

/** creates the handler for conjunction constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrConjunction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create conjunction constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyConjunction,
         consFreeConjunction, consInitConjunction, consExitConjunction, 
         consInitpreConjunction, consExitpreConjunction, consInitsolConjunction, consExitsolConjunction,
         consDeleteConjunction, consTransConjunction, consInitlpConjunction,
         consSepalpConjunction, consSepasolConjunction, consEnfolpConjunction, consEnfopsConjunction, 
         consCheckConjunction, consPropConjunction, consPresolConjunction, consRespropConjunction, consLockConjunction,
         consActiveConjunction, consDeactiveConjunction, 
         consEnableConjunction, consDisableConjunction,
         consPrintConjunction, consCopyConjuction, consParseConjuction,
         conshdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a conjunction constraint */
SCIP_RETCODE SCIPcreateConsConjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nconss,             /**< number of initial constraints in conjunction */
   SCIP_CONS**           conss,              /**< initial constraint in conjunction */
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

   /* find the conjunction constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("conjunction constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, conss, nconss) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, enforce, check, FALSE,
         local, modifiable, dynamic, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** adds constraint to the conjunction of constraints */
SCIP_RETCODE SCIPaddConsElemConjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< conjunction constraint */
   SCIP_CONS*            addcons             /**< additional constraint in conjunction */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(addcons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a conjunction constraint\n");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( consdataAddCons(scip, consdata, addcons) );

   return SCIP_OKAY;
   
}

