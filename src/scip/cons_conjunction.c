/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_conjunction.c,v 1.2 2004/04/27 15:49:57 bzfpfend Exp $"

/**@file   cons_conjunction.c
 * @brief  constraint handler for conjunction constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_conjunction.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "conjunction"
#define CONSHDLR_DESC          "conjunction of constraints"
#define CONSHDLR_SEPAPRIORITY   +000000
#define CONSHDLR_ENFOPRIORITY   +900000
#define CONSHDLR_CHECKPRIORITY  -900000
#define CONSHDLR_SEPAFREQ            -1
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS         TRUE




/*
 * Data structures
 */

/** constraint data for conjunction constraints */
struct ConsData
{
   CONS**           conss;              /**< constraints in conjunction */
   int              consssize;          /**< size of conss array */
   int              nconss;             /**< number of constraints in conjunction */
};




/*
 * Local methods
 */

/** creates conjunction constraint data, captures initial constraints of conjunction */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to constraint data */
   CONS**           conss,              /**< initial constraint in conjunction */
   int              nconss              /**< number of initial constraints in conjunction */
   )
{
   assert(consdata != NULL);
   
   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );
   if( nconss > 0 )
   {
      int c;

      CHECK_OKAY( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->conss, conss, nconss) );
      (*consdata)->consssize = nconss;
      (*consdata)->nconss = nconss;
      for( c = 0; c < nconss; ++c )
      {
         if( SCIPconsIsInitial(conss[c]) || SCIPconsIsChecked(conss[c]) )
         {
            errorMessage("constraints in a conjunction must not be initial or checked\n");
            return SCIP_INVALIDDATA;
         }
         CHECK_OKAY( SCIPcaptureCons(scip, conss[c]) );
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
RETCODE consdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata            /**< pointer to constraint data */
   )
{
   int c;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release constraints */
   for( c = 0; c < (*consdata)->nconss; ++c )
   {
      CHECK_OKAY( SCIPreleaseCons(scip, &(*consdata)->conss[c]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->conss, (*consdata)->consssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** adds constraint to conjunction */
static
RETCODE consdataAddCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< constraint data */
   CONS*            cons                /**< constraint to add to the conjunction */
   )
{
   assert(consdata != NULL);

   if( SCIPconsIsInitial(cons) || SCIPconsIsChecked(cons) )
   {
      errorMessage("constraints in a conjunction must not be initial or checked\n");
      return SCIP_INVALIDDATA;
   }

   /* get memory for additional constraint */
   CHECK_OKAY( SCIPensureBlockMemoryArray(scip, &consdata->conss, &consdata->consssize, consdata->nconss+1) );
   assert(consdata->conss != NULL);
   assert(consdata->nconss < consdata->consssize);

   /* insert constraint in array */
   consdata->conss[consdata->nconss] = cons;
   consdata->nconss++;

   /* capture constraint */
   CHECK_OKAY( SCIPcaptureCons(scip, cons) );

   return SCIP_OKAY;
}

/** adds all constraints in conjunction constraints to the problem; disables unmodifiable conjunction constraints */
static
RETCODE addAllConss(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           conss,              /**< active conjunction constraints */
   int              nconss,             /**< number of active conjunction constraints */
   RESULT*          result              /**< pointer to store the result */
   )
{
   CONSDATA* consdata;
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
            debugMessage("adding constraint <%s> from add conjunction <%s>\n",
               SCIPconsGetName(consdata->conss[i]), SCIPconsGetName(conss[c]));
            CHECK_OKAY( SCIPaddConsLocal(scip, consdata->conss[i]) );
            *result = SCIP_CONSADDED;
         }
      }

      /* disable conjunction constraint, if it is unmodifiable */
      if( !SCIPconsIsModifiable(conss[c]) )
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** checks all constraints in conjunction constraints for feasibility */
static
RETCODE checkAllConss(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           conss,              /**< active conjunction constraints */
   int              nconss,             /**< number of active conjunction constraints */
   SOL*             sol,                /**< solution to check */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result */
   )
{
   CONSDATA* consdata;
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
         CHECK_OKAY( SCIPcheckCons(scip, consdata->conss[i], sol, checkintegrality, checklprows, result) );
      }

      /* disable conjunction constraint, if it is unmodifiable */
      if( !SCIPconsIsModifiable(conss[c]) )
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#define consFreeConjunction NULL


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitConjunction NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitConjunction NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolConjunction NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#define consExitsolConjunction NULL


/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteConjunction)
{  /*lint --e{715}*/
   CHECK_OKAY( consdataFree(scip, consdata) );
   
   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransConjunction)
{  /*lint --e{715}*/
   CONSDATA* sourcedata;
   CONSDATA* targetdata;
   int c;

   /* get constraint data of source constraint */
   sourcedata = SCIPconsGetData(sourcecons);

   /* create constraint data for target constraint */
   CHECK_OKAY( SCIPallocBlockMemory(scip, &targetdata) );
   if( sourcedata->nconss > 0 )
   {
      targetdata->consssize = sourcedata->nconss;
      targetdata->nconss = sourcedata->nconss;
      CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &targetdata->conss, targetdata->consssize) );
      for( c = 0; c < sourcedata->nconss; ++c )
      {
         CHECK_OKAY( SCIPtransformCons(scip, sourcedata->conss[c], &targetdata->conss[c]) );
      }
   }
   else
   {
      targetdata->conss = NULL;
      targetdata->consssize = 0;
      targetdata->nconss = 0;
   }

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                  SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                  SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                  SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
#define consInitlpConjunction NULL


/** separation method of constraint handler */
#define consSepaConjunction NULL


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpConjunction)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* add all constraints to the current node */
   CHECK_OKAY( addAllConss(scip, conss, nconss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsConjunction)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* add all constraints to the current node */
   CHECK_OKAY( addAllConss(scip, conss, nconss, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckConjunction)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* check all constraints of the conjunction */
   CHECK_OKAY( checkAllConss(scip, conss, nconss, sol, checkintegrality, checklprows, result) );

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#define consPropConjunction NULL


/** presolving method of constraint handler */
static
DECL_CONSPRESOL(consPresolConjunction)
{  /*lint --e{715}*/
   CONSDATA* consdata;
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
         CHECK_OKAY( SCIPsetConsChecked(scip, consdata->conss[i]) );

         /* add constraint, if it is not active yet */
         if( !SCIPconsIsActive(consdata->conss[i]) )
         {
            debugMessage("adding constraint <%s> from add conjunction <%s>\n",
               SCIPconsGetName(consdata->conss[i]), SCIPconsGetName(conss[c]));
            CHECK_OKAY( SCIPaddCons(scip, consdata->conss[i]) );
            *result = SCIP_SUCCESS;
         }
         
         /* release constraint from the conjunction constraint */
         CHECK_OKAY( SCIPreleaseCons(scip, &consdata->conss[i]) );
      }

      /* now, the conjunction constraint is empty, since all constraints are added directly to the problem */
      consdata->nconss = 0;

      /* delete conjunction constraint, if it is unmodifiable */
      if( !SCIPconsIsModifiable(conss[c]) )
      {
         CHECK_OKAY( SCIPdelCons(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}


/** conflict variable resolving method of constraint handler */
#define consRescvarConjunction NULL


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockConjunction)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int c;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* lock sub constraints */
   for( c = 0; c < consdata->nconss; ++c )
   {
      CHECK_OKAY( SCIPlockConsVars(scip, consdata->conss[c], nlockspos, nlocksneg) );
   }

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockConjunction)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int c;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* unlock sub constraints */
   for( c = 0; c < consdata->nconss; ++c )
   {
      CHECK_OKAY( SCIPunlockConsVars(scip, consdata->conss[c], nunlockspos, nunlocksneg) );
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




/*
 * constraint specific interface methods
 */

/** creates the handler for conjunction constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrConjunction(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create conjunction constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_NEEDSCONS,
                  consFreeConjunction, consInitConjunction, consExitConjunction, 
                  consInitsolConjunction, consExitsolConjunction,
                  consDeleteConjunction, consTransConjunction, consInitlpConjunction,
                  consSepaConjunction, consEnfolpConjunction, consEnfopsConjunction, consCheckConjunction, 
                  consPropConjunction, consPresolConjunction, consRescvarConjunction,
                  consLockConjunction, consUnlockConjunction,
                  consActiveConjunction, consDeactiveConjunction, 
                  consEnableConjunction, consDisableConjunction,
                  conshdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a conjunction constraint */
RETCODE SCIPcreateConsConjunction(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nconss,             /**< number of initial constraints in conjunction */
   CONS**           conss,              /**< initial constraint in conjunction */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable          /**< is constraint modifiable (subject to column generation)? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   /* find the conjunction constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("conjunction constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, &consdata, conss, nconss) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, enforce, check, FALSE,
                  local, modifiable, FALSE) );

   return SCIP_OKAY;
}

/** adds constraint to the conjunction of constraints */
RETCODE SCIPaddConsElemConjunction(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< conjunction constraint */
   CONS*            addcons             /**< additional constraint in conjunction */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(addcons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not a conjunction constraint\n");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   CHECK_OKAY( consdataAddCons(scip, consdata, addcons) );

   return SCIP_OKAY;
   
}
