/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_and.c
 * @brief  constraint handler for and constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_and.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "and"
#define CONSHDLR_DESC          "and concatenation of constraints"
#define CONSHDLR_SEPAPRIORITY   +000000
#define CONSHDLR_ENFOPRIORITY   +900000
#define CONSHDLR_CHECKPRIORITY  -900000
#define CONSHDLR_SEPAFREQ            -1
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS         TRUE




/*
 * Data structures
 */

/** constraint data for and constraints */
struct ConsData
{
   CONS**           conss;              /**< constraints in concatenation */
   int              consssize;          /**< size of conss array */
   int              nconss;             /**< number of constraints in concatenation */
};

/** constraint handler data */
struct ConshdlrData
{
};




/*
 * Local methods
 */

/** creates and constraint data, captures initial constraints of concatenation */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to constraint data */
   CONS**           conss,              /**< initial constraint in concatenation */
   int              nconss              /**< number of initial constraints in concatenation */
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
            errorMessage("constraints in an and concatenation must not be initial or checked\n");
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

/** frees constraint data and releases all constraints in concatenation */
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

/** adds constraint to concatenation */
static
RETCODE consdataAddCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< constraint data */
   CONS*            cons                /**< constraint to add to the concatenation */
   )
{
   assert(consdata != NULL);

   if( SCIPconsIsInitial(cons) || SCIPconsIsChecked(cons) )
   {
      errorMessage("constraints in an and concatenation must not be initial or checked\n");
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

/** adds all constraints in and constraints to the problem; disables unmodifiable and constraints */
static
RETCODE addAllConss(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           conss,              /**< active and constraints */
   int              nconss,             /**< number of active and constraints */
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
            debugMessage("adding constraint <%s> from add concatenation <%s>\n",
               SCIPconsGetName(consdata->conss[i]), SCIPconsGetName(conss[c]));
            CHECK_OKAY( SCIPaddConsLocal(scip, consdata->conss[i]) );
            *result = SCIP_CONSADDED;
         }
      }

      /* disable and constraint, if it is unmodifiable */
      if( !SCIPconsIsModifiable(conss[c]) )
      {
         CHECK_OKAY( SCIPdisableConsLocal(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** checks all constraints in and constraints for feasibility */
static
RETCODE checkAllConss(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           conss,              /**< active and constraints */
   int              nconss,             /**< number of active and constraints */
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

      /* disable and constraint, if it is unmodifiable */
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
#define consFreeAnd NULL


/** initialization method of constraint handler (called when problem solving starts) */
#define consInitAnd NULL


/** deinitialization method of constraint handler (called when problem solving exits) */
#define consExitAnd NULL


/** solving start notification method of constraint handler (called when presolving was finished) */
#define consSolstartAnd NULL


/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteAnd)
{  /*lint --e{715}*/
   CHECK_OKAY( consdataFree(scip, consdata) );
   
   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransAnd)
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
#define consInitlpAnd NULL


/** separation method of constraint handler */
#define consSepaAnd NULL


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpAnd)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* add all constraints to the current node */
   CHECK_OKAY( addAllConss(scip, conss, nconss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsAnd)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* add all constraints to the current node */
   CHECK_OKAY( addAllConss(scip, conss, nconss, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckAnd)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   /* check all constraints of the concatenation */
   CHECK_OKAY( checkAllConss(scip, conss, nconss, sol, checkintegrality, checklprows, result) );

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#define consPropAnd NULL


/** presolving method of constraint handler */
static
DECL_CONSPRESOL(consPresolAnd)
{  /*lint --e{715}*/
   CONSDATA* consdata;
   int c;
   int i;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* all constraints in an and constraint of the global problem can be added directly to the problem and 
    * removed from the and constraint;
    * an unmodifiable and constraint can be deleted
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
            debugMessage("adding constraint <%s> from add concatenation <%s>\n",
               SCIPconsGetName(consdata->conss[i]), SCIPconsGetName(conss[c]));
            CHECK_OKAY( SCIPaddCons(scip, consdata->conss[i]) );
            *result = SCIP_SUCCESS;
         }
         
         /* release constraint from the and constraint */
         CHECK_OKAY( SCIPreleaseCons(scip, &consdata->conss[i]) );
      }

      /* now, the and constraint is empty, since all constraints are added directly to the problem */
      consdata->nconss = 0;

      /* delete and constraint, if it is unmodifiable */
      if( !SCIPconsIsModifiable(conss[c]) )
      {
         CHECK_OKAY( SCIPdelCons(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}


/** conflict variable resolving method of constraint handler */
#define consRescvarAnd NULL


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockAnd)
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
DECL_CONSUNLOCK(consUnlockAnd)
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
#define consActiveAnd NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveAnd NULL


/** constraint enabling notification method of constraint handler */
#define consEnableAnd NULL


/** constraint disabling notification method of constraint handler */
#define consDisableAnd NULL




/*
 * constraint specific interface methods
 */

/** creates the handler for and constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrAnd(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create and constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_NEEDSCONS,
                  consFreeAnd, consInitAnd, consExitAnd, consSolstartAnd,
                  consDeleteAnd, consTransAnd, consInitlpAnd,
                  consSepaAnd, consEnfolpAnd, consEnfopsAnd, consCheckAnd, 
                  consPropAnd, consPresolAnd, consRescvarAnd,
                  consLockAnd, consUnlockAnd,
                  consActiveAnd, consDeactiveAnd, 
                  consEnableAnd, consDisableAnd,
                  conshdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a and constraint */
RETCODE SCIPcreateConsAnd(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nandconss,          /**< number of initial constraints in concatenation */
   CONS**           andconss,           /**< initial constraint in concatenation */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable          /**< is constraint modifiable (subject to column generation)? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   /* find the and constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("and constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, &consdata, andconss, nandconss) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, enforce, check, FALSE,
                  local, modifiable, FALSE) );

   return SCIP_OKAY;
}

/** adds constraint to the concatenation of an and constraint */
RETCODE SCIPaddConsElemAnd(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< and constraint */
   CONS*            andcons             /**< additional constraint in concatenation */
   )
{
   CONSDATA* consdata;

   assert(cons != NULL);
   assert(andcons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      errorMessage("constraint is not an and constraint\n");
      return SCIP_INVALIDDATA;
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   CHECK_OKAY( consdataAddCons(scip, consdata, andcons) );

   return SCIP_OKAY;
   
}
