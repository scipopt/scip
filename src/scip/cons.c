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

/**@file   cons.c
 * @brief  datastructures and methods for managing constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "cons.h"


/** constraint handler */
struct ConsHdlr
{
   char*            name;               /**< name of constraint handler */
   char*            desc;               /**< description of constraint handler */
   int              sepapriority;       /**< priority of the constraint handler for separation */
   int              enfopriority;       /**< priority of the constraint handler for constraint enforcing */
   int              chckpriority;       /**< priority of the constraint handler for checking infeasibility */
   int              propfreq;           /**< frequency for propagating domains; zero means only preprocessing propagation */
   DECL_CONSFREE((*consfree));          /**< destructor of constraint handler */
   DECL_CONSINIT((*consinit));          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit));          /**< deinitialise constraint handler */
   DECL_CONSDELE((*consdele));          /**< frees specific constraint data */
   DECL_CONSTRAN((*constran));          /**< transforms constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA((*conssepa));          /**< separates cutting planes */
   DECL_CONSENLP((*consenlp));          /**< enforcing constraints for LP solutions */
   DECL_CONSENPS((*consenps));          /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHCK((*conschck));          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop));          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata;       /**< constraint handler data */
   CONS**           sepaconss;          /**< array with active constraints that must be separated during LP processing */
   int              sepaconsssize;      /**< size of sepaconss array */
   int              nsepaconss;         /**< number of active constraints that must be separated during LP processing */
   CONS**           enfoconss;          /**< array with active constraints that must be enforced during node processing */
   int              enfoconsssize;      /**< size of enfoconss array */
   int              nenfoconss;         /**< number of active constraints that must be enforced during node processing */
   CONS**           chckconss;          /**< array with active constraints that must be checked for feasibility */
   int              chckconsssize;      /**< size of chckconss array */
   int              nchckconss;         /**< number of active constraints that must be checked for feasibility */
   CONS**           propconss;          /**< array with active constraints that must be propagated during node processing */
   int              propconsssize;      /**< size of propconss array */
   int              npropconss;         /**< number of active constraints that must be propagated during node processing */
   int              nactiveconss;       /**< total number of active constraints of the handler */
   int              nenabledconss;      /**< total number of enabled constraints of the handler */
   int              lastnsepaconss;     /**< number of already separated constraints after last conshdlrResetSepa() call */
   int              lastnenfoconss;     /**< number of already enforced constraints after last conshdlrResetEnfo() call */
   unsigned int     needscons:1;        /**< should the constraint handler be skipped, if no constraints are available? */
   unsigned int     initialized:1;      /**< is constraint handler initialized? */
};

/** dynamic size attachment for constraint set change data */
struct ConsSetChgDyn
{
   CONSSETCHG**     conssetchg;         /**< pointer to constraint set change data */
   int              addedconsssize;     /**< size of added constraints array */
   int              disabledconsssize;  /**< size of disabled constraints array */
};




/*
 * dynamic memory arrays
 */


/** resizes sepaconss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsureSepaconssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->sepaconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&conshdlr->sepaconss, newsize) );
      conshdlr->sepaconsssize = newsize;
   }
   assert(num <= conshdlr->sepaconsssize);

   return SCIP_OKAY;
}

/** resizes enfoconss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsureEnfoconssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->enfoconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&conshdlr->enfoconss, newsize) );
      conshdlr->enfoconsssize = newsize;
   }
   assert(num <= conshdlr->enfoconsssize);

   return SCIP_OKAY;
}

/** resizes chckconss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsureChckconssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->chckconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&conshdlr->chckconss, newsize) );
      conshdlr->chckconsssize = newsize;
   }
   assert(num <= conshdlr->chckconsssize);

   return SCIP_OKAY;
}

/** resizes propconss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsurePropconssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->propconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&conshdlr->propconss, newsize) );
      conshdlr->propconsssize = newsize;
   }
   assert(num <= conshdlr->propconsssize);

   return SCIP_OKAY;
}


/*
 * Constraint handler methods
 */

/** enables separation, enforcement, and propagation of constraint */
static
RETCODE conshdlrEnableCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(!cons->enabled);
   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->propconsspos == -1);

   debugMessage("enable constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* enable constraint */
   cons->enabled = TRUE;
   conshdlr->nenabledconss++;

   /* add constraint to the separation array */
   if( cons->separate )
   {
      CHECK_OKAY( conshdlrEnsureSepaconssMem(conshdlr, set, conshdlr->nsepaconss+1) );
      cons->sepaconsspos = conshdlr->nsepaconss;
      conshdlr->sepaconss[conshdlr->nsepaconss] = cons;
      conshdlr->nsepaconss++;
   }
      
   /* add constraint to the enforcement array */
   if( cons->enforce )
   {
      CHECK_OKAY( conshdlrEnsureEnfoconssMem(conshdlr, set, conshdlr->nenfoconss+1) );
      cons->enfoconsspos = conshdlr->nenfoconss;
      conshdlr->enfoconss[conshdlr->nenfoconss] = cons;
      conshdlr->nenfoconss++;
   }

   /* add constraint to the propagation array */
   if( cons->propagate )
   {
      CHECK_OKAY( conshdlrEnsurePropconssMem(conshdlr, set, conshdlr->npropconss+1) );
      cons->propconsspos = conshdlr->npropconss;
      conshdlr->propconss[conshdlr->npropconss] = cons;
      conshdlr->npropconss++;
   }

   return SCIP_OKAY;
}

/** disables separation, enforcement, and propagation of constraint */
static
RETCODE conshdlrDisableCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->enabled);
   assert((cons->separate) ^ (cons->sepaconsspos == -1));
   assert((cons->enforce) ^ (cons->enfoconsspos == -1));
   assert((cons->propagate) ^ (cons->propconsspos == -1));

   debugMessage("disable constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* delete constraint from the separation array */
   if( cons->separate )
   {
      delpos = cons->sepaconsspos;
      assert(0 <= delpos && delpos < conshdlr->nsepaconss);
      conshdlr->sepaconss[delpos] = conshdlr->sepaconss[conshdlr->nsepaconss-1];
      conshdlr->sepaconss[delpos]->sepaconsspos = delpos;
      conshdlr->nsepaconss--;
      cons->sepaconsspos = -1;
   }

   /* delete constraint from the enforcement array */
   if( cons->enforce )
   {
      delpos = cons->enfoconsspos;
      assert(0 <= delpos && delpos < conshdlr->nenfoconss);
      conshdlr->enfoconss[delpos] = conshdlr->enfoconss[conshdlr->nenfoconss-1];
      conshdlr->enfoconss[delpos]->enfoconsspos = delpos;
      conshdlr->nenfoconss--;
      cons->enfoconsspos = -1;
   }

   /* delete constraint from the propagation array */
   if( cons->propagate )
   {
      delpos = cons->propconsspos;
      assert(0 <= delpos && delpos < conshdlr->npropconss);
      conshdlr->propconss[delpos] = conshdlr->propconss[conshdlr->npropconss-1];
      conshdlr->propconss[delpos]->propconsspos = delpos;
      conshdlr->npropconss--;
      cons->propconsspos = -1;
   }

   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->propconsspos == -1);

   /* disable constraint */
   cons->enabled = FALSE;
   conshdlr->nenabledconss--;

   return SCIP_OKAY;
}

/** activates and adds constraint to constraint handler's constraint arrays */
static
RETCODE conshdlrAddCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->active);
   assert(!cons->enabled);
   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->chckconsspos == -1);
   assert(cons->propconsspos == -1);

   debugMessage("add constraint <%s> to constraint handler <%s>\n", cons->name, conshdlr->name);

   /* activate constraint */
   cons->active = TRUE;
   conshdlr->nactiveconss++;

   /* add constraint to the check array */
   if( cons->check )
   {
      CHECK_OKAY( conshdlrEnsureChckconssMem(conshdlr, set, conshdlr->nchckconss+1) );
      cons->chckconsspos = conshdlr->nchckconss;
      conshdlr->chckconss[conshdlr->nchckconss] = cons;
      conshdlr->nchckconss++;
   }

   /* enable separation, enforcement, and propagation of constraint */
   CHECK_OKAY( conshdlrEnableCons(conshdlr, set, cons) );

   return SCIP_OKAY;
}

/** deactivates and removes constraint from constraint handler's conss array */
static
RETCODE conshdlrDelCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert((cons->check) ^ (cons->chckconsspos == -1));

   debugMessage("delete constraint <%s> from constraint handler <%s>\n", cons->name, conshdlr->name);

   /* disable constraint */
   if( cons->enabled )
   {
      CHECK_OKAY( conshdlrDisableCons(conshdlr, cons) );
   }
   assert(!cons->enabled);

   /* delete constraint from the check array */
   if( cons->check )
   {
      delpos = cons->chckconsspos;
      assert(0 <= delpos && delpos < conshdlr->nchckconss);
      conshdlr->chckconss[delpos] = conshdlr->chckconss[conshdlr->nchckconss-1];
      conshdlr->chckconss[delpos]->chckconsspos = delpos;
      conshdlr->nchckconss--;
      cons->chckconsspos = -1;
   }

   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->chckconsspos == -1);
   assert(cons->propconsspos == -1);

   /* deactivate constraint */
   cons->active = FALSE;
   conshdlr->nactiveconss--;

   return SCIP_OKAY;
}

DECL_SORTPTRCOMP(SCIPconshdlrCompSepa)  /**< compares two constraint handlers w. r. to their separation priority */
{
   return ((CONSHDLR*)elem2)->sepapriority - ((CONSHDLR*)elem1)->sepapriority;
}

DECL_SORTPTRCOMP(SCIPconshdlrCompEnfo)  /**< compares two constraint handlers w. r. to their enforcing priority */
{
   return ((CONSHDLR*)elem2)->enfopriority - ((CONSHDLR*)elem1)->enfopriority;
}

DECL_SORTPTRCOMP(SCIPconshdlrCompChck)  /**< compares two constraint handlers w. r. to their feasibility check priority */
{
   return ((CONSHDLR*)elem2)->chckpriority - ((CONSHDLR*)elem1)->chckpriority;
}

/** creates a constraint handler */
RETCODE SCIPconshdlrCreate(
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              chckpriority,       /**< priority of the constraint handler for checking infeasibility */
   int              propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE((*consfree)),          /**< destructor of constraint handler */
   DECL_CONSINIT((*consinit)),          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit)),          /**< deinitialise constraint handler */
   DECL_CONSDELE((*consdele)),          /**< free specific constraint data */
   DECL_CONSTRAN((*constran)),          /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA((*conssepa)),          /**< separate cutting planes */
   DECL_CONSENLP((*consenlp)),          /**< enforcing constraints for LP solutions */
   DECL_CONSENPS((*consenps)),          /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHCK((*conschck)),          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop)),          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   assert(conshdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert((propfreq >= 0) ^ (consprop == NULL));

   ALLOC_OKAY( allocMemory(conshdlr) );
   ALLOC_OKAY( duplicateMemoryArray(&(*conshdlr)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*conshdlr)->desc, desc, strlen(desc)+1) );
   (*conshdlr)->sepapriority = sepapriority;
   (*conshdlr)->enfopriority = enfopriority;
   (*conshdlr)->chckpriority = chckpriority;
   (*conshdlr)->propfreq = propfreq;
   (*conshdlr)->consfree = consfree;
   (*conshdlr)->consinit = consinit;
   (*conshdlr)->consexit = consexit;
   (*conshdlr)->consdele = consdele;
   (*conshdlr)->constran = constran;
   (*conshdlr)->conssepa = conssepa;
   (*conshdlr)->consenlp = consenlp;
   (*conshdlr)->consenps = consenps;
   (*conshdlr)->conschck = conschck;
   (*conshdlr)->consprop = consprop;
   (*conshdlr)->conshdlrdata = conshdlrdata;
   (*conshdlr)->sepaconss = NULL;
   (*conshdlr)->sepaconsssize = 0;
   (*conshdlr)->nsepaconss = 0;
   (*conshdlr)->enfoconss = NULL;
   (*conshdlr)->enfoconsssize = 0;
   (*conshdlr)->nenfoconss = 0;
   (*conshdlr)->chckconss = NULL;
   (*conshdlr)->chckconsssize = 0;
   (*conshdlr)->nchckconss = 0;
   (*conshdlr)->propconss = NULL;
   (*conshdlr)->propconsssize = 0;
   (*conshdlr)->npropconss = 0;
   (*conshdlr)->nactiveconss = 0;
   (*conshdlr)->nenabledconss = 0;
   (*conshdlr)->lastnsepaconss = 0;
   (*conshdlr)->lastnenfoconss = 0;
   (*conshdlr)->needscons = needscons;
   (*conshdlr)->initialized = FALSE;

   return SCIP_OKAY;
}

/** calls destructor and frees memory of constraint handler */
RETCODE SCIPconshdlrFree(
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conshdlr != NULL);
   assert(*conshdlr != NULL);
   assert(!(*conshdlr)->initialized);
   assert(scip != NULL);

   /* call destructor of constraint handler */
   if( (*conshdlr)->consfree != NULL )
   {
      CHECK_OKAY( (*conshdlr)->consfree(scip, *conshdlr) );
   }

   freeMemoryArray(&(*conshdlr)->name);
   freeMemoryArray(&(*conshdlr)->desc);
   freeMemoryArrayNull(&(*conshdlr)->sepaconss);
   freeMemoryArrayNull(&(*conshdlr)->enfoconss);
   freeMemoryArrayNull(&(*conshdlr)->chckconss);
   freeMemoryArrayNull(&(*conshdlr)->propconss);
   freeMemory(conshdlr);

   return SCIP_OKAY;
}

/** initializes constraint handler */
RETCODE SCIPconshdlrInit(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conshdlr != NULL);
   assert(scip != NULL);

   if( conshdlr->initialized )
   {
      char s[255];
      sprintf(s, "Constraint handler <%s> already initialized", conshdlr->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( conshdlr->consinit != NULL )
   {
      CHECK_OKAY( conshdlr->consinit(scip, conshdlr) );
   }
   conshdlr->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of constraint handler */
RETCODE SCIPconshdlrExit(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conshdlr != NULL);
   assert(scip != NULL);

   if( !conshdlr->initialized )
   {
      char s[255];
      sprintf(s, "Constraint handler <%s> not initialized", conshdlr->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( conshdlr->consexit != NULL )
   {
      CHECK_OKAY( conshdlr->consexit(scip, conshdlr) );
   }
   conshdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** calls separator method of constraint handler to separate all constraints added after last conshdlrResetSepa() call */
RETCODE SCIPconshdlrSeparate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(0 <= conshdlr->lastnsepaconss && conshdlr->lastnsepaconss <= conshdlr->nsepaconss);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->conssepa != NULL )
   {
      int nconss;

      nconss = conshdlr->nsepaconss - conshdlr->lastnsepaconss;

      if( !conshdlr->needscons || nconss > 0 )
      {
         CONS** conss;
         
         debugMessage("separating constraints %d to %d of %d constraints of handler <%s>\n",
            conshdlr->lastnsepaconss, conshdlr->lastnsepaconss + nconss - 1, conshdlr->nsepaconss, conshdlr->name);

         conss = &(conshdlr->sepaconss[conshdlr->lastnsepaconss]);
         conshdlr->lastnsepaconss = conshdlr->nsepaconss;

         CHECK_OKAY( conshdlr->conssepa(set->scip, conshdlr, conss, nconss, result) );
         debugMessage(" -> separating returned result <%d>\n", *result);

         if( *result != SCIP_SEPARATED
            && *result != SCIP_CONSADDED
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN )
         {
            char s[255];
            sprintf(s, "separation method of constraint handler <%s> returned invalid result <%d>", 
               conshdlr->name, *result);
            errorMessage(s);
            return SCIP_INVALIDRESULT;
         }
      }
   }

   return SCIP_OKAY;
}

/** calls enforcing method of constraint handler for LP solution for all constraints added after last
 *  conshdlrResetEnfo() call
 */
RETCODE SCIPconshdlrEnforceLPSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(0 <= conshdlr->lastnenfoconss && conshdlr->lastnenfoconss <= conshdlr->nenfoconss);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->consenlp != NULL )
   {
      int nconss;

      nconss = conshdlr->nenfoconss - conshdlr->lastnenfoconss;

      if( !conshdlr->needscons || nconss > 0 )
      {
         CONS** conss;
         
         debugMessage("enforcing constraints %d to %d of %d constraints of handler <%s>\n",
            conshdlr->lastnenfoconss, conshdlr->lastnenfoconss + nconss - 1, conshdlr->nenfoconss, conshdlr->name);

         conss = &(conshdlr->enfoconss[conshdlr->lastnenfoconss]);
         conshdlr->lastnenfoconss = conshdlr->nenfoconss;

         CHECK_OKAY( conshdlr->consenlp(set->scip, conshdlr, conss, nconss, result) );
         debugMessage(" -> enforcing returned result <%d>\n", *result);

         if( *result != SCIP_CUTOFF
            && *result != SCIP_BRANCHED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_CONSADDED
            && *result != SCIP_INFEASIBLE
            && *result != SCIP_FEASIBLE )
         {
            char s[255];
            sprintf(s, "enforcing method of constraint handler <%s> for LP solutions returned invalid result <%d>", 
               conshdlr->name, *result);
            errorMessage(s);
            return SCIP_INVALIDRESULT;
         }
      }
   }

   return SCIP_OKAY;
}

/** calls enforcing method of constraint handler for pseudo solution for all constraints added after last
 *  conshdlrResetEnfo() call
 */
RETCODE SCIPconshdlrEnforcePseudoSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(0 <= conshdlr->lastnenfoconss && conshdlr->lastnenfoconss <= conshdlr->nenfoconss);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->consenps != NULL )
   {
      int nconss;

      nconss = conshdlr->nenfoconss - conshdlr->lastnenfoconss;

      if( !conshdlr->needscons || nconss > 0 )
      {
         CONS** conss;
         
         debugMessage("enforcing constraints %d to %d of %d constraints of handler <%s>\n",
            conshdlr->lastnenfoconss, conshdlr->lastnenfoconss + nconss - 1, conshdlr->nenfoconss, conshdlr->name);

         conss = &(conshdlr->enfoconss[conshdlr->lastnenfoconss]);
         conshdlr->lastnenfoconss = conshdlr->nenfoconss;

         CHECK_OKAY( conshdlr->consenps(set->scip, conshdlr, conss, nconss, result) );
         debugMessage(" -> enforcing returned result <%d>\n", *result);

         if( *result != SCIP_CUTOFF
            && *result != SCIP_BRANCHED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_CONSADDED
            && *result != SCIP_INFEASIBLE
            && *result != SCIP_FEASIBLE )
         {
            char s[255];
            sprintf(s, "enforcing method of constraint handler <%s> for pseudo solutions returned invalid result <%d>", 
               conshdlr->name, *result);
            errorMessage(s);
            return SCIP_INVALIDRESULT;
         }
      }
   }

   return SCIP_OKAY;
}

/** calls feasibility check method of constraint handler */
RETCODE SCIPconshdlrCheck(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   SOL*             sol,                /**< primal CIP solution */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->conschck != NULL && (!conshdlr->needscons || conshdlr->nchckconss > 0) )
   {
      debugMessage("checking %d constraints of handler <%s>\n", conshdlr->nchckconss, conshdlr->name);
      CHECK_OKAY( conshdlr->conschck(set->scip, conshdlr, conshdlr->chckconss, conshdlr->nchckconss, 
                     sol, chckintegrality, chcklprows, result) );
      debugMessage(" -> checking returned result <%d>\n", *result);
      if( *result != SCIP_INFEASIBLE
         && *result != SCIP_FEASIBLE )
      {
         char s[255];
         sprintf(s, "feasibility check of constraint handler <%s> returned invalid result <%d>", 
            conshdlr->name, *result);
         errorMessage(s);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** calls propagation method of constraint handler */
RETCODE SCIPconshdlrPropagate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   int              actdepth,           /**< depth of active node; -1 if preprocessing domain propagation */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->consprop != NULL
      && (!conshdlr->needscons || conshdlr->npropconss > 0)
      && (actdepth == -1 || (conshdlr->propfreq > 0 && actdepth % conshdlr->propfreq == 0)) )
   {
      debugMessage("propagating %d constraints of handler <%s>\n", conshdlr->npropconss, conshdlr->name);
      CHECK_OKAY( conshdlr->consprop(set->scip, conshdlr, conshdlr->propconss, conshdlr->npropconss, result) );
      debugMessage(" -> propagation returned result <%d>\n", *result);
      if( *result != SCIP_CUTOFF
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         char s[255];
         sprintf(s, "propagation method of constraint handler <%s> returned invalid result <%d>", 
            conshdlr->name, *result);
         errorMessage(s);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** resets separation to start with first constraint in the next call */
void SCIPconshdlrResetSepa(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   conshdlr->lastnsepaconss = 0;
}

/** resets enforcement to start with first constraint in the next call */
void SCIPconshdlrResetEnfo(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   conshdlr->lastnenfoconss = 0;
}

/** gets name of constraint handler */
const char* SCIPconshdlrGetName(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->name;
}

/** gets user data of constraint handler */
CONSHDLRDATA* SCIPconshdlrGetData(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->conshdlrdata;
}

/** sets user data of constraint handler; user has to free old data in advance! */
void SCIPconshdlrSetData(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   CONSHDLRDATA*    conshdlrdata        /**< new constraint handler user data */
   )
{
   assert(conshdlr != NULL);

   conshdlr->conshdlrdata = conshdlrdata;
}

/** gets number of active constraints of constraint handler */
int SCIPconshdlrGetNActiveConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nactiveconss;
}

/** gets number of enabled constraints of constraint handler */
int SCIPconshdlrGetNEnabledConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nenabledconss;
}

/** gets checking priority of constraint handler */
int SCIPconshdlrGetChckPriority(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->chckpriority;
}

/** gets propagation frequency of constraint handler */
int SCIPconshdlrGetPropFreq(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->propfreq;
}

/** is constraint handler initialized? */
Bool SCIPconshdlrIsInitialized(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->initialized;
}




/*
 * Constraint methods
 */

/** creates and captures a constraint
 *  Warning! If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
RETCODE SCIPconsCreate(
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             original            /**< is constraint belonging to the original problem? */
   )
{
   assert(cons != NULL);
   assert(memhdr != NULL);
   assert(conshdlr != NULL);

   /* create constraint data */
   ALLOC_OKAY( allocBlockMemory(memhdr, cons) );
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, &(*cons)->name, name, strlen(name)+1) );
   (*cons)->conshdlr = conshdlr;
   (*cons)->consdata = consdata;
   (*cons)->node = NULL;
   (*cons)->nuses = 0;
   (*cons)->sepaconsspos = -1;
   (*cons)->enfoconsspos = -1;
   (*cons)->chckconsspos = -1;
   (*cons)->propconsspos = -1;
   (*cons)->arraypos = -1;
   (*cons)->separate = separate;
   (*cons)->enforce = enforce;
   (*cons)->check = check;
   (*cons)->propagate = propagate;
   (*cons)->original = original;
   (*cons)->active = FALSE;
   (*cons)->enabled = FALSE;

   /* capture constraint */
   SCIPconsCapture(*cons);

   return SCIP_OKAY;
}

/** frees a constraint */
RETCODE SCIPconsFree(
   CONS**           cons,               /**< constraint to free */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);
   assert((*cons)->nuses == 0);
   assert((*cons)->conshdlr != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);

   /* free constraint data */
   if( (*cons)->conshdlr->consdele != NULL )
   {
      CHECK_OKAY( (*cons)->conshdlr->consdele(set->scip, (*cons)->conshdlr, &(*cons)->consdata) );
   }
   freeBlockMemoryArray(memhdr, &(*cons)->name, strlen((*cons)->name)+1);
   freeBlockMemory(memhdr, cons);

   return SCIP_OKAY;
}

/** increases usage counter of constraint */
void SCIPconsCapture(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(cons->nuses >= 0);

   debugMessage("capture constraint <%s> with nuses=%d\n", cons->name, cons->nuses);
   cons->nuses++;
}

/** decreases usage counter of constraint, and frees memory if necessary */
RETCODE SCIPconsRelease(
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(memhdr != NULL);
   assert(cons != NULL);
   assert(*cons != NULL);
   assert((*cons)->nuses >= 1);

   debugMessage("release constraint <%s> with nuses=%d\n", (*cons)->name, (*cons)->nuses);
   (*cons)->nuses--;
   if( (*cons)->nuses == 0 )
   {
      CHECK_OKAY( SCIPconsFree(cons, memhdr, set) );
   }
   *cons  = NULL;

   return SCIP_OKAY;
}

/** activates constraint */
RETCODE SCIPconsActivate(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(!cons->active);
   assert(set != NULL);
   
   CHECK_OKAY( conshdlrAddCons(cons->conshdlr, set, cons) );
   assert(cons->active);

   return SCIP_OKAY;
}

/** deactivates constraint */
RETCODE SCIPconsDeactivate(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(cons->active);
   
   CHECK_OKAY( conshdlrDelCons(cons->conshdlr, cons) );
   assert(!cons->active);

   return SCIP_OKAY;
}

/** enables constraint's separation, enforcing, and propagation capabilities */
RETCODE SCIPconsEnable(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->active);
   assert(!cons->enabled);
   assert(set != NULL);
   
   CHECK_OKAY( conshdlrEnableCons(cons->conshdlr, set, cons) );
   assert(cons->enabled);

   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities */
RETCODE SCIPconsDisable(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(cons->active);
   assert(cons->enabled);

   CHECK_OKAY( conshdlrDisableCons(cons->conshdlr, cons) );
   assert(!cons->enabled);

   return SCIP_OKAY;
}

/** copies original constraint into transformed constraint, that is captured */
RETCODE SCIPconsTransform(
   CONS**           transcons,          /**< pointer to store the transformed constraint */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   CONS*            origcons            /**< original constraint */
   )
{
   CONSDATA* consdata;

   assert(transcons != NULL);
   assert(memhdr != NULL);
   assert(origcons != NULL);

   /* transform constraint data */
   consdata = NULL;
   if( origcons->conshdlr->constran != NULL )
   {
      CHECK_OKAY( origcons->conshdlr->constran(set->scip, origcons->conshdlr, origcons->consdata, &consdata) );
   }

   /* create new constraint with transformed data */
   CHECK_OKAY( SCIPconsCreate(transcons, memhdr, origcons->name, origcons->conshdlr, consdata,
                  origcons->separate, origcons->enforce, origcons->check, origcons->propagate, FALSE) );

   return SCIP_OKAY;
}

/** returns the name of the constraint */
const char* SCIPconsGetName(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->name;
}

/** returns the constraint handler of the constraint */
CONSHDLR* SCIPconsGetConsHdlr(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->conshdlr;
}

/** returns the constraint data field of the constraint */
CONSDATA* SCIPconsGetConsData(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->consdata;
}

/** returns TRUE iff constraint is belonging to original problem */
Bool SCIPconsIsOriginal(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->original;
}




/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given constraint */
DECL_HASHGETKEY(SCIPhashGetKeyCons)
{
   CONS* cons = (CONS*)elem;

   assert(cons != NULL);
   return cons->name;
}




/*
 * Constraint set change methods
 */

/** creates empty fixed size constraint set change data */
static
RETCODE conssetchgCreate(
   CONSSETCHG**     conssetchg,         /**< pointer to fixed size constraint set change data */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(conssetchg != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, conssetchg) );
   (*conssetchg)->addedconss = NULL;
   (*conssetchg)->disabledconss = NULL;
   (*conssetchg)->naddedconss = 0;
   (*conssetchg)->ndisabledconss = 0;

   return SCIP_OKAY;
}

/** releases all constraints of the constraint set change data */
static
RETCODE conssetchgRelease(
   CONSSETCHG*      conssetchg,         /**< constraint set change data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   )
{
   CONS* cons;
   int i;
   
   assert(conssetchg != NULL);
   
   /* release constraints */
   for( i = 0; i < conssetchg->naddedconss; ++i )
   {
      cons = conssetchg->addedconss[i];
      if( cons != NULL )
      {
         assert(cons->arraypos == i);
         cons->node = NULL;
         cons->arraypos = -1;
         CHECK_OKAY( SCIPconsRelease(&conssetchg->addedconss[i], memhdr, set) );
      }
   }
   for( i = 0; i < conssetchg->ndisabledconss; ++i )
   {
      if( conssetchg->disabledconss[i] != NULL )
      {
         CHECK_OKAY( SCIPconsRelease(&conssetchg->disabledconss[i], memhdr, set) );
      }
   }

   return SCIP_OKAY;
}

/** frees fixed size constraint set change data and releases all included constraints */
RETCODE SCIPconssetchgFree(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(conssetchg != NULL);
   assert(memhdr != NULL);

   if( *conssetchg != NULL )
   {
      /* release constraints */
      CHECK_OKAY( conssetchgRelease(*conssetchg, memhdr, set) );

      /* free memory */
      freeBlockMemoryArrayNull(memhdr, &(*conssetchg)->addedconss, (*conssetchg)->naddedconss);
      freeBlockMemoryArrayNull(memhdr, &(*conssetchg)->disabledconss, (*conssetchg)->ndisabledconss);
      freeBlockMemory(memhdr, conssetchg);
   }

   return SCIP_OKAY;
}

/** deletes and releases deactivated constraint from the addedconss array of the constraint set change data */
RETCODE SCIPconssetchgDelAddedCons(
   CONSSETCHG*      conssetchg,         /**< constraint set change to delete constraint from */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to delete from addedconss array */
   )
{
   int arraypos;

   assert(conssetchg != NULL);
   assert(cons != NULL);
   assert(!cons->active);
   assert(!cons->enabled);

   arraypos = cons->arraypos;
   assert(0 <= arraypos && arraypos < conssetchg->naddedconss);
   assert(conssetchg->addedconss[arraypos] == cons);

   cons->node = NULL;
   cons->arraypos = -1;
   CHECK_OKAY( SCIPconsRelease(&conssetchg->addedconss[arraypos], memhdr, set) );

   assert(conssetchg->addedconss[arraypos] == NULL);

   return SCIP_OKAY;
}

/** deletes and releases deactivated constraint from the disabledconss array of the constraint set change data */
static
RETCODE conssetchgDelDisabledCons(
   CONSSETCHG*      conssetchg,         /**< constraint set change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   int              arraypos            /**< position of constraint in disabledconss array */
   )
{
   assert(conssetchg != NULL);
   assert(0 <= arraypos && arraypos < conssetchg->ndisabledconss);
   assert(conssetchg->disabledconss[arraypos] != NULL);
   assert(!conssetchg->disabledconss[arraypos]->active);
   assert(!conssetchg->disabledconss[arraypos]->enabled);

   CHECK_OKAY( SCIPconsRelease(&conssetchg->disabledconss[arraypos], memhdr, set) );

   assert(conssetchg->disabledconss[arraypos] == NULL);

   return SCIP_OKAY;
}

/** applies constraint set change */
RETCODE SCIPconssetchgApply(
   CONSSETCHG*      conssetchg,         /**< constraint set change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   )
{
   CONS* cons;
   int i;

   debugMessage("applying constraint set changes at %p\n", conssetchg);

   if( conssetchg == NULL )
      return SCIP_OKAY;

   debugMessage(" -> %d constraint additions, %d constraint disablings\n", 
      conssetchg->naddedconss, conssetchg->ndisabledconss);

   /* apply constraint additions */
   for( i = 0; i < conssetchg->naddedconss; ++i )
   {
      cons = conssetchg->addedconss[i];
      if( cons != NULL )
      {
         CHECK_OKAY( SCIPconsActivate(cons, set) );
      }
   }

   /* apply constraint disablings */
   for( i = 0; i < conssetchg->ndisabledconss; ++i )
   {
      cons = conssetchg->disabledconss[i];

      if( cons != NULL )
      {
         int arraypos;

         arraypos = cons->arraypos;
         assert(!cons->active || arraypos >= 0);

         /* if the constraint is inactive, we can permanently remove it from the disabledconss array
          * if the constraint was added at this node and is no check-constraint, we can deactivate and remove it from
          * both arrays
          */
         if( !cons->active )
         {
            debugMessage("constraint <%s> of handler <%s> was deactivated -> remove it from disabledconss array\n",
               cons->name, cons->conshdlr->name);
            
            /* release and remove constraint from the disabledconss array */
            CHECK_OKAY( SCIPconsRelease(&conssetchg->disabledconss[i], memhdr, set) );
         }
         else if( !cons->check && arraypos < conssetchg->naddedconss && cons == conssetchg->addedconss[arraypos] )
         {
            debugMessage("constraint <%s> of handler <%s> was deactivated at same node -> remove it from both arrays\n",
               cons->name, cons->conshdlr->name);
            
            /* deactivate constraint */
            CHECK_OKAY( SCIPconsDeactivate(conssetchg->addedconss[i]) );

            /* release and remove constraint from the addedconss array */
            CHECK_OKAY( SCIPconssetchgDelAddedCons(conssetchg, memhdr, set, cons) );

            /* release and remove constraint from the disabledconss array */
            CHECK_OKAY( conssetchgDelDisabledCons(conssetchg, memhdr, set, i) );
         }
         else
         {
            CHECK_OKAY( SCIPconsDisable(conssetchg->disabledconss[i]) );
         }
      }
   }

   return SCIP_OKAY;
}

/** undoes constraint set change */
RETCODE SCIPconssetchgUndo(
   CONSSETCHG*      conssetchg,         /**< constraint set change to undo */
   const SET*       set                 /**< global SCIP settings */
   )
{
   int i;

   debugMessage("undoing constraint set changes at %p\n", conssetchg);

   if( conssetchg == NULL )
      return SCIP_OKAY;

   debugMessage(" -> %d constraint additions, %d constraint disablings\n", 
      conssetchg->naddedconss, conssetchg->ndisabledconss);

   /* undo constraint disablings */
   for( i = conssetchg->ndisabledconss-1; i >= 0; --i )
   {
      CHECK_OKAY( SCIPconsEnable(conssetchg->disabledconss[i], set) );
   }

   /* undo constraint additions */
   for( i = conssetchg->naddedconss-1; i >= 0; --i )
   {
      CHECK_OKAY( SCIPconsDeactivate(conssetchg->addedconss[i]) );
   }

   return SCIP_OKAY;
}




/*
 * dynamic size attachment methods for constraint set changes
 */

/** ensures, that addedconss array can store at least num entries */
static
RETCODE ensureAddedconssSize(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamically sized constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   CONSSETCHG** conssetchg;

   assert(conssetchgdyn != NULL);
   assert(conssetchgdyn->conssetchg != NULL);

   conssetchg = conssetchgdyn->conssetchg;
   assert(*conssetchg != NULL || conssetchgdyn->addedconsssize == 0);
   assert(*conssetchg == NULL || (*conssetchg)->naddedconss <= conssetchgdyn->addedconsssize);

   if( num > conssetchgdyn->addedconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      if( *conssetchg == NULL )
      {
         CHECK_OKAY( conssetchgCreate(conssetchg, memhdr) );
      }
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*conssetchg)->addedconss, conssetchgdyn->addedconsssize, newsize) );
      conssetchgdyn->addedconsssize = newsize;
   }
   assert(num <= conssetchgdyn->addedconsssize);

   return SCIP_OKAY;
}

/** ensures, that disabledconss array can store at least num entries */
static
RETCODE ensureDisabledconssSize(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamically sized constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   CONSSETCHG** conssetchg;

   assert(conssetchgdyn != NULL);
   assert(conssetchgdyn->conssetchg != NULL);

   conssetchg = conssetchgdyn->conssetchg;
   assert(*conssetchg != NULL || conssetchgdyn->disabledconsssize == 0);
   assert(*conssetchg == NULL || (*conssetchg)->ndisabledconss <= conssetchgdyn->disabledconsssize);

   if( num > conssetchgdyn->disabledconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      if( *conssetchg == NULL )
      {
         CHECK_OKAY( conssetchgCreate(conssetchg, memhdr) );
      }
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*conssetchg)->disabledconss, conssetchgdyn->disabledconsssize, 
                     newsize) );
      conssetchgdyn->disabledconsssize = newsize;
   }
   assert(num <= conssetchgdyn->disabledconsssize);

   return SCIP_OKAY;
}

/** creates a dynamic size attachment for a constraint set change data structure */
RETCODE SCIPconssetchgdynCreate(
   CONSSETCHGDYN**  conssetchgdyn,      /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(conssetchgdyn != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, conssetchgdyn) );

   (*conssetchgdyn)->conssetchg = NULL;
   (*conssetchgdyn)->addedconsssize = 0;
   (*conssetchgdyn)->disabledconsssize = 0;

   return SCIP_OKAY;
}

/** frees a dynamic size attachment for a constraint set change data structure */
void SCIPconssetchgdynFree(
   CONSSETCHGDYN**  conssetchgdyn,      /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(conssetchgdyn != NULL);
   assert(*conssetchgdyn != NULL);
   assert((*conssetchgdyn)->conssetchg == NULL);
   assert((*conssetchgdyn)->addedconsssize == 0);
   assert((*conssetchgdyn)->disabledconsssize == 0);

   freeBlockMemory(memhdr, conssetchgdyn);
}

/** attaches dynamic size information to constraint set change data */
void SCIPconssetchgdynAttach(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamic size information */
   CONSSETCHG**     conssetchg          /**< pointer to static constraint set change */
   )
{
   assert(conssetchgdyn != NULL);
   assert(conssetchgdyn->conssetchg == NULL);
   assert(conssetchgdyn->addedconsssize == 0);
   assert(conssetchgdyn->disabledconsssize == 0);
   assert(conssetchg != NULL);

   debugMessage("attaching dynamic size information at %p to constraint set change data at %p\n",
      conssetchgdyn, conssetchg);

   conssetchgdyn->conssetchg = conssetchg;
   if( *conssetchg != NULL )
   {
      conssetchgdyn->addedconsssize = (*conssetchg)->naddedconss;
      conssetchgdyn->disabledconsssize = (*conssetchg)->ndisabledconss;
   }
   else
   {
      conssetchgdyn->addedconsssize = 0;
      conssetchgdyn->disabledconsssize = 0;
   }
}

/** detaches dynamic size information and shrinks constraint set change data */
RETCODE SCIPconssetchgdynDetach(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamic size information */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(conssetchgdyn != NULL);
   assert(conssetchgdyn->conssetchg != NULL);
   assert(conssetchgdyn->addedconsssize == 0 || *conssetchgdyn->conssetchg != NULL);
   assert(conssetchgdyn->addedconsssize == 0 || (*conssetchgdyn->conssetchg)->naddedconss <= conssetchgdyn->addedconsssize);
   assert(conssetchgdyn->disabledconsssize == 0 || *conssetchgdyn->conssetchg != NULL);
   assert(conssetchgdyn->disabledconsssize == 0
      || (*conssetchgdyn->conssetchg)->ndisabledconss <= conssetchgdyn->disabledconsssize);

   debugMessage("detaching dynamic size information at %p from constraint set change data at %p\n",
      conssetchgdyn, conssetchgdyn->conssetchg);

   /* shrink static constraint set change data to the size of the used elements */
   if( *conssetchgdyn->conssetchg != NULL )
   {
      if( (*conssetchgdyn->conssetchg)->naddedconss == 0 )
      {
         freeBlockMemoryArrayNull(memhdr, &(*conssetchgdyn->conssetchg)->addedconss, conssetchgdyn->addedconsssize);
      }
      else
      {
         ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*conssetchgdyn->conssetchg)->addedconss,
                        conssetchgdyn->addedconsssize, (*conssetchgdyn->conssetchg)->naddedconss) );
      }
      if( (*conssetchgdyn->conssetchg)->ndisabledconss == 0 )
      {
         freeBlockMemoryArrayNull(memhdr, &(*conssetchgdyn->conssetchg)->disabledconss, conssetchgdyn->disabledconsssize);
      }
      else
      {
         ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &(*conssetchgdyn->conssetchg)->disabledconss,
                        conssetchgdyn->disabledconsssize, (*conssetchgdyn->conssetchg)->ndisabledconss) );
      }
      if( (*conssetchgdyn->conssetchg)->naddedconss == 0 && (*conssetchgdyn->conssetchg)->ndisabledconss == 0 )
      {
         CHECK_OKAY( SCIPconssetchgFree(conssetchgdyn->conssetchg, memhdr, set) );
      }
   }
   else
   {
      assert(conssetchgdyn->addedconsssize == 0);
      assert(conssetchgdyn->disabledconsssize == 0);
   }

   /* detach constraint set change data */
   conssetchgdyn->conssetchg = NULL;
   conssetchgdyn->disabledconsssize = 0;
   conssetchgdyn->addedconsssize = 0;

   return SCIP_OKAY;
}

/** frees attached constraint set change data and detaches dynamic size attachment */
RETCODE SCIPconssetchgdynDiscard(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamically sized constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(conssetchgdyn != NULL);
   assert(conssetchgdyn->conssetchg != NULL);
   assert(conssetchgdyn->addedconsssize == 0 || *conssetchgdyn->conssetchg != NULL);
   assert(conssetchgdyn->addedconsssize == 0 || (*conssetchgdyn->conssetchg)->naddedconss <= conssetchgdyn->addedconsssize);
   assert(conssetchgdyn->disabledconsssize == 0 || *conssetchgdyn->conssetchg != NULL);
   assert(conssetchgdyn->disabledconsssize == 0
      || (*conssetchgdyn->conssetchg)->ndisabledconss <= conssetchgdyn->disabledconsssize);

   debugMessage("discarding dynamic size information at %p from constraint set change data at %p\n", 
      conssetchgdyn, conssetchgdyn->conssetchg);

   /* free static constraint set change data */
   if( *conssetchgdyn->conssetchg != NULL )
   {  
      /* release constraints of the constraint set change data */
      CHECK_OKAY( conssetchgRelease(*conssetchgdyn->conssetchg, memhdr, set) );

      /* clear the constraint set change data */
      freeBlockMemoryArrayNull(memhdr, &(*conssetchgdyn->conssetchg)->addedconss, conssetchgdyn->addedconsssize);
      (*conssetchgdyn->conssetchg)->naddedconss = 0;
      freeBlockMemoryArrayNull(memhdr, &(*conssetchgdyn->conssetchg)->disabledconss, conssetchgdyn->disabledconsssize);
      (*conssetchgdyn->conssetchg)->ndisabledconss = 0;

      /* free the constraint set change data */
      CHECK_OKAY( SCIPconssetchgFree(conssetchgdyn->conssetchg, memhdr, set) );
      conssetchgdyn->addedconsssize = 0;
      conssetchgdyn->disabledconsssize = 0;
   }

   /* detach constraint set change data */
   conssetchgdyn->conssetchg = NULL;
   assert(conssetchgdyn->addedconsssize == 0);
   assert(conssetchgdyn->disabledconsssize == 0);

   return SCIP_OKAY;
}

/** adds constraint addition to constraint set changes, and captures constraint */
RETCODE SCIPconssetchgdynAddAddedCons(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamically sized constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   NODE*            node,               /**< node that the constraint set change belongs to */
   CONS*            cons                /**< added constraint */
   )
{
   CONSSETCHG* conssetchg;
   int naddedconss;

   assert(conssetchgdyn != NULL);
   assert(node != NULL);
   assert(conssetchgdyn->conssetchg == &node->conssetchg);
   assert(cons != NULL);
   assert(cons->node == NULL);
   assert(cons->arraypos == -1);

   naddedconss = (*conssetchgdyn->conssetchg == NULL ? 0 : (*conssetchgdyn->conssetchg)->naddedconss);
   CHECK_OKAY( ensureAddedconssSize(conssetchgdyn, memhdr, set, naddedconss+1) );

   conssetchg = *conssetchgdyn->conssetchg;
   assert(conssetchg != NULL);

   conssetchg->addedconss[conssetchg->naddedconss] = cons;
   cons->node = node;
   cons->arraypos = conssetchg->naddedconss;
   conssetchg->naddedconss++;

   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

/** adds constraint disabling to constraint set changes, and captures constraint */
RETCODE SCIPconssetchgdynAddDisabledCons(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamically sized constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< disabled constraint */
   )
{
   CONSSETCHG* conssetchg;
   int ndisabledconss;

   assert(conssetchgdyn != NULL);
   assert(conssetchgdyn->conssetchg != NULL);
   assert(cons != NULL);

   ndisabledconss = (*conssetchgdyn->conssetchg == NULL ? 0 : (*conssetchgdyn->conssetchg)->ndisabledconss);
   CHECK_OKAY( ensureDisabledconssSize(conssetchgdyn, memhdr, set, ndisabledconss+1) );

   conssetchg = *conssetchgdyn->conssetchg;
   assert(conssetchg != NULL);

   conssetchg->disabledconss[conssetchg->ndisabledconss] = cons;
   conssetchg->ndisabledconss++;

   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

/** gets pointer to constraint set change data the dynamic size information references */
CONSSETCHG** SCIPconssetchgdynGetConssetchgPtr(
   CONSSETCHGDYN*   conssetchgdyn       /**< dynamically sized constraint set change data structure */
   )
{
   assert(conssetchgdyn != NULL);

   return conssetchgdyn->conssetchg;
}


