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
   CONS**           conss;              /**< array with active constraints of this handler; model constraints at first */
   int              consssize;          /**< size of constraints array (available slots in conss array) */
   int              nconss;             /**< number of active constraints (used slots in conss array) */
   int              nmodelconss;        /**< number of model constraints (stored in the first positions in conss array) */
   CONS**           probconss;          /**< array with model constraints of the initial (transformed) problem */
   int              probconsssize;      /**< available slots in probconss array */
   int              nprobconss;         /**< number of initial model constraints */
   unsigned int     needscons:1;        /**< should the constraint handler be skipped, if no constraints are available? */
   unsigned int     initialized:1;      /**< is constraint handler initialized? */
};

/** constraint data structure */
struct Cons
{
   char*            name;               /**< name of the constraint */
   CONSHDLR*        conshdlr;           /**< constraint handler for this constraint */
   CONSDATA*        consdata;           /**< data for this specific constraint */
   int              nuses;              /**< number of times, this constraint is referenced */
   unsigned int     model:1;            /**< TRUE iff constraint is necessary for feasibility */
   unsigned int     original:1;         /**< TRUE iff constraint belongs to original problem */
   unsigned int     active:1;           /**< TRUE iff constraint is active in the active node */
   unsigned int     arraypos:29;        /**< position of constraint in the constraint array of the handler */
};




/*
 * dynamic memory arrays
 */


/** resizes conss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsureConssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->consssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&conshdlr->conss, newsize) );
      conshdlr->consssize = newsize;
   }
   assert(num <= conshdlr->consssize);

   return SCIP_OKAY;
}

/** resizes probconss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsureProbconssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of node slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->probconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&conshdlr->probconss, newsize) );
      conshdlr->probconsssize = newsize;
   }
   assert(num <= conshdlr->probconsssize);

   return SCIP_OKAY;
}


/*
 * Constraint handler methods
 */

/** activates and adds constraint to constraint handler's constraint array */
static
RETCODE conshdlrAddCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nmodelconss >= 0);
   assert(conshdlr->nmodelconss <= conshdlr->nconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->active);
   assert(cons->arraypos == 0);

   debugMessage("add constraint <%s> to constraint handler <%s>\n", cons->name, conshdlr->name);
   CHECK_OKAY( conshdlrEnsureConssMem(conshdlr, set, conshdlr->nconss+1) );

   insertpos = conshdlr->nconss;
   if( cons->model )
   {
      if( conshdlr->nconss > conshdlr->nmodelconss )
      {
         assert(conshdlr->conss[conshdlr->nmodelconss]->arraypos == conshdlr->nmodelconss);
         conshdlr->conss[insertpos] = conshdlr->conss[conshdlr->nmodelconss];
         conshdlr->conss[insertpos]->arraypos = insertpos;
         insertpos = conshdlr->nmodelconss;
      }
      conshdlr->nmodelconss++;
   }
   conshdlr->conss[insertpos] = cons;
   cons->arraypos = insertpos;
   cons->active = TRUE;
   conshdlr->nconss++;
   debugMessage(" -> inserted at position %d (nmodelconss=%d, nconss=%d)\n", 
      insertpos, conshdlr->nmodelconss, conshdlr->nconss);

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
   assert(conshdlr->nconss > 0);
   assert(conshdlr->conss != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->arraypos < conshdlr->nconss);
   assert(cons->model || cons->arraypos >= conshdlr->nmodelconss);
   assert(!cons->model || cons->arraypos < conshdlr->nmodelconss);

   debugMessage("delete constraint <%s> from constraint handler <%s>\n", cons->name, conshdlr->name);
   debugMessage(" -> delete from position %d (nmodelconss=%d, nconss=%d)\n", 
      cons->arraypos, conshdlr->nmodelconss, conshdlr->nconss);

   delpos = cons->arraypos;
   if( cons->model )
   {
      assert(conshdlr->nmodelconss > 0);
      assert(delpos < conshdlr->nmodelconss);
      if( delpos < conshdlr->nmodelconss-1 )
      {
         conshdlr->conss[delpos] = conshdlr->conss[conshdlr->nmodelconss-1];
         conshdlr->conss[delpos]->arraypos = delpos;
         delpos = conshdlr->nmodelconss-1;
      }
      conshdlr->nmodelconss--;
   }
   assert(delpos < conshdlr->nconss);
   if( delpos < conshdlr->nconss-1 )
   {
      conshdlr->conss[delpos] = conshdlr->conss[conshdlr->nconss-1];
      conshdlr->conss[delpos]->arraypos = delpos;
   }
   conshdlr->nconss--;
   cons->arraypos = 0;
   cons->active = FALSE;

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
   (*conshdlr)->needscons = needscons;
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
   (*conshdlr)->conss = NULL;
   (*conshdlr)->consssize = 0;
   (*conshdlr)->nconss = 0;
   (*conshdlr)->nmodelconss = 0;
   (*conshdlr)->probconss = NULL;
   (*conshdlr)->probconsssize = 0;
   (*conshdlr)->nprobconss = 0;
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
      CHECK_OKAY( (*conshdlr)->consfree(*conshdlr, scip) );
   }

   freeMemoryArray(&(*conshdlr)->name);
   freeMemoryArray(&(*conshdlr)->desc);
   freeMemoryArrayNull(&(*conshdlr)->conss);
   freeMemoryArrayNull(&(*conshdlr)->probconss);
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
      CHECK_OKAY( conshdlr->consinit(conshdlr, scip) );
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
      CHECK_OKAY( conshdlr->consexit(conshdlr, scip) );
   }
   conshdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** calls separator method of constraint handler */
RETCODE SCIPconshdlrSeparate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   if( conshdlr->conssepa != NULL && (!conshdlr->needscons || conshdlr->nconss > 0) )
   {
      debugMessage("separate %d constraints of handler <%s>\n", conshdlr->nconss, conshdlr->name);
      CHECK_OKAY( conshdlr->conssepa(conshdlr, set->scip, conshdlr->conss, conshdlr->nconss, result) );
      if( *result != SCIP_SEPARATED
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
   else
      *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** calls enforcing method of constraint handler for LP solutions */
RETCODE SCIPconshdlrEnforceLPSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   if( conshdlr->consenlp != NULL && (!conshdlr->needscons || conshdlr->nmodelconss > 0) )
   {
      debugMessage("enforcing %d constraints of handler <%s> for LP solutions\n", conshdlr->nmodelconss, conshdlr->name);
      CHECK_OKAY( conshdlr->consenlp(conshdlr, set->scip, conshdlr->conss, conshdlr->nmodelconss, result) );
      if( *result != SCIP_CUTOFF
         && *result != SCIP_BRANCHED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_SEPARATED
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
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** calls enforcing method of constraint handler for pseudo solutions */
RETCODE SCIPconshdlrEnforcePseudoSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   if( conshdlr->consenps != NULL && (!conshdlr->needscons || conshdlr->nmodelconss > 0) )
   {
      debugMessage("enforcing %d constraints of handler <%s> for pseudo solutions\n", 
         conshdlr->nmodelconss, conshdlr->name);
      CHECK_OKAY( conshdlr->consenps(conshdlr, set->scip, conshdlr->conss, conshdlr->nmodelconss, result) );
      if( *result != SCIP_CUTOFF
         && *result != SCIP_BRANCHED
         && *result != SCIP_REDUCEDDOM
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
   else
      *result = SCIP_FEASIBLE;

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

   if( conshdlr->conschck != NULL && (!conshdlr->needscons || conshdlr->nprobconss > 0) )
   {
      debugMessage("checking %d constraints of handler <%s>\n", conshdlr->nprobconss, conshdlr->name);
      CHECK_OKAY( conshdlr->conschck(conshdlr, set->scip, conshdlr->probconss, conshdlr->nprobconss, 
                     sol, chckintegrality, chcklprows, result) );
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
   else
      *result = SCIP_FEASIBLE;

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

   if( conshdlr->consprop != NULL
      && (!conshdlr->needscons || conshdlr->nconss > 0)
      && (actdepth == -1 || (conshdlr->propfreq > 0 && actdepth % conshdlr->propfreq == 0)) )
   {
      debugMessage("propagating %d constraints of handler <%s>\n", conshdlr->nconss, conshdlr->name);
      CHECK_OKAY( conshdlr->consprop(conshdlr, set->scip, conshdlr->conss, conshdlr->nconss, result) );
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
   else
      *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** adds constraint to constraint handler's problem constraint array */
RETCODE SCIPconshdlrAddProbCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< model constraint of initial problem to add */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->model);

   debugMessage("add problem constraint <%s> to constraint handler <%s>\n", cons->name, conshdlr->name);
   CHECK_OKAY( conshdlrEnsureProbconssMem(conshdlr, set, conshdlr->nprobconss+1) );

   conshdlr->probconss[conshdlr->nprobconss] = cons;
   conshdlr->nprobconss++;

   return SCIP_OKAY;
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

/** gets constraints array of constraint handler */
CONS** SCIPconshdlrGetConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->conss;
}

/** gets number of constraints in constraints array of constraint handler */
int SCIPconshdlrGetNConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nconss;
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

/** creates and captures a constraint */
RETCODE SCIPconsCreate(
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             model,              /**< is constraint necessary for feasibility? */
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
   (*cons)->nuses = 0;
   (*cons)->model = model;
   (*cons)->original = original;
   (*cons)->active = FALSE;
   (*cons)->arraypos = 0;

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
      CHECK_OKAY( (*cons)->conshdlr->consdele((*cons)->conshdlr, set->scip, &(*cons)->consdata) );
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
      CHECK_OKAY( origcons->conshdlr->constran(origcons->conshdlr, set->scip, origcons->consdata, &consdata) );
   }

   /* create new constraint with transformed data */
   CHECK_OKAY( SCIPconsCreate(transcons, memhdr, origcons->name, origcons->conshdlr, consdata, origcons->model, FALSE) );

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

/** returns TRUE iff constraint is necessary for feasibility */
Bool SCIPconsIsModel(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->model;
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
 * Constraint list methods
 */

/** adds constraint to a list of constraints and captures it */
RETCODE SCIPconslistAdd(
   CONSLIST**       conslist,           /**< constraint list to extend */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   )
{
   CONSLIST* newlist;

   assert(conslist != NULL);
   assert(memhdr != NULL);
   assert(cons != NULL);
   
   ALLOC_OKAY( allocBlockMemory(memhdr, &newlist) );
   newlist->cons = cons;
   newlist->next = *conslist;
   *conslist = newlist;

   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

/** partially unlinks and releases the constraints in the list */
RETCODE SCIPconslistFreePart(
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   CONSLIST*        firstkeep           /**< first constraint list entry to keep */
   )
{
   CONSLIST* next;

   assert(conslist != NULL);
   assert(memhdr != NULL);
   
   while(*conslist != NULL && *conslist != firstkeep)
   {
      CHECK_OKAY( SCIPconsRelease(&(*conslist)->cons, memhdr, set) );
      next = (*conslist)->next;
      freeBlockMemory(memhdr, conslist);
      *conslist = next;
   }
   assert(*conslist == firstkeep); /* firstkeep should be part of conslist */

   return SCIP_OKAY;
}

/** unlinks and releases all the constraints in the list */
RETCODE SCIPconslistFree(
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   )
{
   return SCIPconslistFreePart(conslist, memhdr, set, NULL);
}

