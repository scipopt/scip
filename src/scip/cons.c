/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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
   DECL_CONSFREE((*consfree));          /**< destructor of constraint handler */
   DECL_CONSINIT((*consinit));          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit));          /**< deinitialise constraint handler */
   DECL_CONSDELE((*consdele));          /**< frees specific constraint data */
   DECL_CONSTRAN((*constran));          /**< transforms constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA((*conssepa));          /**< separates cutting planes */
   DECL_CONSENFO((*consenfo));          /**< enforcing constraints */
   DECL_CONSCHCK((*conschck));          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop));          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata;       /**< constraint handler data */
   CONS**           conss;              /**< array with active constraints of this handler; model constraints at first */
   int              consssize;          /**< size of constraints array (available slots in conss array) */
   int              nconss;             /**< number of active constraints (used slots in conss array) */
   int              nmodelconss;        /**< number of model constraints (stored in the first positions in conss array) */
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

static
RETCODE conshdlrEnsureConssMem(         /**< resizes conss array to be able to store at least num constraints */
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
      ALLOC_OKAY( reallocMemoryArray(conshdlr->conss, newsize) );
      conshdlr->consssize = newsize;
   }
   assert(num <= conshdlr->consssize);

   return SCIP_OKAY;
}


/*
 * Constraint handler methods
 */

static
RETCODE conshdlrAddCons(                /**< activates and adds constraint to constraint handler's constraint array */
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

static
RETCODE conshdlrDelCons(                /**< deactivates and removes constraint from constraint handler's conss array */
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

RETCODE SCIPconshdlrCreate(             /**< creates a constraint handler */
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              chckpriority,       /**< priority of the constraint handler for checking infeasibility */
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE((*consfree)),          /**< destructor of constraint handler */
   DECL_CONSINIT((*consinit)),          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit)),          /**< deinitialise constraint handler */
   DECL_CONSDELE((*consdele)),          /**< free specific constraint data */
   DECL_CONSTRAN((*constran)),          /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA((*conssepa)),          /**< separate cutting planes */
   DECL_CONSENFO((*consenfo)),          /**< enforcing constraints */
   DECL_CONSCHCK((*conschck)),          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop)),          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   ALLOC_OKAY( allocMemory(*conshdlr) );
   ALLOC_OKAY( duplicateMemoryArray((*conshdlr)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray((*conshdlr)->desc, desc, strlen(desc)+1) );
   (*conshdlr)->sepapriority = sepapriority;
   (*conshdlr)->enfopriority = enfopriority;
   (*conshdlr)->chckpriority = chckpriority;
   (*conshdlr)->needscons = needscons;
   (*conshdlr)->consfree = consfree;
   (*conshdlr)->consinit = consinit;
   (*conshdlr)->consexit = consexit;
   (*conshdlr)->consdele = consdele;
   (*conshdlr)->constran = constran;
   (*conshdlr)->conssepa = conssepa;
   (*conshdlr)->consenfo = consenfo;
   (*conshdlr)->conschck = conschck;
   (*conshdlr)->consprop = consprop;
   (*conshdlr)->conshdlrdata = conshdlrdata;
   (*conshdlr)->conss = NULL;
   (*conshdlr)->consssize = 0;
   (*conshdlr)->nconss = 0;
   (*conshdlr)->nmodelconss = 0;
   (*conshdlr)->initialized = FALSE;

   return SCIP_OKAY;
}

RETCODE SCIPconshdlrFree(               /**< calls destructor and frees memory of constraint handler */
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conshdlr != NULL);
   assert(*conshdlr != NULL);
   assert(!(*conshdlr)->initialized);

   /* call destructor of constraint handler */
   if( (*conshdlr)->consfree != NULL )
   {
      CHECK_OKAY( (*conshdlr)->consfree(*conshdlr, scip) );
   }

   freeMemoryArray((*conshdlr)->name);
   freeMemoryArray((*conshdlr)->desc);
   freeMemoryArrayNull((*conshdlr)->conss);
   freeMemory(*conshdlr);

   return SCIP_OKAY;
}

RETCODE SCIPconshdlrInit(               /**< initializes constraint handler */
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

   assert(conshdlr->nconss == 0);
   if( conshdlr->consinit != NULL )
   {
      CHECK_OKAY( conshdlr->consinit(conshdlr, scip) );
   }
   conshdlr->initialized = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPconshdlrExit(               /**< calls exit method of constraint handler */
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

RETCODE SCIPconshdlrSeparate(           /**< calls separator method of constraint handler */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(result != NULL);

   if( conshdlr->conssepa != NULL && (!conshdlr->needscons || conshdlr->nconss > 0) )
   {
      debugMessage("separate constraints of handler <%s>\n", conshdlr->name);
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

RETCODE SCIPconshdlrEnforce(            /**< calls enforcing method of constraint handler */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   Bool             lpvalid,            /**< is the LP being processed at the current node? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(result != NULL);

   if( conshdlr->consenfo != NULL && (!conshdlr->needscons || conshdlr->nmodelconss > 0) )
   {
      debugMessage("enforcing constraints of handler <%s>\n", conshdlr->name);
      CHECK_OKAY( conshdlr->consenfo(conshdlr, set->scip, conshdlr->conss, conshdlr->nmodelconss, lpvalid, result) );
      if( *result != SCIP_BRANCHED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_SEPARATED
         && *result != SCIP_INFEASIBLE
         && *result != SCIP_FEASIBLE )
      {
         char s[255];
         sprintf(s, "enforcing method of constraint handler <%s> returned invalid result <%d>", 
            conshdlr->name, *result);
         errorMessage(s);
         return SCIP_INVALIDRESULT;
      }
   }
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

RETCODE SCIPconshdlrCheck(              /**< calls feasibility check method of constraint handler */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   SOL*             sol,                /**< primal CIP solution */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(result != NULL);

   if( conshdlr->conschck != NULL && (!conshdlr->needscons || conshdlr->nmodelconss > 0) )
   {
      debugMessage("checking constraints of handler <%s>\n", conshdlr->name);
      CHECK_OKAY( conshdlr->conschck(conshdlr, set->scip, conshdlr->conss, conshdlr->nmodelconss, sol, result) );
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

RETCODE SCIPconshdlrPropagate(          /**< calls propagation method of constraint handler */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(result != NULL);

   if( conshdlr->consprop != NULL && (!conshdlr->needscons || conshdlr->nconss > 0) )
   {
      debugMessage("propagating constraints of handler <%s>\n", conshdlr->name);
      CHECK_OKAY( conshdlr->consprop(conshdlr, set->scip, conshdlr->conss, conshdlr->nconss, result) );
      if( *result != SCIP_REDUCEDDOM
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

const char* SCIPconshdlrGetName(        /**< gets name of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->name;
}

CONSHDLRDATA* SCIPconshdlrGetData(      /**< gets user data of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->conshdlrdata;
}

void SCIPconshdlrSetData(               /**< sets user data of constraint handler; user has to free old data in advance! */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   CONSHDLRDATA*    conshdlrdata        /**< new constraint handler user data */
   )
{
   assert(conshdlr != NULL);

   conshdlr->conshdlrdata = conshdlrdata;
}

CONS** SCIPconshdlrGetConss(            /**< gets constraints array of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->conss;
}

int SCIPconshdlrGetNConss(              /**< gets number of constraints in constraints array of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nconss;
}

Bool SCIPconshdlrIsInitialized(         /**< is constraint handler initialized? */
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->initialized;
}




/*
 * Constraint methods
 */

RETCODE SCIPconsCreate(                 /**< creates and captures a constraint */
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
   assert(conshdlr->initialized);

   /* create constraint data */
   ALLOC_OKAY( allocBlockMemory(memhdr, *cons) );
   ALLOC_OKAY( duplicateBlockMemoryArray(memhdr, (*cons)->name, name, strlen(name)+1) );
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

RETCODE SCIPconsFree(                   /**< frees a constraint */
   CONS**           cons,               /**< constraint to free */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);
   assert((*cons)->nuses == 0);
   assert((*cons)->conshdlr != NULL);
   assert((*cons)->conshdlr->initialized);
   assert(memhdr != NULL);
   assert(set != NULL);

   /* free constraint data */
   if( (*cons)->conshdlr->consdele != NULL )
   {
      CHECK_OKAY( (*cons)->conshdlr->consdele((*cons)->conshdlr, set->scip, &(*cons)->consdata) );
   }
   freeBlockMemoryArray(memhdr, (*cons)->name, strlen((*cons)->name)+1);
   freeBlockMemory(memhdr, *cons);

   return SCIP_OKAY;
}

void SCIPconsCapture(                   /**< increases usage counter of constraint */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(cons->nuses >= 0);

   debugMessage("capture constraint <%s> with nuses=%d\n", cons->name, cons->nuses);
   cons->nuses++;
}

RETCODE SCIPconsRelease(                /**< decreases usage counter of constraint, and frees memory if necessary */
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

RETCODE SCIPconsActivate(               /**< activates constraint */
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

RETCODE SCIPconsDeactivate(             /**< deactivates constraint */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(cons->active);
   
   CHECK_OKAY( conshdlrDelCons(cons->conshdlr, cons) );
   assert(!cons->active);

   return SCIP_OKAY;
}

RETCODE SCIPconsTransform(              /**< copies original constraint into transformed constraint, that is captured */
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

const char* SCIPconsGetName(            /**< returns the name of the constraint */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->name;
}

CONSHDLR* SCIPconsGetConsHdlr(          /**< returns the constraint handler of the constraint */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->conshdlr;
}

CONSDATA* SCIPconsGetConsdata(          /**< returns the constraint data field of the constraint */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->consdata;
}

Bool SCIPconsIsOriginal(                /**< returns TRUE iff constraint is belonging to original problem */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->original;
}

Bool SCIPconsIsModel(                   /**< returns TRUE iff constraint is necessary for feasibility */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->model;
}


/*
 * Hash functions
 */

DECL_HASHGETKEY(SCIPhashGetKeyCons)     /**< gets the key (i.e. the name) of the given constraint */
{
   CONS* cons = (CONS*)elem;

   assert(cons != NULL);
   return cons->name;
}



/*
 * Constraint list methods
 */

RETCODE SCIPconslistAdd(                /**< adds constraint to a list of constraints and captures it */
   CONSLIST**       conslist,           /**< constraint list to extend */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   )
{
   CONSLIST* newlist;

   assert(conslist != NULL);
   assert(memhdr != NULL);
   assert(cons != NULL);
   
   ALLOC_OKAY( allocBlockMemory(memhdr, newlist) );
   newlist->cons = cons;
   newlist->next = *conslist;
   *conslist = newlist;

   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

RETCODE SCIPconslistFreePart(           /**< partially unlinks and releases the constraints in the list */
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
      freeBlockMemory(memhdr, *conslist);
      *conslist = next;
   }
   assert(*conslist == firstkeep); /* firstkeep should be part of conslist */

   return SCIP_OKAY;
}

RETCODE SCIPconslistFree(               /**< unlinks and releases all the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   )
{
   return SCIPconslistFreePart(conslist, memhdr, set, NULL);
}

