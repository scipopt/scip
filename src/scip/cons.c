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
   const char*      name;               /**< name of constraint handler */
   const char*      desc;               /**< description of constraint handler */
   DECL_CONSINIT((*consinit));          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit));          /**< deinitialise constraint handler */
   DECL_CONSFREE((*consfree));          /**< frees specific constraint data */
   DECL_CONSTRAN((*constran));          /**< transforms constraint data into data belonging to the transformed problem */
   DECL_CONSCHCK((*conschck));          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop));          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata;       /**< constraint handler data */
   unsigned int     initialized:1;      /**< is constraint handler initialized? */
   unsigned int     nconstraints:31;    /**< number of existing constraints of this type */
};

/** constraint data structure */
struct Cons
{
   CONSHDLR*        conshdlr;           /**< constraint handler for this constraint */
   CONSDATA*        consdata;           /**< data for this specific constraint */
   unsigned int     ismodel:1;          /**< TRUE iff constraint is necessary for feasibility */
};



RETCODE SCIPconshdlrCreate(             /**< creates a constraint handler */
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   DECL_CONSINIT((*consinit)),          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit)),          /**< deinitialise constraint handler */
   DECL_CONSFREE((*consfree)),          /**< frees specific constraint data */
   DECL_CONSTRAN((*constran)),          /**< transforms constraint data into data belonging to the transformed problem */
   DECL_CONSCHCK((*conschck)),          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop)),          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   ALLOC_OKAY( allocMemory(*conshdlr) );
   ALLOC_OKAY( duplicateMemoryArray((*conshdlr)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray((*conshdlr)->desc, desc, strlen(desc)+1) );
   (*conshdlr)->consinit = consinit;
   (*conshdlr)->consexit = consexit;
   (*conshdlr)->consfree = consfree;
   (*conshdlr)->constran = constran;
   (*conshdlr)->conschck = conschck;
   (*conshdlr)->consprop = consprop;
   (*conshdlr)->conshdlrdata = conshdlrdata;
   (*conshdlr)->initialized = FALSE;
   (*conshdlr)->nconstraints = 0;

   return SCIP_OKAY;
}

RETCODE SCIPconshdlrFree(               /**< frees memory of constraint handler */
   CONSHDLR**       conshdlr            /**< pointer to constraint handler data structure */
   )
{
   assert(conshdlr != NULL);
   assert(*conshdlr != NULL);
   assert(!(*conshdlr)->initialized);

   freeMemoryArray((*conshdlr)->name);
   freeMemoryArray((*conshdlr)->desc);
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

   assert(conshdlr->nconstraints == 0);
   CHECK_OKAY( (*conshdlr->consinit)(conshdlr, scip) );
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

   CHECK_OKAY( (*conshdlr->consexit)(conshdlr, scip) );
   conshdlr->initialized = FALSE;

   return SCIP_OKAY;
}

const char* SCIPconshdlrGetName(        /**< gets name of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handlert */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->name;
}

Bool SCIPconshdlrIsInitialized(         /**< is constraint handler initialized? */
   CONSHDLR*        conshdlr            /**< constraint handlert */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->initialized;
}


RETCODE SCIPconsCreate(                 /**< creates a constraint */
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   Bool             ismodel,            /**< is constraint necessary for feasibility? */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata            /**< data for this specific constraint */
   )
{
   assert(cons != NULL);
   assert(memhdr != NULL);
   assert(conshdlr != NULL);
   assert(conshdlr->initialized);

   /* create constraint data */
   ALLOC_OKAY( allocBlockMemory(memhdr, *cons) );
   (*cons)->conshdlr = conshdlr;
   (*cons)->consdata = consdata;
   (*cons)->ismodel = ismodel;

   /* increase constraint counter of constraint handler */
   conshdlr->nconstraints++;

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
   assert((*cons)->conshdlr != NULL);
   assert((*cons)->conshdlr->initialized);
   assert((*cons)->conshdlr->nconstraints > 0);   
   assert(memhdr != NULL);
   assert(set != NULL);

   /* decrease constraint counter of constraint handler */
   (*cons)->conshdlr->nconstraints--;

   /* free constraint data */
   CHECK_OKAY( (*cons)->conshdlr->consfree((*cons)->conshdlr, set->scip, &(*cons)->consdata) );
   freeBlockMemory(memhdr, *cons);

   return SCIP_OKAY;
}

RETCODE SCIPconslistAdd(                /**< adds constraint to a list of constraints */
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

   return SCIP_OKAY;
}

RETCODE SCIPconslistFreePart(           /**< partially unlinks and frees the constraints in the list */
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
      CHECK_OKAY( SCIPconsFree(&(*conslist)->cons, memhdr, set) );
      next = (*conslist)->next;
      freeBlockMemory(memhdr, *conslist);
      *conslist = next;
   }
   assert(*conslist == firstkeep); /* firstkeep should be part of conslist */

   return SCIP_OKAY;
}

RETCODE SCIPconslistFree(               /**< unlinks and frees all the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   )
{
   return SCIPconslistFreePart(conslist, memhdr, set, NULL);
}

RETCODE SCIPconsTransform(              /**< copies original constraint into transformed constraint */
   CONS*            origcons,           /**< original constraint */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   CONS**           transcons           /**< pointer to transformed constraint */
   )
{
   CONSDATA* consdata;

   assert(origcons != NULL);
   assert(memhdr != NULL);
   assert(transcons != NULL);

   /* transform constraint data */
   CHECK_OKAY( origcons->conshdlr->constran(origcons->conshdlr, set->scip, origcons->consdata, &consdata) );

   /* create new constraint with transformed data */
   CHECK_OKAY( SCIPconsCreate(transcons, memhdr, origcons->ismodel, origcons->conshdlr, consdata) );

   return SCIP_OKAY;
}
