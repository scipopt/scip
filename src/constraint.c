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

/**@file   constraint.c
 * @brief  datastructures and methods for managing constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "constraint.h"


/** constraint handler */
struct ConsHdlr
{
   char*            name;               /**< name of constraint handler */
   DECL_CONSINIT((*consinit));          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit));          /**< deinitialise constraint handler */
   DECL_CONSFREE((*consfree));          /**< frees specific constraint data */
   DECL_CONSCHCK((*conschck));          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop));          /**< propagate variable domains */
};

/** constraint data structure */
struct Cons
{
   CONSHDLR*        conshdlr;           /**< constraint handler for this constraint */
   CONSDATA*        consdata;           /**< data for this specific constraint */
   int              numuses;            /**< number of times, this constraint is referenced */
   unsigned int     original:1;         /**< TRUE iff constraint belongs to the original problem formulation */
   unsigned int     model:1;            /**< TRUE iff constraint is necessary for feasibility */
};

/** linked list of constraints */
struct ConsList
{
   CONS*            cons;               /**< pointer to constraint data structure */
   CONSLIST*        next;               /**< next list entry */
};



RETCODE SCIPconsCreate(                 /**< creates a constraint */
   CONS**           cons,               /**< pointer to constraint */
   MEM*             mem,                /**< block memory buffers */
   Bool             original,           /**< belongs constraint to the original problem formulation? */
   Bool             model,              /**< is constraint necessary for feasibility? */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata            /**< data for this specific constraint */
   )
{
   assert(cons != NULL);
   assert(mem != NULL);
   assert(conshdlr != NULL);

   ALLOC_OKAY( allocBlockMemory(mem->consmem, *cons) );
   (*cons)->conshdlr = conshdlr;
   (*cons)->consdata = consdata;
   (*cons)->numuses = 0;
   (*cons)->original = original;
   (*cons)->model = model;

   return SCIP_OKAY;
}

void SCIPconsFree(                      /**< frees a constraint */
   CONS**           cons,               /**< constraint to free */
   MEM*             mem                 /**< block memory buffers */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);
   assert((*cons)->conshdlr != NULL);
   assert((*cons)->numuses == 0);
   assert(mem != NULL);

   (*cons)->conshdlr->consfree((*cons)->conshdlr, mem->consmem, (*cons)->consdata);
   freeBlockMemory(mem->consmem, *cons);
}

void SCIPconsCapture(                   /**< increases usage counter of constraint */
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(cons->numuses >= 0);

   cons->numuses++;
}

void SCIPconsRelease(                   /**< decreases usage counter of constraint, and frees memory if necessary */
   CONS**           cons,               /**< pointer to constraint */
   MEM*             mem                 /**< block memory buffers */
   )
{
   assert(mem != NULL);
   assert(cons != NULL);
   assert(*cons != NULL);
   assert((*cons)->numuses >= 1);

   (*cons)->numuses--;
   if( (*cons)->numuses == 0 )
      SCIPconsFree(cons, mem);
}

RETCODE SCIPconslistAdd(                /**< adds constraint to a list of constraints and captures it */
   CONSLIST**       conslist,           /**< constraint list to extend */
   MEM*             mem,                /**< block memory buffers */
   CONS*            cons                /**< constraint to add */
   )
{
   CONSLIST* newlist;

   assert(conslist != NULL);
   assert(mem != NULL);
   assert(cons != NULL);
   
   ALLOC_OKAY( allocBlockMemory(mem->consmem, newlist) );
   newlist->cons = cons;
   newlist->next = *conslist;
   *conslist = newlist;

   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

void SCIPconslistFreePart(              /**< partially unlinks and releases the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEM*             mem,                /**< block memory buffers */
   CONSLIST*        firstkeep           /**< first constraint list entry to keep */
   )
{
   CONSLIST* next;

   assert(conslist != NULL);
   assert(mem != NULL);
   assert(firstkeep != NULL);
   
   while(*conslist != NULL && *conslist != firstkeep)
   {
      SCIPconsRelease(&((*conslist)->cons), mem);
      next = (*conslist)->next;
      freeBlockMemory(mem->consmem, *conslist);
      *conslist = next;
   }
   assert(*conslist == firstkeep); /* firstkeep should be part of conslist */
}

void SCIPconslistFree(                  /**< unlinks and releases all the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEM*             mem                 /**< block memory buffers */
   )
{
   SCIPconslistFreePart(conslist, mem, NULL);
}
