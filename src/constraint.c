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
   MEMHDR*          memhdr,             /**< block memory */
   Bool             original,           /**< belongs constraint to the original problem formulation? */
   Bool             model,              /**< is constraint necessary for feasibility? */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata            /**< data for this specific constraint */
   )
{
   assert(cons != NULL);
   assert(memhdr != NULL);
   assert(conshdlr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, *cons) );
   (*cons)->conshdlr = conshdlr;
   (*cons)->consdata = consdata;
   (*cons)->numuses = 0;
   (*cons)->original = original;
   (*cons)->model = model;

   return SCIP_OKAY;
}

void SCIPconsFree(                      /**< frees a constraint */
   CONS**           cons,               /**< constraint to free */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);
   assert((*cons)->conshdlr != NULL);
   assert((*cons)->numuses == 0);
   assert(memhdr != NULL);

   (*cons)->conshdlr->consfree((*cons)->conshdlr, memhdr, (*cons)->consdata);
   freeBlockMemory(memhdr, *cons);
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
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(memhdr != NULL);
   assert(cons != NULL);
   assert(*cons != NULL);
   assert((*cons)->numuses >= 1);

   (*cons)->numuses--;
   if( (*cons)->numuses == 0 )
      SCIPconsFree(cons, memhdr);
}

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

void SCIPconslistFreePart(              /**< partially unlinks and releases the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr,             /**< block memory */
   CONSLIST*        firstkeep           /**< first constraint list entry to keep */
   )
{
   CONSLIST* next;

   assert(conslist != NULL);
   assert(memhdr != NULL);
   assert(firstkeep != NULL);
   
   while(*conslist != NULL && *conslist != firstkeep)
   {
      SCIPconsRelease(&((*conslist)->cons), memhdr);
      next = (*conslist)->next;
      freeBlockMemory(memhdr, *conslist);
      *conslist = next;
   }
   assert(*conslist == firstkeep); /* firstkeep should be part of conslist */
}

void SCIPconslistFree(                  /**< unlinks and releases all the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   SCIPconslistFreePart(conslist, memhdr, NULL);
}
