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

/**@file   interrupt.c
 * @brief  methods for catching the user CTRL-C interrupt
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <sys/types.h>
#include <stdlib.h>
#include <signal.h>

#include "interrupt.h"
#include "memory.h"


static volatile
Bool                interrupted = FALSE;/**< static variable set to TRUE in case of an CTRL-C interrupt */


/** CTRL-C interrupt data */
struct Interrupt
{
   struct sigaction oldsigaction;       /**< old CTRL-C interrupt handler */
   int              nuses;              /**< number of times, the interrupt is captured */
};


/** interrupt handler for CTRL-C interrupts */
static
void interruptHandler(
   int              signum              /**< interrupt signal number */
   )
{
   interrupted = TRUE;
}

/** creates a CTRL-C interrupt data */
RETCODE SCIPinterruptCreate(
   INTERRUPT**      interrupt           /**< pointer to store the CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);

   ALLOC_OKAY( allocMemory(interrupt) );
   (*interrupt)->nuses = 0;

   return SCIP_OKAY;
}

/** frees a CTRL-C interrupt data */
void SCIPinterruptFree(
   INTERRUPT**      interrupt           /**< pointer to the CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);

   freeMemory(interrupt);
}

/** captures the CTRL-C interrupt to call the SCIP's own interrupt handler */
void SCIPinterruptCapture(
   INTERRUPT*       interrupt,          /**< CTRL-C interrupt data */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(interrupt != NULL);
   assert(interrupt->nuses >= 0);
   assert(set != NULL);

   if( set->catchctrlc )
   {
      if( interrupt->nuses == 0 )
      {
         struct sigaction newaction;

         /* initialize new signal action */
         newaction.sa_handler = interruptHandler;
         newaction.sa_flags = 0;
         sigemptyset(&newaction.sa_mask);
      
         /* set new signal action, and remember old one */
         sigaction(SIGINT, &newaction, &interrupt->oldsigaction);
         interrupted = FALSE;
      }
      interrupt->nuses++;
   }
}

/** releases the CTRL-C interrupt and restores the old interrupt handler */
void SCIPinterruptRelease(
   INTERRUPT*       interrupt,          /**< CTRL-C interrupt data */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(interrupt != NULL);
   assert(!set->catchctrlc || interrupt->nuses >= 1);

   if( set->catchctrlc || interrupt->nuses > 0 )
   {
      interrupt->nuses--;
      if( interrupt->nuses == 0 )
         sigaction(SIGINT, &interrupt->oldsigaction, NULL);
   }
}

/** returns whether the user interrupted by pressing CTRL-C */
Bool SCIPinterrupted(
   void
   )
{
   return interrupted;
}

