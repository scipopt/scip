/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: interrupt.c,v 1.14 2005/01/21 09:16:55 bzfpfend Exp $"

/**@file   interrupt.c
 * @brief  methods and datastructures for catching the user CTRL-C interrupt
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <sys/types.h>
#include <stdlib.h>
#include <signal.h>

#include "def.h"
#include "memory.h"
#include "interrupt.h"


static volatile
int                 ninterrupts = 0;    /**< static variable counting the number of CTRL-C interrupts */


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
{  /*lint --e{715}*/
   ninterrupts++;
   if( ninterrupts >= 5 )
   {
      printf("pressed CTRL-C %d times. forcing termination.\n", ninterrupts);
      exit(1);
   }
   else
      printf("pressed CTRL-C %d times (5 times for forcing termination)\n", ninterrupts);
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
   INTERRUPT*       interrupt           /**< CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);
   assert(interrupt->nuses >= 0);

   if( interrupt->nuses == 0 )
   {
      struct sigaction newaction;
      
      /* initialize new signal action */
      newaction.sa_handler = interruptHandler;
      newaction.sa_flags = 0;
      (void)sigemptyset(&newaction.sa_mask);
      
      /* set new signal action, and remember old one */
      (void)sigaction(SIGINT, &newaction, &interrupt->oldsigaction);
   }
   interrupt->nuses++;
   ninterrupts = 0;
}

/** releases the CTRL-C interrupt and restores the old interrupt handler */
void SCIPinterruptRelease(
   INTERRUPT*       interrupt           /**< CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);
   assert(interrupt->nuses >= 1);

   if( interrupt->nuses > 0 )
   {
      interrupt->nuses--;
      if( interrupt->nuses == 0 )
         (void)sigaction(SIGINT, &interrupt->oldsigaction, NULL);
   }
}

/** returns whether the user interrupted by pressing CTRL-C */
Bool SCIPinterrupted(
   void
   )
{
   return (ninterrupts > 0);
}

