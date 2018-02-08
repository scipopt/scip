/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   interrupt.c
 * @brief  methods and datastructures for catching the user CTRL-C interrupt
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <sys/types.h>
#include <stdlib.h>
#include <signal.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/interrupt.h"


static volatile
int                      ninterrupts = 0;    /**< static variable counting the number of CTRL-C interrupts */
static volatile
int                      nterms = 0;         /**< static variable counting the number of times that the process received a SIGTERM signal */


#ifdef NO_SIGACTION
typedef void (*SigHdlr)(int);

/** CTRL-C interrupt data */
struct SCIP_Interrupt
{
   SigHdlr               oldsiginthdlr;      /**< old CTRL-C interrupt handler */
   SigHdlr               oldsigtermhdlr;     /**< old SIGTERM handler */
   int                   nuses;              /**< number of times, the interrupt is captured */
};

#else

/** CTRL-C interrupt data */
struct SCIP_Interrupt
{
   struct sigaction      oldsigintaction;    /**< old CTRL-C interrupt handler */
   struct sigaction      oldsigtermaction;   /**< old SIGTERM handler (of a parent process) */
   int                   nuses;              /**< number of times, the interrupt is captured */
};
#endif

/** interrupt handler for CTRL-C interrupts */
static
void interruptHandler(
   int                   signum              /**< interrupt signal number */
   )
{  switch (signum) {
      case SIGINT:
         ninterrupts++;
         if( ninterrupts >= 5 )
         {
            printf("pressed CTRL-C %d times. forcing termination.\n", ninterrupts);
            exit(1);
         }
         else
         {
            printf("pressed CTRL-C %d times (5 times for forcing termination)\n", ninterrupts);
         }

         break;
      case SIGTERM:
         nterms++;
         break;
      default:
         break;
   }

}

/** creates a CTRL-C interrupt data */
SCIP_RETCODE SCIPinterruptCreate(
   SCIP_INTERRUPT**      interrupt           /**< pointer to store the CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);

   SCIP_ALLOC( BMSallocMemory(interrupt) );
   (*interrupt)->nuses = 0;

   return SCIP_OKAY;
}

/** frees a CTRL-C interrupt data */
void SCIPinterruptFree(
   SCIP_INTERRUPT**      interrupt           /**< pointer to the CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);

   BMSfreeMemory(interrupt);
}

/** captures the CTRL-C interrupt to call the SCIP's own interrupt handler */
void SCIPinterruptCapture(
   SCIP_INTERRUPT*       interrupt           /**< CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);
   assert(interrupt->nuses >= 0);

   if( interrupt->nuses == 0 )
   {
#ifdef NO_SIGACTION
      interrupt->oldsiginthdlr = signal(SIGINT, interruptHandler);
      interrupt->oldsigtermhdlr = signal(SIGTERM, interruptHandler);
#else
      struct sigaction newaction;

      /* initialize new signal action */
      newaction.sa_handler = interruptHandler;
      newaction.sa_flags = 0;
      (void)sigemptyset(&newaction.sa_mask);

      /* set new signal action, and remember old one */
      (void)sigaction(SIGINT, &newaction, &interrupt->oldsigintaction);

      (void)sigaction(SIGTERM, &newaction, &interrupt->oldsigtermaction);
#endif

      ninterrupts = 0;
      nterms = 0;
   }
   interrupt->nuses++;
}

/** releases the CTRL-C interrupt and restores the old interrupt handler */
void SCIPinterruptRelease(
   SCIP_INTERRUPT*       interrupt           /**< CTRL-C interrupt data */
   )
{
   assert(interrupt != NULL);
   assert(interrupt->nuses >= 1);

   interrupt->nuses--;
   if( interrupt->nuses == 0 )
   {
#ifdef NO_SIGACTION
      (void)signal(SIGINT, interrupt->oldsiginthdlr);
      (void)signal(SIGTERM, interrupt->oldsigtermhdlr);
#else
      (void)sigaction(SIGINT, &interrupt->oldsigintaction, NULL);
      (void)sigaction(SIGTERM, &interrupt->oldsigtermaction, NULL);
#endif
   }
}

/** returns whether the user interrupted by pressing CTRL-C */
SCIP_Bool SCIPinterrupted(
   void
   )
{
   return (ninterrupts > 0);
}

/** returns whether the process has received a SIGTERM */
SCIP_Bool SCIPterminated(
   void
   )
{
   return (nterms > 0);
}

/** resets the number of interrupts to 0 */
void SCIPresetInterrupted(
   void
   )
{
   ninterrupts = 0;
   nterms = 0;
}

