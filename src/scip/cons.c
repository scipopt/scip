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
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons.c,v 1.107 2005/01/18 14:34:28 bzfpfend Exp $"

/**@file   cons.c
 * @brief  methods for constraints and constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "def.h"
#include "set.h"
#include "stat.h"
#include "clock.h"
#include "var.h"
#include "prob.h"
#include "tree.h"
#include "sepastore.h"
#include "scip.h"
#include "cons.h"

#ifndef NDEBUG
#include "struct_cons.h"
#endif




/*
 * dynamic memory arrays
 */


/** resizes conss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsureConssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
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

/** resizes sepaconss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsureSepaconssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
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
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
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

/** resizes checkconss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsureCheckconssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->checkconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&conshdlr->checkconss, newsize) );
      conshdlr->checkconsssize = newsize;
   }
   assert(num <= conshdlr->checkconsssize);

   return SCIP_OKAY;
}

/** resizes propconss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsurePropconssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
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

/** resizes updateconss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsureUpdateconssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);

   if( num > conshdlr->updateconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&conshdlr->updateconss, newsize) );
      conshdlr->updateconsssize = newsize;
   }
   assert(num <= conshdlr->updateconsssize);

   return SCIP_OKAY;
}




/*
 * Constraint handler methods
 */

/** returns whether the constraint's age exceeds the age limit */
static
Bool consExceedsAgelimit(
   CONS*            cons,               /**< constraint to check */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(set != NULL);

   return (set->cons_agelimit >= 0 && cons->age > set->cons_agelimit);
}

/** returns whether the constraint's age exceeds the obsolete age limit */
static
Bool consExceedsObsoleteage(
   CONS*            cons,               /**< constraint to check */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(set != NULL);

   return (set->cons_obsoleteage >= 0 && cons->age > set->cons_obsoleteage);
}

/** marks constraint to be obsolete; it will be moved to the last part of the constraint arrays, such that
 *  it is checked, enforced, separated, and propagated after the useful constraints
 */
static
RETCODE conshdlrMarkConsObsolete(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   CONS*            cons                /**< constraint to be marked obsolete */
   )
{
   CONS* tmpcons;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(!cons->obsolete);

   cons->obsolete = TRUE;
   
   if( cons->active )
   {
      if( cons->check )
      {
         assert(0 <= cons->checkconsspos && cons->checkconsspos < conshdlr->nusefulcheckconss);
         
         /* switch the last useful (non-obsolete) check constraint with this constraint */
         tmpcons = conshdlr->checkconss[conshdlr->nusefulcheckconss-1];
         assert(tmpcons->checkconsspos == conshdlr->nusefulcheckconss-1);
         
         conshdlr->checkconss[conshdlr->nusefulcheckconss-1] = cons;
         conshdlr->checkconss[cons->checkconsspos] = tmpcons;
         tmpcons->checkconsspos = cons->checkconsspos;
         cons->checkconsspos = conshdlr->nusefulcheckconss-1;
         
         conshdlr->nusefulcheckconss--;
      }
   }
   if( cons->enabled )
   {
      if( cons->separate )
      {
         assert(0 <= cons->sepaconsspos && cons->sepaconsspos < conshdlr->nusefulsepaconss);
         
         if( cons->sepaconsspos < conshdlr->lastnusefulsepaconss )
            conshdlr->lastnusefulsepaconss--;

         /* switch the last useful (non-obsolete) sepa constraint with this constraint */
         tmpcons = conshdlr->sepaconss[conshdlr->nusefulsepaconss-1];
         assert(tmpcons->sepaconsspos == conshdlr->nusefulsepaconss-1);
         
         conshdlr->sepaconss[conshdlr->nusefulsepaconss-1] = cons;
         conshdlr->sepaconss[cons->sepaconsspos] = tmpcons;
         tmpcons->sepaconsspos = cons->sepaconsspos;
         cons->sepaconsspos = conshdlr->nusefulsepaconss-1;
         
         conshdlr->nusefulsepaconss--;
      }
      if( cons->enforce )
      {
         assert(0 <= cons->enfoconsspos && cons->enfoconsspos < conshdlr->nusefulenfoconss);
         
         if( cons->enfoconsspos < conshdlr->lastnusefulenfoconss )
            conshdlr->lastnusefulenfoconss--;
         else
         {
            /* the constraint that becomes obsolete is not yet enforced on the current solution:
             * we have to make sure that it will be enforced the next time; this is not done, if the current
             * solution was already enforced and only enforcement on the additional constraints is performed
             * (because in this case, only the new useful constraints are enforced);
             * thus, we have to reset the enforcement counters in order to enforce all constraints again, especially
             * the now obsolete one;
             * this case should occur almost never, because a constraint that was not enforced in the last enforcement
             * is a newly added one, and it is very unlikely that this constraint will become obsolete before the next
             * enforcement call;
             * this reset is not performed for separation and propagation, because they are not vital for correctness
             */
            conshdlr->lastenfolplpcount = -1;
            conshdlr->lastenfolpdomchgcount = -1;
            conshdlr->lastenfopsdomchgcount = -1;
         }

         /* switch the last useful (non-obsolete) enfo constraint with this constraint */
         tmpcons = conshdlr->enfoconss[conshdlr->nusefulenfoconss-1];
         assert(tmpcons->enfoconsspos == conshdlr->nusefulenfoconss-1);
         
         conshdlr->enfoconss[conshdlr->nusefulenfoconss-1] = cons;
         conshdlr->enfoconss[cons->enfoconsspos] = tmpcons;
         tmpcons->enfoconsspos = cons->enfoconsspos;
         cons->enfoconsspos = conshdlr->nusefulenfoconss-1;
         
         conshdlr->nusefulenfoconss--;
      }
      if( cons->propagate && cons->propenabled )
      {
         assert(0 <= cons->propconsspos && cons->propconsspos < conshdlr->nusefulpropconss);
         
         if( cons->propconsspos < conshdlr->lastnusefulpropconss )
            conshdlr->lastnusefulpropconss--;

         /* switch the last useful (non-obsolete) prop constraint with this constraint */
         tmpcons = conshdlr->propconss[conshdlr->nusefulpropconss-1];
         assert(tmpcons->propconsspos == conshdlr->nusefulpropconss-1);
         
         conshdlr->propconss[conshdlr->nusefulpropconss-1] = cons;
         conshdlr->propconss[cons->propconsspos] = tmpcons;
         tmpcons->propconsspos = cons->propconsspos;
         cons->propconsspos = conshdlr->nusefulpropconss-1;
         
         conshdlr->nusefulpropconss--;
      }
   }

   return SCIP_OKAY;
}

/** marks obsolete constraint to be not obsolete anymore;
 *  it will be moved to the first part of the constraint arrays, such that it is checked, enforced, separated,
 *  and propagated before the obsolete constraints
 */
static
RETCODE conshdlrMarkConsUseful(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   CONS*            cons                /**< constraint to be marked obsolete */
   )
{
   CONS* tmpcons;
      
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(cons->obsolete);

   cons->obsolete = FALSE;

   if( cons->active )
   {
      if( cons->check )
      {
         assert(conshdlr->nusefulcheckconss <= cons->checkconsspos && cons->checkconsspos < conshdlr->ncheckconss);
         
         /* switch the first obsolete check constraint with this constraint */
         tmpcons = conshdlr->checkconss[conshdlr->nusefulcheckconss];
         assert(tmpcons->checkconsspos == conshdlr->nusefulcheckconss);
         
         conshdlr->checkconss[conshdlr->nusefulcheckconss] = cons;
         conshdlr->checkconss[cons->checkconsspos] = tmpcons;
         tmpcons->checkconsspos = cons->checkconsspos;
         cons->checkconsspos = conshdlr->nusefulcheckconss;
         
         conshdlr->nusefulcheckconss++;
      }
   }
   if( cons->enabled )
   {
      if( cons->separate )
      {
         assert(conshdlr->nusefulsepaconss <= cons->sepaconsspos && cons->sepaconsspos < conshdlr->nsepaconss);
         
         /* switch the first obsolete sepa constraint with this constraint */
         tmpcons = conshdlr->sepaconss[conshdlr->nusefulsepaconss];
         assert(tmpcons->sepaconsspos == conshdlr->nusefulsepaconss);
         
         conshdlr->sepaconss[conshdlr->nusefulsepaconss] = cons;
         conshdlr->sepaconss[cons->sepaconsspos] = tmpcons;
         tmpcons->sepaconsspos = cons->sepaconsspos;
         cons->sepaconsspos = conshdlr->nusefulsepaconss;
         
         conshdlr->nusefulsepaconss++;
      }
      if( cons->enforce )
      {
         assert(conshdlr->nusefulenfoconss <= cons->enfoconsspos && cons->enfoconsspos < conshdlr->nenfoconss);
         
         /* switch the first obsolete enfo constraint with this constraint */
         tmpcons = conshdlr->enfoconss[conshdlr->nusefulenfoconss];
         assert(tmpcons->enfoconsspos == conshdlr->nusefulenfoconss);
         
         conshdlr->enfoconss[conshdlr->nusefulenfoconss] = cons;
         conshdlr->enfoconss[cons->enfoconsspos] = tmpcons;
         tmpcons->enfoconsspos = cons->enfoconsspos;
         cons->enfoconsspos = conshdlr->nusefulenfoconss;
         
         conshdlr->nusefulenfoconss++;
      }
      if( cons->propagate && cons->propenabled )
      {
         assert(conshdlr->nusefulpropconss <= cons->propconsspos && cons->propconsspos < conshdlr->npropconss);
         
         /* switch the first obsolete prop constraint with this constraint */
         tmpcons = conshdlr->propconss[conshdlr->nusefulpropconss];
         assert(tmpcons->propconsspos == conshdlr->nusefulpropconss);
         
         conshdlr->propconss[conshdlr->nusefulpropconss] = cons;
         conshdlr->propconss[cons->propconsspos] = tmpcons;
         tmpcons->propconsspos = cons->propconsspos;
         cons->propconsspos = conshdlr->nusefulpropconss;
         
         conshdlr->nusefulpropconss++;
      }
   }

   return SCIP_OKAY;
}

/** adds constraint to the sepaconss array of constraint handler */
static
RETCODE conshdlrAddSepacons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->separate);
   assert(cons->sepaconsspos == -1);

   CHECK_OKAY( conshdlrEnsureSepaconssMem(conshdlr, set, conshdlr->nsepaconss+1) );
   insertpos = conshdlr->nsepaconss;
   if( !cons->obsolete )
   {
      if( conshdlr->nusefulsepaconss < conshdlr->nsepaconss )
      {
         conshdlr->sepaconss[conshdlr->nsepaconss] = conshdlr->sepaconss[conshdlr->nusefulsepaconss];
         conshdlr->sepaconss[conshdlr->nsepaconss]->sepaconsspos = conshdlr->nsepaconss;
         insertpos = conshdlr->nusefulsepaconss;
      }
      conshdlr->nusefulsepaconss++;
   }
   conshdlr->sepaconss[insertpos] = cons;
   cons->sepaconsspos = insertpos;
   conshdlr->nsepaconss++;

   return SCIP_OKAY;
}

/** deletes constraint from the sepaconss array of constraint handler */
static
void conshdlrDelSepacons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->separate);
   assert(cons->sepaconsspos != -1);

   delpos = cons->sepaconsspos;
   if( !cons->obsolete )
   {
      assert(0 <= delpos && delpos < conshdlr->nusefulsepaconss);

      if( delpos < conshdlr->lastnusefulsepaconss )
         conshdlr->lastnusefulsepaconss--;

      conshdlr->sepaconss[delpos] = conshdlr->sepaconss[conshdlr->nusefulsepaconss-1];
      conshdlr->sepaconss[delpos]->sepaconsspos = delpos;
      delpos = conshdlr->nusefulsepaconss-1;
      conshdlr->nusefulsepaconss--;
      assert(conshdlr->nusefulsepaconss >= 0);
      assert(conshdlr->lastnusefulsepaconss >= 0);
   }
   assert(conshdlr->nusefulsepaconss <= delpos && delpos < conshdlr->nsepaconss);
   if( delpos < conshdlr->nsepaconss-1 )
   {
      conshdlr->sepaconss[delpos] = conshdlr->sepaconss[conshdlr->nsepaconss-1];
      conshdlr->sepaconss[delpos]->sepaconsspos = delpos;
   }
   conshdlr->nsepaconss--;
   cons->sepaconsspos = -1;
}

/** adds constraint to the enfoconss array of constraint handler */
static
RETCODE conshdlrAddEnfocons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->enforce);
   assert(cons->enfoconsspos == -1);

   CHECK_OKAY( conshdlrEnsureEnfoconssMem(conshdlr, set, conshdlr->nenfoconss+1) );
   insertpos = conshdlr->nenfoconss;
   if( !cons->obsolete )
   {
      if( conshdlr->nusefulenfoconss < conshdlr->nenfoconss )
      {
         conshdlr->enfoconss[conshdlr->nenfoconss] = conshdlr->enfoconss[conshdlr->nusefulenfoconss];
         conshdlr->enfoconss[conshdlr->nenfoconss]->enfoconsspos = conshdlr->nenfoconss;
         insertpos = conshdlr->nusefulenfoconss;
      }
      conshdlr->nusefulenfoconss++;
   }
   else
   {
      /* we have to make sure that even this obsolete constraint is enforced in the next enforcement call;
       * if the same LP or pseudo solution is enforced again, only the newly added useful constraints are
       * enforced; thus, we have to reset the enforcement counters and force all constraints to be 
       * enforced again; this is not needed for separation and propagation, because they are not vital for correctness
       */
      conshdlr->lastenfolplpcount = -1;
      conshdlr->lastenfolpdomchgcount = -1;
      conshdlr->lastenfopsdomchgcount = -1;
   }
   conshdlr->enfoconss[insertpos] = cons;
   cons->enfoconsspos = insertpos;
   conshdlr->nenfoconss++;

   return SCIP_OKAY;
}

/** deletes constraint from the enfoconss array of constraint handler */
static
void conshdlrDelEnfocons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->enforce);
   assert(cons->enfoconsspos != -1);

   delpos = cons->enfoconsspos;
   if( !cons->obsolete )
   {
      assert(0 <= delpos && delpos < conshdlr->nusefulenfoconss);

      if( delpos < conshdlr->lastnusefulenfoconss )
         conshdlr->lastnusefulenfoconss--;

      conshdlr->enfoconss[delpos] = conshdlr->enfoconss[conshdlr->nusefulenfoconss-1];
      conshdlr->enfoconss[delpos]->enfoconsspos = delpos;
      delpos = conshdlr->nusefulenfoconss-1;
      conshdlr->nusefulenfoconss--;

      /* if the constraint that moved to the free position was a newly added constraint and not enforced in the last
       * enforcement, we have to make sure it will be enforced in the next run;
       * this check is not performed for separation and propagation, because they are not vital for correctness
       */
      if( delpos >= conshdlr->lastnusefulenfoconss )
         conshdlr->lastnusefulenfoconss = cons->enfoconsspos;
      conshdlr->lastnusefulenfoconss = MAX(conshdlr->lastnusefulenfoconss, 0);
      assert(conshdlr->nusefulenfoconss >= 0);
      assert(conshdlr->lastnusefulenfoconss >= 0);
   }
   assert(conshdlr->nusefulenfoconss <= delpos && delpos < conshdlr->nenfoconss);
   if( delpos < conshdlr->nenfoconss-1 )
   {
      conshdlr->enfoconss[delpos] = conshdlr->enfoconss[conshdlr->nenfoconss-1];
      conshdlr->enfoconss[delpos]->enfoconsspos = delpos;
   }
   conshdlr->nenfoconss--;
   cons->enfoconsspos = -1;
}

/** adds constraint to the chckconss array of constraint handler */
static
RETCODE conshdlrAddCheckcons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->check);
   assert(cons->checkconsspos == -1);

   CHECK_OKAY( conshdlrEnsureCheckconssMem(conshdlr, set, conshdlr->ncheckconss+1) );
   insertpos = conshdlr->ncheckconss;
   if( !cons->obsolete )
   {
      if( conshdlr->nusefulcheckconss < conshdlr->ncheckconss )
      {
         assert(conshdlr->checkconss[conshdlr->nusefulcheckconss] != NULL);
         conshdlr->checkconss[conshdlr->ncheckconss] = conshdlr->checkconss[conshdlr->nusefulcheckconss];
         conshdlr->checkconss[conshdlr->ncheckconss]->checkconsspos = conshdlr->ncheckconss;
         insertpos = conshdlr->nusefulcheckconss;
      }
      conshdlr->nusefulcheckconss++;
   }
   assert(0 <= insertpos && insertpos <= conshdlr->ncheckconss);
   conshdlr->checkconss[insertpos] = cons;
   cons->checkconsspos = insertpos;
   conshdlr->ncheckconss++;

   return SCIP_OKAY;
}

/** deletes constraint from the chckconss array of constraint handler */
static
void conshdlrDelCheckcons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->check);
   assert(cons->checkconsspos != -1);

   delpos = cons->checkconsspos;
   if( !cons->obsolete )
   {
      assert(0 <= delpos && delpos < conshdlr->nusefulcheckconss);
      conshdlr->checkconss[delpos] = conshdlr->checkconss[conshdlr->nusefulcheckconss-1];
      conshdlr->checkconss[delpos]->checkconsspos = delpos;
      delpos = conshdlr->nusefulcheckconss-1;
      conshdlr->nusefulcheckconss--;
   }
   assert(conshdlr->nusefulcheckconss <= delpos && delpos < conshdlr->ncheckconss);
   if( delpos < conshdlr->ncheckconss-1 )
   {
      conshdlr->checkconss[delpos] = conshdlr->checkconss[conshdlr->ncheckconss-1];
      conshdlr->checkconss[delpos]->checkconsspos = delpos;
   }
   conshdlr->ncheckconss--;
   cons->checkconsspos = -1;
}

/** adds constraint to the propconss array of constraint handler */
static
RETCODE conshdlrAddPropcons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->enabled);
   assert(cons->propagate);
   assert(cons->propenabled);
   assert(cons->propconsspos == -1);

   /* add constraint to the propagation array */
   CHECK_OKAY( conshdlrEnsurePropconssMem(conshdlr, set, conshdlr->npropconss+1) );
   insertpos = conshdlr->npropconss;
   if( !cons->obsolete )
   {
      if( conshdlr->nusefulpropconss < conshdlr->npropconss )
      {
         conshdlr->propconss[conshdlr->npropconss] = conshdlr->propconss[conshdlr->nusefulpropconss];
         conshdlr->propconss[conshdlr->npropconss]->propconsspos = conshdlr->npropconss;
         insertpos = conshdlr->nusefulpropconss;
      }
      conshdlr->nusefulpropconss++;
   }
   conshdlr->propconss[insertpos] = cons;
   cons->propconsspos = insertpos;
   conshdlr->npropconss++;

   return SCIP_OKAY;
}

/** deletes constraint from the propconss array of constraint handler */
static
void conshdlrDelPropcons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->propagate);
   assert(cons->propenabled);
   assert(cons->propconsspos != -1);

   /* delete constraint from the propagation array */
   delpos = cons->propconsspos;
   if( !cons->obsolete )
   {
      assert(0 <= delpos && delpos < conshdlr->nusefulpropconss);

      if( delpos < conshdlr->lastnusefulpropconss )
         conshdlr->lastnusefulpropconss--;

      conshdlr->propconss[delpos] = conshdlr->propconss[conshdlr->nusefulpropconss-1];
      conshdlr->propconss[delpos]->propconsspos = delpos;
      delpos = conshdlr->nusefulpropconss-1;
      conshdlr->nusefulpropconss--;
      assert(conshdlr->nusefulpropconss >= 0);
      assert(conshdlr->lastnusefulpropconss >= 0);
   }
   assert(conshdlr->nusefulpropconss <= delpos && delpos < conshdlr->npropconss);
   if( delpos < conshdlr->npropconss-1 )
   {
      conshdlr->propconss[delpos] = conshdlr->propconss[conshdlr->npropconss-1];
      conshdlr->propconss[delpos]->propconsspos = delpos;
   }
   conshdlr->npropconss--;
   cons->propconsspos = -1;
}

/** enables propagation of constraint */
static
RETCODE conshdlrEnableConsPropagation(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->propenabled);
   assert(cons->propconsspos == -1);

   debugMessage("enable propagation of constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* enable propagation of constraint */
   cons->propenabled = TRUE;

   /* add constraint to the propagation array */
   if( cons->enabled && cons->propagate )
   {
      CHECK_OKAY( conshdlrAddPropcons(conshdlr, set, cons) );
   }

   return SCIP_OKAY;
}

/** disables propagation of constraint */
static
RETCODE conshdlrDisableConsPropagation(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to remove */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->propenabled);
   assert((cons->propagate && cons->enabled) == (cons->propconsspos != -1));

   debugMessage("disable propagation of constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* delete constraint from the propagation array */
   if( cons->propagate && cons->enabled )
   {
      conshdlrDelPropcons(conshdlr, set, cons);
   }
   assert(cons->propconsspos == -1);

   /* disable propagation of constraint */
   cons->propenabled = FALSE;

   return SCIP_OKAY;
}

/** enables separation, enforcement, and propagation of constraint */
static
RETCODE conshdlrEnableCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(stat != NULL);
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
   stat->nenabledconss++;

   /* add constraint to the separation array */
   if( cons->separate )
   {
      CHECK_OKAY( conshdlrAddSepacons(conshdlr, set, cons) );
   }
      
   /* add constraint to the enforcement array */
   if( cons->enforce )
   {
      CHECK_OKAY( conshdlrAddEnfocons(conshdlr, set, cons) );
   }

   /* add constraint to the propagation array */
   if( cons->propagate && cons->propenabled )
   {
      CHECK_OKAY( conshdlrAddPropcons(conshdlr, set, cons) );
   }

   /* call constraint handler's enabling notification method */
   if( conshdlr->consenable != NULL )
   {
      CHECK_OKAY( conshdlr->consenable(set->scip, conshdlr, cons) );
   }

   return SCIP_OKAY;
}

/** disables separation, enforcement, and propagation of constraint */
static
RETCODE conshdlrDisableCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   CONS*            cons                /**< constraint to remove */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(stat != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->enabled);
   assert(cons->separate == (cons->sepaconsspos != -1));
   assert(cons->enforce == (cons->enfoconsspos != -1));
   assert((cons->propagate && cons->propenabled) == (cons->propconsspos != -1));

   debugMessage("disable constraint <%s> at sepa position %d in constraint handler <%s> (%d/%d)\n", 
      cons->name, cons->sepaconsspos, conshdlr->name, conshdlr->nusefulsepaconss, conshdlr->nsepaconss);

   /* call constraint handler's disabling notification method */
   if( conshdlr->consdisable != NULL )
   {
      CHECK_OKAY( conshdlr->consdisable(set->scip, conshdlr, cons) );
   }

   /* delete constraint from the separation array */
   if( cons->separate )
   {
      conshdlrDelSepacons(conshdlr, set, cons);
   }

   /* delete constraint from the enforcement array */
   if( cons->enforce )
   {
      conshdlrDelEnfocons(conshdlr, set, cons);
   }

   /* delete constraint from the propagation array */
   if( cons->propagate && cons->propenabled )
   {
      conshdlrDelPropcons(conshdlr, set, cons);
   }

   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->propconsspos == -1);

   /* disable constraint */
   cons->enabled = FALSE;
   conshdlr->nenabledconss--;
   stat->nenabledconss--;

   return SCIP_OKAY;
}

/** activates and adds constraint to constraint handler's constraint arrays */
static
RETCODE conshdlrActivateCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   CONS*            cons,               /**< constraint to add */
   int              depth               /**< depth in the tree where the activation takes place, or -1 for global problem */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(stat != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->active);
   assert(!cons->enabled);
   assert(cons->consspos == -1);
   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->checkconsspos == -1);
   assert(cons->propconsspos == -1);
   assert(depth >= -1);
   assert(depth == -1 || depth <= SCIPgetDepth(set->scip));

   debugMessage("activate constraint <%s> in constraint handler <%s> (depth %d)\n", cons->name, conshdlr->name, depth);

   /* activate constraint */
   CHECK_OKAY( conshdlrEnsureConssMem(conshdlr, set, conshdlr->nconss+1) );
   cons->active = TRUE;
   cons->activedepth = depth;
   cons->consspos = conshdlr->nconss;
   conshdlr->conss[conshdlr->nconss] = cons;
   conshdlr->nconss++;
   conshdlr->maxnconss = MAX(conshdlr->maxnconss, conshdlr->nconss);
   stat->nactiveconss++;

   /* add constraint to the check array */
   if( cons->check )
   {
      CHECK_OKAY( conshdlrAddCheckcons(conshdlr, set, cons) );
   }

   /* call constraint handler's activation notification method */
   if( conshdlr->consactive != NULL )
   {
      CHECK_OKAY( conshdlr->consactive(set->scip, conshdlr, cons) );
   }

   /* enable separation, enforcement, and propagation of constraint */
   CHECK_OKAY( conshdlrEnableCons(conshdlr, set, stat, cons) );

   return SCIP_OKAY;
}

/** deactivates and removes constraint from constraint handler's conss array */
static
RETCODE conshdlrDeactivateCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(stat != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->consspos != -1);
   assert(cons->check == (cons->checkconsspos != -1));

   debugMessage("deactivate constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* disable constraint */
   if( cons->enabled )
   {
      CHECK_OKAY( conshdlrDisableCons(conshdlr, set, stat, cons) );
   }
   assert(!cons->enabled);

   /* call constraint handler's deactivation notification method */
   if( conshdlr->consdeactive != NULL )
   {
      CHECK_OKAY( conshdlr->consdeactive(set->scip, conshdlr, cons) );
   }

   /* delete constraint from the check array */
   if( cons->check )
   {
      conshdlrDelCheckcons(conshdlr, set, cons);
   }

   /* delete constraint from the conss array */
   delpos = cons->consspos;
   assert(0 <= delpos && delpos < conshdlr->nconss);
   if( delpos < conshdlr->nconss-1 )
   {
      conshdlr->conss[delpos] = conshdlr->conss[conshdlr->nconss-1];
      conshdlr->conss[delpos]->consspos = delpos;
   }
   conshdlr->nconss--;
   cons->consspos = -1;
   cons->active = FALSE;
   cons->activedepth = -2;
   stat->nactiveconss--;

   assert(cons->consspos == -1);
   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->checkconsspos == -1);
   assert(cons->propconsspos == -1);

   return SCIP_OKAY;
}

/** processes all delayed updates of constraints:
 *  recently (de)activated constraints will be (de)activated;
 *  recently en/disabled constraints will be en/disabled;
 *  recent obsolete non-check constraints will be globally deleted;
 *  recent obsolete check constraints will be moved to the last positions in the sepa-, enfo-, check-, and prop-arrays;
 *  recent useful constraints will be moved to the first positions in the sepa-, enfo-, check-, and prop-arrays;
 */
static
RETCODE conshdlrProcessUpdates(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob                /**< problem data */
   )
{
   CONS* cons;
   int i;

   assert(conshdlr != NULL);
   assert(!conshdlr->delayupdates);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);

   debugMessage("processing %d constraints that have to be updated in constraint handler <%s>\n",
      conshdlr->nupdateconss, conshdlr->name);

   for( i = 0; i < conshdlr->nupdateconss; ++i )
   {
      cons = conshdlr->updateconss[i];
      assert(cons != NULL);
      assert(cons->conshdlr == conshdlr);
      assert(cons->update);
      assert(cons->updateactivate || cons->updatedeactivate || cons->updateenable || cons->updatedisable
         || cons->updatepropenable || cons->updatepropdisable || cons->updatedelete || cons->updateobsolete);

      debugMessage(" -> constraint <%s>: activate=%d, deactivate=%d, enable=%d, disable=%d, propenable=%d, propdisable=%d, delete=%d, obsolete=%d (consdata=%p)\n",
         cons->name, cons->updateactivate, cons->updatedeactivate, cons->updateenable, cons->updatedisable,
         cons->updatepropenable, cons->updatepropdisable, cons->updatedelete, cons->updateobsolete, cons->consdata);

      if( cons->updateactivate )
      {
         assert(!cons->active);
         assert(!cons->updatedeactivate);
         assert(!cons->updateenable);
         assert(!cons->updatedisable);
         assert(!cons->updatedelete);
         assert(!cons->updateobsolete);

         /* the activation depth was already stored in SCIPconsActivate() */
         CHECK_OKAY( conshdlrActivateCons(conshdlr, set, stat, cons, cons->activedepth) );
         assert(cons->active);
         cons->updateactivate = FALSE;
      }
      else if( cons->updatedeactivate )
      {
         assert(cons->active);

         CHECK_OKAY( conshdlrDeactivateCons(conshdlr, set, stat, cons) );
         assert(!cons->active);
         cons->updatedeactivate = FALSE;
         cons->updateenable = FALSE;
         cons->updatedisable = FALSE;
         cons->obsolete = consExceedsObsoleteage(cons, set);
         cons->updateobsolete = FALSE;
      }
      else if( cons->updateenable )
      {
         assert(!cons->enabled);
         assert(!cons->updatedisable);

         CHECK_OKAY( conshdlrEnableCons(conshdlr, set, stat, cons) );
         assert(cons->enabled);
         cons->updateenable = FALSE;
      }
      else if( cons->updatedisable )
      {
         assert(cons->enabled);

         CHECK_OKAY( conshdlrDisableCons(conshdlr, set, stat, cons) );
         assert(!cons->enabled);
         cons->updatedisable = FALSE;
      }

      if( cons->updatepropenable )
      {
         assert(!cons->updatepropdisable);
         if( !cons->propenabled )
         {
            CHECK_OKAY( conshdlrEnableConsPropagation(conshdlr, set, cons) );
            assert(cons->propenabled);
         }
         cons->updatepropenable = FALSE;
      }
      else if( cons->updatepropdisable )
      {
         if( cons->propenabled )
         {         
            CHECK_OKAY( conshdlrDisableConsPropagation(conshdlr, set, cons) );
            assert(!cons->propenabled);
         }
         cons->updatepropdisable = FALSE;
      }

      if( cons->updatedelete )
      {
         assert(!cons->check);
         CHECK_OKAY( SCIPconsDelete(cons, memhdr, set, stat, prob) );
         cons->updatedelete = FALSE;
         cons->updateobsolete = FALSE;
      }
      else if( cons->updateobsolete )
      {
         if( !cons->obsolete && consExceedsObsoleteage(cons, set) )
         {
            /* the constraint's status must be switched to obsolete */
            CHECK_OKAY( conshdlrMarkConsObsolete(conshdlr, cons) );
         }
         else if( cons->obsolete && !consExceedsObsoleteage(cons, set) )
         {
            /* the constraint's status must be switched to useful */
            CHECK_OKAY( conshdlrMarkConsUseful(conshdlr, cons) );
         }
         cons->updateobsolete = FALSE;
      }
      assert(!cons->updateactivate);
      assert(!cons->updatedeactivate);
      assert(!cons->updateenable);
      assert(!cons->updatedisable);
      assert(!cons->updatepropenable);
      assert(!cons->updatepropdisable);
      assert(!cons->updatedelete);
      assert(!cons->updateobsolete);
      cons->update = FALSE;

      /* release the constraint */
      CHECK_OKAY( SCIPconsRelease(&conshdlr->updateconss[i], memhdr, set) );
   }

   conshdlr->nupdateconss = 0;

   return SCIP_OKAY;
}

/** marks constraint handler to delay all constraint updates until the next SCIPconshdlrProcessUpdates() call */
static
void conshdlrDelayUpdates(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);
   assert(!conshdlr->delayupdates);
   
   debugMessage("constraint updates of constraint handler <%s> will be delayed\n", conshdlr->name);

   conshdlr->delayupdates = TRUE;
}

/** marks constraint handler to perform all constraint updates immediately;
 *  all delayed constraint updates will be processed
 */
static
RETCODE conshdlrForceUpdates(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob                /**< problem data */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->delayupdates);
   
   debugMessage("constraint updates of constraint handler <%s> will be processed immediately\n", conshdlr->name);
   conshdlr->delayupdates = FALSE;

   CHECK_OKAY( conshdlrProcessUpdates(conshdlr, memhdr, set, stat, prob) );
   assert(conshdlr->nupdateconss == 0);

   return SCIP_OKAY;
}

/** adds constraint to constraint handler's update constraint array and captures it */
static
RETCODE conshdlrAddUpdateCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);

   if( !cons->update )
   {
      debugMessage("constraint <%s> of age %g has to be updated in constraint handler <%s> (consdata=%p)\n",
         cons->name, cons->age, conshdlr->name, cons->consdata);
      
      /* add constraint to the updateconss array */
      CHECK_OKAY( conshdlrEnsureUpdateconssMem(conshdlr, set, conshdlr->nupdateconss+1) );
      conshdlr->updateconss[conshdlr->nupdateconss] = cons;
      conshdlr->nupdateconss++;
      
      /* capture constraint */
      SCIPconsCapture(cons);
      
      cons->update = TRUE;
   }

   return SCIP_OKAY;
}

/** compares two constraint handlers w. r. to their separation priority */
DECL_SORTPTRCOMP(SCIPconshdlrCompSepa)
{  /*lint --e{715}*/
   return ((CONSHDLR*)elem2)->sepapriority - ((CONSHDLR*)elem1)->sepapriority;
}

/** compares two constraint handlers w. r. to their enforcing priority */
DECL_SORTPTRCOMP(SCIPconshdlrCompEnfo)
{  /*lint --e{715}*/
   return ((CONSHDLR*)elem2)->enfopriority - ((CONSHDLR*)elem1)->enfopriority;
}

/** compares two constraint handlers w. r. to their feasibility check priority */
DECL_SORTPTRCOMP(SCIPconshdlrCompCheck)
{  /*lint --e{715}*/
   return ((CONSHDLR*)elem2)->checkpriority - ((CONSHDLR*)elem1)->checkpriority;
}

/** creates a constraint handler */
RETCODE SCIPconshdlrCreate(
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              checkpriority,      /**< priority of the constraint handler for checking infeasibility */
   int              sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int              propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int              eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   int              maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   DECL_CONSINIT    ((*consinit)),      /**< initialize constraint handler */
   DECL_CONSEXIT    ((*consexit)),      /**< deinitialize constraint handler */
   DECL_CONSINITPRE ((*consinitpre)),   /**< presolving initialization method of constraint handler */
   DECL_CONSEXITPRE ((*consexitpre)),   /**< presolving deinitialization method of constraint handler */
   DECL_CONSINITSOL ((*consinitsol)),   /**< solving process initialization method of constraint handler */
   DECL_CONSEXITSOL ((*consexitsol)),   /**< solving process deinitialization method of constraint handler */
   DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSINITLP  ((*consinitlp)),    /**< initialize LP with relaxations of "initial" constraints */
   DECL_CONSSEPA    ((*conssepa)),      /**< separate cutting planes */
   DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method */
   DECL_CONSRESPROP ((*consresprop)),   /**< propagation conflict resolving method */
   DECL_CONSLOCK    ((*conslock)),      /**< variable rounding lock method */
   DECL_CONSACTIVE  ((*consactive)),    /**< activation notification method */
   DECL_CONSDEACTIVE((*consdeactive)),  /**< deactivation notification method */
   DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   DECL_CONSPRINT   ((*consprint)),     /**< constraint display method */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   char paramname[MAXSTRLEN];

   assert(conshdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert((conssepa != NULL) || (sepafreq == -1));
   assert((consprop != NULL) || (propfreq == -1));
   assert(eagerfreq >= -1);

   ALLOC_OKAY( allocMemory(conshdlr) );
   ALLOC_OKAY( duplicateMemoryArray(&(*conshdlr)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*conshdlr)->desc, desc, strlen(desc)+1) );
   (*conshdlr)->sepapriority = sepapriority;
   (*conshdlr)->enfopriority = enfopriority;
   (*conshdlr)->checkpriority = checkpriority;
   (*conshdlr)->sepafreq = sepafreq;
   (*conshdlr)->propfreq = propfreq;
   (*conshdlr)->eagerfreq = eagerfreq;
   (*conshdlr)->maxprerounds = maxprerounds;
   (*conshdlr)->consfree = consfree;
   (*conshdlr)->consinit = consinit;
   (*conshdlr)->consexit = consexit;
   (*conshdlr)->consinitpre = consinitpre;
   (*conshdlr)->consexitpre = consexitpre;
   (*conshdlr)->consinitsol = consinitsol;
   (*conshdlr)->consexitsol = consexitsol;
   (*conshdlr)->consdelete = consdelete;
   (*conshdlr)->constrans = constrans;
   (*conshdlr)->consinitlp = consinitlp;
   (*conshdlr)->conssepa = conssepa;
   (*conshdlr)->consenfolp = consenfolp;
   (*conshdlr)->consenfops = consenfops;
   (*conshdlr)->conscheck = conscheck;
   (*conshdlr)->consprop = consprop;
   (*conshdlr)->conspresol = conspresol;
   (*conshdlr)->consresprop = consresprop;
   (*conshdlr)->conslock = conslock;
   (*conshdlr)->consactive = consactive;
   (*conshdlr)->consdeactive = consdeactive;
   (*conshdlr)->consenable = consenable;
   (*conshdlr)->consdisable = consdisable;
   (*conshdlr)->consprint = consprint;
   (*conshdlr)->conshdlrdata = conshdlrdata;
   (*conshdlr)->conss = NULL;
   (*conshdlr)->consssize = 0;
   (*conshdlr)->nconss = 0;
   (*conshdlr)->maxnconss = 0;
   (*conshdlr)->startnconss = 0;
   (*conshdlr)->sepaconss = NULL;
   (*conshdlr)->sepaconsssize = 0;
   (*conshdlr)->nsepaconss = 0;
   (*conshdlr)->nusefulsepaconss = 0;
   (*conshdlr)->enfoconss = NULL;
   (*conshdlr)->enfoconsssize = 0;
   (*conshdlr)->nenfoconss = 0;
   (*conshdlr)->nusefulenfoconss = 0;
   (*conshdlr)->checkconss = NULL;
   (*conshdlr)->checkconsssize = 0;
   (*conshdlr)->ncheckconss = 0;
   (*conshdlr)->nusefulcheckconss = 0;
   (*conshdlr)->propconss = NULL;
   (*conshdlr)->propconsssize = 0;
   (*conshdlr)->npropconss = 0;
   (*conshdlr)->nusefulpropconss = 0;
   (*conshdlr)->updateconss = NULL;
   (*conshdlr)->updateconsssize = 0;
   (*conshdlr)->nupdateconss = 0;
   (*conshdlr)->nenabledconss = 0;
   (*conshdlr)->lastnusefulpropconss = 0;
   (*conshdlr)->lastnusefulsepaconss = 0;
   (*conshdlr)->lastnusefulenfoconss = 0;

   CHECK_OKAY( SCIPclockCreate(&(*conshdlr)->presoltime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*conshdlr)->sepatime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*conshdlr)->enfolptime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*conshdlr)->enfopstime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*conshdlr)->proptime, SCIP_CLOCKTYPE_DEFAULT) );

   (*conshdlr)->nsepacalls = 0;
   (*conshdlr)->nenfolpcalls = 0;
   (*conshdlr)->nenfopscalls = 0;
   (*conshdlr)->npropcalls = 0;
   (*conshdlr)->ncutoffs = 0;
   (*conshdlr)->ncutsfound = 0;
   (*conshdlr)->nconssfound = 0;
   (*conshdlr)->ndomredsfound = 0;
   (*conshdlr)->nchildren = 0;
   (*conshdlr)->lastpropdomchgcount = -1;
   (*conshdlr)->lastsepalpcount = -1;
   (*conshdlr)->lastenfolplpcount = -1;
   (*conshdlr)->lastenfolpdomchgcount = -1;
   (*conshdlr)->lastenfopsdomchgcount = -1;
   (*conshdlr)->lastnfixedvars = 0;
   (*conshdlr)->lastnaggrvars = 0;
   (*conshdlr)->lastnchgvartypes = 0;
   (*conshdlr)->lastnchgbds = 0;
   (*conshdlr)->lastnaddholes = 0;
   (*conshdlr)->lastndelconss = 0;
   (*conshdlr)->lastnupgdconss = 0;
   (*conshdlr)->lastnchgcoefs = 0;
   (*conshdlr)->lastnchgsides = 0;
   (*conshdlr)->nfixedvars = 0;
   (*conshdlr)->naggrvars = 0;
   (*conshdlr)->nchgvartypes = 0;
   (*conshdlr)->nchgbds = 0;
   (*conshdlr)->naddholes = 0;
   (*conshdlr)->ndelconss = 0;
   (*conshdlr)->nupgdconss = 0;
   (*conshdlr)->nchgcoefs = 0;
   (*conshdlr)->nchgsides = 0;
   (*conshdlr)->needscons = needscons;
   (*conshdlr)->initialized = FALSE;
   (*conshdlr)->delayupdates = FALSE;

   /* add parameters */
   sprintf(paramname, "constraints/%s/sepafreq", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, 
         "frequency for separating cuts (-1: never, 0: only in root node)",
         &(*conshdlr)->sepafreq, sepafreq, -1, INT_MAX, NULL, NULL) );

   sprintf(paramname, "constraints/%s/propfreq", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, 
         "frequency for propagating domains (-1: never, 0: only in root node)",
         &(*conshdlr)->propfreq, propfreq, -1, INT_MAX, NULL, NULL) );

   sprintf(paramname, "constraints/%s/eagerfreq", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, 
         "frequency for using all instead of only the useful constraints in separation, propagation and enforcement (-1: never, 0: only in first evaluation)",
         &(*conshdlr)->eagerfreq, eagerfreq, -1, INT_MAX, NULL, NULL) );

   sprintf(paramname, "constraints/%s/maxprerounds", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, 
         "maximal number of presolving rounds the constraint handler participates in (-1: no limit)",
         &(*conshdlr)->maxprerounds, maxprerounds, -1, INT_MAX, NULL, NULL) );

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
   assert(SCIPgetStage(scip) == SCIP_STAGE_INIT);

   /* call destructor of constraint handler */
   if( (*conshdlr)->consfree != NULL )
   {
      CHECK_OKAY( (*conshdlr)->consfree(scip, *conshdlr) );
   }

   SCIPclockFree(&(*conshdlr)->presoltime);
   SCIPclockFree(&(*conshdlr)->sepatime);
   SCIPclockFree(&(*conshdlr)->enfolptime);
   SCIPclockFree(&(*conshdlr)->enfopstime);
   SCIPclockFree(&(*conshdlr)->proptime);

   freeMemoryArray(&(*conshdlr)->name);
   freeMemoryArray(&(*conshdlr)->desc);
   freeMemoryArrayNull(&(*conshdlr)->conss);
   freeMemoryArrayNull(&(*conshdlr)->sepaconss);
   freeMemoryArrayNull(&(*conshdlr)->enfoconss);
   freeMemoryArrayNull(&(*conshdlr)->checkconss);
   freeMemoryArrayNull(&(*conshdlr)->propconss);
   freeMemoryArrayNull(&(*conshdlr)->updateconss);
   freeMemory(conshdlr);

   return SCIP_OKAY;
}

/** calls init method of constraint handler */
RETCODE SCIPconshdlrInit(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conshdlr != NULL);

   if( conshdlr->initialized )
   {
      errorMessage("constraint handler <%s> already initialized\n", conshdlr->name);
      return SCIP_INVALIDCALL;
   }

   SCIPclockReset(conshdlr->presoltime);
   SCIPclockReset(conshdlr->sepatime);
   SCIPclockReset(conshdlr->enfolptime);
   SCIPclockReset(conshdlr->enfopstime);
   SCIPclockReset(conshdlr->proptime);

   conshdlr->nsepacalls = 0;
   conshdlr->nenfolpcalls = 0;
   conshdlr->nenfopscalls = 0;
   conshdlr->npropcalls = 0;
   conshdlr->ncutoffs = 0;
   conshdlr->ncutsfound = 0;
   conshdlr->nconssfound = 0;
   conshdlr->ndomredsfound = 0;
   conshdlr->nchildren = 0;
   conshdlr->lastpropdomchgcount = -1;
   conshdlr->lastenfolpdomchgcount = -1;
   conshdlr->lastenfopsdomchgcount = -1;
   conshdlr->maxnconss = conshdlr->nconss;
   conshdlr->startnconss = 0;
   conshdlr->lastsepalpcount = -1;
   conshdlr->lastenfolplpcount = -1;
   conshdlr->lastnusefulpropconss = 0;
   conshdlr->lastnusefulsepaconss = 0;
   conshdlr->lastnusefulenfoconss = 0;
   conshdlr->lastnfixedvars = 0;
   conshdlr->lastnaggrvars = 0;
   conshdlr->lastnchgvartypes = 0;
   conshdlr->lastnchgbds = 0;
   conshdlr->lastnaddholes = 0;
   conshdlr->lastndelconss = 0;
   conshdlr->lastnupgdconss = 0;
   conshdlr->lastnchgcoefs = 0;
   conshdlr->lastnchgsides = 0;
   conshdlr->nfixedvars = 0;
   conshdlr->naggrvars = 0;
   conshdlr->nchgvartypes = 0;
   conshdlr->nchgbds = 0;
   conshdlr->naddholes = 0;
   conshdlr->ndelconss = 0;
   conshdlr->nupgdconss = 0;
   conshdlr->nchgcoefs = 0;
   conshdlr->nchgsides = 0;

   /* call initialization method of constraint handler */
   if( conshdlr->consinit != NULL )
   {
      CHECK_OKAY( conshdlr->consinit(scip, conshdlr, conshdlr->conss, conshdlr->nconss) );
   }
   conshdlr->initialized = TRUE;
   assert(!conshdlr->delayupdates);

   return SCIP_OKAY;
}

/** calls exit method of constraint handler */
RETCODE SCIPconshdlrExit(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conshdlr != NULL);

   if( !conshdlr->initialized )
   {
      errorMessage("constraint handler <%s> not initialized\n", conshdlr->name);
      return SCIP_INVALIDCALL;
   }

   /* call deinitialization method of constraint handler */
   if( conshdlr->consexit != NULL )
   {
      CHECK_OKAY( conshdlr->consexit(scip, conshdlr, conshdlr->conss, conshdlr->nconss) );
   }
   conshdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs constraint handler that the presolving process is being started */
RETCODE SCIPconshdlrInitpre(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* call presolving initialization method of constraint handler */
   if( conshdlr->consinitpre != NULL )
   {
      CHECK_OKAY( conshdlr->consinitpre(scip, conshdlr, conshdlr->conss, conshdlr->nconss, result) );

      /* evaluate result */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_UNBOUNDED
         && *result != SCIP_FEASIBLE )
      {
         errorMessage("presolving initialization method of constraint handler <%s> returned invalid result <%d>\n", 
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** informs constraint handler that the presolving is finished */
RETCODE SCIPconshdlrExitpre(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* call presolving deinitialization method of constraint handler */
   if( conshdlr->consexitpre != NULL )
   {
      CHECK_OKAY( conshdlr->consexitpre(scip, conshdlr, conshdlr->conss, conshdlr->nconss, result) );

      /* evaluate result */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_UNBOUNDED
         && *result != SCIP_FEASIBLE )
      {
         errorMessage("presolving deinitialization method of constraint handler <%s> returned invalid result <%d>\n", 
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   /* update statistics */
   conshdlr->maxnconss = conshdlr->nconss;
   conshdlr->startnconss = conshdlr->nconss;

   return SCIP_OKAY;
}

/** informs constraint handler that the branch and bound process is being started */
RETCODE SCIPconshdlrInitsol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conshdlr != NULL);

   /* call solving process initialization method of constraint handler */
   if( conshdlr->consinitsol != NULL )
   {
      CHECK_OKAY( conshdlr->consinitsol(scip, conshdlr, conshdlr->conss, conshdlr->nconss) );
   }

   return SCIP_OKAY;
}

/** informs constraint handler that the branch and bound process data is being freed */
RETCODE SCIPconshdlrExitsol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conshdlr != NULL);

   /* call solving process deinitialization method of constraint handler */
   if( conshdlr->consexitsol != NULL )
   {
      CHECK_OKAY( conshdlr->consexitsol(scip, conshdlr, conshdlr->conss, conshdlr->nconss) );
   }

   return SCIP_OKAY;
}

/** calls LP initialization method of constraint handler to separate all initial constraints */
RETCODE SCIPconshdlrInitLP(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob                /**< problem data */
   )
{
   assert(conshdlr != NULL);

   if( conshdlr->consinitlp != NULL )
   {
      debugMessage("initializing LP with %d constraints of handler <%s>\n", conshdlr->nconss, conshdlr->name);
         
      /* because during constraint processing, constraints of this handler may be activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);
      
      /* call external method */
      CHECK_OKAY( conshdlr->consinitlp(set->scip, conshdlr, conshdlr->conss, conshdlr->nconss) );

      /* perform the cached constraint updates */
      CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, stat, prob) );
   }

   return SCIP_OKAY;
}

/** calls separator method of constraint handler to separate all constraints added after last conshdlrResetSepa() call */
RETCODE SCIPconshdlrSeparate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< problem data */
   SEPASTORE*       sepastore,          /**< separation storage */
   int              depth,              /**< depth of current node */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(conshdlr->lastsepalpcount != stat->lpcount
      || (0 <= conshdlr->lastnusefulsepaconss && conshdlr->lastnusefulsepaconss <= conshdlr->nusefulsepaconss));
   assert(set != NULL);
   assert(stat != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->conssepa != NULL
      && ((depth == 0 && conshdlr->sepafreq == 0) || (conshdlr->sepafreq > 0 && depth % conshdlr->sepafreq == 0)) )
   {
      int nconss;
      int nusefulconss;
      int firstcons;

      /* check, if this LP solution was already separated */
      if( conshdlr->lastsepalpcount == stat->lpcount )
      {
         /* all constraints that were not yet separated on the new LP solution must be useful constraints, which means,
          * that the new constraints are the last constraints of the useful ones
          */
         nconss = conshdlr->nusefulsepaconss - conshdlr->lastnusefulsepaconss;
         nusefulconss = nconss;
         firstcons = conshdlr->lastnusefulsepaconss;
      }
      else
      {
         /* on a new LP solution, we want to separate all constraints */
         nconss = conshdlr->nsepaconss;
         nusefulconss = conshdlr->nusefulsepaconss;
         firstcons = 0;
      }
      assert(firstcons >= 0);
      assert(firstcons + nconss <= conshdlr->nsepaconss);
      assert(nusefulconss <= nconss);

      /* constraint handlers without constraints should only be called once */
      if( nconss > 0 || (!conshdlr->needscons && conshdlr->lastsepalpcount != stat->lpcount) )
      {
         CONS** conss;
         Longint oldndomchgs;
         int oldncutsstored;
         int oldnactiveconss;

         debugMessage("separating constraints %d to %d of %d constraints of handler <%s> (%s LP solution)\n",
            firstcons, firstcons + nconss - 1, conshdlr->nsepaconss, conshdlr->name,
            conshdlr->lastsepalpcount == stat->lpcount ? "old" : "new");

         /* remember the number of processed constraints on the current LP solution */
         conshdlr->lastsepalpcount = stat->lpcount;
         conshdlr->lastnusefulsepaconss = conshdlr->nusefulsepaconss;

         /* get the array of the constraints to be processed */
         conss = &(conshdlr->sepaconss[firstcons]);
         
         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldncutsstored = SCIPsepastoreGetNCutsStored(sepastore);
         oldnactiveconss = stat->nactiveconss;

         /* check, if we want to use eager evaluation */
         if( (conshdlr->eagerfreq == 0 && conshdlr->nsepacalls == 0)
            || (conshdlr->eagerfreq > 0 && conshdlr->nsepacalls % conshdlr->eagerfreq == 0) )
            nusefulconss = nconss;

         /* because during constraint processing, constraints of this handler may be activated, deactivated,
          * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
          * external method; to avoid this, these changes will be buffered and processed after the method call
          */
         conshdlrDelayUpdates(conshdlr);

         /* start timing */
         SCIPclockStart(conshdlr->sepatime, set);

         /* call external method */
         CHECK_OKAY( conshdlr->conssepa(set->scip, conshdlr, conss, nconss, nusefulconss, result) );
         debugMessage(" -> separating returned result <%d>\n", *result);

         /* stop timing */
         SCIPclockStop(conshdlr->sepatime, set);

         /* perform the cached constraint updates */
         CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, stat, prob) );

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN )
         {
            errorMessage("separation method of constraint handler <%s> returned invalid result <%d>\n", 
               conshdlr->name, *result);
            return SCIP_INVALIDRESULT;
         }

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            conshdlr->nsepacalls++;
         if( *result == SCIP_CUTOFF )
            conshdlr->ncutoffs++;
         conshdlr->ncutsfound += SCIPsepastoreGetNCutsStored(sepastore) - oldncutsstored; /*lint !e776*/
         conshdlr->nconssfound += MAX(stat->nactiveconss - oldnactiveconss, 0);
         conshdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
      }
   }

   return SCIP_OKAY;
}

/** calls enforcing method of constraint handler for LP solution for all constraints added after last
 *  conshdlrResetEnfo() call
 */
RETCODE SCIPconshdlrEnforceLPSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   SEPASTORE*       sepastore,          /**< separation storage */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(conshdlr->lastenfolplpcount != stat->lpcount || conshdlr->lastenfolpdomchgcount != stat->domchgcount
      || (0 <= conshdlr->lastnusefulenfoconss && conshdlr->lastnusefulenfoconss <= conshdlr->nusefulenfoconss));
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->consenfolp != NULL )
   {
      int nconss;
      int nusefulconss;
      int firstcons;
      Bool lpchanged;

      /* check, if this LP solution was already enforced */
      if( conshdlr->lastenfolplpcount == stat->lpcount && conshdlr->lastenfolpdomchgcount == stat->domchgcount )
      {
         /* all constraints that were not yet enforced on the new LP solution must be useful constraints, which means,
          * that the new constraints are the last constraints of the useful ones
          */
         nconss = conshdlr->nusefulenfoconss - conshdlr->lastnusefulenfoconss;
         nusefulconss = nconss;
         firstcons = conshdlr->lastnusefulenfoconss;
         lpchanged = FALSE;
      }
      else
      {
         /* on a new LP solution, we want to enforce all constraints */
         nconss = conshdlr->nenfoconss;
         nusefulconss = conshdlr->nusefulenfoconss;
         firstcons = 0;
         lpchanged = TRUE;
      }
      assert(firstcons >= 0);
      assert(firstcons + nconss <= conshdlr->nenfoconss);
      assert(nusefulconss <= nconss);

      /* constraint handlers without constraints should only be called once */
      if( nconss > 0 || (!conshdlr->needscons && lpchanged) )
      {
         CONS** conss;
         Longint oldndomchgs;
         int oldncutsstored;
         int oldnactiveconss;

         debugMessage("enforcing constraints %d to %d of %d constraints of handler <%s> (%s LP solution)\n",
            firstcons, firstcons + nconss - 1, conshdlr->nenfoconss, conshdlr->name, lpchanged ? "new" : "old");

         /* remember the number of processed constraints on the current LP solution */
         conshdlr->lastenfolplpcount = stat->lpcount;
         conshdlr->lastenfolpdomchgcount = stat->domchgcount;
         conshdlr->lastnusefulenfoconss = conshdlr->nusefulenfoconss;

         /* get the array of the constraints to be processed */
         conss = &(conshdlr->enfoconss[firstcons]);

         oldncutsstored = SCIPsepastoreGetNCutsStored(sepastore);
         oldnactiveconss = stat->nactiveconss;
         oldndomchgs = stat->nboundchgs + stat->nholechgs;

         /* check, if we want to use eager evaluation */
         if( (conshdlr->eagerfreq == 0 && conshdlr->nenfolpcalls == 0)
            || (conshdlr->eagerfreq > 0 && conshdlr->nenfolpcalls % conshdlr->eagerfreq == 0) )
            nusefulconss = nconss;

         /* because during constraint processing, constraints of this handler may be activated, deactivated,
          * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
          * external method; to avoid this, these changes will be buffered and processed after the method call
          */
         conshdlrDelayUpdates(conshdlr);

         /* start timing */
         SCIPclockStart(conshdlr->enfolptime, set);

         /* call external method */
         CHECK_OKAY( conshdlr->consenfolp(set->scip, conshdlr, conss, nconss, nusefulconss, result) );
         debugMessage(" -> enforcing returned result <%d>\n", *result);

         /* stop timing */
         SCIPclockStop(conshdlr->enfolptime, set);

         /* perform the cached constraint updates */
         CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, stat, prob) );

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_BRANCHED
            && *result != SCIP_INFEASIBLE
            && *result != SCIP_FEASIBLE )
         {
            errorMessage("enforcing method of constraint handler <%s> for LP solutions returned invalid result <%d>\n", 
               conshdlr->name, *result);
            return SCIP_INVALIDRESULT;
         }

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            conshdlr->nenfolpcalls++;
         if( *result == SCIP_CUTOFF )
            conshdlr->ncutoffs++;
         conshdlr->ncutsfound += SCIPsepastoreGetNCutsStored(sepastore) - oldncutsstored; /*lint !e776*/
         conshdlr->nconssfound += MAX(stat->nactiveconss - oldnactiveconss, 0);
         if( *result != SCIP_BRANCHED )
         {
            assert(tree->nchildren == 0);
            conshdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
         }
         else
            conshdlr->nchildren += tree->nchildren;
      }
   }

   return SCIP_OKAY;
}

/** calls enforcing method of constraint handler for pseudo solution for all constraints added after last
 *  conshdlrResetEnfo() call
 */
RETCODE SCIPconshdlrEnforcePseudoSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   Bool             objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(conshdlr->lastenfopsdomchgcount != stat->domchgcount
      || (0 <= conshdlr->lastnusefulenfoconss && conshdlr->lastnusefulenfoconss <= conshdlr->nusefulenfoconss));
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->consenfops != NULL )
   {
      int nconss;
      int nusefulconss;
      int firstcons;
      Bool pschanged;

      /* check, if this LP solution was already enforced */
      if( conshdlr->lastenfopsdomchgcount == stat->domchgcount )
      {
         /* all constraints that were not yet enforced on the new LP solution must be useful constraints, which means,
          * that the new constraints are the last constraints of the useful ones
          */
         nconss = conshdlr->nusefulenfoconss - conshdlr->lastnusefulenfoconss;
         nusefulconss = nconss;
         firstcons = conshdlr->lastnusefulenfoconss;
         pschanged = FALSE;
      }
      else
      {
         /* on a new pseudo solution, we want to enforce all constraints */
         nconss = conshdlr->nenfoconss;
         nusefulconss = conshdlr->nusefulenfoconss;
         firstcons = 0;
         pschanged = TRUE;
      }
      assert(firstcons >= 0);
      assert(firstcons + nconss <= conshdlr->nenfoconss);
      assert(nusefulconss <= nconss);

      /* constraint handlers without constraints should only be called once */
      if( nconss > 0 || (!conshdlr->needscons && pschanged) )
      {
         CONS** conss;
         Longint oldndomchgs;
         
         debugMessage("enforcing constraints %d to %d of %d constraints of handler <%s> (%s pseudo solution)\n",
            firstcons, firstcons + nconss - 1, conshdlr->nenfoconss, conshdlr->name, pschanged ? "new" : "old");

         /* remember the number of processed constraints on the current pseudo solution */
         conshdlr->lastenfopsdomchgcount = stat->domchgcount;
         conshdlr->lastnusefulenfoconss = conshdlr->nusefulenfoconss;

         /* get the array of the constraints to be processed */
         conss = &(conshdlr->enfoconss[firstcons]);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;

         /* check, if we want to use eager evaluation */
         if( (conshdlr->eagerfreq == 0 && conshdlr->nenfopscalls == 0)
            || (conshdlr->eagerfreq > 0 && conshdlr->nenfopscalls % conshdlr->eagerfreq == 0) )
            nusefulconss = nconss;

         /* because during constraint processing, constraints of this handler may be activated, deactivated,
          * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
          * external method; to avoid this, these changes will be buffered and processed after the method call
          */
         conshdlrDelayUpdates(conshdlr);

         /* start timing */
         SCIPclockStart(conshdlr->enfopstime, set);

         /* call external method */
         CHECK_OKAY( conshdlr->consenfops(set->scip, conshdlr, conss, nconss, nusefulconss, objinfeasible, result) );
         debugMessage(" -> enforcing returned result <%d>\n", *result);

         /* stop timing */
         SCIPclockStop(conshdlr->enfopstime, set);

         /* perform the cached constraint updates */
         CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, stat, prob) );

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_BRANCHED
            && *result != SCIP_SOLVELP
            && *result != SCIP_INFEASIBLE
            && *result != SCIP_FEASIBLE
            && *result != SCIP_DIDNOTRUN )
         {
            errorMessage("enforcing method of constraint handler <%s> for pseudo solutions returned invalid result <%d>\n", 
               conshdlr->name, *result);
            return SCIP_INVALIDRESULT;
         }

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            conshdlr->nenfopscalls++;
         else if( !objinfeasible )
         {
            errorMessage("enforcing method of constraint handler <%s> for pseudo solutions was skipped, even though the solution was not objective-infeasible\n", 
               conshdlr->name);
            return SCIP_INVALIDRESULT;
         }
         if( *result == SCIP_CUTOFF )
            conshdlr->ncutoffs++;
         if( *result != SCIP_BRANCHED )
         {
            assert(tree->nchildren == 0);
            conshdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
         }
         else
            conshdlr->nchildren += tree->nchildren;
      }
   }

   return SCIP_OKAY;
}

/** calls feasibility check method of constraint handler */
RETCODE SCIPconshdlrCheck(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< problem data */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->conscheck != NULL && (!conshdlr->needscons || conshdlr->ncheckconss > 0) )
   {
      debugMessage("checking %d constraints of handler <%s>\n", conshdlr->ncheckconss, conshdlr->name);

      /* because during constraint processing, constraints of this handler may be activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* call external method */
      CHECK_OKAY( conshdlr->conscheck(set->scip, conshdlr, conshdlr->checkconss, conshdlr->ncheckconss, 
                     sol, checkintegrality, checklprows, result) );
      debugMessage(" -> checking returned result <%d>\n", *result);
      
      /* perform the cached constraint updates */
      CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, stat, prob) );

      if( *result != SCIP_INFEASIBLE
         && *result != SCIP_FEASIBLE )
      {
         errorMessage("feasibility check of constraint handler <%s> returned invalid result <%d>\n", 
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** calls propagation method of constraint handler */
RETCODE SCIPconshdlrPropagate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< problem data */
   int              depth,              /**< depth of current node; -1 if preprocessing domain propagation */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(conshdlr->lastpropdomchgcount != stat->domchgcount
      || (0 <= conshdlr->lastnusefulpropconss && conshdlr->lastnusefulpropconss <= conshdlr->nusefulpropconss));
   assert(set != NULL);
   assert(stat != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->consprop != NULL
      && (!conshdlr->needscons || conshdlr->npropconss > 0)
      && (depth == -1 || (conshdlr->propfreq > 0 && depth % conshdlr->propfreq == 0)) )
   {
      int nconss;
      int nusefulconss;
      int firstcons;

      /* check, if the current domains were already propagated */
      if( conshdlr->lastpropdomchgcount == stat->domchgcount )
      {
         /* all constraints that were not yet propagated on the new domains must be useful constraints, which means,
          * that the new constraints are the last constraints of the useful ones
          */
         nconss = conshdlr->nusefulpropconss - conshdlr->lastnusefulpropconss;
         nusefulconss = nconss;
         firstcons = conshdlr->lastnusefulpropconss;
      }
      else
      {
         /* on new domains, we want to proprate all constraints */
         nconss = conshdlr->npropconss;
         nusefulconss = conshdlr->nusefulpropconss;
         firstcons = 0;
      }
      assert(firstcons >= 0);
      assert(firstcons + nconss <= conshdlr->npropconss);
      assert(nusefulconss <= nconss);

      /* constraint handlers without constraints should only be called once */
      if( nconss > 0 || (!conshdlr->needscons && conshdlr->lastpropdomchgcount != stat->domchgcount) )
      {
         CONS** conss;
         Longint oldndomchgs;
         
         debugMessage("propagating constraints %d to %d of %d constraints of handler <%s> (%s pseudo solution, %d useful)\n",
            firstcons, firstcons + nconss - 1, conshdlr->npropconss, conshdlr->name,
            conshdlr->lastpropdomchgcount == stat->domchgcount ? "old" : "new", nusefulconss);

         /* remember the number of processed constraints on the current domains */
         conshdlr->lastpropdomchgcount = stat->domchgcount;
         conshdlr->lastnusefulpropconss = conshdlr->nusefulpropconss;

         /* get the array of the constraints to be processed */
         conss = &(conshdlr->propconss[firstcons]);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;

         /* check, if we want to use eager evaluation */
         if( (conshdlr->eagerfreq == 0 && conshdlr->npropcalls == 0)
            || (conshdlr->eagerfreq > 0 && conshdlr->npropcalls % conshdlr->eagerfreq == 0) )
            nusefulconss = nconss;

         /* because during constraint processing, constraints of this handler may be activated, deactivated,
          * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
          * external method; to avoid this, these changes will be buffered and processed after the method call
          */
         conshdlrDelayUpdates(conshdlr);

         /* start timing */
         SCIPclockStart(conshdlr->proptime, set);

         /* call external method */
         CHECK_OKAY( conshdlr->consprop(set->scip, conshdlr, conss, nconss, nusefulconss, result) );
         debugMessage(" -> propagation returned result <%d>\n", *result);
      
         /* stop timing */
         SCIPclockStop(conshdlr->proptime, set);

         /* perform the cached constraint updates */
         CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, stat, prob) );

         /* check result code of callback method */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN )
         {
            errorMessage("propagation method of constraint handler <%s> returned invalid result <%d>\n", 
               conshdlr->name, *result);
            return SCIP_INVALIDRESULT;
         }

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            conshdlr->npropcalls++;
         if( *result == SCIP_CUTOFF )
            conshdlr->ncutoffs++;
         conshdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
      }
   }

   return SCIP_OKAY;
}

/** calls presolving method of constraint handler */
RETCODE SCIPconshdlrPresolve(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< problem data */
   int              nrounds,            /**< number of presolving rounds already done */
   int*             nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*             naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*             nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*             nchgbds,            /**< pointer to total number of variable bounds tightend of all presolvers */
   int*             naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*             ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*             nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*             nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*             nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nchgvartypes != NULL);
   assert(nchgbds != NULL);
   assert(naddholes != NULL);
   assert(ndelconss != NULL);
   assert(nupgdconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->conspresol != NULL
      && (!conshdlr->needscons || conshdlr->nconss > 0)
      && (conshdlr->maxprerounds == -1 || nrounds < conshdlr->maxprerounds) )
   {
      int nnewfixedvars;
      int nnewaggrvars;
      int nnewchgvartypes;
      int nnewchgbds;
      int nnewholes;
      int nnewdelconss;
      int nnewupgdconss;
      int nnewchgcoefs;
      int nnewchgsides;

      debugMessage("presolving %d constraints of handler <%s>\n", conshdlr->nconss, conshdlr->name);

      /* because during constraint processing, constraints of this handler may be activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* calculate the number of changes since last call */
      nnewfixedvars = *nfixedvars - conshdlr->lastnfixedvars;
      nnewaggrvars = *naggrvars - conshdlr->lastnaggrvars;
      nnewchgvartypes = *nchgvartypes - conshdlr->lastnchgvartypes;
      nnewchgbds = *nchgbds - conshdlr->lastnchgbds;
      nnewholes = *naddholes - conshdlr->lastnaddholes;
      nnewdelconss = *ndelconss - conshdlr->lastndelconss;
      nnewupgdconss = *nupgdconss - conshdlr->lastnupgdconss;
      nnewchgcoefs = *nchgcoefs - conshdlr->lastnchgcoefs;
      nnewchgsides = *nchgsides - conshdlr->lastnchgsides;
      
      /* remember the old number of changes */
      conshdlr->lastnfixedvars = *nfixedvars;
      conshdlr->lastnaggrvars = *naggrvars;
      conshdlr->lastnchgvartypes = *nchgvartypes;
      conshdlr->lastnchgbds = *nchgbds;
      conshdlr->lastnaddholes = *naddholes;
      conshdlr->lastndelconss = *ndelconss;
      conshdlr->lastnupgdconss = *nupgdconss;
      conshdlr->lastnchgcoefs = *nchgcoefs;
      conshdlr->lastnchgsides = *nchgsides;

      /* start timing */
      SCIPclockStart(conshdlr->presoltime, set);

      /* call external method */
      CHECK_OKAY( conshdlr->conspresol(set->scip, conshdlr, conshdlr->conss, conshdlr->nconss, nrounds,
                     nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
                     nnewdelconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
                     nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
                     ndelconss, nupgdconss, nchgcoefs, nchgsides, result) );

      /* stop timing */
      SCIPclockStop(conshdlr->presoltime, set);

      /* count the new changes */
      conshdlr->nfixedvars += *nfixedvars - conshdlr->lastnfixedvars;
      conshdlr->naggrvars += *naggrvars - conshdlr->lastnaggrvars;
      conshdlr->nchgvartypes += *nchgvartypes - conshdlr->lastnchgvartypes;
      conshdlr->nchgbds += *nchgbds - conshdlr->lastnchgbds;
      conshdlr->naddholes += *naddholes - conshdlr->lastnaddholes;
      conshdlr->ndelconss += *ndelconss - conshdlr->lastndelconss;
      conshdlr->nupgdconss += *nupgdconss - conshdlr->lastnupgdconss;
      conshdlr->nchgcoefs += *nchgcoefs - conshdlr->lastnchgcoefs;
      conshdlr->nchgsides += *nchgsides - conshdlr->lastnchgsides;

      /* perform the cached constraint updates */
      CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, stat, prob) );

      /* check result code of callback method */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_UNBOUNDED
         && *result != SCIP_SUCCESS
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         errorMessage("presolving method of constraint handler <%s> returned invalid result <%d>\n", 
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

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

/** gets description of constraint handler */
const char* SCIPconshdlrGetDesc(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->desc;
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

/** gets array with active constraints of constraint handler */
CONS** SCIPconshdlrGetConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->conss;
}

/** gets number of active constraints of constraint handler */
int SCIPconshdlrGetNConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nconss;
}

/** gets number of enabled constraints of constraint handler */
int SCIPconshdlrGetNEnabledConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nenabledconss;
}

/** gets time in seconds used for presolving in this constraint handler */
Real SCIPconshdlrGetPresolTime(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->presoltime);
}

/** gets time in seconds used for separation in this constraint handler */
Real SCIPconshdlrGetSepaTime(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->sepatime);
}

/** gets time in seconds used for LP enforcement in this constraint handler */
Real SCIPconshdlrGetEnfoLPTime(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->enfolptime);
}

/** gets time in seconds used for pseudo enforcement in this constraint handler */
Real SCIPconshdlrGetEnfoPSTime(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->enfopstime);
}

/** gets time in seconds used for propagation in this constraint handler */
Real SCIPconshdlrGetPropTime(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPclockGetTime(conshdlr->proptime);
}

/** gets number of calls to the constraint handler's separation method */
Longint SCIPconshdlrGetNSepaCalls(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nsepacalls;
}

/** gets number of calls to the constraint handler's LP enforcing method */
Longint SCIPconshdlrGetNEnfoLPCalls(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nenfolpcalls;
}

/** gets number of calls to the constraint handler's pseudo enforcing method */
Longint SCIPconshdlrGetNEnfoPSCalls(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nenfopscalls;
}

/** gets number of calls to the constraint handler's propagation method */
Longint SCIPconshdlrGetNPropCalls(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->npropcalls;
}

/** gets total number of times, this constraint handler detected a cutoff */
Longint SCIPconshdlrGetNCutoffs(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ncutoffs;
}

/** gets total number of cuts found by this constraint handler */
Longint SCIPconshdlrGetNCutsFound(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ncutsfound;
}

/** gets total number of additional constraints added by this constraint handler */
Longint SCIPconshdlrGetNConssFound(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nconssfound;
}

/** gets total number of domain reductions found by this constraint handler */
Longint SCIPconshdlrGetNDomredsFound(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ndomredsfound;
}

/** gets number of children created by this constraint handler */
Longint SCIPconshdlrGetNChildren(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nchildren;
}

/** gets maximum number of active constraints of constraint handler existing at the same time */
int SCIPconshdlrGetMaxNConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->maxnconss;
}

/** gets initial number of active constraints of constraint handler */
int SCIPconshdlrGetStartNConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->startnconss;
}

/** gets number of variables fixed in presolving method of constraint handler */
int SCIPconshdlrGetNFixedVars(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nfixedvars;
}

/** gets number of variables aggregated in presolving method of constraint handler */
int SCIPconshdlrGetNAggrVars(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->naggrvars;
}

/** gets number of variable types changed in presolving method of constraint handler */
int SCIPconshdlrGetNVarTypes(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nchgvartypes;
}

/** gets number of bounds changed in presolving method of constraint handler */
int SCIPconshdlrGetNChgBds(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nchgbds;
}

/** gets number of holes added to domains of variables in presolving method of constraint handler */
int SCIPconshdlrGetNAddHoles(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->naddholes;
}

/** gets number of constraints deleted in presolving method of constraint handler */
int SCIPconshdlrGetNDelConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ndelconss;
}

/** gets number of constraints upgraded in presolving method of constraint handler */
int SCIPconshdlrGetNUpgdConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nupgdconss;
}

/** gets number of coefficients changed in presolving method of constraint handler */
int SCIPconshdlrGetNChgCoefs(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nchgcoefs;
}

/** gets number of constraint sides changed in presolving method of constraint handler */
int SCIPconshdlrGetNChgSides(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nchgsides;
}

/** gets separation priority of constraint handler */
int SCIPconshdlrGetSepaPriority(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->sepapriority;
}

/** gets enforcing priority of constraint handler */
int SCIPconshdlrGetEnfoPriority(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->enfopriority;
}

/** gets checking priority of constraint handler */
int SCIPconshdlrGetCheckPriority(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->checkpriority;
}

/** gets separation frequency of constraint handler */
int SCIPconshdlrGetSepaFreq(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->sepafreq;
}

/** gets propagation frequency of constraint handler */
int SCIPconshdlrGetPropFreq(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->propfreq;
}

/** gets frequency of constraint handler for eager evaluations in separation, propagation and enforcement */
int SCIPconshdlrGetEagerFreq(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->eagerfreq;
}

/** needs constraint handler a constraint to be called? */
Bool SCIPconshdlrNeedsCons(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->needscons;
}

/** does the constraint handler perform presolving? */
Bool SCIPconshdlrDoesPresolve(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return (conshdlr->conspresol != NULL);
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
 * Constraint set change methods
 */

/** creates empty constraint set change data */
static
RETCODE conssetchgCreate(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   assert(conssetchg != NULL);
   assert(memhdr != NULL);

   ALLOC_OKAY( allocBlockMemory(memhdr, conssetchg) );
   (*conssetchg)->addedconss = NULL;
   (*conssetchg)->disabledconss = NULL;
   (*conssetchg)->addedconsssize = 0;
   (*conssetchg)->naddedconss = 0;
   (*conssetchg)->disabledconsssize = 0;
   (*conssetchg)->ndisabledconss = 0;

   return SCIP_OKAY;
}

/** releases all constraints of the constraint set change data */
static
RETCODE conssetchgRelease(
   CONSSETCHG*      conssetchg,         /**< constraint set change data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
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
         assert(!cons->active || cons->updatedeactivate);
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

/** frees constraint set change data and releases all included constraints */
RETCODE SCIPconssetchgFree(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(conssetchg != NULL);
   assert(memhdr != NULL);

   if( *conssetchg != NULL )
   {
      /* release constraints */
      CHECK_OKAY( conssetchgRelease(*conssetchg, memhdr, set) );

      /* free memory */
      freeBlockMemoryArrayNull(memhdr, &(*conssetchg)->addedconss, (*conssetchg)->addedconsssize);
      freeBlockMemoryArrayNull(memhdr, &(*conssetchg)->disabledconss, (*conssetchg)->disabledconsssize);
      freeBlockMemory(memhdr, conssetchg);
   }

   return SCIP_OKAY;
}

/** ensures, that addedconss array can store at least num entries */
static
RETCODE conssetchgEnsureAddedconssSize(
   CONSSETCHG*      conssetchg,         /**< constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(conssetchg != NULL);

   if( num > conssetchg->addedconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &conssetchg->addedconss, conssetchg->addedconsssize, newsize) );
      conssetchg->addedconsssize = newsize;
   }
   assert(num <= conssetchg->addedconsssize);

   return SCIP_OKAY;
}

/** ensures, that disabledconss array can store at least num entries */
static
RETCODE conssetchgEnsureDisabledconssSize(
   CONSSETCHG*      conssetchg,         /**< constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(conssetchg != NULL);

   if( num > conssetchg->disabledconsssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocBlockMemoryArray(memhdr, &conssetchg->disabledconss, conssetchg->disabledconsssize, newsize) );
      conssetchg->disabledconsssize = newsize;
   }
   assert(num <= conssetchg->disabledconsssize);

   return SCIP_OKAY;
}

/** adds constraint addition to constraint set changes, and captures constraint; activates constraint if the
 *  constraint set change data is currently active
 */
RETCODE SCIPconssetchgAddAddedCons(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   CONS*            cons,               /**< added constraint */
   int              depth,              /**< depth of constraint set change's node */
   Bool             active              /**< is the constraint set change currently active? */
   )
{
   assert(conssetchg != NULL);
   assert(cons != NULL);

   /* if constraint set change doesn't exist, create it */
   if( *conssetchg == NULL )
   {
      CHECK_OKAY( conssetchgCreate(conssetchg, memhdr) );
   }

   /* add constraint to the addedconss array */
   CHECK_OKAY( conssetchgEnsureAddedconssSize(*conssetchg, memhdr, set, (*conssetchg)->naddedconss+1) );
   (*conssetchg)->addedconss[(*conssetchg)->naddedconss] = cons;
   (*conssetchg)->naddedconss++;

   /* undelete constraint, if it was globally deleted in the past */
   cons->deleted = FALSE;

   /* capture constraint */
   SCIPconsCapture(cons);

   /* activate constraint, if node is active */
   if( active && !SCIPconsIsActive(cons) )
   {
      CHECK_OKAY( SCIPconsActivate(cons, set, stat, depth) );
      assert(SCIPconsIsActive(cons));
         
      /* remember, that this constraint set change data was resposible for the constraint's addition */
      cons->addconssetchg = *conssetchg;
      cons->addarraypos = (*conssetchg)->naddedconss-1;
   }

   return SCIP_OKAY;
}

/** adds constraint disabling to constraint set changes, and captures constraint */
RETCODE SCIPconssetchgAddDisabledCons(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< disabled constraint */
   )
{
   assert(conssetchg != NULL);
   assert(cons != NULL);

   /* if constraint set change doesn't exist, create it */
   if( *conssetchg == NULL )
   {
      CHECK_OKAY( conssetchgCreate(conssetchg, memhdr) );
   }

   /* add constraint to the disabledconss array */
   CHECK_OKAY( conssetchgEnsureDisabledconssSize(*conssetchg, memhdr, set, (*conssetchg)->ndisabledconss+1) );
   (*conssetchg)->disabledconss[(*conssetchg)->ndisabledconss] = cons;
   (*conssetchg)->ndisabledconss++;

   /* capture constraint */
   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

/** deactivates, deletes, and releases constraint from the addedconss array of the constraint set change data */
static
RETCODE conssetchgDelAddedCons(
   CONSSETCHG*      conssetchg,         /**< constraint set change to delete constraint from */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              arraypos            /**< position of constraint in disabledconss array */
   )
{
   CONS* cons;

   assert(conssetchg != NULL);
   assert(conssetchg->addedconss != NULL);
   assert(0 <= arraypos && arraypos < conssetchg->naddedconss);

   cons = conssetchg->addedconss[arraypos];

   debugMessage("delete added constraint <%s> at position %d from constraint set change data\n", cons->name, arraypos);

   /* release constraint */
   CHECK_OKAY( SCIPconsRelease(&conssetchg->addedconss[arraypos], memhdr, set) );

   /* move the last constraint of the addedconss array to the empty slot */
   if( arraypos < conssetchg->naddedconss-1 )
   {
      conssetchg->addedconss[arraypos] = conssetchg->addedconss[conssetchg->naddedconss-1];
      assert(conssetchg->addedconss[arraypos] != NULL);
      if( conssetchg->addedconss[arraypos]->addconssetchg == conssetchg )
      {
         assert(conssetchg->addedconss[arraypos]->addarraypos == conssetchg->naddedconss-1);
         conssetchg->addedconss[arraypos]->addarraypos = arraypos;
      }
   }
   conssetchg->naddedconss--;

   /* remove the link to the constraint set change data */
   if( cons->addconssetchg == conssetchg )
   {
      cons->addconssetchg = NULL;
      cons->addarraypos = -1;
   }

   return SCIP_OKAY;
}

/** deletes and releases deactivated constraint from the disabledconss array of the constraint set change data */
static
RETCODE conssetchgDelDisabledCons(
   CONSSETCHG*      conssetchg,         /**< constraint set change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   int              arraypos            /**< position of constraint in disabledconss array */
   )
{
   assert(conssetchg != NULL);
   assert(0 <= arraypos && arraypos < conssetchg->ndisabledconss);
   assert(conssetchg->disabledconss[arraypos] != NULL);

   debugMessage("delete disabled constraint <%s> at position %d from constraint set change data\n",
      conssetchg->disabledconss[arraypos]->name, arraypos);

   /* release constraint */
   CHECK_OKAY( SCIPconsRelease(&conssetchg->disabledconss[arraypos], memhdr, set) );

   /* move the last constraint of the disabledconss array to the empty slot */
   if( arraypos < conssetchg->ndisabledconss-1 )
   {
      assert(conssetchg->disabledconss[conssetchg->ndisabledconss-1] != NULL);
      conssetchg->disabledconss[arraypos] = conssetchg->disabledconss[conssetchg->ndisabledconss-1];
   }
   conssetchg->ndisabledconss--;

   return SCIP_OKAY;
}

/** applies constraint set change */
RETCODE SCIPconssetchgApply(
   CONSSETCHG*      conssetchg,         /**< constraint set change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   int              depth               /**< depth of constraint set change's node */
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
      assert(cons != NULL);
      assert(!cons->update);

      /* if constraint is already active, or if constraint is globally deleted, it can be removed from addedconss array */
      if( cons->active || cons->deleted )
      {
         CHECK_OKAY( conssetchgDelAddedCons(conssetchg, memhdr, set, i) );
         i--; /* the empty slot is now used by the last constraint, and naddedconss was decreased */
      }
      else
      {
         assert(cons->addconssetchg == NULL);
         assert(cons->addarraypos == -1);

         /* activate constraint */
         CHECK_OKAY( SCIPconsActivate(cons, set, stat, depth) );
         assert(cons->active);
         assert(!cons->update);
         
         /* remember, that this constraint set change data was resposible for the constraint's addition */
         cons->addconssetchg = conssetchg;
         cons->addarraypos = i;
      }
   }

   /* apply constraint disablings */
   for( i = 0; i < conssetchg->ndisabledconss; ++i )
   {
      cons = conssetchg->disabledconss[i];

      if( cons != NULL )
      {
         assert(!cons->update);

         /* if the constraint is disabled, we can permanently remove it from the disabledconss array */
         if( !cons->enabled )
         {
            debugMessage("constraint <%s> of handler <%s> was deactivated -> remove it from disabledconss array\n",
               cons->name, cons->conshdlr->name);
            
            /* release and remove constraint from the disabledconss array */
            CHECK_OKAY( conssetchgDelDisabledCons(conssetchg, memhdr, set, i) );
         }
         else
         {
            assert(cons->addarraypos >= 0);
            assert(!cons->deleted); /* deleted constraints must not be enabled! */
            CHECK_OKAY( SCIPconsDisable(conssetchg->disabledconss[i], set, stat) );
         }
         assert(!cons->update);
         assert(!cons->enabled);
      }
   }

   return SCIP_OKAY;
}

/** undoes constraint set change */
RETCODE SCIPconssetchgUndo(
   CONSSETCHG*      conssetchg,         /**< constraint set change to undo */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   )
{
   CONS* cons;
   int i;

   debugMessage("undoing constraint set changes at %p\n", conssetchg);

   if( conssetchg == NULL )
      return SCIP_OKAY;

   debugMessage(" -> %d constraint additions, %d constraint disablings\n", 
      conssetchg->naddedconss, conssetchg->ndisabledconss);

   /* undo constraint disablings */
   for( i = conssetchg->ndisabledconss-1; i >= 0; --i )
   {
      cons = conssetchg->disabledconss[i];
      if( cons != NULL )
      {
         assert(!cons->update);

         /* If the constraint is inactive, we can permanently remove it from the disabledconss array. It was deactivated
          * in the subtree of the current node but not reactivated on the switching way back to the current node, which
          * means, the deactivation was more global (i.e. valid on a higher level node) than the current node and the
          * disabling at the current node doesn't have any effect anymore.
          * If the constraint is already enabled, we need not to do anything. This may happen on a path A -> B,
          * if the constraint is disabled at node B, and while processing the subtree of B, it is also disabled at
          * the more global node A. Then on the switching path back to A, the node is enabled at node B (which is
          * actually wrong, since it now should be disabled in the whole subtree of A, but we cannot know this), and
          * again enabled at node A (where enabling is ignored). If afterwards, a subnode of B is processed, the
          * switching disables the constraint in node A, and the disabling is then removed from node B.
          */
         if( !cons->active )
         {
            debugMessage("constraint <%s> of handler <%s> was deactivated -> remove it from disabledconss array\n",
               cons->name, cons->conshdlr->name);
            
            /* release and remove constraint from the disabledconss array */
            CHECK_OKAY( conssetchgDelDisabledCons(conssetchg, memhdr, set, i) );
         }
         else if( !cons->enabled )
         {
            assert(cons->addarraypos >= 0);
            assert(!cons->deleted); /* deleted constraints must not be active! */
            CHECK_OKAY( SCIPconsEnable(cons, set, stat) );
         }
         assert(!cons->update);
         assert(!cons->active || cons->enabled);
      }
   }

   /* undo constraint additions */
   for( i = conssetchg->naddedconss-1; i >= 0; --i )
   {
      cons = conssetchg->addedconss[i];
      if( cons != NULL )
      {
         assert(!cons->update);

         /* If the constraint is already deactivated, we need not to do anything. This may happen on a path A -> B,
          * if the constraint is added at node B, and while processing the subtree of B, it is also added at
          * the more global node A. Then on the switching path back to A, the node is deactivated at node B (which is
          * actually wrong, since it now should be active in the whole subtree of A, but we cannot know this), and
          * again deactivated at node A (where deactivation is ignored). If afterwards, a subnode of B is processed, the
          * switching activates the constraint in node A, and the activation is then removed from node B.
          */
         if( cons->active )
         {
            assert(cons->addconssetchg == conssetchg);
            assert(cons->addarraypos == i);
            
            /* deactivate constraint */
            CHECK_OKAY( SCIPconsDeactivate(cons, set, stat) );
         
            /* unlink the constraint and the constraint set change */
            cons->addconssetchg = NULL;
            cons->addarraypos = -1;
         }
         assert(!cons->active);
         assert(!cons->update);
      }
   }

   return SCIP_OKAY;
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
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable,         /**< should the constraint be removed from the LP due to aging or cleanup? */
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
   (*cons)->transorigcons = NULL;
   (*cons)->addconssetchg = NULL;
   (*cons)->addarraypos = -1;
   (*cons)->consspos = -1;
   (*cons)->sepaconsspos = -1;
   (*cons)->enfoconsspos = -1;
   (*cons)->checkconsspos = -1;
   (*cons)->propconsspos = -1;
   (*cons)->activedepth = -2;
   (*cons)->nuses = 0;
   (*cons)->age = 0.0;
   (*cons)->nlockspos = 0;
   (*cons)->nlocksneg = 0;
   (*cons)->initial = initial;
   (*cons)->separate = separate;
   (*cons)->enforce = enforce;
   (*cons)->check = check;
   (*cons)->propagate = propagate;
   (*cons)->propenabled = propagate;
   (*cons)->local = local;
   (*cons)->modifiable = modifiable;
   (*cons)->removeable = removeable;
   (*cons)->original = original;
   (*cons)->active = FALSE;
   (*cons)->enabled = FALSE;
   (*cons)->obsolete = FALSE;
   (*cons)->deleted = FALSE;
   (*cons)->update = FALSE;
   (*cons)->updateactivate = FALSE;
   (*cons)->updatedeactivate = FALSE;
   (*cons)->updateenable = FALSE;
   (*cons)->updatedisable = FALSE;
   (*cons)->updatepropenable = FALSE;
   (*cons)->updatepropdisable = FALSE;
   (*cons)->updatedelete = FALSE;
   (*cons)->updateobsolete = FALSE;
   
   /* capture constraint */
   SCIPconsCapture(*cons);

   return SCIP_OKAY;
}

/** frees constraint data of a constraint, leaving the constraint itself as a zombie constraint; marks the constraint
 *  as deleted
 */
RETCODE SCIPconsFreeData(
   CONS*            cons,               /**< constraint to free */
   MEMHDR*          memhdr,             /**< block memory buffer */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(memhdr != NULL);
   assert(set != NULL);

   /* the constraint data must not be deleted, if the constraint is member of the update queue, because the
    * constraint handler method called in the update queue processing may use the constraint data
    */
   if( !cons->update )
   {
      /* free constraint data */
      if( cons->conshdlr->consdelete != NULL && cons->consdata != NULL )
      {
         CHECK_OKAY( cons->conshdlr->consdelete(set->scip, cons->conshdlr, cons, &cons->consdata) );
      }
      assert(cons->consdata == NULL);
   }
   cons->deleted = TRUE;

   return SCIP_OKAY;
}

/** frees a constraint */
RETCODE SCIPconsFree(
   CONS**           cons,               /**< constraint to free */
   MEMHDR*          memhdr,             /**< block memory buffer */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);
   assert((*cons)->nuses == 0);
   assert((*cons)->conshdlr != NULL);
   assert(!(*cons)->update);
   assert(!(*cons)->original || (*cons)->transorigcons == NULL);
   assert(memhdr != NULL);
   assert(set != NULL);

   /* free constraint data */
   CHECK_OKAY( SCIPconsFreeData(*cons, memhdr, set) );
   assert((*cons)->consdata == NULL);

   /* unlink transformed and original constraint */
   if( (*cons)->transorigcons != NULL )
   {
      assert(!(*cons)->original);
      assert((*cons)->transorigcons->original);
      assert((*cons)->transorigcons->transorigcons == *cons);

      (*cons)->transorigcons->transorigcons = NULL;
   }

   /* free constraint */
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
   SET*             set                 /**< global SCIP settings */
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

/** outputs constraint information to file stream */
RETCODE SCIPconsPrint(
   CONS*            cons,               /**< constraint to print */
   SET*             set,                /**< global SCIP settings */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   if( file == NULL )
      file = stdout;

   fprintf(file, "  [%s] <%s>: ", conshdlr->name, cons->name);
   if( conshdlr->consprint != NULL )
   {
      CHECK_OKAY( conshdlr->consprint(set->scip, conshdlr, cons, file) );
   }
   else
      fprintf(file, "constraint handler <%s> doesn't support printing constraints\n", conshdlr->name);

   return SCIP_OKAY;
}

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was created, or from the problem, if it was a problem constraint
 */
RETCODE SCIPconsDelete(
   CONS*            cons,               /**< constraint to delete */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob                /**< problem data */
   )
{
   assert(cons != NULL);
   assert(!cons->active || cons->updatedeactivate || cons->addarraypos >= 0);

   debugMessage("globally deleting constraint <%s>\n", cons->name);

   /* deactivate constraint, if it is currently active */
   if( cons->active && !cons->updatedeactivate )
   {
      CHECK_OKAY( SCIPconsDeactivate(cons, set, stat) );
   }
   assert(!cons->active || cons->updatedeactivate);
   assert(!cons->enabled || cons->updatedeactivate);

   /* mark constraint deleted */
   cons->deleted = TRUE;

   /* remove formerly active constraint from the conssetchg's/prob's addedconss/conss array */
   if( cons->addarraypos >= 0 )
   {
      if( cons->addconssetchg == NULL )
      {
         /* remove problem constraint from the problem */
         CHECK_OKAY( SCIPprobDelCons(prob, memhdr, set, stat, cons) );
      }
      else
      {
         assert(cons->addconssetchg->addedconss != NULL);
         assert(0 <= cons->addarraypos && cons->addarraypos < cons->addconssetchg->naddedconss);
         assert(cons->addconssetchg->addedconss[cons->addarraypos] == cons);
         
         /* remove constraint from the constraint set change addedconss array */
         CHECK_OKAY( conssetchgDelAddedCons(cons->addconssetchg, memhdr, set, cons->addarraypos) );
      }
   }

   return SCIP_OKAY;
}

/** gets and captures transformed constraint of a given original constraint; if the constraint is not yet transformed,
 *  a new transformed constraint for this constraint is created
 */
RETCODE SCIPconsTransform(
   CONS*            origcons,           /**< original constraint */
   MEMHDR*          memhdr,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   CONS**           transcons           /**< pointer to store the transformed constraint */
   )
{
   assert(origcons != NULL);
   assert(origcons->conshdlr != NULL);
   assert(origcons->original);
   assert(transcons != NULL);

   /* check, if the constraint is already transformed */
   if( origcons->transorigcons != NULL )
   {
      *transcons = origcons->transorigcons;
      SCIPconsCapture(*transcons);
   }
   else
   {
      /* create transformed constraint */
      if( origcons->conshdlr->constrans != NULL )
      {
         /* use constraint handler's own method to transform constraint */
         CHECK_OKAY( origcons->conshdlr->constrans(set->scip, origcons->conshdlr, origcons, transcons) );
      }
      else
      {
         /* create new constraint with empty constraint data */
         CHECK_OKAY( SCIPconsCreate(transcons, memhdr, origcons->name, origcons->conshdlr, NULL, origcons->initial,
               origcons->separate, origcons->enforce, origcons->check, origcons->propagate, 
               origcons->local, origcons->modifiable, origcons->removeable, FALSE) );
      }

      /* link original and transformed constraint */
      origcons->transorigcons = *transcons;
      (*transcons)->transorigcons = origcons;
   }
   assert(*transcons != NULL);

   return SCIP_OKAY;
}

/** gets associated transformed constraint of an original constraint, or NULL if no associated transformed constraint
 *  exists
 */
CONS* SCIPconsGetTransformed(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons->original);

   return cons->transorigcons;
}

/** activates constraint or marks constraint to be activated in next update */
RETCODE SCIPconsActivate(
   CONS*            cons,               /**< constraint */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   int              depth               /**< depth in the tree where the constraint activation takes place, or -1 for global problem */
   )
{
   assert(cons != NULL);
   assert(!cons->active);
   assert(!cons->updateactivate);
   assert(!cons->updatedeactivate);
   assert(!cons->updateenable);
   assert(!cons->updatedisable);
   assert(!cons->updatedelete);
   assert(!cons->updateobsolete);
   assert(cons->activedepth == -2);
   assert(cons->conshdlr != NULL);
   assert(depth == -1 || depth <= SCIPgetDepth(set->scip));

   if( cons->conshdlr->delayupdates )
   {
      debugMessage("delayed activation of constraint <%s> in constraint handler <%s> (depth %d)\n", 
         cons->name, cons->conshdlr->name, depth);
      cons->updateactivate = TRUE;
      cons->activedepth = depth;
      CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      CHECK_OKAY( conshdlrActivateCons(cons->conshdlr, set, stat, cons, depth) );
      assert(cons->active);
   }

   return SCIP_OKAY;
}

/** deactivates constraint or marks constraint to be deactivated in next update */
RETCODE SCIPconsDeactivate(
   CONS*            cons,               /**< constraint */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(cons != NULL);
   assert(cons->active);
   assert(!cons->updateactivate);
   assert(!cons->updatedeactivate);
   assert(cons->activedepth >= -1);
   assert(cons->conshdlr != NULL);

   if( cons->conshdlr->delayupdates )
   {
      debugMessage("delayed deactivation of constraint <%s> in constraint handler <%s>\n", 
         cons->name, cons->conshdlr->name);
      cons->updatedeactivate = TRUE;
      cons->activedepth = -2;
      CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      CHECK_OKAY( conshdlrDeactivateCons(cons->conshdlr, set, stat, cons) );
      assert(!cons->active);
   }

   return SCIP_OKAY;
}

/** enables constraint's separation, enforcing, and propagation capabilities or marks them to be enabled in next update */
RETCODE SCIPconsEnable(
   CONS*            cons,               /**< constraint */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(cons != NULL);
   assert(cons->active);
   assert(!cons->enabled);
   assert(!cons->updateactivate);
   assert(!cons->updateenable);
   assert(!cons->updatedisable);
   assert(cons->conshdlr != NULL);
   
   if( cons->conshdlr->delayupdates )
   {
      cons->updateenable = TRUE;
      CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      CHECK_OKAY( conshdlrEnableCons(cons->conshdlr, set, stat, cons) );
      assert(cons->enabled);
   }

   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities or marks them to be disabled in next update */
RETCODE SCIPconsDisable(
   CONS*            cons,               /**< constraint */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(cons != NULL);
   assert(cons->active);
   assert(cons->enabled);
   assert(!cons->updateactivate);
   assert(!cons->updateenable);
   assert(!cons->updatedisable);
   assert(cons->conshdlr != NULL);

   if( cons->conshdlr->delayupdates )
   {
      cons->updatedisable = TRUE;
      CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      CHECK_OKAY( conshdlrDisableCons(cons->conshdlr, set, stat, cons) );
      assert(!cons->enabled);
   }

   return SCIP_OKAY;
}

/** enables constraint's propagation capabilities or marks them to be enabled in next update */
RETCODE SCIPconsEnablePropagation(
   CONS*            cons,               /**< constraint */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);

   if( cons->updatepropenable || (cons->propenabled && !cons->updatepropdisable) )
      return SCIP_OKAY;

   if( cons->conshdlr->delayupdates )
   {
      cons->updatepropdisable = FALSE;
      cons->updatepropenable = TRUE;
      CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      CHECK_OKAY( conshdlrEnableConsPropagation(cons->conshdlr, set, cons) );
      assert(cons->propenabled);
   }

   return SCIP_OKAY;
}

/** disables constraint's propagation capabilities or marks them to be disabled in next update */
RETCODE SCIPconsDisablePropagation(
   CONS*            cons,               /**< constraint */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);

   if( cons->updatepropdisable || (!cons->propenabled && !cons->updatepropenable) )
      return SCIP_OKAY;

   if( cons->conshdlr->delayupdates )
   {
      cons->updatepropenable = FALSE;
      cons->updatepropdisable = TRUE;
      CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      assert(cons->update);
   }
   else
   {
      CHECK_OKAY( conshdlrDisableConsPropagation(cons->conshdlr, set, cons) );
      assert(!cons->propenabled);
   }

   return SCIP_OKAY;
}

/** adds given value to age of constraint, but age can never become negative;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *  if it's age exceeds the constraint age limit, makes constraint obsolete or marks constraint to be made obsolete
 *  in next update,
 */
RETCODE SCIPconsAddAge(
   CONS*            cons,               /**< constraint */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< problem data */
   Real             deltaage            /**< value to add to the constraint's age */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(!cons->updateactivate);
   assert(set != NULL);

   debugMessage("adding %g to age (%g) of constraint <%s> of handler <%s>\n",
      deltaage, cons->age, cons->name, cons->conshdlr->name);

   cons->age += deltaage;
   cons->age = MAX(cons->age, 0.0);

   if( !cons->check && consExceedsAgelimit(cons, set) )
   {
      if( cons->conshdlr->delayupdates )
      {
         cons->updatedelete = TRUE;
         CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
         assert(cons->update);
      }
      else
      {
         CHECK_OKAY( SCIPconsDelete(cons, memhdr, set, stat, prob) );
      }
   }
   else if( !cons->obsolete && consExceedsObsoleteage(cons, set) )
   {
      if( cons->conshdlr->delayupdates )
      {
         cons->updateobsolete = TRUE;
         CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
         assert(cons->update);
      }
      else
      {
         CHECK_OKAY( conshdlrMarkConsObsolete(cons->conshdlr, cons) );
         assert(cons->obsolete);
      }
   }

   return SCIP_OKAY;
}

/** increases age of constraint by 1.0;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *  if it's age exceeds the constraint age limit, makes constraint obsolete or marks constraint to be made obsolete
 *  in next update,
 */
RETCODE SCIPconsIncAge(
   CONS*            cons,               /**< constraint */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob                /**< problem data */
   )
{
   CHECK_OKAY( SCIPconsAddAge(cons, memhdr, set, stat, prob, 1.0) );

   return SCIP_OKAY;
}

/** resets age of constraint to zero;
 *  should be called
 *   - in constraint separation, if a cut was found for this constraint,
 *   - in constraint enforcing, if the constraint was violated, and
 *   - in constraint propagation, if a domain reduction was deduced;
 *  if it was obsolete, makes constraint useful again or marks constraint to be made useful again in next update
 */
RETCODE SCIPconsResetAge(
   CONS*            cons,               /**< constraint */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(!cons->updateactivate);

   debugMessage("resetting age %g of constraint <%s> of handler <%s>\n", cons->age, cons->name, cons->conshdlr->name);

   cons->age = 0.0;

   if( cons->obsolete )
   {
      if( cons->conshdlr->delayupdates )
      {
         cons->updateobsolete = TRUE;
         CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
         assert(cons->update);
      }
      else
      {
         CHECK_OKAY( conshdlrMarkConsUseful(cons->conshdlr, cons) );
         assert(!cons->obsolete);
      }
   }

   return SCIP_OKAY;
}

/** resolves the given conflicting bound, that was deduced by the given constraint, by putting all "reason" bounds
 *  leading to the deduction into the conflict queue with calls to SCIPaddConflictLb() and SCIPaddConflictUb()
 */
RETCODE SCIPconsResolvePropagation(
   CONS*            cons,               /**< constraint that deduced the assignment */
   SET*             set,                /**< global SCIP settings */
   VAR*             infervar,           /**< variable whose bound was deduced by the constraint */
   int              inferinfo,          /**< user inference information attached to the bound change */
   BOUNDTYPE        inferboundtype,     /**< bound that was deduced (lower or upper bound) */
   BDCHGIDX*        bdchgidx,           /**< bound change index, representing the point of time where change took place */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert((inferboundtype == SCIP_BOUNDTYPE_LOWER
         && SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) > SCIPvarGetLbGlobal(infervar))
      || (inferboundtype == SCIP_BOUNDTYPE_UPPER
         && SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < SCIPvarGetUbGlobal(infervar)));
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);

   if( conshdlr->consresprop != NULL )
   {
      CHECK_OKAY( conshdlr->consresprop(set->scip, conshdlr, cons, infervar, inferinfo, inferboundtype, bdchgidx,
            result) );
      
      /* check result code */
      if( *result != SCIP_SUCCESS && *result != SCIP_DIDNOTFIND )
      {
         errorMessage("propagation conflict resolving method of constraint handler <%s> returned invalid result <%d>\n", 
            conshdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }
   else
   {
      errorMessage("propagation conflict resolving method of constraint handler <%s> is not implemented\n", 
         conshdlr->name);
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** adds given values to lock status of the constraint and updates the rounding locks of the involved variables */
RETCODE SCIPconsAddLocks(
   CONS*            cons,               /**< constraint */
   SET*             set,                /**< global SCIP settings */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   int oldnlockspos;
   int oldnlocksneg;
   int updlockpos;
   int updlockneg;

   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(cons->conshdlr->conslock != NULL);
   assert(cons->nlockspos >= 0);
   assert(cons->nlocksneg >= 0);
   assert(-2 <= nlockspos && nlockspos <= 2);
   assert(-2 <= nlocksneg && nlocksneg <= 2);

   /* update the rounding locks */
   oldnlockspos = cons->nlockspos;
   oldnlocksneg = cons->nlocksneg;
   cons->nlockspos += nlockspos;
   cons->nlocksneg += nlocksneg;
   assert(cons->nlockspos >= 0);
   assert(cons->nlocksneg >= 0);

   /* check, if the constraint switched from unlocked to locked, or from locked to unlocked */
   updlockpos = (int)(cons->nlockspos > 0) - (int)(oldnlockspos > 0);
   updlockneg = (int)(cons->nlocksneg > 0) - (int)(oldnlocksneg > 0);

   /* lock the variables, if the constraint switched from unlocked to locked or from locked to unlocked */
   if( updlockpos != 0 || updlockneg != 0 )
   {
      CHECK_OKAY( cons->conshdlr->conslock(set->scip, cons->conshdlr, cons, updlockpos, updlockneg) );
   }

   return SCIP_OKAY;
}

/** checks single constraint for feasibility of the given solution */
RETCODE SCIPconsCheck(
   CONS*            cons,               /**< constraint to check */
   SET*             set,                /**< global SCIP settings */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   CONSHDLR* conshdlr;

   assert(cons != NULL);
   assert(set != NULL);

   conshdlr = cons->conshdlr;
   assert(conshdlr != NULL);
   
   /* call external method */
   CHECK_OKAY( conshdlr->conscheck(set->scip, conshdlr, &cons, 1, sol, checkintegrality, checklprows, result) );
   debugMessage(" -> checking returned result <%d>\n", *result);
   
   if( *result != SCIP_INFEASIBLE
      && *result != SCIP_FEASIBLE )
   {
      errorMessage("feasibility check of constraint handler <%s> on constraint <%s> returned invalid result <%d>\n", 
         conshdlr->name, cons->name, *result);
      return SCIP_INVALIDRESULT;
   }

   return SCIP_OKAY;
}

/** marks the constraint to be essential for feasibility */
RETCODE SCIPconsSetChecked(
   CONS*            cons,               /**< constraint */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);

   if( !cons->check )
   {
      cons->check = TRUE;

      /* if constraint is a problem constraint, lock variable roundings */
      if( cons->addconssetchg == NULL && cons->addarraypos >= 0 )
      {
         CHECK_OKAY( SCIPconsAddLocks(cons, set, +1, 0) );
      }

      /* if constraint is active, add it to the chckconss array of the constraint handler */
      if( cons->active )
      {
         CHECK_OKAY( conshdlrAddCheckcons(cons->conshdlr, set, cons) );
      }
   }

   return SCIP_OKAY;
}




/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given constraint */
DECL_HASHGETKEY(SCIPhashGetKeyCons)
{  /*lint --e{715}*/
   CONS* cons = (CONS*)elem;

   assert(cons != NULL);
   return cons->name;
}




/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPconsGetName
#undef SCIPconsGetHdlr
#undef SCIPconsGetData
#undef SCIPconsGetNUses
#undef SCIPconsGetActiveDepth
#undef SCIPconsIsActive
#undef SCIPconsIsEnabled
#undef SCIPconsIsPropagationEnabled
#undef SCIPconsIsDeleted
#undef SCIPconsIsObsolete
#undef SCIPconsGetAge
#undef SCIPconsIsInitial
#undef SCIPconsIsSeparated
#undef SCIPconsIsEnforced
#undef SCIPconsIsChecked
#undef SCIPconsIsPropagated
#undef SCIPconsIsGlobal
#undef SCIPconsIsLocal
#undef SCIPconsIsModifiable
#undef SCIPconsIsRemoveable
#undef SCIPconsIsInProb
#undef SCIPconsIsOriginal
#undef SCIPconsIsTransformed
#undef SCIPconsIsLockedPos
#undef SCIPconsIsLockedNeg
#undef SCIPconsIsLocked

/** returns the name of the constraint */
const char* SCIPconsGetName(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->name;
}

/** returns the constraint handler of the constraint */
CONSHDLR* SCIPconsGetHdlr(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->conshdlr;
}

/** returns the constraint data field of the constraint */
CONSDATA* SCIPconsGetData(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->consdata;
}

/** gets number of times, the constraint is currently captured */
int SCIPconsGetNUses(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->nuses;
}

/** for an active constraint, returns the depth in the tree at which the constraint was activated;
 *  for local active constraints, this is the depth at which the constraint is valid
 */
int SCIPconsGetActiveDepth(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsIsActive(cons));

   return cons->activedepth;
}

/** returns TRUE iff constraint is active in the current node */
Bool SCIPconsIsActive(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->updateactivate || (cons->active && !cons->updatedeactivate);
}

/** returns TRUE iff constraint is enabled in the current node */
Bool SCIPconsIsEnabled(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->updateenable || (cons->enabled && !cons->updatedisable);
}

/** returns TRUE iff constraint's propagation is enabled in the current node */
Bool SCIPconsIsPropagationEnabled(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return SCIPconsIsEnabled(cons)
      && (cons->updatepropenable || (cons->propenabled && !cons->updatepropdisable));
}

/** returns TRUE iff constraint is deleted or marked to be deleted */
Bool SCIPconsIsDeleted(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->deleted || cons->updatedelete;
}

/** returns TRUE iff constraint is marked obsolete */
Bool SCIPconsIsObsolete(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->updateobsolete || cons->obsolete;
}

/** gets age of constraint */
Real SCIPconsGetAge(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->age;
}

/** returns TRUE iff the LP relaxation of constraint should be in the initial LP */
Bool SCIPconsIsInitial(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->initial;
}

/** returns TRUE iff constraint should be separated during LP processing */
Bool SCIPconsIsSeparated(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->separate;
}

/** returns TRUE iff constraint should be enforced during node processing */
Bool SCIPconsIsEnforced(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->enforce;
}

/** returns TRUE iff constraint should be checked for feasibility */
Bool SCIPconsIsChecked(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->check;
}

/** returns TRUE iff constraint should be propagated during node processing */
Bool SCIPconsIsPropagated(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->propagate;
}

/** returns TRUE iff constraint is globally valid */
Bool SCIPconsIsGlobal(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return !cons->local;
}

/** returns TRUE iff constraint is only locally valid or not added to any (sub)problem */
Bool SCIPconsIsLocal(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->local;
}

/** returns TRUE iff constraint is modifiable (subject to column generation) */
Bool SCIPconsIsModifiable(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->modifiable;
}

/** returns TRUE iff constraint should be removed from the LP due to aging or cleanup */
Bool SCIPconsIsRemoveable(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->removeable;
}

/** returns TRUE iff constraint belongs to the global problem */
Bool SCIPconsIsInProb(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return (cons->addconssetchg == NULL && cons->addarraypos >= 0);
}

/** returns TRUE iff constraint is belonging to original space */
Bool SCIPconsIsOriginal(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->original;
}

/** returns TRUE iff constraint is belonging to transformed space */
Bool SCIPconsIsTransformed(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return !cons->original;
}

/** returns TRUE iff roundings for variables in constraint are locked */
Bool SCIPconsIsLockedPos(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return (cons->nlockspos > 0);
}

/** returns TRUE iff roundings for variables in constraint's negation are locked */
Bool SCIPconsIsLockedNeg(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return (cons->nlocksneg > 0);
}

/** returns TRUE iff roundings for variables in constraint or in constraint's negation are locked */
Bool SCIPconsIsLocked(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return (cons->nlockspos > 0 || cons->nlocksneg > 0);
}
