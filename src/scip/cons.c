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
#include "clock.h"


/** constraint handler */
struct ConsHdlr
{
   char*            name;               /**< name of constraint handler */
   char*            desc;               /**< description of constraint handler */
   int              sepapriority;       /**< priority of the constraint handler for separation */
   int              enfopriority;       /**< priority of the constraint handler for constraint enforcing */
   int              checkpriority;      /**< priority of the constraint handler for checking infeasibility */
   int              sepafreq;           /**< frequency for separating cuts; zero means to separate only in the root node */
   int              propfreq;           /**< frequency for propagating domains; zero means only preprocessing propagation */
   DECL_CONSFREE    ((*consfree));      /**< destructor of constraint handler */
   DECL_CONSINIT    ((*consinit));      /**< initialise constraint handler */
   DECL_CONSEXIT    ((*consexit));      /**< deinitialise constraint handler */
   DECL_CONSDELETE  ((*consdelete));    /**< free specific constraint data */
   DECL_CONSTRANS   ((*constrans));     /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA    ((*conssepa));      /**< separate cutting planes */
   DECL_CONSENFOLP  ((*consenfolp));    /**< enforcing constraints for LP solutions */
   DECL_CONSENFOPS  ((*consenfops));    /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHECK   ((*conscheck));     /**< check feasibility of primal solution */
   DECL_CONSPROP    ((*consprop));      /**< propagate variable domains */
   DECL_CONSPRESOL  ((*conspresol));    /**< presolving method */
   DECL_CONSENABLE  ((*consenable));    /**< enabling notification method */
   DECL_CONSDISABLE ((*consdisable));   /**< disabling notification method */
   CONSHDLRDATA*    conshdlrdata;       /**< constraint handler data */
   CONS**           conss;              /**< array with all active constraints */
   int              consssize;          /**< size of conss array */
   int              nconss;             /**< total number of active constraints */
   int              maxnconss;          /**< maximal number of active constraints existing at the same time */
   CONS**           sepaconss;          /**< array with active constraints that must be separated during LP processing */
   int              sepaconsssize;      /**< size of sepaconss array */
   int              nsepaconss;         /**< number of active constraints that may be separated during LP processing */
   int              nusefulsepaconss;   /**< number of non-obsolete active constraints that should be separated */
   CONS**           enfoconss;          /**< array with active constraints that must be enforced during node processing */
   int              enfoconsssize;      /**< size of enfoconss array */
   int              nenfoconss;         /**< number of active constraints that must be enforced during node processing */
   int              nusefulenfoconss;   /**< number of non-obsolete active constraints that must be enforced */
   CONS**           checkconss;         /**< array with active constraints that must be checked for feasibility */
   int              checkconsssize;     /**< size of checkconss array */
   int              ncheckconss;        /**< number of active constraints that must be checked for feasibility */
   int              nusefulcheckconss;  /**< number of non-obsolete active constraints that must be checked */
   CONS**           propconss;          /**< array with active constraints that must be propagated during node processing */
   int              propconsssize;      /**< size of propconss array */
   int              npropconss;         /**< number of active constraints that may be propagated during node processing */
   int              nusefulpropconss;   /**< number of non-obsolete active constraints that should be propagated */
   CONS**           updateconss;        /**< array with constraints that changed and have to be update in the handler */
   int              updateconsssize;    /**< size of updateconss array */
   int              nupdateconss;       /**< number of update constraints */
   int              nenabledconss;      /**< total number of enabled constraints of the handler */
   int              lastnsepaconss;     /**< number of already separated constraints after last conshdlrResetSepa() call */
   int              lastnenfoconss;     /**< number of already enforced constraints after last conshdlrResetEnfo() call */
   CLOCK*           presoltime;         /**< time used for presolving of this constraint handler */
   CLOCK*           sepatime;           /**< time used for separation of this constraint handler */
   CLOCK*           enfolptime;         /**< time used for LP enforcement of this constraint handler */
   CLOCK*           enfopstime;         /**< time used for pseudo enforcement of this constraint handler */
   CLOCK*           proptime;           /**< time used for propagation of this constraint handler */
   Longint          nsepacalls;         /**< number of times, the separator was called */
   Longint          nenfolpcalls;       /**< number of times, the LP enforcer was called */
   Longint          nenfopscalls;       /**< number of times, the pseudo enforcer was called */
   Longint          npropcalls;         /**< number of times, the propagator was called */
   Longint          ncutsfound;         /**< total number of cuts found by this constraint handler */
   Longint          nbranchings;        /**< number of times, the constraint handler performed a branching */
   int              lastnfixedvars;     /**< number of variables fixed before the last call to the presolver */
   int              lastnaggrvars;      /**< number of variables aggregated before the last call to the presolver */
   int              lastnchgvartypes;   /**< number of variable type changes before the last call to the presolver */
   int              lastnchgbds;        /**< number of variable bounds tightend before the last call to the presolver */
   int              lastnaddholes;      /**< number of domain holes added before the last call to the presolver */
   int              lastndelconss;      /**< number of deleted constraints before the last call to the presolver */
   int              lastnupgdconss;     /**< number of upgraded constraints before the last call to the presolver */
   int              lastnchgcoefs;      /**< number of changed coefficients before the last call to the presolver */
   int              lastnchgsides;      /**< number of changed left or right hand sides before the last call */
   int              nfixedvars;         /**< total number of variables fixed by this presolver */
   int              naggrvars;          /**< total number of variables aggregated by this presolver */
   int              nchgvartypes;       /**< total number of variable type changes by this presolver */
   int              nchgbds;            /**< total number of variable bounds tightend by this presolver */
   int              naddholes;          /**< total number of domain holes added by this presolver */
   int              ndelconss;          /**< total number of deleted constraints by this presolver */
   int              nupgdconss;         /**< total number of upgraded constraints by this presolver */
   int              nchgcoefs;          /**< total number of changed coefficients by this presolver */
   int              nchgsides;          /**< total number of changed left or right hand sides by this presolver */
   unsigned int     needscons:1;        /**< should the constraint handler be skipped, if no constraints are available? */
   unsigned int     initialized:1;      /**< is constraint handler initialized? */
   unsigned int     delayupdates:1;     /**< must the updates of the constraint arrays be delayed until processUpdate()? */
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


/** resizes conss array to be able to store at least num constraints */
static
RETCODE conshdlrEnsureConssMem(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
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

/** marks constraint to be obsolete; if constraint is not necessary for feasibility, it will be deleted completely;
 *  otherwise, it will be moved to the last part of the constraint arrays, such that it is checked, enforced, separated,
 *  and propagated after the useful constraints
 */
static
RETCODE conshdlrMarkConsObsolete(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   CONS*            cons                /**< constraint to be marked obsolete */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(!cons->obsolete);

   if( !cons->check )
   {
      CHECK_OKAY( SCIPconsDelete(cons, memhdr, set, prob) );
   }
   else
   {
      CONS* tmpcons;

      cons->obsolete = TRUE;

      /* constraint is needed for feasibility -> it should be moved to the last positions in the conss arrays */
      if( cons->active )
      {
         /* switch the last useful (non-obsolete) check constraint with this constraint */
         assert(0 <= cons->checkconsspos && cons->checkconsspos < conshdlr->nusefulcheckconss);

         tmpcons = conshdlr->checkconss[conshdlr->nusefulcheckconss-1];
         assert(tmpcons->checkconsspos == conshdlr->nusefulcheckconss-1);

         conshdlr->checkconss[conshdlr->nusefulcheckconss-1] = cons;
         conshdlr->checkconss[cons->checkconsspos] = tmpcons;
         tmpcons->checkconsspos = cons->checkconsspos;
         cons->checkconsspos = conshdlr->nusefulcheckconss-1;

         conshdlr->nusefulcheckconss--;
      }
      if( cons->enabled )
      {
         if( cons->separate )
         {
            /* switch the last useful (non-obsolete) sepa constraint with this constraint */
            assert(0 <= cons->sepaconsspos && cons->sepaconsspos < conshdlr->nusefulsepaconss);
            
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
            /* switch the last useful (non-obsolete) enfo constraint with this constraint */
            assert(0 <= cons->enfoconsspos && cons->enfoconsspos < conshdlr->nusefulenfoconss);
            
            tmpcons = conshdlr->enfoconss[conshdlr->nusefulenfoconss-1];
            assert(tmpcons->enfoconsspos == conshdlr->nusefulenfoconss-1);
            
            conshdlr->enfoconss[conshdlr->nusefulenfoconss-1] = cons;
            conshdlr->enfoconss[cons->enfoconsspos] = tmpcons;
            tmpcons->enfoconsspos = cons->enfoconsspos;
            cons->enfoconsspos = conshdlr->nusefulenfoconss-1;
            
            conshdlr->nusefulenfoconss--;
         }
         if( cons->propagate )
         {
            /* switch the last useful (non-obsolete) prop constraint with this constraint */
            assert(0 <= cons->propconsspos && cons->propconsspos < conshdlr->nusefulpropconss);
            
            tmpcons = conshdlr->propconss[conshdlr->nusefulpropconss-1];
            assert(tmpcons->propconsspos == conshdlr->nusefulpropconss-1);
            
            conshdlr->propconss[conshdlr->nusefulpropconss-1] = cons;
            conshdlr->propconss[cons->propconsspos] = tmpcons;
            tmpcons->propconsspos = cons->propconsspos;
            cons->propconsspos = conshdlr->nusefulpropconss-1;
            
            conshdlr->nusefulpropconss--;
         }
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
   assert(cons->check);

   cons->obsolete = FALSE;

   if( cons->active )
   {
      /* switch the first obsolete check constraint with this constraint */
      assert(conshdlr->nusefulcheckconss <= cons->checkconsspos && cons->checkconsspos < conshdlr->ncheckconss);
      
      tmpcons = conshdlr->checkconss[conshdlr->nusefulcheckconss];
      assert(tmpcons->checkconsspos == conshdlr->nusefulcheckconss);
      
      conshdlr->checkconss[conshdlr->nusefulcheckconss] = cons;
      conshdlr->checkconss[cons->checkconsspos] = tmpcons;
      tmpcons->checkconsspos = cons->checkconsspos;
      cons->checkconsspos = conshdlr->nusefulcheckconss;
      
      conshdlr->nusefulcheckconss++;
   }
   if( cons->enabled )
   {
      if( cons->separate )
      {
         /* switch the first obsolete sepa constraint with this constraint */
         assert(conshdlr->nusefulsepaconss <= cons->sepaconsspos && cons->sepaconsspos < conshdlr->nsepaconss);
         
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
         /* switch the first obsolete enfo constraint with this constraint */
         assert(conshdlr->nusefulenfoconss <= cons->enfoconsspos && cons->enfoconsspos < conshdlr->nenfoconss);
         
         tmpcons = conshdlr->enfoconss[conshdlr->nusefulenfoconss];
         assert(tmpcons->enfoconsspos == conshdlr->nusefulenfoconss);
         
         conshdlr->enfoconss[conshdlr->nusefulenfoconss] = cons;
         conshdlr->enfoconss[cons->enfoconsspos] = tmpcons;
         tmpcons->enfoconsspos = cons->enfoconsspos;
         cons->enfoconsspos = conshdlr->nusefulenfoconss;
         
         conshdlr->nusefulenfoconss++;
      }
      if( cons->propagate )
      {
         /* switch the first obsolete prop constraint with this constraint */
         assert(conshdlr->nusefulpropconss <= cons->propconsspos && cons->propconsspos < conshdlr->npropconss);
         
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

/** enables separation, enforcement, and propagation of constraint */
static
RETCODE conshdlrEnableCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
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
   }
      
   /* add constraint to the enforcement array */
   if( cons->enforce )
   {
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
      conshdlr->enfoconss[insertpos] = cons;
      cons->enfoconsspos = insertpos;
      conshdlr->nenfoconss++;
   }

   /* add constraint to the propagation array */
   if( cons->propagate )
   {
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
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to remove */
   )
{
   int delpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->enabled);
   assert((cons->separate) ^ (cons->sepaconsspos == -1));
   assert((cons->enforce) ^ (cons->enfoconsspos == -1));
   assert((cons->propagate) ^ (cons->propconsspos == -1));

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
      delpos = cons->sepaconsspos;
      if( !cons->obsolete )
      {
         assert(0 <= delpos && delpos < conshdlr->nusefulsepaconss);
         conshdlr->sepaconss[delpos] = conshdlr->sepaconss[conshdlr->nusefulsepaconss-1];
         conshdlr->sepaconss[delpos]->sepaconsspos = delpos;
         delpos = conshdlr->nusefulsepaconss-1;
         conshdlr->nusefulsepaconss--;
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

   /* delete constraint from the enforcement array */
   if( cons->enforce )
   {
      delpos = cons->enfoconsspos;
      if( !cons->obsolete )
      {
         assert(0 <= delpos && delpos < conshdlr->nusefulenfoconss);
         conshdlr->enfoconss[delpos] = conshdlr->enfoconss[conshdlr->nusefulenfoconss-1];
         conshdlr->enfoconss[delpos]->enfoconsspos = delpos;
         delpos = conshdlr->nusefulenfoconss-1;
         conshdlr->nusefulenfoconss--;
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

   /* delete constraint from the propagation array */
   if( cons->propagate )
   {
      delpos = cons->propconsspos;
      if( !cons->obsolete )
      {
         assert(0 <= delpos && delpos < conshdlr->nusefulpropconss);
         conshdlr->propconss[delpos] = conshdlr->propconss[conshdlr->nusefulpropconss-1];
         conshdlr->propconss[delpos]->propconsspos = delpos;
         delpos = conshdlr->nusefulpropconss-1;
         conshdlr->nusefulpropconss--;
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
RETCODE conshdlrActivateCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   int insertpos;

   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(!cons->active);
   assert(!cons->enabled);
   assert(cons->consspos == -1);
   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->checkconsspos == -1);
   assert(cons->propconsspos == -1);

   debugMessage("activate constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* activate constraint */
   CHECK_OKAY( conshdlrEnsureConssMem(conshdlr, set, conshdlr->nconss+1) );
   cons->active = TRUE;
   cons->consspos = conshdlr->nconss;
   conshdlr->conss[conshdlr->nconss] = cons;
   conshdlr->nconss++;
   conshdlr->maxnconss = MAX(conshdlr->maxnconss, conshdlr->nconss);

   /* add constraint to the check array */
   if( cons->check )
   {
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
   }

   /* enable separation, enforcement, and propagation of constraint */
   CHECK_OKAY( conshdlrEnableCons(conshdlr, set, cons) );

   return SCIP_OKAY;
}

/** deactivates and removes constraint from constraint handler's conss array */
static
RETCODE conshdlrDeactivateCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
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
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);
   assert(cons->active);
   assert(cons->consspos != -1);
   assert((cons->check) ^ (cons->checkconsspos == -1));

   debugMessage("deactivate constraint <%s> in constraint handler <%s>\n", cons->name, conshdlr->name);

   /* disable constraint */
   if( cons->enabled )
   {
      CHECK_OKAY( conshdlrDisableCons(conshdlr, set, cons) );
   }
   assert(!cons->enabled);

   /* delete constraint from the check array */
   if( cons->check )
   {
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

   assert(cons->consspos == -1);
   assert(cons->sepaconsspos == -1);
   assert(cons->enfoconsspos == -1);
   assert(cons->checkconsspos == -1);
   assert(cons->propconsspos == -1);

   return SCIP_OKAY;
}

/** processes all delayed updates of constraints:
 *  newly (de)activated constraints will be (de)activated;
 *  newly en/disabled constraints will be en/disabled;
 *  newly obsolete non-check constraints will be globally deleted;
 *  newly obsolete check constraints will be moved to the last positions in the sepa-, enfo-, check-, and prop-arrays;
 *  newly useful constraints will be moved to the first positions in the sepa-, enfo-, check-, and prop-arrays;
 */
static
RETCODE conshdlrProcessUpdates(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
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
         || cons->updateobsolete);

      debugMessage(" -> constraint <%s>: activate=%d, deactivate=%d, enable=%d, disable=%d, obsolete=%d (consdata=%p)\n",
         cons->name, cons->updateactivate, cons->updatedeactivate, cons->updateenable, cons->updatedisable,
         cons->updateobsolete, cons->consdata);

      if( cons->updateactivate )
      {
         assert(!cons->active);
         assert(!cons->updatedeactivate);
         assert(!cons->updateenable);
         assert(!cons->updatedisable);
         assert(!cons->updateobsolete);

         CHECK_OKAY( conshdlrActivateCons(conshdlr, set, cons) );
         assert(cons->active);
         cons->updateactivate = FALSE;
      }
      else if( cons->updatedeactivate )
      {
         assert(cons->active);

         CHECK_OKAY( conshdlrDeactivateCons(conshdlr, set, cons) );
         assert(!cons->active);
         cons->updatedeactivate = FALSE;
         cons->updateenable = FALSE;
         cons->updatedisable = FALSE;
         cons->obsolete = (cons->age >= set->consagelimit);
         cons->updateobsolete = FALSE;
      }
      else if( cons->updateenable )
      {
         assert(!cons->enabled);
         assert(!cons->updatedisable);

         CHECK_OKAY( conshdlrEnableCons(conshdlr, set, cons) );
         assert(cons->enabled);
         cons->updateenable = FALSE;
      }
      else if( cons->updatedisable )
      {
         assert(cons->enabled);

         CHECK_OKAY( conshdlrDisableCons(conshdlr, set, cons) );
         assert(!cons->enabled);
         cons->updatedisable = FALSE;
      }
      if( cons->updateobsolete )
      {
         if( !cons->obsolete && cons->age >= set->consagelimit )
         {
            /* the constraint's status must be switched to obsolete */
            CHECK_OKAY( conshdlrMarkConsObsolete(conshdlr, memhdr, set, prob, cons) );
         }
         else if( cons->obsolete && cons->age < set->consagelimit )
         {
            /* the constraint's status must be switched to useful */
            CHECK_OKAY( conshdlrMarkConsUseful(conshdlr, cons) );
         }
         cons->updateobsolete = FALSE;
      }
      assert(!cons->updateactivate && !cons->updatedeactivate && !cons->updateenable && !cons->updatedisable
         && !cons->updateobsolete);
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
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->delayupdates);
   
   debugMessage("constraint updates of constraint handler <%s> will be processed immediately\n", conshdlr->name);
   conshdlr->delayupdates = FALSE;

   CHECK_OKAY( conshdlrProcessUpdates(conshdlr, memhdr, set, prob) );
   assert(conshdlr->nupdateconss == 0);

   return SCIP_OKAY;
}

/** adds constraint to constraint handler's update constraint array and captures it */
static
RETCODE conshdlrAddUpdateCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(conshdlr != NULL);
   assert(set != NULL);
   assert(cons != NULL);
   assert(cons->conshdlr == conshdlr);

   if( !cons->update )
   {
      debugMessage("constraint <%s> of age %d has to be updated in constraint handler <%s> (consdata=%p)\n",
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
{
   return ((CONSHDLR*)elem2)->sepapriority - ((CONSHDLR*)elem1)->sepapriority;
}

/** compares two constraint handlers w. r. to their enforcing priority */
DECL_SORTPTRCOMP(SCIPconshdlrCompEnfo)
{
   return ((CONSHDLR*)elem2)->enfopriority - ((CONSHDLR*)elem1)->enfopriority;
}

/** compares two constraint handlers w. r. to their feasibility check priority */
DECL_SORTPTRCOMP(SCIPconshdlrCompCheck)
{
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
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   DECL_CONSINIT    ((*consinit)),      /**< initialise constraint handler */
   DECL_CONSEXIT    ((*consexit)),      /**< deinitialise constraint handler */
   DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA    ((*conssepa)),      /**< separate cutting planes */
   DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method */
   DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   char paramname[MAXSTRLEN];

   assert(conshdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert((sepafreq >= 0) ^ (conssepa == NULL));
   assert((propfreq >= 0) ^ (consprop == NULL));

   ALLOC_OKAY( allocMemory(conshdlr) );
   ALLOC_OKAY( duplicateMemoryArray(&(*conshdlr)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*conshdlr)->desc, desc, strlen(desc)+1) );
   (*conshdlr)->sepapriority = sepapriority;
   (*conshdlr)->enfopriority = enfopriority;
   (*conshdlr)->checkpriority = checkpriority;
   (*conshdlr)->sepafreq = sepafreq;
   (*conshdlr)->propfreq = propfreq;
   (*conshdlr)->consfree = consfree;
   (*conshdlr)->consinit = consinit;
   (*conshdlr)->consexit = consexit;
   (*conshdlr)->consdelete = consdelete;
   (*conshdlr)->constrans = constrans;
   (*conshdlr)->conssepa = conssepa;
   (*conshdlr)->consenfolp = consenfolp;
   (*conshdlr)->consenfops = consenfops;
   (*conshdlr)->conscheck = conscheck;
   (*conshdlr)->consprop = consprop;
   (*conshdlr)->conspresol = conspresol;
   (*conshdlr)->consenable = consenable;
   (*conshdlr)->consdisable = consdisable;
   (*conshdlr)->conshdlrdata = conshdlrdata;
   (*conshdlr)->conss = NULL;
   (*conshdlr)->consssize = 0;
   (*conshdlr)->nconss = 0;
   (*conshdlr)->maxnconss = 0;
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
   (*conshdlr)->lastnsepaconss = 0;
   (*conshdlr)->lastnenfoconss = 0;

   SCIPclockCreate(&(*conshdlr)->presoltime, SCIP_CLOCKTYPE_DEFAULT);
   SCIPclockCreate(&(*conshdlr)->sepatime, SCIP_CLOCKTYPE_DEFAULT);
   SCIPclockCreate(&(*conshdlr)->enfolptime, SCIP_CLOCKTYPE_DEFAULT);
   SCIPclockCreate(&(*conshdlr)->enfopstime, SCIP_CLOCKTYPE_DEFAULT);
   SCIPclockCreate(&(*conshdlr)->proptime, SCIP_CLOCKTYPE_DEFAULT);

   (*conshdlr)->nsepacalls = 0;
   (*conshdlr)->nenfolpcalls = 0;
   (*conshdlr)->nenfopscalls = 0;
   (*conshdlr)->npropcalls = 0;
   (*conshdlr)->ncutsfound = 0;
   (*conshdlr)->nbranchings = 0;
   (*conshdlr)->needscons = needscons;
   (*conshdlr)->initialized = FALSE;
   (*conshdlr)->delayupdates = FALSE;

   /* add parameters */
   sprintf(paramname, "conshdlr/%s/sepafreq", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, 
                  paramname, "frequency for separating cuts (-1: never, 0: only in root node)",
                  &(*conshdlr)->sepafreq, sepafreq, -1, INT_MAX, NULL, NULL) );
   sprintf(paramname, "conshdlr/%s/propfreq", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, 
                  paramname, "frequency for propagating domains (-1: never, 0: only in root node)",
                  &(*conshdlr)->propfreq, propfreq, -1, INT_MAX, NULL, NULL) );

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
   assert(SCIPstage(scip) == SCIP_STAGE_INIT);

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
      char s[MAXSTRLEN];
      sprintf(s, "Constraint handler <%s> already initialized", conshdlr->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( conshdlr->consinit != NULL )
   {
      CHECK_OKAY( conshdlr->consinit(scip, conshdlr) );
      SCIPclockReset(conshdlr->presoltime);
      SCIPclockReset(conshdlr->sepatime);
      SCIPclockReset(conshdlr->enfolptime);
      SCIPclockReset(conshdlr->enfopstime);
      SCIPclockReset(conshdlr->proptime);
      conshdlr->nsepacalls = 0;
      conshdlr->nenfolpcalls = 0;
      conshdlr->nenfopscalls = 0;
      conshdlr->npropcalls = 0;
      conshdlr->ncutsfound = 0;
      conshdlr->maxnconss = conshdlr->nconss;
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
      char s[MAXSTRLEN];
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
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   SEPASTORE*       sepastore,          /**< separation storage */
   int              actdepth,           /**< depth of active node */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(0 <= conshdlr->lastnsepaconss && conshdlr->lastnsepaconss <= conshdlr->nsepaconss);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->conssepa != NULL
      && ((actdepth == 0 && conshdlr->sepafreq == 0) || (conshdlr->sepafreq > 0 && actdepth % conshdlr->sepafreq == 0)) )
   {
      int nconss;
      int nusefulconss;
      int firstcons;

      if( conshdlr->lastnsepaconss > 0 )
      {
         /* all new constraints after the last conshdlrResetSepa() call must be useful constraints, which means, that
          * the new constraints are the last constraints of the useful ones
          */
         nconss = conshdlr->nsepaconss - conshdlr->lastnsepaconss;
         nusefulconss = nconss;
         firstcons = conshdlr->nusefulsepaconss - nconss;
      }
      else
      {
         /* immediately after a conshdlrResetSepa() call, we want to separate all constraints */
         nconss = conshdlr->nsepaconss;
         nusefulconss = conshdlr->nusefulsepaconss;
         firstcons = 0;
      }
      assert(firstcons >= 0);
      assert(firstcons + nconss <= conshdlr->nsepaconss);
      assert(nusefulconss <= nconss);

      if( !conshdlr->needscons || nconss > 0 )
      {
         CONS** conss;
         int oldncutsfound;

         debugMessage("separating constraints %d to %d of %d constraints of handler <%s>\n",
            firstcons, firstcons + nconss - 1, conshdlr->nsepaconss, conshdlr->name);

         conss = &(conshdlr->sepaconss[firstcons]);
         
         /* remember the current total number of found cuts */
         oldncutsfound = SCIPsepastoreGetNCutsFound(sepastore);

         /* because the during constraint processing, constraints of this handler may be activated, deactivated,
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
         CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, prob) );

         /* remember, that these constraints have already been processed */
         conshdlr->lastnsepaconss = conshdlr->nsepaconss;

         /* update the number of found cuts */
         conshdlr->ncutsfound += SCIPsepastoreGetNCutsFound(sepastore) - oldncutsfound;

         if( *result != SCIP_CUTOFF
            && *result != SCIP_SEPARATED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_CONSADDED
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN )
         {
            char s[MAXSTRLEN];
            sprintf(s, "separation method of constraint handler <%s> returned invalid result <%d>", 
               conshdlr->name, *result);
            errorMessage(s);
            return SCIP_INVALIDRESULT;
         }
         if( *result != SCIP_DIDNOTRUN )
            conshdlr->nsepacalls++;
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
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   SEPASTORE*       sepastore,          /**< separation storage */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(0 <= conshdlr->lastnenfoconss && conshdlr->lastnenfoconss <= conshdlr->nenfoconss);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->consenfolp != NULL )
   {
      int nconss;
      int nusefulconss;
      int firstcons;

      if( conshdlr->lastnenfoconss > 0 )
      {
         /* all new constraints after the last conshdlrResetEnfo() call must be useful constraints, which means, that
          * the new constraints are the last constraints of the useful ones
          */
         nconss = conshdlr->nenfoconss - conshdlr->lastnenfoconss;
         nusefulconss = nconss;
         firstcons = conshdlr->nusefulenfoconss - nconss;
      }
      else
      {
         /* immediately after a conshdlrResetEnfo() call, we want to enforce all constraints */
         nconss = conshdlr->nenfoconss;
         nusefulconss = conshdlr->nusefulenfoconss;
         firstcons = 0;
      }
      assert(firstcons >= 0);
      assert(firstcons + nconss <= conshdlr->nenfoconss);
      assert(nusefulconss <= nconss);

      if( !conshdlr->needscons || nconss > 0 )
      {
         CONS** conss;
         int oldncutsfound;

         debugMessage("enforcing constraints %d to %d of %d constraints of handler <%s>\n",
            firstcons, firstcons + nconss - 1, conshdlr->nenfoconss, conshdlr->name);

         conss = &(conshdlr->enfoconss[firstcons]);

         /* remember the current total number of found cuts */
         oldncutsfound = SCIPsepastoreGetNCutsFound(sepastore);

         /* because the during constraint processing, constraints of this handler may be activated, deactivated,
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
         CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, prob) );

         /* remember, that these constraints have already been processed */
         conshdlr->lastnenfoconss = conshdlr->nenfoconss;

         /* update the number of found cuts */
         conshdlr->ncutsfound += SCIPsepastoreGetNCutsFound(sepastore) - oldncutsfound;

         if( *result != SCIP_CUTOFF
            && *result != SCIP_BRANCHED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_CONSADDED
            && *result != SCIP_INFEASIBLE
            && *result != SCIP_FEASIBLE )
         {
            char s[MAXSTRLEN];
            sprintf(s, "enforcing method of constraint handler <%s> for LP solutions returned invalid result <%d>", 
               conshdlr->name, *result);
            errorMessage(s);
            return SCIP_INVALIDRESULT;
         }
         if( *result != SCIP_DIDNOTRUN )
         {
            conshdlr->nenfolpcalls++;
            if( *result == SCIP_BRANCHED )
               conshdlr->nbranchings++;
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
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   Bool             objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(0 <= conshdlr->lastnenfoconss && conshdlr->lastnenfoconss <= conshdlr->nenfoconss);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->consenfops != NULL )
   {
      int nconss;
      int nusefulconss;
      int firstcons;

      if( conshdlr->lastnenfoconss > 0 )
      {
         /* all new constraints after the last conshdlrResetEnfo() call must be useful constraints, which means, that
          * the new constraints are the last constraints of the useful ones
          */
         nconss = conshdlr->nenfoconss - conshdlr->lastnenfoconss;
         nusefulconss = nconss;
         firstcons = conshdlr->nusefulenfoconss - nconss;
      }
      else
      {
         /* immediately after a conshdlrResetEnfo() call, we want to enforce all constraints */
         nconss = conshdlr->nenfoconss;
         nusefulconss = conshdlr->nusefulenfoconss;
         firstcons = 0;
      }
      assert(firstcons >= 0);
      assert(firstcons + nconss <= conshdlr->nenfoconss);
      assert(nusefulconss <= nconss);

      if( !conshdlr->needscons || nconss > 0 )
      {
         CONS** conss;
         
         debugMessage("enforcing constraints %d to %d of %d constraints of handler <%s>\n",
            firstcons, firstcons + nconss - 1, conshdlr->nenfoconss, conshdlr->name);

         conss = &(conshdlr->enfoconss[firstcons]);

         /* because the during constraint processing, constraints of this handler may be activated, deactivated,
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
         CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, prob) );

         /* remember, that these constraints have already been processed */
         conshdlr->lastnenfoconss = conshdlr->nenfoconss;

         if( *result != SCIP_DIDNOTRUN
            && *result != SCIP_CUTOFF
            && *result != SCIP_BRANCHED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_CONSADDED
            && *result != SCIP_SOLVELP
            && *result != SCIP_INFEASIBLE
            && *result != SCIP_FEASIBLE )
         {
            char s[MAXSTRLEN];
            sprintf(s, "enforcing method of constraint handler <%s> for pseudo solutions returned invalid result <%d>", 
               conshdlr->name, *result);
            errorMessage(s);
            return SCIP_INVALIDRESULT;
         }
         if( *result != SCIP_DIDNOTRUN )
         {
            conshdlr->nenfopscalls++;
            if( *result == SCIP_BRANCHED )
               conshdlr->nbranchings++;
         }
         else if( !objinfeasible )
         {
            char s[MAXSTRLEN];
            sprintf(s, "enforcing method of constraint handler <%s> for pseudo solutions was skipped, even though the solution was not objective-infeasible", 
               conshdlr->name);
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
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,    /**< has integrality to be checked? */
   Bool             checklprows,         /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( conshdlr->conscheck != NULL && (!conshdlr->needscons || conshdlr->ncheckconss > 0) )
   {
      debugMessage("checking %d constraints of handler <%s>\n", conshdlr->ncheckconss, conshdlr->name);

      /* because the during constraint processing, constraints of this handler may be activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* call external method */
      CHECK_OKAY( conshdlr->conscheck(set->scip, conshdlr, conshdlr->checkconss, conshdlr->ncheckconss, 
                     sol, checkintegrality, checklprows, result) );
      debugMessage(" -> checking returned result <%d>\n", *result);
      
      /* perform the cached constraint updates */
      CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, prob) );

      if( *result != SCIP_INFEASIBLE
         && *result != SCIP_FEASIBLE )
      {
         char s[MAXSTRLEN];
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
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   int              actdepth,           /**< depth of active node; -1 if preprocessing domain propagation */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conshdlr != NULL);
   assert(conshdlr->nusefulsepaconss <= conshdlr->nsepaconss);
   assert(conshdlr->nusefulenfoconss <= conshdlr->nenfoconss);
   assert(conshdlr->nusefulcheckconss <= conshdlr->ncheckconss);
   assert(conshdlr->nusefulpropconss <= conshdlr->npropconss);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( conshdlr->consprop != NULL
      && (!conshdlr->needscons || conshdlr->npropconss > 0)
      && (actdepth == -1 || (conshdlr->propfreq > 0 && actdepth % conshdlr->propfreq == 0)) )
   {
      debugMessage("propagating %d constraints of handler <%s>\n", conshdlr->npropconss, conshdlr->name);

      /* because during constraint processing, constraints of this handler may be activated, deactivated,
       * enabled, disabled, marked obsolete or useful, which would change the conss array given to the
       * external method; to avoid this, these changes will be buffered and processed after the method call
       */
      conshdlrDelayUpdates(conshdlr);

      /* start timing */
      SCIPclockStart(conshdlr->proptime, set);

      /* call external method */
      CHECK_OKAY( conshdlr->consprop(set->scip, conshdlr, conshdlr->propconss, conshdlr->npropconss, 
                     conshdlr->nusefulpropconss, result) );
      debugMessage(" -> propagation returned result <%d>\n", *result);
      
      /* stop timing */
      SCIPclockStop(conshdlr->proptime, set);

      /* perform the cached constraint updates */
      CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, prob) );

      /* check result code of callback method */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         char s[MAXSTRLEN];
         sprintf(s, "propagation method of constraint handler <%s> returned invalid result <%d>", 
            conshdlr->name, *result);
         errorMessage(s);
         return SCIP_INVALIDRESULT;
      }
      if( *result != SCIP_DIDNOTRUN )
         conshdlr->npropcalls++;
   }

   return SCIP_OKAY;
}

/** calls presolving method of constraint handler */
RETCODE SCIPconshdlrPresolve(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
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
   assert(set->scip != NULL);
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
      && (!conshdlr->needscons || conshdlr->nconss > 0) )
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
      CHECK_OKAY( conshdlrForceUpdates(conshdlr, memhdr, set, prob) );

      /* check result code of callback method */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_UNBOUNDED
         && *result != SCIP_SUCCESS
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         char s[MAXSTRLEN];
         sprintf(s, "presolving method of constraint handler <%s> returned invalid result <%d>", 
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

/** gets total number of cuts found by this constraint handler */
Longint SCIPconshdlrGetNCutsFound(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->ncutsfound;
}

/** gets number of branchings performed by this constraint handler */
Longint SCIPconshdlrGetNBranchings(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->nbranchings;
}

/** gets maximum number of active constraints of constraint handler existing at the same time */
int SCIPconshdlrGetMaxNConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   return conshdlr->maxnconss;
}

/** resets maximum number of active constraints to current number of active constraints */
void SCIPconshdlrResetNMaxNConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(conshdlr != NULL);

   conshdlr->maxnconss = conshdlr->nconss;
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

/** deactivates, deletes, and releases constraint from the addedconss array of the constraint set change data */
static
RETCODE conssetchgDelAddedCons(
   CONSSETCHG*      conssetchg,         /**< constraint set change to delete constraint from */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to delete from addedconss array */
   )
{
   int arraypos;

   assert(conssetchg != NULL);
   assert(cons != NULL);

   debugMessage("delete added constraint <%s> from constraint set change data\n", cons->name);

   assert(!cons->active || cons->updatedeactivate);
   assert(!cons->enabled || cons->updatedeactivate);

   /* release and remove constraint from the addedconss array */
   arraypos = cons->arraypos;
   assert(0 <= arraypos && arraypos < conssetchg->naddedconss);
   assert(conssetchg->addedconss[arraypos] == cons);

   /* mark the constraint to be no longer in the problem */
   cons->node = NULL;
   cons->arraypos = -1;

   /* free constraint data, such that constraint exists only as a zombie constraint from now on */
   CHECK_OKAY( SCIPconsFreeData(cons, memhdr, set) );

   /* release constraint */
   CHECK_OKAY( SCIPconsRelease(&cons, memhdr, set) );

   conssetchg->addedconss[arraypos] = NULL;

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
         assert(cons->active || cons->updateactivate);
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
            CHECK_OKAY( conssetchgDelDisabledCons(conssetchg, memhdr, set, i) );
         }
         else if( !cons->check && arraypos < conssetchg->naddedconss && cons == conssetchg->addedconss[arraypos] )
         {
            debugMessage("constraint <%s> of handler <%s> was added and disabled at same node -> remove from both arrays\n",
               cons->name, cons->conshdlr->name);
            
            /* deactivate the just activated constraint */
            CHECK_OKAY( SCIPconsDeactivate(cons, set) );

            /* release and remove constraint from the addedconss array */
            CHECK_OKAY( conssetchgDelAddedCons(conssetchg, memhdr, set, cons) );

            /* release and remove constraint from the disabledconss array */
            CHECK_OKAY( conssetchgDelDisabledCons(conssetchg, memhdr, set, i) );
         }
         else
         {
            CHECK_OKAY( SCIPconsDisable(conssetchg->disabledconss[i], set) );
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
         CHECK_OKAY( SCIPconsEnable(cons, set) );
      }
   }

   /* undo constraint additions */
   for( i = conssetchg->naddedconss-1; i >= 0; --i )
   {
      cons = conssetchg->addedconss[i];
      if( cons != NULL )
      {
         CHECK_OKAY( SCIPconsDeactivate(cons, set) );
         assert(!cons->active || cons->updatedeactivate);
      }
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
   (*cons)->age = 0;
   (*cons)->consspos = -1;
   (*cons)->sepaconsspos = -1;
   (*cons)->enfoconsspos = -1;
   (*cons)->checkconsspos = -1;
   (*cons)->propconsspos = -1;
   (*cons)->arraypos = -1;
   (*cons)->separate = separate;
   (*cons)->enforce = enforce;
   (*cons)->check = check;
   (*cons)->propagate = propagate;
   (*cons)->original = original;
   (*cons)->active = FALSE;
   (*cons)->enabled = FALSE;
   (*cons)->obsolete = FALSE;
   (*cons)->update = FALSE;
   (*cons)->updateactivate = FALSE;
   (*cons)->updatedeactivate = FALSE;
   (*cons)->updateenable = FALSE;
   (*cons)->updatedisable = FALSE;
   (*cons)->updateobsolete = FALSE;
   
   /* capture constraint */
   SCIPconsCapture(*cons);

   return SCIP_OKAY;
}

/** frees constraint data of a constraint, leaving the constraint itself as a zombie constraint */
RETCODE SCIPconsFreeData(
   CONS*            cons,               /**< constraint to free */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
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
         CHECK_OKAY( cons->conshdlr->consdelete(set->scip, cons->conshdlr, &cons->consdata) );
      }
      assert(cons->consdata == NULL);
   }

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
   assert(!(*cons)->update);
   assert(memhdr != NULL);
   assert(set != NULL);

   /* free constraint data */
   CHECK_OKAY( SCIPconsFreeData(*cons, memhdr, set) );
   assert((*cons)->consdata == NULL);

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

/** globally removes constraint from all subproblems; removes constraint from the addedconss array of the node, where it
 *  was created, or from the problem, if it was a problem constraint;
 *  the constraint data is freed, and if the constraint is no longer used, it is freed completely
 */
RETCODE SCIPconsDelete(
   CONS*            cons,               /**< constraint to delete */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   )
{
   assert(cons != NULL);

   debugMessage("globally deleting constraint <%s>\n", cons->name);

   /* deactivate constraint, if it is currently active */
   if( cons->active && !cons->updatedeactivate )
   {
      CHECK_OKAY( SCIPconsDeactivate(cons, set) );
   }
   assert(!cons->active || cons->updatedeactivate);
   assert(!cons->enabled || cons->updatedeactivate);

   if( cons->node == NULL )
   {
      /* deactivate and remove problem constraint from the problem */
      CHECK_OKAY( SCIPprobDelCons(prob, memhdr, set, cons) );
   }
   else
   {
      /* deactivate and remove constraint from the node's addedconss array */
      CHECK_OKAY( conssetchgDelAddedCons(cons->node->conssetchg, memhdr, set, cons) );
   }

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
   assert(transcons != NULL);
   assert(memhdr != NULL);
   assert(origcons != NULL);

   if( origcons->conshdlr->constrans != NULL )
   {
      /* use constraints own method to transform constraint */
      CHECK_OKAY( origcons->conshdlr->constrans(set->scip, origcons->conshdlr, origcons, transcons) );
   }
   else
   {
      /* create new constraint with empty constraint data */
      CHECK_OKAY( SCIPconsCreate(transcons, memhdr, origcons->name, origcons->conshdlr, NULL,
                     origcons->separate, origcons->enforce, origcons->check, origcons->propagate, FALSE) );
   }
   assert(*transcons != NULL);

   return SCIP_OKAY;
}

/** activates constraint or marks constraint to be activated in next update */
RETCODE SCIPconsActivate(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(!cons->active);
   assert(!cons->updatedeactivate);
   assert(cons->conshdlr != NULL);

   if( cons->conshdlr->delayupdates )
   {
      cons->updateactivate = TRUE;
      CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
   }
   else
   {
      CHECK_OKAY( conshdlrActivateCons(cons->conshdlr, set, cons) );
      assert(cons->active);
   }

   return SCIP_OKAY;
}

/** deactivates constraint or marks constraint to be deactivated in next update */
RETCODE SCIPconsDeactivate(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->active);
   assert(!cons->updatedeactivate);
   assert(cons->conshdlr != NULL);
   
   if( cons->conshdlr->delayupdates )
   {
      cons->updatedeactivate = TRUE;
      CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
   }
   else
   {
      CHECK_OKAY( conshdlrDeactivateCons(cons->conshdlr, set, cons) );
      assert(!cons->active);
   }

   return SCIP_OKAY;
}

/** enables constraint's separation, enforcing, and propagation capabilities or marks them to be enabled in next update */
RETCODE SCIPconsEnable(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->active);
   assert(!cons->enabled);
   assert(cons->conshdlr != NULL);
   
   if( cons->conshdlr->delayupdates )
   {
      cons->updateenable = TRUE;
      CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
   }
   else
   {
      CHECK_OKAY( conshdlrEnableCons(cons->conshdlr, set, cons) );
      assert(cons->enabled);
   }

   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities or marks them to be disabled in next update */
RETCODE SCIPconsDisable(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->active);
   assert(cons->enabled);
   assert(cons->conshdlr != NULL);

   if( cons->conshdlr->delayupdates )
   {
      cons->updatedisable = TRUE;
      CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
   }
   else
   {
      CHECK_OKAY( conshdlrDisableCons(cons->conshdlr, set, cons) );
      assert(!cons->enabled);
   }

   return SCIP_OKAY;
}

/** increases age of constraint; should be called in constraint separation, if no cut was found for this constraint,
 *  in constraint enforcing, if constraint was feasible, and in constraint propagation, if no domain reduction was
 *  deduced;
 *  if it's age exceeds the constraint age limit, makes constraint obsolete or marks constraint to be made obsolete
 *  in next update,
 */
RETCODE SCIPconsIncAge(
   CONS*            cons,               /**< constraint */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);
   assert(set != NULL);

   debugMessage("increasing age (%d) of constraint <%s> of handler <%s>\n",
      cons->age, cons->name, cons->conshdlr->name);

   cons->age++;
   
   if( !cons->obsolete && cons->age >= set->consagelimit )
   {
      if( cons->conshdlr->delayupdates )
      {
         cons->updateobsolete = TRUE;
         CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      }
      else
      {
         CHECK_OKAY( conshdlrMarkConsObsolete(cons->conshdlr, memhdr, set, prob, cons) );
         assert(cons->obsolete);
      }
   }

   return SCIP_OKAY;
}

/** resets age of constraint to zero; should be called in constraint separation, if a cut was found for this constraint,
 *  in constraint enforcing, if the constraint was violated, and in constraint propagation, if a domain reduction was
 *  deduced;
 *  if it was obsolete, makes constraint useful again or marks constraint to be made useful again in next update
 */
RETCODE SCIPconsResetAge(
   CONS*            cons,               /**< constraint */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(cons != NULL);
   assert(cons->conshdlr != NULL);

   debugMessage("resetting age (%d) of constraint <%s> of handler <%s>\n",
      cons->age, cons->name, cons->conshdlr->name);

   cons->age = 0;

   if( cons->obsolete )
   {
      if( cons->conshdlr->delayupdates )
      {
         cons->updateobsolete = TRUE;
         CHECK_OKAY( conshdlrAddUpdateCons(cons->conshdlr, set, cons) );
      }
      else
      {
         CHECK_OKAY( conshdlrMarkConsUseful(cons->conshdlr, cons) );
         assert(!cons->obsolete);
      }
   }

   return SCIP_OKAY;
}

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

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

/** returns TRUE iff constraint is active in the current node */
Bool SCIPconsIsActive(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->updateactivate || (cons->active && !cons->updatedeactivate);
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

/** returns TRUE iff constraint is belonging to original problem */
Bool SCIPconsIsOriginal(
   CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   return cons->original;
}

#endif




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


