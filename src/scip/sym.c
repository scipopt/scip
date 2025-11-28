/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sym.c
 * @ingroup OTHER_CFILES
 * @brief  methods for symmetry handlers
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/clock.h"
#include "scip/set.h"
#include "scip/struct_sym.h"
#include "scip/sym.h"

/** internal method for creating a symmetry handler */
static
SCIP_RETCODE doSymhdlrCreate(
   SCIP_SYMHDLR**        symhdlr,            /**< pointer to symmetry handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of symmetry handler */
   const char*           desc,               /**< description of symmetry handler */
   int                   priority,           /**< priority of the symmetry handler */
   int                   proppriority,       /**< priority of the symmetry handler for propagation */
   int                   sepapriority,       /**< priority of the symmetry handler for separation */
   int                   presolpriority,     /**< priority of the symmetry handler for presolving */
   int                   propfreq,           /**< frequency for calling propagator of symmetry handler */
   int                   sepafreq,           /**< frequency for calling separator of symmetry handler */
   SCIP_Bool             delayprop,          /**< should propagation be delayed, if other sym-propagators found reductions? */
   SCIP_Bool             delaysepa,          /**< should separation be delayed, if other sym-separators found reductions? */
   int                   maxprerounds,       /**< maximal number of presolving rounds the symmetry handler participates in (-1: no limit) */
   SCIP_PROPTIMING       proptiming,         /**< positions in the node solving loop where propagation method of symmetry handlers should be executed */
   SCIP_PRESOLTIMING     presoltiming,       /**< timing mask of the symmetry handler's presolving method */
   SCIP_DECL_SYMHDLRTRYADD((*symtryadd)),    /**< addition method for symmetry method handler plugins */
   SCIP_DECL_SYMHDLRCOPY ((*symcopy)),       /**< copy method of symmetry handler */
   SCIP_DECL_SYMHDLRFREE ((*symfree)),       /**< destructor method of symmetry handler */
   SCIP_DECL_SYMHDLRINIT ((*syminit)),       /**< initialization method of symmetry handler */
   SCIP_DECL_SYMHDLREXIT ((*symexit)),       /**< deinitialization method of symmetry handler */
   SCIP_DECL_SYMHDLRTRANS((*symtrans)),      /**< transformation method of symmetry hanlder */
   SCIP_DECL_SYMHDLRSEPALP((*symsepalp)),    /**< separator for LP solutions */
   SCIP_DECL_SYMHDLRSEPASOL((*symsepasol)),  /**< separator for arbitrary primal solutions */
   SCIP_DECL_SYMHDLRPROP ((*symprop)),       /**< propagation method of symmetry handler */
   SCIP_DECL_SYMHDLRPRESOL((*sympresol)),    /**< presolving method of symmetry handler */
   SCIP_SYMHDLRDATA*     symhdlrdata         /**< symmetry handler data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(symhdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(symsepalp != NULL || symsepasol != NULL || sepafreq == -1);
   assert(symprop != NULL || propfreq == -1);

   SCIP_ALLOC( BMSallocMemory(symhdlr) );
   BMSclearMemory(*symhdlr);

   SCIP_ALLOC( BMSduplicateMemoryArray(&(*symhdlr)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*symhdlr)->desc, desc, strlen(desc)+1) );

   (*symhdlr)->tryaddpriority = priority;
   (*symhdlr)->sepapriority = sepapriority;
   (*symhdlr)->sepafreq = sepafreq;
   (*symhdlr)->proppriority = proppriority;
   (*symhdlr)->propfreq = propfreq;
   (*symhdlr)->presolpriority = presolpriority;
   (*symhdlr)->maxprerounds = maxprerounds;
   (*symhdlr)->delaysepa = delaysepa;
   (*symhdlr)->delayprop = delayprop;
   (*symhdlr)->proptiming = proptiming;
   (*symhdlr)->presoltiming = presoltiming;

   (*symhdlr)->symtryadd = symtryadd;
   (*symhdlr)->symcopy = symcopy;
   (*symhdlr)->symfree = symfree;
   (*symhdlr)->syminit = syminit;
   (*symhdlr)->symexit = symexit;
   (*symhdlr)->symtrans = symtrans;
   (*symhdlr)->symsepalp = symsepalp;
   (*symhdlr)->symsepasol = symsepasol;
   (*symhdlr)->symprop = symprop;
   (*symhdlr)->sympresol = sympresol;
   (*symhdlr)->symhdlrdata = symhdlrdata;

   SCIP_CALL( SCIPclockCreate(&(*symhdlr)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*symhdlr)->presoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*symhdlr)->sepatime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*symhdlr)->proptime, SCIP_CLOCKTYPE_DEFAULT) );

   (*symhdlr)->initialized = FALSE;
   (*symhdlr)->nsepacalls = 0;
   (*symhdlr)->npropcalls = 0;
   (*symhdlr)->ncutoffs = 0;
   (*symhdlr)->ncutsfound = 0;
   (*symhdlr)->ncutsapplied = 0;
   (*symhdlr)->nconssfound = 0;
   (*symhdlr)->ndomredsfound = 0;

   (*symhdlr)->lastnfixedvars = 0;
   (*symhdlr)->lastnaggrvars = 0;
   (*symhdlr)->lastnchgvartypes = 0;
   (*symhdlr)->lastnchgbds = 0;
   (*symhdlr)->lastnaddholes = 0;
   (*symhdlr)->lastndelconss = 0;
   (*symhdlr)->lastnaddconss = 0;
   (*symhdlr)->lastnupgdconss = 0;
   (*symhdlr)->lastnchgcoefs = 0;
   (*symhdlr)->lastnchgsides = 0;
   (*symhdlr)->nfixedvars = 0;
   (*symhdlr)->naggrvars = 0;
   (*symhdlr)->nchgvartypes = 0;
   (*symhdlr)->nchgbds = 0;
   (*symhdlr)->naddholes = 0;
   (*symhdlr)->ndelconss = 0;
   (*symhdlr)->naddconss = 0;
   (*symhdlr)->nupgdconss = 0;
   (*symhdlr)->nchgcoefs = 0;
   (*symhdlr)->nchgsides = 0;
   (*symhdlr)->npresolcalls = 0;
   (*symhdlr)->sepalpwasdelayed = FALSE;
   (*symhdlr)->sepasolwasdelayed = FALSE;
   (*symhdlr)->propwasdelayed = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/sepafreq", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "frequency for separating cuts (-1: never, 0: only in root node)",
         &(*symhdlr)->sepafreq, FALSE, sepafreq, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/propfreq", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "frequency for propagating domains (-1: never, 0: only in root node)",
         &(*symhdlr)->propfreq, FALSE, propfreq, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/proptiming", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing when constraint propagation should be called (%u:BEFORELP, %u:DURINGLPLOOP, %u:AFTERLPLOOP, %u:ALWAYS)", SCIP_PROPTIMING_BEFORELP, SCIP_PROPTIMING_DURINGLPLOOP, SCIP_PROPTIMING_AFTERLPLOOP, SCIP_PROPTIMING_ALWAYS);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         (int*)(&(*symhdlr)->proptiming), TRUE, (int) proptiming, (int) SCIP_PROPTIMING_BEFORELP, (int) SCIP_PROPTIMING_ALWAYS, NULL, NULL) ); /*lint !e713*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/maxprerounds", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "maximal number of presolving rounds the symmetry handler participates in (-1: no limit)",
         &(*symhdlr)->maxprerounds, TRUE, maxprerounds, -1, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/delaysepa", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should separation method be delayed, if other separators found cuts?",
         &(*symhdlr)->delaysepa, TRUE, delaysepa, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/delayprop", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should propagation method be delayed, if other propagators found reductions?",
         &(*symhdlr)->delayprop, TRUE, delayprop, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/presoltiming", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing mask of the symmetry handler's presolving method (%u:FAST, %u:MEDIUM, %u:EXHAUSTIVE, %u:FINAL)",
      SCIP_PRESOLTIMING_FAST, SCIP_PRESOLTIMING_MEDIUM, SCIP_PRESOLTIMING_EXHAUSTIVE, SCIP_PRESOLTIMING_FINAL);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         (int*)&(*symhdlr)->presoltiming, TRUE, (int) presoltiming, (int) SCIP_PRESOLTIMING_FAST, (int) SCIP_PRESOLTIMING_MAX, NULL, NULL) ); /*lint !e740 !e713*/

   return SCIP_OKAY;
}

/** creates a symmetry handler */
SCIP_RETCODE SCIPsymhdlrCreate(
   SCIP_SYMHDLR**        symhdlr,            /**< pointer to symmetry handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of symmetry handler */
   const char*           desc,               /**< description of symmetry handler */
   int                   priority,           /**< priority of the symmetry handler */
   int                   proppriority,       /**< priority of the symmetry handler for propagation */
   int                   sepapriority,       /**< priority of the symmetry handler for separation */
   int                   presolpriority,     /**< priority of the symmetry handler for presolving */
   int                   propfreq,           /**< frequency for calling propagator of symmetry handler */
   int                   sepafreq,           /**< frequency for calling separator of symmetry handler */
   SCIP_Bool             delayprop,          /**< should propagation be delayed, if other sym-propagators found reductions? */
   SCIP_Bool             delaysepa,          /**< should separation be delayed, if other sym-separators found reductions? */
   int                   maxprerounds,       /**< maximal number of presolving rounds the symmetry handler participates in (-1: no limit) */
   SCIP_PROPTIMING       proptiming,         /**< positions in the node solving loop where propagation method of symmetry handlers should be executed */
   SCIP_PRESOLTIMING     presoltiming,       /**< timing mask of the symmetry handler's presolving method */
   SCIP_DECL_SYMHDLRTRYADD((*symtryadd)),    /**< addition method for symmetry method handler plugins */
   SCIP_DECL_SYMHDLRCOPY ((*symcopy)),       /**< copy method of symmetry handler */
   SCIP_DECL_SYMHDLRFREE ((*symfree)),       /**< destructor method of symmetry handler */
   SCIP_DECL_SYMHDLRINIT ((*syminit)),       /**< initialization method of symmetry handler */
   SCIP_DECL_SYMHDLREXIT ((*symexit)),       /**< deinitialization method of symmetry handler */
   SCIP_DECL_SYMHDLRINITSOL((*syminitsol)),  /**< solving process initialization method of symmetry handler */
   SCIP_DECL_SYMHDLREXITSOL((*symexitsol)),  /**< solving process deinitialization method of symmetry handler */
   SCIP_DECL_SYMHDLRTRANS((*symtrans)),      /**< transformation method of symmetry hanlder */
   SCIP_DECL_SYMHDLRSEPALP((*symsepalp)),    /**< separator for LP solutions */
   SCIP_DECL_SYMHDLRSEPASOL((*symsepasol)),  /**< separator for arbitrary primal solutions */
   SCIP_DECL_SYMHDLRPROP ((*symprop)),       /**< propagation method of symmetry handler */
   SCIP_DECL_SYMHDLRPRESOL((*sympresol)),    /**< presolving method of symmetry handler */
   SCIP_SYMHDLRDATA*     symhdlrdata         /**< symmetry handler data */
   )
{
   assert(symhdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(symsepalp != NULL || symsepasol != NULL || sepafreq == -1);
   assert(symprop != NULL || propfreq == -1);

   SCIP_CALL_FINALLY( doSymhdlrCreate(symhdlr, set, messagehdlr, blkmem, name, desc, priority, proppriority,
         sepapriority, presolpriority, propfreq, sepafreq, delayprop, delaysepa, maxprerounds, proptiming,
         presoltiming, symtryadd, symcopy, symfree, syminit, symexit, symtrans, symsepalp, symsepasol,
         symprop, sympresol, symhdlrdata),
         (void) SCIPsymhdlrFree(symhdlr, set) );

   return SCIP_OKAY;
} /*lint !e715*/

/** copies the given symmetry handler to a new scip */
SCIP_RETCODE SCIPsymhdlrCopyInclude(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(symhdlr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( symhdlr->symcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including symmetry handler %s in subscip %p\n",
         SCIPsymhdlrGetName(symhdlr), (void*)set->scip);
      SCIP_CALL( symhdlr->symcopy(set->scip, symhdlr) );
   }
   return SCIP_OKAY;
}

/** calls destructor and frees memory of symmetry handler */
SCIP_RETCODE SCIPsymhdlrFree(
   SCIP_SYMHDLR**        symhdlr,            /**< pointer to symmetry handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(symhdlr != NULL);
   if( *symhdlr == NULL )
      return SCIP_OKAY;
   assert(!(*symhdlr)->initialized);
   assert(set != NULL);

   /* call destructor of symmetry handler */
   if( (*symhdlr)->symfree != NULL )
   {
      SCIP_CALL( (*symhdlr)->symfree(set->scip, *symhdlr) );
   }

   SCIPclockFree(&(*symhdlr)->proptime);
   SCIPclockFree(&(*symhdlr)->sepatime);
   SCIPclockFree(&(*symhdlr)->presoltime);
   SCIPclockFree(&(*symhdlr)->setuptime);

   BMSfreeMemoryArrayNull(&(*symhdlr)->name);
   BMSfreeMemoryArrayNull(&(*symhdlr)->desc);
   BMSfreeMemory(symhdlr);

   return SCIP_OKAY;
}

/** calls exit method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrExit(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(symhdlr != NULL);

   if( !symhdlr->initialized )
   {
      SCIPerrorMessage("symmetry handler <%s> not initialized\n", symhdlr->name);
      return SCIP_INVALIDCALL;
   }

   if( symhdlr->symexit != NULL )
   {
      /* start timing */
      SCIPclockStart(symhdlr->setuptime, set);

      SCIP_CALL( symhdlr->symexit(set->scip, symhdlr) );

      /* stop timing */
      SCIPclockStop(symhdlr->setuptime, set);
   }
   symhdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** calls exit method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrInit(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(symhdlr != NULL);
   assert(set != NULL);

   if( symhdlr->initialized )
   {
      SCIPerrorMessage("symmetry handler <%s> already initialized\n", symhdlr->name);
      return SCIP_INVALIDCALL;
   }

   if( symhdlr->syminit != NULL )
   {
      /* start timing */
      SCIPclockStart(symhdlr->setuptime, set);

      SCIP_CALL( symhdlr->syminit(set->scip, symhdlr) );

      /* stop timing */
      SCIPclockStop(symhdlr->setuptime, set);
   }
   symhdlr->initialized = TRUE;

   return SCIP_OKAY;
}

/** informs symmetry handler that the branch and bound process is being started */
SCIP_RETCODE SCIPsymhdlrInitsol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(symhdlr != NULL);
   assert(set != NULL);

   /* call solving process initialization method of symmetry handler */
   if( symhdlr->syminitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(symhdlr->setuptime, set);

      SCIP_CALL( symhdlr->syminitsol(set->scip, symhdlr) );

      /* stop timing */
      SCIPclockStop(symhdlr->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs symmetry handler rule that the branch and bound process data is being freed */
SCIP_RETCODE SCIPsymhdlrExitsol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             restart             /**< was this exit solve call triggered by a restart? */
   )
{
   assert(symhdlr != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of symmetry handler */
   if( symhdlr->symexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(symhdlr->setuptime, set);

      SCIP_CALL( symhdlr->symexitsol(set->scip, symhdlr, restart) );

      /* stop timing */
      SCIPclockStop(symhdlr->setuptime, set);
   }

   return SCIP_OKAY;
}

/** executes presolving method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrPresol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRESOLTIMING     timing,             /**< current presolving timing */
   int                   nrounds,            /**< number of presolving rounds already done */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(symhdlr != NULL);
   assert(set != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nchgvartypes != NULL);
   assert(nchgbds != NULL);
   assert(naddholes != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);
   assert(nupgdconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* @symtodo how to deal with exact SCIP? */
   if( symhdlr->sympresol == NULL || set->exact_enable )
      return SCIP_OKAY;

   /* check number of presolving rounds */
   if( symhdlr->maxprerounds >= 0 && symhdlr->npresolcalls >= symhdlr->maxprerounds )
      return SCIP_OKAY;

   /* check, if presolver should be delayed */
   if( symhdlr->presoltiming & timing )
   {
      int nnewfixedvars;
      int nnewaggrvars;
      int nnewchgvartypes;
      int nnewchgbds;
      int nnewaddholes;
      int nnewdelconss;
      int nnewaddconss;
      int nnewupgdconss;
      int nnewchgcoefs;
      int nnewchgsides;

      SCIPsetDebugMsg(set, "calling presolving method of symmetry handler <%s>\n", symhdlr->name);

      /* calculate the number of changes since last call */
      nnewfixedvars = *nfixedvars - symhdlr->lastnfixedvars;
      nnewaggrvars = *naggrvars - symhdlr->lastnaggrvars;
      nnewchgvartypes = *nchgvartypes - symhdlr->lastnchgvartypes;
      nnewchgbds = *nchgbds - symhdlr->lastnchgbds;
      nnewaddholes = *naddholes - symhdlr->lastnaddholes;
      nnewdelconss = *ndelconss - symhdlr->lastndelconss;
      nnewaddconss = *naddconss - symhdlr->lastnaddconss;
      nnewupgdconss = *nupgdconss - symhdlr->lastnupgdconss;
      nnewchgcoefs = *nchgcoefs - symhdlr->lastnchgcoefs;
      nnewchgsides = *nchgsides - symhdlr->lastnchgsides;

      /* remember the number of changes prior to the call of the presolver method of the symmetry handler */
      symhdlr->lastnfixedvars = *nfixedvars;
      symhdlr->lastnaggrvars = *naggrvars;
      symhdlr->lastnchgvartypes = *nchgvartypes;
      symhdlr->lastnchgbds = *nchgbds;
      symhdlr->lastnaddholes = *naddholes;
      symhdlr->lastndelconss = *ndelconss;
      symhdlr->lastnaddconss = *naddconss;
      symhdlr->lastnupgdconss = *nupgdconss;
      symhdlr->lastnchgcoefs = *nchgcoefs;
      symhdlr->lastnchgsides = *nchgsides;

      /* start timing */
      SCIPclockStart(symhdlr->presoltime, set);

      /* call external method */
      SCIP_CALL( symhdlr->sympresol(set->scip, symhdlr, nrounds, timing,
            nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewaddholes,
            nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
            nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
            ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

      /* stop timing */
      SCIPclockStop(symhdlr->presoltime, set);

      /* add/count the new changes */
      symhdlr->nfixedvars += *nfixedvars - symhdlr->lastnfixedvars;
      symhdlr->naggrvars += *naggrvars - symhdlr->lastnaggrvars;
      symhdlr->nchgvartypes += *nchgvartypes - symhdlr->lastnchgvartypes;
      symhdlr->nchgbds += *nchgbds - symhdlr->lastnchgbds;
      symhdlr->naddholes += *naddholes - symhdlr->lastnaddholes;
      symhdlr->ndelconss += *ndelconss - symhdlr->lastndelconss;
      symhdlr->naddconss += *naddconss - symhdlr->lastnaddconss;
      symhdlr->nupgdconss += *nupgdconss - symhdlr->lastnupgdconss;
      symhdlr->nchgcoefs += *nchgcoefs - symhdlr->lastnchgcoefs;
      symhdlr->nchgsides += *nchgsides - symhdlr->lastnchgsides;

      /* check result code of callback method */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_UNBOUNDED
         && *result != SCIP_SUCCESS
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         SCIPerrorMessage("symmetry handler <%s> returned invalid result <%d>\n", symhdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }

      /* increase the number of presolving calls, if the symmetry handler tried to find reductions */
      if( *result != SCIP_DIDNOTRUN )
         ++(symhdlr->npresolcalls);
   }

   return SCIP_OKAY;
}

/** calls try-add method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrTryadd(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int**                 symmetries,         /**< array of symmetries */
   int                   nsymmetries,        /**< number of symmetries in symmetries array */
   SYM_SYMTYPE           symtype,            /**< type of symmetry */
   SCIP_VAR**            symvars,            /**< variables on which symmetries act */
   int                   nsymvars,           /**< number of variables in symvars */
   SYM_GRAPH*            symgraph,           /**< symmetry detection graph */
   int                   id,                 /**< identifier of component for which symmetry handling shall be added */
   SCIP_Bool*            success             /**< pointer to store whether symmetry handling method could be added */
   )
{
   assert(symhdlr != NULL);
   assert(set != NULL);
   assert(success != NULL);

   SCIPsetDebugMsg(set, "try to add symmetry handler of type %s for symmetry component %d\n", symhdlr->name, id);

   /* start timing */
   SCIPclockStart(symhdlr->setuptime, set);

   SCIP_CALL( symhdlr->symtryadd(set->scip, symhdlr, symtype, symmetries, nsymmetries,
         symvars, nsymvars, symgraph, id, success) );

   /* end timing */
   SCIPclockStop(symhdlr->setuptime, set);

   if( *success )
      SCIPsetDebugMsg(set, "\t-->success\n");
   else
      SCIPsetDebugMsg(set, "\t-->symmetry handler %s is not applicable\n", symhdlr->name);

   return SCIP_OKAY;
}

/** gets name of symmetry handler */
const char* SCIPsymhdlrGetName(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->name;
}

/** gets user data of symmetry handler */
SCIP_SYMHDLRDATA* SCIPsymhdlrGetData(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->symhdlrdata;
}

/** compares two symmetry handlers w. r. to their try-add priority */
SCIP_DECL_SORTPTRCOMP(SCIPsymhdlrCompTryadd)
{  /*lint --e{715}*/
   return ((SCIP_SYMHDLR*)elem2)->tryaddpriority - ((SCIP_SYMHDLR*)elem1)->tryaddpriority;
}

/** compares two symmetry handlers w.r.t. their separation priority */
SCIP_DECL_SORTPTRCOMP(SCIPsymhdlrCompSepa)
{  /*lint --e{715}*/
   return ((SCIP_SYMHDLR*)elem2)->sepapriority - ((SCIP_SYMHDLR*)elem1)->sepapriority;
}

/** compares two symmetry handlers w.r.t. their propagation priority */
SCIP_DECL_SORTPTRCOMP(SCIPsymhdlrCompProp)
{  /*lint --e{715}*/
   return ((SCIP_SYMHDLR*)elem2)->proppriority - ((SCIP_SYMHDLR*)elem1)->proppriority;
}

/** compares two symmetry handlers w.r.t. their presolving priority */
SCIP_DECL_SORTPTRCOMP(SCIPsymhdlrCompPresol)
{  /*lint --e{715}*/
   return ((SCIP_SYMHDLR*)elem2)->presolpriority - ((SCIP_SYMHDLR*)elem1)->presolpriority;
}
