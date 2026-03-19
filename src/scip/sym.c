/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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
#include "scip/pub_cutpool.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip_cut.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_sym.h"
#include "scip/scip_var.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/struct_cutpool.h"
#include "scip/struct_mem.h"
#include "scip/struct_stat.h"
#include "scip/struct_sym.h"
#include "scip/sym.h"
#include "scip/struct_scip.h"
#include "symmetry/type_symmetry.h"

#define COMPRESSNVARSLB             25000    /**< lower bound on the number of variables above which compression could be performed */

/** initializes statistics and data of symmetry handler*/
static
void initSymhdlr(
   SCIP_SYMHDLR*         symhdlr,            /**< pointer to symmetry handler data structure */
   SCIP_Bool             ininitphase         /**< whether the method is called during SCIP's initialization phase */
   )
{
   assert(symhdlr != NULL);

   if( !symhdlr->initialized )
   {
      SCIPclockReset(symhdlr->setuptime);
      SCIPclockReset(symhdlr->presoltime);
      SCIPclockReset(symhdlr->sepatime);
      SCIPclockReset(symhdlr->proptime);
      SCIPclockReset(symhdlr->sbproptime);
      SCIPclockReset(symhdlr->resproptime);

      symhdlr->initialized = ininitphase;
      symhdlr->nsepacalls = 0;
      symhdlr->npropcalls = 0;
      symhdlr->nrespropcalls = 0;
      symhdlr->ncutoffs = 0;
      symhdlr->ncutsfound = 0;
      symhdlr->ncutsapplied = 0;
      symhdlr->nconssfound = 0;
      symhdlr->ndomredsfound = 0;

      symhdlr->lastnfixedvars = 0;
      symhdlr->lastnaggrvars = 0;
      symhdlr->lastnchgvartypes = 0;
      symhdlr->lastnchgbds = 0;
      symhdlr->lastnaddholes = 0;
      symhdlr->lastndelconss = 0;
      symhdlr->lastnaddconss = 0;
      symhdlr->lastnupgdconss = 0;
      symhdlr->lastnchgcoefs = 0;
      symhdlr->lastnchgsides = 0;
      symhdlr->nfixedvars = 0;
      symhdlr->naggrvars = 0;
      symhdlr->nchgvartypes = 0;
      symhdlr->nchgbds = 0;
      symhdlr->naddholes = 0;
      symhdlr->ndelconss = 0;
      symhdlr->naddconss = 0;
      symhdlr->nupgdconss = 0;
      symhdlr->nchgcoefs = 0;
      symhdlr->nchgsides = 0;
      symhdlr->npresolcalls = 0;
      symhdlr->lastsepanode = -1;
      symhdlr->nsepacallsatnode = 0;
      symhdlr->nseparootcalls = 0;
      symhdlr->ncutsfoundatnode = 0;
      symhdlr->nsepacutoffs = 0;
      symhdlr->nsepadomredsfound = 0;
      symhdlr->sepalpwasdelayed = FALSE;
      symhdlr->sepasolwasdelayed = FALSE;
      symhdlr->propwasdelayed = FALSE;

      symhdlr->symcomps = NULL;
      symhdlr->nsymcomps = 0;
      symhdlr->symcompssize = 0;
   }
}

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
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
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
   SCIP_DECL_SYMHDLRSEPALP((*symsepalp)),    /**< separator for LP solutions */
   SCIP_DECL_SYMHDLRSEPASOL((*symsepasol)),  /**< separator for arbitrary primal solutions */
   SCIP_DECL_SYMHDLRPROP ((*symprop)),       /**< propagation method of symmetry handler */
   SCIP_DECL_SYMHDLRRESPROP((*symresprop)),  /**< propagation conflict resolving method */
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

   SCIP_CALL( SCIPallocBlockMemory(set->scip, symhdlr) );
   BMSclearMemory(*symhdlr);

   SCIP_CALL( SCIPduplicateBlockMemoryArray(set->scip, &(*symhdlr)->name, name, strlen(name) + 1) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(set->scip, &(*symhdlr)->desc, desc, strlen(desc) + 1) );

   (*symhdlr)->tryaddpriority = priority;
   (*symhdlr)->symcopy = symcopy;
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
   (*symhdlr)->symfree = symfree;
   (*symhdlr)->syminit = syminit;
   (*symhdlr)->symexit = symexit;
   (*symhdlr)->syminitsol = syminitsol;
   (*symhdlr)->symexitsol = symexitsol;
   (*symhdlr)->symsepalp = symsepalp;
   (*symhdlr)->symsepasol = symsepasol;
   (*symhdlr)->symprop = symprop;
   (*symhdlr)->sympresol = sympresol;
   (*symhdlr)->symresprop = symresprop;
   (*symhdlr)->symhdlrdata = symhdlrdata;

   SCIP_CALL( SCIPclockCreate(&(*symhdlr)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*symhdlr)->presoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*symhdlr)->sepatime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*symhdlr)->proptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*symhdlr)->sbproptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*symhdlr)->resproptime, SCIP_CLOCKTYPE_DEFAULT) );

   initSymhdlr(*symhdlr, FALSE);

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/tryaddpriority", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, "priority of try-add function",
         &(*symhdlr)->tryaddpriority, TRUE, priority, INT_MIN/4, INT_MAX/4, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/proppriority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of propagator <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*symhdlr)->proppriority, TRUE, proppriority, INT_MIN/4, INT_MAX/4, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/sepafreq", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "frequency for separating cuts (-1: never, 0: only in root node)",
         &(*symhdlr)->sepafreq, FALSE, sepafreq, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/sepapriority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of separator <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*symhdlr)->sepapriority, TRUE, sepapriority, INT_MIN/4, INT_MAX/4,
         NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/maxbounddist", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying separator <%s> (0.0: only on current best node, 1.0: on all nodes)", name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*symhdlr)->maxbounddist, TRUE, maxbounddist, 0.0, 1.0, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/delaysepa", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should separator be delayed, if other separators found cuts?",
         &(*symhdlr)->delaysepa, TRUE, delaysepa, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/expbackoff", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "base for exponential increase of frequency at which separator <%s> is called (1: call at each multiple of frequency)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*symhdlr)->expbackoff, TRUE, 4, 1, 100, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/propfreq", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "frequency for propagating domains (-1: never, 0: only in root node)",
         &(*symhdlr)->propfreq, FALSE, propfreq, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/proptiming", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing when constraint propagation should be called (%u:BEFORELP, %u:DURINGLPLOOP, %u:AFTERLPLOOP, %u:ALWAYS)",
      SCIP_PROPTIMING_BEFORELP, SCIP_PROPTIMING_DURINGLPLOOP, SCIP_PROPTIMING_AFTERLPLOOP, SCIP_PROPTIMING_ALWAYS);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         (int*)(&(*symhdlr)->proptiming), TRUE, (int) proptiming, (int) SCIP_PROPTIMING_BEFORELP, (int) SCIP_PROPTIMING_ALWAYS, NULL, NULL) ); /*lint !e713*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/maxprerounds", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "maximal number of presolving rounds the symmetry handler participates in (-1: no limit)",
         &(*symhdlr)->maxprerounds, TRUE, maxprerounds, -1, INT_MAX, NULL, NULL) );

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
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
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
   SCIP_DECL_SYMHDLRSEPALP((*symsepalp)),    /**< separator for LP solutions */
   SCIP_DECL_SYMHDLRSEPASOL((*symsepasol)),  /**< separator for arbitrary primal solutions */
   SCIP_DECL_SYMHDLRPROP ((*symprop)),       /**< propagation method of symmetry handler */
   SCIP_DECL_SYMHDLRRESPROP((*symresprop)),  /**< propagation conflict resolving method */
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
         sepapriority, presolpriority, propfreq, sepafreq, delayprop, delaysepa, maxbounddist, maxprerounds, proptiming,
         presoltiming, symtryadd, symcopy, symfree, syminit, symexit, syminitsol, symexitsol, symsepalp, symsepasol,
         symprop, symresprop, sympresol, symhdlrdata), (void) SCIPsymhdlrFree(symhdlr, set) );

   return SCIP_OKAY;
} /*lint !e715*/

/** copies the given symmetry handler to a new SCIP */
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

   SCIPclockFree(&(*symhdlr)->resproptime);
   SCIPclockFree(&(*symhdlr)->sbproptime);
   SCIPclockFree(&(*symhdlr)->proptime);
   SCIPclockFree(&(*symhdlr)->sepatime);
   SCIPclockFree(&(*symhdlr)->presoltime);
   SCIPclockFree(&(*symhdlr)->setuptime);

   SCIPfreeBlockMemoryArray(set->scip, &(*symhdlr)->name, strlen((*symhdlr)->name) + 1);
   SCIPfreeBlockMemoryArray(set->scip, &(*symhdlr)->desc, strlen((*symhdlr)->desc) + 1);
   SCIPfreeBlockMemory(set->scip, symhdlr);

   return SCIP_OKAY;
}

/** calls exit method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrExit(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
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

      SCIP_CALL( symhdlr->symexit(set->scip, symhdlr, symhdlr->symcomps, symhdlr->nsymcomps) );

      /* stop timing */
      SCIPclockStop(symhdlr->setuptime, set);
   }
   BMSfreeBlockMemoryArrayNull(blkmem, &symhdlr->symcomps, symhdlr->symcompssize);
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

   /* reset statistics */
   initSymhdlr(symhdlr, TRUE);

   if( symhdlr->syminit != NULL )
   {
      /* start timing */
      SCIPclockStart(symhdlr->setuptime, set);

      SCIP_CALL( symhdlr->syminit(set->scip, symhdlr, symhdlr->symcomps, symhdlr->nsymcomps) );

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

      SCIP_CALL( symhdlr->syminitsol(set->scip, symhdlr, symhdlr->symcomps, symhdlr->nsymcomps) );

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

      SCIP_CALL( symhdlr->symexitsol(set->scip, symhdlr, symhdlr->symcomps, symhdlr->nsymcomps, restart) );

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
      SCIP_CALL( symhdlr->sympresol(set->scip, symhdlr, symhdlr->symcomps, symhdlr->nsymcomps, nrounds, timing,
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

/** calls propagation method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrProp(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             execdelayed,        /**< execute propagator even if it is marked to be delayed */
   SCIP_Bool             instrongbranching,  /**< are we currently doing strong branching? */
   SCIP_PROPTIMING       proptiming,         /**< current point in the node solving process */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(symhdlr != NULL);
   assert(symhdlr->propfreq >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(stat != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   if( symhdlr->symprop == NULL )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( ((depth == 0 && symhdlr->propfreq == 0) || (symhdlr->propfreq > 0 && depth % symhdlr->propfreq == 0))
         && !set->exact_enable )
   {
      if( !symhdlr->delayprop || execdelayed )
      {
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;

         SCIPsetDebugMsg(set, "executing propagator of symmetry handler <%s>\n", symhdlr->name);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;

         /* start timing */
         if( instrongbranching )
            SCIPclockStart(symhdlr->sbproptime, set);
         else
            SCIPclockStart(symhdlr->proptime, set);

         /* call external propagation method */
         SCIP_CALL( symhdlr->symprop(set->scip, symhdlr, symhdlr->symcomps, symhdlr->nsymcomps,
               proptiming, result) );

         /* stop timing */
         if( instrongbranching )
            SCIPclockStop(symhdlr->sbproptime, set);
         else
            SCIPclockStop(symhdlr->proptime, set);

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
            symhdlr->npropcalls++;
         if( *result == SCIP_CUTOFF )
            symhdlr->ncutoffs++;

         /* update domain reductions; therefore remove the domain
          * reduction counts which were generated in probing mode */
         symhdlr->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
         symhdlr->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED
            && *result != SCIP_DELAYNODE )
         {
            SCIPerrorMessage("propagation method of symmetry handler <%s> returned invalid result <%d>\n",
               symhdlr->name, *result);
            return SCIP_INVALIDRESULT;
         }
      }
      else
      {
         SCIPsetDebugMsg(set, "propagator of symmetry handler <%s> was delayed\n", symhdlr->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether propagator was delayed */
      symhdlr->propwasdelayed = (*result == SCIP_DELAYED);
   }
   else
      *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** calls LP separation method of symmetry handler's separator */
SCIP_RETCODE SCIPsymhdlrSepaLP(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   int                   depth,              /**< depth of current node */
   SCIP_Real             bounddist,          /**< current relative distance of local dual bound to global dual bound */
   SCIP_Bool             allowlocal,         /**< should the separator be asked to separate local cuts */
   SCIP_Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(symhdlr != NULL);
   assert(symhdlr->sepafreq >= -1);
   assert(0.0 <= symhdlr->maxbounddist && symhdlr->maxbounddist <= 1.0);
   assert(0.0 <= bounddist && bounddist <= 1.0);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(stat != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   if( symhdlr->symsepalp != NULL && SCIPsetIsLE(set, bounddist, symhdlr->maxbounddist) &&
       ( (depth == 0 && symhdlr->sepafreq != -1) ||
         (symhdlr->sepafreq > 0 && depth % symhdlr->sepafreq == 0 &&
            (symhdlr->expbackoff == 1 || SCIPsetIsIntegral(set, LOG2(depth * (1.0 / symhdlr->sepafreq)) / LOG2((SCIP_Real)symhdlr->expbackoff)))) ||
         symhdlr->sepalpwasdelayed) && !set->exact_enable
     )
   {
      if( (!symhdlr->delaysepa && !symhdlr->sepalpwasdelayed) || execdelayed )
      {
         SCIP_CUTPOOL* cutpool;
         SCIP_CUTPOOL* delayedcutpool;
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;
         int oldncutsfound;
         int oldnactiveconss;
         int ncutsfound;

         SCIPsetDebugMsg(set, "executing separator of symmetry handler <%s> on LP solution\n", symhdlr->name);

         cutpool = SCIPgetGlobalCutpool(set->scip);
         delayedcutpool = SCIPgetDelayedGlobalCutpool(set->scip);
         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;
         oldncutsfound = SCIPsepastoreGetNCuts(sepastore) + SCIPcutpoolGetNCuts(cutpool) + SCIPcutpoolGetNCuts(delayedcutpool);

         oldnactiveconss = stat->nactiveconss;

         /* reset the statistics for current node */
         if( symhdlr->lastsepanode != stat->ntotalnodes )
         {
            symhdlr->nsepacallsatnode = 0;
            symhdlr->ncutsfoundatnode = 0;
         }

         /* start timing */
         SCIPclockStart(symhdlr->sepatime, set);

         /* call external separation method */
         SCIP_CALL( symhdlr->symsepalp(set->scip, symhdlr, symhdlr->symcomps, symhdlr->nsymcomps,
               result, allowlocal, depth) );

         /* stop timing */
         SCIPclockStop(symhdlr->sepatime, set);

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
         {
            symhdlr->nsepacalls++;
            if( depth == 0 )
               symhdlr->nseparootcalls++;
            symhdlr->nsepacallsatnode++;
            symhdlr->lastsepanode = stat->ntotalnodes;
         }
         if( *result == SCIP_CUTOFF )
            symhdlr->nsepacutoffs++;

         ncutsfound = SCIPsepastoreGetNCuts(sepastore) + SCIPcutpoolGetNCuts(cutpool) + SCIPcutpoolGetNCuts(delayedcutpool) - oldncutsfound;
         symhdlr->ncutsfound += ncutsfound;
         symhdlr->ncutsfoundatnode += ncutsfound;
         symhdlr->nconssfound += MAX(stat->nactiveconss - oldnactiveconss, 0); /*lint !e776*/  /* cppcheck-suppress duplicateValueTernary */

         /* update domain reductions; therefore remove the domain
          * reduction counts which were generated in probing mode */
         symhdlr->nsepadomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
         symhdlr->nsepadomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_NEWROUND
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED )
         {
            SCIPerrorMessage("separation method of symmetry handler <%s> returned invalid result <%d>\n",
               symhdlr->name, *result);
            return SCIP_INVALIDRESULT;
         }
      }
      else
      {
         SCIPsetDebugMsg(set, "separator of symmetry handler <%s> was delayed\n", symhdlr->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether separator was delayed */
      symhdlr->sepalpwasdelayed = (*result == SCIP_DELAYED);
   }
   else
      *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** resolves the given conflicting bound, that was deduced by the given propagator, by putting all "reason" bounds
 *  leading to the deduction into the conflict queue with calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar();
 *
 *  @note it is sufficient to explain the relaxed bound change
 */
SCIP_RETCODE SCIPsymhdlrResolvePropagation(
   SCIP_SYMCOMP*         symcomp,            /**< symmetry component responsible for propagation */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             infervar,           /**< variable whose bound was deduced by the constraint */
   int                   inferinfo,          /**< user inference information attached to the bound change */
   SCIP_BOUNDTYPE        inferboundtype,     /**< bound that was deduced (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index, representing the point of time where change took place */
   SCIP_Real             relaxedbd,          /**< the relaxed bound */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_SYMHDLR* symhdlr;

   assert(symcomp != NULL);
   assert((inferboundtype == SCIP_BOUNDTYPE_LOWER
         && SCIPgetVarLbAtIndex(set->scip, infervar, bdchgidx, TRUE) > SCIPvarGetLbGlobal(infervar))
      || (inferboundtype == SCIP_BOUNDTYPE_UPPER
         && SCIPgetVarUbAtIndex(set->scip, infervar, bdchgidx, TRUE) < SCIPvarGetUbGlobal(infervar)));
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   symhdlr = symcomp->symhdlr;
   assert(symhdlr != NULL);

   if( symhdlr->symresprop != NULL )
   {
      /* start timing */
      SCIPclockStart(symhdlr->resproptime, set);

      SCIP_CALL( symhdlr->symresprop(set->scip, symhdlr, symcomp, infervar, inferinfo, inferboundtype, bdchgidx,
            relaxedbd, result) );

      /* stop timing */
      SCIPclockStop(symhdlr->resproptime, set);

      /* update statistic */
      symhdlr->nrespropcalls++;

      /* check result code */
      if( *result != SCIP_SUCCESS && *result != SCIP_DIDNOTFIND )
      {
         SCIPerrorMessage("propagation conflict resolving method of symmetry handler <%s> returned invalid result <%d>\n",
            symhdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }
   else
   {
      SCIPerrorMessage("propagation conflict resolving method of symmetry handler <%s> is not implemented\n", symhdlr->name);
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** calls primal solution separation method of symmetry handler's separator */
SCIP_RETCODE SCIPsymhdlrSepaSol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             allowlocal,         /**< should the separator allow local cuts */
   SCIP_Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(symhdlr != NULL);
   assert(symhdlr->sepafreq >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(stat != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   if( symhdlr->symsepasol != NULL &&
      ( (depth == 0 && symhdlr->sepafreq != -1) ||
         (symhdlr->sepafreq > 0 && depth % symhdlr->sepafreq == 0 &&
            (symhdlr->expbackoff == 1 || SCIPsetIsIntegral(set, LOG2(depth * (1.0 / symhdlr->sepafreq) / LOG2((SCIP_Real)symhdlr->expbackoff))))) ||
         symhdlr->sepasolwasdelayed ) && !set->exact_enable
     )
   {
      if( (!symhdlr->delaysepa && !symhdlr->sepasolwasdelayed) || execdelayed )
      {
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;
         int oldncutsfound;
         int oldnactiveconss;
         int ncutsfound;

         SCIPsetDebugMsg(set, "executing separator of symmetry handler <%s> on solution %p\n",
            symhdlr->name, (void*)sol);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;
         oldncutsfound = SCIPsepastoreGetNCuts(sepastore);
         oldnactiveconss = stat->nactiveconss;

         /* reset the statistics for current node */
         if( symhdlr->lastsepanode != stat->ntotalnodes )
         {
            symhdlr->nsepacallsatnode = 0;
            symhdlr->ncutsfoundatnode = 0;
         }

         /* start timing */
         SCIPclockStart(symhdlr->sepatime, set);

         /* call external separation method */
         SCIP_CALL( symhdlr->symsepasol(set->scip, symhdlr, sol, symhdlr->symcomps, symhdlr->nsymcomps,
               result, allowlocal, depth) );

         /* stop timing */
         SCIPclockStop(symhdlr->sepatime, set);

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
         {
            symhdlr->nsepacalls++;
            if( depth == 0 )
               symhdlr->nseparootcalls++;
            symhdlr->nsepacallsatnode++;
            symhdlr->lastsepanode = stat->ntotalnodes;
         }
         if( *result == SCIP_CUTOFF )
            symhdlr->nsepacutoffs++;

         ncutsfound = SCIPsepastoreGetNCuts(sepastore) - oldncutsfound;

         symhdlr->ncutsfound += ncutsfound;
         symhdlr->ncutsfoundatnode += ncutsfound;
         symhdlr->nconssfound += MAX(stat->nactiveconss - oldnactiveconss, 0); /*lint !e776*/  /* cppcheck-suppress duplicateValueTernary */

         /* update domain reductions; therefore remove the domain
          * reduction counts which were generated in probing mode */
         symhdlr->nsepadomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
         symhdlr->nsepadomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_NEWROUND
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED )
         {
            SCIPerrorMessage("execution method of separation method of symmetry handler <%s> returned invalid result <%d>\n",
               symhdlr->name, *result);
            return SCIP_INVALIDRESULT;
         }
      }
      else
      {
         SCIPsetDebugMsg(set, "separator of symmetry handler <%s> was delayed\n", symhdlr->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether separator was delayed */
      symhdlr->sepasolwasdelayed = (*result == SCIP_DELAYED);
   }
   else
      *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** calls try-add method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrTryAdd(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int**                 perms,              /**< array of (signed) permutations */
   int                   nperms,             /**< number of permutations in perms array */
   SYM_SYMTYPE           symtype,            /**< type of symmetry */
   SCIP_VAR**            permvars,           /**< variables on which permutations act */
   int                   npermvars,          /**< number of variables in permvars */
   SCIP_Real*            vardomcenter,       /**< domain center of variables (or NULL) */
   SCIP_HASHMAP*         permvarmap,         /**< map of variables to indices in permvars array */
   SYM_GRAPH*            symgraph,           /**< symmetry detection graph (or NULL) */
   int                   id,                 /**< identifier of component for which symmetry handling shall be added */
   SCIP_SYMCOMPDATA**    symcompdata,        /**< pointer to store data of symmetry component */
   int*                  naddedconss,        /**< pointer to store number of constraints added by symhdlr */
   SCIP_Bool             allowbdchgs,        /**< whether bound changed permitted (needed for stage EXITPRESOLVE) */
   int*                  nchgbds,            /**< pointer to store number of changed variable bounds */
   SCIP_Bool*            success             /**< pointer to store whether symmetry handling method could be added */
   )
{
   assert(symhdlr != NULL);
   assert(set != NULL);
   assert(permvarmap != NULL);
   assert(symcompdata != NULL);
   assert(naddedconss != NULL);
   assert(nchgbds != NULL);
   assert(success != NULL);

   SCIPsetDebugMsg(set, "try to add method of symmetry handler <%s> for symmetry component %d\n", symhdlr->name, id);

   /* start timing */
   SCIPclockStart(symhdlr->setuptime, set);

   *naddedconss = 0;
   *nchgbds = 0;

   SCIP_CALL( symhdlr->symtryadd(set->scip, symhdlr, symtype, perms, nperms, permvars, npermvars,
         vardomcenter, permvarmap, symgraph, id, symcompdata, naddedconss, allowbdchgs, nchgbds, success) );

   /* end timing */
   SCIPclockStop(symhdlr->setuptime, set);

   if( *success )
      SCIPsetDebugMsg(set, "\t-->success\n");
   else
      SCIPsetDebugMsg(set, "\t-->symmetry handler %s is not applicable\n", symhdlr->name);

   return SCIP_OKAY;
}

/** sets all separation related callbacks of the symmetry handler */
void SCIPsymhdlrSetSepa(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRSEPALP((*symsepalp)),    /**< separate cutting planes for LP solution */
   SCIP_DECL_SYMHDLRSEPASOL((*symsepasol)),  /**< separate cutting planes for arbitrary primal solution */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   sepapriority,       /**< priority of the symmetry handler for separation */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   SCIP_Bool             delaysepa           /**< should separation method be delayed, if other separators found cuts? */
   )
{
   assert(symhdlr != NULL);

   assert(symsepalp != NULL || symsepasol != NULL || sepafreq == -1);

   symhdlr->symsepalp = symsepalp;
   symhdlr->symsepasol = symsepasol;
   symhdlr->sepafreq = sepafreq;
   symhdlr->sepapriority = sepapriority;
   symhdlr->maxbounddist = maxbounddist;
   symhdlr->delaysepa = delaysepa;
}

/** sets both the propagation callback and the propagation frequency of the symmetry handler */
void SCIPsymhdlrSetProp(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRPROP ((*symprop)),       /**< propagate variable domains */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   int                   proppriority,       /**< priority of the symmetry handler for propagation */
   SCIP_PROPTIMING       timingmask          /**< positions in the node solving loop where propagators should be executed */
   )
{
   assert(symhdlr != NULL);

   assert(symprop != NULL || propfreq == -1);

   symhdlr->symprop = symprop;
   symhdlr->propfreq = propfreq;
   symhdlr->proppriority = proppriority;
   symhdlr->delayprop = delayprop;
   symhdlr->proptiming = timingmask;
}

/** sets copy metod of symmetry handler */
void SCIPsymhdlrSetCopy(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRCOPY ((*symcopy))        /**< copy symmetry handler */
   )
{
   assert(symhdlr != NULL);

   symhdlr->symcopy = symcopy;
}

/** sets destructor method of symmetry handler */
void SCIPsymhdlrSetFree(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRFREE ((*symfree))        /**< destructor of symmetry handler */
   )
{
   assert(symhdlr != NULL);

   symhdlr->symfree = symfree;
}

/** sets solving process initialization method of symmetry handler */
void SCIPsymhdlrSetInit(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRINIT ((*syminit))        /**< initialize symmetry handler */
   )
{
   assert(symhdlr != NULL);

   symhdlr->syminit = syminit;
}

/** sets deinitialization method of symmetry handler */
void SCIPsymhdlrSetExit(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLREXIT ((*symexit))        /**< deinitialize symmetry handler */
   )
{
   assert(symhdlr != NULL);

   symhdlr->symexit = symexit;
}

/** sets solving process initialization method of symmetry handler */
void SCIPsymhdlrSetInitsol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRINITSOL((*syminitsol))   /**< solving process initialization method of symmetry handler */
   )
{
   assert(symhdlr != NULL);

   symhdlr->syminitsol = syminitsol;
}

/** sets solving process deinitialization method of symmetry handler */
void SCIPsymhdlrSetExitsol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLREXITSOL((*symexitsol))   /**< solving process deinitialization method of symmetry handler */
   )
{
   assert(symhdlr != NULL);

   symhdlr->symexitsol = symexitsol;
}

/** sets presolving method of symmetry handler */
SCIP_RETCODE SCIPsymhdlrSetPresol(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRPRESOL((*sympresol)),    /**< presolving method of symmetry handler */
   int                   maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   int                   presolpriority,     /**< priority of the symmetry handler for presolving */
   SCIP_PRESOLTIMING     presoltiming        /**< timing mask of the constraint handler's presolving method */
   )
{
   assert(symhdlr != NULL);

   symhdlr->sympresol = sympresol;
   symhdlr->maxprerounds = maxprerounds;

   /* the interface change from delay flags to timings cannot be recognized at compile time: Exit with an appropriate
    * error message
    */
   if( presoltiming < SCIP_PRESOLTIMING_FAST || presoltiming > SCIP_PRESOLTIMING_MAX )
   {
      SCIPmessagePrintError("ERROR: 'PRESOLDELAY'-flag no longer available since SCIP 3.2, use an appropriate "
         "'SCIP_PRESOLTIMING' for <%s> symmetry handler instead.\n", symhdlr->name);

      return SCIP_PARAMETERWRONGVAL;
   }

   symhdlr->presoltiming = presoltiming;
   symhdlr->presolpriority = presolpriority;

   return SCIP_OKAY;
}

/** sets propagation conflict resolving method of symmetry handler */
void SCIPsymhdlrSetResprop(
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRRESPROP((*symresprop))   /**< propagation conflict resolving method */
   )
{
   assert(symhdlr != NULL);

   symhdlr->symresprop = symresprop;
}


/** gets data of symmetry component */
SCIP_SYMCOMPDATA* SCIPsymcompGetData(
   SCIP_SYMCOMP*         symcomp             /**< symmetry component */
   )
{
   assert(symcomp != NULL);

   return symcomp->symcompdata;
}

/** gets symmetry handler of symmetry component */
SCIP_SYMHDLR* SCIPsymcompGetHdlr(
   SCIP_SYMCOMP*         symcomp             /**< symmetry component */
   )
{
   assert(symcomp != NULL);

   return symcomp->symhdlr;
}

/** gets name of symmetry component */
const char* SCIPsymcompGetName(
   SCIP_SYMCOMP*         symcomp             /**< symmetry component */
   )
{
   assert(symcomp != NULL);

   return symcomp->name;
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

/** gets the total number of times, the propagator of the symmetry handler was called */
SCIP_Longint SCIPsymhdlrGetNPropCalls(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->npropcalls;
}

/** gets the total number of times, the propagator of the symmetry handler was called for resolving a propagation */
SCIP_Longint SCIPsymhdlrGetNRespropCalls(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->nrespropcalls;
}

/** gets total number of times, the symmetry handler's propagator detected a cutoff */
SCIP_Longint SCIPsymhdlrGetNCutoffs(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->ncutoffs;
}

/** gets total number of domain reductions found by this symmetry handler's propagator */
SCIP_Longint SCIPsymhdlrGetNDomredsFound(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->ndomredsfound;
}

/** returns the timing mask of the symmetry handler's propagator */
SCIP_PROPTIMING SCIPsymhdlrPropGetTimingmask(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->proptiming;
}

/** gets priority of symmetry handler's separator */
int SCIPsymhdlrSepaGetPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->sepapriority;
}

/** gets priority of symmetry handler's propagator */
int SCIPsymhdlrPropGetPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->proppriority;
}

/** was symmetry handler's separator delayed at the last call? */
SCIP_Bool SCIPsymhdlrSepaWasLPDelayed(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->sepalpwasdelayed;
}

/** was symmetry handler's propagator delayed at the last call? */
SCIP_Bool SCIPsymhdlrPropWasDelayed(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->propwasdelayed;
}

/** gets the total number of times, the separator of the symmetry handler was called */
SCIP_Longint SCIPsymhdlrGetNSepaCalls(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->nsepacalls;
}

/** gets the total number of times the symmetry handler's separator was called */
SCIP_Longint SCIPsymhdlrGetNCutsFound(
      SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->ncutsfound;
}

/** gets the total number of cutting planes added by the symmetry handler's separator */
SCIP_Longint SCIPsymhdlrGetNCutsAdded(
      SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->ncutsapplied;
}

/** gets time in seconds used in the symmetry handler's separator */
SCIP_Real SCIPsymhdlrGetSepaTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return SCIPclockGetTime(symhdlr->sepatime);
}

/** gets time in seconds used in the symmetry handler's propagator */
SCIP_Real SCIPsymhdlrGetPropTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return SCIPclockGetTime(symhdlr->proptime);
}


/** gets time in seconds used in the symmetry handler's propagator during strong branching */
SCIP_Real SCIPsymhdlrGetStrongBranchPropTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return SCIPclockGetTime(symhdlr->sbproptime);
}

/** gets time in seconds used in th symmetry handler's propagator for resolve propagation */
SCIP_Real SCIPsymhdlrGetRespropTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return SCIPclockGetTime(symhdlr->resproptime);
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

/** comparison method for sorting symmetry handlers w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPsymhdlrCompName)
{
   return strcmp(SCIPsymhdlrGetName((SCIP_SYMHDLR*)elem1), SCIPsymhdlrGetName((SCIP_SYMHDLR*)elem2));
}


/** creates a symmetry component */
SCIP_RETCODE SCIPcreateSymmetryComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMCOMP**        symcomp,            /**< pointer to symmetry component */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler active on symmetry component */
   SCIP_SYMCOMPDATA*     symcompdata,        /**< symmetry component data */
   int                   id                  /**< numerical identifier of symmetry component */
   )
{
   int len = 9;

   assert(scip != NULL);
   assert(symcomp != NULL);
   assert(symhdlr != NULL);
   assert(id >= 0);

   SCIP_CALL( SCIPallocBlockMemory(scip, symcomp) );
   (*symcomp)->symhdlr = symhdlr;
   (*symcomp)->symcompdata = symcompdata;

   if( id > 0 )
      len = 8 + (int) log10(id) + 1; /*lint !e747*/
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*symcomp)->name, len) );
   (void) snprintf((*symcomp)->name, SCIP_MAXSTRLEN, "symcomp%d", id);

   return SCIP_OKAY;
}

/** gets presolving priority of symmetry handler */
int SCIPsymhdlrGetPresolPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->presolpriority;
}

/** gets description of symmetry handler */
SCIP_EXPORT
const char* SCIPsymhdlrGetDesc(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->desc;
}

/** does the symmetry handler perform presolving? */
SCIP_Bool SCIPsymhdlrDoesPresolve(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->sympresol != NULL;
}

/** gets time in seconds used for setting up this symmetry handler for new stages */
SCIP_Real SCIPsymhdlrGetSetupTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return SCIPclockGetTime(symhdlr->setuptime);
}

/** gets time in seconds used in this symmetry handler for presolving */
SCIP_Real SCIPsymhdlrGetPresolTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return SCIPclockGetTime(symhdlr->presoltime);
}

/** gets number of times the symmetry handler was called in presolving and tried to find reductions */
SCIP_Longint SCIPsymhdlrGetNPresolCalls(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->npresolcalls;
}

/** gets number of variables fixed during presolving of symmetry handler */
int SCIPsymhdlrGetNFixedVars(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->nfixedvars;
}

/** gets number of variables aggregated during presolving of symmetry handler */
int SCIPsymhdlrGetNAggrVars(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->naggrvars;
}

/** gets number of variable types changed during presolving of symmetry handler */
int SCIPsymhdlrGetNChgVarTypes(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->nchgvartypes;
}

/** gets number of bounds changed during presolving of symmetry handler */
int SCIPsymhdlrGetNChgBds(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->nchgbds;
}

/** gets number of holes added to domains of variables during presolving of symmetry handler */
int SCIPsymhdlrGetNAddHoles(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->naddholes;
}

/** gets number of constraints deleted during presolving of symmetry handler */
int SCIPsymhdlrGetNDelConss(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->ndelconss;
}

/** gets number of constraints added during presolving of symmetry handler */
int SCIPsymhdlrGetNAddConss(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->naddconss;
}

/** gets number of coefficients changed during presolving of symmetry handler */
int SCIPsymhdlrGetNChgCoefs(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->nchgcoefs;
}

/** gets number of constraint sides changed during presolving of symmetry handler */
int SCIPsymhdlrGetNChgSides(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->nchgsides;
}

/** gets priority of separation method of symmetry handler */
int SCIPsymhdlrGetSepaPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->sepapriority;
}

/** gets priority of propagation method of symmetry handler */
int SCIPsymhdlrGetPropPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   )
{
   assert(symhdlr != NULL);

   return symhdlr->proppriority;
}

/** compresses information about permutations in symmetry information data structure */
SCIP_RETCODE SCIPsyminfoCompressPermInfo(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_SYMINFO*         syminfo,            /**< symmetry information */
   SCIP_VAR**            vars,               /**< vars present at time of symmetry computation */
   int                   nvars,              /**< number of vars present at time of symmetry computation */
   SCIP_Real             compressthreshold   /**< if percentage of moved vars is below threshold, compression is done */
   )
{
   SCIP_Real percentagemovedvars;
   int* labelmovedvars;
   int* labeltopermvaridx;
   int nmovedvars = 0;
   int i;
   int p;

   assert(scip != NULL);
   assert(syminfo != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(SCIPisGE(scip, compressthreshold, 0.0));
   assert(SCIPisLE(scip, compressthreshold, 1.0));
   assert(nvars == syminfo->npermvars);

   /* terminate if not enough variables are present */
   if( nvars < COMPRESSNVARSLB )
      return SCIP_OKAY;

   /* detect number of moved vars and label moved vars */
   SCIP_CALL( SCIPallocBufferArray(scip, &labelmovedvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &labeltopermvaridx, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      labelmovedvars[i] = -1;

      for( p = 0; p < syminfo->nperms; ++p )
      {
         if( syminfo->perms[p][i] != i )
         {
            labeltopermvaridx[nmovedvars] = i;
            labelmovedvars[i] = nmovedvars++;
            break;
         }
      }
   }

   /* check whether compression should be performed */
   percentagemovedvars = (SCIP_Real) nmovedvars / (SCIP_Real) nvars;

   if( nmovedvars <= 0 || SCIPisGT(scip, percentagemovedvars, compressthreshold) )
   {
      SCIPfreeBufferArray(scip, &labeltopermvaridx);
      SCIPfreeBufferArray(scip, &labelmovedvars);

      return SCIP_OKAY;
   }

   /* remove variables from permutations that are not affected by any permutation */
   for( p = 0; p < syminfo->nperms; ++p )
   {
      /* iterate over labels and adapt permutation (possibly taking signed permutations into account) */
      for( i = 0; i < nmovedvars; ++i )
      {
         assert(i <= labeltopermvaridx[i]);
         if( syminfo->perms[p][labeltopermvaridx[i]] >= nvars )
         {
            int imgvaridx;

            assert(syminfo->symtype == SYM_SYMTYPE_SIGNPERM);

            imgvaridx = syminfo->perms[p][labeltopermvaridx[i]] - nvars;
            syminfo->perms[p][i] = labelmovedvars[imgvaridx] + nmovedvars;

            assert(0 <= syminfo->perms[p][i] && syminfo->perms[p][i] < 2 * nmovedvars);
         }
         else
            syminfo->perms[p][i] = labelmovedvars[syminfo->perms[p][labeltopermvaridx[i]]];
      }

      if( syminfo->symtype == SYM_SYMTYPE_SIGNPERM )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &syminfo->perms[p], 2 * nvars, 2 * nmovedvars) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &syminfo->perms[p], nvars, nmovedvars) );
      }
   }

   /* remove variables from permvars array that are not affected by any symmetry */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &syminfo->permvars, syminfo->npermvars, nmovedvars) );
   for (i = 0; i < nmovedvars; ++i)
      syminfo->permvars[i] = vars[labeltopermvaridx[i]];
   syminfo->npermvars = nmovedvars;

   SCIPfreeBufferArray(scip, &labeltopermvaridx);
   SCIPfreeBufferArray(scip, &labelmovedvars);

   return SCIP_OKAY;
}

/** displays the cycle of a symmetry */
static
SCIP_RETCODE displayCycleOfSymmetry(
   SCIP*                 scip,               /**< SCIP pointer */
   int*                  perm,               /**< symmetry */
   SYM_SYMTYPE           symtype,            /**< type of symmetry */
   int                   baseidx,            /**< variable index for which cycle is computed */
   SCIP_Bool*            covered,            /**< allocated array to store covered variables */
   int                   nvars,              /**< number of (non-negated) variables in symmetry */
   SCIP_VAR**            vars                /**< variables on which symmetry acts */
   )
{
   char* string;
   int varidx;
   int j;

   assert(scip != NULL);
   assert(perm != NULL);
   assert(0 <= baseidx);
   assert((symtype == SYM_SYMTYPE_PERM && baseidx < nvars) ||
      (symtype == SYM_SYMTYPE_SIGNPERM && baseidx < 2 * nvars));
   assert(covered != NULL);

   /* skip fixed points or elements already covered in previous cycle */
   if( perm[baseidx] == baseidx || covered[baseidx] )
      return SCIP_OKAY;

   varidx = baseidx >= nvars ? baseidx - nvars : baseidx;
   string = (char*) SCIPvarGetName(vars[varidx]);
   SCIPinfoMessage(scip, NULL, "  (%s<%s>", baseidx >= nvars ? "negated " : "", string);
   j = perm[baseidx];
   covered[baseidx] = TRUE;
   while( j != baseidx )
   {
      covered[j] = TRUE;
      varidx = j >= nvars ? j - nvars : j;
      string = (char*) SCIPvarGetName(vars[varidx]);
      SCIPinfoMessage(scip, NULL, ",%s<%s>", j >= nvars ? "negated " : "", string);
      j = perm[j];
   }
   SCIPinfoMessage(scip, NULL, ")\n");

   return SCIP_OKAY;
}

/** displays symmetry information taking components into account */
static
SCIP_RETCODE displaySymmetriesWithComponents(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_SYMINFO*         syminfo             /**< propagator data */
   )
{
   SCIP_Bool* covered;
   SYM_SYMTYPE symtype;
   int* perm;
   int permlen;
   int npermvars;
   int i;
   int p;
   int c;

   assert(scip != NULL);
   assert(syminfo != NULL);
   assert(syminfo->nperms > 0);
   assert(syminfo->permvars != NULL);
   assert(syminfo->npermvars > 0);
   assert(syminfo->ncomponents > 0);

   symtype = syminfo->symtype;
   npermvars = syminfo->npermvars;
   permlen = symtype == SYM_SYMTYPE_PERM ? npermvars : 2 * npermvars;

   SCIP_CALL( SCIPallocClearBufferArray(scip, &covered, permlen) );

   if( symtype == SYM_SYMTYPE_PERM )
      SCIPinfoMessage(scip, NULL, "   Symmetries of different components are displayed as permutations.\n\n");
   else
      SCIPinfoMessage(scip, NULL, "   Symmetries of different components are displayed as signed permutations\n\n");
   for( c = 0; c < syminfo->ncomponents; ++c )
   {
      int cnt;

      SCIPinfoMessage(scip, NULL, "Display symmetries of component %d.\n", c);

      for( p = syminfo->componentbegins[c], cnt = 0; p < syminfo->componentbegins[c + 1]; ++p, ++cnt )
      {
         SCIPinfoMessage(scip, NULL, "Permutation %d:\n", cnt);
         perm = syminfo->perms[syminfo->components[p]];

         for( i = 0; i < permlen; ++i )
         {
            SCIP_CALL( displayCycleOfSymmetry(scip, perm, symtype, i, covered, npermvars, syminfo->permvars) );
         }

         for( i = 0; i < permlen; ++i )
            covered[i] = FALSE;
      }
   }

   SCIPfreeBufferArray(scip, &covered);

   return SCIP_OKAY;
}

/** displays generators of symmetry group, if available */
SCIP_RETCODE SCIPdisplaySymmetryGenerators(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMINFO*         syminfo             /**< symmetry information */
   )
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(syminfo != NULL);

   if( syminfo->nperms == -1 )
   {
      SCIPinfoMessage(scip, NULL, "Cannot display symmetries. Symmetries have not been computed yet.\n");
   }
   else if( syminfo->nperms == 0 )
   {
      SCIPinfoMessage(scip, NULL, "Cannot display symmetries. No symmetries detected.\n");
   }
   else
   {
      SCIP_CALL( displaySymmetriesWithComponents(scip, syminfo) );
   }

   return SCIP_OKAY;
}

/** adds a component to a symmetry handler */
SCIP_RETCODE SCIPaddSymhdlrComponent(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_SYMCOMP*         symcomp             /**< symmetry component */
   )
{
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(symhdlr != NULL);
   assert(symcomp != NULL);

   if( symhdlr->symcompssize == 0 )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &symhdlr->symcomps, 1) );
      symhdlr->symcompssize = 1;
   }
   else if( symhdlr->nsymcomps >= symhdlr->symcompssize )
   {
      int newlen;

      newlen = SCIPcalcMemGrowSize(set->scip, symhdlr->nsymcomps + 1);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &symhdlr->symcomps, symhdlr->symcompssize, newlen) );
      symhdlr->symcompssize = newlen;
   }
   symhdlr->symcomps[(symhdlr->nsymcomps)++] = symcomp;

   return SCIP_OKAY;
}

/** creates new operator node type (used for symmetry detection) and returns its representation
 *
 *  If the operator node already exists, the function terminates with SCIP_INVALIDDATA.
 */
SCIP_RETCODE SCIPcreateSymOpNodeType(
   SCIP*                 scip,               /**< SCIP pointer */
   const char*           opnodename,         /**< name of new operator node type */
   int*                  nodetype            /**< pointer to store the node type */
   )
{
   SCIP_SYMINFO* syminfo;

   assert(scip != NULL);
   assert(nodetype != NULL);

   syminfo = SCIPgetSyminfo(scip);
   if( syminfo == NULL )
   {
      SCIPerrorMessage("Cannot create operator node type, symmetry information has not been created yet.\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPhashmapExists(syminfo->customsymopnodetypes, (void*) opnodename) )
   {
      SCIPerrorMessage("Cannot create operator node type %s, it already exists.\n", opnodename);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPhashmapInsertInt(syminfo->customsymopnodetypes, (void*) opnodename, syminfo->nopnodetypes) );
   *nodetype = syminfo->nopnodetypes++;

   return SCIP_OKAY;
}

/** returns representation of an operator node type.
 *
 *  If the node type does not already exist, a new node type will be created.
 */
SCIP_RETCODE SCIPgetSymOpNodeType(
   SCIP*                 scip,               /**< SCIP pointer */
   const char*           opnodename,         /**< name of new operator node type */
   int*                  nodetype            /**< pointer to store the node type */
   )
{
   SCIP_SYMINFO* syminfo;

   assert(scip != NULL);

   syminfo = SCIPgetSyminfo(scip);
   if( syminfo == NULL )
   {
      SCIPerrorMessage("Cannot return operator node type, symmetry information has not been created yet.\n");
      return SCIP_INVALIDDATA;
   }

   if( !SCIPhashmapExists(syminfo->customsymopnodetypes, (void*) opnodename) )
   {
      SCIP_CALL( SCIPcreateSymOpNodeType(scip, opnodename, nodetype) );
   }
   else
      *nodetype = SCIPhashmapGetImageInt(syminfo->customsymopnodetypes, (void*) opnodename);

   return SCIP_OKAY;
}
