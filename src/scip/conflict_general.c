/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   conflict_general.c
 * @ingroup OTHER_CFILES
 * @brief  methods and datastructures for conflict analysis
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Jakob Witzig
 *
 * SCIP contains two kinds of conflict analysis:
 *    - In graph based conflict analysis, the graph consisting of derived
 *      is analysed. Code and documentation is available in conflict_graphanalysis.h
 *    - In dual proof analysis, an infeasible LP relaxation is analysed.
 *      Using the dual solution, a valid constraint is derived that is violated
 *      by all values in the domain. This constraint is added to the problem
 *      and can then be used for domain propagation.
 *      Code is available in conflict_dualproofanalysis.h
 * This file contains the methods that are shared by both kinds of conflict analysis.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "lpi/lpi.h"
#include "scip/clock.h"
#include "scip/conflict.h"
#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/cons_linear.h"
#include "scip/cuts.h"
#include "scip/history.h"
#include "scip/lp.h"
#include "scip/presolve.h"
#include "scip/prob.h"
#include "scip/prop.h"
#include "scip/pub_conflict.h"
#include "scip/pub_cons.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_paramset.h"
#include "scip/pub_prop.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_mem.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/struct_conflict.h"
#include "scip/struct_lp.h"
#include "scip/struct_prob.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_tree.h"
#include "scip/struct_var.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/visual.h"
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif

/* because calculations might cancel out some values, we stop the infeasibility analysis if a value is bigger than
 * 2^53 = 9007199254740992
 */
#define NUMSTOP 9007199254740992.0

/** returns the current number of conflict sets in the conflict set storage */
int SCIPconflictGetNConflicts(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nconflictsets;
}

/** returns the total number of conflict constraints that were added to the problem */
SCIP_Longint SCIPconflictGetNAppliedConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbconss + conflict->nappliedlocconss;
}

/** returns the total number of literals in conflict constraints that were added to the problem */
SCIP_Longint SCIPconflictGetNAppliedLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbliterals + conflict->nappliedlocliterals;
}

/** returns the total number of global bound changes applied by the conflict analysis */
SCIP_Longint SCIPconflictGetNGlobalChgBds(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nglbchgbds;
}

/** returns the total number of conflict constraints that were added globally to the problem */
SCIP_Longint SCIPconflictGetNAppliedGlobalConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbconss;
}

/** returns the total number of literals in conflict constraints that were added globally to the problem */
SCIP_Longint SCIPconflictGetNAppliedGlobalLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbliterals;
}

/** returns the total number of local bound changes applied by the conflict analysis */
SCIP_Longint SCIPconflictGetNLocalChgBds(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlocchgbds;
}

/** returns the total number of conflict constraints that were added locally to the problem */
SCIP_Longint SCIPconflictGetNAppliedLocalConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedlocconss;
}

/** returns the total number of literals in conflict constraints that were added locally to the problem */
SCIP_Longint SCIPconflictGetNAppliedLocalLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedlocliterals;
}

/** compares two conflict set entries, such that bound changes inferred later are
 *  ordered prior to ones that were inferred earlier
 */
static
SCIP_DECL_SORTPTRCOMP(conflictBdchginfoComp)
{  /*lint --e{715}*/
   SCIP_BDCHGINFO* bdchginfo1;
   SCIP_BDCHGINFO* bdchginfo2;

   bdchginfo1 = (SCIP_BDCHGINFO*)elem1;
   bdchginfo2 = (SCIP_BDCHGINFO*)elem2;
   assert(bdchginfo1 != NULL);
   assert(bdchginfo2 != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo1));
   assert(!SCIPbdchginfoIsRedundant(bdchginfo2));

   if( bdchginfo1 == bdchginfo2 )
      return 0;

   if( !SCIPbdchgidxIsEarlierNonNull(SCIPbdchginfoGetIdx(bdchginfo1), SCIPbdchginfoGetIdx(bdchginfo2)) )
      return -1;
   else
      return +1;
}

/** enables or disables all clocks of \p conflict, depending on the value of the flag */
void SCIPconflictEnableOrDisableClocks(
   SCIP_CONFLICT*        conflict,           /**< the conflict analysis data for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the conflict analysis data be enabled? */
   )
{
   assert(conflict != NULL);

   SCIPclockEnableOrDisable(conflict->boundlpanalyzetime, enable);
   SCIPclockEnableOrDisable(conflict->dIBclock, enable);
   SCIPclockEnableOrDisable(conflict->inflpanalyzetime, enable);
   SCIPclockEnableOrDisable(conflict->propanalyzetime, enable);
   SCIPclockEnableOrDisable(conflict->pseudoanalyzetime, enable);
   SCIPclockEnableOrDisable(conflict->sbanalyzetime, enable);
}

/** creates conflict analysis data for propagation conflicts */
SCIP_RETCODE SCIPconflictCreate(
   SCIP_CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflict != NULL);

   SCIP_ALLOC( BMSallocMemory(conflict) );

   SCIP_CALL( SCIPclockCreate(&(*conflict)->dIBclock, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->propanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->inflpanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->boundlpanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->sbanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->pseudoanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );

   /* enable or disable timing depending on the parameter statistic timing */
   SCIPconflictEnableOrDisableClocks((*conflict), set->time_statistictiming);

   SCIP_CALL( SCIPpqueueCreate(&(*conflict)->bdchgqueue, set->mem_arraygrowinit, set->mem_arraygrowfac,
         conflictBdchginfoComp, NULL) );
   SCIP_CALL( SCIPpqueueCreate(&(*conflict)->forcedbdchgqueue, set->mem_arraygrowinit, set->mem_arraygrowfac,
         conflictBdchginfoComp, NULL) );
   SCIP_CALL( SCIPconflictsetCreate(&(*conflict)->conflictset, blkmem) );
   (*conflict)->conflictsets = NULL;
   (*conflict)->conflictsetscores = NULL;
   (*conflict)->tmpbdchginfos = NULL;
   (*conflict)->conflictsetssize = 0;
   (*conflict)->nconflictsets = 0;
   (*conflict)->proofsets = NULL;
   (*conflict)->proofsetssize = 0;
   (*conflict)->nproofsets = 0;
   (*conflict)->tmpbdchginfossize = 0;
   (*conflict)->ntmpbdchginfos = 0;
   (*conflict)->count = 0;
   (*conflict)->nglbchgbds = 0;
   (*conflict)->nappliedglbconss = 0;
   (*conflict)->nappliedglbliterals = 0;
   (*conflict)->nlocchgbds = 0;
   (*conflict)->nappliedlocconss = 0;
   (*conflict)->nappliedlocliterals = 0;
   (*conflict)->npropcalls = 0;
   (*conflict)->npropsuccess = 0;
   (*conflict)->npropconfconss = 0;
   (*conflict)->npropconfliterals = 0;
   (*conflict)->npropreconvconss = 0;
   (*conflict)->npropreconvliterals = 0;
   (*conflict)->ninflpcalls = 0;
   (*conflict)->ninflpsuccess = 0;
   (*conflict)->ninflpconfconss = 0;
   (*conflict)->ninflpconfliterals = 0;
   (*conflict)->ninflpreconvconss = 0;
   (*conflict)->ninflpreconvliterals = 0;
   (*conflict)->ninflpiterations = 0;
   (*conflict)->nboundlpcalls = 0;
   (*conflict)->nboundlpsuccess = 0;
   (*conflict)->nboundlpconfconss = 0;
   (*conflict)->nboundlpconfliterals = 0;
   (*conflict)->nboundlpreconvconss = 0;
   (*conflict)->nboundlpreconvliterals = 0;
   (*conflict)->nboundlpiterations = 0;
   (*conflict)->nsbcalls = 0;
   (*conflict)->nsbsuccess = 0;
   (*conflict)->nsbconfconss = 0;
   (*conflict)->nsbconfliterals = 0;
   (*conflict)->nsbreconvconss = 0;
   (*conflict)->nsbreconvliterals = 0;
   (*conflict)->nsbiterations = 0;
   (*conflict)->npseudocalls = 0;
   (*conflict)->npseudosuccess = 0;
   (*conflict)->npseudoconfconss = 0;
   (*conflict)->npseudoconfliterals = 0;
   (*conflict)->npseudoreconvconss = 0;
   (*conflict)->npseudoreconvliterals = 0;
   (*conflict)->ndualproofsinfglobal = 0;
   (*conflict)->ndualproofsinflocal = 0;
   (*conflict)->ndualproofsinfsuccess = 0;
   (*conflict)->dualproofsinfnnonzeros = 0;
   (*conflict)->ndualproofsbndglobal = 0;
   (*conflict)->ndualproofsbndlocal = 0;
   (*conflict)->ndualproofsbndsuccess = 0;
   (*conflict)->dualproofsbndnnonzeros = 0;

   SCIP_CALL( SCIPconflictInitProofset((*conflict), blkmem) );

   return SCIP_OKAY;
}

/** frees conflict analysis data for propagation conflicts */
SCIP_RETCODE SCIPconflictFree(
   SCIP_CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(conflict != NULL);
   assert(*conflict != NULL);
   assert((*conflict)->nconflictsets == 0);
   assert((*conflict)->ntmpbdchginfos == 0);

#if defined(SCIP_CONFGRAPH) || defined(SCIP_CONFGRAPH_DOT)
   confgraphFree();
#endif

   SCIPclockFree(&(*conflict)->dIBclock);
   SCIPclockFree(&(*conflict)->propanalyzetime);
   SCIPclockFree(&(*conflict)->inflpanalyzetime);
   SCIPclockFree(&(*conflict)->boundlpanalyzetime);
   SCIPclockFree(&(*conflict)->sbanalyzetime);
   SCIPclockFree(&(*conflict)->pseudoanalyzetime);
   SCIPpqueueFree(&(*conflict)->bdchgqueue);
   SCIPpqueueFree(&(*conflict)->forcedbdchgqueue);
   SCIPconflictsetFree(&(*conflict)->conflictset, blkmem);
   SCIPproofsetFree(&(*conflict)->proofset, blkmem);

   BMSfreeMemoryArrayNull(&(*conflict)->conflictsets);
   BMSfreeMemoryArrayNull(&(*conflict)->conflictsetscores);
   BMSfreeMemoryArrayNull(&(*conflict)->proofsets);
   BMSfreeMemoryArrayNull(&(*conflict)->tmpbdchginfos);
   BMSfreeMemory(conflict);

   return SCIP_OKAY;
}

/** returns the conflict lower bound if the variable is present in the current conflict set; otherwise the global lower
 *  bound
 */
SCIP_Real SCIPconflictGetVarLb(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   if( var->conflictlbcount == conflict->count )
   {
      assert(EPSGE(var->conflictlb, var->conflictrelaxedlb, 1e-09));
      return var->conflictrelaxedlb;
   }

   return SCIPvarGetLbGlobal(var);
}

/** returns the conflict upper bound if the variable is present in the current conflict set; otherwise the global upper
 *  bound
 */
SCIP_Real SCIPconflictGetVarUb(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   if( var->conflictubcount == conflict->count )
   {
      assert(EPSLE(var->conflictub, var->conflictrelaxedub, 1e-09));
      return var->conflictrelaxedub;
   }

   return SCIPvarGetUbGlobal(var);
}

/** gets time in seconds used for preprocessing global conflict constraint before appliance */
SCIP_Real SCIPconflictGetGlobalApplTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->dIBclock);
}

/** gets time in seconds used for analyzing propagation conflicts */
SCIP_Real SCIPconflictGetPropTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->propanalyzetime);
}

/** gets number of calls to propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropcalls;
}

/** gets number of calls to propagation conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNPropSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropsuccess;
}

/** gets number of conflict constraints detected in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconfconss;
}

/** gets total number of literals in conflict constraints created in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconfliterals;
}

/** gets number of reconvergence constraints detected in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropreconvconss;
}

/** gets total number of literals in reconvergence constraints created in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropreconvliterals;
}

/** gets time in seconds used for analyzing infeasible LP conflicts */
SCIP_Real SCIPconflictGetInfeasibleLPTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->inflpanalyzetime);
}

/** gets number of calls to infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpcalls;
}

/** gets number of calls to infeasible LP conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNInfeasibleLPSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpsuccess;
}

/** gets number of conflict constraints detected in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpconfconss;
}

/** gets total number of literals in conflict constraints created in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpconfliterals;
}

/** gets number of reconvergence constraints detected in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpreconvconss;
}

/** gets total number of literals in reconvergence constraints created in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpreconvliterals;
}

/** gets number of LP iterations in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpiterations;
}

/** gets time in seconds used for analyzing bound exceeding LP conflicts */
SCIP_Real SCIPconflictGetBoundexceedingLPTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->boundlpanalyzetime);
}

/** gets number of calls to bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpcalls;
}

/** gets number of calls to bound exceeding LP conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNBoundexceedingLPSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpsuccess;
}

/** gets number of conflict constraints detected in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpconfconss;
}

/** gets total number of literals in conflict constraints created in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpconfliterals;
}

/** gets number of reconvergence constraints detected in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpreconvconss;
}

/** gets total number of literals in reconvergence constraints created in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpreconvliterals;
}

/** gets number of LP iterations in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpiterations;
}

/** gets time in seconds used for analyzing infeasible strong branching conflicts */
SCIP_Real SCIPconflictGetStrongbranchTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->sbanalyzetime);
}

/** gets number of successful calls to dual proof analysis derived from infeasible LPs */
SCIP_Longint SCIPconflictGetNDualproofsInfSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ndualproofsinfsuccess;
}

/** gets number of globally valid dual proof constraints derived from infeasible LPs */
SCIP_Longint SCIPconflictGetNDualproofsInfGlobal(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ndualproofsinfglobal;
}

/** gets number of locally valid dual proof constraints derived from infeasible LPs */
SCIP_Longint SCIPconflictGetNDualproofsInfLocal(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ndualproofsinflocal;
}

/** gets average length of dual proof constraints derived from infeasible LPs */
SCIP_Longint SCIPconflictGetNDualproofsInfNonzeros(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->dualproofsinfnnonzeros;
}

/** gets number of successfully analyzed dual proofs derived from bound exceeding LPs */
SCIP_Longint SCIPconflictGetNDualproofsBndSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ndualproofsbndsuccess;
}

/** gets number of globally applied dual proofs derived from bound exceeding LPs */
SCIP_Longint SCIPconflictGetNDualproofsBndGlobal(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ndualproofsbndglobal;
}

/** gets number of locally applied dual proofs derived from bound exceeding LPs */
SCIP_Longint SCIPconflictGetNDualproofsBndLocal(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ndualproofsbndlocal;
}

/** gets average length of dual proofs derived from bound exceeding LPs */
SCIP_Longint SCIPconflictGetNDualproofsBndNonzeros(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->dualproofsbndnnonzeros;
}

/** gets number of calls to infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbcalls;
}

/** gets number of calls to infeasible strong branching conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNStrongbranchSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbsuccess;
}

/** gets number of conflict constraints detected in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconfconss;
}

/** gets total number of literals in conflict constraints created in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconfliterals;
}

/** gets number of reconvergence constraints detected in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbreconvconss;
}

/** gets total number of literals in reconvergence constraints created in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbreconvliterals;
}

/** gets number of LP iterations in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbiterations;
}

/** adds a weighted LP row to an aggregation row */
static
SCIP_RETCODE addRowToAggrRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             weight,             /**< weight for scaling */
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(set != NULL);
   assert(row != NULL);
   assert(weight != 0.0);

   /* add minimal value to dual row's left hand side: y_i < 0 -> lhs, y_i > 0 -> rhs */
   if( weight < 0.0 )
   {
      assert(!SCIPsetIsInfinity(set, -row->lhs));
      SCIP_CALL( SCIPaggrRowAddRow(set->scip, aggrrow, row, weight, -1) );
   }
   else
   {
      assert(!SCIPsetIsInfinity(set, row->rhs));
      SCIP_CALL( SCIPaggrRowAddRow(set->scip, aggrrow, row, weight, +1) );
   }
   SCIPsetDebugMsg(set, " -> add %s row <%s>[%g,%g](lp depth: %d): dual=%g -> dualrhs=%g\n",
      row->local ? "local" : "global",
      SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant,
      row->lpdepth, weight, SCIPaggrRowGetRhs(aggrrow));

   return SCIP_OKAY;
}

/** checks validity of an LP row and a corresponding weight */
static
SCIP_Bool checkDualFeasibility(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             weight,             /**< weight for scaling */
   SCIP_Bool*            zerocontribution    /**< pointer to store whether every row entry is zero within tolerances */
   )
{
   SCIP_Bool valid = TRUE;

   *zerocontribution = TRUE;

   /* dual solution values of 0.0 are always valid */
   if( REALABS(weight) > QUAD_EPSILON )
   {
      *zerocontribution = FALSE;

      /* check dual feasibility */
      if( (SCIPsetIsInfinity(set, -row->lhs) && weight > 0.0) || (SCIPsetIsInfinity(set, row->rhs) && weight < 0.0) )
      {
         int i;

         /* ignore slight numerical violations if the contribution of every component of the row is close to zero */
         if( weight > 0.0 )
            *zerocontribution = SCIPsetIsDualfeasZero(set, row->rhs * weight);
         else
            *zerocontribution = SCIPsetIsDualfeasZero(set, row->lhs * weight);

         for( i = 0; i < row->len && *zerocontribution; i++ )
         {
            if( !SCIPsetIsDualfeasZero(set, weight * row->vals[i]) )
               *zerocontribution = FALSE;
         }

         if( !(*zerocontribution) )
         {
            SCIPsetDebugMsg(set, " -> invalid dual solution value %g for row <%s>: lhs=%g, rhs=%g\n",
               weight, SCIProwGetName(row), row->lhs, row->rhs);

            valid = FALSE;
         }
      }
   }

   return valid;
}

/** calculates the minimal activity of a given aggregation row */
SCIP_Real SCIPaggrRowGetMinActivity(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables (or NULL for global bounds) */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables (or NULL for global bounds) */
   SCIP_Bool*            infdelta            /**< pointer to store whether at least one variable contributes with an infinite value */
   )
{
   SCIP_VAR** vars;
   SCIP_Real QUAD(minact);
   int* inds;
   int nnz;
   int i;

   vars = SCIPprobGetVars(transprob);
   assert(vars != NULL);

   nnz = SCIPaggrRowGetNNz(aggrrow);
   inds = SCIPaggrRowGetInds(aggrrow);

   QUAD_ASSIGN(minact, 0.0);

   if( infdelta != NULL )
      *infdelta = FALSE;

   for( i = 0; i < nnz; i++ )
   {
      SCIP_Real val;
      SCIP_Real QUAD(delta);
      int v = inds[i];

      assert(SCIPvarGetProbindex(vars[v]) == v);

      val = SCIPaggrRowGetProbvarValue(aggrrow, v);

      if( val > 0.0 )
      {
         SCIP_Real bnd = (curvarlbs == NULL ? SCIPvarGetLbGlobal(vars[v]) : curvarlbs[v]);

         if( SCIPsetIsInfinity(set, -bnd) )
         {
            if( infdelta != NULL )
               *infdelta = TRUE;

            return -SCIPsetInfinity(set);
         }

         SCIPquadprecProdDD(delta, val, bnd);
      }
      else
      {
         SCIP_Real bnd = (curvarubs == NULL ? SCIPvarGetUbGlobal(vars[v]) : curvarubs[v]);

         if( SCIPsetIsInfinity(set, bnd) )
         {
            if( infdelta != NULL )
               *infdelta = TRUE;

            return -SCIPsetInfinity(set);
         }

         SCIPquadprecProdDD(delta, val, bnd);
      }

      /* update minimal activity */
      SCIPquadprecSumQQ(minact, minact, delta);
   }

   /* check whether the minimal activity is infinite */
   if( SCIPsetIsInfinity(set, QUAD_TO_DBL(minact)) )
      return SCIPsetInfinity(set);
   if( SCIPsetIsInfinity(set, -QUAD_TO_DBL(minact)) )
      return -SCIPsetInfinity(set);

   return QUAD_TO_DBL(minact);
}

/** sort local rows by increasing depth and number of nonzeros as tie-breaker */
static
SCIP_RETCODE sortLocalRows(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row */
   SCIP_ROW**            rows,               /**< array of local rows */
   int*                  rowinds,            /**< array of row indices */
   int*                  rowdepth,           /**< array of LP depths */
   int                   nrows               /**< number of local rows */
   )
{
   int* rownnz;
   int i;

   assert(aggrrow != NULL);
   assert(rows != NULL);
   assert(nrows > 0);
   assert(rowinds != NULL);
   assert(rowdepth != NULL);

   /* sort row indices by increasing depth */
   SCIPsortIntInt(rowdepth, rowinds, nrows);
   assert(rowdepth[0] <= rowdepth[nrows-1]);

   SCIP_CALL( SCIPsetAllocBufferArray(set, &rownnz, nrows) );

   /* get number of nonzero entries for every row */
   for( i = 0; i < nrows; i++ )
   {
      SCIP_ROW* row = rows[rowinds[i]];
      assert(row != NULL);

      rownnz[i] = row->len;
   }

   /* since SCIP has no stable sorting function we sort each bucket separately */
   for( i = 0; i < nrows; i++ )
   {
      int j = i;
      int d = rowdepth[i];

      /* search for the next row with a greater depth */
      while( j+1 < nrows && rowdepth[j+1] == d )
         j++;

      /* the bucket has size one */
      if( j == i )
         continue;

      assert(j-i+1 <= nrows);

      /* sort row indices by increasing number of nonzero elements */
      SCIPsortIntIntInt(&rownnz[i], &rowdepth[i], &rowinds[i], j-i+1);
      assert(rownnz[i] <= rownnz[j]);

      i = j;
   } /*lint --e{850} i is modified in the body of the for loop */

#ifndef NDEBUG
   for( i = 0; i < nrows-1; i++ )
      assert(rowdepth[i] < rowdepth[i+1] || (rowdepth[i] == rowdepth[i+1] && rownnz[i] <= rownnz[i+1]));
#endif

   SCIPsetFreeBufferArray(set, &rownnz);

   return SCIP_OKAY;
}

/** adds locally valid rows to the proof constraint */
static
SCIP_RETCODE addLocalRows(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_AGGRROW*         proofrow,           /**< aggregated row representing the proof */
   SCIP_ROW**            rows,               /**< array if locally valid rows */
   SCIP_Real*            dualsols,           /**< dual solution vector */
   int*                  localrowinds,       /**< array of row indecies */
   int*                  localrowdepth,      /**< array of row depths */
   int                   nlocalrows,         /**< number of local rows stored in rows array */
   SCIP_Real*            proofact,           /**< pointer to store the activity of the proof constraint */
   int*                  validdepth,         /**< pointer to store the depth where the proof constraint is valid */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_Bool*            valid               /**< pointer store whether the proof constraint is valid */
   )
{
   SCIP_Bool infdelta;
   int i;

   assert(set != NULL);
   assert(lp != NULL);

   *validdepth = 0;

   if( !set->conf_uselocalrows )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "add local rows to dual proof:\n");

   /* check whether the proof is already valid, e.g., violated within the local bounds */
   *proofact = SCIPaggrRowGetMinActivity(set, transprob, proofrow, curvarlbs, curvarubs, &infdelta);

   /* we stop if the minimal activity is infinite but all variables have a finite activity delta (bad numerics) */
   if( !infdelta && SCIPsetIsInfinity(set, REALABS(*proofact)) )
   {
      *valid = FALSE;
      return SCIP_OKAY;
   }

   /* break if the proof is valid w.r.t local bounds
    * note: it can happen that the proof contains a variable with an infinite activity delta.
    *       here, we don't break immediately because we might be able to fix it by adding local rows
    */
   if( !infdelta && SCIPsetIsGT(set, *proofact, SCIPaggrRowGetRhs(proofrow)) )
   {
      *valid = TRUE;
      return SCIP_OKAY;
   }

   /* sort local rows by depth */
   SCIP_CALL( sortLocalRows(set, proofrow, rows, localrowinds, localrowdepth, nlocalrows) );

   /* add successively local rows */
   for( i = 0; i < nlocalrows; ++i )
   {
      SCIP_ROW* row;
      int r;

      r = localrowinds[i];
      row = rows[r];

      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->vals != NULL);
      assert(row == lp->lpirows[r]);
      assert(row->local);
      assert(row->lpdepth == localrowdepth[i]);

      /* ignore dual solution values of 0.0 (in this case: y_i == 0) */
      if( REALABS(dualsols[r]) > 0.0 )
      {
#ifndef NDEBUG
         SCIP_Bool zerocontribution;

         /* check dual feasibility */
         *valid = checkDualFeasibility(set, row, dualsols[r], &zerocontribution);
         assert(*valid);
         assert(!zerocontribution);
#endif

         if( SCIPsetIsDualfeasZero(set, dualsols[r]) )
            continue;

         /* add row to dual proof */
         SCIP_CALL( addRowToAggrRow(set, row, -dualsols[r], proofrow) );

         /* update depth where the proof is valid */
         if( *validdepth < localrowdepth[i] )
            *validdepth = localrowdepth[i];

         /* get the new minimal activity */
         *proofact = SCIPaggrRowGetMinActivity(set, transprob, proofrow, curvarlbs, curvarubs, &infdelta);

         /* we stop if the minimal activity is infinite but all variables have a finite activity delta (bad numerics) */
         if( !infdelta && SCIPsetIsInfinity(set, REALABS(*proofact)) )
         {
            *valid = FALSE;
            goto TERMINATE;
         }

         /* break if the proof is valid w.r.t local bounds */
         if( !infdelta && SCIPsetIsGT(set, *proofact, SCIPaggrRowGetRhs(proofrow)) )
         {
            *valid = TRUE;
            break;
         }
      }
   }

   /* remove all nearly zero coefficients */
   SCIPaggrRowRemoveZeros(set->scip, proofrow, TRUE, valid);

  TERMINATE:
   if( !(*valid) )
   {
      SCIPsetDebugMsg(set, " -> proof is not valid: %g <= %g\n", *proofact, SCIPaggrRowGetRhs(proofrow));
      SCIPsetDebugMsg(set, " -> stop due to numerical troubles\n");
   }
   else
   {
      *proofact = SCIPaggrRowGetMinActivity(set, transprob, proofrow, curvarlbs, curvarubs, &infdelta);

      /* we stop if the minimal activity is infinite but all variables have a finite activity delta (bad numerics) */
      if( !infdelta && SCIPsetIsInfinity(set, REALABS(*proofact)) )
      {
         *valid = FALSE;
         SCIPsetDebugMsg(set, " -> proof is not valid: %g <= %g [infdelta: %u]\n", *proofact, SCIPaggrRowGetRhs(proofrow), infdelta);
      }
      else if( infdelta || SCIPsetIsLE(set, *proofact, SCIPaggrRowGetRhs(proofrow)) )
      {
         *valid = FALSE;
         SCIPsetDebugMsg(set, " -> proof is not valid: %g <= %g [infdelta: %u]\n", *proofact, SCIPaggrRowGetRhs(proofrow), infdelta);
      }
   }

   return SCIP_OKAY;
}

/** calculates a Farkas proof from the current dual LP solution */
SCIP_RETCODE SCIPgetFarkasProof(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPI*             lpi,                /**< LPI data */
   SCIP_TREE*            tree,               /**< tree data */
   SCIP_AGGRROW*         farkasrow,          /**< aggregated row representing the proof */
   SCIP_Real*            farkasact,          /**< maximal activity of the proof constraint */
   int*                  validdepth,         /**< pointer to store the valid depth of the proof constraint */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_Bool*            valid               /**< pointer store whether the proof constraint is valid */
   )
{
   SCIP_ROW** rows;
   SCIP_Real* dualfarkas;
   SCIP_ROW* row;
   int* localrowinds;
   int* localrowdepth;
   SCIP_Bool infdelta;
   int nlocalrows;
   int nrows;
   int r;

   assert(set != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(valid != NULL);

   assert(SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsDualFeasible(lpi));
   assert(SCIPlpiIsPrimalInfeasible(lpi) || !SCIPlpDivingObjChanged(lp));

   /* get LP rows and problem variables */
   rows = SCIPlpGetRows(lp);
   nrows = SCIPlpGetNRows(lp);
   assert(nrows == 0 || rows != NULL);
   assert(nrows == lp->nlpirows);

   /* it can happen that infeasibility is detetected within LP presolve. in that case, the LP solver may not be able to
    * to return the dual ray.
    */
   if( !SCIPlpiHasDualRay(lpi) )
   {
      *valid = FALSE;
      return SCIP_OKAY;
   }

   assert(farkasrow != NULL);

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualfarkas, nrows) );
   BMSclearMemoryArray(dualfarkas, nrows);

   /* get dual Farkas values of rows */
   SCIP_CALL( SCIPlpiGetDualfarkas(lpi, dualfarkas) );

   localrowinds = NULL;
   localrowdepth = NULL;
   nlocalrows = 0;

   /* calculate the Farkas row */
   (*valid) = TRUE;
   (*validdepth) = 0;

   for( r = 0; r < nrows; ++r )
   {
      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->vals != NULL);
      assert(row == lp->lpirows[r]);

      /* ignore dual ray values of 0.0 (in this case: y_i == z_i == 0) */
      if( REALABS(dualfarkas[r]) > 0.0 )
      {
         SCIP_Bool zerocontribution;

         /* check dual feasibility */
         *valid = checkDualFeasibility(set, row, dualfarkas[r], &zerocontribution);

         if( !(*valid) )
            goto TERMINATE;

         if( zerocontribution )
            continue;

         if( SCIPsetIsDualfeasZero(set, dualfarkas[r]) )
            continue;

         if( !row->local )
         {
            SCIP_CALL( addRowToAggrRow(set, row, -dualfarkas[r], farkasrow) );

            /* due to numerical reasons we want to stop */
            if( REALABS(SCIPaggrRowGetRhs(farkasrow)) > NUMSTOP )
            {
               (*valid) = FALSE;
               goto TERMINATE;
            }
         }
         else
         {
            int lpdepth = SCIProwGetLPDepth(row);

            if( nlocalrows == 0 && lpdepth < SCIPtreeGetFocusDepth(tree) )
            {
               SCIP_CALL( SCIPsetAllocBufferArray(set, &localrowinds, nrows-r) );
               SCIP_CALL( SCIPsetAllocBufferArray(set, &localrowdepth, nrows-r) );
            }

            if( lpdepth < SCIPtreeGetFocusDepth(tree) )
            {
               assert(localrowinds != NULL);
               assert(localrowdepth != NULL);

               localrowinds[nlocalrows] = r;
               localrowdepth[nlocalrows++] = lpdepth;
            }
         }
      }
   }

   /* remove all coefficients that are too close to zero */
   SCIPaggrRowRemoveZeros(set->scip, farkasrow, TRUE, valid);

   if( !(*valid) )
      goto TERMINATE;

   infdelta = FALSE;

   /* calculate the current Farkas activity, always using the best bound w.r.t. the Farkas coefficient */
   *farkasact = SCIPaggrRowGetMinActivity(set, prob, farkasrow, curvarlbs, curvarubs, &infdelta);

   SCIPsetDebugMsg(set, " -> farkasact=%g farkasrhs=%g [infdelta: %u], \n",
      (*farkasact), SCIPaggrRowGetRhs(farkasrow), infdelta);

   /* The constructed proof is not valid, this can happen due to numerical reasons,
    * e.g., we only consider rows r with !SCIPsetIsZero(set, dualfarkas[r]),
    * or because of local rows were ignored so far.
    * Due to the latter case, it might happen at least one variable contributes
    * with an infinite value to the activity (see: https://git.zib.de/integer/scip/issues/2743)
    */
   if( infdelta || SCIPsetIsFeasLE(set, *farkasact, SCIPaggrRowGetRhs(farkasrow)))
   {
      /* add contribution of local rows */
      if( nlocalrows > 0 && set->conf_uselocalrows > 0 )
      {
         SCIP_CALL( addLocalRows(set, prob, lp, farkasrow, rows, dualfarkas, localrowinds, localrowdepth,
               nlocalrows, farkasact, validdepth, curvarlbs, curvarubs, valid) );
      }
      else
      {
         (*valid) = FALSE;
         SCIPsetDebugMsg(set, " -> proof is not valid to due infinite activity delta\n");
      }
   }

  TERMINATE:

   SCIPfreeBufferArrayNull(set->scip, &localrowdepth);
   SCIPfreeBufferArrayNull(set->scip, &localrowinds);
   SCIPsetFreeBufferArray(set, &dualfarkas);

   return SCIP_OKAY;
}

/** calculates a dual proof from the current dual LP solution */
SCIP_RETCODE SCIPgetDualProof(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPI*             lpi,                /**< LPI data */
   SCIP_TREE*            tree,               /**< tree data */
   SCIP_AGGRROW*         farkasrow,          /**< aggregated row representing the proof */
   SCIP_Real*            farkasact,          /**< maximal activity of the proof constraint */
   int*                  validdepth,         /**< pointer to store the valid depth of the proof constraint */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_Bool*            valid               /**< pointer store whether the proof constraint is valid */
   )
{
   SCIP_RETCODE retcode;
   SCIP_ROW** rows;
   SCIP_ROW* row;
   SCIP_Real* primsols;
   SCIP_Real* dualsols;
   SCIP_Real* redcosts;
   int* localrowinds;
   int* localrowdepth;
   SCIP_Real maxabsdualsol;
   SCIP_Bool infdelta;
   int nlocalrows;
   int nrows;
   int ncols;
   int r;

   assert(set != NULL);
   assert(transprob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(valid != NULL);

   *validdepth = 0;
   *valid = TRUE;

   localrowinds = NULL;
   localrowdepth = NULL;
   nlocalrows = 0;

   /* get LP rows and problem variables */
   rows = SCIPlpGetRows(lp);
   nrows = SCIPlpGetNRows(lp);
   ncols = SCIPlpGetNCols(lp);
   assert(nrows == 0 || rows != NULL);
   assert(nrows == lp->nlpirows);

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &primsols, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualsols, nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &redcosts, ncols) );

   /* get solution from LPI */
   retcode = SCIPlpiGetSol(lpi, NULL, primsols, dualsols, NULL, redcosts);
   if( retcode == SCIP_LPERROR ) /* on an error in the LP solver, just abort the conflict analysis */
   {
      (*valid) = FALSE;
      goto TERMINATE;
   }
   SCIP_CALL( retcode );
#ifdef SCIP_DEBUG
   {
      SCIP_Real objval;
      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      SCIPsetDebugMsg(set, " -> LP objval: %g\n", objval);
   }
#endif

   /* check whether the dual solution is numerically stable */
   maxabsdualsol = 0;
   for( r = 0; r < nrows; r++ )
   {
      SCIP_Real absdualsol = REALABS(dualsols[r]);

      if( absdualsol > maxabsdualsol )
         maxabsdualsol = absdualsol;
   }

   /* don't consider dual solution with maxabsdualsol > 1e+07, this would almost cancel out the objective constraint */
   if( maxabsdualsol > 1e+07 )
   {
      (*valid) = FALSE;
      goto TERMINATE;
   }

   /* clear the proof */
   SCIPaggrRowClear(farkasrow);

   /* Let y be the dual solution and r be the reduced cost vector. Let z be defined as
    *    z_i := y_i if i is a global row,
    *    z_i := 0   if i is a local row.
    * Define the set X := {x | lhs <= Ax <= rhs, lb <= x <= ub, c^Tx <= c*}, with c* being the current primal bound.
    * Then the following inequalities are valid for all x \in X:
    *                                 - c* <= -c^Tx
    *   <=>                     z^TAx - c* <= (z^TA - c^T) x
    *   <=>                     z^TAx - c* <= (y^TA - c^T - (y-z)^TA) x
    *   <=>                     z^TAx - c* <= (-r^T - (y-z)^TA) x         (dual feasibility of (y,r): y^TA + r^T == c^T)
    * Because lhs <= Ax <= rhs and lb <= x <= ub, the inequality can be relaxed to give
    *     min{z^Tq | lhs <= q <= rhs} - c* <= max{(-r^T - (y-z)^TA) x | lb <= x <= ub}, or X = {}.
    *
    * The resulting dual row is:  z^T{lhs,rhs} - c* <= (-r^T - (y-z)^TA){lb,ub},
    * where lhs, rhs, lb, and ub are selected in order to maximize the feasibility of the row.
    */

   /* add the objective function to the aggregation row with respect to the current cutoff bound
    *
    * for an integral objective the right-hand side is reduced by the cutoff bound delta to cut off up to the next
    * possible objective value below the cutoff bound
    */
   SCIP_CALL( SCIPaggrRowAddObjectiveFunction(set->scip, farkasrow, lp->cutoffbound - (SCIPprobIsObjIntegral(transprob) ? SCIPsetCutoffbounddelta(set) : 0.0), 1.0) );

   /* dual row: z^T{lhs,rhs} - c* <= (-r^T - (y-z)^TA){lb,ub}
    * process rows: add z^T{lhs,rhs} to the dual row's left hand side, and -(y-z)^TA to the dual row's coefficients
    */
   for( r = 0; r < nrows; ++r )
   {
      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->vals != NULL);
      assert(row == lp->lpirows[r]);

      /* ignore dual solution values of 0.0 (in this case: y_i == z_i == 0) */
      if( REALABS(dualsols[r]) > 0.0 )
      {
         SCIP_Bool zerocontribution;

         /* check dual feasibility */
         *valid = checkDualFeasibility(set, row, dualsols[r], &zerocontribution);

         if( !(*valid) )
            goto TERMINATE;

         if( zerocontribution )
            continue;

         if( SCIPsetIsDualfeasZero(set, dualsols[r]) )
            continue;

         /* skip local row */
         if( !row->local )
         {
            SCIP_CALL( addRowToAggrRow(set, row, -dualsols[r], farkasrow) );

            /* due to numerical reasons we want to stop */
            if( REALABS(SCIPaggrRowGetRhs(farkasrow)) > NUMSTOP )
            {
               (*valid) = FALSE;
               goto TERMINATE;
            }
         }
         else
         {
            int lpdepth = SCIProwGetLPDepth(row);

            if( nlocalrows == 0 && lpdepth < SCIPtreeGetFocusDepth(tree) )
            {
               SCIP_CALL( SCIPsetAllocBufferArray(set, &localrowinds, nrows-r) );
               SCIP_CALL( SCIPsetAllocBufferArray(set, &localrowdepth, nrows-r) );
            }

            if( lpdepth < SCIPtreeGetFocusDepth(tree) )
            {
               assert(localrowinds != NULL);
               assert(localrowdepth != NULL);

               localrowinds[nlocalrows] = r;
               localrowdepth[nlocalrows++] = lpdepth;
            }
         }
      }
   }

   /* remove all nearly zero coefficients */
   SCIPaggrRowRemoveZeros(set->scip, farkasrow, TRUE, valid);

   if( !(*valid) )
      goto TERMINATE;

   infdelta = FALSE;

   /* check validity of the proof */
   *farkasact = SCIPaggrRowGetMinActivity(set, transprob, farkasrow, curvarlbs, curvarubs, &infdelta);

   SCIPsetDebugMsg(set, " -> farkasact=%g farkasrhs=%g [infdelta: %u], \n",
      (*farkasact), SCIPaggrRowGetRhs(farkasrow), infdelta);

   /* The constructed proof is not valid, this can happen due to numerical reasons,
    * e.g., we only consider rows r with !SCIPsetIsZero(set, dualsol[r]),
    * or because of local rows were ignored so far.
    * Due to the latter case, it might happen at least one variable contributes
    * with an infinite value to the activity (see: https://git.zib.de/integer/scip/issues/2743)
    */
   if( infdelta || SCIPsetIsFeasLE(set, *farkasact, SCIPaggrRowGetRhs(farkasrow)))
   {
      /* add contribution of local rows */
      if( nlocalrows > 0 && set->conf_uselocalrows > 0 )
      {
         SCIP_CALL( addLocalRows(set, transprob, lp, farkasrow, rows, dualsols, localrowinds, localrowdepth,
               nlocalrows, farkasact, validdepth, curvarlbs, curvarubs, valid) );
      }
      else
      {
         (*valid) = FALSE;
         SCIPsetDebugMsg(set, " -> proof is not valid to due infinite activity delta\n");
      }
   }

  TERMINATE:

   SCIPfreeBufferArrayNull(set->scip, &localrowdepth);
   SCIPfreeBufferArrayNull(set->scip, &localrowinds);
   SCIPsetFreeBufferArray(set, &redcosts);
   SCIPsetFreeBufferArray(set, &dualsols);
   SCIPsetFreeBufferArray(set, &primsols);

   return SCIP_OKAY;
}


/*
 * pseudo solution conflict analysis
 */

/** analyzes a pseudo solution with objective value exceeding the current cutoff to find out the bound changes on
 *  variables that were responsible for the objective value degradation;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for pseudo solution conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyzePseudo(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real* curvarlbs;
   SCIP_Real* curvarubs;
   int* lbchginfoposs;
   int* ubchginfoposs;
   SCIP_Real* pseudocoefs;
   SCIP_Real pseudolhs;
   SCIP_Real pseudoact;
   int nvars;
   int v;

   assert(conflict != NULL);
   assert(conflict->nconflictsets == 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(lp != NULL);
   assert(!SCIPsetIsInfinity(set, -SCIPlpGetPseudoObjval(lp, set, transprob)));
   assert(!SCIPsetIsInfinity(set, lp->cutoffbound));

   if( success != NULL )
      *success = FALSE;

   /* check, if pseudo solution conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_usepseudo )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "analyzing pseudo solution (obj: %g) that exceeds objective limit (%g)\n",
      SCIPlpGetPseudoObjval(lp, set, transprob), lp->cutoffbound);

   conflict->conflictset->conflicttype = SCIP_CONFTYPE_BNDEXCEEDING;
   conflict->conflictset->usescutoffbound = TRUE;

   /* start timing */
   SCIPclockStart(conflict->pseudoanalyzetime, set);
   conflict->npseudocalls++;

   vars = transprob->vars;
   nvars = transprob->nvars;
   assert(nvars == 0 || vars != NULL);

   /* The current primal bound c* gives an upper bound for the current pseudo objective value:
    *   min{c^T x | lb <= x <= ub} <= c*.
    * We have to transform this row into a >= inequality in order to use methods above:
    *                          -c* <= max{-c^T x | lb <= x <= ub}.
    * In the local subproblem, this row is violated. We want to undo bound changes while still keeping the
    * row violated.
    */

   /* get temporary memory for remembering variables' current bounds and corresponding bound change information
    * positions in variable's bound change information arrays
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lbchginfoposs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ubchginfoposs, nvars) );

   /* get temporary memory for infeasibility proof coefficients */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &pseudocoefs, nvars) );

   /* for an integral objective use the cutoff bound reduced by the cutoff bound delta to cut off up to the next better
    * objective value
    */
   pseudolhs = -(lp->cutoffbound - (SCIPprobIsObjIntegral(transprob) ? SCIPsetCutoffbounddelta(set) : 0.0));

   /* store the objective values as infeasibility proof coefficients, and recalculate the pseudo activity */
   pseudoact = 0.0;
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      pseudocoefs[v] = -SCIPvarGetObj(var);
      curvarlbs[v] = SCIPvarGetLbLocal(var);
      curvarubs[v] = SCIPvarGetUbLocal(var);
      lbchginfoposs[v] = var->nlbchginfos-1;
      ubchginfoposs[v] = var->nubchginfos-1;

      if( SCIPsetIsZero(set, pseudocoefs[v]) )
      {
         pseudocoefs[v] = 0.0;
         continue;
      }

      if( pseudocoefs[v] > 0.0 )
         pseudoact += pseudocoefs[v] * curvarubs[v];
      else
         pseudoact += pseudocoefs[v] * curvarlbs[v];
   }
   assert(SCIPsetIsFeasEQ(set, pseudoact, -SCIPlpGetPseudoObjval(lp, set, transprob)));
   SCIPsetDebugMsg(set, "  -> recalculated pseudo infeasibility proof:  %g <= %g\n", pseudolhs, pseudoact);

   /* check, if the pseudo row is still violated (after recalculation of pseudo activity) */
   if( SCIPsetIsFeasGT(set, pseudolhs, pseudoact) )
   {
      int nconss;
      int nliterals;
      int nreconvconss;
      int nreconvliterals;

      /* undo bound changes without destroying the infeasibility proof */
      SCIP_CALL( SCIPundoBdchgsProof(set, transprob, SCIPtreeGetCurrentDepth(tree), pseudocoefs, pseudolhs, &pseudoact,
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs, NULL, NULL, NULL, lp->lpi) );

      /* analyze conflict on remaining bound changes */
      SCIP_CALL( SCIPconflictAnalyzeRemainingBdchgs(conflict, blkmem, set, stat, transprob, tree, FALSE, \
            lbchginfoposs, ubchginfoposs, &nconss, &nliterals, &nreconvconss, &nreconvliterals) );
      conflict->npseudosuccess += (nconss > 0 ? 1 : 0);
      conflict->npseudoconfconss += nconss;
      conflict->npseudoconfliterals += nliterals;
      conflict->npseudoreconvconss += nreconvconss;
      conflict->npseudoreconvliterals += nreconvliterals;
      if( success != NULL )
         *success = (nconss > 0);
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &pseudocoefs);
   SCIPsetFreeBufferArray(set, &ubchginfoposs);
   SCIPsetFreeBufferArray(set, &lbchginfoposs);
   SCIPsetFreeBufferArray(set, &curvarubs);
   SCIPsetFreeBufferArray(set, &curvarlbs);

   /* flush conflict set storage */
   SCIP_CALL( SCIPconflictFlushConss(conflict, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cliquetable) );

   /* stop timing */
   SCIPclockStop(conflict->pseudoanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing pseudo solution conflicts */
SCIP_Real SCIPconflictGetPseudoTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->pseudoanalyzetime);
}

/** gets number of calls to pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudocalls;
}

/** gets number of calls to pseudo solution conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNPseudoSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudosuccess;
}

/** gets number of conflict constraints detected in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconfconss;
}

/** gets total number of literals in conflict constraints created in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconfliterals;
}

/** gets number of reconvergence constraints detected in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoreconvconss;
}

/** gets total number of literals in reconvergence constraints created in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoreconvliterals;
}

/** actually performs analysis of infeasible LP */
static
SCIP_RETCODE conflictAnalyzeLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   SCIP_Bool*            dualproofsuccess,   /**< pointer to store success result of dual proof analysis */
   int*                  iterations,         /**< pointer to store the total number of LP iterations used */
   int*                  nconss,             /**< pointer to store the number of generated conflict constraints */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict constraints */
   int*                  nreconvconss,       /**< pointer to store the number of generated reconvergence constraints */
   int*                  nreconvliterals,    /**< pointer to store the number of literals generated reconvergence constraints */
   SCIP_Bool             marklpunsolved      /**< whether LP should be marked unsolved after analysis (needed for strong branching) */
   )
{
   SCIP_VAR** vars;
   SCIP_AGGRROW* farkasrow;
   SCIP_LPI* lpi;
   SCIP_Bool valid;
   SCIP_Bool globalinfeasible;
   int* lbchginfoposs;
   int* ubchginfoposs;
   int validdepth;
   int nvars;
   int v;
   SCIP_Real* curvarlbs;
   SCIP_Real* curvarubs;
   SCIP_Real farkasactivity;

   assert(conflict != NULL);
   assert(conflict->nconflictsets == 0);
   assert(set != NULL);
   assert(SCIPprobAllColsInLP(transprob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(iterations != NULL);
   assert(nconss != NULL);
   assert(nliterals != NULL);
   assert(nreconvconss != NULL);
   assert(nreconvliterals != NULL);

   *iterations = 0;
   *nconss = 0;
   *nliterals = 0;
   *nreconvconss = 0;
   *nreconvliterals = 0;

   vars = transprob->vars;
   nvars = transprob->nvars;

   valid = TRUE;
   validdepth = 0;

   /* get LP solver interface */
   lpi = SCIPlpGetLPI(lp);
   assert(SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsDualFeasible(lpi));
   assert(SCIPlpiIsPrimalInfeasible(lpi) || !SCIPlpDivingObjChanged(lp));

   if( !SCIPlpiIsPrimalInfeasible(lpi) )
   {
      SCIP_Real objval;

      assert(!SCIPlpDivingObjChanged(lp));

      /* make sure, a dual feasible solution exists, that exceeds the objective limit;
       * With FASTMIP setting, CPLEX does not apply the final pivot to reach the dual solution exceeding the objective
       * limit. Therefore, we have to either turn off FASTMIP and resolve the problem or continue solving it without
       * objective limit for at least one iteration. It seems that the strategy to continue with FASTMIP for one
       * additional simplex iteration yields better results.
       */
      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      if( objval < lp->lpiobjlim )
      {
         SCIP_RETCODE retcode;

         /* temporarily disable objective limit and install an iteration limit */
         SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_OBJLIM, SCIPlpiInfinity(lpi)) );
         SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, 1) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve LP */
         retcode = SCIPlpiSolveDual(lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         valid = (retcode != SCIP_LPERROR);
         if( valid )
         {
            int iter;

            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lpi, &iter) );
            (*iterations) += iter;
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            SCIPsetDebugMsg(set, " -> resolved objlim exceeding LP in %d iterations (total: %" SCIP_LONGINT_FORMAT ") (infeasible:%u, objlim: %u, optimal:%u)\n",
               iter, stat->nconflictlpiterations, SCIPlpiIsPrimalInfeasible(lpi), SCIPlpiIsObjlimExc(lpi), SCIPlpiIsOptimal(lpi));
            valid = (SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsDualFeasible(lpi));
         }

         /* reinstall old objective and iteration limits in LP solver */
         SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_OBJLIM, lp->lpiobjlim) );
         SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, lp->lpiitlim) );

         /* abort, if the LP produced an error */
         if( !valid )
            return SCIP_OKAY;
      }
   }
   assert(SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsDualFeasible(lpi));

   if( !SCIPlpiIsPrimalInfeasible(lpi) )
   {
      SCIP_Real objval;

      assert(!SCIPlpDivingObjChanged(lp));

      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      if( objval < lp->lpiobjlim )
      {
         SCIPsetDebugMsg(set, " -> LP does not exceed the cutoff bound: obj=%g, cutoff=%g\n", objval, lp->lpiobjlim);
         return SCIP_OKAY;
      }
      else
      {
         SCIPsetDebugMsg(set, " -> LP exceeds the cutoff bound: obj=%g, cutoff=%g\n", objval, lp->lpiobjlim);
      }
   }

   assert(valid);

   SCIP_CALL( SCIPaggrRowCreate(set->scip, &farkasrow) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lbchginfoposs, transprob->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ubchginfoposs, transprob->nvars) );

   farkasactivity = 0.0;

   /* get temporary memory for remembering variables' current bounds and corresponding bound change information
    * positions in variable's bound change information arrays
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );

   /* get current bounds and current positions in lb/ubchginfos arrays of variables */
   valid = TRUE;
   for( v = 0; v < nvars && valid; ++v )
   {
      SCIP_VAR* var;

      var = vars[v];

      curvarlbs[v] = SCIPvarGetLbLP(var, set);
      curvarubs[v] = SCIPvarGetUbLP(var, set);
      lbchginfoposs[v] = var->nlbchginfos-1;
      ubchginfoposs[v] = var->nubchginfos-1;
      assert(diving || SCIPsetIsEQ(set, curvarlbs[v], SCIPvarGetLbLocal(var)));
      assert(diving || SCIPsetIsEQ(set, curvarubs[v], SCIPvarGetUbLocal(var)));

      /* check, if last bound changes were due to strong branching or diving */
      if( diving )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
         if( SCIPsetIsGT(set, curvarlbs[v], lb) )
            lbchginfoposs[v] = var->nlbchginfos;
         else if( SCIPsetIsLT(set, curvarlbs[v], lb) )
         {
            /* the bound in the diving LP was relaxed -> the LP is not a subproblem of the current node -> abort! */
            /**@todo we could still analyze such a conflict, but we would have to take care with our data structures */
            valid = FALSE;
         }
         if( SCIPsetIsLT(set, curvarubs[v], ub) )
            ubchginfoposs[v] = var->nubchginfos;
         else if( SCIPsetIsGT(set, curvarubs[v], ub) )
         {
            /* the bound in the diving LP was relaxed -> the LP is not a subproblem of the current node -> abort! */
            /**@todo we could still analyze such a conflict, but we would have to take care with our data structures */
            valid = FALSE;
         }
      }
   }

   if( !valid )
      goto TERMINATE;

   /* the LP is prooven to be infeasible */
   if( SCIPlpiIsPrimalInfeasible(lpi) )
   {
      SCIP_CALL( SCIPgetFarkasProof(set, transprob, lp, lpi, tree, farkasrow, &farkasactivity, &validdepth,
         curvarlbs, curvarubs, &valid) );
   }
   /* the LP is dual feasible and/or exceeds the current incumbant solution */
   else
   {
      assert(SCIPlpiIsDualFeasible(lpi) || SCIPlpiIsObjlimExc(lpi));
      SCIP_CALL( SCIPgetDualProof(set, transprob, lp, lpi, tree, farkasrow, &farkasactivity, &validdepth,
         curvarlbs, curvarubs, &valid) );
   }

   if( !valid || validdepth >= SCIPtreeGetCurrentDepth(tree) )
      goto TERMINATE;

   globalinfeasible = FALSE;

   /* start dual proof analysis */
   if( ((set->conf_useinflp == 'b' || set->conf_useinflp == 'd') && conflict->conflictset->conflicttype == SCIP_CONFTYPE_INFEASLP)
      || ((set->conf_useboundlp == 'b' || set->conf_useboundlp == 'd') && conflict->conflictset->conflicttype == SCIP_CONFTYPE_BNDEXCEEDING) )
   {
      /* start dual proof analysis */
      SCIP_CALL( SCIPconflictAnalyzeDualProof(conflict, set, stat, blkmem, origprob, transprob, tree, reopt, lp, farkasrow, \
         validdepth, curvarlbs, curvarubs, TRUE, &globalinfeasible, dualproofsuccess) );
   }

   assert(valid);

   /* todo: in theory, we could apply conflict graph analysis for locally valid proofs, too, but this needs to be implemented */
   if( !globalinfeasible && validdepth <= SCIPtreeGetEffectiveRootDepth(tree)
      && (((set->conf_useinflp == 'b' || set->conf_useinflp == 'c') && conflict->conflictset->conflicttype == SCIP_CONFTYPE_INFEASLP)
      || ((set->conf_useboundlp == 'b' || set->conf_useboundlp == 'c') && conflict->conflictset->conflicttype == SCIP_CONFTYPE_BNDEXCEEDING)) )
   {
      SCIP_Real* farkascoefs;
      SCIP_Real farkaslhs;
      int* inds;
      int nnz;

#ifdef SCIP_DEBUG
      {
         SCIP_Real objlim;
         SCIPsetDebugMsg(set, "analyzing conflict on infeasible LP (infeasible: %u, objlimexc: %u, optimal:%u) in depth %d (diving: %u)\n",
               SCIPlpiIsPrimalInfeasible(lpi), SCIPlpiIsObjlimExc(lpi), SCIPlpiIsOptimal(lpi), SCIPtreeGetCurrentDepth(tree), diving);

         SCIP_CALL( SCIPlpiGetRealpar(lpi, SCIP_LPPAR_OBJLIM, &objlim) );
         SCIPsetDebugMsg(set, " -> objective limit in LP solver: %g (in LP: %g)\n", objlim, lp->lpiobjlim);
      }
#endif

      SCIP_CALL( SCIPsetAllocBufferArray(set, &farkascoefs, SCIPprobGetNVars(transprob)) );
      BMSclearMemoryArray(farkascoefs, SCIPprobGetNVars(transprob));

      farkaslhs = -SCIPaggrRowGetRhs(farkasrow);
      farkasactivity = -farkasactivity;

      inds = SCIPaggrRowGetInds(farkasrow);
      nnz = SCIPaggrRowGetNNz(farkasrow);

      for( v = 0; v < nnz; v++ )
      {
         int i = inds[v];

         assert(SCIPvarGetProbindex(vars[i]) == inds[v]);

         farkascoefs[i] = -SCIPaggrRowGetProbvarValue(farkasrow, i);
      }

      SCIP_CALL( SCIPrunBoundHeuristic(conflict, set, stat, origprob, transprob, tree, reopt, lp, lpi, blkmem, farkascoefs,
            &farkaslhs, &farkasactivity, curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs, iterations, marklpunsolved,
            dualproofsuccess, &valid) );

      SCIPsetFreeBufferArray(set, &farkascoefs);

      if( !valid )
         goto FLUSHPROOFSETS;

      /* analyze the conflict starting with remaining bound changes */
      SCIP_CALL( SCIPconflictAnalyzeRemainingBdchgs(conflict, blkmem, set, stat, transprob, tree, diving, \
            lbchginfoposs, ubchginfoposs, nconss, nliterals, nreconvconss, nreconvliterals) );

      /* flush conflict set storage */
      SCIP_CALL( SCIPconflictFlushConss(conflict, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, \
            eventqueue, cliquetable) );
   }

  FLUSHPROOFSETS:
   /* flush proof set */
   if( SCIPproofsetGetNVars(conflict->proofset) > 0 || conflict->nproofsets > 0 )
   {
      SCIP_CALL( SCIPconflictFlushProofset(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, \
            branchcand, eventqueue, cliquetable) );
   }

  TERMINATE:
   SCIPsetFreeBufferArray(set, &curvarubs);
   SCIPsetFreeBufferArray(set, &curvarlbs);
   SCIPsetFreeBufferArray(set, &ubchginfoposs);
   SCIPsetFreeBufferArray(set, &lbchginfoposs);
   SCIPaggrRowFree(set->scip, &farkasrow);

   return SCIP_OKAY;
}


/*
 * infeasible strong branching conflict analysis
 */

/** analyses infeasible strong branching sub problems for conflicts */
SCIP_RETCODE SCIPconflictAnalyzeStrongbranch(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_COL*             col,                /**< LP column with at least one infeasible strong branching subproblem */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict          /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   )
{
   int* cstat;
   int* rstat;
   SCIP_RETCODE retcode;
   SCIP_Bool resolve;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real newlb;
   SCIP_Real newub;
   SCIP_Bool dualraysuccess;
   int iter;
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPprobAllColsInLP(transprob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */
   assert(col != NULL);
   assert((col->sbdownvalid && SCIPsetIsGE(set, col->sbdown, lp->cutoffbound)
         && SCIPsetFeasCeil(set, col->primsol-1.0) >= col->lb - 0.5)
      || (col->sbupvalid && SCIPsetIsGE(set, col->sbup, lp->cutoffbound)
         && SCIPsetFeasFloor(set, col->primsol+1.0) <= col->ub + 0.5));
   assert(SCIPtreeGetCurrentDepth(tree) > 0);

   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_usesb )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   /* inform the LPI that strong branch is (temporarily) finished */
   SCIP_CALL( SCIPlpiEndStrongbranch(lp->lpi) );

   /* start timing */
   SCIPclockStart(conflict->sbanalyzetime, set);

   /* get temporary memory for storing current LP basis */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cstat, lp->nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rstat, lp->nlpirows) );

   /* get current LP basis */
   SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );

   /* remember old bounds */
   oldlb = col->lb;
   oldub = col->ub;

   resolve = FALSE;

   /* is down branch infeasible? */
   if( col->sbdownvalid && SCIPsetIsGE(set, col->sbdown, lp->cutoffbound) )
   {
      newub = SCIPsetFeasCeil(set, col->primsol-1.0);
      if( newub >= col->lb - 0.5 )
      {
         SCIPsetDebugMsg(set, "analyzing conflict on infeasible downwards strongbranch for variable <%s>[%g,%g] in depth %d\n",
            SCIPvarGetName(SCIPcolGetVar(col)), SCIPvarGetLbLocal(SCIPcolGetVar(col)), SCIPvarGetUbLocal(SCIPcolGetVar(col)),
            SCIPtreeGetCurrentDepth(tree));

         conflict->conflictset->conflicttype = SCIP_CONFTYPE_INFEASLP;
         conflict->nsbcalls++;

         /* change the upper bound */
         col->ub = newub;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve the LP */
         retcode = SCIPlpiSolveDual(lp->lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         if( retcode != SCIP_LPERROR )
         {
            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lp->lpi, &iter) );
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            conflict->nsbiterations += iter;
            SCIPsetDebugMsg(set, " -> resolved downwards strong branching LP in %d iterations\n", iter);

            /* perform conflict analysis on infeasible LP; last parameter guarantees status 'solved' on return */
            SCIP_CALL( conflictAnalyzeLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt, \
                  lp, branchcand, eventqueue, cliquetable, TRUE, &dualraysuccess, &iter, &nconss, &nliterals, \
                  &nreconvconss, &nreconvliterals, FALSE) );
            conflict->nsbsuccess += ((nconss > 0 || dualraysuccess) ? 1 : 0);
            conflict->nsbiterations += iter;
            conflict->nsbconfconss += nconss;
            conflict->nsbconfliterals += nliterals;
            conflict->nsbreconvconss += nreconvconss;
            conflict->nsbreconvliterals += nreconvliterals;
            if( downconflict != NULL )
               *downconflict = (nconss > 0);
         }

         /* reset the upper bound */
         col->ub = oldub;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* reset LP basis */
         SCIP_CALL( SCIPlpiSetBase(lp->lpi, cstat, rstat) );

         /* mark the LP to be resolved at the end */
         resolve = TRUE;
      }
   }

   /* is up branch infeasible? */
   if( col->sbupvalid && SCIPsetIsGE(set, col->sbup, lp->cutoffbound) )
   {
      newlb = SCIPsetFeasFloor(set, col->primsol+1.0);
      if( newlb <= col->ub + 0.5 )
      {
         SCIPsetDebugMsg(set, "analyzing conflict on infeasible upwards strongbranch for variable <%s>[%g,%g] in depth %d\n",
            SCIPvarGetName(SCIPcolGetVar(col)), SCIPvarGetLbLocal(SCIPcolGetVar(col)), SCIPvarGetUbLocal(SCIPcolGetVar(col)),
            SCIPtreeGetCurrentDepth(tree));

         conflict->conflictset->conflicttype = SCIP_CONFTYPE_INFEASLP;
         conflict->nsbcalls++;

         /* change the lower bound */
         col->lb = newlb;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve the LP */
         retcode = SCIPlpiSolveDual(lp->lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         if( retcode != SCIP_LPERROR )
         {
            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lp->lpi, &iter) );
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            conflict->nsbiterations += iter;
            SCIPsetDebugMsg(set, " -> resolved upwards strong branching LP in %d iterations\n", iter);

            /* perform conflict analysis on infeasible LP; last parameter guarantees status 'solved' on return */
            SCIP_CALL( conflictAnalyzeLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt, \
                  lp, branchcand, eventqueue, cliquetable, TRUE, &dualraysuccess, &iter, &nconss, &nliterals, \
                  &nreconvconss, &nreconvliterals, FALSE) );
            conflict->nsbsuccess += ((nconss > 0 || dualraysuccess) ? 1 : 0);
            conflict->nsbiterations += iter;
            conflict->nsbconfconss += nconss;
            conflict->nsbconfliterals += nliterals;
            conflict->nsbreconvconss += nreconvconss;
            conflict->nsbreconvliterals += nreconvliterals;
            if( upconflict != NULL )
               *upconflict = (nconss > 0);
         }

         /* reset the lower bound */
         col->lb = oldlb;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* reset LP basis */
         SCIP_CALL( SCIPlpiSetBase(lp->lpi, cstat, rstat) );

         /* mark the LP to be resolved at the end */
         resolve = TRUE;
      }
   }

   /* free temporary memory for storing current LP basis */
   SCIPsetFreeBufferArray(set, &rstat);
   SCIPsetFreeBufferArray(set, &cstat);

   assert(lp->flushed);

   /* resolve LP if something has changed in order to synchronize LPI and LP */
   if ( resolve )
   {
      /* start LP timer */
      SCIPclockStart(stat->conflictlptime, set);

      /* resolve the LP */
      SCIP_CALL( SCIPlpiSolveDual(lp->lpi) );

      /* stop LP timer */
      SCIPclockStop(stat->conflictlptime, set);
   }

   /* stop timing */
   SCIPclockStop(conflict->sbanalyzetime, set);

   /* inform the LPI that strong branch starts (again) */
   SCIP_CALL( SCIPlpiStartStrongbranch(lp->lpi) );

   return SCIP_OKAY;
}

/** analyzes an infeasible LP to find out the bound changes on variables that were responsible for the infeasibility;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible LP conflict analysis
 */
static
SCIP_RETCODE conflictAnalyzeInfeasibleLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_Bool dualraysuccess = FALSE;
   SCIP_Longint olddualproofsuccess;
   int iterations;
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(SCIPprobAllColsInLP(transprob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */

   assert(success == NULL || *success == FALSE);

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->conf_enable || set->conf_useinflp == 'o' )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "analyzing conflict on infeasible LP in depth %d (solstat: %d, objchanged: %u)\n",
      SCIPtreeGetCurrentDepth(tree), SCIPlpGetSolstat(lp), SCIPlpDivingObjChanged(lp));

   /* start timing */
   SCIPclockStart(conflict->inflpanalyzetime, set);
   conflict->ninflpcalls++;

   conflict->conflictset->conflicttype = SCIP_CONFTYPE_INFEASLP;

   olddualproofsuccess = conflict->ndualproofsinfsuccess;

   /* perform conflict analysis */
   SCIP_CALL( conflictAnalyzeLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, \
         cliquetable, SCIPlpDiving(lp), &dualraysuccess, &iterations, &nconss, &nliterals, &nreconvconss, &nreconvliterals, TRUE) );
   conflict->ninflpsuccess += ((nconss > 0 || conflict->ndualproofsinfsuccess > olddualproofsuccess) ? 1 : 0);
   conflict->ninflpiterations += iterations;
   conflict->ninflpconfconss += nconss;
   conflict->ninflpconfliterals += nliterals;
   conflict->ninflpreconvconss += nreconvconss;
   conflict->ninflpreconvliterals += nreconvliterals;
   if( success != NULL )
      *success = (nconss > 0 || conflict->ndualproofsinfsuccess > olddualproofsuccess);

   /* stop timing */
   SCIPclockStop(conflict->inflpanalyzetime, set);

   return SCIP_OKAY;
}

/** analyzes a bound exceeding LP to find out the bound changes on variables that were responsible for exceeding the
 *  primal bound;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for bound exceeding LP conflict analysis
 */
static
SCIP_RETCODE conflictAnalyzeBoundexceedingLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_Bool dualraysuccess;
   SCIP_Longint oldnsuccess;
   int iterations;
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(!SCIPlpDivingObjChanged(lp));
   assert(SCIPprobAllColsInLP(transprob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */

   assert(success == NULL || *success == FALSE);

   /* check, if bound exceeding LP conflict analysis is enabled */
   if( !set->conf_enable || set->conf_useboundlp == 'o')
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "analyzing conflict on bound exceeding LP in depth %d (solstat: %d)\n",
      SCIPtreeGetCurrentDepth(tree), SCIPlpGetSolstat(lp));

   /* start timing */
   SCIPclockStart(conflict->boundlpanalyzetime, set);
   conflict->nboundlpcalls++;

   /* mark the conflict to depend on the cutoff bound */
   conflict->conflictset->conflicttype = SCIP_CONFTYPE_BNDEXCEEDING;
   conflict->conflictset->usescutoffbound = TRUE;

   oldnsuccess = conflict->ndualproofsbndsuccess + conflict->ndualproofsinfsuccess;

   /* perform conflict analysis */
   SCIP_CALL( conflictAnalyzeLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, \
         cliquetable, SCIPlpDiving(lp), &dualraysuccess, &iterations, &nconss, &nliterals, &nreconvconss, &nreconvliterals, TRUE) );
   conflict->nboundlpsuccess += ((nconss > 0 || conflict->ndualproofsbndsuccess + conflict->ndualproofsinfsuccess > oldnsuccess) ? 1 : 0);
   conflict->nboundlpiterations += iterations;
   conflict->nboundlpconfconss += nconss;
   conflict->nboundlpconfliterals += nliterals;
   conflict->nboundlpreconvconss += nreconvconss;
   conflict->nboundlpreconvliterals += nreconvliterals;
   if( success != NULL )
      *success = (nconss > 0 || conflict->ndualproofsbndsuccess + conflict->ndualproofsinfsuccess > oldnsuccess);

   /* stop timing */
   SCIPclockStop(conflict->boundlpanalyzetime, set);

   return SCIP_OKAY;
}

/** analyzes an infeasible or bound exceeding LP to find out the bound changes on variables that were responsible for the
 *  infeasibility or for exceeding the primal bound;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible or bound exceeding LP conflict analysis;
 *  may only be called if SCIPprobAllColsInLP()
 */
SCIP_RETCODE SCIPconflictAnalyzeLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_LPSOLVALS storedsolvals;
   SCIP_COLSOLVALS* storedcolsolvals;
   SCIP_ROWSOLVALS* storedrowsolvals;
   int c;
   int r;

   if( success != NULL )
      *success = FALSE;

   /* check if the conflict analysis is applicable */
   if( !set->conf_enable || (set->conf_useinflp == 'o' && set->conf_useboundlp == 'o') )
      return SCIP_OKAY;

   /* in rare cases, it might happen that the solution stati of the LP and the LPI are out of sync; in particular this
    * happens when a new incumbent which cuts off the current node is found during the LP solving loop; in this case the
    * LP has status objlimit, but if diving has been used, the LPI only has the basis information, but is not solved
    *
    * @todo: alternatively, solve the LPI
    */
   if( !SCIPlpiWasSolved(SCIPlpGetLPI(lp)) )
      return SCIP_OKAY;

   /* LP conflict analysis is only valid, if all variables are known */
   assert( SCIPprobAllColsInLP(transprob, set, lp) );
   assert( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT
      || (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && set->lp_disablecutoff == 1) );

   /* save status */
   storedsolvals.lpsolstat = lp->lpsolstat;
   storedsolvals.lpobjval = lp->lpobjval;
   storedsolvals.primalfeasible = lp->primalfeasible;
   storedsolvals.primalchecked = lp->primalchecked;
   storedsolvals.dualfeasible = lp->dualfeasible;
   storedsolvals.dualchecked = lp->dualchecked;
   storedsolvals.solisbasic = lp->solisbasic;
   storedsolvals.lpissolved = lp->solved;

   /* store solution values */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &storedcolsolvals, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &storedrowsolvals, lp->nrows) );
   for (c = 0; c < lp->ncols; ++c)
   {
      SCIP_COL* col;

      col = lp->cols[c];
      assert( col != NULL );

      storedcolsolvals[c].primsol = col->primsol;
      storedcolsolvals[c].redcost = col->redcost;
      storedcolsolvals[c].basisstatus = col->basisstatus; /*lint !e641 !e732*/
   }
   for (r = 0; r < lp->nrows; ++r)
   {
      SCIP_ROW* row;

      row = lp->rows[r];
      assert( row != NULL );

      if ( lp->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE )
         storedrowsolvals[r].dualsol = row->dualfarkas;
      else
      {
         assert( lp->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT ||
            (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && set->lp_disablecutoff == 1) );
         storedrowsolvals[r].dualsol = row->dualsol;
      }
      storedrowsolvals[r].activity = row->activity;
      storedrowsolvals[r].basisstatus = row->basisstatus; /*lint !e641 !e732*/
   }

   /* check, if the LP was infeasible or bound exceeding */
   if( SCIPlpiIsPrimalInfeasible(SCIPlpGetLPI(lp)) )
   {
      SCIP_CALL( conflictAnalyzeInfeasibleLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, \
            reopt, lp, branchcand, eventqueue, cliquetable, success) );
   }
   else
   {
      SCIP_CALL( conflictAnalyzeBoundexceedingLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, \
            reopt, lp, branchcand, eventqueue, cliquetable, success) );
   }

   /* possibly restore solution values */
   if( lp->flushed && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_NOTSOLVED )
   {
      /* restore status */
      lp->lpsolstat = storedsolvals.lpsolstat;
      lp->lpobjval = storedsolvals.lpobjval;
      lp->primalfeasible = storedsolvals.primalfeasible;
      lp->primalchecked = storedsolvals.primalchecked;
      lp->dualfeasible = storedsolvals.dualfeasible;
      lp->dualchecked = storedsolvals.dualchecked;
      lp->solisbasic = storedsolvals.solisbasic;
      lp->solved = storedsolvals.lpissolved;

      for (c = 0; c < lp->ncols; ++c)
      {
         SCIP_COL* col;

         col = lp->cols[c];
         assert( col != NULL );
         col->primsol = storedcolsolvals[c].primsol;
         col->redcost = storedcolsolvals[c].redcost;
         col->basisstatus = storedcolsolvals[c].basisstatus; /*lint !e641 !e732*/
      }
      for (r = 0; r < lp->nrows; ++r)
      {
         SCIP_ROW* row;

         row = lp->rows[r];
         assert( row != NULL );

         if ( lp->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE )
            row->dualfarkas = storedrowsolvals[r].dualsol;
         else
         {
            assert( lp->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT );
            row->dualsol = storedrowsolvals[r].dualsol;
         }
         row->activity = storedrowsolvals[r].activity;
         row->basisstatus = storedrowsolvals[r].basisstatus; /*lint !e641 !e732*/
      }
   }
   SCIPsetFreeBufferArray(set, &storedrowsolvals);
   SCIPsetFreeBufferArray(set, &storedcolsolvals);

   return SCIP_OKAY;
}
