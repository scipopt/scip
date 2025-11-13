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

/**@file   struct_sym.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for symmetry method handlers
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SYM_H__
#define __SCIP_STRUCT_SYM_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_sym.h"

#ifdef __cplusplus
extern "C" {
#endif

/** symmetry group */
struct SCIP_Sym
{
   SCIP_SYMHDLR*         symhdlr;            /**< symmetry handler for this symmetry group */
};

/** symmetry handler */
struct SCIP_Symhdlr
{
   SCIP_Longint          nsepacalls;         /**< number of times, the separator was called */
   SCIP_Longint          npropcalls;         /**< number of times, the propagator was called */
   SCIP_Longint          ncutoffs;           /**< number of cutoffs found so far by this symmetry handler */
   SCIP_Longint          ncutsfound;         /**< number of cuts found by this symmetry handler */
   SCIP_Longint          ncutsapplied;       /**< number of cuts found by this symmetry handler applied to lp */
   SCIP_Longint          nconssfound;        /**< number of additional constraints added by this symmetry handler */
   SCIP_Longint          ndomredsfound;      /**< number of domain reductions found so far by this symmetry handler */
   char*                 name;               /**< name of symmetry handler */
   char*                 desc;               /**< description of symmetry handler */
   SCIP_DECL_SYMTRYADD   ((*symtryadd));     /**< try to add symmetry handler */
   SCIP_DECL_SYMCOPY     ((*symcopy));       /**< copy method of symmetry handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_SYMFREE     ((*symfree));       /**< destructor of symmetry handler */
   SCIP_DECL_SYMINIT     ((*syminit));       /**< initialize symmetry handler */
   SCIP_DECL_SYMEXIT     ((*symexit));       /**< deinitialize symmetry handler */
   SCIP_DECL_SYMTRANS    ((*symtrans));      /**< transform symmetry handler data into data belonging to the transformed problem */
   SCIP_DECL_SYMSEPALP   ((*symsepalp));     /**< separate cutting planes for LP solution */
   SCIP_DECL_SYMSEPASOL  ((*symsepasol));    /**< separate cutting planes for arbitrary primal solution */
   SCIP_DECL_SYMPROP     ((*symprop));       /**< propagate variable domains */
   SCIP_DECL_SYMPRESOL   ((*sympresol));     /**< presolving method */
   SCIP_SYMHDLRDATA*     symhdlrdata;        /**< symmetry handler data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this symmetry handler for the next stages */
   SCIP_CLOCK*           presoltime;         /**< time used for presolving of this symmetry handler */
   SCIP_CLOCK*           sepatime;           /**< time used for separation of this symmetry handler */
   SCIP_CLOCK*           proptime;           /**< time used for propagation of this symmetry handler */
   int                   tryaddpriority;     /**< priority of the symmetry handler for trying to add its symmetry handling methods */
   int                   sepapriority;       /**< priority of the symmetry handler for separation */
   int                   sepafreq;           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   propfreq;           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int                   maxprerounds;       /**< maximal number of presolving rounds the symmetry handler participates in (-1: no limit) */
   int                   nfixedvars;         /**< total number of variables fixed by this presolver */
   int                   naggrvars;          /**< total number of variables aggregated by this presolver */
   int                   nchgvartypes;       /**< total number of variable type changes by this presolver */
   int                   nchgbds;            /**< total number of variable bounds tightened by this presolver */
   int                   naddholes;          /**< total number of domain holes added by this presolver */
   int                   ndelconss;          /**< total number of deleted constraints by this presolver */
   int                   naddconss;          /**< total number of added constraints by this presolver */
   int                   nupgdconss;         /**< total number of upgraded constraints by this presolver */
   int                   nchgcoefs;          /**< total number of changed coefficients by this presolver */
   int                   nchgsides;          /**< total number of changed left or right hand sides by this presolver */
   int                   npresolcalls;       /**< number of times the symmetry handler was called in presolving and tried to find reductions */
   SCIP_Bool             delaysepa;          /**< should separation method be delayed, if other separators found cuts? */
   SCIP_Bool             delayprop;          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_Bool             sepalpwasdelayed;   /**< was the LP separation method delayed at the last call? */
   SCIP_Bool             sepasolwasdelayed;  /**< was the SOL separation method delayed at the last call? */
   SCIP_Bool             propwasdelayed;     /**< was the propagation method delayed at the last call? */
   SCIP_PROPTIMING       proptiming;         /**< positions in the node solving loop where propagation method of symmetry handlers should be executed */
   SCIP_PRESOLTIMING     presoltiming;       /**< timing mask of the symmetry handler's presolving method */
};

#ifdef __cplusplus
}
#endif

#endif
