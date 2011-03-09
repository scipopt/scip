/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_prop.h
 * @brief  datastructures for propagators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_PROP_H__
#define __SCIP_STRUCT_PROP_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_prop.h"

#ifdef __cplusplus
extern "C" {
#endif

/** propagators data */
struct SCIP_Prop
{
   SCIP_Longint          ncalls;             /**< number of times, this propagator was called */
   SCIP_Longint          nrespropcalls;      /**< number of times, the resolve propagtion was called */
   SCIP_Longint          ncutoffs;           /**< number of cutoffs found so far by this constraint handler */
   SCIP_Longint          ndomredsfound;      /**< number of domain reductions found so far by this constraint handler */
   char*                 name;               /**< name of propagator */
   char*                 desc;               /**< description of propagator */
   SCIP_DECL_PROPCOPY    ((*propcopy));      /**< copy method of propagator or NULL if you don't want to copy your plugin into subscips */
   SCIP_DECL_PROPFREE    ((*propfree));      /**< destructor of propagator */
   SCIP_DECL_PROPINIT    ((*propinit));      /**< initialize propagator */
   SCIP_DECL_PROPEXIT    ((*propexit));      /**< deinitialize propagator */
   SCIP_DECL_PROPINITSOL ((*propinitsol));   /**< solving process initialization method of propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol));   /**< solving process deinitialization method of propagator */
   SCIP_DECL_PROPEXEC    ((*propexec));      /**< execution method of propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop));   /**< propagation conflict resolving method */
   SCIP_PROPDATA*        propdata;           /**< propagators local data */
   SCIP_CLOCK*           proptime;           /**< propagation time */
   SCIP_CLOCK*           resproptime;        /**< time used for resolve propagation of this constraint handler */
   int                   priority;           /**< priority of the propagator */
   int                   freq;               /**< frequency for calling propagator */
   SCIP_Bool             delay;              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_Bool             wasdelayed;         /**< was the propagator delayed at the last call? */
   SCIP_Bool             initialized;        /**< is propagator initialized? */
};

#ifdef __cplusplus
}
#endif

#endif
