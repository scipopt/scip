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
#pragma ident "@(#) $Id: struct_prop.h,v 1.5 2005/02/08 14:22:30 bzfpfend Exp $"

/**@file   struct_prop.h
 * @brief  datastructures for propagators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_PROP_H__
#define __STRUCT_PROP_H__


#include "def.h"
#include "type_clock.h"
#include "type_prop.h"


/** propagators data */
struct Prop
{
   Longint          ncalls;             /**< number of times, this propagator was called */
   Longint          ncutoffs;           /**< number of cutoffs found so far by this constraint handler */
   Longint          ndomredsfound;      /**< number of domain reductions found so far by this constraint handler */
   char*            name;               /**< name of propagator */
   char*            desc;               /**< description of propagator */
   DECL_PROPFREE    ((*propfree));      /**< destructor of propagator */
   DECL_PROPINIT    ((*propinit));      /**< initialize propagator */
   DECL_PROPEXIT    ((*propexit));      /**< deinitialize propagator */
   DECL_PROPINITSOL ((*propinitsol));   /**< solving process initialization method of propagator */
   DECL_PROPEXITSOL ((*propexitsol));   /**< solving process deinitialization method of propagator */
   DECL_PROPEXEC    ((*propexec));      /**< execution method of propagator */
   DECL_PROPRESPROP ((*propresprop));   /**< propagation conflict resolving method */
   PROPDATA*        propdata;           /**< propagators local data */
   CLOCK*           clock;              /**< propagation time */
   int              priority;           /**< priority of the propagator */
   int              freq;               /**< frequency for calling propagator */
   Bool             delay;              /**< should propagator be delayed, if other propagators found reductions? */
   Bool             wasdelayed;         /**< was the propagator delayed at the last call? */
   Bool             initialized;        /**< is propagator initialized? */
};


#endif
