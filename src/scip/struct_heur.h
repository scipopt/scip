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
#pragma ident "@(#) $Id: struct_heur.h,v 1.10 2005/02/07 14:08:28 bzfpfend Exp $"

/**@file   struct_heur.h
 * @brief  datastructures for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_HEUR_H__
#define __STRUCT_HEUR_H__


#include "def.h"
#include "type_clock.h"
#include "type_heur.h"


/** primal heuristics data */
struct Heur
{
   Longint          ncalls;             /**< number of times, this heuristic was called */
   Longint          nsolsfound;         /**< number of feasible primal solutions found so far by this heuristic */
   char*            name;               /**< name of primal heuristic */
   char*            desc;               /**< description of primal heuristic */
   DECL_HEURFREE    ((*heurfree));      /**< destructor of primal heuristic */
   DECL_HEURINIT    ((*heurinit));      /**< initialize primal heuristic */
   DECL_HEUREXIT    ((*heurexit));      /**< deinitialize primal heuristic */
   DECL_HEURINITSOL ((*heurinitsol));   /**< solving process initialization method of primal heuristic */
   DECL_HEUREXITSOL ((*heurexitsol));   /**< solving process deinitialization method of primal heuristic */
   DECL_HEUREXEC    ((*heurexec));      /**< execution method of primal heuristic */
   HEURDATA*        heurdata;           /**< primal heuristics local data */
   CLOCK*           clock;              /**< heuristic execution time */
   int              priority;           /**< priority of the primal heuristic */
   int              freq;               /**< frequency for calling primal heuristic */
   int              freqofs;            /**< frequency offset for calling primal heuristic */
   int              maxdepth;           /**< maximal depth level to call heuristic at (-1: no limit) */
   int              delaypos;           /**< position in the delayed heuristics queue, or -1 if not delayed */
   Bool             pseudonodes;        /**< call heuristic at nodes where only a pseudo solution exist? */
   Bool             duringplunging;     /**< call heuristic during plunging? */
   Bool             initialized;        /**< is primal heuristic initialized? */
   char             dispchar;           /**< display character of primal heuristic */
};


#endif
