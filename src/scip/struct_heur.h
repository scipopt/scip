/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_heur.h,v 1.2 2003/12/04 15:11:31 bzfpfend Exp $"

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
   char*            name;               /**< name of primal heuristic */
   char*            desc;               /**< description of primal heuristic */
   char             dispchar;           /**< display character of primal heuristic */
   int              priority;           /**< priority of the primal heuristic */
   int              freq;               /**< frequency for calling primal heuristic */
   DECL_HEURFREE    ((*heurfree));      /**< destructor of primal heuristic */
   DECL_HEURINIT    ((*heurinit));      /**< initialize primal heuristic */
   DECL_HEUREXIT    ((*heurexit));      /**< deinitialize primal heuristic */
   DECL_HEUREXEC    ((*heurexec));      /**< execution method of primal heuristic */
   HEURDATA*        heurdata;           /**< primal heuristics local data */
   CLOCK*           clock;              /**< heuristic execution time */
   Longint          ncalls;             /**< number of times, this heuristic was called */
   Longint          nsolsfound;         /**< number of feasible primal solutions found so far by this heuristic */
   Bool             pseudonodes;        /**< call heuristic at nodes where only a pseudo solution exist? */
   Bool             initialized;        /**< is primal heuristic initialized? */
};


#endif
