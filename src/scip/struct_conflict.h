/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_conflict.h,v 1.8 2004/08/25 15:40:07 bzfpfend Exp $"

/**@file   struct_conflict.h
 * @brief  datastructures for conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_CONFLICT_H__
#define __STRUCT_CONFLICT_H__


#include "def.h"
#include "type_clock.h"
#include "type_lpi.h"
#include "type_misc.h"
#include "type_var.h"
#include "type_conflict.h"



/** conflict handler */
struct Conflicthdlr
{
   char*            name;               /**< name of conflict handler */
   char*            desc;               /**< description of conflict handler */
   DECL_CONFLICTFREE((*conflictfree));  /**< destructor of conflict handler */
   DECL_CONFLICTINIT((*conflictinit));  /**< initialize conflict handler */
   DECL_CONFLICTEXIT((*conflictexit));  /**< deinitialize conflict handler */
   DECL_CONFLICTEXEC((*conflictexec));  /**< conflict processing method of conflict handler */
   CONFLICTHDLRDATA* conflicthdlrdata;  /**< conflict handler data */
   int              priority;           /**< priority of the conflict handler */
   Bool             initialized;        /**< is conflict handler initialized? */
};

/** conflict analysis data structure */
struct Conflict
{
   Longint          npropcalls;         /**< number of calls to propagation conflict analysis */
   Longint          npropconflicts;     /**< number of valid conflicts detected in propagation conflict analysis */
   Longint          nlpcalls;           /**< number of calls to infeasible LP conflict analysis */
   Longint          nlpconflicts;       /**< number of valid conflicts detected in infeasible LP conflict analysis */
   Longint          nlpiterations;      /**< total number of LP iterations used in LP conflict analysis */
   Longint          nsbcalls;           /**< number of calls to infeasible strong branching conflict analysis */
   Longint          nsbconflicts;       /**< number of valid conflicts detected in infeasible strong branching analysis */
   Longint          nsbiterations;      /**< total number of LP iterations used in infeasible strong branching analysis */
   Longint          npseudocalls;       /**< number of calls to pseudo solution conflict analysis */
   Longint          npseudoconflicts;   /**< number of valid conflicts detected in pseudo solution conflict analysis */
   CLOCK*           propanalyzetime;    /**< time used for propagation conflict analysis */
   CLOCK*           lpanalyzetime;      /**< time used for infeasible LP conflict analysis */
   CLOCK*           sbanalyzetime;      /**< time used for infeasible LP conflict analysis */
   CLOCK*           pseudoanalyzetime;  /**< time used for pseudo solution conflict analysis */
   LPI*             lpi;                /**< LP problem object for the alternative polyhedron */
   PQUEUE*          binbdchgqueue;      /**< unprocessed conflict bound changes on binary variables */
   PQUEUE*          nonbinbdchgqueue;   /**< unprocessed conflict bound changes on non-binary variables */
   VAR**            conflictvars;       /**< variables resembling the conflict clause */
   int              conflictvarssize;   /**< size of conflictvars array */
   int              nconflictvars;      /**< number of variables in the conflict set (used slots of conflictvars array) */
   int              count;              /**< conflict set counter to label conflict variables with */
};


#endif
