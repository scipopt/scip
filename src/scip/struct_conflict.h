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
#pragma ident "@(#) $Id: struct_conflict.h,v 1.2 2004/01/16 11:25:04 bzfpfend Exp $"

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
   int              priority;           /**< priority of the conflict handler */
   DECL_CONFLICTFREE((*conflictfree));  /**< destructor of conflict handler */
   DECL_CONFLICTINIT((*conflictinit));  /**< initialize conflict handler */
   DECL_CONFLICTEXIT((*conflictexit));  /**< deinitialize conflict handler */
   DECL_CONFLICTEXEC((*conflictexec));  /**< conflict processing method of conflict handler */
   CONFLICTHDLRDATA* conflicthdlrdata;  /**< conflict handler data */
   Bool             initialized;        /**< is conflict handler initialized? */
};

/** conflict analysis data structure for propagation conflicts */
struct Conflict
{
   CLOCK*           analyzetime;        /**< time used for propagation conflict analysis */
   PQUEUE*          varqueue;           /**< unprocessed conflict variables */
   VAR**            conflictvars;       /**< variables resembling the conflict clause */
   int              conflictvarssize;   /**< size of conflictvars array */
   int              nconflictvars;      /**< number of variables in the conflict set (used slots of conflictvars array) */
   Longint          ncalls;             /**< number of calls to propagation conflict analysis */
   Longint          nconflicts;         /**< number of valid conflicts detected in propagation conflict analysis */
};

/** conflict analysis data structure for infeasible LP conflicts */
struct LPConflict
{
   CLOCK*           analyzetime;        /**< time used for infeasible LP conflict analysis */
   LPI*             lpi;                /**< LP problem object for the alternative polyhedron */
   Longint          ncalls;             /**< number of calls to infeasible LP conflict analysis */
   Longint          nconflicts;         /**< number of valid conflicts detected in infeasible LP conflict analysis */
   Longint          nlpiterations;      /**< total number of iterations used for solving alternative LPs */
};


#endif
