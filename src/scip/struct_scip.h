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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_scip.h,v 1.17 2005/08/09 16:27:07 bzfpfend Exp $"

/**@file   struct_scip.h
 * @brief  SCIP main data structure
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SCIP_H__
#define __SCIP_STRUCT_SCIP_H__


#include "scip/def.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_clock.h"
#include "scip/type_interrupt.h"
#include "scip/type_mem.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_implics.h"
#include "scip/type_prob.h"
#include "scip/type_primal.h"
#include "scip/type_tree.h"
#include "scip/type_pricestore.h"
#include "scip/type_sepastore.h"
#include "scip/type_cutpool.h"
#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_dialog.h"


/** SCIP main data structure */
struct Scip
{
   /* INIT */
   MEM*             mem;                /**< block memory buffers */
   SET*             set;                /**< global SCIP settings */
   INTERRUPT*       interrupt;          /**< CTRL-C interrupt data */
   DIALOGHDLR*      dialoghdlr;         /**< dialog handler for user interface */
   CLOCK*           totaltime;          /**< total SCIP running time */

   /* PROBLEM */
   STAT*            stat;               /**< dynamic problem statistics */
   PROB*            origprob;           /**< original problem data */

   /* TRANSFORMED */
   EVENTFILTER*     eventfilter;        /**< event filter for global (not variable dependent) events */
   EVENTQUEUE*      eventqueue;         /**< event queue to cache events and process them later (bound change events) */
   BRANCHCAND*      branchcand;         /**< storage for branching candidates */
   LP*              lp;                 /**< LP data */
   PRIMAL*          primal;             /**< primal data and solution storage */
   TREE*            tree;               /**< branch and bound tree */
   CONFLICT*        conflict;           /**< conflict analysis data */
   CLIQUETABLE*     cliquetable;        /**< collection of cliques */
   PROB*            transprob;          /**< transformed problem after presolve */

   /* SOLVING */
   PRICESTORE*      pricestore;         /**< storage for priced variables */
   SEPASTORE*       sepastore;          /**< storage for separated cuts */
   CUTPOOL*         cutpool;            /**< global cut pool */
};


#endif
