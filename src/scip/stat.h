/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   stat.h
 * @brief  problem statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STAT_H__
#define __STAT_H__


typedef struct Stat STAT;               /**< problem and runtime specific statistics */


#include "def.h"
#include "retcode.h"


/** problem and runtime specific statistics */
struct Stat
{
   int              nvaridx;            /**< number of used variable indices */
   int              ncolidx;            /**< number of used column indices */
   int              nrowidx;            /**< number of used row indices */
   int              marked_nvaridx;     /**< number of used variable indices before solving started */
   int              marked_ncolidx;     /**< number of used column indices before solving started */
   int              marked_nrowidx;     /**< number of used row indices before solving started */
   int              lpcount;            /**< internal LP counter, where all SCIPlpSolve() calls are counted */
   int              nlps;               /**< number of LPs solved (primal + dual) with at least 1 iteration */
   int              nprimallps;         /**< number of primal LPs solved */
   int              nduallps;           /**< number of dual LPs solved */
   int              nlpiterations;      /**< number of simplex iterations (primal + dual) */
   int              nprimallpiterations;/**< number of iterations in primal simplex */
   int              nduallpiterations;  /**< number of iterations in dual simplex */
   int              nstrongbranch;      /**< number of strong branching calls */
   int              nseparounds;        /**< number of separation rounds performed in actual node */
   Longint          nnodes;             /**< number of nodes processed (including active node) */
   Longint          nboundchanges;      /**< number of times a variable's bound has been changed */
   Longint          lastdispnode;       /**< last node for which an information line was displayed */
   int              ndisplines;         /**< number of displayed information lines */
   int              maxdepth;           /**< maximal depth of all processed nodes */
   int              plungedepth;        /**< actual plunging depth (successive times, a child was selected as next node) */
};


/** creates problem statistics data */
extern
RETCODE SCIPstatCreate(
   STAT**           stat                /**< pointer to problem statistics data */
   );

/** frees problem statistics data */
extern
RETCODE SCIPstatFree(
   STAT**           stat                /**< pointer to problem statistics data */
   );

/** marks statistics to be able to reset them when solving process is freed */
extern
void SCIPstatMark(
   STAT*            stat                /**< problem statistics data */
   );

/** reset statistics to the data before solving started */
extern
void SCIPstatReset(
   STAT*            stat                /**< problem statistics data */
   );

#endif
