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

/**@file   sepa.h
 * @brief  methods and datastructures for separating cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SEPA_H__
#define __SEPA_H__


typedef struct Sepa SEPA;               /**< storage for separated variables */


#include "def.h"
#include "retcode.h"
#include "set.h"
#include "mem.h"
#include "lp.h"
#include "tree.h"


/** creates separation storage */
extern
RETCODE SCIPsepaCreate(
   SEPA**           sepa                /**< pointer to store separation storage */
   );

/** frees separation storage */
extern
RETCODE SCIPsepaFree(
   SEPA**           sepa                /**< pointer to store separation storage */
   );

/** adds cut to separation storage and captures it */
extern
RETCODE SCIPsepaAddCut(
   SEPA*            sepa,               /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   ROW*             cut,                /**< separated cut */
   Real             score,              /**< separation score of cut (the larger, the better the cut) */
   Bool             root                /**< are we at the root node? */
   );

/** adds cuts to the LP and clears separation storage */
extern
RETCODE SCIPsepaApplyCuts(
   SEPA*            sepa,               /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< LP data */
   );

/** get number of cuts in the separation storage */
extern
int SCIPsepaGetNCuts(
   SEPA*            sepa                /**< separation storage */
   );



#endif
