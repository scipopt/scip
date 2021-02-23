/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dpterms_core.c
 * @brief  Core of dynamic programming solver for Steiner tree (sub-) problems with small number of terminals
 * @author Daniel Rehfeldt
 *
 * Contains core methods.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "scip/rbtree.h"
#include "dpterms.h"
#include "dptermsinterns.h"
#include "stpbitset.h"
#include "stpvector.h"
#include "stpprioqueue.h"

/*
 * Local methods
 */


/*
 * Interface methods
 */


/** solves problem */
SCIP_RETCODE dpterms_coreSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< (compressed) graph */
   DPSOLVER*             dpsolver            /**< solver */
)
{
   // solve

   return SCIP_OKAY;
}
