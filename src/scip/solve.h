/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   solve.h
 * @brief  main solving loop and node processing
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SOLVE_H__
#define __SOLVE_H__


#include "def.h"
#include "retcode.h"
#include "set.h"
#include "mem.h"
#include "stat.h"
#include "prob.h"
#include "tree.h"
#include "lp.h"
#include "price.h"
#include "sepa.h"
#include "primal.h"


extern
RETCODE SCIPsolveCIP(                   /**< main solving loop */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   PRICE*           price,              /**< pricing storage */
   SEPA*            sepa,               /**< separation storage */
   PRIMAL*          primal              /**< primal data */
   );


#endif
