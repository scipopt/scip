/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   solex.h
 * @brief  internal methods for storing exact primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SOLEX_H__
#define __SCIP_SOLEX_H__


#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_var.h"
#include "scip/type_solex.h"
#include "scip/type_heur.h"
#include "scip/pub_solex.h"
#ifdef WITH_EXACTSOLVE
#include "gmp.h"


#ifdef __cplusplus
extern "C" {
#endif

/** creates exact primal CIP solution, initialized to zero */
extern
SCIP_RETCODE SCIPsolexCreate(
   SCIP_SOLEX**          sol,                /**< pointer to exact primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** frees exact primal CIP solution */
extern
SCIP_RETCODE SCIPsolexFree(
   SCIP_SOLEX**          sol,                /**< pointer to exact primal CIP solution */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** sets value of variable in exact primal CIP solution */
extern
SCIP_RETCODE SCIPsolexSetVal(
   SCIP_SOLEX*           sol,                /**< exact primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to add to solution */
   const mpq_t           obj,                /**< objective value of variable */
   const mpq_t           lb,                 /**< global lower bound of variable */
   const mpq_t           val                 /**< solution value of variable */
   );

/** returns value of variable in exact primal CIP solution */
extern
void SCIPsolexGetVal(
   SCIP_SOLEX*           sol,                /**< exact primal CIP solution */
   SCIP_VAR*             var,                /**< variable to get value for */
   mpq_t                 val                 /**< pointer to store value of variable */
   );

/** gets objective value of exact primal CIP solution in transformed problem */
extern
void SCIPsolexGetObj(
   SCIP_SOLEX*           sol,                /**< primal CIP solution */
   mpq_t                 obj                 /**< pointer to store objective value of solution */
   );

/** returns whether the given exact solutions in transformed space are equal */
extern
SCIP_Bool SCIPsolexsAreEqual(
   SCIP_SOLEX*           sol1,               /**< first exact primal CIP solution */
   SCIP_SOLEX*           sol2,               /**< second exact primal CIP solution */
   SCIP_PROB*            prob                /**< transformed problem data */
   );

/** outputs non-zero elements of exact solution to file stream */
extern
SCIP_RETCODE SCIPsolexPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOLEX*           sol,                /**< primal CIP solution */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

#ifdef __cplusplus
}
#endif

#endif

#endif
