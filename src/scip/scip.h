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

/**@file   scip.h
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_H__
#define __SCIP_H__


typedef struct Scip SCIP;               /**< SCIP main data structure */


#include <stdio.h>

#include "def.h"
#include "retcode.h"
#include "constraint.h"
#include "lp.h"


#define CHECK_SCIP(x) { RETCODE _retcode_; \
                        if( (_retcode_ = (x)) != SCIP_OKAY ) \
                          SCIPerror(stderr, _retcode_, __FILE__, __LINE__); \
                      }

extern
void SCIPerror(                         /**< prints error message and aborts program execution */
   FILE*            errout,             /**< file stream to write error message */
   RETCODE          retcode,            /**< SCIP return code causing the error */
   const char*      filename,           /**< source code file name */
   int              line                /**< source line */
   );

extern
RETCODE SCIPcreate(                     /**< creates and initializes SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   );

extern
RETCODE SCIPfree(                       /**< frees SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   );

extern
RETCODE SCIPcreateProb(                 /**< creates empty problem and initializes all solving data structures */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< problem name */
   );

extern
RETCODE SCIPfreeProb(                   /**< frees problem and branch-and-bound data structures */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
RETCODE SCIPsolve(                      /**< solves problem */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
RETCODE SCIPfreeSolve(                  /**< frees all solution process data, only original problem is kept */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
RETCODE SCIPcreateCol(                  /**< create variable with empty column */
   SCIP*            scip,               /**< SCIP data structure */
   COL**            col,                /**< pointer to column object */
   const char*      name,               /**< name of column */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   COLTYPE          coltype             /**< type of variable */
   );

extern
RETCODE SCIPcreateRow(                  /**< creates an LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row,                /**< pointer to row */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             rhs,                /**< right hand side of row */
   Real             lhs,                /**< left hand side of row (for ranged rows) */
   Real             epsilon,            /**< maximal normed violation of row */
   ROWTYPE          rowtype             /**< type of row */
   );

extern
RETCODE SCIPcaptureRow(                 /**< increases usage counter of LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< row to capture */
   );

extern
RETCODE SCIPcreateConstraint(           /**< creates a constraint of the given constraint handler */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             original,           /**< belongs constraint to the original problem formulation? */
   Bool             model               /**< is constraint necessary for feasibility? */
   );

extern
RETCODE SCIPaddConstraint(              /**< adds constraint to the problem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   );

extern
RETCODE SCIPaddLocalConstraint(         /**< adds local constraint to the actual subproblem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   );



/*
 * debug methods
 */

#ifndef NDEBUG

extern
void SCIPdebugMemory(                   /**< prints output about used memory */
   SCIP*            scip                /**< SCIP data structure */
   );

#endif


#endif
