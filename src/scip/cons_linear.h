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

/**@file   cons_linear.h
 * @brief  constraint handler for linear constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_LINEAR_H__
#define __CONS_LINEAR_H__


#include "scip.h"
#include "constraint.h"


#define CONSHDLR_NAME "linear"

DECL_CONSINIT(SCIPconsInit_Linear);
DECL_CONSEXIT(SCIPconsExit_Linear);
DECL_CONSFREE(SCIPconsFree_Linear);
DECL_CONSCHCK(SCIPconsChck_Linear);
DECL_CONSPROP(SCIPconsProp_Linear);


/** Constraint Handler for linear constraints. This can be made static, because the
 *  constraint handler doesn't need any specific data.
 */
static
CONSHDLR ConsHdlr_Linear = {
   CONSHDLR_NAME,
   SCIPconsInit_Linear,
   SCIPconsExit_Linear,
   SCIPconsFree_Linear,
   SCIPconsChck_Linear,
   SCIPconsProp_Linear,
   NULL
};


extern
RETCODE SCIPconsCreate_Linear(          /**< creates a linear constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             rhs,                /**< right hand side of row */
   Real             lhs,                /**< left hand side of row (for ranged rows) */
   Real             epsilon,            /**< maximal normed violation of row */
   ROWTYPE          rowtype,            /**< type of row */
   Bool             original,           /**< belongs constraint to the original problem formulation? */
   Bool             model               /**< is constraint necessary for feasibility? */
   );

#endif
