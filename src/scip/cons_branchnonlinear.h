/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_branchnonlinear.h,v 1.2 2009/07/28 10:05:26 bzfviger Exp $"

/**@file   cons_branchnonlinear.h
 * @brief  constraint handler for branching on variables in nonlinear (nonconvex) constraints
 * @author Stefan Vigerske
 * 
 * We cannot put this into a branching rule plugin, since these rules are only designed for branchings on integer variables. 
 * This is ok, since the SCIP design bases on a definition of CIP where the CIP becomes an LP after fixing all integer variables.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_BRANCHNONLINEAR_H__
#define __SCIP_CONS_BRANCHNONLINEAR_H__


#include "scip/scip.h"


/** creates the branching rule for nonlinear variables and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrBranchNonlinear(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** Updates or initializes the infeasibility of a variable.
 * If called the first time for some variable, then this variable is added to the list of branching candidates.
 */
extern
SCIP_RETCODE SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(
   SCIP*               scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*      conshdlr,           /**< constraint handler */
   SCIP_VAR*           var,                /**< variable */
   SCIP_Real           varinfeasibility    /**< infeasibility of variable */
   );

#endif
