/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_countsols.h,v 1.3 2008/01/30 14:19:23 bzfheinz Exp $"

/**@file   cons_countsols.h
 * @brief  constraint handler for counting feasible solutions
 * @author Stefan Heinz
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_COUNTSOLS_H__
#define __SCIP_CONS_COUNTSOLS_H__

#include "scip/scip.h"

/** dialog execution method for the count command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecCount);

/** creates the handler for countsol constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrCountsols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* execute counting */
extern
SCIP_RETCODE SCIPcount(
   SCIP*                 scip                /**< SCIP data structure */
   );

#if 0
/* returns TRUE if the counting process was correct; otherwise FALSE */
extern
SCIP_Bool SCIPisCountValid(
   SCIP*                 scip                /**< SCIP data structure */
   ); 
#endif

/** returns number of feasible solutions found as SCIP_Longint; if the number does not fit into 
 *  a SCIP_Longint the valid flag is set to FALSE */
extern
SCIP_Longint SCIPgetNCountedSols(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_Bool*            valid              /**< pointer to store if the return value is valid */             
   ); 

/** returns number of counted solutions as string */
extern
void SCIPgetNCountedSolsstr(
   SCIP*                 scip,              /**< SCIP data structure */
   char**                buffer,             /**< buffer to store the number for counted solutions */
   int                   buffersize,         /**< buffer size */
   int*                  requiredsize        /**< pointer to store the required size */
   );

/** returns number of counted feasible unimodular subtrees */
extern
SCIP_Longint SCIPgetNCountedFeasUS(
   SCIP*                 scip                /**< SCIP data structure */
   ); 

#endif
