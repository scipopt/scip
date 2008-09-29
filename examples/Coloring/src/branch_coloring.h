/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch_coloring.h,v 1.2 2008/09/29 19:49:57 bzfheinz Exp $"

/**@file   branch_coloring.h
 * @brief  coloring branching rule
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_COLORING_H__
#define __SCIP_BRANCH_COLORING_H__


#include "scip/scip.h"
#include "probdata_coloring.h"
#include "cons_storeGraph.h"
#include "scip/cons_linear.h"


/** creates the coloring branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchruleColoring(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
