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
#pragma ident "@(#) $Id: prop_rootredcost.h,v 1.5 2009/07/31 11:37:16 bzfwinkm Exp $"

/**@file   prop_rootredcost.h
 * @brief  reduced cost strengthening at the root node
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_ROOTREDCOST_H__
#define __SCIP_PROP_ROOTREDCOST_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the root node reduced cost strengthening propagator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePropRootredcost(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
