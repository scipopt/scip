/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_objpscostdiving.h,v 1.2 2005/01/18 09:26:46 bzfpfend Exp $"

/**@file   heur_objpscostdiving.h
 * @brief  LP diving heuristic that changes variable's objective value instead of bounds, using pseudo cost values as guide
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __HEUR_OBJPSCOSTDIVING_H__
#define __HEUR_OBJPSCOSTDIVING_H__


#include "scip.h"


/** creates the objpscostdiving heuristic and includes it in SCIP */
extern
RETCODE SCIPincludeHeurObjpscostdiving(
   SCIP*            scip                /**< SCIP data structure */
   );

#endif
