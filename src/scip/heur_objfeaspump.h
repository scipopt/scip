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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_objfeaspump.h,v 1.1 2005/01/31 15:21:03 bzfberth Exp $"

/**@file   heur_objfeaspump.h
 * @brief  feasibility pump heuristic by Fischetti, Glover and Lodi 
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __HEUR_OBJFEASPUMP_H__
#define __HEUR_OBJFEASPUMP_H__


#include "scip.h"


/** creates the objfeaspump heuristic and includes it in SCIP */
extern
RETCODE SCIPincludeHeurObjfeaspump(
   SCIP*            scip                /**< SCIP data structure */
   );

#endif
