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
#pragma ident "@(#) $Id: heur.h,v 1.33 2005/01/21 09:16:53 bzfpfend Exp $"

/**@file   heur.h
 * @brief  internal methods for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __HEUR_H__
#define __HEUR_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_result.h"
#include "type_set.h"
#include "type_primal.h"
#include "type_scip.h"
#include "type_heur.h"
#include "pub_heur.h"



/** creates a primal heuristic */
extern
RETCODE SCIPheurCreate(
   HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of primal heuristic */
   const char*      desc,               /**< description of primal heuristic */
   char             dispchar,           /**< display character of primal heuristic */
   int              priority,           /**< priority of the primal heuristic */
   int              freq,               /**< frequency for calling primal heuristic */
   int              freqofs,            /**< frequency offset for calling primal heuristic */
   int              maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   Bool             pseudonodes,        /**< call heuristic at nodes where only a pseudo solution exist? */
   Bool             duringplunging,     /**< call heuristic during plunging? */
   DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   DECL_HEURINIT    ((*heurinit)),      /**< initialize primal heuristic */
   DECL_HEUREXIT    ((*heurexit)),      /**< deinitialize primal heuristic */
   DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   HEURDATA*        heurdata            /**< primal heuristic data */
   );

/** calls destructor and frees memory of primal heuristic */
extern
RETCODE SCIPheurFree(
   HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes primal heuristic */
extern
RETCODE SCIPheurInit(
   HEUR*            heur,               /**< primal heuristic */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls exit method of primal heuristic */
extern
RETCODE SCIPheurExit(
   HEUR*            heur,               /**< primal heuristic */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes solution process data of primal heuristic */
extern
RETCODE SCIPheurInitsol(
   HEUR*            heur,               /**< primal heuristic */
   SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of primal heuristic */
extern
RETCODE SCIPheurExec(
   HEUR*            heur,               /**< primal heuristic */
   SET*             set,                /**< global SCIP settings */
   PRIMAL*          primal,             /**< primal data */
   int              depth,              /**< depth of current node */
   int              lpforkdepth,        /**< depth of the last node with solved LP */
   Bool             currentnodehaslp,   /**< is LP being processed in the current node? */
   Bool             plunging,           /**< is the next node to be processed a child or sibling? */
   int*             ndelayedheurs,      /**< pointer to count the number of delayed heuristics */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of primal heuristic */
extern
void SCIPheurSetPriority(
   HEUR*            heur,               /**< primal heuristic */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the primal heuristic */
   );


#endif
