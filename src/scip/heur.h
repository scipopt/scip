/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur.h
 * @brief  methods and datastructures for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __HEUR_H__
#define __HEUR_H__


typedef struct Heur HEUR;               /**< primal heuristic */
typedef struct HeurData HEURDATA;       /**< locally defined primal heuristic data */



/** destructor of primal heuristic to free user data (called when SCIP is exiting)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    heur            : the primal heuristic itself
 */
#define DECL_HEURFREE(x) RETCODE x (SCIP* scip, HEUR* heur)

/** initialization method of primal heuristic (called when problem solving starts)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    heur            : the primal heuristic itself
 */
#define DECL_HEURINIT(x) RETCODE x (SCIP* scip, HEUR* heur)

/** deinitialization method of primal heuristic (called when problem solving exits)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    heur            : the primal heuristic itself
 */
#define DECL_HEUREXIT(x) RETCODE x (SCIP* scip, HEUR* heur)

/** execution method of primal heuristic
 *
 *  Searches for feasible primal solutions. The method is called in the node processing loop.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    heur            : the primal heuristic itself
 *    result          : pointer to store the result of the separation call
 *
 *  possible return values for *result:
 *    SCIP_FOUNDSOL   : at least one feasible primal solution was found
 *    SCIP_DIDNOTFIND : the heuristic searched, but didn't found a feasible solution
 *    SCIP_DIDNOTRUN  : the heuristic was skipped
 */
#define DECL_HEUREXEC(x) RETCODE x (SCIP* scip, HEUR* heur, RESULT* result)




#include "scip.h"
#include "def.h"
#include "retcode.h"



/** creates a primal heuristic */
extern
RETCODE SCIPheurCreate(
   HEUR**           heur,               /**< pointer to primal heuristic data structure */
   const char*      name,               /**< name of primal heuristic */
   const char*      desc,               /**< description of primal heuristic */
   char             dispchar,           /**< display character of primal heuristic */
   int              priority,           /**< priority of the primal heuristic */
   int              freq,               /**< frequency for calling primal heuristic */
   Bool             pseudonodes,        /**< call heuristic at nodes where only a pseudo solution exist? */
   DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   DECL_HEURINIT    ((*heurinit)),      /**< initialise primal heuristic */
   DECL_HEUREXIT    ((*heurexit)),      /**< deinitialise primal heuristic */
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

/** calls execution method of primal heuristic */
extern
RETCODE SCIPheurExec(
   HEUR*            heur,               /**< primal heuristic */
   const SET*       set,                /**< global SCIP settings */
   int              actdepth,           /**< depth of active node */
   Bool             actnodehaslp,       /**< is LP being processed in the active node? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** gets name of primal heuristic */
extern
const char* SCIPheurGetName(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets display character of primal heuristic */
extern
char SCIPheurGetDispchar(
   HEUR*            heur                /**< primal heuristic */
   );

/** gets user data of primal heuristic */
extern
HEURDATA* SCIPheurGetData(
   HEUR*            heur                /**< primal heuristic */
   );

/** sets user data of primal heuristic; user has to free old data in advance! */
extern
void SCIPheurSetData(
   HEUR*            heur,               /**< primal heuristic */
   HEURDATA*        heurdata            /**< new primal heuristic user data */
   );

/** gets frequency of primal heuristic */
extern
int SCIPheurGetFreq(
   HEUR*            heur                /**< primal heuristic */
   );

/** is primal heuristic initialized? */
extern
Bool SCIPheurIsInitialized(
   HEUR*            heur                /**< primal heuristic */
   );


#endif
