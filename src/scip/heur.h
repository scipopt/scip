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
 *    heur            : the primal heuristic itself
 *    scip            : SCIP main data structure
 */
#define DECL_HEURFREE(x) RETCODE x (HEUR* heur, SCIP* scip)

/** initialization method of primal heuristic (called when problem solving starts)
 *
 *  input:
 *    heur            : the primal heuristic itself
 *    scip            : SCIP main data structure
 */
#define DECL_HEURINIT(x) RETCODE x (HEUR* heur, SCIP* scip)

/** deinitialization method of primal heuristic (called when problem solving exits)
 *
 *  input:
 *    heur            : the primal heuristic itself
 *    scip            : SCIP main data structure
 */
#define DECL_HEUREXIT(x) RETCODE x (HEUR* heur, SCIP* scip)

/** execution method of primal heuristic
 *
 *  Searches for feasible primal solutions. The method is called in the node processing loop.
 *
 *  input:
 *    heur            : the primal heuristic itself
 *    scip            : SCIP main data structure
 *    result          : pointer to store the result of the separation call
 *
 *  possible return values for *result:
 *    SCIP_FOUNDSOL   : at least one feasible primal solution was found
 *    SCIP_DIDNOTFIND : the heuristic searched, but didn't found a feasible solution
 *    SCIP_DIDNOTRUN  : the heuristic was skipped
 */
#define DECL_HEUREXEC(x) RETCODE x (HEUR* heur, SCIP* scip, RESULT* result)




#include "scip.h"
#include "def.h"
#include "retcode.h"



extern
RETCODE SCIPheurCreate(                 /**< creates a primal heuristic */
   HEUR**           heur,               /**< pointer to primal heuristic data structure */
   const char*      name,               /**< name of primal heuristic */
   const char*      desc,               /**< description of primal heuristic */
   char             dispchar,           /**< display character of primal heuristic */
   int              priority,           /**< priority of the primal heuristic */
   int              freq,               /**< frequency for calling primal heuristic */
   DECL_HEURFREE((*heurfree)),          /**< destructor of primal heuristic */
   DECL_HEURINIT((*heurinit)),          /**< initialise primal heuristic */
   DECL_HEUREXIT((*heurexit)),          /**< deinitialise primal heuristic */
   DECL_HEUREXEC((*heurexec)),          /**< execution method of primal heuristic */
   HEURDATA*        heurdata            /**< primal heuristic data */
   );

extern
RETCODE SCIPheurFree(                   /**< calls destructor and frees memory of primal heuristic */
   HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPheurInit(                   /**< initializes primal heuristic */
   HEUR*            heur,               /**< primal heuristic */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPheurExit(                   /**< calls exit method of primal heuristic */
   HEUR*            heur,               /**< primal heuristic */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPheurExec(                   /**< calls execution method of primal heuristic */
   HEUR*            heur,               /**< primal heuristic */
   const SET*       set,                /**< global SCIP settings */
   int              actdepth,           /**< depth of active node */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
const char* SCIPheurGetName(            /**< gets name of primal heuristic */
   HEUR*            heur                /**< primal heuristic */
   );

extern
char SCIPheurGetDispchar(               /**< gets display character of primal heuristic */
   HEUR*            heur                /**< primal heuristic */
   );

extern
HEURDATA* SCIPheurGetData(              /**< gets user data of primal heuristic */
   HEUR*            heur                /**< primal heuristic */
   );

extern
void SCIPheurSetData(                   /**< sets user data of primal heuristic; user has to free old data in advance! */
   HEUR*            heur,               /**< primal heuristic */
   HEURDATA*        heurdata            /**< new primal heuristic user data */
   );

extern
int SCIPheurGetFreq(                    /**< gets frequency of primal heuristic */
   HEUR*            heur                /**< primal heuristic */
   );

extern
Bool SCIPheurIsInitialized(             /**< is primal heuristic initialized? */
   HEUR*            heur                /**< primal heuristic */
   );


#endif
