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

/**@file   sepa.h
 * @brief  methods and datastructures for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SEPA_H__
#define __SEPA_H__


typedef struct Sepa SEPA;               /**< separator */
typedef struct SepaData SEPADATA;       /**< locally defined separator data */



/** destructor of separator to free user data (called when SCIP is exiting)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    sepa            : the separator itself
 */
#define DECL_SEPAFREE(x) RETCODE x (SCIP* scip, SEPA* sepa)

/** initialization method of separator (called when problem solving starts)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    sepa            : the separator itself
 */
#define DECL_SEPAINIT(x) RETCODE x (SCIP* scip, SEPA* sepa)

/** deinitialization method of separator (called when problem solving exits)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    sepa            : the separator itself
 */
#define DECL_SEPAEXIT(x) RETCODE x (SCIP* scip, SEPA* sepa)

/** execution method of separator
 *
 *  Searches for cutting planes. The method is called in the LP solving loop.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    sepa            : the separator itself
 *    result          : pointer to store the result of the separation call
 *
 *  possible return values for *result:
 *    SCIP_CUTOFF     : at least one unmodifiable row is infeasible in the variable's bounds -> node is infeasible
 *    SCIP_SEPARATED  : a cutting plane was generated
 *    SCIP_REDUCEDDOM : no cutting plane was generated, but at least one domain was reduced
 *    SCIP_CONSADDED  : no cutting plane or domain reductions, but at least one additional constraint was generated
 *    SCIP_DIDNOTFIND : the separator searched, but didn't found a feasible cutting plane
 *    SCIP_DIDNOTRUN  : the separator was skipped
 */
#define DECL_SEPAEXEC(x) RETCODE x (SCIP* scip, SEPA* sepa, RESULT* result)




#include "scip.h"
#include "def.h"
#include "retcode.h"
#include "sepastore.h"



/** creates a separator */
extern
RETCODE SCIPsepaCreate(
   SEPA**           sepa,               /**< pointer to separator data structure */
   const char*      name,               /**< name of separator */
   const char*      desc,               /**< description of separator */
   int              priority,           /**< priority of separator */
   int              freq,               /**< frequency for calling separator */
   DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   DECL_SEPAINIT    ((*sepainit)),      /**< initialise separator */
   DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialise separator */
   DECL_SEPAEXEC    ((*sepaexec)),      /**< execution method of separator */
   SEPADATA*        sepadata            /**< separator data */
   );

/** calls destructor and frees memory of separator */
extern
RETCODE SCIPsepaFree(
   SEPA**           sepa,               /**< pointer to separator data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes separator */
extern
RETCODE SCIPsepaInit(
   SEPA*            sepa,               /**< separator */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls exit method of separator */
extern
RETCODE SCIPsepaExit(
   SEPA*            sepa,               /**< separator */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls execution method of separator */
extern
RETCODE SCIPsepaExec(
   SEPA*            sepa,               /**< separator */
   const SET*       set,                /**< global SCIP settings */
   SEPASTORE*       sepastore,          /**< separation storage */
   int              actdepth,           /**< depth of active node */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** gets name of separator */
extern
const char* SCIPsepaGetName(
   SEPA*            sepa                /**< separator */
   );

/** gets user data of separator */
extern
SEPADATA* SCIPsepaGetData(
   SEPA*            sepa                /**< separator */
   );

/** sets user data of separator; user has to free old data in advance! */
extern
void SCIPsepaSetData(
   SEPA*            sepa,               /**< separator */
   SEPADATA*        sepadata            /**< new separator user data */
   );

/** gets frequency of separator */
extern
int SCIPsepaGetFreq(
   SEPA*            sepa                /**< separator */
   );

/** gets time in seconds used in this separator */
extern
Real SCIPsepaGetTime(
   SEPA*            sepa                /**< separator */
   );

/** gets the number of times, the separator was called and tried to find a solution */
extern
Longint SCIPsepaGetNCalls(
   SEPA*            sepa                /**< separator */
   );

/** gets the number of cutting planes found by this separator */
extern
Longint SCIPsepaGetNCutsFound(
   SEPA*            sepa                /**< separator */
   );

/** is separator initialized? */
extern
Bool SCIPsepaIsInitialized(
   SEPA*            sepa                /**< separator */
   );


#endif
