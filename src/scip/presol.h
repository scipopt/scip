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

/**@file   presol.h
 * @brief  datastructures and methods for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRESOL_H__
#define __PRESOL_H__


typedef struct Presol PRESOL;           /**< presolver data structure */
typedef struct PresolData PRESOLDATA;   /**< presolver specific data */


/** destructor of presolver to free user data (called when SCIP is exiting)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    presol          : the presolver itself
 */
#define DECL_PRESOLFREE(x) RETCODE x (SCIP* scip, PRESOL* presol)

/** initialization method of presolver (called when problem solving starts)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    presol          : the presolver itself
 */
#define DECL_PRESOLINIT(x) RETCODE x (SCIP* scip, PRESOL* presol)

/** deinitialization method of presolver (called when problem solving exits)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    presol          : the presolver itself
 */
#define DECL_PRESOLEXIT(x) RETCODE x (SCIP* scip, PRESOL* presol)

/** presolving execution method
 *
 *  input:
 *    scip            : SCIP main data structure
 *    presol          : the presolver itself
 *    nrounds         : number of presolving rounds already done
 *    nnewfixedvars   : number of variables fixed since the last call to the presolver
 *    nnewaggrvars    : number of variables aggregated since the last call to the presolver
 *    nnewchgvartypes : number of variable type changes since the last call to the presolver
 *    nnewchgbds      : number of variable bounds tightend since the last call to the presolver
 *    nnewholes       : number of domain holes added since the last call to the presolver
 *    nnewdelconss    : number of deleted constraints since the last call to the presolver
 *    nnewupgdconss   : number of upgraded constraints since the last call to the presolver
 *    nnewchgcoefs    : number of changed coefficients since the last call to the presolver
 *    nnewchgsides    : number of changed left or right hand sides since the last call to the presolver
 *
 *  input/output:
 *    nfixedvars      : pointer to total number of variables fixed of all presolvers
 *    naggrvars       : pointer to total number of variables aggregated of all presolvers
 *    nchgvartypes    : pointer to total number of variable type changes of all presolvers
 *    nchgbds         : pointer to total number of variable bounds tightend of all presolvers
 *    naddholes       : pointer to total number of domain holes added of all presolvers
 *    ndelconss       : pointer to total number of deleted constraints of all presolvers
 *    nupgdconss      : pointer to total number of upgraded constraints of all presolvers
 *    nchgcoefs       : pointer to total number of changed coefficients of all presolvers
 *    nchgsides       : pointer to total number of changed left/right hand sides of all presolvers
 */
#define DECL_PRESOLEXEC(x) RETCODE x (SCIP* scip, PRESOL* presol, int nrounds,              \
   int nnewfixedvars, int nnewaggrvars, int nnewchgvartypes, int nnewchgbds, int nnewholes, \
   int nnewdelconss, int nnewupgdconss, int nnewchgcoefs, int nnewchgsides,                 \
   int* nfixedvars, int* naggrvars, int* nchgvartypes, int* nchgbds, int* naddholes,        \
   int* ndelconss, int* nupgdconss, int* nchgcoefs, int* nchgsides)



#include "scip.h"
#include "retcode.h"



/** creates a presolver */
extern
RETCODE SCIPpresolCreate(
   PRESOL**         presol,             /**< pointer to store presolver */
   const char*      name,               /**< name of presolver */
   const char*      desc,               /**< description of presolver */
   int              priority,           /**< priority of the presolver */
   DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver */
   DECL_PRESOLINIT  ((*presolinit)),    /**< initialise presolver */
   DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialise presolver */
   DECL_PRESOLEXEC  ((*presolexec)),    /**< presolving execution method */
   PRESOLDATA*      presoldata          /**< presolver data */
   );

/** frees memory of presolver */
extern
RETCODE SCIPpresolFree(
   PRESOL**         presol,             /**< pointer to presolver data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes presolver */
extern
RETCODE SCIPpresolInit(
   PRESOL*          presol,             /**< presolver */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** deinitializes presolver */
extern
RETCODE SCIPpresolExit(
   PRESOL*          presol,             /**< presolver */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** executes presolver */
extern
RETCODE SCIPpresolExec(
   PRESOL*          presol,             /**< presolver */
   SCIP*            scip,               /**< SCIP data structure */   
   int              nrounds,            /**< number of presolving rounds already done */
   int*             nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*             naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*             nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*             nchgbds,            /**< pointer to total number of variable bounds tightend of all presolvers */
   int*             naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*             ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*             nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*             nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*             nchgsides           /**< pointer to total number of changed left/right hand sides of all presolvers */
   );

/** gets name of presolver */
extern
const char* SCIPpresolGetName(
   PRESOL*          presol              /**< presolver */
   );

/** gets priority of presolver */
extern
int SCIPpresolGetPriority(
   PRESOL*          presol              /**< presolver */
   );

/** gets user data of presolver */
extern
PRESOLDATA* SCIPpresolGetData(
   PRESOL*          presol              /**< presolver */
   );

/** sets user data of presolver; user has to free old data in advance! */
extern
void SCIPpresolSetData(
   PRESOL*          presol,             /**< presolver */
   PRESOLDATA*      presoldata          /**< new presolver user data */
   );

/** is presolver initialized? */
extern
Bool SCIPpresolIsInitialized(
   PRESOL*          presol              /**< presolver */
   );


#endif
