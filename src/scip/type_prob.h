/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_prob.h,v 1.4 2004/04/27 15:50:06 bzfpfend Exp $"

/**@file   type_prob.h
 * @brief  type definitions for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_PROB_H__
#define __TYPE_PROB_H__


/** objective sense: minimization or maximization */
enum Objsense
{
   SCIP_OBJSENSE_MAXIMIZE = -1,         /**< maximization of objective function */
   SCIP_OBJSENSE_MINIMIZE = +1          /**< minimization of objective function (the default) */
};
typedef enum Objsense OBJSENSE;

typedef struct Prob PROB;               /**< main problem to solve */
typedef struct ProbData PROBDATA;       /**< user problem data set by the reader */


/** frees user data of original problem (called when the original problem is freed)
 *
 *  This method should free the user data of the original problem.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    probdata        : pointer to the user problem data to free
 */
#define DECL_PROBDELORIG(x) RETCODE x (SCIP* scip, PROBDATA** probdata)

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed)
 *
 *  Because the original problem and the user data of the original problem should not be
 *  modified during the solving process, a transformed problem is created as a copy of
 *  the original problem. If the user problem data is never modified during the solving
 *  process anyways, it is enough to simple copy the user data's pointer. This is the
 *  default implementation, which is used when a NULL is given as PROBTRANS method.
 *  If the user data may be modified during the solving process (e.g. during preprocessing),
 *  the PROBTRANS method must be given and has to copy the user problem data to a different
 *  memory location.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    sourcedata      : source problem data to transform
 *    targetdata      : pointer to store created transformed problem data
 */
#define DECL_PROBTRANS(x) RETCODE x (SCIP* scip, PROBDATA* sourcedata, PROBDATA** targetdata)

/** frees user data of transformed problem (called when the transformed problem is freed)
 *
 *  This method has to be implemented, if the PROBTRANS method is not a simple pointer
 *  copy operation like in the default PROBTRANS implementation. It should free the
 *  user data of the transformed problem, that was created in the PROBTRANS method.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    probdata        : pointer to the user problem data to free
 */
#define DECL_PROBDELTRANS(x) RETCODE x (SCIP* scip, PROBDATA** probdata)

/** solving process initialization method of transformed data (called before the branch and bound process begins)
 *
 *  This method is called before the branch and bound process begins and can be used to initialize user problem
 *  data that depends for example on the number of active problem variables, because these are now fixed.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    probdata        : user problem data
 */
#define DECL_PROBINITSOL(x) RETCODE x (SCIP* scip, PROBDATA* probdata)

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed)
 *
 *  This method is called before the branch and bound data is freed and should be used to free all data that
 *  was allocated in the solving process initialization method. The user has to make sure, that all LP rows associated
 *  to the transformed user problem data are released.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    probdata        : user problem data
 */
#define DECL_PROBEXITSOL(x) RETCODE x (SCIP* scip, PROBDATA* probdata)




#include "def.h"
#include "type_retcode.h"
#include "type_scip.h"


#endif
