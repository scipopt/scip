/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_prob.h,v 1.2 2003/12/08 11:51:05 bzfpfend Exp $"

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
 *  (called when problem solving starts)
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



#include "def.h"
#include "type_retcode.h"
#include "type_scip.h"


#endif
