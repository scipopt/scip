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
#pragma ident "@(#) $Id: type_prob.h,v 1.1 2003/12/01 14:41:37 bzfpfend Exp $"

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


/** frees user problem data
 *
 *  input:
 *    scip            : SCIP main data structure
 *    probdata        : pointer to the user problem data to free
 */
#define DECL_PROBDELETE(x) RETCODE x (SCIP* scip, PROBDATA** probdata)

/** transforms user problem data into data belonging to the transformed problem
 *
 *  input:
 *    scip            : SCIP main data structure
 *    sourcedata      : source problem data to transform
 *    targetdata      : pointer to store created transformed problem data
 */
#define DECL_PROBTRANS(x) RETCODE x (SCIP* scip, PROBDATA* sourcedata, PROBDATA** targetdata)



#include "def.h"
#include "type_retcode.h"
#include "type_scip.h"


#endif
