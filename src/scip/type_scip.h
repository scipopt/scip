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
#pragma ident "@(#) $Id: type_scip.h,v 1.1 2003/12/01 14:41:37 bzfpfend Exp $"

/**@file   type_scip.h
 * @brief  type definitions for SCIP's main datastructure
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_SCIP_H__
#define __TYPE_SCIP_H__


/** SCIP operation stage */
enum Stage
{
   SCIP_STAGE_INIT       = 0,           /**< SCIP datastructures are initialized, no problem exists */
   SCIP_STAGE_PROBLEM    = 1,           /**< the problem is being created and modified */
   SCIP_STAGE_INITSOLVE  = 2,           /**< the solving process data is being initialized */
   SCIP_STAGE_PRESOLVING = 3,           /**< the problem is being presolved */
   SCIP_STAGE_PRESOLVED  = 4,           /**< the problem was presolved */
   SCIP_STAGE_SOLVING    = 5,           /**< the problem is being solved */
   SCIP_STAGE_SOLVED     = 6,           /**< the problem was solved */
   SCIP_STAGE_FREESOLVE  = 7            /**< the solving process data is being freed */
};
typedef enum Stage STAGE;

typedef struct Scip SCIP;               /**< SCIP main data structure */


#endif
