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
#pragma ident "@(#) $Id: type_vbc.h,v 1.1 2004/03/22 16:03:31 bzfpfend Exp $"

/**@file   type_vbc.h
 * @brief  type definitions for VBC Tool output
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_VBC_H__
#define __TYPE_VBC_H__


/** node colors in VBC output */
enum VBCColor
{
   SCIP_VBCCOLOR_UNSOLVED = 3,          /**< color for newly created, unsolved nodes */
   SCIP_VBCCOLOR_SOLVED   = 2,          /**< color for solved nodes */
   SCIP_VBCCOLOR_SOLUTION = 5           /**< color for solved nodes, where a solution has been found */
};


typedef struct Vbc VBC;                 /**< VBC Tool data structure */


#endif
