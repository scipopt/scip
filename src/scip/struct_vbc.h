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
#pragma ident "@(#) $Id: struct_vbc.h,v 1.4 2005/01/21 09:17:10 bzfpfend Exp $"

/**@file   struct_vbc.h
 * @brief  datastructures for VBC Tool output
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_VBC_H__
#define __STRUCT_VBC_H__

#include <stdio.h>

#include "def.h"
#include "type_misc.h"


/** VBC Tool data structure */
struct Vbc
{
   FILE*            file;               /**< file to store VBC information */
   HASHMAP*         nodenum;            /**< hash map for mapping nodes to node numbers */
   Longint          timestep;           /**< time step counter for non real time output */
   Bool             userealtime;        /**< should the real solving time be used instead of a time step counter? */
};


#endif
