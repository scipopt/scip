/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_vbc.h
 * @brief  datastructures for VBC Tool output
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_VBC_H__
#define __SCIP_STRUCT_VBC_H__

#include <stdio.h>

#include "scip/def.h"
#include "scip/type_vbc.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** VBC Tool data structure */
struct SCIP_Vbc
{
   FILE*                 file;               /**< file to store VBC information */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< message handler to use */
   SCIP_HASHMAP*         nodenum;            /**< hash map for mapping nodes to node numbers */
   SCIP_Longint          timestep;           /**< time step counter for non real time output */
   SCIP_NODE*            lastnode;           /**< last node that was colored */
   SCIP_VBCCOLOR         lastcolor;          /**< last color that was used */
   SCIP_Bool             userealtime;        /**< should the real solving time be used instead of a time step counter? */
};

#ifdef __cplusplus
}
#endif

#endif
