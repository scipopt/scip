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
#pragma ident "@(#) $Id: struct_sepa.h,v 1.3 2004/02/04 17:27:46 bzfpfend Exp $"

/**@file   struct_sepa.h
 * @brief  datastructures for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_SEPA_H__
#define __STRUCT_SEPA_H__


#include "def.h"
#include "type_clock.h"
#include "type_sepa.h"


/** separators data */
struct Sepa
{
   char*            name;               /**< name of separator */
   char*            desc;               /**< description of separator */
   int              priority;           /**< priority of the separator */
   int              freq;               /**< frequency for calling separator */
   DECL_SEPAFREE    ((*sepafree));      /**< destructor of separator */
   DECL_SEPAINIT    ((*sepainit));      /**< initialize separator */
   DECL_SEPAEXIT    ((*sepaexit));      /**< deinitialize separator */
   DECL_SEPAEXEC    ((*sepaexec));      /**< execution method of separator */
   SEPADATA*        sepadata;           /**< separators local data */
   CLOCK*           clock;              /**< separation time */
   Longint          lastsepanode;       /**< last node where this separator was called */
   Longint          ncalls;             /**< number of times, this separator was called */
   Longint          ncutsfound;         /**< number of cutting planes found so far by this separator */
   int              ncallsatnode;       /**< number of times, this separator was called at the current node */
   int              ncutsfoundatnode;   /**< number of cutting planes found at the current node */
   Bool             initialized;        /**< is separator initialized? */
};


#endif
