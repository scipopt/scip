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
#pragma ident "@(#) $Id: struct_presol.h,v 1.1 2003/12/01 14:41:33 bzfpfend Exp $"

/**@file   struct_presol.h
 * @brief  datastructures for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_PRESOL_H__
#define __STRUCT_PRESOL_H__


#include "def.h"
#include "type_clock.h"
#include "type_presol.h"


/** presolver */
struct Presol
{
   char*            name;               /**< name of presolver */
   char*            desc;               /**< description of presolver */
   int              priority;           /**< priority of the presolver */
   DECL_PRESOLFREE  ((*presolfree));    /**< destructor of presolver */
   DECL_PRESOLINIT  ((*presolinit));    /**< initialize presolver */
   DECL_PRESOLEXIT  ((*presolexit));    /**< deinitialize presolver */
   DECL_PRESOLEXEC  ((*presolexec));    /**< presolver execution method */
   PRESOLDATA*      presoldata;         /**< presolver data */
   CLOCK*           clock;              /**< presolving time */
   unsigned int     initialized:1;      /**< is presolver initialized? */
   int              lastnfixedvars;     /**< number of variables fixed before the last call to the presolver */
   int              lastnaggrvars;      /**< number of variables aggregated before the last call to the presolver */
   int              lastnchgvartypes;   /**< number of variable type changes before the last call to the presolver */
   int              lastnchgbds;        /**< number of variable bounds tightend before the last call to the presolver */
   int              lastnaddholes;      /**< number of domain holes added before the last call to the presolver */
   int              lastndelconss;      /**< number of deleted constraints before the last call to the presolver */
   int              lastnupgdconss;     /**< number of upgraded constraints before the last call to the presolver */
   int              lastnchgcoefs;      /**< number of changed coefficients before the last call to the presolver */
   int              lastnchgsides;      /**< number of changed left or right hand sides before the last call */
   int              nfixedvars;         /**< total number of variables fixed by this presolver */
   int              naggrvars;          /**< total number of variables aggregated by this presolver */
   int              nchgvartypes;       /**< total number of variable type changes by this presolver */
   int              nchgbds;            /**< total number of variable bounds tightend by this presolver */
   int              naddholes;          /**< total number of domain holes added by this presolver */
   int              ndelconss;          /**< total number of deleted constraints by this presolver */
   int              nupgdconss;         /**< total number of upgraded constraints by this presolver */
   int              nchgcoefs;          /**< total number of changed coefficients by this presolver */
   int              nchgsides;          /**< total number of changed left or right hand sides by this presolver */
};


#endif
