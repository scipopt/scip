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
#pragma ident "@(#) $Id: struct_presol.h,v 1.9 2005/02/07 18:12:02 bzfpfend Exp $"

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
   DECL_PRESOLFREE  ((*presolfree));    /**< destructor of presolver to free user data (called when SCIP is exiting) */
   DECL_PRESOLINIT  ((*presolinit));    /**< initialization method of presolver (called after problem was transformed) */
   DECL_PRESOLEXIT  ((*presolexit));    /**< deinitialization method of presolver (called before transformed problem is freed) */
   DECL_PRESOLINITPRE((*presolinitpre));/**< presolving initialization method of presolver (called when presolving is about to begin) */
   DECL_PRESOLEXITPRE((*presolexitpre));/**< presolving deinitialization method of presolver (called after presolving has been finished) */
   DECL_PRESOLEXEC  ((*presolexec));    /**< execution method of presolver */
   PRESOLDATA*      presoldata;         /**< presolver data */
   CLOCK*           clock;              /**< presolving time */
   int              priority;           /**< priority of the presolver */
   int              maxrounds;          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
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
   Bool             delay;              /**< should presolver be delayed, if other presolvers found reductions? */
   Bool             wasdelayed;         /**< was the presolver delayed at the last call? */
   Bool             initialized;        /**< is presolver initialized? */
};


#endif
