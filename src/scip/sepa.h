/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa.h
 * @brief  methods and datastructures for separating cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SEPA_H__
#define __SEPA_H__


typedef struct Sepa SEPA;             /**< storage for sepad variables */


#include "def.h"
#include "retcode.h"


extern
RETCODE SCIPsepaCreate(                 /**< creates separation storage */
   SEPA**           sepa                /**< pointer to store separation storage */
   );

extern
RETCODE SCIPsepaFree(                   /**< frees separation storage */
   SEPA**           sepa                /**< pointer to store separation storage */
   );

extern
RETCODE SCIPsepaAddCut(                 /**< adds cut to separation storage and captures it */
   SEPA*            sepa,               /**< separation storage */
   const SET*       set,                /**< global SCIP settings */
   ROW*             cut,                /**< separated cut */
   Real             score,              /**< separation score of cut (the larger, the better the cut) */
   Bool             pool                /**< should the cut be used in the global cut pool? Cut must be global valid! */
   );

extern
RETCODE SCIPsepaApplyCuts(              /**< adds cuts to the LP and clears separation storage */
   SEPA*            sepa,               /**< separation storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp                  /**< LP data */
   );

extern
int SCIPsepaGetNCuts(                   /**< get number of cuts in the separation storage */
   SEPA*            sepa                /**< separation storage */
   );



#endif
