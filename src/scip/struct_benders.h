/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_benders.h
 * @ingroup INTERNALAPI
 * @brief  data structures required for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_BENDERS_H__
#define __SCIP_STRUCT_BENDERS_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_benders.h"

#ifdef __cplusplus
extern "C" {
#endif

/** variable benderss data */
struct SCIP_Benders
{
   char*                 name;               /**< name of Benders' decomposition */
   char*                 desc;               /**< description of Benders' decomposition */
   SCIP_DECL_BENDERSCOPY ((*benderscopy));   /**< copy method of benders or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree));   /**< destructor of variable benders */
   SCIP_DECL_BENDERSINIT ((*bendersinit));   /**< initialize variable benders */
   SCIP_DECL_BENDERSEXIT ((*bendersexit));   /**< deinitialize variable benders */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol));/**< solving process initialization method of variable benders */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol));/**< solving process deinitialization method of variable benders */
   SCIP_DECL_BENDERSGETMASTERVAR((*bendersgetmastervar));/**< returns the master variable for the given subproblem variable*/
   SCIP_DECL_BENDERSEXEC ((*bendersexec));   /**< executes the solve method for the Benders' decomposition subproblems */
   SCIP_BENDERSDATA*     bendersdata;        /**< variable benderss local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this benders for the next stages */
   SCIP_CLOCK*           bendersclock;       /**< benders execution time */
   int                   priority;           /**< priority of the Benders' decomposition */
   int                   ncalls;             /**< number of times, this benders was called */
   int                   noptcutsfound;      /**< number of optimality cuts found by the Benders' decomposition */
   int                   nfeascutsfound;     /**< number of feasibility cuts found by the Benders' decomposition */
   SCIP_Bool             initialized;        /**< is Benders' decomposition initialized? */

   /* the subproblem information */
   int                   nsubproblems;       /**< number of subproblems */
};

#ifdef __cplusplus
}
#endif

#endif
