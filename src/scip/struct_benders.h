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
#include "scip/type_benderscut.h"

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
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre));/**< presolving initialization method for Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre));/**< presolving deinitialization method for Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol));/**< solving process initialization method of variable benders */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol));/**< solving process deinitialization method of variable benders */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)); /**< returns the corresponding variable from the master or subproblem */
   SCIP_DECL_BENDERSEXEC ((*bendersexec));   /**< executes the solve method for Benders' decomposition */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub));/**< creates the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub));/**< the solving method for the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve));/**< called after the subproblems are solved. */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub));/**< the freeing method for the Benders' decomposition subproblems */
   SCIP_BENDERSDATA*     bendersdata;        /**< variable benderss local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this benders for the next stages */
   SCIP_CLOCK*           bendersclock;       /**< benders execution time */
   int                   priority;           /**< priority of the Benders' decomposition */
   int                   ncalls;             /**< number of times, this benders was called */
   int                   noptcutsfound;      /**< number of optimality cuts found by the Benders' decomposition */
   int                   nfeascutsfound;     /**< number of feasibility cuts found by the Benders' decomposition */
   SCIP_Bool             initialized;        /**< is Benders' decomposition initialized? */
   SCIP_Bool             cutlp;              /**< should Benders' cuts be generated for LP solutions? */
   SCIP_Bool             cutpseudo;          /**< should Benders' cuts be generated for pseudo solutions? */
   SCIP_Bool             cutrelax;           /**< should Benders' cuts be generated for relaxation solutions? */

   /* data for the Magnanti-Wong cut strengthening */
   SCIP_Bool             usemagnantiwong;    /**< Should the Magnanti-Wong cut strengthening technique be used? */
   SCIP_Bool             computerelint;      /**< Should the relative interior point be computed? */
   int                   maxlpiterfactor;    /**< the factor for the maximum number of lp iterations. */
   SCIP_SOL*             relintsol;          /**< the relative interior point used for the Magnanti-Wong technique. */
   SCIP_SOL*             currentsol;         /**< the current solution used to fix variables in the subproblem. */
   SCIP_VAR**            mwauxiliaryvars;    /**< the auxiliary variables for the magnanti-wong method */
   SCIP_Real             updatefactor;       /**< the factor used to update the core point between iterations */
   SCIP_Bool             coreptupdated;      /**< indicates whether the core point has been updated. */

   /* the subproblem information */
   SCIP**                subproblems;        /**< the Benders' decomposition subproblems */
   SCIP_VAR**            auxiliaryvars;      /**< the auxiliary variables for the Benders' optimality cuts */
   SCIP_Real*            subprobobjval;      /**< the objective value of the subproblem in the current iteration */
   SCIP_Real*            bestsubprobobjval;  /**< the best objective value of the subproblem */
   int                   addedsubprobs;      /**< subproblems added to the Benders' decomposition data */
   int                   nsubproblems;       /**< number of subproblems */
   SCIP_Bool*            subprobislp;        /**< is the subproblem formulated as an LP? */

   /* Bender's cut information */
   SCIP_BENDERSCUT**     benderscuts;        /**< the available Benders' cut algorithms */
   int                   nbenderscuts;       /**< the number of Benders' cut algorithms */
   int                   benderscutssize;    /**< the size of the Benders' cuts algorithms array */
   SCIP_Bool             benderscutssorted;  /**< are the Benders' cuts algorithms sorted by priority */
   SCIP_Bool             benderscutsnamessorted;/**< are the Benders' cuts algorithms sorted by name */
};

#ifdef __cplusplus
}
#endif

#endif
