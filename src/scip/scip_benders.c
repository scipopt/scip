/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_benders.c
 * @brief  public methods for Benders decomposition
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Robert Lion Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif


#include "lpi/lpi.h"
#include "nlpi/exprinterpret.h"
#include "nlpi/nlpi.h"
#include "scip/benders.h"
#include "scip/benderscut.h"
#include "scip/branch.h"
#include "scip/branch_nodereopt.h"
#include "scip/clock.h"
#include "scip/compr.h"
#include "scip/concsolver.h"
#include "scip/concurrent.h"
#include "scip/conflict.h"
#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/cons_linear.h"
#include "scip/cutpool.h"
#include "scip/cuts.h"
#include "scip/debug.h"
#include "scip/def.h"
#include "scip/dialog.h"
#include "scip/dialog_default.h"
#include "scip/disp.h"
#include "scip/event.h"
#include "scip/heur.h"
#include "scip/heur_ofins.h"
#include "scip/heur_reoptsols.h"
#include "scip/heur_trivialnegation.h"
#include "scip/heuristics.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/interrupt.h"
#include "scip/lp.h"
#include "scip/mem.h"
#include "scip/message_default.h"
#include "scip/misc.h"
#include "scip/nlp.h"
#include "scip/nodesel.h"
#include "scip/paramset.h"
#include "scip/presol.h"
#include "scip/presolve.h"
#include "scip/pricer.h"
#include "scip/pricestore.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/prop.h"
#include "scip/reader.h"
#include "scip/relax.h"
#include "scip/reopt.h"
#include "scip/retcode.h"
#include "scip/scipbuildflags.h"
#include "scip/scipcoreplugins.h"
#include "scip/scipgithash.h"
#include "scip/sepa.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/syncstore.h"
#include "scip/table.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/visual.h"
#include "xml/xml.h"

#include "scip/scip_benders.h"
#include "scip/scip_cons.h"
#include "scip/scip_lp.h"

#include "scip/pub_benders.h"
#include "scip/pub_message.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** creates a Benders' decomposition and includes it in SCIP
 *
 *  To use the Benders' decomposition for solving a problem, it first has to be activated with a call to SCIPactivateBenders().
 *  This should be done during the problem creation stage.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note method has all Benders' decomposition callbacks as arguments and is thus changed every time a new callback is
 *        added in future releases; consider using SCIPincludeBendersBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   SCIP_Bool             cutlp,              /**< should Benders' cuts be generated for LP solutions */
   SCIP_Bool             cutpseudo,          /**< should Benders' cuts be generated for pseudo solutions */
   SCIP_Bool             cutrelax,           /**< should Benders' cuts be generated for relaxation solutions */
   SCIP_Bool             shareauxvars,       /**< should this Benders' use the highest priority Benders aux vars */
   SCIP_DECL_BENDERSCOPY ((*benderscopy)),   /**< copy method of Benders' decomposition or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree)),   /**< destructor of Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit)),   /**< initialize Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit)),   /**< deinitialize Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre)),/**< presolving initialization method for Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre)),/**< presolving deinitialization method for Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol)),/**< solving process initialization method of Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol)),/**< solving process deinitialization method of Benders' decomposition */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub)),/**< creates a Benders' decomposition subproblem */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve)),/**< the execution method of the Benders' decomposition algorithm */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex)),/**< the solving method for convex Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub)),/**< the solving method for the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve)),/**< called after the subproblems are solved. */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub)),/**< the freeing method for the Benders' decomposition subproblems */
   SCIP_BENDERSDATA*     bendersdata         /**< Benders' decomposition data */
   )
{
   SCIP_BENDERS* benders;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeBenders", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether pricer is already present */
   if( SCIPfindBenders(scip, name) != NULL )
   {
      SCIPerrorMessage("benders <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   /* Checking whether the benderssolvesub and the bendersfreesub are both implemented or both are not implemented */
   if( (benderssolvesubconvex == NULL && benderssolvesub == NULL && bendersfreesub != NULL)
      || ((benderssolvesubconvex != NULL || benderssolvesub != NULL) && bendersfreesub == NULL) )
   {
      SCIPerrorMessage("Benders' decomposition <%s> requires that if bendersFreesub%s is "
         "implemented at least one of bendersSolvesubconvex%s or bendersSolvesub%s are implemented, "
         "or if bendersFreesub%s is not implemented, then none are implented.\n", name, name, name, name, name);
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPbendersCreate(&benders, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         cutlp, cutpseudo, cutrelax, shareauxvars, benderscopy, bendersfree, bendersinit, bendersexit, bendersinitpre,
         bendersexitpre, bendersinitsol, bendersexitsol, bendersgetvar, benderscreatesub, benderspresubsolve,
         benderssolvesubconvex, benderssolvesub, benderspostsolve, bendersfreesub, bendersdata) );
   SCIP_CALL( SCIPsetIncludeBenders(scip->set, benders) );

   return SCIP_OKAY;
}

/** creates a Benders' decomposition and includes it in SCIP with all non-fundamental callbacks set to NULL
 *
 *  If needed, the non-fundamental callbacks can be added afterwards via setter functions SCIPsetBendersCopy(),
 *  SCIPsetBendersFree(), SCIPsetBendersInity(), SCIPsetBendersExit(), SCIPsetBendersInitsol(), SCIPsetBendersExitsol(),
 *  SCIPsetBendersFarkas().
 *
 *  To use the Benders' decomposition for solving a problem, it first has to be activated with a call to SCIPactivateBenders().
 *  This should be done during the problem creation stage.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeBenders() instead
 */
SCIP_RETCODE SCIPincludeBendersBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS**        bendersptr,         /**< reference to a benders, or NULL */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   SCIP_Bool             cutlp,              /**< should Benders' cuts be generated for LP solutions */
   SCIP_Bool             cutpseudo,          /**< should Benders' cuts be generated for pseudo solutions */
   SCIP_Bool             cutrelax,           /**< should Benders' cuts be generated for relaxation solutions */
   SCIP_Bool             shareauxvars,       /**< should this Benders' use the highest priority Benders aux vars */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub)),/**< creates a Benders' decomposition subproblem */
   SCIP_BENDERSDATA*     bendersdata         /**< Benders' decomposition data */
   )
{
   SCIP_BENDERS* benders;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeBendersBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether Benders' decomposition is already present */
   if( SCIPfindBenders(scip, name) != NULL )
   {
      SCIPerrorMessage("benders <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbendersCreate(&benders, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         cutlp, cutpseudo, cutrelax, shareauxvars, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, bendersgetvar,
         benderscreatesub, NULL, NULL, NULL, NULL, NULL, bendersdata) );
   SCIP_CALL( SCIPsetIncludeBenders(scip->set, benders) );

   if( bendersptr != NULL )
      *bendersptr = benders;

   return SCIP_OKAY;
}

/** sets copy method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSCOPY((*benderscopy))     /**< copy method of Benders' decomposition or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetCopy(benders, benderscopy);

   return SCIP_OKAY;
}

/** sets destructor method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSFREE((*bendersfree))     /**< destructor of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetFree(benders, bendersfree);

   return SCIP_OKAY;
}

/** sets initialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit))    /**< initialize Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetInit(benders, bendersinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit))    /**< deinitialize Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetExit(benders, bendersexit);

   return SCIP_OKAY;
}

/** sets presolving initialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersInitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre))/**< presolving initialization method of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersInitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetInitpre(benders, bendersinitpre);

   return SCIP_OKAY;
}

/** sets presolving deinitialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersExitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre))/**< presolving deinitialization method of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersExitpre", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetExitpre(benders, bendersexitpre);

   return SCIP_OKAY;
}

/** sets solving process initialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol))/**< solving process initialization method of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetInitsol(benders, bendersinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol))/**< solving process deinitialization method of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetExitsol(benders, bendersexitsol);

   return SCIP_OKAY;
}

/** sets the method called prior to solving the subproblems for benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersPresubsolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve))/**< method called prior to solving the subproblems */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersPresubsolve", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetPresubsolve(benders, benderspresubsolve);

   return SCIP_OKAY;
}

/** sets the subproblem solving and freeing methods for Benders' decomposition
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersSolveAndFreesub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex)),/**< the solving method for convex Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub)),/**< solving method for a Benders' decomposition subproblem */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub))/**< the subproblem freeing method for Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersSolveAndFreesub", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   /* Checking whether the benderssolvesub and the bendersfreesub are both implemented or both are not implemented */
   if( (benderssolvesubconvex == NULL && benderssolvesub == NULL && bendersfreesub != NULL)
      || ((benderssolvesubconvex != NULL || benderssolvesub != NULL) && bendersfreesub == NULL) )
   {
      SCIPerrorMessage("Benders' decomposition <%s> requires that if bendersFreesub%s is "
         "implemented at least one of bendersSolvesubconvex%s or bendersSolvesub%s are implemented, "
         "or if bendersFreesub%s is not implemented, then none are implented.\n", SCIPbendersGetName(benders),
         SCIPbendersGetName(benders), SCIPbendersGetName(benders), SCIPbendersGetName(benders),
         SCIPbendersGetName(benders));
      return SCIP_INVALIDCALL;
   }

   SCIPbendersSetSolvesubconvex(benders, benderssolvesubconvex);
   SCIPbendersSetSolvesub(benders, benderssolvesub);
   SCIPbendersSetFreesub(benders, bendersfreesub);

   return SCIP_OKAY;
}

/** sets the post solving methods for benders
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBendersPostsolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve))/**< solving process deinitialization method of Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBendersPostsolve", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benders != NULL);

   SCIPbendersSetPostsolve(benders, benderspostsolve);

   return SCIP_OKAY;
}

/** returns the Benders' decomposition of the given name, or NULL if not existing */
SCIP_BENDERS* SCIPfindBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of Benders' decomposition */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindBenders(scip->set, name);
}

/** returns the array of currently available Benders' decomposition; active Benders' decomposition are in the first
 * slots of the array
 */
SCIP_BENDERS** SCIPgetBenders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortBenders(scip->set);

   return scip->set->benders;
}

/** returns the number of currently available Benders' decomposition */
int SCIPgetNBenders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nbenders;
}

/** returns the number of currently active Benders' decomposition */
int SCIPgetNActiveBenders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nactivebenders;
}

/** activates the Benders' decomposition to be used for the current problem
 *
 *  This method should be called during the problem creation stage for all pricers that are necessary to solve
 *  the problem model.
 *
 *  @note The Benders' decompositions are automatically deactivated when the problem is freed.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPactivateBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPactivateBenders", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersActivate(benders, scip->set, nsubproblems) );

   return SCIP_OKAY;
}

/** deactivates the Benders' decomposition
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPdeactivateBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdeactivateBenders", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   SCIPbendersDeactivate(benders, scip->set);

   return SCIP_OKAY;
}

/** sets the priority of a Benders' decomposition */
void SCIPsetBendersPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   priority            /**< new priority of the Benders' decomposition */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);

   SCIPbendersSetPriority(benders, scip->set, priority);
}

/** calls the exec method of Benders' decomposition to solve the subproblems
 *
 *  The checkint flag indicates whether integer feasibility can be assumed. If it is not assumed, i.e. checkint ==
 *  FALSE, then only the convex relaxations of the subproblems are solved. If integer feasibility is assumed, i.e.
 *  checkint == TRUE, then the convex relaxations and the full CIP are solved to generate Benders' cuts and check
 *  solution feasibility.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPsolveBendersSubproblems(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL */
   SCIP_RESULT*          result,             /**< result of the pricing process */
   SCIP_Bool*            infeasible,         /**< is the master problem infeasible with respect to the Benders' cuts? */
   SCIP_Bool*            auxviol,            /**< set to TRUE only if the solution is feasible but the aux vars are violated */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint            /**< should the integer solution be checked by the subproblems */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveBendersSubproblems", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersExec(benders, scip->set, sol, result, infeasible, auxviol, type, checkint) );

   return SCIP_OKAY;
}

/** returns the master problem variable for the given subproblem variable
 *
 *  This function is used as part of the cut generation process.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPgetBendersMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_VAR*             var,                /**< the subproblem variable */
   SCIP_VAR**            mappedvar           /**< pointer to store the master variable that var is mapped to */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);
   assert(var != NULL);
   assert(mappedvar != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetBendersMasterVar", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersGetVar(benders, scip->set, var, mappedvar, -1) );

   return SCIP_OKAY;
}

/** returns the subproblem problem variable for the given master variable
 *
 *  This function is used as part of the cut generation process.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPgetBendersSubproblemVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_VAR*             var,                /**< the master variable */
   SCIP_VAR**            mappedvar,          /**< pointer to store the subproblem variable that var is mapped to */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);
   assert(var != NULL);
   assert(mappedvar != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetBendersSubproblemVar", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersGetVar(benders, scip->set, var, mappedvar, probnumber) );

   return SCIP_OKAY;
}

/** returns the number of subproblems that are stored in the given Benders' decomposition
 *
 *  @return the number of subproblems in the Benders' decomposition
 */
int SCIPgetBendersNSubproblems(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);

   return SCIPbendersGetNSubproblems(benders);
}

/** registers the Benders' decomposition subproblem with the Benders' decomposition struct
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 */
SCIP_RETCODE SCIPaddBendersSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP*                 subproblem          /**< Benders' decomposition subproblem */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);
   assert(subproblem != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddBendersSubproblem", FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersAddSubproblem(benders, subproblem) );

   return SCIP_OKAY;
}

/** calls the generic subproblem setup method for a Benders' decomposition subproblem
 *
 *  This is called if the user requires to solve the Benders' decomposition subproblem separately from the main Benders'
 *  solving loop. This could be in the case of enhancement techniques.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPsetupBendersSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   SCIP_SOL*             sol,                /**< primal solution used to setup the problem, NULL for LP solution */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetupBendersSubproblem", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersSetupSubproblem(benders, scip->set, sol, probnumber) );

   return SCIP_OKAY;
}

/** calls the solving method for a single Benders' decomposition subproblem
 *
 *  The method either calls the users solve subproblem method or calls the generic method. In the case of the generic
 *  method, the user must set up the subproblem prior to calling this method.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPsolveBendersSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL for the current LP/Pseudo solution */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_Bool             solvecip,           /**< directly solve the CIP subproblem */
   SCIP_Real*            objective           /**< the objective function value of the subproblem, can be NULL */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveBendersSubproblem", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

  SCIP_CALL( SCIPbendersSolveSubproblem(benders, scip->set, sol, probnumber, infeasible, type, solvecip, objective) );

   return SCIP_OKAY;
}

/** frees the subproblem after calling the solve subproblem method
 *
 *  This will either call the user defined free
 *  subproblem callback for Benders' decomposition or the default freeing methods. In the default case, if the
 *  subproblem is an LP, then SCIPendProbing is called. If the subproblem is a MIP, then SCIPfreeTransform is called.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPfreeBendersSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeBendersSubproblem", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersFreeSubproblem(benders, scip->set, probnumber) );

   return SCIP_OKAY;
}

/** checks the optimality of a Benders' decomposition subproblem by comparing the objective function value against the
 *  value of the corresponding auxiliary variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if requested subproblem is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcheckBendersSubproblemOptimality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL for the current LP solution */
   int                   probnumber,         /**< the number of the pricing problem */
   SCIP_Bool*            optimal             /**< flag to indicate whether the current subproblem is optimal for the master */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* check stages for both, SCIP and the requested subproblem data structure */
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcheckBendersSubproblemOptimality", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(SCIPbendersSubproblem(benders, probnumber), "SCIPcheckBendersSubproblemOptimality",
         FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersCheckSubproblemOptimality(benders, scip->set, sol, probnumber, optimal) );

   return SCIP_OKAY;
}

/** returns the value of the auxiliary variable for a given subproblem
 *
 *  @return the value of the auxiliary variable for the given subproblem
 */
SCIP_Real SCIPgetBendersAuxiliaryVarVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL for the current LP solution */
   int                   probnumber          /**< the number of the pricing problem */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return SCIPbendersGetAuxiliaryVarVal(benders, scip->set, sol, probnumber);
}

/** solves an independent subproblem to identify its lower bound and updates the lower bound of the corresponding
 *  auxiliary variable
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPcomputeBendersSubproblemLowerbound(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem to be evaluated */
   SCIP_Real*            lowerbound,         /**< the lowerbound for the subproblem */
   SCIP_Bool*            infeasible          /**< was the subproblem found to be infeasible? */
   )
{
   assert(scip != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcomputeBendersSubproblemLowerbound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersComputeSubproblemLowerbound(benders, scip->set, probnumber, lowerbound, infeasible) );

   return SCIP_OKAY;
}

/** Merges a subproblem into the master problem. This process just adds a copy of the subproblem variables and
 *  constraints to the master problem, but keeps the subproblem stored in the Benders' decomposition data structure.
 *  The reason for keeping the subproblem available is for when it is queried for solutions after the problem is solved.
 *
 *  Once the subproblem is merged into the master problem, then the subproblem is flagged as disabled. This means that
 *  it will not be solved in the subsequent subproblem solving loops.
 *
 *  The associated auxiliary variables are kept in the master problem. The objective function of the merged subproblem
 *  is added as an underestimator constraint.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPmergeBendersSubproblemIntoMaster(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of subproblem variables corresponding
                                              *   to the newly created master variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of subproblem constraints to the
                                                  corresponding newly created constraints, or NULL */
   int                   probnumber          /**< the number of the subproblem that will be merged into the master problem*/
   )
{
   assert(scip != NULL);
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPgetBendersNSubproblems(scip, benders));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPmergeBendersSubproblemIntoMaster", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbendersMergeSubproblemIntoMaster(benders, scip->set, varmap, consmap, probnumber) );

   return SCIP_OKAY;
}

/** creates a Benders' cut algorithms and includes it in the associated Benders' decomposition
 *
 *  This should be called from the SCIPincludeBendersXyz for the associated Benders' decomposition. It is only possible
 *  to include a Benders' cut algorithm if a Benders' decomposition has already been included
 *  This should be done during the problem creation stage.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note method has all Benders' decomposition callbacks as arguments and is thus changed every time a new callback is
 *        added in future releases; consider using SCIPincludeBendersBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeBenderscut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   const char*           name,               /**< name of Benders' decomposition cuts */
   const char*           desc,               /**< description of Benders' decomposition cuts */
   int                   priority,           /**< priority of the Benders' decomposition cuts */
   SCIP_Bool             islpcut,            /**< indicates whether the cut is generated from the LP solution */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy)),/**< copy method of Benders' decomposition cuts or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree)),/**< destructor of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit)),/**< initialize Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit)),/**< deinitialize Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol)),/**< solving process initialization method of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol)),/**< solving process deinitialization method of Benders' decomposition cuts */
   SCIP_DECL_BENDERSCUTEXEC((*benderscutexec)),/**< execution method of Benders' decomposition cuts */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< Benders' decomposition cuts data */
   )
{
   SCIP_BENDERSCUT* benderscut;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeBenderscut", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether pricer is already present */
   if( SCIPfindBenderscut(benders, name) != NULL )
   {
      SCIPerrorMessage("benders <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbenderscutCreate(benders, &benderscut, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         islpcut, benderscutcopy, benderscutfree, benderscutinit, benderscutexit,
         benderscutinitsol, benderscutexitsol, benderscutexec, benderscutdata) );
   SCIP_CALL( SCIPbendersIncludeBenderscut(benders, scip->set, benderscut) );

   return SCIP_OKAY;
}

/** creates a Benders' cut and includes it an associated Benders' decomposition with all non-fundamental callbacks set to NULL
 *
 *  If needed, the non-fundamental callbacks can be added afterwards via setter functions SCIPsetBenderscutCopy(),
 *  SCIPsetBenderscutFree(), SCIPsetBenderscutInit(), SCIPsetBenderscutExit(), SCIPsetBenderscutInitsol(),
 *  SCIPsetBenderscutExitsol().
 *
 *  This should be done during the problem creation stage.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeBenders() instead
 */
SCIP_RETCODE SCIPincludeBenderscutBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSCUT**     benderscutptr,      /**< reference to a Benders' decomposition cut, or NULL */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   SCIP_Bool             islpcut,            /**< indicates whether the cut is generated from the LP solution */
   SCIP_DECL_BENDERSCUTEXEC((*benderscutexec)),/**< the execution method of the Benders' cut algorithm */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< Benders' cut data */
   )
{
   SCIP_BENDERSCUT* benderscut;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeBenderscutBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether Benders' decomposition cut is already present */
   if( SCIPfindBenderscut(benders, name) != NULL )
   {
      SCIPerrorMessage("Benders' cut <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbenderscutCreate(benders, &benderscut, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc,
         priority, islpcut, NULL, NULL, NULL, NULL, NULL, NULL, benderscutexec, benderscutdata) );
   SCIP_CALL( SCIPbendersIncludeBenderscut(benders, scip->set, benderscut) );

   if( benderscutptr != NULL )
      *benderscutptr = benderscut;

   return SCIP_OKAY;
}

/** sets copy method of Benders' decomposition cut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy))/**< copy method of benderscut or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetCopy(benderscut, benderscutcopy);

   return SCIP_OKAY;
}

/** sets destructor method of benderscut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree))/**< destructor of benderscut */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetFree(benderscut, benderscutfree);

   return SCIP_OKAY;
}

/** sets initialization method of benderscut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit))/**< initialize benderscut */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetInit(benderscut, benderscutinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of benderscut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit))/**< deinitialize benderscut */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetExit(benderscut, benderscutexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of benderscut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol))/**< solving process initialization method of benderscut */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetInitsol(benderscut, benderscutinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of benderscut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol))/**< solving process deinitialization method of benderscut */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(benderscut != NULL);

   SCIPbenderscutSetExitsol(benderscut, benderscutexitsol);

   return SCIP_OKAY;
}

/** sets the priority of a Benders' decomposition cut algorithm
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetBenderscutPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< benderscut */
   int                   priority            /**< new priority of the Benders' decomposition */
   )
{
   SCIP_BENDERS** benders;
   int nbenders;
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(benderscut != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetBenderscutPriority", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPbenderscutSetPriority(benderscut, priority);

   /* DIRTY: This is not a good fix */
   /* Changing the priority of one Benders' cut in a Benders' decomposition requires all Benders' cuts to be set to
    * unsorted. This is a fix that is not very nice, but it does the job */
   benders = SCIPgetBenders(scip);
   nbenders = SCIPgetNBenders(scip);
   for( i = 0; i < nbenders; i++ )
      SCIPbendersSetBenderscutsSorted(benders[i], FALSE);

   return SCIP_OKAY;
}

/** adds the generated constraint to the Benders' decomposition cut storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPstoreBenderscutCons(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_CONS*            cons                /**< the constraint to be added to the Benders' cut storage */
   )
{
   assert(scip != NULL);
   assert(benderscut != NULL);
   assert(cons != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPstoreBenderscutCons", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbenderscutStoreCons(benderscut, scip->set, cons) );

   /* capturing the stored constraint */
   SCIP_CALL( SCIPcaptureCons(scip, cons) );

   return SCIP_OKAY;
}

/** adds the generated cuts to the Benders' cut storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPstoreBenderscutCut(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cuts */
   SCIP_ROW*             cut                 /**< the cut to be added to the Benders' cut storage */
   )
{
   assert(scip != NULL);
   assert(benderscut != NULL);
   assert(cut != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPstoreBenderscutCut", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbenderscutStoreCut(benderscut, scip->set, cut) );

   /* capturing the row */
   SCIP_CALL( SCIPcaptureRow(scip, cut) );

   return SCIP_OKAY;
}
