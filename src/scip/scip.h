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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip.h
 * @ingroup PUBLICCOREAPI
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_H__
#define __SCIP_SCIP_H__


#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_clock.h"
#include "scip/type_misc.h"
#include "scip/type_timing.h"
#include "scip/type_paramset.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_nlp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_scip.h"

#include "scip/type_bandit.h"
#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_cons.h"
#include "scip/type_dialog.h"
#include "scip/type_disp.h"
#include "scip/type_heur.h"
#include "scip/type_compr.h"
#include "scip/type_history.h"
#include "scip/type_nodesel.h"
#include "scip/type_presol.h"
#include "scip/type_pricer.h"
#include "scip/type_reader.h"
#include "scip/type_relax.h"
#include "scip/type_sepa.h"
#include "scip/type_table.h"
#include "scip/type_prop.h"
#include "nlpi/type_nlpi.h"
#include "scip/type_concsolver.h"
#include "scip/type_syncstore.h"
#include "scip/type_benders.h"
#include "scip/type_benderscut.h"

/* include public interfaces, s.t. the user only needs to include scip.h */
#include "scip/pub_bandit.h"
#include "scip/pub_benders.h"
#include "scip/pub_benderscut.h"
#include "scip/pub_branch.h"
#include "scip/pub_conflict.h"
#include "scip/pub_cons.h"
#include "scip/pub_cutpool.h"
#include "scip/pub_dialog.h"
#include "scip/pub_disp.h"
#include "scip/pub_event.h"
#include "scip/pub_fileio.h"
#include "scip/pub_heur.h"
#include "scip/pub_compr.h"
#include "scip/pub_history.h"
#include "scip/pub_implics.h"
#include "scip/pub_lp.h"
#include "scip/pub_nlp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_nodesel.h"
#include "scip/pub_paramset.h"
#include "scip/pub_presol.h"
#include "scip/pub_pricer.h"
#include "scip/pub_reader.h"
#include "scip/pub_relax.h"
#include "scip/pub_sepa.h"
#include "scip/pub_prop.h"
#include "scip/pub_sol.h"
#include "scip/pub_table.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "lpi/lpi.h"
#include "nlpi/pub_expr.h"

/* include global presolving, cuts, and heuristics methods */
#include "scip/presolve.h"
#include "scip/cuts.h"
#include "scip/heuristics.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/var.h"
#include "scip/cons.h"
#include "scip/solve.h"
#include "scip/debug.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * miscellaneous methods
 */





/*
 * general SCIP methods
 */


/*
 * debug solution methods
 */



/*
 * message output methods
 */





/*
 * SCIP copy methods
 */


/*
 * parameter settings
 */



/*
 * SCIP user functionality methods: managing plugins
 */








/**@addtogroup PublicSeparatorMethods
 *
 * @{
 */

/** creates a separator and includes it in SCIP.
 *
 *  @note method has all separator callbacks as arguments and is thus changed every time a new
 *        callback is added
 *        in future releases; consider using SCIPincludeSepaBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
EXTERN
SCIP_RETCODE SCIPincludeSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of separator */
   const char*           desc,               /**< description of separator */
   int                   priority,           /**< priority of separator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling separator */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   SCIP_Bool             usessubscip,        /**< does the separator use a secondary SCIP instance? */
   SCIP_Bool             delay,              /**< should separator be delayed, if other separators found cuts? */
   SCIP_DECL_SEPACOPY    ((*sepacopy)),      /**< copy method of separator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   SCIP_DECL_SEPAINIT    ((*sepainit)),      /**< initialize separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialize separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol)),   /**< solving process initialization method of separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol)),   /**< solving process deinitialization method of separator */
   SCIP_DECL_SEPAEXECLP  ((*sepaexeclp)),    /**< LP solution separation method of separator */
   SCIP_DECL_SEPAEXECSOL ((*sepaexecsol)),   /**< arbitrary primal solution separation method of separator */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   );

/** creates a separator and includes it in SCIP with its most fundamental callbacks. All non-fundamental
 *  (or optional) callbacks as, e.g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetSepaInit(), SCIPsetSepaFree(),
 *  SCIPsetSepaInitsol(), SCIPsetSepaExitsol(), SCIPsetSepaCopy(), SCIPsetExit().
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeSepa() instead
 */
EXTERN
SCIP_RETCODE SCIPincludeSepaBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA**           sepa,               /**< reference to a separator, or NULL */
   const char*           name,               /**< name of separator */
   const char*           desc,               /**< description of separator */
   int                   priority,           /**< priority of separator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling separator */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   SCIP_Bool             usessubscip,        /**< does the separator use a secondary SCIP instance? */
   SCIP_Bool             delay,              /**< should separator be delayed, if other separators found cuts? */
   SCIP_DECL_SEPAEXECLP  ((*sepaexeclp)),    /**< LP solution separation method of separator */
   SCIP_DECL_SEPAEXECSOL ((*sepaexecsol)),   /**< arbitrary primal solution separation method of separator */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   );

/** sets copy method of separator */
EXTERN
SCIP_RETCODE SCIPsetSepaCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPACOPY    ((*sepacopy))       /**< copy method of separator or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of separator */
EXTERN
SCIP_RETCODE SCIPsetSepaFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAFREE    ((*sepafree))       /**< destructor of separator */
   );

/** sets initialization method of separator */
EXTERN
SCIP_RETCODE SCIPsetSepaInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAINIT    ((*sepainit))       /**< initialize separator */
   );

/** sets deinitialization method of separator */
EXTERN
SCIP_RETCODE SCIPsetSepaExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit))       /**< deinitialize separator */
   );

/** sets solving process initialization method of separator */
EXTERN
SCIP_RETCODE SCIPsetSepaInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol))    /**< solving process initialization method of separator */
   );

/** sets solving process deinitialization method of separator */
EXTERN
SCIP_RETCODE SCIPsetSepaExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol))    /**< solving process deinitialization method of separator */
   );

/** returns the separator of the given name, or NULL if not existing */
EXTERN
SCIP_SEPA* SCIPfindSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of separator */
   );

/** returns the array of currently available separators */
EXTERN
SCIP_SEPA** SCIPgetSepas(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available separators */
EXTERN
int SCIPgetNSepas(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a separator */
EXTERN
SCIP_RETCODE SCIPsetSepaPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   int                   priority            /**< new priority of the separator */
   );

/** gets value of minimal efficacy for a cut to enter the LP
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @return value of "separating/minefficacyroot" if at root node, otherwise value of "separating/minefficacy"
 */
EXTERN
SCIP_Real SCIPgetSepaMinEfficacy(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef NDEBUG
#define SCIPgetSepaMinEfficacy(scip)         (SCIPtreeGetCurrentDepth((scip)->tree) == 0 ? (scip)->set->sepa_minefficacyroot : (scip)->set->sepa_minefficacy)
#endif

/* @} */











/**@addtogroup PublicTableMethods
 *
 * @{
 */

/** creates a statistics table and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeTable(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of statistics table */
   const char*           desc,               /**< description of statistics table */
   SCIP_Bool             active,             /**< should the table be activated by default? */
   SCIP_DECL_TABLECOPY   ((*tablecopy)),     /**< copy method of statistics table or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_TABLEFREE   ((*tablefree)),     /**< destructor of statistics table */
   SCIP_DECL_TABLEINIT   ((*tableinit)),     /**< initialize statistics table */
   SCIP_DECL_TABLEEXIT   ((*tableexit)),     /**< deinitialize statistics table */
   SCIP_DECL_TABLEINITSOL ((*tableinitsol)), /**< solving process initialization method of statistics table */
   SCIP_DECL_TABLEEXITSOL ((*tableexitsol)), /**< solving process deinitialization method of statistics table */
   SCIP_DECL_TABLEOUTPUT ((*tableoutput)),   /**< output method */
   SCIP_TABLEDATA*       tabledata,          /**< statistics table data */
   int                   position,           /**< position of statistics table */
   SCIP_STAGE            earlieststage       /**< output of the statistics table is only printed from this stage onwards */
   );

/** returns the statistics table of the given name, or NULL if not existing */
EXTERN
SCIP_TABLE* SCIPfindTable(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of statistics table */
   );

/** returns the array of currently available statistics tables */
EXTERN
SCIP_TABLE** SCIPgetTables(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available statistics tables */
EXTERN
int SCIPgetNTables(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* @} */



/*
 * user interactive dialog methods
 */


/*
 * global problem methods
 */





/*
 * local subproblem methods
 */





/*
 * solve methods
 */

/**@addtogroup PublicSolveMethods
 *
 * @{
 */

/** initializes solving data structures and transforms problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
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
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post When calling this method in the \ref SCIP_STAGE_PROBLEM stage, the \SCIP stage is changed to \ref
 *        SCIP_STAGE_TRANSFORMED; otherwise, the stage is not changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPtransformProb(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms and presolves problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages:
 *        - \ref SCIP_STAGE_PRESOLVING if the presolving process was interrupted
 *        - \ref SCIP_STAGE_PRESOLVED if the presolving process was finished and did not solve the problem
 *        - \ref SCIP_STAGE_SOLVED if the problem was solved during presolving
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPpresolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms, presolves, and solves problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages depending on if and when the solution
 *        process was interrupted:

 *        - \ref SCIP_STAGE_PRESOLVING if the solution process was interrupted during presolving
 *        - \ref SCIP_STAGE_SOLVING if the solution process was interrupted during the tree search
 *        - \ref SCIP_STAGE_SOLVED if the solving process was not interrupted
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPsolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms, presolves, and solves problem using additional solvers which emphasize on
 *  finding solutions.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *        - \ref SCIP_STAGE_PRESOLVING if the solution process was interrupted during presolving
 *        - \ref SCIP_STAGE_SOLVING if the solution process was interrupted during the tree search
 *        - \ref SCIP_STAGE_SOLVED if the solving process was not interrupted
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @deprecated Please use SCIPsolveConcurrent() instead.
 */
EXTERN
SCIP_RETCODE SCIPsolveParallel(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms, presolves, and solves problem using additional solvers which emphasize on
 *  finding solutions.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *        - \ref SCIP_STAGE_PRESOLVING if the solution process was interrupted during presolving
 *        - \ref SCIP_STAGE_SOLVING if the solution process was interrupted during the tree search
 *        - \ref SCIP_STAGE_SOLVED if the solving process was not interrupted
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPsolveConcurrent(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** frees branch and bound tree and all solution process data; statistics, presolving data and transformed problem is
 *  preserved
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post If this method is called in \SCIP stage \ref SCIP_STAGE_INIT or \ref SCIP_STAGE_PROBLEM, the stage of
 *        \SCIP is not changed; otherwise, the \SCIP stage is changed to \ref SCIP_STAGE_TRANSFORMED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPfreeSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             restart             /**< should certain data be preserved for improved restarting? */
   );

/** frees branch and bound tree and all solution process data; statistics, presolving data and transformed problem is
 *  preserved
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post If this method is called in \SCIP stage \ref SCIP_STAGE_INIT or \ref SCIP_STAGE_PROBLEM, the stage of
 *        \SCIP is not changed; otherwise, the \SCIP stage is changed to \ref SCIP_STAGE_PRESOLVED.
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPfreeReoptSolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** frees all solution process data including presolving and transformed problem, only original problem is kept
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages:
 *        - \ref SCIP_STAGE_INIT if the method was called from \SCIP stage \ref SCIP_STAGE_INIT
 *        - \ref SCIP_STAGE_PROBLEM if the method was called from any other of the allowed stages
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPfreeTransform(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** informs \SCIP that the solving process should be interrupted as soon as possible (e.g., after the current node has
 *   been solved)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note the \SCIP stage does not get changed
 */
EXTERN
SCIP_RETCODE SCIPinterruptSolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** informs SCIP that the solving process should be restarted as soon as possible (e.g., after the current node has
 *  been solved)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note the \SCIP stage does not get changed
 */
EXTERN
SCIP_RETCODE SCIPrestartSolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** include specific heuristics and branching rules for reoptimization
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
EXTERN
SCIP_RETCODE SCIPenableReoptimization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             enable              /**< enable reoptimization (TRUE) or disable it (FALSE) */
   );

/** returns whether reoptimization is enabled or not */
EXTERN
SCIP_Bool SCIPisReoptEnabled(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the stored solutions corresponding to a given run */
EXTERN
SCIP_RETCODE SCIPgetReoptSolsRun(
   SCIP*                 scip,               /**< SCIP data structue */
   int                   run,                /**< number of the run */
   SCIP_SOL**            sols,               /**< array to store solutions */
   int                   allocmem,           /**< allocated size of the array */
   int*                  nsols               /**< number of solutions */
   );

/** mark all stored solutions as not updated */
EXTERN
void SCIPresetReoptSolMarks(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** check if the reoptimization process should be restarted
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcheckReoptRestart(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< current node of the branch and bound tree (or NULL) */
   SCIP_Bool*            restart             /**< pointer to store of the reoptimitation process should be restarted */
   );

/** save bound change based on dual information in the reoptimization tree
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPaddReoptDualBndchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             newbound,           /**< new bound of the variable */
   SCIP_Real             oldbound            /**< old bound of the variable */
   );

/** returns the optimal solution of the last iteration or NULL of none exists */
EXTERN
SCIP_SOL* SCIPgetReoptLastOptSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the objective coefficent of a given variable in a previous iteration
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetReoptOldObjCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable */
   int                   run,                /**< number of the run */
   SCIP_Real*            objcoef             /**< pointer to store the objective coefficient */
   );

/** returns whether we are in the restarting phase
 *
 *  @return TRUE, if we are in the restarting phase; FALSE, otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
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
EXTERN
SCIP_Bool SCIPisInRestart(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */




/*
 * variable methods
 */

/**@addtogroup PublicVariableMethods
 *
 *@{
 */

/** creates and captures problem variable; if variable is of integral type, fractional bounds are automatically rounded;
 *  an integer variable with bounds zero and one is automatically converted into a binary variable;
 *
 *  @warning When doing column generation and the original problem is a maximization problem, notice that SCIP will
 *           transform the problem into a minimization problem by multiplying the objective function by -1.  Thus, the
 *           original objective function value of variables created during the solving process has to be multiplied by
 *           -1, too.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note the variable gets captured, hence at one point you have to release it using the method SCIPreleaseVar()
 */
EXTERN
SCIP_RETCODE SCIPcreateVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable, or NULL */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data, or NULL */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable, or NULL */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable, or NULL */
   );

/** creates and captures problem variable with optional callbacks and variable data set to NULL, which can be set
 *  afterwards using SCIPvarSetDelorigData(), SCIPvarSetTransData(),
 *  SCIPvarSetDeltransData(), SCIPvarSetCopy(), and SCIPvarSetData(); sets variable flags initial=TRUE
 *  and removable = FALSE, which can be adjusted by using SCIPvarSetInitial() and SCIPvarSetRemovable(), resp.;
 *  if variable is of integral type, fractional bounds are automatically rounded;
 *  an integer variable with bounds zero and one is automatically converted into a binary variable;
 *
 *  @warning When doing column generation and the original problem is a maximization problem, notice that SCIP will
 *           transform the problem into a minimization problem by multiplying the objective function by -1.  Thus, the
 *           original objective function value of variables created during the solving process has to be multiplied by
 *           -1, too.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note the variable gets captured, hence at one point you have to release it using the method SCIPreleaseVar()
 */
EXTERN
SCIP_RETCODE SCIPcreateVarBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype             /**< type of variable */
   );

/** outputs the variable name to the file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPwriteVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR*             var,                /**< variable to output */
   SCIP_Bool             type                /**< should the variable type be also posted */
   );

/** print the given list of variables to output stream separated by the given delimiter character;
 *
 *  i. e. the variables x1, x2, ..., xn with given delimiter ',' are written as: \<x1\>, \<x2\>, ..., \<xn\>;
 *
 *  the method SCIPparseVarsList() can parse such a string
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  @note The printing process is done via the message handler system.
 */
EXTERN
SCIP_RETCODE SCIPwriteVarsList(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR**            vars,               /**< variable array to output */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             type,               /**< should the variable type be also posted */
   char                  delimiter           /**< character which is used for delimitation */
   );

/** print the given variables and coefficients as linear sum in the following form
 *  c1 \<x1\> + c2 \<x2\>   ... + cn \<xn\>
 *
 *  This string can be parsed by the method SCIPparseVarsLinearsum().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  @note The printing process is done via the message handler system.
 */
EXTERN
SCIP_RETCODE SCIPwriteVarsLinearsum(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR**            vars,               /**< variable array to output */
   SCIP_Real*            vals,               /**< array of coefficients or NULL if all coefficients are 1.0 */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             type                /**< should the variable type be also posted */
   );

/** print the given monomials as polynomial in the following form
 *  c1 \<x11\>^e11 \<x12\>^e12 ... \<x1n\>^e1n + c2 \<x21\>^e21 \<x22\>^e22 ... + ... + cn \<xn1\>^en1 ...
 *
 *  This string can be parsed by the method SCIPparseVarsPolynomial().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  @note The printing process is done via the message handler system.
 */
EXTERN
SCIP_RETCODE SCIPwriteVarsPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR***           monomialvars,       /**< arrays with variables for each monomial */
   SCIP_Real**           monomialexps,       /**< arrays with variable exponents, or NULL if always 1.0 */
   SCIP_Real*            monomialcoefs,      /**< array with monomial coefficients */
   int*                  monomialnvars,      /**< array with number of variables for each monomial */
   int                   nmonomials,         /**< number of monomials */
   SCIP_Bool             type                /**< should the variable type be also posted */
   );

/** parses variable information (in cip format) out of a string; if the parsing process was successful a variable is
 *  created and captured; if variable is of integral type, fractional bounds are automatically rounded; an integer
 *  variable with bounds zero and one is automatically converted into a binary variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPparseVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to store the problem variable */
   const char*           str,                /**< string to parse */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   SCIP_VARDATA*         vardata,            /**< user data for this specific variable */
   char**                endptr,             /**< pointer to store the final string position if successful */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   );

/** parses the given string for a variable name and stores the variable in the corresponding pointer if such a variable
 *  exits and returns the position where the parsing stopped
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPparseVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            var,                /**< pointer to store the problem variable, or NULL if it does not exit */
   char**                endptr              /**< pointer to store the final string position if successful */
   );

/** parse the given string as variable list (here ',' is the delimiter)) (\<x1\>, \<x2\>, ..., \<xn\>) (see
 *  SCIPwriteVarsList() ); if it was successful, the pointer success is set to TRUE
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note The pointer success in only set to FALSE in the case that a variable with a parsed variable name does not exist.
 *
 *  @note If the number of (parsed) variables is greater than the available slots in the variable array, nothing happens
 *        except that the required size is stored in the corresponding integer; the reason for this approach is that we
 *        cannot reallocate memory, since we do not know how the memory has been allocated (e.g., by a C++ 'new' or SCIP
 *        memory functions).
 */
EXTERN
SCIP_RETCODE SCIPparseVarsList(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            vars,               /**< array to store the parsed variable */
   int*                  nvars,              /**< pointer to store number of parsed variables */
   int                   varssize,           /**< size of the variable array */
   int*                  requiredsize,       /**< pointer to store the required array size for the active variables */
   char**                endptr,             /**< pointer to store the final string position if successful */
   char                  delimiter,          /**< character which is used for delimitation */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successful or not */
   );

/** parse the given string as linear sum of variables and coefficients (c1 \<x1\> + c2 \<x2\> + ... + cn \<xn\>)
 *  (see SCIPwriteVarsLinearsum() ); if it was successful, the pointer success is set to TRUE
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note The pointer success in only set to FALSE in the case that a variable with a parsed variable name does not exist.
 *
 *  @note If the number of (parsed) variables is greater than the available slots in the variable array, nothing happens
 *        except that the required size is stored in the corresponding integer; the reason for this approach is that we
 *        cannot reallocate memory, since we do not know how the memory has been allocated (e.g., by a C++ 'new' or SCIP
 *        memory functions).
 */
EXTERN
SCIP_RETCODE SCIPparseVarsLinearsum(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            vars,               /**< array to store the parsed variables */
   SCIP_Real*            vals,               /**< array to store the parsed coefficients */
   int*                  nvars,              /**< pointer to store number of parsed variables */
   int                   varssize,           /**< size of the variable array */
   int*                  requiredsize,       /**< pointer to store the required array size for the active variables */
   char**                endptr,             /**< pointer to store the final string position if successful */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successful or not */
   );

/** parse the given string as polynomial of variables and coefficients
 *  (c1 \<x11\>^e11 \<x12\>^e12 ... \<x1n\>^e1n + c2 \<x21\>^e21 \<x22\>^e22 ... + ... + cn \<xn1\>^en1 ...)
 *  (see SCIPwriteVarsPolynomial()); if it was successful, the pointer success is set to TRUE
 *
 *  The user has to call SCIPfreeParseVarsPolynomialData(scip, monomialvars, monomialexps,
 *  monomialcoefs, monomialnvars, *nmonomials) short after SCIPparseVarsPolynomial to free all the
 *  allocated memory again.  Do not keep the arrays created by SCIPparseVarsPolynomial around, since
 *  they use buffer memory that is intended for short term use only.
 *
 *  Parsing is stopped at the end of string (indicated by the \\0-character) or when no more monomials
 *  are recognized.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPparseVarsPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR****          monomialvars,       /**< pointer to store arrays with variables for each monomial */
   SCIP_Real***          monomialexps,       /**< pointer to store arrays with variable exponents */
   SCIP_Real**           monomialcoefs,      /**< pointer to store array with monomial coefficients */
   int**                 monomialnvars,      /**< pointer to store array with number of variables for each monomial */
   int*                  nmonomials,         /**< pointer to store number of parsed monomials */
   char**                endptr,             /**< pointer to store the final string position if successful */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successful or not */
   );

/** frees memory allocated when parsing a polynomial from a string
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
void SCIPfreeParseVarsPolynomialData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR****          monomialvars,       /**< pointer to store arrays with variables for each monomial */
   SCIP_Real***          monomialexps,       /**< pointer to store arrays with variable exponents */
   SCIP_Real**           monomialcoefs,      /**< pointer to store array with monomial coefficients */
   int**                 monomialnvars,      /**< pointer to store array with number of variables for each monomial */
   int                   nmonomials          /**< pointer to store number of parsed monomials */
   );

/** increases usage counter of variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_RETCODE SCIPcaptureVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to capture */
   );

/** decreases usage counter of variable, if the usage pointer reaches zero the variable gets freed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  @note the pointer of the variable will be NULLed
 */
EXTERN
SCIP_RETCODE SCIPreleaseVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var                 /**< pointer to variable */
   );

/** changes the name of a variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_PROBLEM
 *
 *  @note to get the current name of a variable, use SCIPvarGetName() from pub_var.h
 */
EXTERN
SCIP_RETCODE SCIPchgVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable */
   const char*           name                /**< new name of constraint */
   );

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPtransformVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get/create transformed variable for */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable */
   );

/** gets and captures transformed variables for an array of variables;
 *  if a variable of the array is not yet transformed, a new transformed variable for this variable is created;
 *  it is possible to call this method with vars == transvars
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPtransformVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get/create transformed variables for */
   SCIP_VAR**            vars,               /**< array with variables to get/create transformed variables for */
   SCIP_VAR**            transvars           /**< array to store the transformed variables */
   );

/** gets corresponding transformed variable of a given variable;
 *  returns NULL as transvar, if transformed variable is not yet existing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPgetTransformedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get transformed variable for */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable */
   );

/** gets corresponding transformed variables for an array of variables;
 *  stores NULL in a transvars slot, if the transformed variable is not yet existing;
 *  it is possible to call this method with vars == transvars, but remember that variables that are not
 *  yet transformed will be replaced with NULL
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPgetTransformedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get transformed variables for */
   SCIP_VAR**            vars,               /**< array with variables to get transformed variables for */
   SCIP_VAR**            transvars           /**< array to store the transformed variables */
   );

/** gets negated variable x' = lb + ub - x of variable x; negated variable is created, if not yet existing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPgetNegatedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get negated variable for */
   SCIP_VAR**            negvar              /**< pointer to store the negated variable */
   );

/** gets negated variables x' = lb + ub - x of variables x; negated variables are created, if not yet existing;
 *  in difference to \ref SCIPcreateVar, the negated variable must not be released (unless captured explicitly)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPgetNegatedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get negated variables for */
   SCIP_VAR**            vars,               /**< array of variables to get negated variables for */
   SCIP_VAR**            negvars             /**< array to store the negated variables */
   );

/** gets a binary variable that is equal to the given binary variable, and that is either active, fixed, or
 *  multi-aggregated, or the negated variable of an active, fixed, or multi-aggregated variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_RETCODE SCIPgetBinvarRepresentative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to get binary representative for */
   SCIP_VAR**            repvar,             /**< pointer to store the binary representative */
   SCIP_Bool*            negated             /**< pointer to store whether the negation of an active variable was returned */
   );

/** gets binary variables that are equal to the given binary variables, and which are either active, fixed, or
 *  multi-aggregated, or the negated variables of active, fixed, or multi-aggregated variables
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_RETCODE SCIPgetBinvarRepresentatives(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of binary variables to get representatives for */
   SCIP_VAR**            vars,               /**< binary variables to get binary representatives for */
   SCIP_VAR**            repvars,            /**< array to store the binary representatives */
   SCIP_Bool*            negated             /**< array to store whether the negation of an active variable was returned */
   );

/** flattens aggregation graph of multi-aggregated variable in order to avoid exponential recursion later on
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPflattenVarAggregationGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** Transforms a given linear sum of variables, that is a_1*x_1 + ... + a_n*x_n + c into a corresponding linear sum of
 *  active variables, that is b_1*y_1 + ... + b_m*y_m + d.
 *
 *  If the number of needed active variables is greater than the available slots in the variable array, nothing happens
 *  except that the required size is stored in the corresponding variable (requiredsize). Otherwise, the active variable
 *  representation is stored in the variable array, scalar array and constant.
 *
 *  The reason for this approach is that we cannot reallocate memory, since we do not know how the memory has been
 *  allocated (e.g., by a C++ 'new' or SCIP functions).
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
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
 *
 *  @note The resulting linear sum is stored into the given variable array, scalar array, and constant. That means the
 *        given entries are overwritten.
 *
 *  @note That method can be used to convert a single variables into variable space of active variables. Therefore call
 *        the method with the linear sum 1.0*x + 0.0.
 */
EXTERN
SCIP_RETCODE SCIPgetProbvarLinearSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array x_1, ..., x_n in the linear sum which will be
                                              *   overwritten by the variable array y_1, ..., y_m in the linear sum
                                              *   w.r.t. active variables */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in linear sum which will be overwritten to the
                                              *   scalars b_1, ..., b_m in the linear sum of the active variables  */
   int*                  nvars,              /**< pointer to number of variables in the linear sum which will be
                                              *   overwritten by the number of variables in the linear sum corresponding
                                              *   to the active variables */
   int                   varssize,           /**< available slots in vars and scalars array which is needed to check if
                                              *   the array are large enough for the linear sum w.r.t. active
                                              *   variables */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c which
                                              *   will chnage to constant d in the linear sum b_1*y_1 + ... + b_m*y_m +
                                              *   d w.r.t. the active variables */
   int*                  requiredsize,       /**< pointer to store the required array size for the linear sum w.r.t. the
                                              *   active variables */
   SCIP_Bool             mergemultiples      /**< should multiple occurrences of a var be replaced by a single coeff? */
   );

/** transforms given variable, scalar and constant to the corresponding active, fixed, or
 *  multi-aggregated variable, scalar and constant; if the variable resolves to a fixed variable,
 *  "scalar" will be 0.0 and the value of the sum will be stored in "constant"; a multi-aggregation
 *  with only one active variable (this can happen due to fixings after the multi-aggregation),
 *  is treated like an aggregation; if the multi-aggregation constant is infinite, "scalar" will be 0.0
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
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
EXTERN
SCIP_RETCODE SCIPgetProbvarSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   SCIP_Real*            scalar,             /**< pointer to scalar a in sum a*x + c */
   SCIP_Real*            constant            /**< pointer to constant c in sum a*x + c */
   );

/** return for given variables all their active counterparts; all active variables will be pairwise different
 *  @note It does not hold that the first output variable is the active variable for the first input variable.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
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
EXTERN
SCIP_RETCODE SCIPgetActiveVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array with given variables and as output all active
                                              *   variables, if enough slots exist */
   int*                  nvars,              /**< number of given variables, and as output number of active variables,
                                              *   if enough slots exist */
   int                   varssize,           /**< available slots in vars array */
   int*                  requiredsize        /**< pointer to store the required array size for the active variables */
   );

/** returns the reduced costs of the variable in the current node's LP relaxation;
 *  the current node has to have a feasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 *
 *  @note The return value of this method should be used carefully if the dual feasibility check was explictely disabled.
 */
EXTERN
SCIP_Real SCIPgetVarRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get reduced costs, should be a column in current node LP */
   );

/** returns the implied reduced costs of the variable in the current node's LP relaxation;
 *  the current node has to have a feasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 *
 *  @note The return value of this method should be used carefully if the dual feasibility check was explictely disabled.
 */
EXTERN
SCIP_Real SCIPgetVarImplRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get reduced costs, should be a column in current node LP */
   SCIP_Bool             varfixing           /**< FALSE if for x == 0, TRUE for x == 1 */
   );

/** returns the Farkas coefficient of the variable in the current node's LP relaxation;
 *  the current node has to have an infeasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_Real SCIPgetVarFarkasCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get reduced costs, should be a column in current node LP */
   );

/** returns lower bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
EXTERN
SCIP_Real SCIPgetVarLbAtIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   );

/** returns upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
EXTERN
SCIP_Real SCIPgetVarUbAtIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   );

/** returns lower or upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
EXTERN
SCIP_Real SCIPgetVarBdAtIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   );

/** returns whether the binary variable was fixed at the time given by the bound change index */
EXTERN
SCIP_Bool SCIPgetVarWasFixedAtIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   );

/** gets solution value for variable in current node
 *
 *  @return solution value for variable in current node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_Real SCIPgetVarSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get solution value for */
   );

/** gets solution values of multiple variables in current node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetVarSols(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get solution value for */
   SCIP_VAR**            vars,               /**< array with variables to get value for */
   SCIP_Real*            vals                /**< array to store solution values of variables */
   );

/** sets the solution value of all variables in the global relaxation solution to zero
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPclearRelaxSolVals(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the value of the given variable in the global relaxation solution;
 *  this solution can be filled by the relaxation handlers  and can be used by heuristics and for separation;
 *  You can use SCIPclearRelaxSolVals() to set all values to zero, initially;
 *  after setting all solution values, you have to call SCIPmarkRelaxSolValid()
 *  to inform SCIP that the stored solution is valid
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This method incrementally updates the objective value of the relaxation solution. If the whole solution
 *        should be updated, using SCIPsetRelaxSolVals() instead or calling SCIPclearRelaxSolVals() before setting
 *        the first value to reset the solution and the objective value to 0 may help the numerics.
 */
EXTERN
SCIP_RETCODE SCIPsetRelaxSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to set value for */
   SCIP_Real             val                 /**< solution value of variable */
   );

/** sets the values of the given variables in the global relaxation solution and informs SCIP about the validity
 *  and whether the solution can be enforced via linear cuts;
 *  this solution can be filled by the relaxation handlers  and can be used by heuristics and for separation;
 *  the solution is automatically cleared, s.t. all other variables get value 0.0
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPsetRelaxSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to set relaxation solution value for */
   SCIP_VAR**            vars,               /**< array with variables to set value for */
   SCIP_Real*            vals,               /**< array with solution values of variables */
   SCIP_Bool             includeslp          /**< does the relaxator contain all cuts in the LP? */
   );

/** sets the values of the variables in the global relaxation solution to the values in the given primal solution
 *  and informs SCIP about the validity and whether the solution can be enforced via linear cuts;
 *  the relaxation solution can be filled by the relaxation handlers and might be used by heuristics and for separation
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPsetRelaxSolValsSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal relaxation solution */
   SCIP_Bool             includeslp          /**< does the relaxator contain all cuts in the LP? */
   );

/** returns whether the relaxation solution is valid
 *
 *  @return TRUE, if the relaxation solution is valid; FALSE, otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_Bool SCIPisRelaxSolValid(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** informs SCIP that the relaxation solution is valid and whether the relaxation can be enforced through linear cuts
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPmarkRelaxSolValid(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             includeslp          /**< does the relaxator contain all cuts in the LP? */
   );

/** informs SCIP, that the relaxation solution is invalid
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPmarkRelaxSolInvalid(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the relaxation solution value of the given variable
 *
 *  @return the relaxation solution value of the given variable
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_Real SCIPgetRelaxSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get value for */
   );

/** gets the relaxation solution objective value
 *
 *  @return the objective value of the relaxation solution
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_Real SCIPgetRelaxSolObj(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** determine which branching direction should be evaluated first by strong branching
 *
 *  @return TRUE iff strong branching should first evaluate the down child
 *
 */
EXTERN
SCIP_Bool SCIPisStrongbranchDownFirst(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to determine the branching direction on */
   );

/** start strong branching - call before any strong branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note if propagation is enabled, strong branching is not done directly on the LP, but probing nodes are created
 *        which allow to perform propagation but also creates some overhead
 */
EXTERN
SCIP_RETCODE SCIPstartStrongbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             enablepropagation   /**< should propagation be done before solving the strong branching LP? */
   );

/** end strong branching - call after any strong branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPendStrongbranch(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets strong branching information on column variable with fractional value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetVarStrongbranchFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get strong branching values for */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict,         /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   );

/** gets strong branching information with previous domain propagation on column variable
 *
 *  Before calling this method, the strong branching mode must have been activated by calling SCIPstartStrongbranch();
 *  after strong branching was done for all candidate variables, the strong branching mode must be ended by
 *  SCIPendStrongbranch(). Since this method applies domain propagation before strongbranching, propagation has to be be
 *  enabled in the SCIPstartStrongbranch() call.
 *
 *  Before solving the strong branching LP, domain propagation can be performed. The number of propagation rounds
 *  can be specified by the parameter @p maxproprounds.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @warning When using this method, LP banching candidates and solution values must be copied beforehand, because
 *           they are updated w.r.t. the strong branching LP solution.
 */
EXTERN
SCIP_RETCODE SCIPgetVarStrongbranchWithPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get strong branching values for */
   SCIP_Real             solval,             /**< value of the variable in the current LP solution */
   SCIP_Real             lpobjval,           /**< LP objective value of the current LP solution */
   int                   itlim,              /**< iteration limit for strong branchings */
   int                   maxproprounds,      /**< maximum number of propagation rounds (-1: no limit, -2: parameter
                                              *   settings) */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Longint*         ndomredsdown,       /**< pointer to store the number of domain reductions down, or NULL */
   SCIP_Longint*         ndomredsup,         /**< pointer to store the number of domain reductions up, or NULL */
   SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict,         /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Real*            newlbs,             /**< array to store valid lower bounds for all active variables, or NULL */
   SCIP_Real*            newubs              /**< array to store valid upper bounds for all active variables, or NULL */
   );

/** gets strong branching information on column variable x with integral LP solution value (val); that is, the down branch
 *  is (val -1.0) and the up brach ins (val +1.0)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note If the integral LP solution value is the lower or upper bound of the variable, the corresponding branch will be
 *        marked as infeasible. That is, the valid pointer and the infeasible pointer are set to TRUE.
 */
EXTERN
SCIP_RETCODE SCIPgetVarStrongbranchInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get strong branching values for */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict,         /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   );

/** gets strong branching information on column variables with fractional values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetVarsStrongbranchesFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables to get strong branching values for */
   int                   nvars,              /**< number of variables */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching variables down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching variables up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< array to store whether the downward branches are infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< array to store whether the upward branches are infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< array to store whether conflict constraints were created for
                                              *   infeasible downward branches, or NULL */
   SCIP_Bool*            upconflict,         /**< array to store whether conflict constraints were created for
                                              *   infeasible upward branches, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   );

/** gets strong branching information on column variables with integral values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetVarsStrongbranchesInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables to get strong branching values for */
   int                   nvars,              /**< number of variables */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching variables down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching variables up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< array to store whether the downward branches are infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< array to store whether the upward branches are infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< array to store whether conflict constraints were created for
                                              *   infeasible downward branches, or NULL */
   SCIP_Bool*            upconflict,         /**< array to store whether conflict constraints were created for
                                              *   infeasible upward branches, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   );

/** get LP solution status of last strong branching call (currently only works for strong branching with propagation) */
EXTERN
SCIP_LPSOLSTAT SCIPgetLastStrongbranchLPSolStat(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        branchdir           /**< branching direction for which LP solution status is requested */
   );

/** gets strong branching information on COLUMN variable of the last SCIPgetVarStrongbranch() call;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given variable;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPgetVarStrongbranchLast(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get last strong branching values for */
   SCIP_Real*            down,               /**< stores dual bound after branching column down, or NULL */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up, or NULL */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Real*            solval,             /**< stores LP solution value of variable at last strong branching call, or NULL */
   SCIP_Real*            lpobjval            /**< stores LP objective value at last strong branching call, or NULL */
   );

/** sets strong branching information for a column variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPsetVarStrongbranchData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to set last strong branching values for */
   SCIP_Real             lpobjval,           /**< objective value of the current LP */
   SCIP_Real             primsol,            /**< primal solution value of the column in the current LP */
   SCIP_Real             down,               /**< dual bound after branching column down */
   SCIP_Real             up,                 /**< dual bound after branching column up */
   SCIP_Bool             downvalid,          /**< is the returned down value a valid dual bound? */
   SCIP_Bool             upvalid,            /**< is the returned up value a valid dual bound? */
   SCIP_Longint          iter,               /**< total number of strong branching iterations */
   int                   itlim               /**< iteration limit applied to the strong branching call */
   );

/** rounds the current solution and tries it afterwards; if feasible, adds it to storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPtryStrongbranchLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            foundsol,           /**< stores whether solution was feasible and good enough to keep */
   SCIP_Bool*            cutoff              /**< stores whether solution was cutoff due to exceeding the cutoffbound */
   );

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given variable, or -1 if strong branching was never applied to the variable in current run
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_Longint SCIPgetVarStrongbranchNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get last strong branching node for */
   );

/** if strong branching was already applied on the variable at the current node, returns the number of LPs solved after
 *  the LP where the strong branching on this variable was applied;
 *  if strong branching was not yet applied on the variable at the current node, returns INT_MAX
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_Longint SCIPgetVarStrongbranchLPAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get strong branching LP age for */
   );

/** gets number of times, strong branching was applied in current run on the given variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
int SCIPgetVarNStrongbranchs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get last strong branching node for */
   );

/** adds given values to lock numbers of type @p locktype of variable for rounding
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
EXTERN
SCIP_RETCODE SCIPaddVarLocksType(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_LOCKTYPE         locktype,           /**< type of the variable locks */
   int                   nlocksdown,         /**< modification in number of rounding down locks */
   int                   nlocksup            /**< modification in number of rounding up locks */
   );


/** adds given values to lock numbers of variable for rounding
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note This method will always add variable locks of type model
 */
EXTERN
SCIP_RETCODE SCIPaddVarLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   nlocksdown,         /**< modification in number of rounding down locks */
   int                   nlocksup            /**< modification in number of rounding up locks */
   );


/** add locks of type @p locktype of variable with respect to the lock status of the constraint and its negation;
 *  this method should be called whenever the lock status of a variable in a constraint changes, for example if
 *  the coefficient of the variable changed its sign or if the left or right hand sides of the constraint were
 *  added or removed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
EXTERN
SCIP_RETCODE SCIPlockVarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             lockdown,           /**< should the rounding be locked in downwards direction? */
   SCIP_Bool             lockup              /**< should the rounding be locked in upwards direction? */
   );

/** remove locks of type @p locktype of variable with respect to the lock status of the constraint and its negation;
 *  this method should be called whenever the lock status of a variable in a constraint changes, for example if
 *  the coefficient of the variable changed its sign or if the left or right hand sides of the constraint were
 *  added or removed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
EXTERN
SCIP_RETCODE SCIPunlockVarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             lockdown,           /**< should the rounding be unlocked in downwards direction? */
   SCIP_Bool             lockup              /**< should the rounding be unlocked in upwards direction? */
   );

/** changes variable's objective value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 */
EXTERN
SCIP_RETCODE SCIPchgVarObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             newobj              /**< new objective value */
   );

/** adds value to variable's objective value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 */
EXTERN
SCIP_RETCODE SCIPaddVarObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             addobj              /**< additional objective value */
   );

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) lower bound value;
 *  does not change the bounds of the variable
 *
 *  @return adjusted lower bound for the given variable; the bound of the variable is not changed
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Real SCIPadjustedVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to adjust the bound for */
   SCIP_Real             lb                  /**< lower bound value to adjust */
   );

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) upper bound value;
 *  does not change the bounds of the variable
 *
 *  @return adjusted upper bound for the given variable; the bound of the variable is not changed
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Real SCIPadjustedVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to adjust the bound for */
   SCIP_Real             ub                  /**< upper bound value to adjust */
   );

/** depending on SCIP's stage, changes lower bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPchgVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** depending on SCIP's stage, changes upper bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPchgVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** changes lower bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPchgVarLbNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to change bound at, or NULL for current node */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** changes upper bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPchgVarUbNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to change bound at, or NULL for current node */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** changes global lower bound of variable; if possible, adjust bound to integral value; also tightens the local bound,
 *  if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPchgVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** changes global upper bound of variable; if possible, adjust bound to integral value; also tightens the local bound,
 *  if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPchgVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** changes lazy lower bound of the variable, this is only possible if the variable is not in the LP yet
 *
 *  lazy bounds are bounds, that are enforced by constraints and the objective function; hence, these bounds do not need
 *  to be put into the LP explicitly.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note lazy bounds are useful for branch-and-price since the corresponding variable bounds are not part of the LP
 */
EXTERN
SCIP_RETCODE SCIPchgVarLbLazy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             lazylb              /**< the lazy lower bound to be set */
   );

/** changes lazy upper bound of the variable, this is only possible if the variable is not in the LP yet
 *
 *  lazy bounds are bounds, that are enforced by constraints and the objective function; hence, these bounds do not need
 *  to be put into the LP explicitly.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note lazy bounds are useful for branch-and-price since the corresponding variable bounds are not part of the LP
 */
EXTERN
SCIP_RETCODE SCIPchgVarUbLazy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             lazyub              /**< the lazy lower bound to be set */
   );

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  doesn't store any inference information in the bound change, such that in conflict analysis, this change
 *  is treated like a branching decision
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPtightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  doesn't store any inference information in the bound change, such that in conflict analysis, this change
 *  is treated like a branching decision
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPtightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** fixes variable in preprocessing or in the current node, if the new bound is tighter (w.r.t. bound strengthening
 *  epsilon) than the current bound; if possible, adjusts bound to integral value; the given inference constraint is
 *  stored, such that the conflict analysis is able to find out the reason for the deduction of the bound change
 *
 *  @note In presolving stage when not in probing mode the variable will be fixed directly, otherwise this method
 *        changes first the lowerbound by calling SCIPinferVarLbCons and second the upperbound by calling
 *        SCIPinferVarUbCons
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPinferVarFixCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             fixedval,           /**< new value for fixation */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPinferVarLbCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPinferVarUbCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** depending on SCIP's stage, fixes binary variable in the problem, in preprocessing, or in current node;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason for the
 *  deduction of the fixing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPinferBinvarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to fix */
   SCIP_Bool             fixedval,           /**< value to fix binary variable to */
   SCIP_CONS*            infercons,          /**< constraint that deduced the fixing */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the fixing tightened the local bounds, or NULL */
   );

/** fixes variable in preprocessing or in the current node, if the new bound is tighter (w.r.t. bound strengthening
 *  epsilon) than the current bound; if possible, adjusts bound to integral value; the given inference constraint is
 *  stored, such that the conflict analysis is able to find out the reason for the deduction of the bound change
 *
 *  @note In presolving stage when not in probing mode the variable will be fixed directly, otherwise this method
 *        changes first the lowerbound by calling SCIPinferVarLbProp and second the upperbound by calling
 *        SCIPinferVarUbProp
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPinferVarFixProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             fixedval,           /**< new value for fixation */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPinferVarLbProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPinferVarUbProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** depending on SCIP's stage, fixes binary variable in the problem, in preprocessing, or in current node;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason for the
 *  deduction of the fixing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPinferBinvarProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to fix */
   SCIP_Bool             fixedval,           /**< value to fix binary variable to */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the fixing */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the fixing tightened the local bounds, or NULL */
   );

/** changes global lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current global bound; if possible, adjusts bound to integral value;
 *  also tightens the local bound, if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPtightenVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** changes global upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current global bound; if possible, adjusts bound to integral value;
 *  also tightens the local bound, if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
EXTERN
SCIP_RETCODE SCIPtightenVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** for a multi-aggregated variable, returns the global lower bound computed by adding the global bounds from all aggregation variables
 *
 *  This global bound may be tighter than the one given by SCIPvarGetLbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 *  calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetLbGlobal.
 *
 *  @return the global lower bound computed by adding the global bounds from all aggregation variables
 */
EXTERN
SCIP_Real SCIPcomputeVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

/** for a multi-aggregated variable, returns the global upper bound computed by adding the global bounds from all aggregation variables
 *
 *  This global bound may be tighter than the one given by SCIPvarGetUbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 *  calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetUbGlobal.
 *
 *  @return the global upper bound computed by adding the global bounds from all aggregation variables
 */
EXTERN
SCIP_Real SCIPcomputeVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

/** for a multi-aggregated variable, returns the local lower bound computed by adding the local bounds from all aggregation variables
 *
 *  This local bound may be tighter than the one given by SCIPvarGetLbLocal, since the latter is not updated if bounds of aggregation variables are changing
 *  calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetLbLocal.
 *
 *  @return the local lower bound computed by adding the global bounds from all aggregation variables
 */
EXTERN
SCIP_Real SCIPcomputeVarLbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

/** for a multi-aggregated variable, returns the local upper bound computed by adding the local bounds from all aggregation variables
 *
 *  This local bound may be tighter than the one given by SCIPvarGetUbLocal, since the latter is not updated if bounds of aggregation variables are changing
 *  calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetUbLocal.
 *
 *  @return the local upper bound computed by adding the global bounds from all aggregation variables
 */
EXTERN
SCIP_Real SCIPcomputeVarUbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

/** for a multi-aggregated variable, gives the global lower bound computed by adding the global bounds from all
 *  aggregation variables, this global bound may be tighter than the one given by SCIPvarGetLbGlobal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
EXTERN
SCIP_Real SCIPgetVarMultaggrLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

/** for a multi-aggregated variable, gives the global upper bound computed by adding the global bounds from all
 *  aggregation variables, this upper bound may be tighter than the one given by SCIPvarGetUbGlobal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
EXTERN
SCIP_Real SCIPgetVarMultaggrUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

/** for a multi-aggregated variable, gives the local lower bound computed by adding the local bounds from all
 *  aggregation variables, this lower bound may be tighter than the one given by SCIPvarGetLbLocal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
EXTERN
SCIP_Real SCIPgetVarMultaggrLbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

/** for a multi-aggregated variable, gives the local upper bound computed by adding the local bounds from all
 *  aggregation variables, this upper bound may be tighter than the one given by SCIPvarGetUbLocal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
EXTERN
SCIP_Real SCIPgetVarMultaggrUbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPcomputeVarLbGlobal(scip, var)  (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR ? SCIPgetVarMultaggrLbGlobal(scip, var) : SCIPvarGetLbGlobal(var))
#define SCIPcomputeVarUbGlobal(scip, var)  (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR ? SCIPgetVarMultaggrUbGlobal(scip, var) : SCIPvarGetUbGlobal(var))
#define SCIPcomputeVarLbLocal(scip, var)   (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR ? SCIPgetVarMultaggrLbLocal(scip, var)  : SCIPvarGetLbLocal(var))
#define SCIPcomputeVarUbLocal(scip, var)   (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR ? SCIPgetVarMultaggrUbLocal(scip, var)  : SCIPvarGetUbLocal(var))

#endif

/** returns solution value and index of variable lower bound that is closest to the variable's value in the given primal
 *  solution or current LP solution if no primal solution is given; returns an index of -1 if no variable lower bound is
 *  available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetVarClosestVlb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_Real*            closestvlb,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvlbidx       /**< pointer to store the index of the closest variable lower bound */
   );

/** returns solution value and index of variable upper bound that is closest to the variable's value in the given primal solution;
 *  or current LP solution if no primal solution is given; returns an index of -1 if no variable upper bound is available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetVarClosestVub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_Real*            closestvub,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvubidx       /**< pointer to store the index of the closest variable lower bound */
   );

/** informs variable x about a globally valid variable lower bound x >= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  if z is non-continuous and 1/b not too small, the corresponding valid upper/lower bound
 *  z <= (x-d)/b or z >= (x-d)/b (depending on the sign of of b) is added, too;
 *  improves the global bounds of the variable and the vlb variable if possible
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPaddVarVlb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_VAR*             vlbvar,             /**< variable z    in x >= b*z + d */
   SCIP_Real             vlbcoef,            /**< coefficient b in x >= b*z + d */
   SCIP_Real             vlbconstant,        /**< constant d    in x >= b*z + d */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   );


/** informs variable x about a globally valid variable upper bound x <= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  if z is non-continuous and 1/b not too small, the corresponding valid lower/upper bound
 *  z >= (x-d)/b or z <= (x-d)/b (depending on the sign of of b) is added, too;
 *  improves the global bounds of the variable and the vlb variable if possible
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPaddVarVub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_VAR*             vubvar,             /**< variable z    in x <= b*z + d */
   SCIP_Real             vubcoef,            /**< coefficient b in x <= b*z + d */
   SCIP_Real             vubconstant,        /**< constant d    in x <= b*z + d */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   );

/** informs binary variable x about a globally valid implication:  x == 0 or x == 1  ==>  y <= b  or  y >= b;
 *  also adds the corresponding implication or variable bound to the implied variable;
 *  if the implication is conflicting, the variable is fixed to the opposite value;
 *  if the variable is already fixed to the given value, the implication is performed immediately;
 *  if the implication is redundant with respect to the variables' global bounds, it is ignored
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPaddVarImplication(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER)
                                              *                          or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   );

/** adds a clique information to SCIP, stating that at most one of the given binary variables can be set to 1;
 *  if a variable appears twice in the same clique, the corresponding implications are performed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPaddClique(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   SCIP_Bool*            values,             /**< values of the variables in the clique; NULL to use TRUE for all vars */
   int                   nvars,              /**< number of variables in the clique */
   SCIP_Bool             isequation,         /**< is the clique an equation or an inequality? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   );

/** calculates a partition of the given set of binary variables into cliques; takes into account independent clique components
 *
 *  The algorithm performs the following steps:
 *  - recomputes connected components of the clique table, if necessary
 *  - computes a clique partition for every connected component greedily.
 *  - relabels the resulting clique partition such that it satisfies the description below
 *
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same clique;
 *  the first variable is always assigned to clique 0, and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  for each clique at most 1 variables can be set to TRUE in a feasible solution;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcalcCliquePartition(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int const             nvars,              /**< number of variables in the clique */
   int*const             cliquepartition,    /**< array of length nvars to store the clique partition */
   int*const             ncliques            /**< pointer to store the number of cliques actually contained in the partition */
   );

/** calculates a partition of the given set of binary variables into negated cliques;
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same negated clique;
 *  the first variable is always assigned to clique 0 and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  for each clique with n_c variables at least n_c-1 variables can be set to TRUE in a feasible solution;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcalcNegatedCliquePartition(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int const             nvars,              /**< number of variables in the clique */
   int*const             cliquepartition,    /**< array of length nvars to store the clique partition */
   int*const             ncliques            /**< pointer to store the number of cliques actually contained in the partition */
   );

/** force SCIP to clean up all cliques; cliques do not get automatically cleaned up after presolving. Use
 *  this method to prevent inactive variables in cliques when retrieved via SCIPgetCliques()
 *
 *  @return SCIP_OKAY if everything worked, otherwise a suitable error code is passed
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_RETCODE SCIPcleanupCliques(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            infeasible          /**< pointer to store if cleanup detected infeasibility */
   );

/** gets the number of cliques in the clique table
 *
 *  @return number of cliques in the clique table
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
int SCIPgetNCliques(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the number of cliques created so far by the cliquetable
 *
 *  @return number of cliques created so far by the cliquetable
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
int SCIPgetNCliquesCreated(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the array of cliques in the clique table
 *
 *  @return array of cliques in the clique table
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_CLIQUE** SCIPgetCliques(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether there is a clique that contains both given variable/value pairs;
 *  the variables must be active binary variables;
 *  if regardimplics is FALSE, only the cliques in the clique table are looked at;
 *  if regardimplics is TRUE, both the cliques and the implications of the implication graph are regarded
 *
 *  @return TRUE, if there is a clique that contains both variable/clique pairs; FALSE, otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note a variable with it's negated variable are NOT! in a clique
 *  @note a variable with itself are in a clique
 */
EXTERN
SCIP_Bool SCIPhaveVarsCommonClique(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_Bool             value1,             /**< value of first variable */
   SCIP_VAR*             var2,               /**< second variable */
   SCIP_Bool             value2,             /**< value of second variable */
   SCIP_Bool             regardimplics       /**< should the implication graph also be searched for a clique? */
   );

/** writes the clique graph to a gml file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note there can be duplicated arcs in the output file
 *
 *  If @p writenodeweights is true, only nodes corresponding to variables that have a fractional value and only edges
 *  between such nodes are written.
 */
EXTERN
SCIP_RETCODE SCIPwriteCliqueGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           fname,              /**< name of file */
   SCIP_Bool             writenodeweights    /**< should we write weights of nodes? */
   );

/** Removes (irrelevant) variable from all its global structures, i.e. cliques, implications and variable bounds.
 *  This is an advanced method which should be used with care.
 *
 *  @return SCIP_OKAY if everything worked, otherwise a suitable error code is passed
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_RETCODE SCIPremoveVarFromGlobalStructures(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to remove from global structures */
   );

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPchgVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             branchfactor        /**< factor to weigh variable's branching score with */
   );

/** scales the branch factor of the variable with the given value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPscaleVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             scale               /**< factor to scale variable's branching factor with */
   );

/** adds the given value to the branch factor of the variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPaddVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             addfactor           /**< value to add to the branch factor of the variable */
   );

/** sets the branch priority of the variable; variables with higher branch priority are always preferred to variables
 *  with lower priority in selection of branching variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 * @note the default branching priority is 0
 */
EXTERN
SCIP_RETCODE SCIPchgVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< branch priority of the variable */
   );

/** changes the branch priority of the variable to the given value, if it is larger than the current priority
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPupdateVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< new branch priority of the variable, if it is larger than current priority */
   );

/** adds the given value to the branch priority of the variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPaddVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   addpriority         /**< value to add to the branch priority of the variable */
   );

/** sets the branch direction of the variable (-1: prefer downwards branch, 0: automatic selection, +1: prefer upwards
 *  branch)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPchgVarBranchDirection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        branchdirection     /**< preferred branch direction of the variable (downwards, upwards, auto) */
   );

/** changes type of variable in the problem;
 *
 *  @warning This type change might change the variable array returned from SCIPgetVars() and SCIPgetVarsData();
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *
 *  @note If SCIP is already beyond the SCIP_STAGE_PROBLEM and a original variable is passed, the variable type of the
 *        corresponding transformed variable is changed; the type of the original variable does not change
 *
 *  @note If the type changes from a continuous variable to a non-continuous variable the bounds of the variable get
 *        adjusted w.r.t. to integrality information
 */
EXTERN
SCIP_RETCODE SCIPchgVarType(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_VARTYPE          vartype,            /**< new type of variable */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected (, due to
                                              *   integrality condition of the new variable type) */
   );

/** in problem creation and solving stage, both bounds of the variable are set to the given value;
 *  in presolving stage, the variable is converted into a fixed variable, and bounds are changed respectively;
 *  conversion into a fixed variable changes the vars array returned from SCIPgetVars() and SCIPgetVarsData(),
 *  and also renders arrays returned from the SCIPvarGetImpl...() methods invalid
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPfixVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to fix */
   SCIP_Real             fixedval,           /**< value to fix variable to */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            fixed               /**< pointer to store whether the fixing was performed (variable was unfixed) */
   );

/** From a given equality a*x + b*y == c, aggregates one of the variables and removes it from the set of
 *  active problem variables. This changes the vars array returned from SCIPgetVars() and SCIPgetVarsData(),
 *  and also renders the arrays returned from the SCIPvarGetImpl...() methods for the two variables invalid.
 *  In the first step, the equality is transformed into an equality with active problem variables
 *  a'*x' + b'*y' == c'. If x' == y', this leads to the detection of redundancy if a' == -b' and c' == 0,
 *  of infeasibility, if a' == -b' and c' != 0, or to a variable fixing x' == c'/(a'+b') (and possible
 *  infeasibility) otherwise.
 *  In the second step, the variable to be aggregated is chosen among x' and y', prefering a less strict variable
 *  type as aggregation variable (i.e. continuous variables are preferred over implicit integers, implicit integers
 *  over integers, and integers over binaries). If none of the variables is continuous, it is tried to find an integer
 *  aggregation (i.e. integral coefficients a'' and b'', such that a''*x' + b''*y' == c''). This can lead to
 *  the detection of infeasibility (e.g. if c'' is fractional), or to a rejection of the aggregation (denoted by
 *  aggregated == FALSE), if the resulting integer coefficients are too large and thus numerically instable.
 *
 *  The output flags have the following meaning:
 *  - infeasible: the problem is infeasible
 *  - redundant:  the equality can be deleted from the constraint set
 *  - aggregated: the aggregation was successfully performed (the variables were not aggregated before)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_PRESOLVING
 */
EXTERN
SCIP_RETCODE SCIPaggregateVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             varx,               /**< variable x in equality a*x + b*y == c */
   SCIP_VAR*             vary,               /**< variable y in equality a*x + b*y == c */
   SCIP_Real             scalarx,            /**< multiplier a in equality a*x + b*y == c */
   SCIP_Real             scalary,            /**< multiplier b in equality a*x + b*y == c */
   SCIP_Real             rhs,                /**< right hand side c in equality a*x + b*y == c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            redundant,          /**< pointer to store whether the equality is (now) redundant */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   );

/** converts variable into multi-aggregated variable; this changes the variable array returned from
 *  SCIPgetVars() and SCIPgetVarsData();
 *
 *  @warning The integrality condition is not checked anymore on the multi-aggregated variable. You must not
 *           multi-aggregate an integer variable without being sure, that integrality on the aggregation variables
 *           implies integrality on the aggregated variable.
 *
 *  The output flags have the following meaning:
 *  - infeasible: the problem is infeasible
 *  - aggregated: the aggregation was successfully performed (the variables were not aggregated before)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_PRESOLVING
 */
EXTERN
SCIP_RETCODE SCIPmultiaggregateVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable x to aggregate */
   int                   naggvars,           /**< number n of variables in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_VAR**            aggvars,            /**< variables y_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real*            scalars,            /**< multipliers a_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real             constant,           /**< constant shift c in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   );

/** returns whether aggregation of variables is not allowed */
EXTERN
SCIP_Bool SCIPdoNotAggr(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether multi-aggregation is disabled */
EXTERN
SCIP_Bool SCIPdoNotMultaggr(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether variable is not allowed to be multi-aggregated */
EXTERN
SCIP_Bool SCIPdoNotMultaggrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable x to aggregate */
   );

/** returns whether dual reductions are allowed during propagation and presolving
 *
 *  @note A reduction is called dual, if it may discard feasible solutions, but leaves at least one optimal solution
 *        intact. Often such reductions are based on analyzing the objective function, reduced costs, and/or dual LPs.
 */
EXTERN
SCIP_Bool SCIPallowDualReds(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether propagation w.r.t. current objective is allowed */
EXTERN
SCIP_Bool SCIPallowObjProp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** modifies an initial seed value with the global shift of random seeds */
EXTERN
unsigned int SCIPinitializeRandomSeed(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          initialseedvalue    /**< initial seed value to be modified */
   );

/** marks the variable that it must not be multi-aggregated
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *
 *  @note There exists no "unmark" method since it has to be ensured that if a plugin requires that a variable is not
 *        multi-aggregated that this is will be the case.
 */
EXTERN
SCIP_RETCODE SCIPmarkDoNotMultaggrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to delete */
   );

/** enables the collection of statistics for a variable
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPenableVarHistory(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** disables the collection of any statistic for a variable
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPdisableVarHistory(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of "solvaldelta" in the
 *  variable's solution value and resulting change of "objdelta" in the in the LP's objective value;
 *  the update is ignored, if the objective value difference is infinite
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPupdateVarPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   SCIP_Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   SCIP_Real             weight              /**< weight in (0,1] of this update in pseudo cost sum */
   );

/** gets the variable's pseudo cost value for the given change of the variable's LP value
 *
 *  @return the variable's pseudo cost value for the given change of the variable's LP value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarPseudocostVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the variable's pseudo cost value for the given change of the variable's LP value,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's pseudo cost value for the given change of the variable's LP value,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarPseudocostValCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the variable's pseudo cost value for the given direction
 *
 *  @return the variable's pseudo cost value for the given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the variable's pseudo cost value for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's pseudo cost value for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarPseudocostCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction
 *
 *  @return the variable's (possible fractional) number of pseudo cost updates for the given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarPseudocostCount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarPseudocostCountCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** get pseudo cost variance of the variable, either for entire solve or only for current branch and bound run
 *
 *  @return returns the (corrected) variance of pseudo code information collected so far.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarPseudocostVariance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Bool             onlycurrentrun      /**< only for pseudo costs of current branch and bound run */
   );

/** calculates a confidence bound for this variable under the assumption of normally distributed pseudo costs
 *
 *  The confidence bound \f$ \theta \geq 0\f$ denotes the interval borders \f$ [X - \theta, \ X + \theta]\f$, which contains
 *  the true pseudo costs of the variable, i.e., the expected value of the normal distribution, with a probability
 *  of 2 * clevel - 1.
 *
 *  @return value of confidence bound for this variable
 */
EXTERN
SCIP_Real SCIPcalculatePscostConfidenceBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable in question */
   SCIP_BRANCHDIR        dir,                /**< the branching direction for the confidence bound */
   SCIP_Bool             onlycurrentrun,     /**< should only the current run be taken into account */
   SCIP_CONFIDENCELEVEL  clevel              /**< confidence level for the interval */
   );

/** check if variable pseudo-costs have a significant difference in location. The significance depends on
 *  the choice of \p clevel and on the kind of tested hypothesis. The one-sided hypothesis, which
 *  should be rejected, is that fracy * mu_y >= fracx * mu_x, where mu_y and mu_x denote the
 *  unknown location means of the underlying pseudo-cost distributions of x and y.
 *
 *  This method is applied best if variable x has a better pseudo-cost score than y. The method hypothesizes that y were actually
 *  better than x (despite the current information), meaning that y can be expected to yield branching
 *  decisions as least as good as x in the long run. If the method returns TRUE, the current history information is
 *  sufficient to safely rely on the alternative hypothesis that x yields indeed a better branching score (on average)
 *  than y.
 *
 *  @note The order of x and y matters for the one-sided hypothesis
 *
 *  @note set \p onesided to FALSE if you are not sure which variable is better. The hypothesis tested then reads
 *        fracy * mu_y == fracx * mu_x vs the alternative hypothesis fracy * mu_y != fracx * mu_x.
 *
 *  @return TRUE if the hypothesis can be safely rejected at the given confidence level
 */
EXTERN
SCIP_Bool SCIPsignificantVarPscostDifference(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             varx,               /**< variable x */
   SCIP_Real             fracx,              /**< the fractionality of variable x */
   SCIP_VAR*             vary,               /**< variable y */
   SCIP_Real             fracy,              /**< the fractionality of variable y */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_CONFIDENCELEVEL  clevel,             /**< confidence level for rejecting hypothesis */
   SCIP_Bool             onesided            /**< should a one-sided hypothesis y >= x be tested? */
   );

/** tests at a given confidence level whether the variable pseudo-costs only have a small probability to
 *  exceed a \p threshold. This is useful to determine if past observations provide enough evidence
 *  to skip an expensive strong-branching step if there is already a candidate that has been proven to yield an improvement
 *  of at least \p threshold.
 *
 *  @note use \p clevel to adjust the level of confidence. For SCIP_CONFIDENCELEVEL_MIN, the method returns TRUE if
 *        the estimated probability to exceed \p threshold is less than 25 %.
 *
 *  @see  SCIP_Confidencelevel for a list of available levels. The used probability limits refer to the one-sided levels
 *        of confidence.
 *
 *  @return TRUE if the variable pseudo-cost probabilistic model is likely to be smaller than \p threshold
 *          at the given confidence level \p clevel.
 */
EXTERN
SCIP_Bool SCIPpscostThresholdProbabilityTest(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable x */
   SCIP_Real             frac,               /**< the fractionality of variable x */
   SCIP_Real             threshold,          /**< the threshold to test against */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_CONFIDENCELEVEL  clevel              /**< confidence level for rejecting hypothesis */
   );

/** check if the current pseudo cost relative error in a direction violates the given threshold. The Relative
 *  Error is calculated at a specific confidence level
 *
 *  @return TRUE if relative error in variable pseudo costs is smaller than \p threshold
 */
EXTERN
SCIP_Bool SCIPisVarPscostRelerrorReliable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable in question */
   SCIP_Real             threshold,          /**< threshold for relative errors to be considered reliable (enough) */
   SCIP_CONFIDENCELEVEL  clevel              /**< a given confidence level */
   );

/** gets the variable's pseudo cost score value for the given LP solution value
 *
 *  @return the variable's pseudo cost score value for the given LP solution value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarPseudocostScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solval              /**< variable's LP solution value */
   );

/** gets the variable's pseudo cost score value for the given LP solution value,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's pseudo cost score value for the given LP solution value,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarPseudocostScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solval              /**< variable's LP solution value */
   );

/** returns the variable's VSIDS value
 *
 *  @return the variable's VSIDS value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarVSIDS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's VSIDS value only using conflicts of the current run
 *
 *  @return the variable's VSIDS value only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarVSIDSCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's conflict score value
 *
 *  @return the variable's conflict score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarConflictScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's conflict score value only using conflicts of the current run
 *
 *  @return the variable's conflict score value only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarConflictScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's conflict length score
 *
 *  @return the variable's conflict length score
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarConflictlengthScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's conflict length score only using conflicts of the current run
 *
 *  @return the variable's conflict length score only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarConflictlengthScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's average conflict length
 *
 *  @return the variable's average conflict length
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgConflictlength(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's average  conflict length only using conflicts of the current run
 *
 *  @return the variable's average conflict length only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgConflictlengthCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of inferences found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of inferences found after branching on the variable in given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgInferences(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of inferences found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of inferences found after branching on the variable in given direction in the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgInferencesCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's average inference score value
 *
 *  @return the variable's average inference score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgInferenceScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's average inference score value only using inferences of the current run
 *
 *  @return the variable's average inference score value only using inferences of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgInferenceScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** initializes the upwards and downwards pseudocosts, conflict scores, conflict lengths, inference scores, cutoff scores
 *  of a variable to the given values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPinitVarBranchStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which should be initialized */
   SCIP_Real             downpscost,         /**< value to which pseudocosts for downwards branching should be initialized */
   SCIP_Real             uppscost,           /**< value to which pseudocosts for upwards branching should be initialized */
   SCIP_Real             downvsids,          /**< value to which VSIDS score for downwards branching should be initialized */
   SCIP_Real             upvsids,            /**< value to which VSIDS score for upwards branching should be initialized */
   SCIP_Real             downconflen,        /**< value to which conflict length score for downwards branching should be initialized */
   SCIP_Real             upconflen,          /**< value to which conflict length score for upwards branching should be initialized */
   SCIP_Real             downinfer,          /**< value to which inference counter for downwards branching should be initialized */
   SCIP_Real             upinfer,            /**< value to which inference counter for upwards branching should be initialized */
   SCIP_Real             downcutoff,         /**< value to which cutoff counter for downwards branching should be initialized */
   SCIP_Real             upcutoff            /**< value to which cutoff counter for upwards branching should be initialized */
   );

/** initializes the upwards and downwards conflict scores, conflict lengths, inference scores, cutoff scores of a
 *  variable w.r.t. a value by the given values (SCIP_VALUEHISTORY)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPinitVarValueBranchStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which should be initialized */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_Real             downvsids,          /**< value to which VSIDS score for downwards branching should be initialized */
   SCIP_Real             upvsids,            /**< value to which VSIDS score for upwards branching should be initialized */
   SCIP_Real             downconflen,        /**< value to which conflict length score for downwards branching should be initialized */
   SCIP_Real             upconflen,          /**< value to which conflict length score for upwards branching should be initialized */
   SCIP_Real             downinfer,          /**< value to which inference counter for downwards branching should be initialized */
   SCIP_Real             upinfer,            /**< value to which inference counter for upwards branching should be initialized */
   SCIP_Real             downcutoff,         /**< value to which cutoff counter for downwards branching should be initialized */
   SCIP_Real             upcutoff            /**< value to which cutoff counter for upwards branching should be initialized */
   );

/** returns the average number of cutoffs found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of cutoffs found after branching on the variable in given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of cutoffs found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of cutoffs found after branching on the variable in given direction in the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgCutoffsCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's average cutoff score value
 *
 *  @return the variable's average cutoff score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgCutoffScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's average cutoff score value, only using cutoffs of the current run
 *
 *  @return the variable's average cutoff score value, only using cutoffs of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgCutoffScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor
 *
 *  @return the variable's average inference/cutoff score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgInferenceCutoffScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   );

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor, only using inferences and cutoffs of the current run
 *
 *  @return the variable's average inference/cutoff score value, only using inferences and cutoffs of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetVarAvgInferenceCutoffScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   );

/** outputs variable information to file stream via the message system
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
EXTERN
SCIP_RETCODE SCIPprintVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/**@} */




/*
 * conflict analysis methods
 */





/*
 * constraint methods
 */





/*
 * LP methods
 */




/*
 * LP column methods
 */





/*
 * LP row methods
 */



/*
 * NLP methods
 */



/*
 * NLP diving methods
 */



/*
 * NLP nonlinear row methods
 */






/*
 * nonlinear methods
 */





/*
 * cutting plane methods
 */





/*
 * LP diving methods
 */





/*
 * probing methods
 */


/*
 * branching methods
 */





/*
 * primal solutions
 */

/**@addtogroup PublicSolutionMethods
 *
 * @{
 */

/** creates a primal solution, initialized to zero
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcreateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcreateLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current NLP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcreateNLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current relaxation solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcreateRelaxSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcreatePseudoSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current LP or pseudo solution, depending on whether the LP was solved
 *  at the current node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcreateCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a partial primal solution, initialized to unknown values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
EXTERN
SCIP_RETCODE SCIPcreatePartialSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to unknown values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcreateUnknownSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution living in the original problem space, initialized to zero;
 *  a solution in original space allows to set original variables to values that would be invalid in the
 *  transformed problem due to preprocessing fixings or aggregations
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcreateOrigSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a copy of a primal solution; note that a copy of a linked solution is also linked and needs to be unlinked
 *  if it should stay unaffected from changes in the LP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcreateSolCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_SOL*             sourcesol           /**< primal CIP solution to copy */
   );

/** creates a copy of a solution in the original primal solution space
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPcreateSolCopyOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_SOL*             sourcesol           /**< primal CIP solution to copy */
   );

/** creates a copy of a primal solution, thereby replacing infinite fixings of variables by finite values;
 *  the copy is always defined in the original variable space;
 *  success indicates whether the objective value of the solution was changed by removing infinite values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_RETCODE SCIPcreateFiniteSolCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_SOL*             sourcesol,          /**< primal CIP solution to copy */
   SCIP_Bool*            success             /**< does the finite solution have the same objective value? */
   );

/** frees primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPfreeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol                 /**< pointer to the solution */
   );

/** links a primal solution to the current LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPlinkLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the current NLP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPlinkNLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the current relaxation solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPlinkRelaxSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPlinkPseudoSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the current LP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPlinkCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** clears a primal solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPclearSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** stores solution values of variables in solution's own array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
EXTERN
SCIP_RETCODE SCIPunlinkSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** sets value of variable in primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPsetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_Real             val                 /**< solution value of variable */
   );

/** sets values of multiple variables in primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPsetSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   int                   nvars,              /**< number of variables to set solution value for */
   SCIP_VAR**            vars,               /**< array with variables to add to solution */
   SCIP_Real*            vals                /**< array with solution values of variables */
   );

/** increases value of variable in primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPincSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_VAR*             var,                /**< variable to increase solution value for */
   SCIP_Real             incval              /**< increment for solution value of variable */
   );

/** returns value of variable in primal CIP solution, or in current LP/pseudo solution
 *
 *  @return value of variable in primal CIP solution, or in current LP/pseudo solution
 *
 *  @pre In case the solution pointer @p sol is @b NULL, that means it is asked for the LP or pseudo solution, this method
 *       can only be called if @p scip is in the solving stage \ref SCIP_STAGE_SOLVING. In any other case, this method
 *       can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Real SCIPgetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   SCIP_VAR*             var                 /**< variable to get value for */
   );

/** gets values of multiple variables in primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPgetSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables to get solution value for */
   SCIP_VAR**            vars,               /**< array with variables to get value for */
   SCIP_Real*            vals                /**< array to store solution values of variables */
   );

/** returns objective value of primal CIP solution w.r.t. original problem, or current LP/pseudo objective value
 *
 *  @return objective value of primal CIP solution w.r.t. original problem, or current LP/pseudo objective value
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Real SCIPgetSolOrigObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   );

/** returns transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value
 *
 *  @return transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Real SCIPgetSolTransObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   );

/** recomputes the objective value of an original solution, e.g., when transferring solutions
 *  from the solution pool (objective coefficients might have changed in the meantime)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 */
EXTERN
SCIP_RETCODE SCIPrecomputeSolObj(
   SCIP*                 scip,
   SCIP_SOL*             sol
   );

/** maps original space objective value into transformed objective value
 *
 *  @return transformed objective value
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPtransformObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             obj                 /**< original space objective value to transform */
   );

/** maps transformed objective value into original space
 *
 *  @return objective value into original space
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPretransformObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             obj                 /**< transformed objective value to retransform in original space */
   );

/** gets clock time, when this solution was found
 *
 *  @return clock time, when this solution was found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Real SCIPgetSolTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** gets branch and bound run number, where this solution was found
 *
 *  @return branch and bound run number, where this solution was found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
int SCIPgetSolRunnum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** gets node number of the specific branch and bound run, where this solution was found
 *
 *  @return node number of the specific branch and bound run, where this solution was found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Longint SCIPgetSolNodenum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** gets heuristic, that found this solution (or NULL if it's from the tree)
 *
 *  @return heuristic, that found this solution (or NULL if it's from the tree)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_HEUR* SCIPgetSolHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** returns whether two given solutions are exactly equal
 *
 *  @return returns whether two given solutions are exactly equal
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Bool SCIPareSolsEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol1,               /**< first primal CIP solution */
   SCIP_SOL*             sol2                /**< second primal CIP solution */
   );

/** adjusts solution values of implicit integer variables in handed solution. Solution objective value is not
 *  deteriorated by this method.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPadjustImplicitSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             uselprows           /**< should LP row information be considered for none-objective variables */
   );

/** outputs non-zero variables of solution in original problem space to the given file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre In case the solution pointer @p sol is NULL (askinking for the current LP/pseudo solution), this method can be
 *       called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @pre In case the solution pointer @p sol is @b not NULL, this method can be called if @p scip is in one of the
 *       following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_RETCODE SCIPprintSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** outputs non-zero variables of solution in transformed problem space to file stream
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
 */
EXTERN
SCIP_RETCODE SCIPprintTransSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** outputs discrete variables of solution in original problem space to the given file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_RETCODE SCIPprintMIPStart(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** returns dual solution value of a constraint */
EXTERN
SCIP_RETCODE SCIPgetDualSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which the dual solution should be returned */
   SCIP_Real*            dualsolval,         /**< pointer to store the dual solution value */
   SCIP_Bool*            boundconstraint     /**< pointer to store whether the constraint is a bound constraint (or NULL) */
   );

/** check whether the dual solution is available
 *
 * @note This is used when calling \ref SCIPprintDualSol()
 *
 * @return is dual solution available?
 *
 * @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Bool SCIPisDualSolAvailable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             printreason         /**< print warning message if dualsol is not available? */
   );

/** outputs dual solution from LP solver to file stream
 *
 *  @note This only works if no presolving has been performed, which can be checked by calling method \ref
 *        SCIPhasPerformedPresolve().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPprintDualSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );


/** outputs non-zero variables of solution representing a ray in original problem space to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_RETCODE SCIPprintRay(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution representing ray */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** gets number of feasible primal solutions stored in the solution storage in case the problem is transformed;
 *  in case the problem stage is SCIP_STAGE_PROBLEM, the number of solution in the original solution candidate
 *  storage is returned
 *
 *  @return number of feasible primal solutions stored in the solution storage in case the problem is transformed; or
 *          number of solution in the original solution candidate storage if the problem stage is SCIP_STAGE_PROBLEM
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
int SCIPgetNSols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets array of feasible primal solutions stored in the solution storage in case the problem is transformed; in case
 *  if the problem stage is in SCIP_STAGE_PROBLEM, it returns the number array of solution candidate stored
 *
 *  @return array of feasible primal solutions
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_SOL** SCIPgetSols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets best feasible primal solution found so far if the problem is transformed; in case the problem is in
 *  SCIP_STAGE_PROBLEM it returns the best solution candidate, or NULL if no solution has been found or the candidate
 *  store is empty;
 *
 *  @return best feasible primal solution so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_SOL* SCIPgetBestSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** outputs best feasible primal solution found so far to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREE
 */
EXTERN
SCIP_RETCODE SCIPprintBestSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** outputs best feasible primal solution found so far in transformed variables to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREE
 */
EXTERN
SCIP_RETCODE SCIPprintBestTransSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** try to round given solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIProundSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Bool*            success             /**< pointer to store whether rounding was successful */
   );

/** retransforms solution to original problem space
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPretransformSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** reads a given solution file, problem has to be transformed in advance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPreadSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of the input file */
   );

/** reads a given solution file and store the solution values in the given solution pointer
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPreadSolFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of the input file */
   SCIP_SOL*             sol,                /**< solution pointer */
   SCIP_Bool             xml,                /**< true, iff the given solution in written in XML */
   SCIP_Bool*            partial,            /**< pointer to store if the solution is partial */
   SCIP_Bool*            error               /**< pointer store if an error occured */
   );

/** adds feasible primal solution to solution storage by copying it
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note Do not call during propagation, use heur_trysol instead.
 */
EXTERN
SCIP_RETCODE SCIPaddSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds primal solution to solution storage, frees the solution afterwards
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note Do not call during propagation, use heur_trysol instead.
 */
EXTERN
SCIP_RETCODE SCIPaddSolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds current LP/pseudo solution to solution storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPaddCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** checks solution for feasibility; if possible, adds it to storage by copying
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
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note Do not call during propagation, use heur_trysol instead.
 */
EXTERN
SCIP_RETCODE SCIPtrySol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards
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
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note Do not call during propagation, use heur_trysol instead.
 */
EXTERN
SCIP_RETCODE SCIPtrySolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violation be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   );

/** checks current LP/pseudo solution for feasibility; if possible, adds it to storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPtryCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violation be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

/** returns all partial solutions
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_SOL** SCIPgetPartialSols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns number of partial solutions
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNPartialSols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** checks solution for feasibility without adding it to the solution store
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPcheckSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violation be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            feasible            /**< stores whether given solution is feasible */
   );

/** checks solution for feasibility in original problem without adding it to the solution store;
 *  this method is used to double check a solution in order to validate the presolving process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPcheckSolOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            feasible,           /**< stores whether given solution is feasible */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool             completely          /**< should all violations be checked? */
   );

/** update integrality violation of a solution */
EXTERN
void SCIPupdateSolIntegralityViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol             /**< absolute violation */
   );

/** update bound violation of a solution */
EXTERN
void SCIPupdateSolBoundViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation */
   SCIP_Real             relviol             /**< relative violation */
   );

/** update LP row violation of a solution */
EXTERN
void SCIPupdateSolLPRowViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation */
   SCIP_Real             relviol             /**< relative violation */
   );

/** update constraint violation of a solution */
EXTERN
void SCIPupdateSolConsViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation */
   SCIP_Real             relviol             /**< relative violation */
   );

/** update LP row and constraint violations of a solution */
EXTERN
void SCIPupdateSolLPConsViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation */
   SCIP_Real             relviol             /**< relative violation */
   );

/** allow violation updates */
EXTERN
void SCIPactivateSolViolationUpdates(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** disallow violation updates */
EXTERN
void SCIPdeactivateSolViolationUpdates(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return whether a primal ray is stored that proves unboundedness of the LP relaxation
 *
 *  @return return whether a primal ray is stored that proves unboundedness of the LP relaxation
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Bool SCIPhasPrimalRay(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets value of given variable in primal ray causing unboundedness of the LP relaxation;
 *  should only be called if such a ray is stored (check with SCIPhasPrimalRay())
 *
 *  @return value of given variable in primal ray causing unboundedness of the LP relaxation
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetPrimalRayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get value for */
   );

/** updates the primal ray thats proves unboundedness
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPupdatePrimalRay(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             primalray           /**< the new primal ray */
   );

/**@} */




/*
 * event methods
 */





/*
 * tree methods
 */

/**@addtogroup PublicTreeMethods
 *
 * @{
 */

/** gets focus node in the tree
 *
 *  if we are in probing/diving mode this method returns the node in the tree where the probing/diving mode was started.
 *
 *  @return the current node of the search tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_NODE* SCIPgetFocusNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current node in the tree
 *
 *  @return the current node of the search tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_NODE* SCIPgetCurrentNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets depth of current node, or -1 if no current node exists; in probing, the current node is the last probing node,
 *  such that the depth includes the probing path
 *
 *  @return the depth of current node, or -1 if no current node exists; in probing, the current node is the last probing node,
 *  such that the depth includes the probing path
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
 */
EXTERN
int SCIPgetDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets depth of the focus node, or -1 if no focus node exists; the focus node is the currently processed node in the
 *  branching tree, excluding the nodes of the probing path
 *
 *  @return the depth of the focus node, or -1 if no focus node exists; the focus node is the currently processed node in the
 *  branching tree, excluding the nodes of the probing path
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
 */
EXTERN
int SCIPgetFocusDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current plunging depth (successive times, a child was selected as next node)
 *
 *  @return the current plunging depth (successive times, a child was selected as next node)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetPlungeDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the root node of the tree
 *
 *  @return the root node of the search tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_NODE* SCIPgetRootNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the effective root depth, i.e., the depth of the deepest node which is part of all paths from the root node
 *  to the unprocessed nodes.
 *
 *  @return effective root depth
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetEffectiveRootDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the current node is already solved and only propagated again
 *
 *  @return TRUE is returned if \SCIP performance repropagation, otherwise FALSE.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_Bool SCIPinRepropagation(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets children of focus node along with the number of children
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          children,           /**< pointer to store children array, or NULL if not needed */
   int*                  nchildren           /**< pointer to store number of children, or NULL if not needed */
   );

/** gets number of children of focus node
 *
 *  @return number of children of the focus node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetNChildren(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets siblings of focus node along with the number of siblings
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetSiblings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          siblings,           /**< pointer to store siblings array, or NULL if not needed */
   int*                  nsiblings           /**< pointer to store number of siblings, or NULL if not needed */
   );

/** gets number of siblings of focus node
 *
 *  @return the number of siblings of focus node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetNSiblings(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets leaves of the tree along with the number of leaves
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetLeaves(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          leaves,             /**< pointer to store leaves array, or NULL if not needed */
   int*                  nleaves             /**< pointer to store number of leaves, or NULL if not needed */
   );

/** gets number of leaves in the tree
 *
 *  @return the number of leaves in the tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetNLeaves(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of nodes left in the tree (children + siblings + leaves)
 *
 *  @return the number of nodes left in the tree (children + siblings + leaves)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
int SCIPgetNNodesLeft(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best child of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @return the best child of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_NODE* SCIPgetPrioChild(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @return the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_NODE* SCIPgetPrioSibling(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best child of the focus node w.r.t. the node selection strategy
 *
 *  @return the best child of the focus node w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_NODE* SCIPgetBestChild(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best sibling of the focus node w.r.t. the node selection strategy
 *
 *  @return the best sibling of the focus node w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_NODE* SCIPgetBestSibling(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best leaf from the node queue w.r.t. the node selection strategy
 *
 *  @return the best leaf from the node queue w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_NODE* SCIPgetBestLeaf(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy
 *
 *  @return the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_NODE* SCIPgetBestNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the node with smallest lower bound from the tree (child, sibling, or leaf)
 *
 *  @return the node with smallest lower bound from the tree (child, sibling, or leaf)
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_NODE* SCIPgetBestboundNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** access to all data of open nodes (leaves, children, and siblings)
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPgetOpenNodesData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          leaves,             /**< pointer to store the leaves, or NULL if not needed */
   SCIP_NODE***          children,           /**< pointer to store the children, or NULL if not needed */
   SCIP_NODE***          siblings,           /**< pointer to store the siblings, or NULL if not needed */
   int*                  nleaves,            /**< pointer to store the number of leaves, or NULL */
   int*                  nchildren,          /**< pointer to store the number of children, or NULL */
   int*                  nsiblings           /**< pointer to store the number of siblings, or NULL */
   );

/** cuts off node and whole sub tree from branch and bound tree
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPcutoffNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node that should be cut off */
   );

/** marks the given node to be propagated again the next time a node of its subtree is processed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPrepropagateNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node that should be propagated again */
   );

/** returns depth of first node in active path that is marked being cutoff
 *
 *  @return depth of first node in active path that is marked being cutoff
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetCutoffdepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns depth of first node in active path that has to be propagated again
 *
 *  @return depth of first node in active path that has to be propagated again
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetRepropdepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** prints all branching decisions on variables from the root to the given node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPprintNodeRootPath(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** sets whether the LP should be solved at the focus node
 *
 *  @note In order to have an effect, this method needs to be called after a node is focused but before the LP is
 *        solved.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
void SCIPsetFocusnodeLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             solvelp             /**< should the LP be solved? */
   );

/**@} */


/*
 * parallel interface methods
 */




/*
 * reoptimization methods
 */

/**@addtogroup PublicReoptimizationMethods
 *
 * @{
 */

/** return the ids of child nodes stored in the reoptimization tree
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPgetReoptChildIDs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int*         ids,                /**< array to store the ids of child nodes */
   int                   mem,                /**< allocated memory */
   int*                  nids                /**< number of child nodes */
   );

/** return the ids of all leave nodes store in the reoptimization tree induced by the given node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPgetReoptLeaveIDs(
   SCIP*                 scip,               /**< SCIP data strcuture */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int*         ids,                /**< array of ids */
   int                   mem,                /**< allocated memory */
   int*                  nids                /**< number of child nodes */
   );

/** returns the number of nodes in the reoptimization tree induced by @p node; if @p node == NULL, the method
 *  returns the number of nodes of the whole reoptimization tree.
 */
EXTERN
int SCIPgetNReoptnodes(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** returns the number of leave nodes of the subtree induced by @p node; if @p node == NULL, the method
 *  returns the number of leaf nodes of the whole reoptimization tree.
 */
EXTERN
int SCIPgetNReoptLeaves(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** gets the node of the reoptimization tree corresponding to the unique @p id */
SCIP_REOPTNODE* SCIPgetReoptnode(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          id                  /**< unique id */
   );

/** add a variable bound change to a given reoptnode
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPaddReoptnodeBndchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR*             var,                /**< variable pointer */
   SCIP_Real             bound,              /**< variable bound to add */
   SCIP_BOUNDTYPE        boundtype           /**< bound type of the variable value */
   );

/** set the @p representation as the new search frontier
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 */
EXTERN
SCIP_RETCODE SCIPsetReoptCompression(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE**      representation,     /**< array of representatives */
   int                   nrepresentatives,   /**< number of representatives */
   SCIP_Bool*            success             /**< pointer to store the result */
   );

/** add stored constraint to a reoptimization node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 */
EXTERN
SCIP_RETCODE SCIPaddReoptnodeCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of variable bounds */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of variable boundtypes */
   SCIP_Real             lhs,                /**< lhs of the constraint */
   SCIP_Real             rhs,                /**< rhs of the constraint */
   int                   nvars,              /**< number of variables */
   REOPT_CONSTYPE        constype,           /**< type of the constraint */
   SCIP_Bool             linear              /**< the given constraint has a linear representation */
   );

/** return the branching path stored in the reoptree at ID id */
EXTERN
void SCIPgetReoptnodePath(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of variable bounds */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of bound types */
   int                   mem,                /**< allocated memory */
   int*                  nvars,              /**< number of variables */
   int*                  nafterdualvars      /**< number of variables directly after the first based on dual information */
   );

/** initialize a set of empty reoptimization nodes
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 */
EXTERN
SCIP_RETCODE SCIPinitRepresentation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE**      representatives,    /**< array of representatives */
   int                   nrepresentatives    /**< number of representatives */
   );

/** reset a set of initialized reoptimization nodes
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 */
SCIP_RETCODE SCIPresetRepresentation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE**      representatives,    /**< array of representatives */
   int                   nrepresentatives    /**< number of representatives */
   );

/** free a set of initialized reoptimization nodes
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 */
EXTERN
SCIP_RETCODE SCIPfreeRepresentation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE**      representatives,    /**< array of representatives */
   int                   nrepresentatives    /**< number of representatives */
   );

/** reactivate the given @p reoptnode and split them into several nodes if necessary
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPapplyReopt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node to reactivate */
   unsigned int          id,                 /**< unique id of the reoptimization node */
   SCIP_Real             estimate,           /**< estimate of the child nodes that should be created */
   SCIP_NODE**           childnodes,         /**< array to store the created child nodes */
   int*                  ncreatedchilds,     /**< pointer to store number of created child nodes */
   int*                  naddedconss,        /**< pointer to store number of generated constraints */
   int                   childnodessize,     /**< available size of childnodes array */
   SCIP_Bool*            success             /**< pointer store the result*/
   );

/** remove the stored information about bound changes based in dual information
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPresetReoptnodeDualcons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** splits the root into several nodes and moves the child nodes of the root to one of the created nodes
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPsplitReoptRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  ncreatedchilds,     /**< pointer to store the number of created nodes */
   int*                  naddedconss         /**< pointer to store the number added constraints */
   );

/** returns if a node should be reoptimized */
EXTERN
SCIP_Bool SCIPreoptimizeNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** deletes the given reoptimization node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPdeleteReoptnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE**      reoptnode           /**< node of the reoptimization tree */
   );

/** return the similarity between two objective functions */
EXTERN
SCIP_Real SCIPgetReoptSimilarity(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   run1,               /**< number of run */
   int                   run2                /**< number of run */
   );

/** check the changes of the variable coefficient in the objective function */
EXTERN
void SCIPgetVarCoefChg(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   varidx,             /**< index of variable */
   SCIP_Bool*            negated,            /**< coefficient changed the sign */
   SCIP_Bool*            entering,           /**< coefficient gets non-zero coefficient */
   SCIP_Bool*            leaving             /**< coefficient gets zero coefficient */
   );

/* @} */

/*
 * statistic methods
 */

/**@addtogroup PublicSolvingStatsMethods
 *
 * @{
 */

/** gets number of branch and bound runs performed, including the current run
 *
 *  @return the number of branch and bound runs performed, including the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
int SCIPgetNRuns(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of reoptimization runs performed, including the current run
 *
 *  @return the number of reoptimization runs performed, including the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
int SCIPgetNReoptRuns(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** add given number to the number of processed nodes in current run and in all runs, including the focus node
 *
 *  @return the number of processed nodes in current run, including the focus node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
void SCIPaddNNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint          nnodes              /**< number of processed nodes to add to the statistics */
   );

/** gets number of processed nodes in current run, including the focus node
 *
 *  @return the number of processed nodes in current run, including the focus node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Longint SCIPgetNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of processed nodes in all runs, including the focus node
 *
 *  @return the total number of processed nodes in all runs, including the focus node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Longint SCIPgetNTotalNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of leaf nodes processed with feasible relaxation solution
 *
 * @return number of leaf nodes processed with feasible relaxation solution
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Longint SCIPgetNFeasibleLeaves(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of infeasible leaf nodes processed
 *
 * @return number of infeasible leaf nodes processed
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Longint SCIPgetNInfeasibleLeaves(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of processed leaf nodes that hit LP objective limit
 *
 * @return number of processed leaf nodes that hit LP objective limit
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Longint SCIPgetNObjlimLeaves(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of times a selected node was from a cut off subtree
 *
 *  @return number of times a selected node was from a cut off subtree
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Longint SCIPgetNDelayedCutoffs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far
 *
 *  @return the total number of LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_Longint SCIPgetNLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of iterations used so far in primal and dual simplex and barrier algorithm
 *
 *  @return the total number of iterations used so far in primal and dual simplex and barrier algorithm
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of active non-zeros in the current transformed problem
 *
 *  @return the number of active non-zeros in the current transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_Longint SCIPgetNNZs(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** gets total number of iterations used so far in primal and dual simplex and barrier algorithm for the root node
 *
 *  @return the total number of iterations used so far in primal and dual simplex and barrier algorithm for the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNRootLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of iterations used in primal and dual simplex and barrier algorithm for the first LP at the root
 *  node
 *
 *  @return the total number of iterations used in primal and dual simplex and barrier algorithm for the first root LP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNRootFirstLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of primal LPs solved so far
 *
 *  @return the total number of primal LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNPrimalLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of iterations used so far in primal simplex
 *
 *  @return total number of iterations used so far in primal simplex
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNPrimalLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of dual LPs solved so far
 *
 *  @return the total number of dual LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNDualLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of iterations used so far in dual simplex
 *
 *  @return the total number of iterations used so far in dual simplex
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNDualLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of barrier LPs solved so far
 *
 *  @return the total number of barrier LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNBarrierLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of iterations used so far in barrier algorithm
 *
 *  @return the total number of iterations used so far in barrier algorithm
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNBarrierLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far that were resolved from an advanced start basis
 *
 *  @return the total number of LPs solved so far that were resolved from an advanced start basis
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in primal and dual simplex calls where an advanced start basis
 *  was available
 *
 *  @return the total number of simplex iterations used so far in primal and dual simplex calls where an advanced start
 *          basis was available
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of primal LPs solved so far that were resolved from an advanced start basis
 *
 *  @return the total number of primal LPs solved so far that were resolved from an advanced start basis
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNPrimalResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in primal simplex calls where an advanced start basis
 *  was available
 *
 *  @return the total number of simplex iterations used so far in primal simplex calls where an advanced start
 *          basis was available
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNPrimalResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of dual LPs solved so far that were resolved from an advanced start basis
 *
 *  @return the total number of dual LPs solved so far that were resolved from an advanced start basis
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNDualResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in dual simplex calls where an advanced start basis
 *  was available
 *
 *  @return the total number of simplex iterations used so far in dual simplex calls where an advanced start
 *          basis was available
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNDualResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far for node relaxations
 *
 *  @return the total number of LPs solved so far for node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNNodeLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far for node relaxations
 *
 *  @return the total number of simplex iterations used so far for node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNNodeLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far for initial LP in node relaxations
 *
 *  @return the total number of LPs solved so far for initial LP in node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNNodeInitLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far for initial LP in node relaxations
 *
 *  @return the total number of simplex iterations used so far for initial LP in node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNNodeInitLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far during diving and probing
 *
 *  @return total number of LPs solved so far during diving and probing
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNDivingLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far during diving and probing
 *
 *  @return the total number of simplex iterations used so far during diving and probing
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNDivingLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of times, strong branching was called (each call represents solving two LPs)
 *
 *  @return the total number of times, strong branching was called (each call represents solving two LPs)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNStrongbranchs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in strong branching
 *
 *  @return the total number of simplex iterations used so far in strong branching
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNStrongbranchLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of times, strong branching was called at the root node (each call represents solving two LPs)
 *
 *  @return the total number of times, strong branching was called at the root node (each call represents solving two LPs)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNRootStrongbranchs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in strong branching at the root node
 *
 *  @return the total number of simplex iterations used so far in strong branching at the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Longint SCIPgetNRootStrongbranchLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of pricing rounds performed so far at the current node
 *
 *  @return the number of pricing rounds performed so far at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetNPriceRounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get current number of variables in the pricing store
 *
 *  @return the current number of variables in the pricing store
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
int SCIPgetNPricevars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of pricing variables found so far
 *
 *  @return the total number of pricing variables found so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
int SCIPgetNPricevarsFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of pricing variables applied to the LPs
 *
 *  @return the total number of pricing variables applied to the LPs
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
int SCIPgetNPricevarsApplied(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of separation rounds performed so far at the current node
 *
 *  @return the number of separation rounds performed so far at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetNSepaRounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of cuts found so far
 *
 *  @return the total number of cuts found so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
int SCIPgetNCutsFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get number of cuts found so far in current separation round
 *
 *  @return the number of cuts found so far in current separation round
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
int SCIPgetNCutsFoundRound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of cuts applied to the LPs
 *
 *  @return the total number of cuts applied to the LPs
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
int SCIPgetNCutsApplied(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of constraints found in conflict analysis (conflict and reconvergence constraints)
 *
 *  @return the total number of constraints found in conflict analysis (conflict and reconvergence constraints)
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
 */
EXTERN
SCIP_Longint SCIPgetNConflictConssFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get number of conflict constraints found so far at the current node
 *
 *  @return the number of conflict constraints found so far at the current node
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
 */
EXTERN
int SCIPgetNConflictConssFoundNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of conflict constraints added to the problem
 *
 *  @return the total number of conflict constraints added to the problem
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
 */
EXTERN
SCIP_Longint SCIPgetNConflictConssApplied(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets maximal depth of all processed nodes in current branch and bound run (excluding probing nodes)
 *
 *  @return the maximal depth of all processed nodes in current branch and bound run (excluding probing nodes)
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
 */
EXTERN
int SCIPgetMaxDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets maximal depth of all processed nodes over all branch and bound runs
 *
 *  @return the maximal depth of all processed nodes over all branch and bound runs
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
 */
EXTERN
int SCIPgetMaxTotalDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of backtracks, i.e. number of times, the new node was selected from the leaves queue
 *
 *  @return the total number of backtracks, i.e. number of times, the new node was selected from the leaves queue
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
 */
EXTERN
SCIP_Longint SCIPgetNBacktracks(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of active constraints at the current node
 *
 *  @return the total number of active constraints at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetNActiveConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of enabled constraints at the current node
 *
 *  @return the total number of enabled constraints at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
int SCIPgetNEnabledConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets average dual bound of all unprocessed nodes for original problem */
EXTERN
SCIP_Real SCIPgetAvgDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets average lower (dual) bound of all unprocessed nodes in transformed problem
 *
 *  @return the average lower (dual) bound of all unprocessed nodes in transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets global dual bound
 *
 *  @return the global dual bound
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets global lower (dual) bound in transformed problem
 *
 *  @return the global lower (dual) bound in transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets dual bound of the root node for the original problem
 *
 *  @return the dual bound of the root node for the original problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetDualboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets lower (dual) bound in transformed problem of the root node
 *
 *  @return the lower (dual) bound in transformed problem of the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetLowerboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets dual bound for the original problem obtained by the first LP solve at the root node
 *
 *  @return the dual bound for the original problem of the first LP solve at the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetFirstLPDualboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets lower (dual) bound in transformed problem obtained by the first LP solve at the root node
 *
 *  @return the lower (dual) bound in transformed problem obtained by first LP solve at the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetFirstLPLowerboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** the primal bound of the very first solution */
EXTERN
SCIP_Real SCIPgetFirstPrimalBound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets global primal bound (objective value of best solution or user objective limit) for the original problem
 *
 *  @return the global primal bound (objective value of best solution or user objective limit) for the original problem
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
 */
EXTERN
SCIP_Real SCIPgetPrimalbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets global upper (primal) bound in transformed problem (objective value of best solution or user objective limit)
 *
 *  @return the global upper (primal) bound in transformed problem (objective value of best solution or user objective limit)
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
 */
EXTERN
SCIP_Real SCIPgetUpperbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets global cutoff bound in transformed problem: a sub problem with lower bound larger than the cutoff
 *  cannot contain a better feasible solution; usually, this bound is equal to the upper bound, but if the
 *  objective value is always integral, the cutoff bound is (nearly) one less than the upper bound;
 *  additionally, due to objective function domain propagation, the cutoff bound can be further reduced
 *
 *  @return global cutoff bound in transformed problem
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
 */
EXTERN
SCIP_Real SCIPgetCutoffbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** updates the cutoff bound
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @note using this method in the solving stage can lead to an erroneous SCIP solving status; in particular,
 *        if a solution not respecting the cutoff bound was found before installing a cutoff bound which
 *        renders the remaining problem infeasible, this solution may be reported as optimal
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note the given cutoff bound has to better or equal to known one (SCIPgetCutoffbound())
 */
EXTERN
SCIP_RETCODE SCIPupdateCutoffbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             cutoffbound         /**< new cutoff bound */
   );

/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 *
 *  @return TRUE if the current primal bound is justified with a feasible primal solution, otherwise FALSE
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
 */
EXTERN
SCIP_Bool SCIPisPrimalboundSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current gap |(primalbound - dualbound)/min(|primalbound|,|dualbound|)| if both bounds have same sign,
 *  or infinity, if they have opposite sign
 *
 *  @return the current gap |(primalbound - dualbound)/min(|primalbound|,|dualbound|)| if both bounds have same sign,
 *  or infinity, if they have opposite sign
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetGap(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current gap |(upperbound - lowerbound)/min(|upperbound|,|lowerbound|)| in transformed problem if both bounds
 *  have same sign, or infinity, if they have opposite sign
 *
 *  @return current gap |(upperbound - lowerbound)/min(|upperbound|,|lowerbound|)| in transformed problem if both bounds
 *  have same sign, or infinity, if they have opposite sign
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetTransGap(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of feasible primal solutions found so far
 *
 *  @return the number of feasible primal solutions found so far
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
 */
EXTERN
SCIP_Longint SCIPgetNSolsFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of feasible primal solutions respecting the objective limit found so far
 *
 *  @return the number of feasible primal solutions respecting the objective limit found so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
EXTERN
SCIP_Longint SCIPgetNLimSolsFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of feasible primal solutions found so far, that improved the primal bound at the time they were found
 *
 *  @return the number of feasible primal solutions found so far, that improved the primal bound at the time they were found
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
 */
EXTERN
SCIP_Longint SCIPgetNBestSolsFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average pseudo cost value for the given direction over all variables
 *
 *  @return the average pseudo cost value for the given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the average pseudo cost value for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average pseudo cost value for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgPseudocostCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the average number of pseudo cost updates for the given direction over all variables
 *
 *  @return the average number of pseudo cost updates for the given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgPseudocostCount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the average number of pseudo cost updates for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average number of pseudo cost updates for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgPseudocostCountCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the number of pseudo cost updates for the given direction over all variables
 *
 *  @return the number of pseudo cost updates for the given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetPseudocostCount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Bool             onlycurrentrun      /**< use only history of current run? */
   );

/** gets the average pseudo cost score value over all variables, assuming a fractionality of 0.5
 *
 *  @return the average pseudo cost score value over all variables, assuming a fractionality of 0.5
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgPseudocostScore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the variance of pseudo costs for all variables in the requested direction
 *
 *  @return the variance of pseudo costs for all variables in the requested direction
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetPseudocostVariance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        branchdir,          /**< the branching direction, up or down */
   SCIP_Bool             onlycurrentrun      /**< use only history of current run? */
   );

/** gets the average pseudo cost score value over all variables, assuming a fractionality of 0.5,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average pseudo cost score value over all variables, assuming a fractionality of 0.5,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgPseudocostScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average conflict score value over all variables */
EXTERN
SCIP_Real SCIPgetAvgConflictScore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average conflict score value over all variables, only using the conflict information of the current run
 *
 *  @return the average conflict score value over all variables, only using the conflict information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgConflictScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average inference score value over all variables
 *
 *  @return the average inference score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgConflictlengthScore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average conflictlength score value over all variables, only using the conflictlength information of the
 *  current run
 *
 *  @return the average conflictlength score value over all variables, only using the conflictlength information of the
 *          current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgConflictlengthScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the average number of inferences found after branching in given direction over all variables
 *
 *  @return the average number of inferences found after branching in given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgInferences(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of inferences found after branching in given direction over all variables,
 *  only using the inference information of the current run
 *
 *  @return the average number of inferences found after branching in given direction over all variables,
 *          only using the inference information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgInferencesCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the average inference score value over all variables
 *
 *  @return the average inference score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgInferenceScore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average inference score value over all variables, only using the inference information of the
 *  current run
 *
 *  @return the average inference score value over all variables, only using the inference information of the
 *          current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgInferenceScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the average number of cutoffs found after branching in given direction over all variables
 *
 *  @return the average number of cutoffs found after branching in given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of cutoffs found after branching in given direction over all variables,
 *  only using the cutoff information of the current run
 *
 *  @return the average number of cutoffs found after branching in given direction over all variables,
 *          only using the cutoff information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgCutoffsCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the average cutoff score value over all variables
 *
 *  @return the average cutoff score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgCutoffScore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average cutoff score value over all variables, only using the cutoff information of the current run
 *
 *  @return the average cutoff score value over all variables, only using the cutoff information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetAvgCutoffScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets deterministic time number of LPs solved so far
 *
 *  @return the total number of LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_Real SCIPgetDeterministicTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** outputs original problem to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
EXTERN
SCIP_RETCODE SCIPprintOrigProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format (or NULL for default CIP format)*/
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   );

/** outputs transformed problem of the current node to file stream
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
EXTERN
SCIP_RETCODE SCIPprintTransProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format (or NULL for default CIP format)*/
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   );

/** outputs status statistics
 *
 *  @note If limits have been changed between the solution and the call to this function, the status is recomputed and
 *        thus may to correspond to the original status.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintStatusStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs timing statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintTimingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs statistics for original problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintOrigProblemStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** outputs statistics for transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintTransProblemStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs presolver statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintPresolverStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs constraint statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintConstraintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs constraint timing statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintConstraintTimingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs propagator statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintPropagatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs conflict statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintConflictStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs separator statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintSeparatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs pricer statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintPricerStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs branching rule statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintBranchruleStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs heuristics statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintHeuristicStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs compression statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintCompressionStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs LP statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs NLP statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintNLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs relaxator statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintRelaxatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs tree statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintTreeStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs root statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintRootStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs solution statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintSolutionStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs concurrent solver statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintConcsolverStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs Benders' decomposition statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
void SCIPprintBendersStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   );

/** outputs solving statistics
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @note If limits have been changed between the solution and the call to this function, the status is recomputed and
 *        thus may to correspond to the original status.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPprintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** outputs reoptimization statistics
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPprintReoptStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** outputs history statistics about branchings on variables
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
EXTERN
SCIP_RETCODE SCIPprintBranchingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** outputs node information display line
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
EXTERN
SCIP_RETCODE SCIPprintDisplayLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VERBLEVEL        verblevel,          /**< minimal verbosity level to actually display the information line */
   SCIP_Bool             endline             /**< should the line be terminated with a newline symbol? */
   );

/** gets total number of implications between variables that are stored in the implication graph
 *
 *  @return the total number of implications between variables that are stored in the implication graph
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
EXTERN
int SCIPgetNImplications(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** stores conflict graph of binary variables' implications into a file, which can be used as input for the DOT tool
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
 *
 *  @deprecated because binary implications are now stored as cliques
 */
EXTERN
SCIP_DEPRECATED
SCIP_RETCODE SCIPwriteImplicationConflictGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name, or NULL for stdout */
   );


/** update statistical information when a new solution was found */
EXTERN
void SCIPstoreSolutionGap(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */

/*
 * timing methods
 */


/**@addtogroup PublicTimingMethods
 *
 * @{
 */

/** gets current time of day in seconds (standard time zone)
 *
 *  @return the current time of day in seconds (standard time zone).
 */
EXTERN
SCIP_Real SCIPgetTimeOfDay(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a clock using the default clock type
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPcreateClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** creates a clock counting the CPU user seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPcreateCPUClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** creates a clock counting the wall clock seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPcreateWallClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** frees a clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPfreeClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** resets the time measurement of a clock to zero and completely stops the clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPresetClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** starts the time measurement of a clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPstartClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** stops the time measurement of a clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPstopClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** enables or disables all statistic clocks of SCIP concerning plugin statistics,
 *  LP execution time, strong branching time, etc.
 *
 *  Method reads the value of the parameter timing/statistictiming. In order to disable statistic timing,
 *  set the parameter to FALSE.
 *
 *  @note: The (pre-)solving time clocks which are relevant for the output during (pre-)solving
 *         are not affected by this method
 *
 *  @see: For completely disabling all timing of SCIP, consider setting the parameter timing/enabled to FALSE
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPenableOrDisableStatisticTiming(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** starts the current solving time
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPstartSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** stops the current solving time in seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPstopSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the measured time of a clock in seconds
 *
 *  @return the measured time of a clock in seconds.
 */
EXTERN
SCIP_Real SCIPgetClockTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** sets the measured time of a clock to the given value in seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPsetClockTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_Real             sec                 /**< time in seconds to set the clock's timer to */
   );

/** gets the current total SCIP time in seconds, possibly accumulated over several problems.
 *
 *  @return the current total SCIP time in seconds, ie. the total time since the SCIP instance has been created
 */
EXTERN
SCIP_Real SCIPgetTotalTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the current solving time in seconds
 *
 *  @return the current solving time in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_Real SCIPgetSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the current reading time in seconds
 *
 *  @return the current reading time in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_Real SCIPgetReadingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the current presolving time in seconds
 *
 *  @return the current presolving time in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_Real SCIPgetPresolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the time need to solve the first LP in the root node
 *
 *  @return the solving time for the first LP in the root node in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_Real SCIPgetFirstLPTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */

/*
 * numeric values and comparisons
 */





/*
 * memory management
 */




/*
 * dynamic arrays
 */




/**@addtogroup PublicValidationMethods
 *
 * @{
 */

/** validate the result of the solve
 *
 *  the validation includes
 *
 *  - checking the feasibility of the incumbent solution in the original problem (using SCIPcheckSolOrig())
 *
 *  - checking if the objective bounds computed by SCIP agree with external primal and dual reference bounds.
 *
 *  All external reference bounds the original problem space and the original objective sense.
 *
 *  For infeasible problems, +/-SCIPinfinity() should be passed as reference bounds depending on the objective sense
 *  of the original problem.
 */
EXTERN
SCIP_RETCODE SCIPvalidateSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             primalreference,    /**< external primal reference value for the problem, or SCIP_UNKNOWN */
   SCIP_Real             dualreference,      /**< external dual reference value for the problem, or SCIP_UNKNOWN */
   SCIP_Real             reftol,             /**< relative tolerance for acceptable violation of reference values */
   SCIP_Bool             quiet,              /**< TRUE if no status line should be printed */
   SCIP_Bool*            feasible,           /**< pointer to store if the best solution is feasible in the original problem,
                                               *  or NULL */
   SCIP_Bool*            primalboundcheck,   /**< pointer to store if the primal bound respects the given dual reference
                                               *  value, or NULL */
   SCIP_Bool*            dualboundcheck      /**< pointer to store if the dual bound respects the given primal reference
                                               *  value, or NULL */
   );

/**@} */


#ifdef __cplusplus
}
#endif

#endif
