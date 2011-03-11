/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip.h
 * @ingroup PUBLICMETHODS
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

#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_cons.h"
#include "scip/type_dialog.h"
#include "scip/type_disp.h"
#include "scip/type_heur.h"
#include "scip/type_nodesel.h"
#include "scip/type_presol.h"
#include "scip/type_pricer.h"
#include "scip/type_reader.h"
#include "scip/type_relax.h"
#include "scip/type_sepa.h"
#include "scip/type_prop.h"
#include "nlpi/type_nlpi.h"

/* include public interfaces, s.t. the user only needs to include scip.h */
#include "scip/pub_branch.h"
#include "scip/pub_conflict.h"
#include "scip/pub_cons.h"
#include "scip/pub_cutpool.h"
#include "scip/pub_dialog.h"
#include "scip/pub_disp.h"
#include "scip/pub_event.h"
#include "scip/pub_fileio.h"
#include "scip/pub_heur.h"
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
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/lpi.h"
#include "nlpi/pub_expr.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/var.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * miscellaneous methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** returns SCIP version number */
extern
SCIP_Real SCIPversion(
   void
   );

/** returns SCIP major version */
extern
int SCIPmajorVersion(
   void
   );

/** returns SCIP minor version */
extern
int SCIPminorVersion(
   void
   );

/** returns SCIP technical version */
extern
int SCIPtechVersion(
   void
   );

/** returns SCIP sub version number */
extern
int SCIPsubversion(
   void
   );

/** prints a version information line to a file stream via the message handler system
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
extern
void SCIPprintVersion(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** prints error message for the given SCIP return code via the error prints method */
extern
void SCIPprintError(
   SCIP_RETCODE          retcode             /**< SCIP return code causing the error */
   );

/**@} */




/*
 * general SCIP methods
 */

/**@name General SCIP Methods */
/**@{ */

/** creates and initializes SCIP data structures
 *
 *  @note The SCIP default message handler is installed. Use the method SCIPsetMessagehdlr() or SCIPsetMessagehdlrFree()
 *        to installed your own message or SCIPsetMessagehdlrLogfile() and SCIPsetMessagehdlrQuiet() to write into a log
 *        file and turn off/on the display output.
 */
extern
SCIP_RETCODE SCIPcreate(
   SCIP**                scip                /**< pointer to SCIP data structure */
   );

/** frees SCIP data structures */
extern
SCIP_RETCODE SCIPfree(
   SCIP**                scip                /**< pointer to SCIP data structure */
   );

/** returns current stage of SCIP */
extern
SCIP_STAGE SCIPgetStage(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** outputs SCIP stage and solution status if applicable via the message handler
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
extern
SCIP_RETCODE SCIPprintStage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** gets solution status */
extern
SCIP_STATUS SCIPgetStatus(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** outputs solution status */
extern
SCIP_RETCODE SCIPprintStatus(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** returns whether the current stage belongs to the transformed problem space */
extern
SCIP_Bool SCIPisTransformed(
   SCIP*                 scip                /**< SCIP data structure */
   );
/** returns whether the solution process should be probably correct */
extern
SCIP_Bool SCIPisExactSolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the presolving process would be finished given no more presolving reductions are found in this
 *  presolving round
 *
 *  Checks whether the number of presolving rounds is not exceeded and the presolving reductions found in the current
 *  presolving round suffice to trigger another presolving round.
 *
 *  @note if subsequent presolvers find more reductions, presolving might continue even if the method returns FALSE
 *  @note does not check whether infeasibility or unboundedness was already detected in presolving (which would result
 *        in presolving being stopped although the method returns TRUE)
 */
extern
SCIP_Bool SCIPisPresolveFinished(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the user pressed CTRL-C to interrupt the solving process */
extern
SCIP_Bool SCIPpressedCtrlC(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the solving process should be / was stopped before proving optimality;
 *  if the solving process should be / was stopped, the status returned by SCIPgetStatus() yields
 *  the reason for the premature abort
 */
extern
SCIP_Bool SCIPisStopped(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**@} */



/*
 * message output methods
 */

/**@name Message Output Methods */
/**@{ */

/** Installs the given message handler, such that all messages are passed to this handler. A messages handler can be
 *  created via SCIPmessagehdlrCreate().
 *
 *  @note The currently installed messages handler gets not freed. That has to be done by the user using
 *        SCIPmessagehdlrFree() or use SCIPsetMessagehdlrFree().
 */
extern
SCIP_RETCODE SCIPsetMessagehdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   );

/** Installs the given message handler, such that all messages are passed to this handler. A messages handler can be
 *  created via SCIPmessagehdlrCreate(). The currently installed messages handler gets freed.
 */
extern
SCIP_RETCODE SCIPsetMessagehdlrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   );

/** returns the currently installed message handler, or NULL if messages are currently suppressed */
extern
SCIP_MESSAGEHDLR* SCIPgetMessagehdlr(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the log file name for the currently installed message handler */
extern
void SCIPsetMessagehdlrLogfile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of log file, or NULL (stdout) */
   );

/** sets the currently installed message handler to be quiet (or not) */
extern
void SCIPsetMessagehdlrQuiet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   );

/** prints a warning message via the message handler */
extern
void SCIPwarningMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a dialog message that requests user interaction or is a direct response to a user interactive command */
extern
void SCIPdialogMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message */
extern
void SCIPinfoMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message depending on the verbosity level */
extern
void SCIPverbMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** returns the current message verbosity level */
extern
SCIP_VERBLEVEL SCIPgetVerbLevel(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**@} */




/*
 * SCIP copy methods
 */

/**@name Copy Methods */
/**@{ */

/** copies plugins from sourcescip to targetscip; in case that a constraint handler which does not need constraints
 *  cannot be copied, valid will return FALSE. All plugins can declare that, if their copy process failed, the 
 *  copied SCIP instance might not represent the same problem semantics as the original. 
 *  Note that in this case dual reductions might be invalid.
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
extern
SCIP_RETCODE SCIPcopyPlugins(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_Bool             copyreaders,        /**< should the file readers be copied */
   SCIP_Bool             copypricers,        /**< should the variable pricers be copied */
   SCIP_Bool             copyconshdlrs,      /**< should the constraint handlers be copied */
   SCIP_Bool             copyconflicthdlrs,  /**< should the conflict handlers be copied */
   SCIP_Bool             copypresolvers,     /**< should the presolvers be copied */
   SCIP_Bool             copyrelaxators,     /**< should the relaxation handler be copied */
   SCIP_Bool             copyseparators,     /**< should the separators be copied */
   SCIP_Bool             copypropagators,    /**< should the propagators be copied */
   SCIP_Bool             copyheuristics,     /**< should the heuristics be copied */
   SCIP_Bool             copyeventhdlrs,     /**< should the event handlers be copied */
   SCIP_Bool             copynodeselectors,  /**< should the node selectors be copied */
   SCIP_Bool             copybranchrules,    /**< should the branchrules be copied */
   SCIP_Bool             copydisplays,       /**< should the display columns be copied */
   SCIP_Bool             copydialogs,        /**< should the dialogs be copied */
   SCIP_Bool             copynlpis,          /**< should the NLPIs be copied */
   SCIP_Bool*            valid               /**< pointer to store whether plugins, in particular all constraint
                                              *   handlers which do not need constraints were validly copied */
   );

/** create a problem by copying the problem data of the source SCIP
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
extern
SCIP_RETCODE SCIPcopyProb(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   const char*           name                /**< problem name */
   );

/** returns copy of the source variable; if there already is a copy of the source variable in the variable hash map,
 *  it is just returned as target variable; elsewise a new variable will be created and added to the target SCIP; this
 *  created variable is added to the variable hash map and returned as target variable
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *  @note if a new variable was created, this variable will be added to the target-SCIP, but it is not captured
 */
extern
SCIP_RETCODE SCIPgetVarCopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_VAR*             sourcevar,          /**< source variable */
   SCIP_VAR**            targetvar,          /**< pointer to store the target variable */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< should global or local bounds be used? */
   SCIP_Bool*            success             /**< pointer to store whether the copying was successful or not */
   );

/** copies all active variables from source-SCIP and adds these variable to the target-SCIP; the mapping between these
 *  variables are stored in the variable hashmap, target-SCIP has to be in problem creation stage, fixed and aggregated
 *  variables do not get copied
 *
 *  @note the variables are added to the target-SCIP but not captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
extern
SCIP_RETCODE SCIPcopyVars(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables  to the corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global              /**< should global or local bounds be used? */
   );

/** returns copy of the source constraint; if there already is a copy of the source constraint in the constraint hash
 *  map, it is just returned as target constraint; elsewise a new constraint will be created and captured in the target
 *  SCIP; this created constraint is added to the constraint hash map and returned as target constraint; the variable
 *  map is used to map the variables of the source SCIP to the variables of the target SCIP
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution may
 *           be declared feasible even if it violates this particular constraint.  This constellation should only be
 *           used, if no LP or pseudo solution can violate the constraint -- e.g. if a local constraint is redundant due
 *           to the variable's local bounds.
 *
 *  @note The constraint is not added to the target SCIP
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
extern
SCIP_RETCODE SCIPgetConsCopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_CONS*            sourcecons,         /**< source constraint of the source SCIP */
   SCIP_CONS**           targetcons,         /**< pointer to store the created target constraint */
   SCIP_CONSHDLR*        sourceconshdlr,     /**< source constraint handler for this constraint */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to the corresponding
                                              *   variables of the target SCIP, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   const char*           name,               /**< name of constraint, or NULL if the name of the source constraint should be used */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            success             /**< pointer to store whether the copying was successful or not */
   );

/** copies constraints from the source-SCIP and adds these to the target-SCIP; for mapping the
 *  variables between the source and the target SCIP a hash map can be given; if the variable hash
 *  map is NULL or necessary variable mapping is missing, the required variables are created in the
 *  target-SCIP and added to the hash map, if not NULL; all variables which are created are added to
 *  the target-SCIP but not (user) captured; if the constraint hash map is not NULL the mapping
 *  between the constraints of the source and target-SCIP is stored
 *
 *  @note the constraints are added to the target-SCIP but are not (user) captured in the target SCIP
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
extern
SCIP_RETCODE SCIPcopyConss(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to the corresponding
                                              *   variables of the target SCIP, must not be NULL! */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance? 
                                              *   If TRUE, the modifiable flag of constraints will be copied. */
   SCIP_Bool*            valid               /**< pointer to store whether all constraints were validly copied */
   );

/** convert all active cuts from cutpool of sourcescip to linear constraints in targetscip, sourcescip and targetscip
 *  could be the same
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
extern
SCIP_RETCODE SCIPconvertCutsToConss(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   int*                  ncutsadded          /**< pointer to store number of added cuts */
   );

/** copies all active cuts from cutpool of sourcescip to constraints in targetscip
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
extern
SCIP_RETCODE SCIPcopyCuts(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global              /**< create a global or a local copy? */
   );

/** copies parameter settings from sourcescip to targetscip
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
extern
SCIP_RETCODE SCIPcopyParamSettings(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   );

/** gets depth of current scip instance (increased by each copy call) */
extern
int SCIPgetSubscipDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** copies source SCIP to target SCIP; the copying process is done in the following order:
 *  1) copy the plugins
 *  2) copy the settings
 *  3) create problem data in target-SCIP and copy the problem data of the source-SCIP
 *  4) copy all active variables
 *  5) copy all constraints
 *
 *  @note all variables and constraints which are created in the target-SCIP are not (user) captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
extern
SCIP_RETCODE SCIPcopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   const char*           suffix,             /**< suffix which will be added to the names of the source SCIP */          
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance? If TRUE, pricer
                                              *   plugins will be copied and activated, and the modifiable flag of
                                              *   constraints will be respected. If FALSE, valid will be set to FALSE, when
                                              *   there are pricers present */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid or not */
   );

/**@} */




/*
 * parameter settings
 */

/**@name Parameter Methods */
/**@{ */

/** creates a SCIP_Bool parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPaddBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Bool             defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPaddIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   int*                  valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   int                   defaultvalue,       /**< default value of the parameter */
   int                   minvalue,           /**< minimum value for parameter */
   int                   maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a SCIP_Longint parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPaddLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Longint          defaultvalue,       /**< default value of the parameter */
   SCIP_Longint          minvalue,           /**< minimum value for parameter */
   SCIP_Longint          maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a SCIP_Real parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPaddRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Real             defaultvalue,       /**< default value of the parameter */
   SCIP_Real             minvalue,           /**< minimum value for parameter */
   SCIP_Real             maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPaddCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char*                 valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   char                  defaultvalue,       /**< default value of the parameter */
   const char*           allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPaddStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char**                valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   const char*           defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** gets the value of an existing SCIP_Bool parameter */
extern
SCIP_RETCODE SCIPgetBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Int parameter */
extern
SCIP_RETCODE SCIPgetIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   int*                  value               /**< pointer to store the parameter */
   );

/** gets the value of an existing SCIP_Longint parameter */
extern
SCIP_RETCODE SCIPgetLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint*         value               /**< pointer to store the parameter */
   );

/** gets the value of an existing SCIP_Real parameter */
extern
SCIP_RETCODE SCIPgetRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Real*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Char parameter */
extern
SCIP_RETCODE SCIPgetCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   char*                 value               /**< pointer to store the parameter */
   );

/** gets the value of an existing String parameter */
extern
SCIP_RETCODE SCIPgetStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   char**                value               /**< pointer to store the parameter */
   );

/** changes the value of an existing parameter */
extern
SCIP_RETCODE SCIPsetParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   void*                 value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Bool parameter */
extern
SCIP_RETCODE SCIPchgBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Bool parameter */
extern
SCIP_RETCODE SCIPsetBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   );

/** changes the value of an existing Int parameter */
extern
SCIP_RETCODE SCIPchgIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   int                   value               /**< new value of the parameter */
   );

/** changes the value of an existing Int parameter */
extern
SCIP_RETCODE SCIPsetIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   int                   value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Longint parameter */
extern
SCIP_RETCODE SCIPchgLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Longint parameter */
extern
SCIP_RETCODE SCIPsetLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Real parameter */
extern
SCIP_RETCODE SCIPchgRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Real             value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Real parameter */
extern
SCIP_RETCODE SCIPsetRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Real             value               /**< new value of the parameter */
   );

/** changes the value of an existing Char parameter */
extern
SCIP_RETCODE SCIPchgCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   char                  value               /**< new value of the parameter */
   );

/** changes the value of an existing Char parameter */
extern
SCIP_RETCODE SCIPsetCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   char                  value               /**< new value of the parameter */
   );

/** changes the value of an existing String parameter */
extern
SCIP_RETCODE SCIPchgStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param,              /**< parameter */
   const char*           value               /**< new value of the parameter */
   );

/** changes the value of an existing String parameter */
extern
SCIP_RETCODE SCIPsetStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           value               /**< new value of the parameter */
   );

/** reads parameters from a file */
extern
SCIP_RETCODE SCIPreadParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   );

/** writes all parameters in the parameter set to a file */
extern
SCIP_RETCODE SCIPwriteParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< file name, or NULL for stdout */
   SCIP_Bool             comments,           /**< should parameter descriptions be written as comments? */
   SCIP_Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
   );

/** resets a single parameter to its default value */
extern
SCIP_RETCODE SCIPresetParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the parameter */
   );

/** resets all parameters to their default values */
extern
SCIP_RETCODE SCIPresetParams(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets parameters to 
 *  - SCIP_PARAMEMPHASIS_DEFAULT to use default values (see also SCIPresetParams())
 *  - SCIP_PARAMEMPHASIS_COUNTER to get feasible and "fast" counting process
 *  - SCIP_PARAMEMPHASIS_CPSOLVER to get CP like search (e.g. no LP relaxation)
 *  - SCIP_PARAMEMPHASIS_EASYCIP to solve easy problems fast
 *  - SCIP_PARAMEMPHASIS_FEASIBILITY to detect feasibility fast 
 *  - SCIP_PARAMEMPHASIS_HARDLP to be capable to handle hard LPs
 *  - SCIP_PARAMEMPHASIS_OPTIMALITY to prove optimality fast
 */
extern
SCIP_RETCODE SCIPsetEmphasis(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMEMPHASIS    paramemphasis,      /**< parameter emphasis */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets parameters to deactivate separators and heuristics that use auxiliary SCIP instances; should be called for
 *  auxiliary SCIP instances to avoid recursion
 */
SCIP_RETCODE SCIPsetSubscipsOff(
   SCIP*                 scip,               /**< (auxiliary) SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets heuristic parameters values to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all heuristic parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for heuristic is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the heuristic are called more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all heuristics
 */
extern
SCIP_RETCODE SCIPsetHeuristics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets presolving parameters to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all presolving parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for presolving is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the presolving is more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all presolving
 */
extern 
SCIP_RETCODE SCIPsetPresolving(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets separating parameters to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all separating parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for separating is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the separating is done more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all separating
 */
extern 
SCIP_RETCODE SCIPsetSeparating(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );
   
/** returns the array of all available SCIP parameters */
extern
SCIP_PARAM** SCIPgetParams(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the total number of all available SCIP parameters */
extern
int SCIPgetNParams(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */




/*
 * SCIP user functionality methods: managing plugins
 */

/**@name SCIP User Functionality Methods: Managing Plugins */
/**@{ */

/** creates a reader and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeReader(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of reader */
   const char*           desc,               /**< description of reader */
   const char*           extension,          /**< file extension that reader processes */
   SCIP_DECL_READERCOPY  ((*readercopy)),    /**< copy method of reader or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_READERFREE  ((*readerfree)),    /**< destructor of reader */
   SCIP_DECL_READERREAD  ((*readerread)),    /**< read method */
   SCIP_DECL_READERWRITE ((*readerwrite)),   /**< write method */
   SCIP_READERDATA*      readerdata          /**< reader data */
   );

/** returns the reader of the given name, or NULL if not existing */
extern
SCIP_READER* SCIPfindReader(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of reader */
   );

/** returns the array of currently available readers */
extern
SCIP_READER** SCIPgetReaders(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available readers */
extern
int SCIPgetNReaders(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a variable pricer and includes it in SCIP
 *  To use the variable pricer for solving a problem, it first has to be activated with a call to SCIPactivatePricer().
 *  This should be done during the problem creation stage.
 */
extern
SCIP_RETCODE SCIPincludePricer(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of variable pricer */
   const char*           desc,               /**< description of variable pricer */
   int                   priority,           /**< priority of the variable pricer */
   SCIP_Bool             delay,              /**< should the pricer be delayed until no other pricers or already existing
                                              *   problem variables with negative reduced costs are found?
                                              *   if this is set to FALSE it may happen that the pricer produces columns
                                              *   that already exist in the problem (which are also priced in by the
                                              *   default problem variable pricing in the same round)
                                              */
   SCIP_DECL_PRICERCOPY  ((*pricercopy)),    /**< copy method of variable pricer or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRICERFREE  ((*pricerfree)),    /**< destructor of variable pricer */
   SCIP_DECL_PRICERINIT  ((*pricerinit)),    /**< initialize variable pricer */
   SCIP_DECL_PRICEREXIT  ((*pricerexit)),    /**< deinitialize variable pricer */
   SCIP_DECL_PRICERINITSOL((*pricerinitsol)),/**< solving process initialization method of variable pricer */
   SCIP_DECL_PRICEREXITSOL((*pricerexitsol)),/**< solving process deinitialization method of variable pricer */
   SCIP_DECL_PRICERREDCOST((*pricerredcost)),/**< reduced cost pricing method of variable pricer for feasible LPs */
   SCIP_DECL_PRICERFARKAS((*pricerfarkas)),  /**< Farkas pricing method of variable pricer for infeasible LPs */
   SCIP_PRICERDATA*      pricerdata          /**< variable pricer data */
   );

/** returns the variable pricer of the given name, or NULL if not existing */
extern
SCIP_PRICER* SCIPfindPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of variable pricer */
   );

/** returns the array of currently available variable pricers; active pricers are in the first slots of the array */
extern
SCIP_PRICER** SCIPgetPricers(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available variable pricers */
extern
int SCIPgetNPricers(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently active variable pricers, that are used in the LP solving loop */
extern
int SCIPgetNActivePricers(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a variable pricer */
extern
SCIP_RETCODE SCIPsetPricerPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer,             /**< variable pricer */
   int                   priority            /**< new priority of the variable pricer */
   );

/** activates pricer to be used for the current problem
 *  This method should be called during the problem creation stage for all pricers that are necessary to solve
 *  the problem model.
 *  The pricers are automatically deactivated when the problem is freed.
 */
extern
SCIP_RETCODE SCIPactivatePricer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** deactivates pricer */
extern
SCIP_RETCODE SCIPdeactivatePricer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** creates a constraint handler and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of constraint handler */
   const char*           desc,               /**< description of constraint handler */
   int                   sepapriority,       /**< priority of the constraint handler for separation */
   int                   enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int                   chckpriority,       /**< priority of the constraint handler for checking feasibility (and propagation) */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int                   eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   int                   maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   SCIP_Bool             delaysepa,          /**< should separation method be delayed, if other separators found cuts? */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_Bool             delaypresol,        /**< should presolving method be delayed, if other presolvers found reductions? */
   SCIP_Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagators should be executed */
   SCIP_DECL_CONSHDLRCOPY((*conshdlrcopy)),  /**< copy method of constraint handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   SCIP_DECL_CONSINIT    ((*consinit)),      /**< initialize constraint handler */
   SCIP_DECL_CONSEXIT    ((*consexit)),      /**< deinitialize constraint handler */
   SCIP_DECL_CONSINITPRE ((*consinitpre)),   /**< presolving initialization method of constraint handler */
   SCIP_DECL_CONSEXITPRE ((*consexitpre)),   /**< presolving deinitialization method of constraint handler */
   SCIP_DECL_CONSINITSOL ((*consinitsol)),   /**< solving process initialization method of constraint handler */
   SCIP_DECL_CONSEXITSOL ((*consexitsol)),   /**< solving process deinitialization method of constraint handler */
   SCIP_DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   SCIP_DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   SCIP_DECL_CONSINITLP  ((*consinitlp)),    /**< initialize LP with relaxations of "initial" constraints */
   SCIP_DECL_CONSSEPALP  ((*conssepalp)),    /**< separate cutting planes for LP solution */
   SCIP_DECL_CONSSEPASOL ((*conssepasol)),   /**< separate cutting planes for arbitrary primal solution */
   SCIP_DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   SCIP_DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   SCIP_DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   SCIP_DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   SCIP_DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method */
   SCIP_DECL_CONSRESPROP ((*consresprop)),   /**< propagation conflict resolving method */
   SCIP_DECL_CONSLOCK    ((*conslock)),      /**< variable rounding lock method */
   SCIP_DECL_CONSACTIVE  ((*consactive)),    /**< activation notification method */
   SCIP_DECL_CONSDEACTIVE((*consdeactive)),  /**< deactivation notification method */
   SCIP_DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   SCIP_DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   SCIP_DECL_CONSDELVARS ((*consdelvars)),   /**< variable deletion method */
   SCIP_DECL_CONSPRINT   ((*consprint)),     /**< constraint display method */
   SCIP_DECL_CONSCOPY    ((*conscopy)),      /**< constraint copying method */
   SCIP_DECL_CONSPARSE   ((*consparse)),     /**< constraint parsing method */
   SCIP_DECL_CONSGETVARS ((*consgetvars)),   /**< constraint get variables method */
   SCIP_DECL_CONSGETNVARS((*consgetnvars)),  /**< constraint get number of variable method */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   );

/** returns the constraint handler of the given name, or NULL if not existing */
extern
SCIP_CONSHDLR* SCIPfindConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint handler */
   );

/** returns the array of currently available constraint handlers */
extern
SCIP_CONSHDLR** SCIPgetConshdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available constraint handlers */
extern
int SCIPgetNConshdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a conflict handler and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConflicthdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of conflict handler */
   const char*           desc,               /**< description of conflict handler */
   int                   priority,           /**< priority of the conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy)),  /**< copy method of conflict handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol)),/**< solving process initialization method of conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol)),/**< solving process deinitialization method of conflict handler */
   SCIP_DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   );

/** returns the conflict handler of the given name, or NULL if not existing */
extern
SCIP_CONFLICTHDLR* SCIPfindConflicthdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of conflict handler */
   );

/** returns the array of currently available conflict handlers */
extern
SCIP_CONFLICTHDLR** SCIPgetConflicthdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available conflict handlers */
extern
int SCIPgetNConflicthdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a conflict handler */
extern
SCIP_RETCODE SCIPsetConflicthdlrPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   int                   priority            /**< new priority of the conflict handler */
   );

/** creates a presolver and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePresol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of presolver */
   const char*           desc,               /**< description of presolver */
   int                   priority,           /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
   int                   maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
   SCIP_Bool             delay,              /**< should presolver be delayed, if other presolvers found reductions? */
   SCIP_DECL_PRESOLCOPY  ((*presolcopy)),    /**< copy method of presolver or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver to free user data (called when SCIP is exiting) */
   SCIP_DECL_PRESOLINIT  ((*presolinit)),    /**< initialization method of presolver (called after problem was transformed) */
   SCIP_DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialization method of presolver (called before transformed problem is freed) */
   SCIP_DECL_PRESOLINITPRE((*presolinitpre)),/**< presolving initialization method of presolver (called when presolving is about to begin) */
   SCIP_DECL_PRESOLEXITPRE((*presolexitpre)),/**< presolving deinitialization method of presolver (called after presolving has been finished) */
   SCIP_DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   );

/** returns the presolver of the given name, or NULL if not existing */
extern
SCIP_PRESOL* SCIPfindPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of presolver */
   );

/** returns the array of currently available presolvers */
extern
SCIP_PRESOL** SCIPgetPresols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available presolvers */
extern
int SCIPgetNPresols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a presolver */
extern
SCIP_RETCODE SCIPsetPresolPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   int                   priority            /**< new priority of the presolver */
   );

/** creates a relaxation handler and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of relaxation handler */
   const char*           desc,               /**< description of relaxation handler */
   int                   priority,           /**< priority of the relaxation handler (negative: after LP, non-negative: before LP) */
   int                   freq,               /**< frequency for calling relaxation handler */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy)),     /**< copy method of relaxation handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxation handler */
   SCIP_DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxation handler */
   SCIP_DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxation handler */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol)),  /**< solving process initialization method of relaxation handler */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol)),  /**< solving process deinitialization method of relaxation handler */
   SCIP_DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxation handler */
   SCIP_RELAXDATA*       relaxdata           /**< relaxation handler data */
   );

/** returns the relaxation handler of the given name, or NULL if not existing */
extern
SCIP_RELAX* SCIPfindRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of relaxation handler*/
   );

/** returns the array of currently available relaxation handlers  */
extern
SCIP_RELAX** SCIPgetRelaxs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available relaxation handlers  */
extern
int SCIPgetNRelaxs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a relaxation handler*/
extern
SCIP_RETCODE SCIPsetRelaxPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   int                   priority            /**< new priority of the relaxation handler */
   );

/** creates a separator and includes it in SCIP */
extern
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

/** returns the separator of the given name, or NULL if not existing */
extern
SCIP_SEPA* SCIPfindSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of separator */
   );

/** returns the array of currently available separators */
extern
SCIP_SEPA** SCIPgetSepas(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available separators */
extern
int SCIPgetNSepas(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a separator */
extern
SCIP_RETCODE SCIPsetSepaPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   int                   priority            /**< new priority of the separator */
   );

/** creates a propagator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of propagator */
   const char*           desc,               /**< description of propagator */
   int                   priority,           /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling propagator */
   SCIP_Bool             delay,              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagators should be executed */
   int                   presolpriority,     /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
   SCIP_Bool             presoldelay,        /**< should presolving be delayed, if other presolvers found reductions? */
   SCIP_DECL_PROPCOPY    ((*propcopy)),      /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PROPFREE    ((*propfree)),      /**< destructor of propagator */
   SCIP_DECL_PROPINIT    ((*propinit)),      /**< initialize propagator */
   SCIP_DECL_PROPEXIT    ((*propexit)),      /**< deinitialize propagator */
   SCIP_DECL_PROPINITPRE ((*propinitpre)),   /**< presolving initialization method of propagator */
   SCIP_DECL_PROPEXITPRE ((*propexitpre)),   /**< presolving deinitialization method of propagator */
   SCIP_DECL_PROPINITSOL ((*propinitsol)),   /**< solving process initialization method of propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol)),   /**< solving process deinitialization method of propagator */
   SCIP_DECL_PROPPRESOL  ((*proppresol)),    /**< presolving method */
   SCIP_DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop)),   /**< propagation conflict resolving method */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   );

/** returns the propagator of the given name, or NULL if not existing */
extern
SCIP_PROP* SCIPfindProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of propagator */
   );

/** returns the array of currently available propagators */
extern
SCIP_PROP** SCIPgetProps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available propagators */
extern
int SCIPgetNProps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a propagator */
extern
SCIP_RETCODE SCIPsetPropPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   int                   priority            /**< new priority of the propagator */
   );

/** sets the presolving priority of a propagator */
extern
SCIP_RETCODE SCIPsetPropPresolPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   int                   presolpriority      /**< new presol priority of the propagator */
   );


/** creates a primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of primal heuristic */
   const char*           desc,               /**< description of primal heuristic */
   char                  dispchar,           /**< display character of primal heuristic */
   int                   priority,           /**< priority of the primal heuristic */
   int                   freq,               /**< frequency for calling primal heuristic */
   int                   freqofs,            /**< frequency offset for calling primal heuristic */
   int                   maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   unsigned int          timingmask,         /**< positions in the node solving loop where heuristic should be executed;
                                              *   see definition of SCIP_HeurTiming for possible values */
   SCIP_Bool             usessubscip,        /**< does the heuristic use a secondary SCIP instance? */
   SCIP_DECL_HEURCOPY    ((*heurcopy)),      /**< copy method of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit)),      /**< initialize primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit)),      /**< deinitialize primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol)),   /**< solving process initialization method of primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol)),   /**< solving process deinitialization method of primal heuristic */
   SCIP_DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   );

/** returns the primal heuristic of the given name, or NULL if not existing */
extern
SCIP_HEUR* SCIPfindHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of primal heuristic */
   );

/** returns the array of currently available primal heuristics */
extern
SCIP_HEUR** SCIPgetHeurs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available primal heuristics */
extern
int SCIPgetNHeurs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a primal heuristic */
extern
SCIP_RETCODE SCIPsetHeurPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< primal heuristic */
   int                   priority            /**< new priority of the primal heuristic */
   );

/** creates an event handler and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of event handler */
   const char*           desc,               /**< description of event handler */
   SCIP_DECL_EVENTCOPY   ((*eventcopy)),     /**< copy method of event handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_EVENTFREE   ((*eventfree)),     /**< destructor of event handler */
   SCIP_DECL_EVENTINIT   ((*eventinit)),     /**< initialize event handler */
   SCIP_DECL_EVENTEXIT   ((*eventexit)),     /**< deinitialize event handler */
   SCIP_DECL_EVENTINITSOL((*eventinitsol)),  /**< solving process initialization method of event handler */
   SCIP_DECL_EVENTEXITSOL((*eventexitsol)),  /**< solving process deinitialization method of event handler */
   SCIP_DECL_EVENTDELETE ((*eventdelete)),   /**< free specific event data */
   SCIP_DECL_EVENTEXEC   ((*eventexec)),     /**< execute event handler */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   );

/** returns the event handler of the given name, or NULL if not existing */
extern
SCIP_EVENTHDLR* SCIPfindEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of event handler */
   );

/** returns the array of currently available event handlers */
extern
SCIP_EVENTHDLR** SCIPgetEventhdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available event handlers */
extern
int SCIPgetNEventhdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a node selector and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of node selector */
   const char*           desc,               /**< description of node selector */
   int                   stdpriority,        /**< priority of the node selector in standard mode */
   int                   memsavepriority,    /**< priority of the node selector in memory saving mode */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy)),   /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_NODESELFREE ((*nodeselfree)),   /**< destructor of node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit)),   /**< initialize node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit)),   /**< deinitialize node selector */
   SCIP_DECL_NODESELINITSOL((*nodeselinitsol)),/**< solving process initialization method of node selector */
   SCIP_DECL_NODESELEXITSOL((*nodeselexitsol)),/**< solving process deinitialization method of node selector */
   SCIP_DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   SCIP_DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   );

/** returns the node selector of the given name, or NULL if not existing */
extern
SCIP_NODESEL* SCIPfindNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of node selector */
   );

/** returns the array of currently available node selectors */
extern
SCIP_NODESEL** SCIPgetNodesels(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available node selectors */
extern
int SCIPgetNNodesels(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a node selector in standard mode */
extern
SCIP_RETCODE SCIPsetNodeselStdPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   int                   priority            /**< new standard priority of the node selector */
   );

/** sets the priority of a node selector in memory saving mode */
extern
SCIP_RETCODE SCIPsetNodeselMemsavePriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   int                   priority            /**< new memory saving priority of the node selector */
   );

/** returns the currently used node selector */
extern
SCIP_NODESEL* SCIPgetNodesel(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of branching rule */
   const char*           desc,               /**< description of branching rule */
   int                   priority,           /**< priority of the branching rule */
   int                   maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying branching rule
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_DECL_BRANCHCOPY  ((*branchcopy)),    /**< copy method of branching rule or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   SCIP_DECL_BRANCHINIT  ((*branchinit)),    /**< initialize branching rule */
   SCIP_DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialize branching rule */
   SCIP_DECL_BRANCHINITSOL((*branchinitsol)),/**< solving process initialization method of branching rule */
   SCIP_DECL_BRANCHEXITSOL((*branchexitsol)),/**< solving process deinitialization method of branching rule */
   SCIP_DECL_BRANCHEXECLP((*branchexeclp)),  /**< branching execution method for fractional LP solutions */
   SCIP_DECL_BRANCHEXECEXT((*branchexecext)),/**< branching execution method for relaxation solutions */
   SCIP_DECL_BRANCHEXECPS((*branchexecps)),  /**< branching execution method for not completely fixed pseudo solutions */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   );

/** returns the branching rule of the given name, or NULL if not existing */
extern
SCIP_BRANCHRULE* SCIPfindBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of branching rule */
   );

/** returns the array of currently available branching rules */
extern
SCIP_BRANCHRULE** SCIPgetBranchrules(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available branching rules */
extern
int SCIPgetNBranchrules(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a branching rule */
extern
SCIP_RETCODE SCIPsetBranchrulePriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   int                   priority            /**< new priority of the branching rule */
   );

/** sets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
extern
SCIP_RETCODE SCIPsetBranchruleMaxdepth(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   int                   maxdepth            /**< new maxdepth of the branching rule */
   );

/** sets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
extern
SCIP_RETCODE SCIPsetBranchruleMaxbounddist(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_Real             maxbounddist        /**< new maxbounddist of the branching rule */
   );

/** creates a display column and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of display column */
   const char*           desc,               /**< description of display column */
   const char*           header,             /**< head line of display column */
   SCIP_DISPSTATUS       dispstatus,         /**< display activation status of display column */
   SCIP_DECL_DISPCOPY    ((*dispcopy)),      /**< copy method of display column or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DISPFREE    ((*dispfree)),      /**< destructor of display column */
   SCIP_DECL_DISPINIT    ((*dispinit)),      /**< initialize display column */
   SCIP_DECL_DISPEXIT    ((*dispexit)),      /**< deinitialize display column */
   SCIP_DECL_DISPINITSOL ((*dispinitsol)),   /**< solving process initialization method of display column */
   SCIP_DECL_DISPEXITSOL ((*dispexitsol)),   /**< solving process deinitialization method of display column */
   SCIP_DECL_DISPOUTPUT  ((*dispoutput)),    /**< output method */
   SCIP_DISPDATA*        dispdata,           /**< display column data */
   int                   width,              /**< width of display column (no. of chars used) */
   int                   priority,           /**< priority of display column */
   int                   position,           /**< relative position of display column */
   SCIP_Bool             stripline           /**< should the column be separated with a line from its right neighbor? */
   );

/** returns the display column of the given name, or NULL if not existing */
extern
SCIP_DISP* SCIPfindDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of display column */
   );

/** returns the array of currently available display columns */
extern
SCIP_DISP** SCIPgetDisps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available display columns */
extern
int SCIPgetNDisps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** automatically selects display columns for being shown w.r.t. the display width parameter */
extern
SCIP_RETCODE SCIPautoselectDisps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes an NLPI in SCIP */
extern
SCIP_RETCODE SCIPincludeNlpi(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi                /**< NLPI data structure */
   );

/** returns the NLPI of the given name, or NULL if not existing */
extern
SCIP_NLPI* SCIPfindNlpi(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of NLPI */
   );

/** returns the array of currently available NLPIs (sorted by priority) */
extern
SCIP_NLPI** SCIPgetNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available NLPIs */
extern
int SCIPgetNNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of an NLPI */
extern
SCIP_RETCODE SCIPsetNlpiPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< NLPI */
   int                   priority            /**< new priority of the NLPI */
   );

/** includes information about an external code linked into the SCIP library */
extern
SCIP_RETCODE SCIPincludeExternalCodeInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of external code */
   const char*           description         /**< description of external code, or NULL */
   );

/** returns an array of names of currently included external codes */
extern
char** SCIPgetExternalCodeNames(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns an array of the descriptions of currently included external codes
 *
 *  @note some descriptions may be NULL
 */
extern
char** SCIPgetExternalCodeDescriptions(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently included information on external codes */
extern
int SCIPgetNExternalCodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** prints information on external codes to a file stream via the message handler system
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
extern
void SCIPprintExternalCodes(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/**@} */




/*
 * user interactive dialog methods
 */

/**@name User Interactive Dialog Methods */
/**@{ */

/** creates and includes dialog */
extern
SCIP_RETCODE SCIPincludeDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         dialog,             /**< pointer to store the dialog */
   SCIP_DECL_DIALOGCOPY  ((*dialogcopy)),    /**< copy method of dialog or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DIALOGEXEC  ((*dialogexec)),    /**< execution method of dialog */
   SCIP_DECL_DIALOGDESC  ((*dialogdesc)),    /**< description output method of dialog, or NULL */
   SCIP_DECL_DIALOGFREE  ((*dialogfree)),    /**< destructor of dialog to free user data, or NULL */
   const char*           name,               /**< name of dialog: command name appearing in parent's dialog menu */
   const char*           desc,               /**< description of dialog used if description output method is NULL */
   SCIP_Bool             issubmenu,          /**< is the dialog a submenu? */
   SCIP_DIALOGDATA*      dialogdata          /**< user defined dialog data */
   );

/** returns if the dialog already exists */
extern
SCIP_Bool SCIPexistsDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** captures a dialog */
extern
SCIP_RETCODE SCIPcaptureDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** releases a dialog */
extern
SCIP_RETCODE SCIPreleaseDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         dialog              /**< pointer to the dialog */
   );

/** makes given dialog the root dialog of SCIP's interactive user shell; captures dialog and releases former root dialog */
extern
SCIP_RETCODE SCIPsetRootDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog to be the root */
   );

/** returns the root dialog of SCIP's interactive user shell */
extern
SCIP_DIALOG* SCIPgetRootDialog(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds a sub dialog to the given dialog as menu entry and captures it */
extern
SCIP_RETCODE SCIPaddDialogEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog,             /**< dialog to extend, or NULL for root dialog */
   SCIP_DIALOG*          subdialog           /**< subdialog to add as menu entry in dialog */
   );

/** adds a single line of input which is treated as if the user entered the command line */
extern
SCIP_RETCODE SCIPaddDialogInputLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           inputline           /**< input line to add */
   );

/** adds a single line of input to the command history which can be accessed with the cursor keys */
extern
SCIP_RETCODE SCIPaddDialogHistoryLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           inputline           /**< input line to add */
   );

/** starts interactive mode of SCIP by executing the root dialog */
extern
SCIP_RETCODE SCIPstartInteraction(
   SCIP*                 scip                /**< SCIP data structure */
   );
   
/**@} */




/*
 * global problem methods
 */

/**@name Global Problem Methods */
/**@{ */

/** creates empty problem and initializes all solving data structures (the objective sense is set to MINIMIZE)
 *  If the problem type requires the use of variable pricers, these pricers should be added to the problem with calls
 *  to SCIPactivatePricer(). These pricers are automatically deactivated, when the problem is freed.
 */
extern
SCIP_RETCODE SCIPcreateProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */
   SCIP_DECL_PROBDELORIG ((*probdelorig)),   /**< frees user data of original problem */
   SCIP_DECL_PROBTRANS   ((*probtrans)),     /**< creates user data of transformed problem by transforming original user data */
   SCIP_DECL_PROBDELTRANS((*probdeltrans)),  /**< frees user data of transformed problem */
   SCIP_DECL_PROBINITSOL ((*probinitsol)),   /**< solving process initialization method of transformed data */
   SCIP_DECL_PROBEXITSOL ((*probexitsol)),   /**< solving process deinitialization method of transformed data */
   SCIP_DECL_PROBCOPY    ((*probcopy)),      /**< copies user data if you want to copy it to a subscip, or NULL */
   SCIP_PROBDATA*        probdata            /**< user problem data set by the reader */
   );

/** reads problem from file and initializes all solving data structures */
extern
SCIP_RETCODE SCIPreadProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< problem file name */
   const char*           extension           /**< extension of the desired file reader, 
                                              *   or NULL if file extension should be used */
   );

/** writes original problem to file  */
extern
SCIP_RETCODE SCIPwriteOrigProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< output file (or NULL for standard output) */
   const char*           extension,          /**< extension of the desired file reader, 
                                              *   or NULL if file extension should be used */
   SCIP_Bool             genericnames        /**< use generic variable and constraint names? */
   );

/** writes transformed problem which are valid in the current node to file */
extern
SCIP_RETCODE SCIPwriteTransProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< output file (or NULL for standard output) */
   const char*           extension,          /**< extension of the desired file reader, 
                                              *   or NULL if file extension should be used */
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   );

/** frees problem and solution process data */
extern
SCIP_RETCODE SCIPfreeProb(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** permutes parts of the problem data structure */
extern
SCIP_RETCODE SCIPpermuteProb(
   SCIP*                 scip,              /**< SCIP data structure */
   unsigned int          randseed,          /**< seed value for random generator */
   SCIP_Bool             permuteconss,      /**< should the list of constraints in each constraint handler be permuted? */
   SCIP_Bool             permutebinvars,    /**< should the list of binary variables be permuted? */
   SCIP_Bool             permuteintvars,    /**< should the list of integer variables be permuted? */
   SCIP_Bool             permuteimplvars,   /**< should the list of implicit integer variables be permuted? */
   SCIP_Bool             permutecontvars    /**< should the list of continuous integer variables be permuted? */
   );

/** gets user problem data */
extern
SCIP_PROBDATA* SCIPgetProbData(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets user problem data */
extern
SCIP_RETCODE SCIPsetProbData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< user problem data to use */
   );

/** gets name of the current problem instance */
extern
const char* SCIPgetProbName(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets name of the current problem instance */
extern
SCIP_RETCODE SCIPsetProbName(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name to be set */
   );

/** gets objective sense of original problem */
extern
SCIP_OBJSENSE SCIPgetObjsense(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets objective sense of problem */
extern
SCIP_RETCODE SCIPsetObjsense(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OBJSENSE         objsense            /**< new objective sense */
   );

/** adds offset of objective function */
extern
SCIP_RETCODE SCIPaddObjoffset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             addval              /**< value to add to objective offset */
   );

/** returns the objective offset of the original problem */
extern
SCIP_Real SCIPgetOrigObjoffset(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the objective scale of the original problem */
extern
SCIP_Real SCIPgetOrigObjscale(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the objective offset of the transformed problem */
extern
SCIP_Real SCIPgetTransObjoffset(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the objective scale of the transformed problem */
extern
SCIP_Real SCIPgetTransObjscale(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets limit on objective function, such that only solutions better than this limit are accepted */
extern
SCIP_RETCODE SCIPsetObjlimit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             objlimit            /**< new primal objective limit */
   );

/** gets current limit on objective function */
extern
SCIP_Real SCIPgetObjlimit(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** informs SCIP, that the objective value is always integral in every feasible solution */
extern
SCIP_RETCODE SCIPsetObjIntegral(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the objective value is known to be integral in every feasible solution */
extern
SCIP_Bool SCIPisObjIntegral(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the Euclidean norm of the objective function vector (available only for transformed problem) */
extern
SCIP_Real SCIPgetObjNorm(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds variable to the problem */
extern
SCIP_RETCODE SCIPaddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to add */
   );

/** adds variable to the problem and uses it as pricing candidate to enter the LP */
extern
SCIP_RETCODE SCIPaddPricedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             score               /**< pricing score of variable (the larger, the better the variable) */
   );

/** removes variable from the problem */
extern
SCIP_RETCODE SCIPdelVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to delete */
   SCIP_Bool*            deleted             /**< pointer to store whether variable was successfully marked to be deleted */
   );

/** gets variables of the problem along with the numbers of different variable types; data may become invalid after
 *  calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 */
extern
SCIP_RETCODE SCIPgetVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*                  nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*                  nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*                  nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*                  nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*                  ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   );

/** gets array with active problem variables
 *
 *  @warning If your are using the methods which add or change bound of variables (e.g., SCIPchgVarType(), SCIPfixVar(),
 *           SCIPaggregateVars(), and SCIPmultiaggregateVar()), it can happen that the internal variable array (which is
 *           accessed via this method) gets resized and/or resorted. This can invalid the data pointer which is returned
 *           by this method.
 */
extern
SCIP_VAR** SCIPgetVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of active problem variables */
extern
int SCIPgetNVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of binary active problem variables */
extern
int SCIPgetNBinVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of integer active problem variables */
extern
int SCIPgetNIntVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of implicit integer active problem variables */
extern
int SCIPgetNImplVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of continuous active problem variables */
extern
int SCIPgetNContVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of active problem variables with a non-zero objective coefficient
 *
 *  @note In case of the original problem the number of variables is counted. In case of the transformed problem the
 *        number of variables is just return since it is stored
 */
extern
int SCIPgetNObjVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets array with fixed and aggregated problem variables; data may become invalid after
 *  calls to SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 */
extern
SCIP_VAR** SCIPgetFixedVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of fixed or aggregated problem variables */
extern
int SCIPgetNFixedVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets variables of the original problem along with the numbers of different variable types; data may become invalid
 *  after a call to SCIPchgVarType()
 */
extern
SCIP_RETCODE SCIPgetOrigVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*                  nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*                  nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*                  nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*                  nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*                  ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   );

/** gets array with original problem variables; data may become invalid after
 *  a call to SCIPchgVarType()
 */
extern
SCIP_VAR** SCIPgetOrigVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of original problem variables */
extern
int SCIPgetNOrigVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of binary original problem variables */
extern
int SCIPgetNOrigBinVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of integer original problem variables */
extern
int SCIPgetNOrigIntVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of implicit integer original problem variables */
extern
int SCIPgetNOrigImplVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of continuous original problem variables */
extern
int SCIPgetNOrigContVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of all problem variables created during creation and solving of problem;
 *  this includes also variables that were deleted in the meantime
 */
extern
int SCIPgetNTotalVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets variables of the original or transformed problem along with the numbers of different variable types;
 *  the returned problem space (original or transformed) corresponds to the given solution;
 *  data may become invalid after calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and 
 *  SCIPmultiaggregateVar()
 */
extern
SCIP_RETCODE SCIPgetSolVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that selects the problem space, NULL for current solution */
   SCIP_VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*                  nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*                  nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*                  nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*                  nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*                  ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   );

/** returns variable of given name in the problem, or NULL if not existing */
extern
SCIP_VAR* SCIPfindVar(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of variable to find */
   );

/** returns TRUE iff all potential variables exist in the problem, and FALSE, if there may be additional variables,
 *  that will be added in pricing and improve the objective value
 */
extern
SCIP_Bool SCIPallVarsInProb(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds constraint to the problem; if constraint is only valid locally, it is added to the local subproblem of the
 *  current node (and all of its subnodes); otherwise it is added to the global problem;
 *  if a local constraint is added at the root node, it is automatically upgraded into a global constraint
 */
extern
SCIP_RETCODE SCIPaddCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to add */
   );

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was added, or from the problem, if it was a problem constraint
 */
extern
SCIP_RETCODE SCIPdelCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to delete */
   );

/** returns original constraint of given name in the problem, or NULL if not existing */
extern
SCIP_CONS* SCIPfindOrigCons(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint to find */
   );

/** returns constraint of given name in the problem, or NULL if not existing */
extern
SCIP_CONS* SCIPfindCons(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint to find */
   );

/** gets number of upgraded constraints */
extern
int SCIPgetNUpgrConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of globally valid constraints currently in the problem */
extern
int SCIPgetNConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets array of globally valid constraints currently in the problem
 *
 *  @warning If your are using the method SCIPaddCons(), it can happen that the internal constraint array (which is
 *           accessed via this method) gets resized. This can invalid the data pointer which is returned by this method.
 */
extern
SCIP_CONS** SCIPgetConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of constraints in the original problem */
extern
int SCIPgetNOrigConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets array of constraints in the original problem */
extern
SCIP_CONS** SCIPgetOrigConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */




/*
 * local subproblem methods
 */

/**@name Local Subproblem Methods */
/**@{ */

/** adds constraint to the given node (and all of its subnodes), even if it is a global constraint;
 *  It is sometimes desirable to add the constraint to a more local node (i.e., a node of larger depth) even if
 *  the constraint is also valid higher in the tree, for example, if one wants to produce a constraint which is
 *  only active in a small part of the tree although it is valid in a larger part.
 *  In this case, one should pass the more global node where the constraint is valid as "validnode".
 *  Note that the same constraint cannot be added twice to the branching tree with different "validnode" parameters.
 *  If the constraint is valid at the same node as it is inserted (the usual case), one should pass NULL as "validnode".
 *  If the "validnode" is the root node, it is automatically upgraded into a global constraint, but still only added to
 *  the given node. If a local constraint is added to the root node, it is added to the global problem instead.
 */
extern
SCIP_RETCODE SCIPaddConsNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to add constraint to */
   SCIP_CONS*            cons,               /**< constraint to add */
   SCIP_NODE*            validnode           /**< node at which the constraint is valid, or NULL */
   );

/** adds constraint locally to the current node (and all of its subnodes), even if it is a global constraint;
 *  It is sometimes desirable to add the constraint to a more local node (i.e., a node of larger depth) even if
 *  the constraint is also valid higher in the tree, for example, if one wants to produce a constraint which is
 *  only active in a small part of the tree although it is valid in a larger part.
 *
 *  If the constraint is valid at the same node as it is inserted (the usual case), one should pass NULL as "validnode".
 *  If the "validnode" is the root node, it is automatically upgraded into a global constraint, but still only added to
 *  the given node. If a local constraint is added to the root node, it is added to the global problem instead.
 *
 *  @note The same constraint cannot be added twice to the branching tree with different "validnode" parameters. This is
 *        the case due internal data structures and performance issues. In such a case you should try to realize your
 *        issue using the method SCIPdisableCons() and SCIPenableCons() and control these via the event system of SCIP.
 */
extern
SCIP_RETCODE SCIPaddConsLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to add */
   SCIP_NODE*            validnode           /**< node at which the constraint is valid, or NULL */
   );

/** disables constraint's separation, enforcing, and propagation capabilities at the given node (and all subnodes);
 *  if the method is called at the root node, the constraint is globally deleted from the problem;
 *  the constraint deletion is being remembered at the given node, s.t. after leaving the node's subtree, the constraint
 *  is automatically enabled again, and after entering the node's subtree, it is automatically disabled;
 *  this may improve performance because redundant checks on this constraint are avoided, but it consumes memory;
 *  alternatively, use SCIPdisableCons()
 */
extern
SCIP_RETCODE SCIPdelConsNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to disable constraint in */
   SCIP_CONS*            cons                /**< constraint to locally delete */
   );

/** disables constraint's separation, enforcing, and propagation capabilities at the current node (and all subnodes);
 *  if the method is called during problem modification or at the root node, the constraint is globally deleted from
 *  the problem;
 *  the constraint deletion is being remembered at the current node, s.t. after leaving the current subtree, the
 *  constraint is automatically enabled again, and after reentering the current node's subtree, it is automatically
 *  disabled again;
 *  this may improve performance because redundant checks on this constraint are avoided, but it consumes memory;
 *  alternatively, use SCIPdisableCons()
 */
extern
SCIP_RETCODE SCIPdelConsLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to locally delete */
   );

/** gets estimate of best primal solution w.r.t. original problem contained in current subtree */
extern
SCIP_Real SCIPgetLocalOrigEstimate(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets estimate of best primal solution w.r.t. transformed problem contained in current subtree */
extern
SCIP_Real SCIPgetLocalTransEstimate(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets dual bound of current node */
extern
SCIP_Real SCIPgetLocalDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets lower bound of current node in transformed problem */
extern
SCIP_Real SCIPgetLocalLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets dual bound of given node */
extern
SCIP_Real SCIPgetNodeDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node to get dual bound for */
   );

/** gets lower bound of given node in transformed problem */
extern
SCIP_Real SCIPgetNodeLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node to get dual bound for */
   );

/** if given value is tighter (larger for minimization, smaller for maximization) than the current node's dual bound,
 *  sets the current node's dual bound to the new value
 */
extern
SCIP_RETCODE SCIPupdateLocalDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   );

/** if given value is larger than the current node's lower bound (in transformed problem), sets the current node's
 *  lower bound to the new value
 */
extern
SCIP_RETCODE SCIPupdateLocalLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   );

/** if given value is tighter (larger for minimization, smaller for maximization) than the node's dual bound,
 *  sets the node's dual bound to the new value
 */
extern
SCIP_RETCODE SCIPupdateNodeDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to update dual bound for */
   SCIP_Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   );

/** if given value is larger than the node's lower bound (in transformed problem), sets the node's lower bound
 *  to the new value
 */
extern
SCIP_RETCODE SCIPupdateNodeLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to update lower bound for */
   SCIP_Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   );

/** change the node selection priority of the given child */
extern 
SCIP_RETCODE SCIPchgChildPrio(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            child,              /**< child to update the node selection priority */
   SCIP_Real             priority            /**< node selection priority value */
   );

/**@} */




/*
 * solve methods
 */

/**@name Solve Methods */
/**@{ */

/** initializes solving data structures and transforms problem */
extern
SCIP_RETCODE SCIPtransformProb(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms and presolves problem */
extern
SCIP_RETCODE SCIPpresolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** transforms, presolves, and solves problem */
extern
SCIP_RETCODE SCIPsolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** frees branch and bound tree and all solution process data; statistics, presolving data and transformed problem is
 *  preserved
 */
extern
SCIP_RETCODE SCIPfreeSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             restart             /**< should certain data be preserved for improved restarting? */
   );

/** frees all solution process data including presolving and transformed problem, only original problem is kept */
extern
SCIP_RETCODE SCIPfreeTransform(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** interrupts solving process as soon as possible (e.g., after the current node has been solved) */
extern
SCIP_RETCODE SCIPinterruptSolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** restarts solving process as soon as possible (e.g., after the current node has been solved) */
extern
SCIP_RETCODE SCIPrestartSolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** whether we are in the restarting phase */
extern
SCIP_Bool SCIPisInRestart(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */




/*
 * variable methods
 */

/**@name Variable Methods */
/**@{ */

/** creates and captures problem variable; if variable is of integral type, fractional bounds are automatically rounded;
 *  an integer variable with bounds zero and one is automatically converted into a binary variable;
 *
 *  @warning When doing column generation and the original problem is a maximization problem, notice that SCIP will
 *           transform the problem into a minimization problem by multiplying the objective function by -1.  Thus, the
 *           original objective function value of variables created during the solving process has to be multiplied by
 *           -1, too.
 *
 *  @note the variable gets captured, hence at one point you have to release it using the method SCIPreleaseVar()
 */
extern
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

/** outputs the variable name to the file stream */
extern
SCIP_RETCODE SCIPwriteVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR*             var,                /**< variable array to output */
   SCIP_Bool             type                /**< should the variable type be also posted */
   );

/** print the given list of variables to output stream separated by the given delimiter character;
 *
 *  i. e. the variables x1, x2, ..., xn with given delimiter ',' are written as: \<x1\>, \<x2\>, ..., \<xn\>;
 *
 *  the method SCIPparseVarsList() can parse such a string
 *
 *  @note The printing process is done via the message handler system.
 */
extern
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
 *  @note The printing process is done via the message handler system.
 */
extern
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
 *  @note The printing process is done via the message handler system.
 */
extern
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
 */
extern
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
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   );

/** parses the given string for a variable name and stores the variable in the corresponding pointer if such a variable
 *  exits and returns the position where the parsing stopped 
 */
extern
SCIP_RETCODE SCIPparseVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            var,                /**< pointer to store the problem variable, or NULL if it does not exit */
   char**                endptr              /**< pointer to store the final string position if successfully */
   );

/** parse the given string as variable list (here ',' is the delimiter)) (\<x1\>, \<x2\>, ..., \<xn\>) (see
 *  SCIPwriteVarsList() ); if it was successful, the pointer success is set to TRUE
 *
 *  @note the pointer success in only set to FALSE in the case that a variable with a parsed variable name does not exist
 *
 *  @note If the number of (parsed) variables is greater than the available slots in the variable array, nothing happens
 *        except that the required size is stored in the corresponding integer; the reason for this approach is that we
 *        cannot reallocate memory, since we do not know how the memory has been allocated (e.g., by a C++ 'new' or SCIP
 *        memory functions).
 */
extern
SCIP_RETCODE SCIPparseVarsList(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            vars,               /**< array to store the parsed variable */
   int*                  nvars,              /**< pointer to store number of parsed variables */
   int                   varssize,           /**< size of the variable array */
   int*                  requiredsize,       /**< pointer to store the required array size for the active variables */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   char                  delimiter,          /**< character which is used for delimitation */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successfully or not */
   );

/** parse the given string as linear sum of variables and coefficients (c1 \<x1\> + c2 \<x2\> + ... + cn \<xn\>)
 *  (see SCIPwriteVarsLinearsum() ); if it was successful, the pointer success is set to TRUE
 *
 *  @note the pointer success in only set to FALSE in the case that a variable with a parsed variable name does not exist
 *
 *  @note If the number of (parsed) variables is greater than the available slots in the variable array, nothing happens
 *        except that the required size is stored in the corresponding integer; the reason for this approach is that we
 *        cannot reallocate memory, since we do not know how the memory has been allocated (e.g., by a C++ 'new' or SCIP
 *        memory functions).
 */
extern
SCIP_RETCODE SCIPparseVarsLinearsum(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            vars,               /**< array to store the parsed variables */
   SCIP_Real*            vals,               /**< array to store the parsed coefficients */
   int*                  nvars,              /**< pointer to store number of parsed variables */
   int                   varssize,           /**< size of the variable array */
   int*                  requiredsize,       /**< pointer to store the required array size for the active variables */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successfully or not */
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
 */
extern
SCIP_RETCODE SCIPparseVarsPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR****          monomialvars,       /**< pointer to store arrays with variables for each monomial */
   SCIP_Real***          monomialexps,       /**< pointer to store arrays with variable exponents */
   SCIP_Real**           monomialcoefs,      /**< pointer to store array with monomial coefficients */
   int**                 monomialnvars,      /**< pointer to store array with number of variables for each monomial */
   int*                  nmonomials,         /**< pointer to store number of parsed monomials */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successfully or not */
   );

/** frees memory allocated when parsing a polynomial from a string */
extern
void SCIPfreeParseVarsPolynomialData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR****          monomialvars,       /**< pointer to store arrays with variables for each monomial */
   SCIP_Real***          monomialexps,       /**< pointer to store arrays with variable exponents */
   SCIP_Real**           monomialcoefs,      /**< pointer to store array with monomial coefficients */
   int**                 monomialnvars,      /**< pointer to store array with number of variables for each monomial */
   int                   nmonomials          /**< pointer to store number of parsed monomials */
   );

/** increases usage counter of variable */
extern
SCIP_RETCODE SCIPcaptureVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to capture */
   );

/** decreases usage counter of variable, if the usage pointer reaches zero the variable gets freed
 *
 *  @note the pointer of the variable will be NULLed
 */
extern
SCIP_RETCODE SCIPreleaseVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var                 /**< pointer to variable */
   );

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 */
extern
SCIP_RETCODE SCIPtransformVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get/create transformed variable for */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable */
   );

/** gets and captures transformed variables for an array of variables;
 *  if a variable of the array is not yet transformed, a new transformed variable for this variable is created;
 *  it is possible to call this method with vars == transvars
 */
extern
SCIP_RETCODE SCIPtransformVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get/create transformed variables for */
   SCIP_VAR**            vars,               /**< array with variables to get/create transformed variables for */
   SCIP_VAR**            transvars           /**< array to store the transformed variables */
   );

/** gets corresponding transformed variable of a given variable;
 *  returns NULL as transvar, if transformed variable is not yet existing
 */
extern
SCIP_RETCODE SCIPgetTransformedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get transformed variable for */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable */
   );

/** gets corresponding transformed variables for an array of variables;
 *  stores NULL in a transvars slot, if the transformed variable is not yet existing;
 *  it is possible to call this method with vars == transvars, but remember that variables that are not
 *  yet transformed will be replaced with NULL
 */
extern
SCIP_RETCODE SCIPgetTransformedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get transformed variables for */
   SCIP_VAR**            vars,               /**< array with variables to get transformed variables for */
   SCIP_VAR**            transvars           /**< array to store the transformed variables */
   );

/** gets negated variable x' = lb + ub - x of variable x; negated variable is created, if not yet existing */
extern
SCIP_RETCODE SCIPgetNegatedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get negated variable for */
   SCIP_VAR**            negvar              /**< pointer to store the negated variable */
   );

/** gets negated variables x' = lb + ub - x of variables x; negated variables are created, if not yet existing */
extern
SCIP_RETCODE SCIPgetNegatedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get negated variables for */
   SCIP_VAR**            vars,               /**< array of variables to get negated variables for */
   SCIP_VAR**            negvars             /**< array to store the negated variables */
   );

/** gets a binary variable that is equal to the given binary variable, and that is either active, fixed, or 
 *  multi-aggregated, or the negated variable of an active, fixed, or multi-aggregated variable
 */
extern
SCIP_RETCODE SCIPgetBinvarRepresentative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to get binary representative for */
   SCIP_VAR**            repvar,             /**< pointer to store the binary representative */
   SCIP_Bool*            negated             /**< pointer to store whether the negation of an active variable was returned */
   );

/** gets binary variables that are equal to the given binary variables, and which are either active, fixed, or 
 *  multi-aggregated, or the negated variables of active, fixed, or multi-aggregated variables
 */
extern
SCIP_RETCODE SCIPgetBinvarRepresentatives(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of binary variables to get representatives for */
   SCIP_VAR**            vars,               /**< binary variables to get binary representatives for */
   SCIP_VAR**            repvars,            /**< array to store the binary representatives */
   SCIP_Bool*            negated             /**< array to store whether the negation of an active variable was returned */
   );

/** flattens aggregation graph of multi-aggregated variable in order to avoid exponential recursion later on */
extern
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
 *  @note The resulting linear sum is stored into the given variable array, scalar array, and constant. That means the
 *        given entries are overwritten.
 *
 *  @note That method can be used to convert a single variables into variable space of active variables. Therefore call
 *        the method with the linear sum 1.0*x + 0.0.
 */
extern
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

/** return for given variables all their active counterparts; all active variables will be pairwise different
 *  @note It does not hold that the first output variable is the active variable for the first input variable.
 */
extern
SCIP_RETCODE SCIPgetActiveVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array with given variables and as output all active
					      *   variables, if enough slots exist
					      */
   int*                  nvars,              /**< number of given variables, and as output number of active variables,
					      *   if enough slots exist
					      */
   int                   varssize,           /**< available slots in vars array */
   int*                  requiredsize        /**< pointer to store the required array size for the active variables */
   );

/** returns the reduced costs of the variable in the current node's LP relaxation;
 *  the current node has to have a feasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 */
extern
SCIP_Real SCIPgetVarRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get reduced costs, should be a column in current node LP */
   );

/** returns the Farkas coefficient of the variable in the current node's LP relaxation;
 *  the current node has to have an infeasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 */
extern
SCIP_Real SCIPgetVarFarkasCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get reduced costs, should be a column in current node LP */
   );

/** gets solution value for variable in current node */
extern
SCIP_Real SCIPgetVarSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get solution value for */
   );

/** gets solution values of multiple variables in current node */
extern
SCIP_RETCODE SCIPgetVarSols(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get solution value for */
   SCIP_VAR**            vars,               /**< array with variables to get value for */
   SCIP_Real*            vals                /**< array to store solution values of variables */
   );

/** sets the solution value of all variables in the global relaxation solution to zero */
extern
SCIP_RETCODE SCIPclearRelaxSolVals(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the value of the given variable in the global relaxation solution;
 *  this solution can be filled by the relaxation handlers and can be used by heuristics and for separation;
 *  You can use SCIPclearRelaxSolVals() to set all values to zero, initially;
 *  after setting all solution values, you have to call SCIPmarkRelaxSolValid() 
 *  to inform SCIP that the stored solution is valid
 */
extern
SCIP_RETCODE SCIPsetRelaxSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to set value for */
   SCIP_Real             val                 /**< solution value of variable */
   );

/** sets the values of the given variables in the global relaxation solution;
 *  this solution can be filled by the relaxation handlers  and can be used by heuristics and for separation;
 *  the solution is automatically cleared, s.t. all other variables get value 0.0
 */
extern
SCIP_RETCODE SCIPsetRelaxSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to set relaxation solution value for */
   SCIP_VAR**            vars,               /**< array with variables to set value for */
   SCIP_Real*            vals                /**< array with solution values of variables */
   );

/** sets the values of the variables in the global relaxation solution to the values 
 *  in the given primal solution; the relaxation solution can be filled by the relaxation hanlders
 *  and might be used by heuristics and for separation
 */
extern
SCIP_RETCODE SCIPsetRelaxSolValsSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal relaxation solution */ 
   );

/** returns whether the relaxation solution is valid */
extern
SCIP_Bool SCIPisRelaxSolValid(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** informs SCIP, that the relaxation solution is valid */
extern
SCIP_RETCODE SCIPmarkRelaxSolValid(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** informs SCIP, that the relaxation solution is invalid */
extern
SCIP_RETCODE SCIPmarkRelaxSolInvalid(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the relaxation solution value of the given variable */
extern
SCIP_Real SCIPgetRelaxSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get value for */
   );

/** gets the relaxation solution objective value */
extern
SCIP_Real SCIPgetRelaxSolObj(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** start strong branching - call before any strong branching */
extern
SCIP_RETCODE SCIPstartStrongbranch(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** end strong branching - call after any strong branching */
extern
SCIP_RETCODE SCIPendStrongbranch(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets strong branching information on column variable with fractional value */
extern
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

/** gets strong branching information on column variable with integral value */
extern
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

/** gets strong branching information on column variables with fractional values */
extern
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

/** gets strong branching information on column variables with integral values */
extern
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

/** gets strong branching information on COLUMN variable of the last SCIPgetVarStrongbranch() call;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given variable;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 */
extern
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

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given variable, or -1 if strong branching was never applied to the variable in current run
 */
extern
SCIP_Longint SCIPgetVarStrongbranchNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get last strong branching node for */
   );

/** if strong branching was already applied on the variable at the current node, returns the number of LPs solved after
 *  the LP where the strong branching on this variable was applied;
 *  if strong branching was not yet applied on the variable at the current node, returns INT_MAX
 */
extern
SCIP_Longint SCIPgetVarStrongbranchLPAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get strong branching LP age for */
   );

/** gets number of times, strong branching was applied in current run on the given variable */
extern
int SCIPgetVarNStrongbranchs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get last strong branching node for */
   );

/** adds given values to lock numbers of variable for rounding */
extern
SCIP_RETCODE SCIPaddVarLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   nlocksdown,         /**< modification in number of rounding down locks */
   int                   nlocksup            /**< modification in number of rounding up locks */
   );

/** locks rounding of variable with respect to the lock status of the constraint and its negation;
 *  this method should be called whenever the lock status of a variable in a constraint changes, for example if
 *  the coefficient of the variable changed its sign or if the left or right hand sides of the constraint were
 *  added or removed
 */
extern
SCIP_RETCODE SCIPlockVarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             lockdown,           /**< should the rounding be locked in downwards direction? */
   SCIP_Bool             lockup              /**< should the rounding be locked in upwards direction? */
   );

/** unlocks rounding of variable with respect to the lock status of the constraint and its negation;
 *  this method should be called whenever the lock status of a variable in a constraint changes, for example if
 *  the coefficient of the variable changed its sign or if the left or right hand sides of the constraint were
 *  added or removed
 */
extern
SCIP_RETCODE SCIPunlockVarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             lockdown,           /**< should the rounding be locked in downwards direction? */
   SCIP_Bool             lockup              /**< should the rounding be locked in upwards direction? */
   );

/** changes variable's objective value */
extern
SCIP_RETCODE SCIPchgVarObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             newobj              /**< new objective value */
   );

/** adds value to variable's objective value */
extern
SCIP_RETCODE SCIPaddVarObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             addobj              /**< additional objective value */
   );

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) lower bound value;
 *  does not change the bounds of the variable
 */
extern
SCIP_Real SCIPadjustedVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to adjust the bound for */
   SCIP_Real             lb                  /**< lower bound value to adjust */
   );

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) upper bound value;
 *  does not change the bounds of the variable
 */
extern
SCIP_Real SCIPadjustedVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to adjust the bound for */
   SCIP_Real             ub                  /**< upper bound value to adjust */
   );

/** depending on SCIP's stage, changes lower bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
SCIP_RETCODE SCIPchgVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** depending on SCIP's stage, changes upper bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
SCIP_RETCODE SCIPchgVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** changes lower bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 */
extern
SCIP_RETCODE SCIPchgVarLbNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to change bound at, or NULL for current node */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** changes upper bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 */
extern
SCIP_RETCODE SCIPchgVarUbNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to change bound at, or NULL for current node */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** changes global lower bound of variable; if possible, adjust bound to integral value; also tightens the local bound,
 *  if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
SCIP_RETCODE SCIPchgVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** changes global upper bound of variable; if possible, adjust bound to integral value; also tightens the local bound,
 *  if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
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
 *  @note lazy bounds are useful for branch-and-price since the the corresponding variable bounds are not part of the LP
 */
extern
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
 *  @note lazy bounds are useful for branch-and-price since the the corresponding variable bounds are not part of the LP
 */
extern
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
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
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
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
SCIP_RETCODE SCIPtightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
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
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
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
 */
extern
SCIP_RETCODE SCIPinferBinvarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to fix */
   SCIP_Bool             fixedval,           /**< value to fix binary variable to */
   SCIP_CONS*            infercons,          /**< constraint that deduced the fixing */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the fixing tightened the local bounds, or NULL */
   );

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
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
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
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
 */
extern
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
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
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
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
extern
SCIP_RETCODE SCIPtightenVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

#ifndef NDEBUG

/** for a multi-aggregated variable, returns the global lower bound computed by adding the global bounds from all aggregation variables
 * this global bound may be tighter than the one given by SCIPvarGetLbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetLbGlobal
 */
extern
SCIP_Real SCIPcomputeVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

/** for a multi-aggregated variable, returns the global upper bound computed by adding the global bounds from all aggregation variables
 * this global bound may be tighter than the one given by SCIPvarGetUbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetUbGlobal
 */
extern
SCIP_Real SCIPcomputeVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

/** for a multi-aggregated variable, returns the local lower bound computed by adding the local bounds from all aggregation variables
 * this local bound may be tighter than the one given by SCIPvarGetLbLocal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetLbLocal
 */
extern
SCIP_Real SCIPcomputeVarLbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

/** for a multi-aggregated variable, returns the local upper bound computed by adding the local bounds from all aggregation variables
 * this local bound may be tighter than the one given by SCIPvarGetUbLocal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetUbLocal
 */
extern
SCIP_Real SCIPcomputeVarUbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   );

#else

#define SCIPcomputeVarLbGlobal(scip, var)  (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrLbGlobal(var, (scip)->set) : SCIPvarGetLbGlobal(var))
#define SCIPcomputeVarUbGlobal(scip, var)  (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrUbGlobal(var, (scip)->set) : SCIPvarGetUbGlobal(var))
#define SCIPcomputeVarLbLocal(scip, var)   (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrLbLocal(var, (scip)->set)  : SCIPvarGetLbLocal(var))
#define SCIPcomputeVarUbLocal(scip, var)   (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR ? SCIPvarGetMultaggrUbLocal(var, (scip)->set)  : SCIPvarGetUbLocal(var))

#endif

/** returns solution value and index of variable lower bound that is closest to the variable's value in the given primal solution
 *  or current LP solution if no primal solution is given; returns an index of -1 if no variable upper bound is available
 */
extern
SCIP_RETCODE SCIPgetVarClosestVlb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_Real*            closestvlb,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvlbidx       /**< pointer to store the index of the closest variable lower bound */
   );

/** returns solution value and index of variable upper bound that is closest to the variable's value in the given primal solution;
 *  or current LP solution if no primal solution is given; returns an index of -1 if no variable upper bound is available
 */
extern
SCIP_RETCODE SCIPgetVarClosestVub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_Real*            closestvub,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvubidx       /**< pointer to store the index of the closest variable lower bound */
   );

/** informs variable x about a globally valid variable lower bound x >= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  improves the global bounds of the variable and the vlb variable if possible
 */
extern
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
 *  improves the global bounds of the variable and the vlb variable if possible
 */
extern
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
 */
extern
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
 */
extern
SCIP_RETCODE SCIPaddClique(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   SCIP_Bool*            values,             /**< values of the variables in the clique; NULL to use TRUE for all vars */
   int                   nvars,              /**< number of variables in the clique */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   );

/** calculates a partition of the given set of binary variables into cliques;
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same clique;
 *  the first variable is always assigned to clique 0, and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  for each clique at most 1 variables can be set to TRUE in a feasible solution;
 */
extern
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
 */
extern
SCIP_RETCODE SCIPcalcNegatedCliquePartition(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int const             nvars,              /**< number of variables in the clique */
   int*const             cliquepartition,    /**< array of length nvars to store the clique partition */
   int*const             ncliques            /**< pointer to store the number of cliques actually contained in the partition */
   );

/** gets the number of cliques in the clique table */
extern
int SCIPgetNCliques(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the array of cliques in the clique table */
extern
SCIP_CLIQUE** SCIPgetCliques(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
extern
SCIP_RETCODE SCIPchgVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             branchfactor        /**< factor to weigh variable's branching score with */
   );

/** scales the branch factor of the variable with the given value */
extern
SCIP_RETCODE SCIPscaleVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             scale               /**< factor to scale variable's branching factor with */
   );

/** adds the given value to the branch factor of the variable */
extern
SCIP_RETCODE SCIPaddVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             addfactor           /**< value to add to the branch factor of the variable */
   );

/** sets the branch priority of the variable; variables with higher branch priority are always preferred to variables
 *  with lower priority in selection of branching variable
 *
 * @note the default branching priority is 0
 */
extern
SCIP_RETCODE SCIPchgVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< branch priority of the variable */
   );

/** changes the branch priority of the variable to the given value, if it is larger than the current priority */
extern
SCIP_RETCODE SCIPupdateVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< new branch priority of the variable, if it is larger than current priority */
   );

/** adds the given value to the branch priority of the variable */
extern
SCIP_RETCODE SCIPaddVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   addpriority         /**< value to add to the branch priority of the variable */
   );

/** sets the branch direction of the variable (-1: prefer downwards branch, 0: automatic selection, +1: prefer upwards
 *  branch)
 */
extern
SCIP_RETCODE SCIPchgVarBranchDirection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        branchdirection     /**< preferred branch direction of the variable (downwards, upwards, auto) */
   );

/** changes type of variable in the problem;
 *
 *  @warning This type changes might change the variable array returned from SCIPgetVars() and SCIPgetVarsData();
 *
 *  @note If SCIP is already beyond the SCIP_STAGE_PROBLEM and a original variable is passed; the variable type of the
 *        corresponding transformed variable is changed; the type of the original variable does not change
 *
 *  @note If the type changes from a continuous variable to a non-continuous variable the bound of the variable get
 *        adjusted w.r.t. to integrality information
 */
extern
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
 */
extern
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
 *  In the second step, the variable to be aggregated is chosen among x' and y', preferring a less strict variable
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
 */
extern
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
 */
extern
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
extern
SCIP_Bool SCIPdoNotAggr(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether variable is not allowed to be multi-aggregated */
extern
SCIP_Bool SCIPdoNotMultaggrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable x to aggregate */
   );

/** marks the variable to not to be multi-aggregated */
extern
SCIP_RETCODE SCIPmarkDoNotMultaggrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to delete */
   );

/** enables the collection of statistics for a variable */
extern
void SCIPenableVarHistory(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** diables the collection of any statistic for a variable */
extern
void SCIPdisableVarHistory(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of "solvaldelta" in the
 *  variable's solution value and resulting change of "objdelta" in the in the LP's objective value;
 *  the update is ignored, if the objective value difference is infinite
 */
extern
SCIP_RETCODE SCIPupdateVarPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   SCIP_Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   SCIP_Real             weight              /**< weight in (0,1] of this update in pseudo cost sum */
   );

/** gets the variable's pseudo cost value for the given change of the variable's LP value */
extern
SCIP_Real SCIPgetVarPseudocostVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the variable's pseudo cost value for the given change of the variable's LP value,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPgetVarPseudocostValCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the variable's pseudo cost value for the given direction */
extern
SCIP_Real SCIPgetVarPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the variable's pseudo cost value for the given direction,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPgetVarPseudocostCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
extern
SCIP_Real SCIPgetVarPseudocostCount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPgetVarPseudocostCountCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the variable's pseudo cost score value for the given LP solution value */
extern
SCIP_Real SCIPgetVarPseudocostScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solval              /**< variable's LP solution value */
   );

/** gets the variable's pseudo cost score value for the given LP solution value,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPgetVarPseudocostScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solval              /**< variable's LP solution value */
   );

/** returns the variable's conflict score value */
extern
SCIP_Real SCIPgetVarVSIDS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's conflict score value only using conflicts of the current run */
extern
SCIP_Real SCIPgetVarVSIDSCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's conflict score value */
extern
SCIP_Real SCIPgetVarConflictScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's conflict score value only using conflicts of the current run */
extern
SCIP_Real SCIPgetVarConflictScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's conflict length score */
extern
SCIP_Real SCIPgetVarConflictlengthScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's conflict length score only using conflicts of the current run */
extern
SCIP_Real SCIPgetVarConflictlengthScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's average conflict length */
extern
SCIP_Real SCIPgetVarAvgConflictlength(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's average conflict length only using conflicts of the current run */
extern
SCIP_Real SCIPgetVarAvgConflictlengthCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of inferences found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 */
extern
SCIP_Real SCIPgetVarAvgInferences(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of inferences found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 */
extern
SCIP_Real SCIPgetVarAvgInferencesCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's average inference score value */
extern
SCIP_Real SCIPgetVarAvgInferenceScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's average inference score value only using inferences of the current run */
extern
SCIP_Real SCIPgetVarAvgInferenceScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** initializes the upwards and downwards pseudocosts, conflict scores, conflict lengths, inference scores, cutoff scores
 *  of a variable to the given values 
 */
extern
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

/** returns the average number of cutoffs found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 */
extern
SCIP_Real SCIPgetVarAvgCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of cutoffs found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 */
extern
SCIP_Real SCIPgetVarAvgCutoffsCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's average cutoff score value */
extern
SCIP_Real SCIPgetVarAvgCutoffScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's average cutoff score value, only using cutoffs of the current run */
extern
SCIP_Real SCIPgetVarAvgCutoffScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor
 */
extern
SCIP_Real SCIPgetVarAvgInferenceCutoffScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   );

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor, only using inferences and cutoffs of the current run
 */
extern
SCIP_Real SCIPgetVarAvgInferenceCutoffScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   );

/** outputs variable information to file stream via the message system
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
extern
SCIP_RETCODE SCIPprintVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/**@} */




/*
 * conflict analysis methods
 */

/**@name Conflict Analysis Methods */
/**@{ */

/** return TRUE if conflict analysis is applicable; In case the function return FALSE there is no need to initialize the
 *  conflict analysis since it will not be applied
 */
extern
SCIP_Bool SCIPisConflictAnalysisApplicable(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** initializes the conflict analysis by clearing the conflict candidate queue; this method must be called before you
 *  enter the conflict variables by calling SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar();
 */
extern
SCIP_RETCODE SCIPinitConflictAnalysis(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds lower bound of variable at the time of the given bound change index to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictLb() should be called for each lower bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictLb() should be called
 *      for each lower bound, whose current assignment lead to the deduction of the given conflict bound.
 */
extern
SCIP_RETCODE SCIPaddConflictLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose lower bound should be added to conflict candidate queue */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   );

/** adds lower bound of variable at the time of the given bound change index to the conflict analysis' candidate storage
 *  with the additional information of a relaxed lower bound; this relaxed lower bound is the one which would be enough
 *  to explain a certain bound change;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictRelaxedLb() should be called for each (relaxed)
 *      lower bound that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictRelexedLb() should be
 *      called for each (relaxed) lower bound, whose current assignment lead to the deduction of the given conflict
 *      bound.
 */
extern
SCIP_RETCODE SCIPaddConflictRelaxedLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose lower bound should be added to conflict candidate queue */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   SCIP_Real             relaxedlb           /**< the relaxed lower bound */
   );

/** adds upper bound of variable at the time of the given bound change index to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictUb() should be called for each upper bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictUb() should be called
 *      for each upper bound, whose current assignment lead to the deduction of the given conflict bound.
 */
extern
SCIP_RETCODE SCIPaddConflictUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   );

/** adds upper bound of variable at the time of the given bound change index to the conflict analysis' candidate storage
 *  with the additional information of a relaxed upper bound; this relaxed upper bound is the one which would be enough
 *  to explain a certain bound change;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictRelaxedUb() should be called for each (relaxed)
 *      upper bound that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictRelaxedUb() should be
 *      called for each (relaxed) upper bound, whose current assignment lead to the deduction of the given conflict
 *      bound.
 */
extern
SCIP_RETCODE SCIPaddConflictRelaxedUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   SCIP_Real             relaxedub           /**< the relaxed upper bound */
   );

/** adds lower or upper bound of variable at the time of the given bound change index to the conflict analysis' candidate
 *  storage; this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictBd() should be called for each bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictBd() should be called
 *      for each bound, whose current assignment lead to the deduction of the given conflict bound.
 */
extern
SCIP_RETCODE SCIPaddConflictBd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the conflicting bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   );

/** adds lower or upper bound of variable at the time of the given bound change index to the conflict analysis'
 *  candidate storage; with the additional information of a relaxed upper bound; this relaxed upper bound is the one
 *  which would be enough to explain a certain bound change;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictRelaxedBd() should be called for each (relaxed)
 *      bound that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictRelaxedBd() should be
 *      called for each (relaxed) bound, whose current assignment lead to the deduction of the given conflict bound.
 */
extern
SCIP_RETCODE SCIPaddConflictRelaxedBd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the conflicting bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   SCIP_Real             relaxedbd           /**< the relaxed bound */
   );

/** adds changed bound of fixed binary variable to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictBinvar() should be called for each fixed binary
 *      variable that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictBinvar() should be called
 *      for each binary variable, whose current fixing lead to the deduction of the given conflict bound.
 */
extern
SCIP_RETCODE SCIPaddConflictBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< binary variable whose changed bound should be added to conflict queue */
   );

/** checks if the given variable is already part of the current conflict set or queued for resolving with the same or
 *  even stronger bound
 */
extern
SCIP_RETCODE SCIPisConflictVarUsed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the conflicting bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   SCIP_Bool*            used                /**< pointer to store if the variable is already used */
   );

/** returns the conflict lower bound if the variable is present in the current conflict set; otherwise SCIP_INFINITY */
extern
SCIP_Real SCIPgetConflictVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the conflict upper bound if the variable is present in the current conflict set; otherwise minus
 *  SCIP_INFINITY
 */
extern
SCIP_Real SCIPgetConflictVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the relaxed conflict lower bound if the variable is present in the current conflict set; otherwise
 *  SCIP_INFINITY
 */
extern
SCIP_Real SCIPgetConflictVarRelaxedLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the relaxed conflict upper bound if the variable is present in the current conflict set; otherwise
 *  minus SCIP_INFINITY
 */
extern
SCIP_Real SCIPgetConflictVarRelaxedUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** analyzes conflict bounds that were added after a call to SCIPinitConflictAnalysis() with calls to
 *  SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(),
 *  SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar(); on success, calls the conflict
 *  handlers to create a conflict constraint out of the resulting conflict set; the given valid depth must be a depth
 *  level, at which the conflict set defined by calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), and SCIPaddConflictBinvar() is
 *  valid for the whole subtree; if the conflict was found by a violated constraint, use SCIPanalyzeConflictCons()
 *  instead of SCIPanalyzeConflict() to make sure, that the correct valid depth is used
 */
extern
SCIP_RETCODE SCIPanalyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   );

/** analyzes conflict bounds that were added with calls to SCIPaddConflictLb(), SCIPaddConflictUb(),
 *  SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or
 *  SCIPaddConflictBinvar(); on success, calls the conflict handlers to create a conflict constraint out of the
 *  resulting conflict set; the given constraint must be the constraint that detected the conflict, i.e. the constraint
 *  that is infeasible in the local bounds of the initial conflict set (defined by calls to SCIPaddConflictLb(),
 *  SCIPaddConflictUb(), SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(),
 *  SCIPaddConflictRelaxedBd(), and SCIPaddConflictBinvar())
 */
extern
SCIP_RETCODE SCIPanalyzeConflictCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that detected the conflict */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   );

/**@} */




/*
 * constraint methods
 */

/**@name Constraint Methods */
/**@{ */

/** creates and captures a constraint of the given constraint handler
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution may
 *           be declared feasible even if it violates this particular constraint.  This constellation should only be
 *           used, if no LP or pseudo solution can violate the constraint -- e.g. if a local constraint is redundant due
 *           to the variable's local bounds.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
extern
SCIP_RETCODE SCIPcreateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   SCIP_CONSDATA*        consdata,           /**< data for this specific constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** parses constraint information (in cip format) out of a string; if the parsing process was successful a constraint is
 *  creates and captures;
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution may
 *           be declared feasible even if it violates this particular constraint.  This constellation should only be
 *           used, if no LP or pseudo solution can violate the constraint -- e.g. if a local constraint is redundant due
 *           to the variable's local bounds.
 */
extern
SCIP_RETCODE SCIPparseCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to constraint */
   const char*           str,                /**< name of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   );

/** increases usage counter of constraint */
extern
SCIP_RETCODE SCIPcaptureCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to capture */
   );

/** decreases usage counter of constraint, if the usage pointer reaches zero the constraint gets freed
 *
 *  @note the pointer of the constraint will be NULLed
 */
extern
SCIP_RETCODE SCIPreleaseCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons                /**< pointer to constraint */
   );

/** sets the initial flag of the given constraint */
extern
SCIP_RETCODE SCIPsetConsInitial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             initial             /**< new value */
   );

/** sets the separate flag of the given constraint */
extern
SCIP_RETCODE SCIPsetConsSeparated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             separate            /**< new value */
   );

/** sets the enforce flag of the given constraint */
extern
SCIP_RETCODE SCIPsetConsEnforced(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             enforce             /**< new value */
   );

/** sets the check flag of the given constraint */
extern
SCIP_RETCODE SCIPsetConsChecked(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             check               /**< new value */
   );

/** sets the propagate flag of the given constraint */
extern
SCIP_RETCODE SCIPsetConsPropagated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             propagate           /**< new value */
   );

/** sets the local flag of the given constraint */
extern
SCIP_RETCODE SCIPsetConsLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             local               /**< new value */
   );

/** sets the modifiable flag of the given constraint */
extern
SCIP_RETCODE SCIPsetConsModifiable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             modifiable          /**< new value */
   );

/** sets the dynamic flag of the given constraint */
extern
SCIP_RETCODE SCIPsetConsDynamic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             dynamic             /**< new value */
   );

/** sets the removable flag of the given constraint */
extern
SCIP_RETCODE SCIPsetConsRemovable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             removable           /**< new value */
   );

/** sets the stickingatnode flag of the given constraint */
extern
SCIP_RETCODE SCIPsetConsStickingAtNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             stickingatnode      /**< new value */
   );

/** gets and captures transformed constraint of a given constraint; if the constraint is not yet transformed,
 *  a new transformed constraint for this constraint is created
 */
extern
SCIP_RETCODE SCIPtransformCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get/create transformed constraint for */
   SCIP_CONS**           transcons           /**< pointer to store the transformed constraint */
   );

/** gets and captures transformed constraints for an array of constraints;
 *  if a constraint in the array is not yet transformed, a new transformed constraint for this constraint is created;
 *  it is possible to call this method with conss == transconss
 */
extern
SCIP_RETCODE SCIPtransformConss(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nconss,             /**< number of constraints to get/create transformed constraints for */
   SCIP_CONS**           conss,              /**< array with constraints to get/create transformed constraints for */
   SCIP_CONS**           transconss          /**< array to store the transformed constraints */
   );

/** gets corresponding transformed constraint of a given constraint;
 *  returns NULL as transcons, if transformed constraint is not yet existing
 */
extern
SCIP_RETCODE SCIPgetTransformedCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get the transformed constraint for */
   SCIP_CONS**           transcons           /**< pointer to store the transformed constraint */
   );

/** gets corresponding transformed constraints for an array of constraints;
 *  stores NULL in a transconss slot, if the transformed constraint is not yet existing;
 *  it is possible to call this method with conss == transconss, but remember that constraints that are not
 *  yet transformed will be replaced with NULL
 */
extern
SCIP_RETCODE SCIPgetTransformedConss(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nconss,             /**< number of constraints to get the transformed constraints for */
   SCIP_CONS**           conss,              /**< constraints to get the transformed constraints for */
   SCIP_CONS**           transconss          /**< array to store the transformed constraints */
   );

/** adds given value to age of constraint, but age can never become negative;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 */
extern
SCIP_RETCODE SCIPaddConsAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             deltaage            /**< value to add to the constraint's age */
   );

/** increases age of constraint by 1.0;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 */
extern
SCIP_RETCODE SCIPincConsAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** resets age of constraint to zero;
 *  should be called
 *   - in constraint separation, if a cut was found for this constraint,
 *   - in constraint enforcing, if the constraint was violated, and
 *   - in constraint propagation, if a domain reduction was deduced;
 */
extern
SCIP_RETCODE SCIPresetConsAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** enables constraint's separation, propagation, and enforcing capabilities */
extern
SCIP_RETCODE SCIPenableCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** disables constraint's separation, propagation, and enforcing capabilities, s.t. the constraint is not propagated,
 *  separated, and enforced anymore until it is enabled again with a call to SCIPenableCons();
 *  in contrast to SCIPdelConsLocal() and SCIPdelConsNode(), the disabling is not associated to a node in the tree and
 *  does not consume memory; therefore, the constraint is neither automatically enabled on leaving the node nor
 *  automatically disabled again on entering the node again;
 *  note that the constraints enforcing capabilities are necessary for the solution's feasibility, if the constraint
 *  is a model constraint; that means, you must be sure that the constraint cannot be violated in the current subtree,
 *  and you have to enable it again manually by calling SCIPenableCons(), if this subtree is left (e.g. by using
 *  an appropriate event handler that watches the corresponding variables' domain changes)
 */
extern
SCIP_RETCODE SCIPdisableCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** enables constraint's separation capabilities */
extern
SCIP_RETCODE SCIPenableConsSeparation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** disables constraint's separation capabilities s.t. the constraint is not propagated anymore until the separation
 *  is enabled again with a call to SCIPenableConsSeparation(); in contrast to SCIPdelConsLocal() and SCIPdelConsNode(),
 *  the disabling is not associated to a node in the tree and does not consume memory; therefore, the constraint
 *  is neither automatically enabled on leaving the node nor automatically disabled again on entering the node again
 */
extern
SCIP_RETCODE SCIPdisableConsSeparation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** enables constraint's propagation capabilities */
extern
SCIP_RETCODE SCIPenableConsPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** disables constraint's propagation capabilities s.t. the constraint is not propagated anymore until the propagation
 *  is enabled again with a call to SCIPenableConsPropagation(); in contrast to SCIPdelConsLocal() and SCIPdelConsNode(),
 *  the disabling is not associated to a node in the tree and does not consume memory; therefore, the constraint
 *  is neither automatically enabled on leaving the node nor automatically disabled again on entering the node again
 */
extern
SCIP_RETCODE SCIPdisableConsPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** adds given values to lock status of the constraint and updates the rounding locks of the involved variables */
extern
SCIP_RETCODE SCIPaddConsLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nlockspos,          /**< increase in number of rounding locks for constraint */
   int                   nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   );

/** checks single constraint for feasibility of the given solution */
extern
SCIP_RETCODE SCIPcheckCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows (both local and global) to be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** enforces single constraint for a given pseudo solution */
extern
SCIP_RETCODE SCIPenfopsCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to enforce */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_Bool             objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** enforces single constraint for a given LP solution */
extern
SCIP_RETCODE SCIPenfolpCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to enforce */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls LP initialization method for single constraint */
extern
SCIP_RETCODE SCIPinitlpCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to initialize */
   );

/** calls separation method of single constraint for LP solution */
extern
SCIP_RETCODE SCIPsepalpCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to separate */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   );

/** calls separation method of single constraint for given primal solution */
extern
SCIP_RETCODE SCIPsepasolCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to separate */
   SCIP_SOL*             sol,                /**< primal solution that should be separated*/
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   );

/** calls domain propagation method of single constraint */
extern
SCIP_RETCODE SCIPpropCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to propagate */
   SCIP_PROPTIMING       proptiming,         /**< current point in the node solving loop */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** resolves propagation conflict of single constraint */
extern
SCIP_RETCODE SCIPrespropCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to resolve conflict for */
   SCIP_VAR*             infervar,           /**< the conflict variable whose bound change has to be resolved */
   int                   inferinfo,          /**< the user information passed to the corresponding SCIPinferVarLbCons() or SCIPinferVarUbCons() call */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** presolves of single constraint */
extern
SCIP_RETCODE SCIPpresolCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to presolve */
   int                   nrounds,            /**< number of presolving rounds already done */
   int                   nnewfixedvars,      /**< number of variables fixed since the last call to the presolving method */
   int                   nnewaggrvars,       /**< number of variables aggregated since the last call to the presolving method */
   int                   nnewchgvartypes,    /**< number of variable type changes since the last call to the presolving method */
   int                   nnewchgbds,         /**< number of variable bounds tightened since the last call to the presolving method */
   int                   nnewholes,          /**< number of domain holes added since the last call to the presolving method */
   int                   nnewdelconss,       /**< number of deleted constraints since the last call to the presolving method */
   int                   nnewaddconss,       /**< number of added constraints since the last call to the presolving method */
   int                   nnewupgdconss,      /**< number of upgraded constraints since the last call to the presolving method */
   int                   nnewchgcoefs,       /**< number of changed coefficients since the last call to the presolving method */
   int                   nnewchgsides,       /**< number of changed left or right hand sides since the last call to the presolving method */
   int*                  nfixedvars,         /**< pointer to count total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to count total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to count total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to count total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to count total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to count total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to count total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to count total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to count total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls constraint activation notification method of single constraint */
extern
SCIP_RETCODE SCIPactiveCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to notify */
   );

/** calls constraint deactivation notification method of single constraint */
extern
SCIP_RETCODE SCIPdeactiveCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to notify */
   );

/** outputs constraint information to file stream via the message handler system
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
extern
SCIP_RETCODE SCIPprintCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** method to collect the variables of a constraint
 *
 *  If the number of variables is greater than the available slots in the variable array, nothing happens except that
 *  the success point is set to FALSE. With the method SCIPgetConsNVars() it is possible to get the number of variables
 *  a constraint has in its scope.
 *
 *  @note The success pointer indicates if all variables were copied into the vars arrray.
 *
 *  @note It might be that a constraint handler does not support this functionality, in that case the success pointer is
 *        set to FALSE.
 */
extern
SCIP_RETCODE SCIPgetConsVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which the variables are wanted */
   SCIP_VAR**            vars,               /**< array to store the involved variable of the constraint */
   int                   varssize,           /**< available slots in vars array which is needed to check if the array is large enough */
   SCIP_Bool*            success             /**< pointer to store whether the variables are successfully copied */
   );

/** method to collect the number of variables of a constraint
 *
 *  @note The success pointer indicates if the contraint handler was able to return the number of variables
 *
 *  @note It might be that a constraint handler does not support this functionality, in that case the success pointer is
 *        set to FALSE
 */
extern
SCIP_RETCODE SCIPgetConsNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which the number of variables is wanted */
   int*                  nvars,              /**< pointer to store the number of variables */
   SCIP_Bool*            success             /**< pointer to store whether the constraint successfully returned the number of variables */
   );

/**@} */




/*
 * LP methods
 */

/**@name LP Methods */
/**@{ */

/** returns, whether the LP was or is to be solved in the current node */
extern
SCIP_Bool SCIPhasCurrentNodeLP(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns, whether the LP of the current node is already constructed */
extern
SCIP_Bool SCIPisLPConstructed(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** makes sure that the LP of the current node is loaded and may be accessed through the LP information methods */
extern
SCIP_RETCODE SCIPconstructLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   );

/** makes sure that the LP of the current node is flushed */
extern
SCIP_RETCODE SCIPflushLP(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets solution status of current LP */
extern
SCIP_LPSOLSTAT SCIPgetLPSolstat(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the current lp is a relaxation of the current problem and its optimal objective value is a local lower bound */
extern
SCIP_Bool SCIPisLPRelax(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets objective value of current LP (which is the sum of column and loose objective value) */
extern
SCIP_Real SCIPgetLPObjval(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets part of objective value of current LP that results from COLUMN variables only */
extern
SCIP_Real SCIPgetLPColumnObjval(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets part of objective value of current LP that results from LOOSE variables only */
extern
SCIP_Real SCIPgetLPLooseObjval(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the global pseudo objective value; that is all variables set to their best  (w.r.t. the objective
 *  function) global bound
 */
extern
SCIP_Real SCIPgetGlobalPseudoObjval(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the pseudo objective value for the current search node; that is all variables set to their best (w.r.t. the
 *  objective function) local bound
 */
extern
SCIP_Real SCIPgetPseudoObjval(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the root lp is a relaxation of the problem and its optimal objective value is a global lower bound */
extern
SCIP_Bool SCIPisRootLPRelax(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the objective value of the root node LP; returns SCIP_INVALID if the root node LP was not (yet) solved */
extern
SCIP_Real SCIPgetLPRootObjval(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets part of the objective value of the root node LP that results from COLUMN variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 */
extern
SCIP_Real SCIPgetLPRootColumnObjval(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets part of the objective value of the root node LP that results from LOOSE variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 */
extern
SCIP_Real SCIPgetLPRootLooseObjval(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current LP columns along with the current number of LP columns */
extern
SCIP_RETCODE SCIPgetLPColsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL***           cols,               /**< pointer to store the array of LP columns, or NULL */
   int*                  ncols               /**< pointer to store the number of LP columns, or NULL */
   );

/** gets current LP columns */
extern
SCIP_COL** SCIPgetLPCols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current number of LP columns */
extern
int SCIPgetNLPCols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current LP rows along with the current number of LP rows */
extern
SCIP_RETCODE SCIPgetLPRowsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW***           rows,               /**< pointer to store the array of LP rows, or NULL */
   int*                  nrows               /**< pointer to store the number of LP rows, or NULL */
   );

/** gets current LP rows */
extern
SCIP_ROW** SCIPgetLPRows(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current number of LP rows */
extern
int SCIPgetNLPRows(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing
 */
extern
SCIP_Bool SCIPallColsInLP(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the current LP solution is basic, i.e. is defined by a valid simplex basis */
extern
SCIP_Bool SCIPisLPSolBasic(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1 */
extern
SCIP_RETCODE SCIPgetLPBasisInd(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  basisind            /**< pointer to store the basis indices */
   );

/** gets a row from the inverse basis matrix B^-1 */
extern
SCIP_RETCODE SCIPgetLPBInvRow(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   );

/** gets a column from the inverse basis matrix B^-1 */
extern
SCIP_RETCODE SCIPgetLPBInvCol(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP
                                              *   returned by SCIPcolGetLPPos(); you have to call SCIPgetBasisInd()
                                              *   to get the array which links the B^-1 column numbers to the row and
                                              *   column numbers of the LP! c must be between 0 and nrows-1, since the
                                              *   basis has the size nrows * nrows */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   );

/** gets a row from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A) */
extern
SCIP_RETCODE SCIPgetLPBInvARow(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   r,                  /**< row number */
   SCIP_Real*            binvrow,            /**< row in B^-1 from prior call to SCIPgetLPBInvRow(), or NULL */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   );

/** gets a column from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A),
 *  i.e., it computes B^-1 * A_c with A_c being the c'th column of A
 */
extern
SCIP_RETCODE SCIPgetLPBInvACol(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   c,                  /**< column number which can be accessed by SCIPcolGetLPPos() */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   );

/** calculates a weighted sum of all LP rows; for negative weights, the left and right hand side of the corresponding
 *  LP row are swapped in the summation
 */
extern
SCIP_RETCODE SCIPsumLPRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_REALARRAY*       sumcoef,            /**< array to store sum coefficients indexed by variables' probindex */
   SCIP_Real*            sumlhs,             /**< pointer to store the left hand side of the row summation */
   SCIP_Real*            sumrhs              /**< pointer to store the right hand side of the row summation */
   );

/* calculates a MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
extern
SCIP_RETCODE SCIPcalcMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,  
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables; 
                                              *   NULL for using closest bound for all variables */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mksetcoefs,         /**< array to store mixed knapsack set coefficients: size nvars; or NULL */
   SCIP_Bool*            mksetcoefsvalid,    /**< pointer to store whether mixed knapsack set coefficients are valid; or NULL */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   SCIP_Bool*            cutislocal          /**< pointer to store whether the returned cut is only valid locally */
   );


/* calculates a strong CG cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
extern
SCIP_RETCODE SCIPcalcStrongCG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce strong CG cut for */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mircoef,            /**< array to store strong CG coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the strong CG row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid strong CG cut */
   SCIP_Bool*            cutislocal          /**< pointer to store whether the returned cut is only valid locally */
   );

/** writes current LP to a file */
extern
SCIP_RETCODE SCIPwriteLP(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   );

/** writes MIP relaxation of the current branch-and-bound node to a file */
extern
SCIP_RETCODE SCIPwriteMIP(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< file name */
   SCIP_Bool             genericnames,       /**< should generic names like x_i and row_j be used in order to avoid
                                              *   troubles with reserved strings? */
   SCIP_Bool             origobj             /**< should the original objective function be used? */
   );

/** gets the LP interface of SCIP;
 *  with the LPI you can use all of the methods defined in scip/lpi.h;
 *
 *  @warning You have to make sure, that the full internal state of the LPI does not change or is recovered completely
 *           after the end of the method that uses the LPI. In particular, if you manipulate the LP or its solution
 *           (e.g. by calling one of the SCIPlpiAdd...() or one of the SCIPlpiSolve...() methods), you have to check in
 *           advance with SCIPlpiWasSolved() whether the LP is currently solved. If this is the case, you have to make
 *           sure, the internal solution status is recovered completely at the end of your method. This can be achieved
 *           by getting the LPI state before applying any LPI manipulations with SCIPlpiGetState() and restoring it
 *           afterwards with SCIPlpiSetState() and SCIPlpiFreeState(). Additionally you have to resolve the LP with the
 *           appropriate SCIPlpiSolve...() call in order to reinstall the internal solution status.
 *
 *  @warning Make also sure, that all parameter values that you have changed are set back to their original values.
 */
extern
SCIP_RETCODE SCIPgetLPI(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LPI**            lpi                 /**< pointer to store the LP interface */
   );

/** Displays quality information about the current LP solution. An LP solution need to be available. Information printed
 *  is subject to what the LP solver supports
 *
 *  @note The printing process is done via the message handler system.
 */
extern
SCIP_RETCODE SCIPprintLPSolutionQuality(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** compute relative interior point to current LP
 * @see SCIPlpComputeRelIntPoint
 */
extern
SCIP_RETCODE SCIPcomputeLPRelIntPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             relaxrows,          /**< should the rows be relaxed */
   SCIP_Bool             inclobjcutoff,      /**< should a row for the objective cutoff be included */
   char                  normtype,           /**< which norm to use: 'o'ne-norm or 's'upremum-norm */
   SCIP_SOL**            point               /**< relative interior point on exit */
   );

/**@} */



/*
 * LP column methods
 */

/**@name LP Column Methods */
/**@{ */

/** returns the reduced costs of a column in the last (feasible) LP */
extern
SCIP_Real SCIPgetColRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL*             col                 /**< LP column */
   );

/** returns the Farkas coefficient of a column in the last (infeasible) LP */
extern
SCIP_Real SCIPgetColFarkasCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL*             col                 /**< LP column */
   );

/**@} */




/*
 * LP row methods
 */

/**@name LP Row Methods */
/**@{ */

/** creates and captures an LP row */
extern
SCIP_RETCODE SCIPcreateRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   const char*           name,               /**< name of row */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COL**            cols,               /**< array with columns of row entries */
   SCIP_Real*            vals,               /**< array with coefficients of row entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   );

/** creates and captures an LP row without any coefficients */
extern
SCIP_RETCODE SCIPcreateEmptyRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   const char*           name,               /**< name of row */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   );

/** increases usage counter of LP row */
extern
SCIP_RETCODE SCIPcaptureRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row to capture */
   );

/** decreases usage counter of LP row, and frees memory if necessary */
extern
SCIP_RETCODE SCIPreleaseRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row                 /**< pointer to LP row */
   );

/** changes left hand side of LP row */
extern
SCIP_RETCODE SCIPchgRowLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             lhs                 /**< new left hand side */
   );

/** changes right hand side of LP row */
extern
SCIP_RETCODE SCIPchgRowRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             rhs                 /**< new right hand side */
   );

/** informs row, that all subsequent additions of variables to the row should be cached and not directly applied;
 *  after all additions were applied, SCIPflushRowExtensions() must be called;
 *  while the caching of row extensions is activated, information methods of the row give invalid results;
 *  caching should be used, if a row is build with SCIPaddVarToRow() calls variable by variable to increase
 *  the performance
 */
extern
SCIP_RETCODE SCIPcacheRowExtensions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** flushes all cached row extensions after a call of SCIPcacheRowExtensions() and merges coefficients with
 *  equal columns into a single coefficient
 */
extern
SCIP_RETCODE SCIPflushRowExtensions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** resolves variable to columns and adds them with the coefficient to the row */
extern
SCIP_RETCODE SCIPaddVarToRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             val                 /**< value of coefficient */
   );

/** resolves variables to columns and adds them with the coefficients to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 */
extern
SCIP_RETCODE SCIPaddVarsToRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Real*            vals                /**< values of coefficients */
   );

/** resolves variables to columns and adds them with the same single coefficient to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 */
extern
SCIP_RETCODE SCIPaddVarsToRowSameCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Real             val                 /**< unique value of all coefficients */
   );

/** tries to find a value, such that all row coefficients, if scaled with this value become integral */
extern
SCIP_RETCODE SCIPcalcRowIntegralScalar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   );

/** tries to scale row, s.t. all coefficients (of integer variables) become integral */
extern
SCIP_RETCODE SCIPmakeRowIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal value to scale row with */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Bool*            success             /**< stores whether row could be made rational */
   );

/** returns minimal absolute value of row vector's non-zero coefficients */
extern
SCIP_Real SCIPgetRowMinCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns maximal absolute value of row vector's non-zero coefficients */
extern
SCIP_Real SCIPgetRowMaxCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the minimal activity of a row w.r.t. the column's bounds */
extern
SCIP_Real SCIPgetRowMinActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the maximal activity of a row w.r.t. the column's bounds */
extern
SCIP_Real SCIPgetRowMaxActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** recalculates the activity of a row in the last LP solution */
extern
SCIP_RETCODE SCIPrecalcRowLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the activity of a row in the last LP solution */
extern
SCIP_Real SCIPgetRowLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the feasibility of a row in the last LP solution: negative value means infeasibility */
extern
SCIP_Real SCIPgetRowLPFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** recalculates the activity of a row for the current pseudo solution */
extern
SCIP_RETCODE SCIPrecalcRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the activity of a row for the current pseudo solution */
extern
SCIP_Real SCIPgetRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the feasibility of a row for the current pseudo solution: negative value means infeasibility */
extern
SCIP_Real SCIPgetRowPseudoFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** recalculates the activity of a row in the last LP or pseudo solution */
extern
SCIP_RETCODE SCIPrecalcRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the activity of a row in the last LP or pseudo solution */
extern
SCIP_Real SCIPgetRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the feasibility of a row in the last LP or pseudo solution */
extern
SCIP_Real SCIPgetRowFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the activity of a row for the given primal solution */
extern
SCIP_Real SCIPgetRowSolActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** returns the feasibility of a row for the given primal solution */
extern
SCIP_Real SCIPgetRowSolFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** output row to file stream via the message handler system */
extern
SCIP_RETCODE SCIPprintRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/**@} */


/*
 * NLP methods
 */

/**@name NLP Methods */
/**@{ */

/** marks that there are constraints that are representable by nonlinear constraints involving continuous variables
 * This method should be called by a constraint handler if it has constraints that have a representation as nonlinear rows
 * and that are still nonlinear after fixing all discrete variables in the CIP.
 * 
 * The function should be called before the branch-and-bound process is initialized, e.g., when presolve is exiting.
 * 
 */ 
extern
void SCIPmarkContinuousNonlinearitiesPresent(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** marks that there are constraints that are representable by nonlinear constraints
 * This method should be called by a constraint handler if it has constraints that have a representation as nonlinear rows.
 * 
 * The function should be called before the branch-and-bound process is initialized, e.g., when presolve is exiting.
 * 
 * Calling SCIPmarkContinuousNonlinearitiesPresent makes a call to SCIPmarkNonlinearitiesPresent dispensable.
 */ 
extern
void SCIPmarkNonlinearitiesPresent(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether constraints representable as nonlinear rows are present that involve continuous nonlinear variables
 * @see SCIPmarkContinuousNonlinearitiesPresent
 */
extern
SCIP_Bool SCIPhasContinuousNonlinearitiesPresent(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether constraints representable as nonlinear rows are present
 * @see SCIPmarkNonlinearitiesPresent
 */
extern
SCIP_Bool SCIPhasNonlinearitiesPresent(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns, whether an NLP has been constructed */
extern
SCIP_Bool SCIPisNLPConstructed(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current NLP variables along with the current number of NLP variables */
extern
SCIP_RETCODE SCIPgetNLPVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to store the array of NLP variables, or NULL */
   int*                  nvars               /**< pointer to store the number of NLP variables, or NULL */
   );

/** gets array with variables of the NLP */
extern
SCIP_VAR** SCIPgetNLPVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current number of variables in NLP */
extern
int SCIPgetNNLPVars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** computes for each variables the number of NLP rows in which the variable appears in a nonlinear var */
extern
SCIP_RETCODE SCIPgetNLPVarsNonlinearity(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  nlcount             /**< an array of length at least SCIPnlpGetNVars() to store nonlinearity counts of variables */
   );

/** gives dual solution values associated with lower bounds of NLP variables */
extern
SCIP_Real* SCIPgetNLPVarsLbDualsol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gives dual solution values associated with upper bounds of NLP variables */
extern
SCIP_Real* SCIPgetNLPVarsUbDualsol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current NLP nonlinear rows along with the current number of NLP nonlinear rows */
extern
SCIP_RETCODE SCIPgetNLPNlRowsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW***         nlrows,             /**< pointer to store the array of NLP nonlinear rows, or NULL */
   int*                  nnlrows             /**< pointer to store the number of NLP nonlinear rows, or NULL */
   );

/** gets array with nonlinear rows of the NLP */
extern
SCIP_NLROW** SCIPgetNLPNlRows(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current number of nonlinear rows in NLP */
extern
int SCIPgetNNLPNlRows(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds a nonlinear row to the NLP
 * row is captured by NLP */
extern
SCIP_RETCODE SCIPaddNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< nonlinear row to add to NLP */
   );

/** makes sure that the NLP of the current node is flushed */
extern
SCIP_RETCODE SCIPflushNLP(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets or clears initial primal guess for NLP solution (start point for NLP solver) */
extern
SCIP_RETCODE SCIPsetNLPInitialGuess(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            initialguess        /**< values of initial guess (corresponding to variables from SCIPgetNLPVarsData), or NULL to use no start point */
   );

/** sets initial primal guess for NLP solution (start point for NLP solver) */
extern
SCIP_RETCODE SCIPsetNLPInitialGuessSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< solution which values should be taken as initial guess, or NULL for LP solution */
   );

/** solves the current NLP */
extern
SCIP_RETCODE SCIPsolveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets solution status of current NLP */
extern
SCIP_NLPSOLSTAT SCIPgetNLPSolstat(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets termination status of last NLP solve */
extern
SCIP_NLPTERMSTAT SCIPgetNLPTermstat(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gives statistics (number of iterations, solving time, ...) of last NLP solve */
extern
SCIP_RETCODE SCIPgetNLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   );

/** gets objective value of current NLP */
extern
SCIP_Real SCIPgetNLPObjval(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** indicates whether a feasible solution for the current NLP is available
 * thus, returns whether the solution status <= feasible  */
extern
SCIP_Bool SCIPhasNLPSolution(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets fractional variables of last NLP solution along with solution values and fractionalities */
extern
SCIP_RETCODE SCIPgetNLPFracVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           fracvars,           /**< pointer to store the array of NLP fractional variables, or NULL */
   SCIP_Real**           fracvarssol,        /**< pointer to store the array of NLP fractional variables solution values, or NULL */
   SCIP_Real**           fracvarsfrac,       /**< pointer to store the array of NLP fractional variables fractionalities, or NULL */
   int*                  nfracvars,          /**< pointer to store the number of NLP fractional variables , or NULL */
   int*                  npriofracvars       /**< pointer to store the number of NLP fractional variables with maximal branching priority, or NULL */
   );

/** gets integer parameter of NLP */
extern
SCIP_RETCODE SCIPgetNLPIntPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   );

/** sets integer parameter of NLP */
extern
SCIP_RETCODE SCIPsetNLPIntPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   );

/** gets floating point parameter of NLP */
extern
SCIP_RETCODE SCIPgetNLPRealPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   );

/** sets floating point parameter of NLP */
extern
SCIP_RETCODE SCIPsetNLPRealPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   );

/** gets string parameter of NLP */
extern
SCIP_RETCODE SCIPgetNLPStringPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char**          sval                /**< pointer to store the parameter value */
   );

/** sets string parameter of NLP */
extern
SCIP_RETCODE SCIPsetNLPStringPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char*           sval                /**< parameter value */
   );

/** writes current NLP to a file */
extern
SCIP_RETCODE SCIPwriteNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   );

/** gets the NLP interface and problem used by the SCIP NLP;
 *  with the NLPI and its problem you can use all of the methods defined in nlpi/nlpi.h;
 *
 *  @warning You have to make sure, that the full internal state of the NLPI does not change or is recovered completely
 *           after the end of the method that uses the NLPI. In particular, if you manipulate the NLP or its solution
 *           (e.g. by calling one of the SCIPnlpiAdd...() or the SCIPnlpiSolve() method), you have to check in advance
 *           whether the NLP is currently solved.  If this is the case, you have to make sure, the internal solution
 *           status is recovered completely at the end of your method. Additionally you have to resolve the NLP with
 *           SCIPnlpiSolve() in order to reinstall the internal solution status.
 *
 *  @warning Make also sure, that all parameter values that you have changed are set back to their original values.
 */
extern
SCIP_RETCODE SCIPgetNLPI(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI**           nlpi,               /**< pointer to store the NLP solver interface */
   SCIP_NLPIPROBLEM**    nlpiproblem         /**< pointer to store the NLP solver interface problem */
   );

/**@} */


/*
 * NLP diving methods
 */

/**@name NLP Diving Methods */
/**@{ */

/** initiates NLP diving
 * making methods SCIPchgVarObjDiveNLP(), SCIPchgVarBoundsDiveNLP(), SCIPchgVarsBoundsDiveNLP(), and SCIPsolveDiveNLP() available */
extern
SCIP_RETCODE SCIPstartDiveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** ends NLP diving
 * resets changes made by SCIPchgVarObjDiveNLP(), SCIPchgVarBoundsDiveNLP(), and SCIPchgVarsBoundsDiveNLP() */
extern
SCIP_RETCODE SCIPendDiveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** changes cutoffbound in current dive */
SCIP_RETCODE SCIPchgCutoffboundDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newcutoffbound      /**< new cutoffbound */
   );

/** changes linear objective coefficient of a variable in diving NLP */
extern
SCIP_RETCODE SCIPchgVarObjDiveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             coef                /**< new value for coefficient */
   );

/** changes bounds of a variable in diving NLP */
extern
SCIP_RETCODE SCIPchgVarBoundsDiveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which bounds to change */
   SCIP_Real             lb,                 /**< new lower bound */
   SCIP_Real             ub                  /**< new upper bound */
   );

/** changes bounds of a set of variables in diving NLP */
extern
SCIP_RETCODE SCIPchgVarsBoundsDiveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables which bounds to changes */
   SCIP_VAR**            vars,               /**< variables which bounds to change */
   SCIP_Real*            lbs,                /**< new lower bounds */
   SCIP_Real*            ubs                 /**< new upper bounds */
   );

/** solves diving NLP */
extern
SCIP_RETCODE SCIPsolveDiveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */


/*
 * NLP nonlinear row methods
 */

/**@name NLP Nonlinear Row Methods */
/**@{ */

/** creates and captures an NLP row */
extern
SCIP_RETCODE SCIPcreateNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   const char*           name,               /**< name of nonlinear row */
   SCIP_Real             constant,           /**< constant */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< linear variables, or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< linear coefficients, or NULL if nlinvars == 0 */
   int                   nquadvars,          /**< number of variables in quadratic term */
   SCIP_VAR**            quadvars,           /**< variables in quadratic terms, or NULL if nquadvars == 0 */
   int                   nquadelems,         /**< number of elements in quadratic term */
   SCIP_QUADELEM*        quadelems,          /**< elements (i.e., monomials) in quadratic term, or NULL if nquadelems == 0 */
   SCIP_EXPRTREE*        expression,         /**< nonlinear expression, or NULL */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   );

/** creates and captures an NLP nonlinear row without any coefficients */
extern
SCIP_RETCODE SCIPcreateEmptyNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   const char*           name,               /**< name of nonlinear row */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   );

/** creates and captures an NLP row from a linear row */
extern
SCIP_RETCODE SCIPcreateNlRowFromRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   SCIP_ROW*             row                 /**< the linear row to copy */
   );

/** increases usage counter of NLP nonlinear row */
extern
SCIP_RETCODE SCIPcaptureNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< nonlinear row to capture */
   );

/** decreases usage counter of NLP nonlinear row, and frees memory if necessary */
extern
SCIP_RETCODE SCIPreleaseNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow               /**< nonlinear row to release */
   );

/** changes left hand side of NLP nonlinear row */
extern
SCIP_RETCODE SCIPchgNlRowLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real             lhs                 /**< new left hand side */
   );

/** changes right hand side of NLP nonlinear row */
extern
SCIP_RETCODE SCIPchgNlRowRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real             rhs                 /**< new right hand side */
   );

/** changes constant of NLP nonlinear row */
extern
SCIP_RETCODE SCIPchgNlRowConstant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_Real             constant            /**< new value for constant */
   );

/** adds variable with a linear coefficient to the nonlinear row */
extern
SCIP_RETCODE SCIPaddLinearCoefToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             val                 /**< value of coefficient in linear part of row */
   );

/** adds variables with linear coefficients to the row */
extern
SCIP_RETCODE SCIPaddLinearCoefsToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Real*            vals                /**< values of coefficients in linear part of row */
   );

/** changes linear coefficient of a variables in a row
 * setting the coefficient to 0.0 means that it is removed from the row
 * the variable does not need to exists before */
extern
SCIP_RETCODE SCIPchgNlRowLinearCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< new value of coefficient */
   );

/** adds quadratic variable to the nonlinear row
 * after adding a quadratic variable, it can be used to add quadratic elements */
extern
SCIP_RETCODE SCIPaddQuadVarToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** adds quadratic variables to the nonlinear row
 * after adding quadratic variables, they can be used to add quadratic elements */
extern
SCIP_RETCODE SCIPaddQuadVarsToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   nvars,              /**< number of problem variables */
   SCIP_VAR**            vars                /**< problem variables */
   );

/** add a quadratic element to the nonlinear row
 * variable indices of the quadratic element need to be relative to quadratic variables array of row */
extern
SCIP_RETCODE SCIPaddQuadElementToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_QUADELEM         quadelem            /**< quadratic element */
   );

/** adds quadratic elements to the nonlinear row
 * variable indices of the quadratic elements need to be relative to quadratic variables array of row */
extern
SCIP_RETCODE SCIPaddQuadElementsToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   nquadelems,         /**< number of quadratic elements */
   SCIP_QUADELEM*        quadelems           /**< quadratic elements */
   );

/** changes coefficient in quadratic part of a row
 * setting the coefficient in the quadelement to 0.0 means that it is removed from the row
 * the element does not need to exists before */
extern
SCIP_RETCODE SCIPchgNlRowQuadElement(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_QUADELEM         quadelement         /**< new quadratic element, or update for existing one */
   );

/** sets or deletes expression tree in the nonlinear row */
extern
SCIP_RETCODE SCIPsetNlRowExprtree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_EXPRTREE*        exprtree            /**< expression tree, or NULL */
   );

/** sets a parameter of expression tree in the nonlinear row */
extern
SCIP_RETCODE SCIPsetNlRowExprtreeParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   paramidx,           /**< index of parameter in expression tree */
   SCIP_Real             paramval            /**< new value of parameter in expression tree */
   );

/** sets parameters of expression tree in the nonlinear row */
extern
SCIP_RETCODE SCIPsetNlRowExprtreeParams(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_Real*            paramvals           /**< new values of parameter in expression tree */
   );

/** recalculates the activity of a nonlinear row in the last NLP solution */
extern
SCIP_RETCODE SCIPrecalcNlRowNLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< NLP nonlinear row */
   );

/** returns the activity of a nonlinear row in the last NLP solution */
extern
SCIP_RETCODE SCIPgetNlRowNLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            activity            /**< buffer to store activity value */
   );

/** gives the feasibility of a nonlinear row in the last NLP solution: negative value means infeasibility */
extern
SCIP_RETCODE SCIPgetNlRowNLPFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   );

/** recalculates the activity of a nonlinear row for the current pseudo solution */
extern
SCIP_RETCODE SCIPrecalcNlRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< NLP nonlinear row */
   );

/** gives the activity of a nonlinear row for the current pseudo solution */
extern
SCIP_RETCODE SCIPgetNlRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            pseudoactivity      /**< buffer to store pseudo activity value */
   );

/** gives the feasibility of a nonlinear row for the current pseudo solution: negative value means infeasibility */
extern
SCIP_RETCODE SCIPgetNlRowPseudoFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            pseudofeasibility   /**< buffer to store pseudo feasibility value */
   );

/** recalculates the activity of a nonlinear row in the last NLP or pseudo solution */
extern
SCIP_RETCODE SCIPrecalcNlRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< NLP nonlinear row */
   );

/** gives the activity of a nonlinear row in the last NLP or pseudo solution */
extern
SCIP_RETCODE SCIPgetNlRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            activity            /**< buffer to store activity value */
   );

/** gives the feasibility of a nonlinear row in the last NLP or pseudo solution */
extern
SCIP_RETCODE SCIPgetNlRowFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   );

/** gives the activity of a nonlinear row for the given primal solution or NLP solution or pseudo solution */
extern
SCIP_RETCODE SCIPgetNlRowSolActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for NLP solution of pseudo solution */
   SCIP_Real*            activity            /**< buffer to store activity value */
   );

/** gives the feasibility of a nonlinear row for the given primal solution */
extern
SCIP_RETCODE SCIPgetNlRowSolFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   );

/** gives the minimal and maximal activity of a nonlinear row w.r.t. the variable's bounds */
extern
SCIP_RETCODE SCIPgetNlRowActivityBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_Real*            minactivity,        /**< buffer to store minimal activity, or NULL */
   SCIP_Real*            maxactivity         /**< buffer to store maximal activity, or NULL */
   );

/** output nonlinear row to file stream */
extern
SCIP_RETCODE SCIPprintNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/**@} */

/**@name Expression tree methods */
/**@{ */

/** replaced array of variables in expression tree by corresponding transformed variables */
extern
SCIP_RETCODE SCIPgetExprtreeTransformedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** evaluates an expression tree for a primal solution or LP solution */
extern
SCIP_RETCODE SCIPevalExprtreeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_SOL*             sol,                /**< a solution, or NULL for current LP solution */
   SCIP_Real*            val                 /**< buffer to store value */
   );

/** evaluates an expression tree w.r.t. current global bounds */
extern
SCIP_RETCODE SCIPevalExprtreeGlobalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        val                 /**< buffer to store result */
   );

/** evaluates an expression tree w.r.t. current local bounds */
extern
SCIP_RETCODE SCIPevalExprtreeLocalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        val                 /**< buffer to store result */
   );

/**@} */


/*
 * cutting plane methods
 */

/**@name Cutting Plane Methods */
/**@{ */

/** returns efficacy of the cut with respect to the given primal solution or the current LP solution:
 *  e = -feasibility/norm
 */
extern
SCIP_Real SCIPgetCutEfficacy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for current LP solution */
   SCIP_ROW*             cut                 /**< separated cut */
   );

/** returns whether the cut's efficacy with respect to the given primal solution or the current LP solution is greater
 *  than the minimal cut efficacy
 */
extern
SCIP_Bool SCIPisCutEfficacious(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for current LP solution */
   SCIP_ROW*             cut                 /**< separated cut */
   );

/** checks, if the given cut's efficacy is larger than the minimal cut efficacy */
extern
SCIP_Bool SCIPisEfficacious(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             efficacy            /**< efficacy of the cut */
   );

/** calculates the efficacy norm of the given vector, which depends on the "separating/efficacynorm" parameter */
extern
SCIP_Real SCIPgetVectorEfficacyNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            vals,               /**< array of values */
   int                   nvals               /**< number of values */
   );

/** adds cut to separation storage */
extern
SCIP_RETCODE SCIPaddCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that was separated, or NULL for LP solution */
   SCIP_ROW*             cut,                /**< separated cut */
   SCIP_Bool             forcecut            /**< should the cut be forced to enter the LP? */
   );

/** if not already existing, adds row to global cut pool */
extern
SCIP_RETCODE SCIPaddPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** removes the row from the global cut pool */
extern
SCIP_RETCODE SCIPdelPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** gets current cuts in the global cut pool */
extern
SCIP_CUT** SCIPgetPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current number of rows in the global cut pool */
extern
int SCIPgetNPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the global cut pool used by SCIP */
extern
SCIP_CUTPOOL* SCIPgetGlobalCutpool(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a cut pool */
extern
SCIP_RETCODE SCIPcreateCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   int                   agelimit            /**< maximum age a cut can reach before it is deleted from the pool */
   );

/** frees a cut pool */
extern
SCIP_RETCODE SCIPfreeCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL**        cutpool             /**< pointer to store cut pool */
   );

/** if not already existing, adds row to a cut pool and captures it */
extern
SCIP_RETCODE SCIPaddRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** adds row to a cut pool and captures it; doesn't check for multiple cuts */
extern
SCIP_RETCODE SCIPaddNewRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** removes the LP row from a cut pool */
extern
SCIP_RETCODE SCIPdelRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< row to remove */
   );

/** separates cuts from a cut pool */
extern
SCIP_RETCODE SCIPseparateCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   );

/** separates the given primal solution or the current LP solution by calling the separators and constraint handlers'
 *  separation methods;
 *  the generated cuts are stored in the separation storage and can be accessed with the methods SCIPgetCuts() and
 *  SCIPgetNCuts();
 *  after evaluating the cuts, you have to call SCIPclearCuts() in order to remove the cuts from the
 *  separation storage;
 *  it is possible to call SCIPseparateSol() multiple times with different solutions and evaluate the found cuts
 *  afterwards
 */
extern
SCIP_RETCODE SCIPseparateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_Bool             pretendroot,        /**< should the cut separators be called as if we are at the root node? */
   SCIP_Bool             onlydelayed,        /**< should only separators be called that were delayed in the previous round? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a separator was delayed */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   );

/** gets the array of cuts currently stored in the separation storage */
extern
SCIP_ROW** SCIPgetCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get current number of cuts in the separation storage */
extern
int SCIPgetNCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** clears the separation storage */
extern
SCIP_RETCODE SCIPclearCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** removes cuts that are inefficacious w.r.t. the current LP solution from separation storage without adding the cuts to the LP */
extern
SCIP_RETCODE SCIPremoveInefficaciousCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */




/*
 * LP diving methods
 */

/**@name LP Diving Methods */
/**@{ */

/** initiates LP diving, making methods SCIPchgVarObjDive(), SCIPchgVarLbDive(), and SCIPchgVarUbDive() available
 *
 *  @note diving is allowed even if the current LP is not flushed, not solved, or not solved to optimality; be aware
 *  that solving the (first) diving LP may take longer than expect and that the latter two cases could stem from
 *  numerical troubles during the last LP solve; because of this, most users will want to call this method only if
 *  SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL
 */
extern
SCIP_RETCODE SCIPstartDive(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** quits LP diving and resets bounds and objective values of columns to the current node's values */
extern
SCIP_RETCODE SCIPendDive(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** changes variable's objective value in current dive */
extern
SCIP_RETCODE SCIPchgVarObjDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             newobj              /**< new objective value */
   );

/** changes variable's lower bound in current dive */
extern
SCIP_RETCODE SCIPchgVarLbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** changes variable's upper bound in current dive */
extern
SCIP_RETCODE SCIPchgVarUbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** gets variable's objective value in current dive */
extern
SCIP_Real SCIPgetVarObjDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get the bound for */
   );

/** gets variable's lower bound in current dive */
extern
SCIP_Real SCIPgetVarLbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get the bound for */
   );

/** gets variable's upper bound in current dive */
extern
SCIP_Real SCIPgetVarUbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get the bound for */
   );

/** solves the LP of the current dive; no separation or pricing is applied
 *
 *  @note be aware that the LP solve may take longer than expected if SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL,
 *  compare the explanation of SCIPstartDive()
 */
extern
SCIP_RETCODE SCIPsolveDiveLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   );

/** returns the number of the node in the current branch and bound run, where the last LP was solved in diving
 *  or probing mode
 */
extern
SCIP_Longint SCIPgetLastDivenode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether we are in diving mode */
extern
SCIP_Bool SCIPinDive(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */




/*
 * probing methods
 */

/**@name Probing Methods */
/**@{ */

/** returns whether we are in probing mode */
extern
SCIP_Bool SCIPinProbing(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** initiates probing, making methods SCIPnewProbingNode(), SCIPbacktrackProbing(), SCIPchgVarLbProbing(), 
 *  SCIPchgVarUbProbing(), SCIPfixVarProbing(), SCIPpropagateProbing(), and SCIPsolveProbingLP() available
 */
extern
SCIP_RETCODE SCIPstartProbing(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a new probing sub node, whose changes can be undone by backtracking to a higher node in the probing path
 *  with a call to SCIPbacktrackProbing();
 *  using a sub node for each set of probing bound changes can improve conflict analysis
 */
extern
SCIP_RETCODE SCIPnewProbingNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the current probing depth, i.e. the number of probing sub nodes existing in the probing path */
extern
int SCIPgetProbingDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** undoes all changes to the problem applied in probing up to the given probing depth;
 *  the changes of the probing node of the given probing depth are the last ones that remain active;
 *  changes that were applied before calling SCIPnewProbingNode() cannot be undone
 */
extern
SCIP_RETCODE SCIPbacktrackProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated */
   );

/** quits probing and resets bounds and constraints to the focus node's environment */
extern
SCIP_RETCODE SCIPendProbing(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** injects a change of variable's lower bound into current probing node; the same can also be achieved with a call to
 *  SCIPchgVarLb(), but in this case, the bound change would be treated like a deduction instead of a branching decision
 */
extern
SCIP_RETCODE SCIPchgVarLbProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** injects a change of variable's upper bound into current probing node; the same can also be achieved with a call to
 *  SCIPchgVarUb(), but in this case, the bound change would be treated like a deduction instead of a branching decision
 */
extern
SCIP_RETCODE SCIPchgVarUbProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** injects a change of variable's bounds into current probing node to fix the variable to the specified value;
 *  the same can also be achieved with a call to SCIPfixVar(), but in this case, the bound changes would be treated
 *  like deductions instead of branching decisions
 */
extern
SCIP_RETCODE SCIPfixVarProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             fixedval            /**< value to fix variable to */
   );

/** applies domain propagation on the probing sub problem, that was changed after SCIPstartProbing() was called;
 *  the propagated domains of the variables can be accessed with the usual bound accessing calls SCIPvarGetLbLocal()
 *  and SCIPvarGetUbLocal(); the propagation is only valid locally, i.e. the local bounds as well as the changed
 *  bounds due to SCIPchgVarLbProbing(), SCIPchgVarUbProbing(), and SCIPfixVarProbing() are used for propagation
 */
extern
SCIP_RETCODE SCIPpropagateProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the probing node can be cut off */
   SCIP_Longint*         ndomredsfound       /**< pointer to store the number of domain reductions found, or NULL */
   );

/** applies domain propagation on the probing sub problem, that was changed after SCIPstartProbing() was called;
 *  only propagations of the binary variables fixed at the current probing node that are triggered by the implication
 *  graph and the clique table are applied;
 *  the propagated domains of the variables can be accessed with the usual bound accessing calls SCIPvarGetLbLocal()
 *  and SCIPvarGetUbLocal(); the propagation is only valid locally, i.e. the local bounds as well as the changed
 *  bounds due to SCIPchgVarLbProbing(), SCIPchgVarUbProbing(), and SCIPfixVarProbing() are used for propagation
 */
extern
SCIP_RETCODE SCIPpropagateProbingImplications(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing node can be cut off */
   );

/** solves the LP at the current probing node (cannot be applied at preprocessing stage);
 *  no separation or pricing is applied
 */
extern
SCIP_RETCODE SCIPsolveProbingLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   );

/** solves the LP at the current probing node (cannot be applied at preprocessing stage) and applies pricing
 *  until the LP is solved to optimality; no separation is applied
 */
extern
SCIP_RETCODE SCIPsolveProbingLPWithPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             pretendroot,        /**< should the pricers be called as if we are at the root node? */
   SCIP_Bool             displayinfo,        /**< should info lines be displayed after each pricing round? */
   int                   maxpricerounds,     /**< maximal number of pricing rounds (-1: no limit);
                                              *   a finite limit means that the LP might not be solved to optimality! */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   );

/**@} */




/*
 * branching methods
 */

/**@name Branching Methods */
/**@{ */

/** gets branching candidates for LP solution branching (fractional variables) along with solution values,
 *  fractionalities, and number of branching candidates;
 *  branching rules should always select the branching candidate among the first npriolpcands of the candidate
 *  list
 */
extern
SCIP_RETCODE SCIPgetLPBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   SCIP_Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   SCIP_Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*                  nlpcands,           /**< pointer to store the number of LP branching candidates, or NULL */
   int*                  npriolpcands        /**< pointer to store the number of candidates with maximal priority, or NULL */
   );

/** gets number of branching candidates for LP solution branching (number of fractional variables) */
extern
int SCIPgetNLPBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of branching candidates with maximal priority for LP solution branching */
extern
int SCIPgetNPrioLPBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets external branching candidates along with solution values, scores, and number of branching candidates;
 *  these branching candidates can be used by relaxations or nonlinear constraint handlers
 *  branching rules should always select the branching candidate among the first nprioexterncands of the candidate
 *  list
 */
extern
SCIP_RETCODE SCIPgetExternBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           externcands,        /**< pointer to store the array of extern branching candidates, or NULL */
   SCIP_Real**           externcandssol,     /**< pointer to store the array of extern candidate solution values, or NULL */
   SCIP_Real**           externcandsscore,   /**< pointer to store the array of extern candidate scores, or NULL */
   int*                  nexterncands,       /**< pointer to store the number of extern branching candidates, or NULL */
   int*                  nprioexterncands,   /**< pointer to store the number of candidates with maximal priority, or NULL */
   int*                  nprioexternbins,    /**< pointer to store the number of binary candidates with maximal priority, or NULL */
   int*                  nprioexternints,    /**< pointer to store the number of integer candidates with maximal priority, or NULL */
   int*                  nprioexternimpls    /**< pointer to store the number of implicit integer candidates with maximal priority, 
                                              *   or NULL */
   );

/** gets number of external branching candidates */
extern
int SCIPgetNExternBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of external branching candidates with maximal branch priority */
extern
int SCIPgetNPrioExternBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of binary external branching candidates with maximal branch priority */
extern
int SCIPgetNPrioExternBranchBins(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of integer external branching candidates with maximal branch priority */
extern
int SCIPgetNPrioExternBranchInts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of implicit integer external branching candidates with maximal branch priority */
extern
int SCIPgetNPrioExternBranchImpls(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of continuous external branching candidates with maximal branch priority */
extern
int SCIPgetNPrioExternBranchConts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** insert variable, its score and its solution value into the external branching candidate storage
 * the relative difference of the current lower and upper bounds of a continuous variable must be at least epsilon
 */
extern
SCIP_RETCODE SCIPaddExternBranchCand(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to insert */
   SCIP_Real             score,              /**< score of external candidate, e.g. infeasibility */
   SCIP_Real             solval              /**< value of the variable in the current solution */
   );

/** removes all external candidates from the storage for external branching */
extern
void SCIPclearExternBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** checks whether the given variable is contained in the candidate storage for external branching */
extern
SCIP_Bool SCIPcontainsExternBranchCand(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to look for */
   );

/** gets branching candidates for pseudo solution branching (non-fixed variables) along with the number of candidates */
extern
SCIP_RETCODE SCIPgetPseudoBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*                  npseudocands,       /**< pointer to store the number of pseudo branching candidates, or NULL */
   int*                  npriopseudocands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   );

/** gets branching candidates for pseudo solution branching (non-fixed variables) */
extern
int SCIPgetNPseudoBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPgetNPrioPseudoBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of binary branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPgetNPrioPseudoBranchBins(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of integer branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPgetNPrioPseudoBranchInts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of implicit integer branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPgetNPrioPseudoBranchImpls(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** calculates the branching score out of the gain predictions for a binary branching */
extern
SCIP_Real SCIPgetBranchScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   SCIP_Real             downgain,           /**< prediction of objective gain for rounding downwards */
   SCIP_Real             upgain              /**< prediction of objective gain for rounding upwards */
   );

/** calculates the branching score out of the gain predictions for a branching with arbitrary many children */
extern
SCIP_Real SCIPgetBranchScoreMultiple(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   int                   nchildren,          /**< number of children that the branching will create */
   SCIP_Real*            gains               /**< prediction of objective gain for each child */
   );

/** computes a branching point for a continuous or discrete variable
 * @see SCIPbranchGetBranchingPoint
 */
extern
SCIP_Real SCIPgetBranchingPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable, of which the branching point should be computed */
   SCIP_Real             suggestion          /**< suggestion for branching point, or SCIP_INVALID if no suggestion */
   );

/** calculates the node selection priority for moving the given variable's LP value to the given target value;
 *  this node selection priority can be given to the SCIPcreateChild() call
 */
extern
SCIP_Real SCIPcalcNodeselPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable on which the branching is applied */
   SCIP_BRANCHDIR        branchdir,          /**< type of branching that was performed: upwards, downwards, or fixed;
                                              *   fixed should only be used, when both bounds changed
                                              */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   );

/** calculates an estimate for the objective of the best feasible solution contained in the subtree after applying the given 
 *  branching; this estimate can be given to the SCIPcreateChild() call
 */
extern
SCIP_Real SCIPcalcChildEstimate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable on which the branching is applied */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   );

/** creates a child node of the focus node */
extern
SCIP_RETCODE SCIPcreateChild(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE**           node,               /**< pointer to node data structure */
   SCIP_Real             nodeselprio,        /**< node selection priority of new node */
   SCIP_Real             estimate            /**< estimate for (transformed) objective value of best feasible solution in subtree */
   );

/** branches on a non-continuous variable v using the current LP or pseudo solution;
 *  if solution value x' is fractional, two child nodes will be created
 *  (x <= floor(x'), x >= ceil(x')),
 *  if solution value is integral, the x' is equal to lower or upper bound of the branching
 *  variable and the bounds of v are finite, then two child nodes will be created
 *  (x <= x", x >= x"+1 with x" = floor((lb + ub)/2)),
 *  otherwise (up to) three child nodes will be created
 *  (x <= x'-1, x == x', x >= x'+1)
 */
extern
SCIP_RETCODE SCIPbranchVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           eqchild,            /**< pointer to return the middle child with variable fixed, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   );

/** branches on a variable x using a given value x'; 
 *  for continuous variables with relative domain width larger epsilon, x' must not be one of the bounds;
 *  two child nodes (x <= x', x >= x') are created;
 *  for integer variables, if solution value x' is fractional, two child nodes are created
 *  (x <= floor(x'), x >= ceil(x')),
 *  if x' is integral, three child nodes are created
 *  (x <= x'-1, x == x', x >= x'+1)
 */
extern
SCIP_RETCODE SCIPbranchVarVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             val,                /**< value to branch on */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           eqchild,            /**< pointer to return the middle child with variable fixed, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   );

/** n-ary branching on a variable x using a given value
 * Branches on variable x such that up to n/2 children are created on each side of the usual branching value.
 * The branching value is selected as in SCIPbranchVarVal().
 * The parameters minwidth and widthfactor determine the domain width of the branching variable in the child nodes.
 * If n is odd, one child with domain width 'width' and having the branching value in the middle is created.
 * Otherwise, two children with domain width 'width' and being left and right of the branching value are created.
 * Next further nodes to the left and right are created, where width is multiplied by widthfactor with increasing distance from the first nodes.
 * The initial width is calculated such that n/2 nodes are created to the left and to the right of the branching value.
 * If this value is below minwidth, the initial width is set to minwidth, which may result in creating less than n nodes.
 *
 * Giving a large value for widthfactor results in creating children with small domain when close to the branching value
 * and large domain when closer to the current variable bounds. That is, setting widthfactor to a very large value and n to 3
 * results in a ternary branching where the branching variable is mostly fixed in the middle child.
 * Setting widthfactor to 1.0 results in children where the branching variable always has the same domain width
 * (except for one child if the branching value is not in the middle).
 */
extern
SCIP_RETCODE SCIPbranchVarValNary(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             val,                /**< value to branch on */
   int                   n,                  /**< attempted number of children to be created, must be >= 2 */
   SCIP_Real             minwidth,           /**< minimal domain width in children */
   SCIP_Real             widthfactor,        /**< multiplier for children domain width with increasing distance from val, must be >= 1.0 */
   int*                  nchildren           /**< buffer to store number of created children, or NULL */
   );

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 */
extern
SCIP_RETCODE SCIPbranchLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   );

/** calls branching rules to branch on a external candidates; if no such candidates exist, the result is SCIP_DIDNOTRUN */
extern
SCIP_RETCODE SCIPbranchExtern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   );

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN */
extern
SCIP_RETCODE SCIPbranchPseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   );

/**@} */




/*
 * primal solutions
 */

/**@name Primal Solution Methods */
/**@{ */

/** creates a primal solution, initialized to zero */
extern
SCIP_RETCODE SCIPcreateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current LP solution */
extern
SCIP_RETCODE SCIPcreateLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current NLP solution */
extern
SCIP_RETCODE SCIPcreateNLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current relaxation solution */
extern
SCIP_RETCODE SCIPcreateRelaxSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current pseudo solution */
extern
SCIP_RETCODE SCIPcreatePseudoSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current LP or pseudo solution, depending on whether the LP was solved
 *  at the current node
 */
extern
SCIP_RETCODE SCIPcreateCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to unknown values */
extern
SCIP_RETCODE SCIPcreateUnknownSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution living in the original problem space, initialized to zero;
 *  a solution in original space allows to set original variables to values, that would be invalid in the
 *  transformed problem due to preprocessing fixings or aggregations
 */
extern
SCIP_RETCODE SCIPcreateOrigSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a copy of a primal solution; note that a copy of a linked solution is also linked and needs to be unlinked
 *  if it should stay unaffected from changes in the LP or pseudo solution
 */
extern
SCIP_RETCODE SCIPcreateSolCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_SOL*             sourcesol           /**< primal CIP solution to copy */
   );

/** frees primal CIP solution */
extern
SCIP_RETCODE SCIPfreeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol                 /**< pointer to the solution */
   );

/** links a primal solution to the current LP solution */
extern
SCIP_RETCODE SCIPlinkLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the current NLP solution */
extern
SCIP_RETCODE SCIPlinkNLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the current relaxation solution */
extern
SCIP_RETCODE SCIPlinkRelaxSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the current pseudo solution */
extern
SCIP_RETCODE SCIPlinkPseudoSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the current LP or pseudo solution */
extern
SCIP_RETCODE SCIPlinkCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** clears a primal solution */
extern
SCIP_RETCODE SCIPclearSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** stores solution values of variables in solution's own array */
extern
SCIP_RETCODE SCIPunlinkSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** sets value of variable in primal CIP solution */
extern
SCIP_RETCODE SCIPsetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_Real             val                 /**< solution value of variable */
   );

/** sets values of multiple variables in primal CIP solution */
extern
SCIP_RETCODE SCIPsetSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   int                   nvars,              /**< number of variables to set solution value for */
   SCIP_VAR**            vars,               /**< array with variables to add to solution */
   SCIP_Real*            vals                /**< array with solution values of variables */
   );

/** increases value of variable in primal CIP solution */
extern
SCIP_RETCODE SCIPincSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_VAR*             var,                /**< variable to increase solution value for */
   SCIP_Real             incval              /**< increment for solution value of variable */
   );

/** returns value of variable in primal CIP solution, or in current LP/pseudo solution */
extern
SCIP_Real SCIPgetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   SCIP_VAR*             var                 /**< variable to get value for */
   );

/** gets values of multiple variables in primal CIP solution */
extern
SCIP_RETCODE SCIPgetSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables to get solution value for */
   SCIP_VAR**            vars,               /**< array with variables to get value for */
   SCIP_Real*            vals                /**< array to store solution values of variables */
   );

/** returns objective value of primal CIP solution w.r.t. original problem, or current LP/pseudo objective value */
extern
SCIP_Real SCIPgetSolOrigObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   );

/** returns transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value */
extern
SCIP_Real SCIPgetSolTransObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   );

/** maps original space objective value into transformed objective value */
extern
SCIP_Real SCIPtransformObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             obj                 /**< original space objective value to transform */
   );

/** maps transformed objective value into original space */
extern
SCIP_Real SCIPretransformObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             obj                 /**< transformed objective value to retransform in original space */
   );

/** gets clock time, when this solution was found */
extern
SCIP_Real SCIPgetSolTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** gets branch and bound run number, where this solution was found */
extern
int SCIPgetSolRunnum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** gets node number of the specific branch and bound run, where this solution was found */
extern
SCIP_Longint SCIPgetSolNodenum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
extern
SCIP_HEUR* SCIPgetSolHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   );

/** returns whether two given solutions are exactly equal */
extern
SCIP_Bool SCIPareSolsEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol1,               /**< first primal CIP solution */
   SCIP_SOL*             sol2                /**< second primal CIP solution */
   );

/** outputs non-zero variables of solution in original problem space to file stream */
extern
SCIP_RETCODE SCIPprintSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** outputs non-zero variables of solution in transformed problem space to file stream */
extern
SCIP_RETCODE SCIPprintTransSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** outputs non-zero variables of solution representing a ray in original problem space to file stream */
extern
SCIP_RETCODE SCIPprintRay(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution representing ray */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** gets number of feasible primal solutions stored in the solution storage in case the problem is transformed; in case
 *  if the problem stage is SCIP_STAGE_PROBLEM, it returns the number solution in the original solution candidate
 *  storage
 */
extern
int SCIPgetNSols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets array of feasible primal solutions stored in the solution storage in case the problem is transformed; in case
 *  if the problem stage is in SCIP_STAGE_PROBLEM, it returns the number array of solution candidate stored
 */
extern
SCIP_SOL** SCIPgetSols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets best feasible primal solution found so far if the problem is transformed; in case the problem is in
 *  SCIP_STAGE_PROBLEM it returns the best solution candidate, or NULL if no solution has been found or the candidate
 *  store is empty;
 */
extern
SCIP_SOL* SCIPgetBestSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** outputs best feasible primal solution found so far to file stream */
extern
SCIP_RETCODE SCIPprintBestSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** outputs best feasible primal solution found so far in transformed variables to file stream */
extern
SCIP_RETCODE SCIPprintBestTransSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** try to round given solution */
extern
SCIP_RETCODE SCIProundSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Bool*            success             /**< pointer to store whether rounding was successful */
   );

/** retransforms solution to original problem space */
extern
SCIP_RETCODE SCIPretransformSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** reads a given solution file, problem has to be transformed in advance */
extern
SCIP_RETCODE SCIPreadSol(
   SCIP*                 scip,              /**< SCIP data structure */
   const char*           filename           /**< name of the input file */
   );

/** adds feasible primal solution to solution storage by copying it */
extern
SCIP_RETCODE SCIPaddSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds primal solution to solution storage, frees the solution afterwards */
extern
SCIP_RETCODE SCIPaddSolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds current LP/pseudo solution to solution storage */
extern
SCIP_RETCODE SCIPaddCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** checks solution for feasibility; if possible, adds it to storage by copying */
extern
SCIP_RETCODE SCIPtrySol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< should all reasons of violations be printed? */
   SCIP_Bool             checkbounds,        /**< should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows (both local and global) to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards */
extern
SCIP_RETCODE SCIPtrySolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool             printreason,        /**< should all reasons of violations be printed? */
   SCIP_Bool             checkbounds,        /**< should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows (both local and global) to be checked? */
   SCIP_Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   );

/** checks current LP/pseudo solution for feasibility; if possible, adds it to storage */
extern
SCIP_RETCODE SCIPtryCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_Bool             printreason,        /**< should all reasons of violations be printed? */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows (both local and global) to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

/** checks solution for feasibility without adding it to the solution store */
extern
SCIP_RETCODE SCIPcheckSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< should all reasons of violations be printed? */
   SCIP_Bool             checkbounds,        /**< should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows (both local and global) to be checked? */
   SCIP_Bool*            feasible            /**< stores whether given solution is feasible */
   );

/** checks solution for feasibility in original problem without adding it to the solution store;
 *  this method is used to double check a solution in order to validate the presolving process
 */
extern
SCIP_RETCODE SCIPcheckSolOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            feasible,           /**< stores whether given solution is feasible */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool             completely          /**< should all violations be checked? */
   );

/** return whether a primal ray is stored that proves unboundedness of the LP relaxation */
extern
SCIP_Bool SCIPhasPrimalRay(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets value of given variable in primal ray causing unboundedness of the LP relaxation; 
 *  should only be called if such a ray is stored (check with SCIPhasPrimalRay()) */
extern
SCIP_Real SCIPgetPrimalRayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get value for */
   );

/**@} */




/*
 * event methods
 */

/**@name Event Methods
 *
 * Events can only be caught during the operation on the transformed problem.
 * Events on variables can only be caught for transformed variables.
 * If you want to catch an event for an original variable, you have to get the corresponding transformed variable
 * with a call to SCIPgetTransformedVar() and catch the event on the transformed variable.
 */
/**@{ */

/** catches a global (not variable or row dependent) event */
extern
SCIP_RETCODE SCIPcatchEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   );

/** drops a global event (stops to track event) */
extern
SCIP_RETCODE SCIPdropEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int                   filterpos           /**< position of event filter entry returned by SCIPcatchEvent(), or -1 */
   );

/** catches an objective value or domain change event on the given transformed variable */
extern
SCIP_RETCODE SCIPcatchVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< transformed variable to catch event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   );

/** drops an objective value or domain change event (stops to track event) on the given transformed variable */
extern
SCIP_RETCODE SCIPdropVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< transformed variable to drop event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int                   filterpos           /**< position of event filter entry returned by SCIPcatchVarEvent(), or -1 */
   );

/** catches a row coefficient, constant, or side change event on the given row */
extern
SCIP_RETCODE SCIPcatchRowEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< linear row to catch event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   );

/** drops a row coefficient, constant, or side change event (stops to track event) on the given row */
extern
SCIP_RETCODE SCIPdropRowEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< linear row to drop event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int                   filterpos           /**< position of event filter entry returned by SCIPcatchVarEvent(), or -1 */
   );

/**@} */




/*
 * tree methods
 */

/**@name Tree Methods */
/**@{ */

/** gets current node in the tree */
extern
SCIP_NODE* SCIPgetCurrentNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the root node of the tree */
extern
SCIP_NODE* SCIPgetRootNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the current node is already solved and only propagated again */
extern
SCIP_Bool SCIPinRepropagation(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets children of focus node along with the number of children */
extern
SCIP_RETCODE SCIPgetChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          children,           /**< pointer to store children array, or NULL if not needed */
   int*                  nchildren           /**< pointer to store number of children, or NULL if not needed */
   );

/** gets number of children of focus node */
extern
int SCIPgetNChildren(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets siblings of focus node along with the number of siblings */
extern
SCIP_RETCODE SCIPgetSiblings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          siblings,           /**< pointer to store siblings array, or NULL if not needed */
   int*                  nsiblings           /**< pointer to store number of siblings, or NULL if not needed */
   );

/** gets number of siblings of focus node */
extern
int SCIPgetNSiblings(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets leaves of the tree along with the number of leaves */
extern
SCIP_RETCODE SCIPgetLeaves(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          leaves,             /**< pointer to store leaves array, or NULL if not needed */
   int*                  nleaves             /**< pointer to store number of leaves, or NULL if not needed */
   );

/** gets number of leaves in the tree */
extern
int SCIPgetNLeaves(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best child of the focus node w.r.t. the node selection priority assigned by the branching rule */
extern
SCIP_NODE* SCIPgetPrioChild(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule */
extern
SCIP_NODE* SCIPgetPrioSibling(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best child of the focus node w.r.t. the node selection strategy */
extern
SCIP_NODE* SCIPgetBestChild(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best sibling of the focus node w.r.t. the node selection strategy */
extern
SCIP_NODE* SCIPgetBestSibling(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best leaf from the node queue w.r.t. the node selection strategy */
extern
SCIP_NODE* SCIPgetBestLeaf(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy */
extern
SCIP_NODE* SCIPgetBestNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the node with smallest lower bound from the tree (child, sibling, or leaf) */
extern
SCIP_NODE* SCIPgetBestboundNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** cuts off node and whole sub tree from branch and bound tree */
extern
SCIP_RETCODE SCIPcutoffNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node that should be cut off */
   );

/** marks the given node to be propagated again the next time a node of its subtree is processed */
extern
SCIP_RETCODE SCIPrepropagateNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node that should be propagated again */
   );

/** returns depth of first node in active path that is marked being cutoff */
extern
int SCIPgetCutoffdepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns depth of first node in active path that has to be propagated again */
extern
int SCIPgetRepropdepth(
   SCIP*                 scip                /**< SCIP data structure */
   );


/* prints all branching decisions on variables from the root to the given node */
extern
SCIP_RETCODE SCIPprintNodeRootPath(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/**@} */


/*
 * statistic methods
 */

/**@name Statistic Methods */
/**@{ */

/** gets number of branch and bound runs performed, including the current run */
extern
int SCIPgetNRuns(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of processed nodes in current run, including the focus node */
extern
SCIP_Longint SCIPgetNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of processed nodes in all runs, including the focus node */
extern
SCIP_Longint SCIPgetNTotalNodes(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of nodes left in the tree (children + siblings + leaves) */
extern
int SCIPgetNNodesLeft(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far */
extern
SCIP_Longint SCIPgetNLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of iterations used so far in primal and dual simplex and barrier algorithm */
extern
SCIP_Longint SCIPgetNLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of iterations used so far in primal and dual simplex and barrier algorithm for the root node */
extern
SCIP_Longint SCIPgetNRootLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of primal LPs solved so far */
extern
SCIP_Longint SCIPgetNPrimalLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of iterations used so far in primal simplex */
extern
SCIP_Longint SCIPgetNPrimalLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of dual LPs solved so far */
extern
SCIP_Longint SCIPgetNDualLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of iterations used so far in dual simplex */
extern
SCIP_Longint SCIPgetNDualLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of barrier LPs solved so far */
extern
SCIP_Longint SCIPgetNBarrierLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of iterations used so far in barrier algorithm */
extern
SCIP_Longint SCIPgetNBarrierLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far that were resolved from an advanced start basis */
extern
SCIP_Longint SCIPgetNResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in primal and dual simplex calls where an advanced start basis
 *  was available
 */
extern
SCIP_Longint SCIPgetNResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of primal LPs solved so far that were resolved from an advanced start basis */
extern
SCIP_Longint SCIPgetNPrimalResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in primal simplex calls where an advanced start basis
 *  was available
 */
extern
SCIP_Longint SCIPgetNPrimalResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of dual LPs solved so far that were resolved from an advanced start basis */
extern
SCIP_Longint SCIPgetNDualResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in dual simplex calls where an advanced start basis
 *  was available
 */
extern
SCIP_Longint SCIPgetNDualResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far for node relaxations */
extern
SCIP_Longint SCIPgetNNodeLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far for node relaxations */
extern
SCIP_Longint SCIPgetNNodeLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far for initial LP in node relaxations */
extern
SCIP_Longint SCIPgetNNodeInitLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far for initial LP in node relaxations */
extern
SCIP_Longint SCIPgetNNodeInitLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far during diving and probing */
extern
SCIP_Longint SCIPgetNDivingLPs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far during diving and probing */
extern
SCIP_Longint SCIPgetNDivingLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of times, strong branching was called (each call represents solving two LPs) */
extern
SCIP_Longint SCIPgetNStrongbranchs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in strong branching */
extern
SCIP_Longint SCIPgetNStrongbranchLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of times, strong branching was called at the root node (each call represents solving two LPs) */
extern
SCIP_Longint SCIPgetNRootStrongbranchs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in strong branching at the root node */
extern
SCIP_Longint SCIPgetNRootStrongbranchLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of pricing rounds performed so far at the current node */
extern
int SCIPgetNPriceRounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get current number of variables in the pricing store */
extern
int SCIPgetNPricevars(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of pricing variables found so far */
extern
int SCIPgetNPricevarsFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of pricing variables applied to the LPs */
extern
int SCIPgetNPricevarsApplied(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of separation rounds performed so far at the current node */
extern
int SCIPgetNSepaRounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of cuts found so far */
extern
int SCIPgetNCutsFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get number of cuts found so far in current separation round */
extern
int SCIPgetNCutsFoundRound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of cuts applied to the LPs */
extern
int SCIPgetNCutsApplied(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of constraints found in conflict analysis (conflict and reconvergence constraints) */
extern
SCIP_Longint SCIPgetNConflictConssFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get number of conflict constraints found so far at the current node */
extern
int SCIPgetNConflictConssFoundNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get total number of conflict constraints added to the problem */
extern
SCIP_Longint SCIPgetNConflictConssApplied(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets depth of current node, or -1 if no current node exists; in probing, the current node is the last probing node,
 *  such that the depth includes the probing path
 */
extern
int SCIPgetDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets depth of the focus node, or -1 if no focus node exists; the focus node is the currently processed node in the
 *  branching tree, excluding the nodes of the probing path
 */
extern
int SCIPgetFocusDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets maximal depth of all processed nodes in current branch and bound run (excluding probing nodes) */
extern
int SCIPgetMaxDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets maximal depth of all processed nodes over all branch and bound runs */
extern
int SCIPgetMaxTotalDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of backtracks, i.e. number of times, the new node was selected from the leaves queue */
extern
SCIP_Longint SCIPgetNBacktracks(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current plunging depth (successive times, a child was selected as next node) */
extern
int SCIPgetPlungeDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of active constraints at the current node */
extern
int SCIPgetNActiveConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total number of enabled constraints at the current node */
extern
int SCIPgetNEnabledConss(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets average dual bound of all unprocessed nodes for original problem */
extern
SCIP_Real SCIPgetAvgDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets average lower (dual) bound of all unprocessed nodes in transformed problem */
extern
SCIP_Real SCIPgetAvgLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets global dual bound for original problem */
extern
SCIP_Real SCIPgetDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets global lower (dual) bound in transformed problem */
extern
SCIP_Real SCIPgetLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets dual bound of the root node for the original problem */
extern
SCIP_Real SCIPgetDualboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets lower (dual) bound in transformed problem of the root node */
extern
SCIP_Real SCIPgetLowerboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets global primal bound (objective value of best solution or user objective limit) for the original problem */
extern
SCIP_Real SCIPgetPrimalbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets global upper (primal) bound in transformed problem (objective value of best solution or user objective limit) */
extern
SCIP_Real SCIPgetUpperbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets global cutoff bound in transformed problem: a sub problem with lower bound larger than the cutoff
 *  cannot contain a better feasible solution; usually, this bound is equal to the upper bound, but if the
 *  objective value is always integral, the cutoff bound is (nearly) one less than the upper bound;
 *  additionally, due to objective function domain propagation, the cutoff bound can be further reduced
 */
extern
SCIP_Real SCIPgetCutoffbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** updates the cutoff bound
 *
 *  @note the given cutoff bound has to better or equal to known one (SCIPgetCutoffbound())
 */
extern
SCIP_RETCODE SCIPupdateCutoffbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             cutoffbound         /**< new cutoff bound */
   );

/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 */
extern
SCIP_Bool SCIPisPrimalboundSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current gap |(primalbound - dualbound)/min(|primalbound|,|dualbound|)| if both bounds have same sign,
 *  or infinity, if they have opposite sign
 */
extern
SCIP_Real SCIPgetGap(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current gap |(upperbound - lowerbound)/min(|upperbound|,|lowerbound|)| in transformed problem if both bounds
 *  have same sign, or infinity, if they have opposite sign
 */
extern
SCIP_Real SCIPgetTransGap(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of feasible primal solutions found so far */
extern
SCIP_Longint SCIPgetNSolsFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of feasible primal solutions found so far, that improved the primal bound at the time they were found */
extern
SCIP_Longint SCIPgetNBestSolsFound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average pseudo cost value for the given direction over all variables */
extern
SCIP_Real SCIPgetAvgPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the average pseudo cost value for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPgetAvgPseudocostCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the average number of pseudo cost updates for the given direction over all variables */
extern
SCIP_Real SCIPgetAvgPseudocostCount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the average number of pseudo cost updates for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPgetAvgPseudocostCountCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the average pseudo cost score value over all variables, assuming a fractionality of 0.5 */
extern
SCIP_Real SCIPgetAvgPseudocostScore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average pseudo cost score value over all variables, assuming a fractionality of 0.5,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPgetAvgPseudocostScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average conflict score value over all variables */
extern
SCIP_Real SCIPgetAvgConflictScore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average conflict score value over all variables, only using the pseudo cost information of the current run */
extern
SCIP_Real SCIPgetAvgConflictScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average inference score value over all variables */
extern
SCIP_Real SCIPgetAvgConflictlengthScore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average conflictlength score value over all variables, only using the pseudo cost information of the current run */
extern
SCIP_Real SCIPgetAvgConflictlengthScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the average number of inferences found after branching in given direction over all variables */
extern
SCIP_Real SCIPgetAvgInferences(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of inferences found after branching in given direction over all variables,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPgetAvgInferencesCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the average inference score value over all variables */
extern
SCIP_Real SCIPgetAvgInferenceScore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average inference score value over all variables, only using the inference information information of the
 *  current run
 */
extern
SCIP_Real SCIPgetAvgInferenceScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the average number of cutoffs found after branching in given direction over all variables */
extern
SCIP_Real SCIPgetAvgCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of cutoffs found after branching in given direction over all variables,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPgetAvgCutoffsCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the average cutoff score value over all variables */
extern
SCIP_Real SCIPgetAvgCutoffScore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the average cutoff score value over all variables, only using the pseudo cost information of the current run */
extern
SCIP_Real SCIPgetAvgCutoffScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** outputs original problem to file stream */
extern
SCIP_RETCODE SCIPprintOrigProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format (or NULL for default CIP format)*/
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   );

/** outputs transformed problem of the current node to file stream */
extern
SCIP_RETCODE SCIPprintTransProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format (or NULL for default CIP format)*/
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   );

/** outputs solving statistics */
extern
SCIP_RETCODE SCIPprintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** outputs history statistics about branchings on variables */
extern
SCIP_RETCODE SCIPprintBranchingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** outputs node information display line */
extern
SCIP_RETCODE SCIPprintDisplayLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VERBLEVEL        verblevel           /**< minimal verbosity level to actually display the information line */
   );

/** gets total number of implications between variables that are stored in the implication graph */
extern
int SCIPgetNImplications(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** stores conflict graph of binary variables' implications into a file, which can be used as input for the DOT tool */
extern
SCIP_RETCODE SCIPwriteImplicationConflictGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name, or NULL for stdout */
   );

/**@} */




/*
 * timing methods
 */

/**@name Timing Methods */
/**@{ */

/** gets current time of day in seconds (standard time zone) */
extern
SCIP_Real SCIPgetTimeOfDay(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a clock using the default clock type */
extern
SCIP_RETCODE SCIPcreateClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** creates a clock counting the CPU user seconds */
extern
SCIP_RETCODE SCIPcreateCPUClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** creates a clock counting the wall clock seconds */
extern
SCIP_RETCODE SCIPcreateWallClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** frees a clock */
extern
SCIP_RETCODE SCIPfreeClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** resets the time measurement of a clock to zero and completely stops the clock */
extern
SCIP_RETCODE SCIPresetClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** starts the time measurement of a clock */
extern
SCIP_RETCODE SCIPstartClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** stops the time measurement of a clock */
extern
SCIP_RETCODE SCIPstopClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** starts the current solving time */
extern
SCIP_RETCODE SCIPstartSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** stops the current solving time in seconds */
extern
SCIP_RETCODE SCIPstopSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the measured time of a clock in seconds */
extern
SCIP_Real SCIPgetClockTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** sets the measured time of a clock to the given value in seconds */
extern
SCIP_RETCODE SCIPsetClockTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_Real             sec                 /**< time in seconds to set the clock's timer to */
   );

/** gets the current total SCIP time in seconds */
extern
SCIP_Real SCIPgetTotalTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the current solving time in seconds */
extern
SCIP_Real SCIPgetSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the current reading time in seconds */
extern
SCIP_Real SCIPgetReadingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the current presolving time in seconds */
extern
SCIP_Real SCIPgetPresolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */




/*
 * numeric values and comparisons
 */

/**@name Numerical Methods */
/**@{ */


/** returns value treated as zero */
extern
SCIP_Real SCIPepsilon(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns value treated as zero for sums of floating point values */
extern
SCIP_Real SCIPsumepsilon(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns feasibility tolerance for constraints */
extern
SCIP_Real SCIPfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns feasibility tolerance for reduced costs */
extern
SCIP_Real SCIPdualfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns convergence tolerance used in barrier algorithm */
extern
SCIP_Real SCIPbarrierconvtol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return the cutoff bound delta */
extern
SCIP_Real SCIPcutoffbounddelta(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the feasibility tolerance for constraints */
extern
SCIP_RETCODE SCIPchgFeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             feastol             /**< new feasibility tolerance for constraints */
   );

/** sets the feasibility tolerance for reduced costs */
extern
SCIP_RETCODE SCIPchgDualfeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dualfeastol         /**< new feasibility tolerance for reduced costs */
   );

/** sets the convergence tolerance used in barrier algorithm */
extern
SCIP_RETCODE SCIPchgBarrierconvtol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             barrierconvtol      /**< new convergence tolerance used in barrier algorithm */
   );

/** marks that some limit parameter was changed */
extern
void SCIPmarkLimitChanged(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

/** returns value treated as infinity */
extern
SCIP_Real SCIPinfinity(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** checks, if values are in range of epsilon */
extern
SCIP_Bool SCIPisEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) lower than val2 */
extern
SCIP_Bool SCIPisLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) greater than val2 */
extern
SCIP_Bool SCIPisLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) greater than val2 */
extern
SCIP_Bool SCIPisGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) lower than val2 */
extern
SCIP_Bool SCIPisGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is (positive) infinite */
extern
SCIP_Bool SCIPisInfinity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to be compared against infinity */
   );

/** checks, if value is in range epsilon of 0.0 */
extern
SCIP_Bool SCIPisZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is greater than epsilon */
extern
SCIP_Bool SCIPisPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is lower than -epsilon */
extern
SCIP_Bool SCIPisNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is integral within epsilon */
extern
SCIP_Bool SCIPisIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks whether the product val * scalar is integral in epsilon scaled by scalar */
extern
SCIP_Bool SCIPisScalingIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val,                /**< unscaled value to check for scaled integrality */
   SCIP_Real             scalar              /**< value to scale val with for checking for integrality */
   );

/** checks, if given fractional part is smaller than epsilon */
extern
SCIP_Bool SCIPisFracIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value + epsilon down to the next integer */
extern
SCIP_Real SCIPfloor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value - epsilon up to the next integer */
extern
SCIP_Real SCIPceil(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value to the nearest integer with epsilon tolerance */
extern
SCIP_Real SCIPround(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** returns fractional part of value, i.e. x - floor(x) in epsilon tolerance */
extern
SCIP_Real SCIPfrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to return fractional part for */
   );

/** checks, if values are in range of sumepsilon */
extern
SCIP_Bool SCIPisSumEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) lower than val2 */
extern
SCIP_Bool SCIPisSumLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
extern
SCIP_Bool SCIPisSumLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) greater than val2 */
extern
SCIP_Bool SCIPisSumGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
extern
SCIP_Bool SCIPisSumGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range sumepsilon of 0.0 */
extern
SCIP_Bool SCIPisSumZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is greater than sumepsilon */
extern
SCIP_Bool SCIPisSumPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is lower than -sumepsilon */
extern
SCIP_Bool SCIPisSumNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if relative difference of values is in range of feasibility tolerance */
extern
SCIP_Bool SCIPisFeasEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference val1 and val2 is lower than feasibility tolerance */
extern
SCIP_Bool SCIPisFeasLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than feasibility tolerance */
extern
SCIP_Bool SCIPisFeasLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than feastol */
extern
SCIP_Bool SCIPisFeasGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -feastol */
extern
SCIP_Bool SCIPisFeasGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range feasibility tolerance of 0.0 */
extern
SCIP_Bool SCIPisFeasZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is greater than feasibility tolerance */
extern
SCIP_Bool SCIPisFeasPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is lower than -feasibility tolerance */
extern
SCIP_Bool SCIPisFeasNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is integral within the LP feasibility bounds */
extern
SCIP_Bool SCIPisFeasIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if given fractional part is smaller than feastol */
extern
SCIP_Bool SCIPisFeasFracIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value + feasibility tolerance down to the next integer */
extern
SCIP_Real SCIPfeasFloor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value - feasibility tolerance up to the next integer */
extern
SCIP_Real SCIPfeasCeil(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value - feasibility tolerance up to the next integer in feasibility tolerance */
extern
SCIP_Real SCIPfeasRound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** returns fractional part of value, i.e. x - floor(x) */
extern
SCIP_Real SCIPfeasFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if the given new lower bound is tighter (w.r.t. bound strengthening epsilon) than the old one */
extern
SCIP_Bool SCIPisLbBetter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   );

/** checks, if the given new upper bound is tighter (w.r.t. bound strengthening epsilon) than the old one */
SCIP_Bool SCIPisUbBetter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newub,              /**< new upper bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   );

/** checks, if relative difference of values is in range of epsilon */
extern
SCIP_Bool SCIPisRelEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than epsilon */
extern
SCIP_Bool SCIPisRelLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
extern
SCIP_Bool SCIPisRelLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than epsilon */
extern
SCIP_Bool SCIPisRelGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
extern
SCIP_Bool SCIPisRelGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of values is in range of sumepsilon */
extern
SCIP_Bool SCIPisSumRelEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
extern
SCIP_Bool SCIPisSumRelLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
extern
SCIP_Bool SCIPisSumRelLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than sumepsilon */
extern
SCIP_Bool SCIPisSumRelGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -sumepsilon */
extern
SCIP_Bool SCIPisSumRelGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** Checks, if an iteratively updated value is reliable or should be recomputed from scratch.
 *  This is useful, if the value, e.g., the activity of a linear constraint or the pseudo objective value, gets a high
 *  absolute value during the optimization process which is later reduced significantly. In this case, the last digits
 *  were canceled out when increasing the value and are random after decreasing it.
 *  We do not consider the cancellations which can occur during increasing the absolute value because they just cannot
 *  be expressed using fixed precision floating point arithmetic, anymore.
 *  In order to get more reliable values, the idea is to always store the last reliable value, where increasing the
 *  absolute of the value is viewed as preserving reliability. Then, after each update, the new absolute value can be
 *  compared against the last reliable one with this method, checking whether it was decreased by a factor of at least
 *  "lp/recompfac" and should be recomputed.
 */
extern
SCIP_Bool SCIPisUpdateUnreliable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newvalue,           /**< new value after update */
   SCIP_Real             oldvalue            /**< old value, i.e., last reliable value */
   );

/** checks, if value is huge and should be handled separately (e.g., in activity computation) */
extern
SCIP_Bool SCIPisHugeValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to be checked whether it is huge */
   );

/** returns the minimum value that is regarded as huge and should be handled separately (e.g., in activity computation) */
extern
SCIP_Bool SCIPgetHugeValue(
   SCIP*                 scip                /**< SCIP data structure */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPinfinity(scip)                        SCIPsetInfinity((scip)->set)
#define SCIPisInfinity(scip, val)                 SCIPsetIsInfinity((scip)->set, val)        
#define SCIPisEQ(scip, val1, val2)                SCIPsetIsEQ((scip)->set, val1, val2)       
#define SCIPisLT(scip, val1, val2)                SCIPsetIsLT((scip)->set, val1, val2)       
#define SCIPisLE(scip, val1, val2)                SCIPsetIsLE((scip)->set, val1, val2)       
#define SCIPisGT(scip, val1, val2)                SCIPsetIsGT((scip)->set, val1, val2)       
#define SCIPisGE(scip, val1, val2)                SCIPsetIsGE((scip)->set, val1, val2)       
#define SCIPisZero(scip, val)                     SCIPsetIsZero((scip)->set, val)            
#define SCIPisPositive(scip, val)                 SCIPsetIsPositive((scip)->set, val)        
#define SCIPisNegative(scip, val)                 SCIPsetIsNegative((scip)->set, val)        
#define SCIPisIntegral(scip, val)                 SCIPsetIsIntegral((scip)->set, val)
#define SCIPisScalingIntegral(scip, val, scalar)  SCIPsetIsScalingIntegral((scip)->set, val, scalar)
#define SCIPisFracIntegral(scip, val)             SCIPsetIsFracIntegral((scip)->set, val)
#define SCIPfloor(scip, val)                      SCIPsetFloor((scip)->set, val)
#define SCIPceil(scip, val)                       SCIPsetCeil((scip)->set, val)
#define SCIPround(scip, val)                      SCIPsetRound((scip)->set, val)
#define SCIPfrac(scip, val)                       SCIPsetFrac((scip)->set, val)
                                                                                
#define SCIPisSumEQ(scip, val1, val2)             SCIPsetIsSumEQ((scip)->set, val1, val2)    
#define SCIPisSumLT(scip, val1, val2)             SCIPsetIsSumLT((scip)->set, val1, val2)    
#define SCIPisSumLE(scip, val1, val2)             SCIPsetIsSumLE((scip)->set, val1, val2)    
#define SCIPisSumGT(scip, val1, val2)             SCIPsetIsSumGT((scip)->set, val1, val2)    
#define SCIPisSumGE(scip, val1, val2)             SCIPsetIsSumGE((scip)->set, val1, val2)    
#define SCIPisSumZero(scip, val)                  SCIPsetIsSumZero((scip)->set, val)         
#define SCIPisSumPositive(scip, val)              SCIPsetIsSumPositive((scip)->set, val)     
#define SCIPisSumNegative(scip, val)              SCIPsetIsSumNegative((scip)->set, val)     
                                                                                
#define SCIPisFeasEQ(scip, val1, val2)            SCIPsetIsFeasEQ((scip)->set, val1, val2)   
#define SCIPisFeasLT(scip, val1, val2)            SCIPsetIsFeasLT((scip)->set, val1, val2)   
#define SCIPisFeasLE(scip, val1, val2)            SCIPsetIsFeasLE((scip)->set, val1, val2)   
#define SCIPisFeasGT(scip, val1, val2)            SCIPsetIsFeasGT((scip)->set, val1, val2)   
#define SCIPisFeasGE(scip, val1, val2)            SCIPsetIsFeasGE((scip)->set, val1, val2)   
#define SCIPisFeasZero(scip, val)                 SCIPsetIsFeasZero((scip)->set, val)        
#define SCIPisFeasPositive(scip, val)             SCIPsetIsFeasPositive((scip)->set, val)    
#define SCIPisFeasNegative(scip, val)             SCIPsetIsFeasNegative((scip)->set, val)    
#define SCIPisFeasIntegral(scip, val)             SCIPsetIsFeasIntegral((scip)->set, val)
#define SCIPisFeasFracIntegral(scip, val)         SCIPsetIsFeasFracIntegral((scip)->set, val)
#define SCIPfeasFloor(scip, val)                  SCIPsetFeasFloor((scip)->set, val)
#define SCIPfeasCeil(scip, val)                   SCIPsetFeasCeil((scip)->set, val)
#define SCIPfeasRound(scip, val)                  SCIPsetFeasRound((scip)->set, val)
#define SCIPfeasFrac(scip, val)                   SCIPsetFeasFrac((scip)->set, val)

                                                                                
#define SCIPisLbBetter(scip, newlb, oldlb, oldub) SCIPsetIsLbBetter(scip->set, newlb, oldlb, oldub)
#define SCIPisUbBetter(scip, newub, oldlb, oldub) SCIPsetIsUbBetter(scip->set, newub, oldlb, oldub)
                                                                                
#define SCIPisRelEQ(scip, val1, val2)             SCIPsetIsRelEQ((scip)->set, val1, val2)    
#define SCIPisRelLT(scip, val1, val2)             SCIPsetIsRelLT((scip)->set, val1, val2)    
#define SCIPisRelLE(scip, val1, val2)             SCIPsetIsRelLE((scip)->set, val1, val2)    
#define SCIPisRelGT(scip, val1, val2)             SCIPsetIsRelGT((scip)->set, val1, val2)    
#define SCIPisRelGE(scip, val1, val2)             SCIPsetIsRelGE((scip)->set, val1, val2)    
                                                                                
#define SCIPisSumRelEQ(scip, val1, val2)          SCIPsetIsSumRelEQ((scip)->set, val1, val2) 
#define SCIPisSumRelLT(scip, val1, val2)          SCIPsetIsSumRelLT((scip)->set, val1, val2) 
#define SCIPisSumRelLE(scip, val1, val2)          SCIPsetIsSumRelLE((scip)->set, val1, val2) 
#define SCIPisSumRelGT(scip, val1, val2)          SCIPsetIsSumRelGT((scip)->set, val1, val2) 
#define SCIPisSumRelGE(scip, val1, val2)          SCIPsetIsSumRelGE((scip)->set, val1, val2) 

#define SCIPisUpdateUnreliable(scip, newval, oldval) SCIPsetIsUpdateUnreliable((scip)->set, newval, oldval)
#define SCIPisHugeValue(scip, val) SCIPsetIsHugeValue((scip)->set, val)
#define SCIPgetHugeValue(scip) SCIPsetGetHugeValue((scip)->set)
#endif

/** outputs a real number, or "+infinity", or "-infinity" to a file */
extern
void SCIPprintReal(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Real             val,                /**< value to print */
   int                   width,              /**< width of the field */
   int                   precision           /**< number of significant digits printed */
   );

/**@} */




/*
 * memory management
 */

/**@name Memory Management */
/**@{ */

#define SCIPallocMemory(scip,ptr)               ( (BMSallocMemory((ptr)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocMemoryArray(scip,ptr,num)      ( (BMSallocMemoryArray((ptr), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocClearMemoryArray(scip,ptr,num) ( (BMSallocClearMemoryArray((ptr), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocMemorySize(scip,ptr,size)      ( (BMSallocMemorySize((ptr), (size)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocMemoryArray(scip,ptr,newnum) ( (BMSreallocMemoryArray((ptr), (newnum)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocMemorySize(scip,ptr,newsize) ( (BMSreallocMemorySize((ptr), (newsize)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateMemory(scip, ptr, source)  ( (BMSduplicateMemory((ptr), (source)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateMemoryArray(scip, ptr, source, num) \
                                                     ( (BMSduplicateMemoryArray((ptr), (source), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPfreeMemory(scip,ptr)                BMSfreeMemory(ptr)
#define SCIPfreeMemoryNull(scip,ptr)            BMSfreeMemoryNull(ptr)
#define SCIPfreeMemoryArray(scip,ptr)           BMSfreeMemoryArray(ptr)
#define SCIPfreeMemoryArrayNull(scip,ptr)       BMSfreeMemoryArrayNull(ptr)
#define SCIPfreeMemorySize(scip,ptr)            BMSfreeMemorySize(ptr)
#define SCIPfreeMemorySizeNull(scip,ptr)        BMSfreeMemorySizeNull(ptr)

#define SCIPallocBlockMemory(scip,ptr)          ( (BMSallocBlockMemory(SCIPblkmem(scip), (ptr)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocBlockMemoryArray(scip,ptr,num) ( (BMSallocBlockMemoryArray(SCIPblkmem(scip), (ptr), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocBlockMemorySize(scip,ptr,size) ( (BMSallocBlockMemorySize(SCIPblkmem(scip), (ptr), (size)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocBlockMemoryArray(scip,ptr,oldnum,newnum) \
                                                     ( (BMSreallocBlockMemoryArray(SCIPblkmem(scip), (ptr), (oldnum), (newnum)) \
                                                       == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocBlockMemorySize(scip,ptr,oldsize,newsize) \
                                                     ( (BMSreallocBlockMemorySize(SCIPblkmem(scip), (ptr), (oldsize), (newsize)) \
                                                       == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateBlockMemory(scip, ptr, source) \
                                                     ( (BMSduplicateBlockMemory(SCIPblkmem(scip), (ptr), (source)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateBlockMemoryArray(scip, ptr, source, num) \
                                                     ( (BMSduplicateBlockMemoryArray(SCIPblkmem(scip), (ptr), (source), (num)) \
                                                       == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPensureBlockMemoryArray(scip,ptr,arraysizeptr,minsize) \
                                                     ( (SCIPensureBlockMemoryArray_call((scip), (void**)(ptr), sizeof(**(ptr)), \
                                                        (arraysizeptr), (minsize))) )
#define SCIPfreeBlockMemory(scip,ptr)           BMSfreeBlockMemory(SCIPblkmem(scip), (ptr))
#define SCIPfreeBlockMemoryNull(scip,ptr)       BMSfreeBlockMemoryNull(SCIPblkmem(scip), (ptr))
#define SCIPfreeBlockMemoryArray(scip,ptr,num)  BMSfreeBlockMemoryArray(SCIPblkmem(scip), (ptr), (num))
#define SCIPfreeBlockMemoryArrayNull(scip,ptr,num) \
                                                     BMSfreeBlockMemoryArrayNull(SCIPblkmem(scip), (ptr), (num))
#define SCIPfreeBlockMemorySize(scip,ptr,size)  BMSfreeBlockMemorySize(SCIPblkmem(scip), (ptr), (size))
#define SCIPfreeBlockMemorySizeNull(scip,ptr,size) \
                                                     BMSfreeBlockMemorySizeNull(SCIPblkmem(scip), (ptr), (size))

#define SCIPallocBuffer(scip,ptr)               SCIPallocBufferSize(scip, (void**)(ptr), (int)sizeof(**(ptr)))
#define SCIPallocBufferArray(scip,ptr,num)      SCIPallocBufferSize(scip, (void**)(ptr), (num)*(int)sizeof(**(ptr)))
#define SCIPreallocBufferArray(scip,ptr,num)    SCIPreallocBufferSize(scip, (void**)(ptr), (num)*(int)sizeof(**(ptr)))
#define SCIPduplicateBuffer(scip,ptr,source)    SCIPduplicateBufferSize(scip, (void**)(ptr), source, (int)sizeof(**(ptr)))
#define SCIPduplicateBufferArray(scip,ptr,source,num) \
                                                     SCIPduplicateBufferSize(scip, (void**)(ptr), source, (num)*(int)sizeof(**(ptr)))
#define SCIPfreeBuffer(scip,ptr)                SCIPfreeBufferSize(scip, (void**)(ptr), 0)
#define SCIPfreeBufferNull(scip,ptr)            { if( *(ptr) != NULL ) SCIPfreeBuffer(scip, ptr); }
#define SCIPfreeBufferArray(scip,ptr)           SCIPfreeBufferSize(scip, (void**)(ptr), 0)
#define SCIPfreeBufferArrayNull(scip,ptr)       { if( *(ptr) != NULL ) SCIPfreeBufferArray(scip, ptr); }


/** returns block memory to use at the current time */
extern
BMS_BLKMEM* SCIPblkmem(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the total number of bytes used in block memory */
extern
SCIP_Longint SCIPgetMemUsed(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** calculate memory size for dynamically allocated arrays */
extern
int SCIPcalcMemGrowSize(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   num                 /**< minimum number of entries to store */
   );

/** extends a dynamically allocated block memory array to be able to store at least the given number of elements;
 *  use SCIPensureBlockMemoryArray() define to call this method!
 */
extern
SCIP_RETCODE SCIPensureBlockMemoryArray_call(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                arrayptr,           /**< pointer to dynamically sized array */
   size_t                elemsize,           /**< size in bytes of each element in array */
   int*                  arraysize,          /**< pointer to current array size */
   int                   minsize             /**< required minimal array size */
   );

/** gets a memory buffer with at least the given size */
extern
SCIP_RETCODE SCIPallocBufferSize(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                ptr,                /**< pointer to store the buffer */
   int                   size                /**< required size in bytes of buffer */
   );

/** allocates a memory buffer with at least the given size and copies the given memory into the buffer */
extern
SCIP_RETCODE SCIPduplicateBufferSize(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                ptr,                /**< pointer to store the buffer */
   const void*           source,             /**< memory block to copy into the buffer */
   int                   size                /**< required size in bytes of buffer */
   );

/** reallocates a memory buffer to at least the given size */
extern
SCIP_RETCODE SCIPreallocBufferSize(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                ptr,                /**< pointer to the buffer */
   int                   size                /**< required size in bytes of buffer */
   );

/** frees a memory buffer */
extern
void SCIPfreeBufferSize(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                ptr,                /**< pointer to the buffer */
   int                   dummysize           /**< used to get a safer define for SCIPfreeBuffer() and SCIPfreeBufferArray() */
   );

/** prints output about used memory */
extern
void SCIPprintMemoryDiagnostic(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */




/*
 * dynamic arrays
 */

/**@name Dynamic Arrays */
/**@{ */

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** creates a dynamic array of real values */
extern
SCIP_RETCODE SCIPcreateRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY**      realarray           /**< pointer to store the real array */
   );

/** frees a dynamic array of real values */
extern
SCIP_RETCODE SCIPfreeRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY**      realarray           /**< pointer to the real array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
SCIP_RETCODE SCIPextendRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic real array */
extern
SCIP_RETCODE SCIPclearRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   );

/** gets value of entry in dynamic array */
extern
SCIP_Real SCIPgetRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
SCIP_RETCODE SCIPsetRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx,                /**< array index to set value for */
   SCIP_Real             val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
SCIP_RETCODE SCIPincRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx,                /**< array index to increase value for */
   SCIP_Real             incval              /**< value to increase array index */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPgetRealarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPgetRealarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   );

/** creates a dynamic array of int values */
extern
SCIP_RETCODE SCIPcreateIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY**       intarray            /**< pointer to store the int array */
   );

/** frees a dynamic array of int values */
extern
SCIP_RETCODE SCIPfreeIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY**       intarray            /**< pointer to the int array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
SCIP_RETCODE SCIPextendIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic int array */
extern
SCIP_RETCODE SCIPclearIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   );

/** gets value of entry in dynamic array */
extern
int SCIPgetIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
SCIP_RETCODE SCIPsetIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx,                /**< array index to set value for */
   int                   val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
SCIP_RETCODE SCIPincIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx,                /**< array index to increase value for */
   int                   incval              /**< value to increase array index */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPgetIntarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPgetIntarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   );

/** creates a dynamic array of bool values */
extern
SCIP_RETCODE SCIPcreateBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY**      boolarray           /**< pointer to store the bool array */
   );

/** frees a dynamic array of bool values */
extern
SCIP_RETCODE SCIPfreeBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY**      boolarray           /**< pointer to the bool array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
SCIP_RETCODE SCIPextendBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic bool array */
extern
SCIP_RETCODE SCIPclearBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** gets value of entry in dynamic array */
extern
SCIP_Bool SCIPgetBoolarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
SCIP_RETCODE SCIPsetBoolarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx,                /**< array index to set value for */
   SCIP_Bool             val                 /**< value to set array index to */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPgetBoolarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPgetBoolarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** creates a dynamic array of pointers */
extern
SCIP_RETCODE SCIPcreatePtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY**       ptrarray            /**< pointer to store the int array */
   );

/** frees a dynamic array of pointers */
extern
SCIP_RETCODE SCIPfreePtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY**       ptrarray            /**< pointer to the int array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
SCIP_RETCODE SCIPextendPtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic pointer array */
extern
SCIP_RETCODE SCIPclearPtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic int array */
   );

/** gets value of entry in dynamic array */
extern
void* SCIPgetPtrarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
SCIP_RETCODE SCIPsetPtrarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   idx,                /**< array index to set value for */
   void*                 val                 /**< value to set array index to */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPgetPtrarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPgetPtrarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPcreateRealarray(scip, realarray)                    SCIPrealarrayCreate(realarray, SCIPblkmem(scip))
#define SCIPfreeRealarray(scip, realarray)                      SCIPrealarrayFree(realarray)
#define SCIPextendRealarray(scip, realarray, minidx, maxidx)    SCIPrealarrayExtend(realarray, (scip)->set, minidx, maxidx)
#define SCIPclearRealarray(scip, realarray)                     SCIPrealarrayClear(realarray)
#define SCIPgetRealarrayVal(scip, realarray, idx)               SCIPrealarrayGetVal(realarray, idx)
#define SCIPsetRealarrayVal(scip, realarray, idx, val)          SCIPrealarraySetVal(realarray, (scip)->set, idx, val)
#define SCIPincRealarrayVal(scip, realarray, idx, incval)       SCIPrealarrayIncVal(realarray, scip->set, idx, incval)
#define SCIPgetRealarrayMinIdx(scip, realarray)                 SCIPrealarrayGetMinIdx(realarray)
#define SCIPgetRealarrayMaxIdx(scip, realarray)                 SCIPrealarrayGetMaxIdx(realarray)

#define SCIPcreateIntarray(scip, intarray)                      SCIPintarrayCreate(intarray, SCIPblkmem(scip))
#define SCIPfreeIntarray(scip, intarray)                        SCIPintarrayFree(intarray)
#define SCIPextendIntarray(scip, intarray, minidx, maxidx)      SCIPintarrayExtend(intarray, (scip)->set, minidx, maxidx)
#define SCIPclearIntarray(scip, intarray)                       SCIPintarrayClear(intarray)
#define SCIPgetIntarrayVal(scip, intarray, idx)                 SCIPintarrayGetVal(intarray, idx)
#define SCIPsetIntarrayVal(scip, intarray, idx, val)            SCIPintarraySetVal(intarray, (scip)->set, idx, val)
#define SCIPincIntarrayVal(scip, intarray, idx, incval)         SCIPintarrayIncVal(intarray, scip->set, idx, incval)
#define SCIPgetIntarrayMinIdx(scip, intarray)                   SCIPintarrayGetMinIdx(intarray)
#define SCIPgetIntarrayMaxIdx(scip, intarray)                   SCIPintarrayGetMaxIdx(intarray)

#define SCIPcreateBoolarray(scip, boolarray)                    SCIPboolarrayCreate(boolarray, SCIPblkmem(scip))
#define SCIPfreeBoolarray(scip, boolarray)                      SCIPboolarrayFree(boolarray)
#define SCIPextendBoolarray(scip, boolarray, minidx, maxidx)    SCIPboolarrayExtend(boolarray, (scip)->set, minidx, maxidx)
#define SCIPclearBoolarray(scip, boolarray)                     SCIPboolarrayClear(boolarray)
#define SCIPgetBoolarrayVal(scip, boolarray, idx)               SCIPboolarrayGetVal(boolarray, idx)
#define SCIPsetBoolarrayVal(scip, boolarray, idx, val)          SCIPboolarraySetVal(boolarray, (scip)->set, idx, val)
#define SCIPgetBoolarrayMinIdx(scip, boolarray)                 SCIPboolarrayGetMinIdx(boolarray)
#define SCIPgetBoolarrayMaxIdx(scip, boolarray)                 SCIPboolarrayGetMaxIdx(boolarray)

#define SCIPcreatePtrarray(scip, ptrarray)                      SCIPptrarrayCreate(ptrarray, SCIPblkmem(scip))
#define SCIPfreePtrarray(scip, ptrarray)                        SCIPptrarrayFree(ptrarray)
#define SCIPextendPtrarray(scip, ptrarray, minidx, maxidx)      SCIPptrarrayExtend(ptrarray, (scip)->set, minidx, maxidx)
#define SCIPclearPtrarray(scip, ptrarray)                       SCIPptrarrayClear(ptrarray)
#define SCIPgetPtrarrayVal(scip, ptrarray, idx)                 SCIPptrarrayGetVal(ptrarray, idx)
#define SCIPsetPtrarrayVal(scip, ptrarray, idx, val)            SCIPptrarraySetVal(ptrarray, (scip)->set, idx, val)
#define SCIPgetPtrarrayMinIdx(scip, ptrarray)                   SCIPptrarrayGetMinIdx(ptrarray)
#define SCIPgetPtrarrayMaxIdx(scip, ptrarray)                   SCIPptrarrayGetMaxIdx(ptrarray)

#endif

/**@} */
#ifdef __cplusplus
}
#endif

#endif
