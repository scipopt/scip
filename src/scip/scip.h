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
#pragma ident "@(#) $Id: scip.h,v 1.143 2004/06/30 14:17:01 bzfpfend Exp $"

/**@file   scip.h
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_H__
#define __SCIP_H__


#include <stdio.h>

#include "def.h"
#include "message.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_result.h"
#include "type_clock.h"
#include "type_misc.h"
#include "type_paramset.h"
#include "type_event.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_prob.h"
#include "type_tree.h"
#include "type_scip.h"

#include "type_branch.h"
#include "type_conflict.h"
#include "type_cons.h"
#include "type_dialog.h"
#include "type_disp.h"
#include "type_heur.h"
#include "type_nodesel.h"
#include "type_presol.h"
#include "type_pricer.h"
#include "type_reader.h"
#include "type_sepa.h"

/* include public interfaces, s.t. the user only needs to include scip.h */
#include "pub_branch.h"
#include "pub_conflict.h"
#include "pub_cons.h"
#include "pub_cutpool.h"
#include "pub_dialog.h"
#include "pub_disp.h"
#include "pub_event.h"
#include "pub_heur.h"
#include "pub_lp.h"
#include "pub_misc.h"
#include "pub_nodesel.h"
#include "pub_paramset.h"
#include "pub_presol.h"
#include "pub_pricer.h"
#include "pub_reader.h"
#include "pub_sepa.h"
#include "pub_sol.h"
#include "pub_tree.h"
#include "pub_var.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "struct_scip.h"
#include "set.h"
#endif




/*
 * miscellaneous methods
 */

/**@name Miscellaneos Methods */
/**@{ */

/** returns scip version number */
extern
Real SCIPversion(
   void
   );

/** prints a version information line to a file stream */
extern
void SCIPprintVersion(
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** prints error message for the given SCIP return code */
extern
void SCIPprintError(
   RETCODE          retcode,            /**< SCIP return code causing the error */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/**@} */




/*
 * general SCIP methods
 */

/**@name General SCIP Methods */
/**@{ */

/** creates and initializes SCIP data structures */
extern
RETCODE SCIPcreate(
   SCIP**           scip                /**< pointer to SCIP data structure */
   );

/** frees SCIP data structures */
extern
RETCODE SCIPfree(
   SCIP**           scip                /**< pointer to SCIP data structure */
   );

/** prints a message depending on the verbosity level */
extern
void SCIPmessage(
   SCIP*            scip,               /**< SCIP data structure */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   );

/** returns current stage of SCIP */
extern
STAGE SCIPstage(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns whether the current stage belongs to the transformed problem space */
extern
Bool SCIPisTransformed(
   SCIP*            scip                /**< SCIP data structure */
   );
/** returns whether the solution process should be provably correct */
extern
Bool SCIPisExactSolve(
   SCIP*            scip                /**< SCIP data structure */
   );

/**@} */




/*
 * parameter settings
 */

/**@name Parameter Methods */
/**@{ */

/** creates a Bool parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddBoolParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Bool             defaultvalue,       /**< default value of the parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddIntParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   int*             valueptr,           /**< pointer to store the current parameter value, or NULL */
   int              defaultvalue,       /**< default value of the parameter */
   int              minvalue,           /**< minimum value for parameter */
   int              maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a Longint parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddLongintParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   Longint          defaultvalue,       /**< default value of the parameter */
   Longint          minvalue,           /**< minimum value for parameter */
   Longint          maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a Real parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddRealParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Real             defaultvalue,       /**< default value of the parameter */
   Real             minvalue,           /**< minimum value for parameter */
   Real             maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddCharParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   char*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   char             defaultvalue,       /**< default value of the parameter */
   const char*      allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddStringParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   char**           valueptr,           /**< pointer to store the current parameter value, or NULL */
   const char*      defaultvalue,       /**< default value of the parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** gets the value of an existing Bool parameter */
extern
RETCODE SCIPgetBoolParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Bool*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Int parameter */
extern
RETCODE SCIPgetIntParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   int*             value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Longint parameter */
extern
RETCODE SCIPgetLongintParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Longint*         value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Real parameter */
extern
RETCODE SCIPgetRealParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Real*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Char parameter */
extern
RETCODE SCIPgetCharParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   char*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing String parameter */
extern
RETCODE SCIPgetStringParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   char**           value               /**< pointer to store the parameter */
   );

/** changes the value of an existing Bool parameter */
extern
RETCODE SCIPsetBoolParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Bool             value               /**< new value of the parameter */
   );

/** changes the value of an existing Int parameter */
extern
RETCODE SCIPsetIntParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   int              value               /**< new value of the parameter */
   );

/** changes the value of an existing Longint parameter */
extern
RETCODE SCIPsetLongintParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Longint          value               /**< new value of the parameter */
   );

/** changes the value of an existing Real parameter */
extern
RETCODE SCIPsetRealParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Real             value               /**< new value of the parameter */
   );

/** changes the value of an existing Char parameter */
extern
RETCODE SCIPsetCharParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   char             value               /**< new value of the parameter */
   );

/** changes the value of an existing String parameter */
extern
RETCODE SCIPsetStringParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      value               /**< new value of the parameter */
   );

/** reads parameters from a file */
extern
RETCODE SCIPreadParams(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< file name */
   );

/** writes all parameters in the parameter set to a file */
extern
RETCODE SCIPwriteParams(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename,           /**< file name, or NULL for stdout */
   Bool             comments,           /**< should parameter descriptions be written as comments? */
   Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
   );

/** returns the array of all available SCIP parameters */
extern
PARAM** SCIPgetParams(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the total number of all available SCIP parameters */
extern
int SCIPgetNParams(
   SCIP*            scip                /**< SCIP data structure */
   );

/**@} */




/*
 * SCIP user functionality methods: managing plugins
 */

/**@name SCIP User Functionality Methods: Managing Plugins */
/**@{ */

/** creates a reader and includes it in SCIP */
extern
RETCODE SCIPincludeReader(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
   DECL_READERFREE  ((*readerfree)),    /**< destructor of reader */
   DECL_READERREAD  ((*readerread)),    /**< read method */
   READERDATA*      readerdata          /**< reader data */
   );

/** returns the reader of the given name, or NULL if not existing */
extern
READER* SCIPfindReader(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of constraint handler */
   );

/** returns the array of currently available readers */
extern
READER** SCIPgetReaders(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available readers */
extern
int SCIPgetNReaders(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates a variable pricer and includes it in SCIP
 *  To use the variable pricer for solving a problem, it first has to be activated with a call to SCIPactivatePricer().
 *  This should be done during the problem creation stage.
 */
extern
RETCODE SCIPincludePricer(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of variable pricer */
   const char*      desc,               /**< description of variable pricer */
   int              priority,           /**< priority of the variable pricer */
   DECL_PRICERFREE  ((*pricerfree)),    /**< destructor of variable pricer */
   DECL_PRICERINIT  ((*pricerinit)),    /**< initialize variable pricer */
   DECL_PRICEREXIT  ((*pricerexit)),    /**< deinitialize variable pricer */
   DECL_PRICERREDCOST((*pricerredcost)),/**< reduced cost pricing method of variable pricer for feasible LPs */
   DECL_PRICERFARKAS((*pricerfarkas)),  /**< farkas pricing method of variable pricer for infeasible LPs */
   PRICERDATA*      pricerdata          /**< variable pricer data */
   );

/** returns the variable pricer of the given name, or NULL if not existing */
extern
PRICER* SCIPfindPricer(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of variable pricer */
   );

/** returns the array of currently available variable pricers; active pricers are in the first slots of the array */
extern
PRICER** SCIPgetPricers(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available variable pricers */
extern
int SCIPgetNPricers(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently active variable pricers, that are used in the LP solving loop */
extern
int SCIPgetNActivePricers(
   SCIP*            scip                /**< SCIP data structure */
   );

/** sets the priority of a variable pricer */
extern
RETCODE SCIPsetPricerPriority(
   SCIP*            scip,               /**< SCIP data structure */
   PRICER*          pricer,             /**< variable pricer */
   int              priority            /**< new priority of the variable pricer */
   );

/** activates pricer to be used for the current problem
 *  This method should be called during the problem creation stage for all pricers that are necessary to solve
 *  the problem model.
 *  The pricers are automatically deactivated when the problem is freed.
 */
extern
RETCODE SCIPactivatePricer(
   SCIP*            scip,               /**< SCIP data structure */
   PRICER*          pricer              /**< variable pricer */
   );

/** creates a constraint handler and includes it in SCIP */
extern
RETCODE SCIPincludeConshdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              chckpriority,       /**< priority of the constraint handler for checking infeasibility */
   int              sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int              propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int              eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   DECL_CONSINIT    ((*consinit)),      /**< initialize constraint handler */
   DECL_CONSEXIT    ((*consexit)),      /**< deinitialize constraint handler */
   DECL_CONSINITPRE ((*consinitpre)),   /**< presolving initialization method of constraint handler */
   DECL_CONSEXITPRE ((*consexitpre)),   /**< presolving deinitialization method of constraint handler */
   DECL_CONSINITSOL ((*consinitsol)),   /**< solving process initialization method of constraint handler */
   DECL_CONSEXITSOL ((*consexitsol)),   /**< solving process deinitialization method of constraint handler */
   DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSINITLP  ((*consinitlp)),    /**< initialize LP with relaxations of "initial" constraints */
   DECL_CONSSEPA    ((*conssepa)),      /**< separate cutting planes */
   DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method */
   DECL_CONSRESCVAR ((*consrescvar)),   /**< conflict variable resolving method */
   DECL_CONSLOCK    ((*conslock)),      /**< variable rounding lock method */
   DECL_CONSUNLOCK  ((*consunlock)),    /**< variable rounding unlock method */
   DECL_CONSACTIVE  ((*consactive)),    /**< activation notification method */
   DECL_CONSDEACTIVE((*consdeactive)),  /**< deactivation notification method */
   DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   DECL_CONSPRINT   ((*consprint)),     /**< constraint display method */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   );

/** returns the constraint handler of the given name, or NULL if not existing */
extern
CONSHDLR* SCIPfindConshdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of constraint handler */
   );

/** returns the array of currently available constraint handlers */
extern
CONSHDLR** SCIPgetConshdlrs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available constraint handlers */
extern
int SCIPgetNConshdlrs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates a conflict handler and includes it in SCIP */
extern
RETCODE SCIPincludeConflicthdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of conflict handler */
   const char*      desc,               /**< description of conflict handler */
   int              priority,           /**< priority of the conflict handler */
   DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   );

/** returns the conflict handler of the given name, or NULL if not existing */
extern
CONFLICTHDLR* SCIPfindConflicthdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of conflict handler */
   );

/** returns the array of currently available conflict handlers */
extern
CONFLICTHDLR** SCIPgetConflicthdlrs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available conflict handlers */
extern
int SCIPgetNConflicthdlrs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** sets the priority of a conflict handler */
extern
RETCODE SCIPsetConflicthdlrPriority(
   SCIP*            scip,               /**< SCIP data structure */
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   int              priority            /**< new priority of the conflict handler */
   );

/** creates a presolver and includes it in SCIP */
extern
RETCODE SCIPincludePresol(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of presolver */
   const char*      desc,               /**< description of presolver */
   int              priority,           /**< priority of the presolver */
   DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver to free user data (called when SCIP is exiting) */
   DECL_PRESOLINIT  ((*presolinit)),    /**< initialization method of presolver (called after problem was transformed) */
   DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialization method of presolver (called before transformed problem is freed) */
   DECL_PRESOLINITPRE((*presolinitpre)),/**< presolving initialization method of presolver (called when presolving is about to begin) */
   DECL_PRESOLEXITPRE((*presolexitpre)),/**< presolving deinitialization method of presolver (called after presolving has been finished) */
   DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   PRESOLDATA*      presoldata          /**< presolver data */
   );

/** returns the presolver of the given name, or NULL if not existing */
extern
PRESOL* SCIPfindPresol(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of presolver */
   );

/** returns the array of currently available presolvers */
extern
PRESOL** SCIPgetPresols(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available presolvers */
extern
int SCIPgetNPresols(
   SCIP*            scip                /**< SCIP data structure */
   );

/** sets the priority of a presolver */
extern
RETCODE SCIPsetPresolPriority(
   SCIP*            scip,               /**< SCIP data structure */
   PRESOL*          presol,             /**< presolver */
   int              priority            /**< new priority of the presolver */
   );

/** creates a separator and includes it in SCIP */
extern
RETCODE SCIPincludeSepa(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of separator */
   const char*      desc,               /**< description of separator */
   int              priority,           /**< priority of the separator */
   int              freq,               /**< frequency for calling separator */
   DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   DECL_SEPAINIT    ((*sepainit)),      /**< initialize separator */
   DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialize separator */
   DECL_SEPAEXEC    ((*sepaexec)),      /**< execution method of separator */
   SEPADATA*        sepadata            /**< separator data */
   );

/** returns the separator of the given name, or NULL if not existing */
extern
SEPA* SCIPfindSepa(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of separator */
   );

/** returns the array of currently available separators */
extern
SEPA** SCIPgetSepas(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available separators */
extern
int SCIPgetNSepas(
   SCIP*            scip                /**< SCIP data structure */
   );

/** sets the priority of a separator */
extern
RETCODE SCIPsetSepaPriority(
   SCIP*            scip,               /**< SCIP data structure */
   SEPA*            sepa,               /**< primal sepaistic */
   int              priority            /**< new priority of the separator */
   );

/** creates a primal heuristic and includes it in SCIP */
extern
RETCODE SCIPincludeHeur(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of primal heuristic */
   const char*      desc,               /**< description of primal heuristic */
   char             dispchar,           /**< display character of primal heuristic */
   int              priority,           /**< priority of the primal heuristic */
   int              freq,               /**< frequency for calling primal heuristic */
   int              freqofs,            /**< frequency offset for calling primal heuristic */
   int              maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   Bool             pseudonodes,        /**< call heuristic at nodes where only a pseudo solution exist? */
   Bool             duringplunging,     /**< call heuristic during plunging? */
   DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   DECL_HEURINIT    ((*heurinit)),      /**< initialize primal heuristic */
   DECL_HEUREXIT    ((*heurexit)),      /**< deinitialize primal heuristic */
   DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   HEURDATA*        heurdata            /**< primal heuristic data */
   );

/** returns the primal heuristic of the given name, or NULL if not existing */
extern
HEUR* SCIPfindHeur(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of primal heuristic */
   );

/** returns the array of currently available primal heuristics */
extern
HEUR** SCIPgetHeurs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available primal heuristics */
extern
int SCIPgetNHeurs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** sets the priority of a primal heuristic */
extern
RETCODE SCIPsetHeurPriority(
   SCIP*            scip,               /**< SCIP data structure */
   HEUR*            heur,               /**< primal heuristic */
   int              priority            /**< new priority of the primal heuristic */
   );

/** creates an event handler and includes it in SCIP */
extern
RETCODE SCIPincludeEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of event handler */
   const char*      desc,               /**< description of event handler */
   DECL_EVENTFREE   ((*eventfree)),     /**< destructor of event handler */
   DECL_EVENTINIT   ((*eventinit)),     /**< initialize event handler */
   DECL_EVENTEXIT   ((*eventexit)),     /**< deinitialize event handler */
   DECL_EVENTDELETE ((*eventdelete)),   /**< free specific event data */
   DECL_EVENTEXEC   ((*eventexec)),     /**< execute event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   );

/** returns the event handler of the given name, or NULL if not existing */
extern
EVENTHDLR* SCIPfindEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   );

/** returns the array of currently available event handlers */
extern
EVENTHDLR** SCIPgetEventhdlrs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available event handlers */
extern
int SCIPgetNEventhdlrs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates a node selector and includes it in SCIP */
extern
RETCODE SCIPincludeNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of node selector */
   const char*      desc,               /**< description of node selector */
   int              stdpriority,        /**< priority of the node selector in standard mode */
   int              memsavepriority,    /**< priority of the node selector in memory saving mode */
   Bool             lowestboundfirst,   /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
   DECL_NODESELFREE ((*nodeselfree)),   /**< destructor of node selector */
   DECL_NODESELINIT ((*nodeselinit)),   /**< initialize node selector */
   DECL_NODESELEXIT ((*nodeselexit)),   /**< deinitialize node selector */
   DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   NODESELDATA*     nodeseldata         /**< node selector data */
   );

/** returns the node selector of the given name, or NULL if not existing */
extern
NODESEL* SCIPfindNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   );

/** returns the array of currently available node selectors */
extern
NODESEL** SCIPgetNodesels(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available node selectors */
extern
int SCIPgetNNodesels(
   SCIP*            scip                /**< SCIP data structure */
   );

/** sets the priority of a node selector in standard mode */
extern
RETCODE SCIPsetNodeselStdPriority(
   SCIP*            scip,               /**< SCIP data structure */
   NODESEL*         nodesel,            /**< node selector */
   int              priority            /**< new standard priority of the node selector */
   );

/** sets the priority of a node selector in memory saving mode */
extern
RETCODE SCIPsetNodeselMemsavePriority(
   SCIP*            scip,               /**< SCIP data structure */
   NODESEL*         nodesel,            /**< node selector */
   int              priority            /**< new memory saving priority of the node selector */
   );

/** returns the currently used node selector */
extern
NODESEL* SCIPgetNodesel(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates a branching rule and includes it in SCIP */
extern
RETCODE SCIPincludeBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   int              maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   DECL_BRANCHINIT  ((*branchinit)),    /**< initialize branching rule */
   DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialize branching rule */
   DECL_BRANCHEXECLP((*branchexeclp)),  /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXECPS((*branchexecps)),  /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   );

/** returns the branching rule of the given name, or NULL if not existing */
extern
BRANCHRULE* SCIPfindBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   );

/** returns the array of currently available branching rules */
extern
BRANCHRULE** SCIPgetBranchrules(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available branching rules */
extern
int SCIPgetNBranchrules(
   SCIP*            scip                /**< SCIP data structure */
   );

/** sets the priority of a branching rule */
extern
RETCODE SCIPsetBranchrulePriority(
   SCIP*            scip,               /**< SCIP data structure */
   BRANCHRULE*      branchrule,         /**< branching rule */
   int              priority            /**< new priority of the branching rule */
   );

/** sets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
extern
RETCODE SCIPsetBranchruleMaxdepth(
   SCIP*            scip,               /**< SCIP data structure */
   BRANCHRULE*      branchrule,         /**< branching rule */
   int              maxdepth            /**< new maxdepth of the branching rule */
   );

/** creates a display column and includes it in SCIP */
extern
RETCODE SCIPincludeDisp(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of display column */
   const char*      desc,               /**< description of display column */
   const char*      header,             /**< head line of display column */
   DISPSTATUS       dispstatus,         /**< display activation status of display column */
   DECL_DISPFREE    ((*dispfree)),      /**< destructor of display column */
   DECL_DISPINIT    ((*dispinit)),      /**< initialize display column */
   DECL_DISPEXIT    ((*dispexit)),      /**< deinitialize display column */
   DECL_DISPOUTPUT  ((*dispoutput)),    /**< output method */
   DISPDATA*        dispdata,           /**< display column data */
   int              width,              /**< width of display column (no. of chars used) */
   int              priority,           /**< priority of display column */
   int              position,           /**< relative position of display column */
   Bool             stripline           /**< should the column be separated with a line from its right neighbour? */
   );

/** returns the display column of the given name, or NULL if not existing */
extern
DISP* SCIPfindDisp(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   );

/** returns the array of currently available display columns */
extern
DISP** SCIPgetDisps(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the number of currently available display columns */
extern
int SCIPgetNDisps(
   SCIP*            scip                /**< SCIP data structure */
   );

/** automatically selects display columns for being shown w.r.t. the display width parameter */
extern
RETCODE SCIPautoselectDisps(
   SCIP*            scip                /**< SCIP data structure */
   );

/**@} */




/*
 * user interactive dialog methods
 */

/**@name User Interactive Dialog Methods */
/**@{ */

/** creates and captures a dialog */
extern
RETCODE SCIPcreateDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG**         dialog,             /**< pointer to store the dialog */
   DECL_DIALOGEXEC  ((*dialogexec)),    /**< execution method of dialog */
   DECL_DIALOGDESC  ((*dialogdesc)),    /**< description output method of dialog, or NULL */
   const char*      name,               /**< name of dialog: command name appearing in parent's dialog menu */
   const char*      desc,               /**< description of dialog used if description output method is NULL */
   Bool             issubmenu,          /**< is the dialog a submenu? */
   DIALOGDATA*      dialogdata          /**< user defined dialog data */
   );

/** captures a dialog */
extern
RETCODE SCIPcaptureDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          dialog              /**< dialog */
   );

/** releases a dialog */
extern
RETCODE SCIPreleaseDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG**         dialog              /**< pointer to the dialog */
   );

/** makes given dialog the root dialog of SCIP's interactive user shell; captures dialog and releases former root dialog */
extern
RETCODE SCIPsetRootDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          dialog              /**< dialog to be the root */
   );

/** returns the root dialog of SCIP's interactive user shell */
extern
DIALOG* SCIPgetRootDialog(
   SCIP*            scip                /**< SCIP data structure */
   );

/** adds a sub dialog to the given dialog as menu entry and captures it */
extern
RETCODE SCIPaddDialogEntry(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          dialog,             /**< dialog to extend, or NULL for root dialog */
   DIALOG*          subdialog           /**< subdialog to add as menu entry in dialog */
   );

/** starts interactive mode of SCIP by executing the root dialog */
extern
RETCODE SCIPstartInteraction(
   SCIP*            scip                /**< SCIP data structure */
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
RETCODE SCIPcreateProb(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< problem name */
   DECL_PROBDELORIG ((*probdelorig)),   /**< frees user data of original problem */
   DECL_PROBTRANS   ((*probtrans)),     /**< creates user data of transformed problem by transforming original user data */
   DECL_PROBDELTRANS((*probdeltrans)),  /**< frees user data of transformed problem */
   DECL_PROBINITSOL ((*probinitsol)),   /**< solving process initialization method of transformed data */
   DECL_PROBEXITSOL ((*probexitsol)),   /**< solving process deinitialization method of transformed data */
   PROBDATA*        probdata            /**< user problem data set by the reader */
   );

/** reads problem from file and initializes all solving data structures */
extern
RETCODE SCIPreadProb(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< problem file name */
   );

/** frees problem and branch and bound data structures */
extern
RETCODE SCIPfreeProb(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets user problem data */
extern
PROBDATA* SCIPgetProbData(
   SCIP*            scip                /**< SCIP data structure */
   );

/** sets user problem data */
extern
RETCODE SCIPsetProbData(
   SCIP*            scip,               /**< SCIP data structure */
   PROBDATA*        probdata            /**< user problem data to use */
   );

/** sets objective sense of problem */
extern
RETCODE SCIPsetObjsense(
   SCIP*            scip,               /**< SCIP data structure */
   OBJSENSE         objsense            /**< new objective sense */
   );

/** sets limit on objective function, such that only solutions better than this limit are accepted */
extern
RETCODE SCIPsetObjlimit(
   SCIP*            scip,               /**< SCIP data structure */
   Real             objlimit            /**< new primal objective limit */
   );

/** gets current limit on objective function */
extern
Real SCIPgetObjlimit(
   SCIP*            scip                /**< SCIP data structure */
   );

/** informs SCIP, that the objective value is always integral in every feasible solution */
extern
RETCODE SCIPsetObjIntegral(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns whether the objective value is known to be integral in every feasible solution */
extern
Bool SCIPisObjIntegral(
   SCIP*            scip                /**< SCIP data structure */
   );

/** adds variable to the problem */
extern
RETCODE SCIPaddVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to add */
   );

/** adds variable to the problem and uses it as pricing candidate to enter the LP */
extern
RETCODE SCIPaddPricedVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to add */
   Real             score               /**< pricing score of variable (the larger, the better the variable) */
   );

/** gets variables of the problem along with the numbers of different variable types; data may become invalid after
 *  calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 */
extern
RETCODE SCIPgetVarsData(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*             nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*             nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*             nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*             nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*             ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   );

/** gets array with active problem variables; data may become invalid after
 *  calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 */
extern
VAR** SCIPgetVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of active problem variables */
extern
int SCIPgetNVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of binary active problem variables */
extern
int SCIPgetNBinVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of integer active problem variables */
extern
int SCIPgetNIntVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of implicit integer active problem variables */
extern
int SCIPgetNImplVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of continuous active problem variables */
extern
int SCIPgetNContVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets variables of the original problem along with the numbers of different variable types; data may become invalid
 *  after a call to SCIPchgVarType()
 */
extern
RETCODE SCIPgetOrigVarsData(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*             nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*             nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*             nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*             nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*             ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   );

/** gets array with original problem variables; data may become invalid after
 *  a call to SCIPchgVarType()
 */
extern
VAR** SCIPgetOrigVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of original problem variables */
extern
int SCIPgetNOrigVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of binary original problem variables */
extern
int SCIPgetNOrigBinVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of integer original problem variables */
extern
int SCIPgetNOrigIntVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of implicit integer original problem variables */
extern
int SCIPgetNOrigImplVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of continuous original problem variables */
extern
int SCIPgetNOrigContVars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns variable of given name in the problem, or NULL if not existing */
extern
VAR* SCIPfindVar(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of variable to find */
   );

/** returns TRUE iff all potential variables exist in the problem, and FALSE, if there may be additional variables,
 *  that will be added in pricing and improve the objective value
 */
extern
Bool SCIPallVarsInProb(
   SCIP*            scip                /**< SCIP data structure */
   );

/** adds constraint to the problem; if constraint is only valid locally, it is added to the local subproblem of the
 *  active node (and all of its subnodes); otherwise it is added to the global problem;
 *  if a local constraint is added at the root node, it is automatically upgraded into a global constraint
 */
extern
RETCODE SCIPaddCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   );

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was added, or from the problem, if it was a problem constraint
 */
extern
RETCODE SCIPdelCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to delete */
   );

/** returns constraint of given name in the problem, or NULL if not existing */
extern
CONS* SCIPfindCons(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of constraint to find */
   );

/**@} */




/*
 * local subproblem methods
 */

/**@name Local Subproblem Methods */
/**@{ */

/** adds constraint to the given node (and all of its subnodes), even if it is a global constraint;
 *  if a local constraint is added to the root node, it is automatically upgraded into a global constraint
 */
extern
RETCODE SCIPaddConsNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to add constraint to */
   CONS*            cons                /**< constraint to add */
   );

/** adds constraint locally to the active node (and all of its subnodes), even if it is a global constraint;
 *  if a local constraint is added at the root node, it is automatically upgraded into a global constraint
 */
extern
RETCODE SCIPaddConsLocal(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   );

/** disables constraint's separation, enforcing, and propagation capabilities at the given node (and all subnodes) */
extern
RETCODE SCIPdisableConsNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to disable constraint in */
   CONS*            cons                /**< constraint to disable */
   );

/** disables constraint's separation, enforcing, and propagation capabilities at the active node (and all subnodes);
 *  if the method is called during problem modification or presolving, the constraint is globally deleted from the problem
 */
extern
RETCODE SCIPdisableConsLocal(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to disable */
   );

/** gets dual bound of active node */
extern
Real SCIPgetLocalDualbound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets lower bound of active node in transformed problem */
extern
Real SCIPgetLocalLowerbound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets dual bound of given node */
extern
Real SCIPgetNodeDualbound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node                /**< node to get dual bound for */
   );

/** gets lower bound of given node in transformed problem */
extern
Real SCIPgetNodeLowerbound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node                /**< node to get dual bound for */
   );

/** if given value is tighter (larger for minimization, smaller for maximization) than the active node's dual bound,
 *  sets the active node's dual bound to the new value
 */
extern
RETCODE SCIPupdateLocalDualbound(
   SCIP*            scip,               /**< SCIP data structure */
   Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   );

/** if given value is larger than the active node's lower bound (in transformed problem), sets the active node's
 *  lower bound to the new value
 */
extern
RETCODE SCIPupdateLocalLowerbound(
   SCIP*            scip,               /**< SCIP data structure */
   Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   );

/** if given value is tighter (larger for minimization, smaller for maximization) than the node's dual bound,
 *  sets the node's dual bound to the new value
 */
extern
RETCODE SCIPupdateNodeDualbound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to update dual bound for */
   Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   );

/** if given value is larger than the node's lower bound (in transformed problem), sets the node's lower bound
 *  to the new value
 */
extern
RETCODE SCIPupdateNodeLowerbound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to update lower bound for */
   Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   );

/**@} */




/*
 * solve methods
 */

/**@name Solve Methods */
/**@{ */

/** presolves problem */
extern
RETCODE SCIPpresolve(
   SCIP*            scip                /**< SCIP data structure */
   );

/** presolves and solves problem */
extern
RETCODE SCIPsolve(
   SCIP*            scip                /**< SCIP data structure */
   );

/** frees branch and bound tree and all solution process data; statistics, presolving data and transformed problem is
 *  preserved
 */
extern
RETCODE SCIPfreeSolve(
   SCIP*            scip                /**< SCIP data structure */
   );

/** frees all solution process data including presolving and transformed problem, only original problem is kept */
extern
RETCODE SCIPfreeTransform(
   SCIP*            scip                /**< SCIP data structure */
   );

/**@} */




/*
 * variable methods
 */

/**@name Variable Methods */
/**@{ */

/** creates and captures problem variable; if variable is of integral type, fractional bounds are automatically rounded; 
 *  an integer variable with bounds zero and one is automatically converted into a binary variable
 */
extern
RETCODE SCIPcreateVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var,                /**< pointer to variable object */
   const char*      name,               /**< name of variable, or NULL for automatic name creation */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             initial,            /**< should var's column be present in the initial root LP? */
   Bool             removeable,         /**< is var's column removeable from the LP (due to aging or cleanup)? */
   DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   VARDATA*         vardata             /**< user data for this specific variable */
   );

/** increases usage counter of variable */
extern
RETCODE SCIPcaptureVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to capture */
   );

/** decreases usage counter of variable, and frees memory if necessary */
extern
RETCODE SCIPreleaseVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var                 /**< pointer to variable */
   );

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 */
extern
RETCODE SCIPtransformVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get/create transformed variable for */
   VAR**            transvar            /**< pointer to store the transformed variable */
   );

/** gets and captures transformed variables for an array of variables;
 *  if a variable of the array is not yet transformed, a new transformed variable for this variable is created;
 *  it is possible to call this method with vars == transvars
 */
extern
RETCODE SCIPtransformVars(
   SCIP*            scip,               /**< SCIP data structure */
   int              nvars,              /**< number of variables to get/create transformed variables for */
   VAR**            vars,               /**< array with variables to get/create transformed variables for */
   VAR**            transvars           /**< array to store the transformed variables */
   );

/** gets corresponding transformed variable of a given variable;
 *  returns NULL as transvar, if transformed variable is not yet existing
 */
extern
RETCODE SCIPgetTransformedVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get transformed variable for */
   VAR**            transvar            /**< pointer to store the transformed variable */
   );

/** gets corresponding transformed variables for an array of variables;
 *  stores NULL in a transvars slot, if the transfored variable is not yet existing;
 *  it is possible to call this method with vars == transvars, but remember that variables that are not
 *  yet transformed will be replaced with NULL
 */
extern
RETCODE SCIPgetTransformedVars(
   SCIP*            scip,               /**< SCIP data structure */
   int              nvars,              /**< number of variables to get transformed variables for */
   VAR**            vars,               /**< array with variables to get transformed variables for */
   VAR**            transvars           /**< array to store the transformed variables */
   );

/** gets negated variable x' = lb + ub - x of variable x; negated variable is created, if not yet existing */
extern
RETCODE SCIPgetNegatedVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get negated variable for */
   VAR**            negvar              /**< pointer to store the negated variable */
   );

/** gets a binary variable that is equal to the given binary variable, and that is either active or the negated
 *  variable of an active binary variable; if the given variable is fixed, NULL is returned as representative,
 *  and *negated is TRUE iff the variable is fixed to TRUE
 */
extern
RETCODE SCIPgetBinvarRepresentative(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< binary variable to get binary representative for */
   VAR**            repvar,             /**< pointer to store the binary representative */
   Bool*            negated             /**< pointer to store whether the negation of an active variable was returned */
   );

/** gets solution value for variable in active node */
extern
Real SCIPgetVarSol(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get solution value for */
   );

/** gets solution values of multiple variables in active node */
extern
RETCODE SCIPgetVarSols(
   SCIP*            scip,               /**< SCIP data structure */
   int              nvars,              /**< number of variables to get solution value for */
   VAR**            vars,               /**< array with variables to get value for */
   Real*            vals                /**< array to store solution values of variables */
   );

/** gets strong branching information on COLUMN variable */
extern
RETCODE SCIPgetVarStrongbranch(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get strong branching values for */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching column down */
   Real*            up,                 /**< stores dual bound after branching column up */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   );

/** gets strong branching information on COLUMN variable of the last SCIPgetVarStrongbranch() call;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given variable;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 */
extern
RETCODE SCIPgetVarStrongbranchLast(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get last strong branching values for */
   Real*            down,               /**< stores dual bound after branching column down, or NULL */
   Real*            up,                 /**< stores dual bound after branching column up, or NULL */
   Real*            solval              /**< stores LP solution value of variable at last strong branching call, or NULL */
   );

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given column, or -1 if strong branching was never applied to the column in current run
 */
extern
Longint SCIPgetVarStrongbranchNode(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get last strong branching node for */
   );

/** changes variable's objective value */
extern
RETCODE SCIPchgVarObj(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the objective value for */
   Real             newobj              /**< new objective value */
   );

/** adds value to variable's objective value */
extern
RETCODE SCIPaddVarObj(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the objective value for */
   Real             addobj              /**< additional objective value */
   );

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) lower bound value;
 *  does not change the bounds of the variable
 */
extern
Real SCIPadjustedVarLb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to adjust the bound for */
   Real             lb                  /**< lower bound value to adjust */
   );

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) lower bound value;
 *  does not change the bounds of the variable
 */
extern
Real SCIPadjustedVarUb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to adjust the bound for */
   Real             ub                  /**< upper bound value to adjust */
   );

/** depending on SCIP's stage, changes lower bound of variable in the problem, in preprocessing, or in active node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 */
extern
RETCODE SCIPchgVarLb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

/** depending on SCIP's stage, changes upper bound of variable in the problem, in preprocessing, or in active node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 */
extern
RETCODE SCIPchgVarUb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

/** changes lower bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 */
extern
RETCODE SCIPchgVarLbNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

/** changes upper bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 */
extern
RETCODE SCIPchgVarUbNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

/** changes lower bound of variable in preprocessing or in the active node, if the new bound is tighter than the
 *  current bound; if possible, adjusts bound to integral value; doesn't store any inference information in the
 *  bound change, such that in conflict analysis, this change is treated like a branching decision
 */
extern
RETCODE SCIPtightenVarLb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new value for bound */
   Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** changes upper bound of variable in preprocessing or in the active node, if the new bound is tighter than the
 *  current bound; if possible, adjusts bound to integral value; doesn't store any inference information in the
 *  bound change, such that in conflict analysis, this change is treated like a branching decision
 */
extern
RETCODE SCIPtightenVarUb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new value for bound */
   Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** fixes binary variable to given value; in problem creation or preprocessing stage, the variable is converted
 *  into a fixed variable, and the given inference constraint is ignored; in solving stage, the variable is fixed
 *  locally at the given node, and the given inference constraint is stored, such that the conflict analysis is
 *  able to find out the reason for the deduction of the variable fixing
 */
extern
RETCODE SCIPinferBinVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< binary variable, that is deduced to a fixed value */
   Bool             fixedval,           /**< value to fix binary variable to */
   CONS*            infercons,          /**< constraint that deduced the fixing */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   );

/** informs variable x about a globally valid variable lower bound x >= b*z + d with integer variable z */
extern
RETCODE SCIPaddVarVlb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   VAR*             vlbvar,             /**< variable z    in x >= b*z + d */
   Real             vlbcoef,            /**< coefficient b in x >= b*z + d */
   Real             vlbconstant         /**< constant d    in x >= b*z + d */
   );

/** informs variable x about a globally valid variable upper bound x <= b*z + d with integer variable z */
extern
RETCODE SCIPaddVarVub(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   VAR*             vubvar,             /**< variable z    in x <= b*z + d */
   Real             vubcoef,            /**< coefficient b in x <= b*z + d */
   Real             vubconstant         /**< constant d    in x <= b*z + d */
   );

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
extern
RETCODE SCIPchgVarBranchFactor(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             branchfactor        /**< factor to weigh variable's branching score with */
   );

/** scales the branch factor of the variable with the given value */
extern
RETCODE SCIPscaleVarBranchFactor(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             scale               /**< factor to scale variable's branching factor with */
   );

/** adds the given value to the branch factor of the variable */
extern
RETCODE SCIPaddVarBranchFactor(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             addfactor           /**< value to add to the branch factor of the variable */
   );

/** sets the branch priority of the variable; variables with higher branch priority are always prefered to variables
 *  with lower priority in selection of branching variable
 */
extern
RETCODE SCIPchgVarBranchPriority(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   int              branchpriority      /**< branch priority of the variable */
   );

/** changes the branch priority of the variable to the given value, if it is larger than the current priority */
extern
RETCODE SCIPupdateVarBranchPriority(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   int              branchpriority      /**< new branch priority of the variable, if it is larger than current priority */
   );

/** adds the given value to the branch priority of the variable */
extern
RETCODE SCIPaddVarBranchPriority(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   int              addpriority         /**< value to add to the branch priority of the variable */
   );

/** changes type of variable in the problem; this changes the vars array returned from
 *  SCIPgetVars() and SCIPgetVarsData()
 */
extern
RETCODE SCIPchgVarType(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   VARTYPE          vartype             /**< new type of variable */
   );

/** in problem creation and solving stage, both bounds of the variable are set to the given value;
 *  in presolving stage, the variable is converted into a fixed variable, and bounds are changed respectively;
 *  conversion into a fixed variable changes the vars array returned from SCIPgetVars() and SCIPgetVarsData()
 */
extern
RETCODE SCIPfixVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to fix */
   Real             fixedval,           /**< value to fix variable at */
   Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   Bool*            fixed               /**< pointer to store whether the fixing was performed (variable was unfixed) */
   );

/** From a given equality a*x + b*y == c, aggregates one of the variables and removes it from the set of
 *  active problem variables. This changes the vars array returned from SCIPgetVars() and SCIPgetVarsData().
 *  In the first step, the equality is transformed into an equality with active problem variables
 *  a'*x' + b'*y' == c'. If x' == y', this leads to the detection of redundancy if a' == -b' and c' == 0,
 *  of infeasibility, if a' == -b' and c' != 0, or to a variable fixing x' == c'/(a'+b') (and possible
 *  infeasibility) otherwise.
 *  In the second step, the variable to be aggregated is chosen among x' and y', prefering a less strict variable
 *  type as aggregation variable (i.e. continuous variables are prefered over implicit integers, implicit integers
 *  over integers, and integers over binaries). If none of the variables is continuous, it is tried to find an integer
 *  aggregation (i.e. integral coefficients a'' and b'', such that a''*x' + b''*y' == c''). This can lead to
 *  the detection of infeasibility (e.g. if c'' is fractional), or to a rejection of the aggregation (denoted by
 *  aggregated == FALSE), if the resulting integer coefficients are too large and thus numerically instable.
 *
 *  The output flags have the following meaning:
 *  - infeasible: the problem is infeasible
 *  - redundant:  the equality can be deleted from the constraint set
 *  - aggregated: the aggregation was successfully performed (aggregated implies redundant)
 */
extern
RETCODE SCIPaggregateVars(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             varx,               /**< variable x in equality a*x + b*y == c */
   VAR*             vary,               /**< variable y in equality a*x + b*y == c */
   Real             scalarx,            /**< multiplier a in equality a*x + b*y == c */
   Real             scalary,            /**< multiplier b in equality a*x + b*y == c */
   Real             rhs,                /**< right hand side c in equality a*x + b*y == c */
   Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   Bool*            redundant,          /**< pointer to store whether the equality is (now) redundant */
   Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   );

/** converts variable into multi-aggregated variable; this changes the vars array returned from
 *  SCIPgetVars() and SCIPgetVarsData()
 */
extern
RETCODE SCIPmultiaggregateVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable x to aggregate */
   int              naggvars,           /**< number n of variables in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   VAR**            aggvars,            /**< variables y_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Real*            scalars,            /**< multipliers a_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Real             constant,           /**< constant shift c in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   );

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of "solvaldelta" in the
 *  variable's solution value and resulting change of "objdelta" in the in the LP's objective value
 */
extern
RETCODE SCIPupdateVarPseudocost(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   Real             weight              /**< weight in (0,1] of this update in pseudo cost sum */
   );

/** gets the variable's pseudo cost value for the given direction */
extern
Real SCIPgetVarPseudocost(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
extern
Real SCIPgetVarPseudocostCount(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   int              dir                 /**< branching direction: 0 (down), or 1 (up) */
   );

/** gets the variable's pseudo cost score value for the given LP solution value */
extern
Real SCIPgetVarPseudocostScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             solval              /**< variable's LP solution value */
   );

/** returns the average number of inferences found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 */
extern
Real SCIPgetVarAvgInferences(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** returns the variable's average inference score value */
extern
Real SCIPgetVarAvgInferenceScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< problem variable */
   );

/** returns the average number of cutoffs found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 */
extern
Real SCIPgetVarAvgCutoffs(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** returns the variable's average cutoff score value */
extern
Real SCIPgetVarAvgCutoffScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< problem variable */
   );

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor
 */
extern
Real SCIPgetVarAvgInferenceCutoffScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   );

/** gets user data for given variable */
extern
VARDATA* SCIPgetVarData(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< problem variable */
   );

/**@} */




/*
 * conflict analysis methods
 */

/**@name Conflict Analysis Methods */
/**@{ */

/** initializes the conflict analysis by clearing the conflict variable candidate queue */
extern
RETCODE SCIPinitConflictAnalysis(
   SCIP*            scip                /**< SCIP data structure */
   );

/** adds currently fixed binary variable to the conflict analysis' candidate storage; this method should be called in
 *  one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictVar() should be called for each variable,
 *      whose current assignment lead to the conflict (i.e. the infeasibility of a globally valid constraint).
 *   2. In the conflict variable resolution method of a constraint handler, SCIPaddConflictVar() should be called
 *      for each variable, whose current assignment lead to the deduction of the given conflict variable.
 */
extern
RETCODE SCIPaddConflictVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< conflict variable to add to conflict candidate queue */
   );

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set; the conflict analysis
 *  should only be called if a globally valid constraint was violated -- otherwise, the resulting conflict
 *  constraint wouldn't be globally valid
 */
extern
RETCODE SCIPanalyzeConflict(
   SCIP*            scip,               /**< SCIP data structure */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   );

/**@} */




/*
 * constraint methods
 */

/**@name Constraint Methods */
/**@{ */

/** creates and captures a constraint of the given constraint handler
 *  Warning! If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
extern
RETCODE SCIPcreateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to constraint */
   const char*      name,               /**< name of constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   );

/** increases usage counter of constraint */
extern
RETCODE SCIPcaptureCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to capture */
   );

/** decreases usage counter of constraint, and frees memory if necessary */
extern
RETCODE SCIPreleaseCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons                /**< pointer to constraint */
   );

/** gets and captures transformed constraint of a given constraint; if the constraint is not yet transformed,
 *  a new transformed constraint for this constraint is created
 */
extern
RETCODE SCIPtransformCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint to get/create transformed constraint for */
   CONS**           transcons           /**< pointer to store the transformed constraint */
   );

/** gets and captures transformed constraints for an array of constraints;
 *  if a constraint in the array is not yet transformed, a new transformed constraint for this constraint is created;
 *  it is possible to call this method with conss == transconss
 */
extern
RETCODE SCIPtransformConss(
   SCIP*            scip,               /**< SCIP data structure */
   int              nconss,             /**< number of constraints to get/create transformed constraints for */
   CONS**           conss,              /**< array with constraints to get/create transformed constraints for */
   CONS**           transconss          /**< array to store the transformed constraints */
   );

/** gets corresponding transformed constraint of a given constraint;
 *  returns NULL as transcons, if transformed constraint is not yet existing
 */
extern
RETCODE SCIPgetTransformedCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint to get the transformed constraint for */
   CONS**           transcons           /**< pointer to store the transformed constraint */
   );

/** gets corresponding transformed constraints for an array of constraints;
 *  stores NULL in a transconss slot, if the transformed constraint is not yet existing;
 *  it is possible to call this method with conss == transconss, but remember that constraints that are not
 *  yet transformed will be replaced with NULL
 */
extern
RETCODE SCIPgetTransformedConss(
   SCIP*            scip,               /**< SCIP data structure */
   int              nconss,             /**< number of constraints to get the transformed constraints for */
   CONS**           conss,              /**< constraints to get the transformed constraints for */
   CONS**           transconss          /**< array to store the transformed constraints */
   );

/** adds given value to age of constraint, but age can never become negative;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 */
extern
RETCODE SCIPaddConsAge(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint */
   Real             deltaage            /**< value to add to the constraint's age */
   );

/** increases age of constraint by 1.0;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 */
extern
RETCODE SCIPincConsAge(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint */
   );

/** resets age of constraint to zero;
 *  should be called
 *   - in constraint separation, if a cut was found for this constraint,
 *   - in constraint enforcing, if the constraint was violated, and
 *   - in constraint propagation, if a domain reduction was deduced;
 */
extern
RETCODE SCIPresetConsAge(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint */
   );

/** locks rounding of variables involved in the costraint */
extern
RETCODE SCIPlockConsVars(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   );

/** unlocks rounding of variables involved in the costraint */
extern
RETCODE SCIPunlockConsVars(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint */
   int              nunlockspos,        /**< decrease in number of rounding locks for constraint */
   int              nunlocksneg         /**< decrease in number of rounding locks for constraint's negation */
   );

/** checks single constraint for feasibility of the given solution */
extern
RETCODE SCIPcheckCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint to check */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** marks the constraint to be essential for feasibility */
extern
RETCODE SCIPsetConsChecked(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint */
   );

/**@} */




/*
 * LP methods
 */

/**@name LP Methods */
/**@{ */

/** checks, whether the LP was solved in the active node */
extern
Bool SCIPhasActNodeLP(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets solution status of current LP */
extern
LPSOLSTAT SCIPgetLPSolstat(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets objective value of current LP */
extern
Real SCIPgetLPObjval(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets pseudo objective value of the current LP */
extern
Real SCIPgetPseudoObjval(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets current LP columns along with the current number of LP columns */
extern
RETCODE SCIPgetLPColsData(
   SCIP*            scip,               /**< SCIP data structure */
   COL***           cols,               /**< pointer to store the array of LP columns, or NULL */
   int*             ncols               /**< pointer to store the number of LP columns, or NULL */
   );

/** gets current LP columns */
extern
COL** SCIPgetLPCols(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets current number of LP columns */
extern
int SCIPgetNLPCols(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets current LP rows along with the current number of LP rows */
extern
RETCODE SCIPgetLPRowsData(
   SCIP*            scip,               /**< SCIP data structure */
   ROW***           rows,               /**< pointer to store the array of LP rows, or NULL */
   int*             nrows               /**< pointer to store the number of LP rows, or NULL */
   );

/** gets current LP rows */
extern
ROW** SCIPgetLPRows(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets current number of LP rows */
extern
int SCIPgetNLPRows(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing
 */
extern
Bool SCIPallColsInLP(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1 */
extern
RETCODE SCIPgetLPBasisInd(
   SCIP*            scip,               /**< SCIP data structure */
   int*             basisind            /**< pointer to store the basis indices */
   );

/** gets a row from the inverse basis matrix B^-1 */
extern
RETCODE SCIPgetLPBInvRow(
   SCIP*            scip,               /**< SCIP data structure */
   int              r,                  /**< row number */
   Real*            coef                /**< pointer to store the coefficients of the row */
   );

/** gets a row from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A) */
extern
RETCODE SCIPgetLPBInvARow(
   SCIP*            scip,               /**< SCIP data structure */
   int              r,                  /**< row number */
   Real*            binvrow,            /**< row in B^-1 from prior call to SCIPgetLPBInvRow(), or NULL */
   Real*            coef                /**< pointer to store the coefficients of the row */
   );

/** calculates a weighted sum of all LP rows; for negative weights, the left and right hand side of the corresponding
 *  LP row are swapped in the summation
 */
extern
RETCODE SCIPsumLPRows(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            weights,            /**< row weights in row summation */
   REALARRAY*       sumcoef,            /**< array to store sum coefficients indexed by variables' probindex */
   Real*            sumlhs,             /**< pointer to store the left hand side of the row summation */
   Real*            sumrhs              /**< pointer to store the right hand side of the row summation */
   );

/* calculates a MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
extern
RETCODE SCIPcalcMIR(
   SCIP*            scip,               /**< SCIP data structure */
   Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   Bool*            success             /**< pointer to store whether the returned coefficients are a valid MIR cut */
   );

/** writes current LP to a file */
extern
RETCODE SCIPwriteLP(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      fname               /**< file name */
   );

/**@} */




/*
 * LP diving methods
 */

/**@name LP Diving Methods */
/**@{ */

/** initiates LP diving, making methods SCIPchgVarObjDive(), SCIPchgVarLbDive(), and SCIPchgVarUbDive() available */
extern
RETCODE SCIPstartDive(
   SCIP*            scip                /**< SCIP data structure */
   );

/** quits LP diving and resets bounds and objective values of columns to the current node's values */
extern
RETCODE SCIPendDive(
   SCIP*            scip                /**< SCIP data structure */
   );

/** changes variable's objective value in current dive */
extern
RETCODE SCIPchgVarObjDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the objective value for */
   Real             newobj              /**< new objective value */
   );

/** changes variable's lower bound in current dive */
extern
RETCODE SCIPchgVarLbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

/** changes variable's upper bound in current dive */
extern
RETCODE SCIPchgVarUbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

/** gets variable's objective value in current dive */
extern
Real SCIPgetVarObjDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get the bound for */
   );

/** gets variable's lower bound in current dive */
extern
Real SCIPgetVarLbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get the bound for */
   );

/** gets variable's upper bound in current dive */
extern
Real SCIPgetVarUbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get the bound for */
   );

/** solves the LP of the current dive */
extern
RETCODE SCIPsolveDiveLP(
   SCIP*            scip,               /**< SCIP data structure */
   int              itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   );

/** returns the number of the node in the current branch and bound run, where the last LP diving was applied */
extern
Longint SCIPgetLastDivenode(
   SCIP*            scip                /**< SCIP data structure */
   );

/**@} */




/*
 * LP row methods
 */

/**@name LP Row Methods */
/**@{ */

/** creates and captures an LP row */
extern
RETCODE SCIPcreateRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row,                /**< pointer to row */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            cols,               /**< array with columns of row entries */
   Real*            vals,               /**< array with coefficients of row entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             local,              /**< is row only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   );

/** creates and captures an LP row without any coefficients */
extern
RETCODE SCIPcreateEmptyRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row,                /**< pointer to row */
   const char*      name,               /**< name of row */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             local,              /**< is row only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   );

/** increases usage counter of LP row */
extern
RETCODE SCIPcaptureRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< row to capture */
   );

/** decreases usage counter of LP row, and frees memory if necessary */
extern
RETCODE SCIPreleaseRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row                 /**< pointer to LP row */
   );

/** changes left hand side of LP row */
extern
RETCODE SCIPchgRowLhs(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             lhs                 /**< new left hand side */
   );

/** changes right hand side of LP row */
extern
RETCODE SCIPchgRowRhs(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             rhs                 /**< new right hand side */
   );

/** informs row, that all subsequent additions of variables to the row should be cached and not directly applied;
 *  after all additions were applied, SCIPflushRowExtensions() must be called;
 *  while the caching of row extensions is activated, information methods of the row give invalid results;
 *  caching should be used, if a row is build with SCIPaddVarToRow() calls variable by variable to increase
 *  the performance
 */
extern
RETCODE SCIPcacheRowExtensions(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** flushes all cached row extensions after a call of SCIPcacheRowExtensions() and merges coefficients with
 *  equal columns into a single coefficient
 */
extern
RETCODE SCIPflushRowExtensions(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** resolves variable to columns and adds them with the coefficient to the row */
extern
RETCODE SCIPaddVarToRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   VAR*             var,                /**< problem variable */
   Real             val                 /**< value of coefficient */
   );

/** resolves variables to columns and adds them with the coefficients to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 */
extern
RETCODE SCIPaddVarsToRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   int              nvars,              /**< number of variables to add to the row */
   VAR**            vars,               /**< problem variables to add */
   Real*            vals                /**< values of coefficients */
   );

/** resolves variables to columns and adds them with the same single coefficient to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 */
extern
RETCODE SCIPaddVarsToRowSameCoef(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   int              nvars,              /**< number of variables to add to the row */
   VAR**            vars,               /**< problem variables to add */
   Real             val                 /**< unique value of all coefficients */
   );

/** tries to find a rational representation of the row and multiplies coefficients with common denominator */
extern
RETCODE SCIPmakeRowRational(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   Real             maxscale,           /**< maximal value to scale row with */
   Bool*            success             /**< stores whether row could be made rational */
   );

/** returns the minimal activity of a row w.r.t. the column's bounds */
extern
Real SCIPgetRowMinActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the maximal activity of a row w.r.t. the column's bounds */
extern
Real SCIPgetRowMaxActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** recalculates the activity of a row in the last LP solution */
extern
RETCODE SCIPrecalcRowLPActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the activity of a row in the last LP solution */
extern
Real SCIPgetRowLPActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the feasibility of a row in the last LP solution: negative value means infeasibility */
extern
Real SCIPgetRowLPFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** recalculates the activity of a row for the current pseudo solution */
extern
RETCODE SCIPrecalcRowPseudoActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the activity of a row for the current pseudo solution */
extern
Real SCIPgetRowPseudoActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the feasibility of a row for the current pseudo solution: negative value means infeasibility */
extern
Real SCIPgetRowPseudoFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** recalculates the activity of a row in the last LP or pseudo solution */
extern
RETCODE SCIPrecalcRowActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the activity of a row in the last LP or pseudo solution */
extern
Real SCIPgetRowActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the feasibility of a row in the last LP or pseudo solution */
extern
Real SCIPgetRowFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the activity of a row for the given primal solution */
extern
Real SCIPgetRowSolActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   SOL*             sol                 /**< primal CIP solution */
   );

/** returns the feasibility of a row for the given primal solution */
extern
Real SCIPgetRowSolFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   SOL*             sol                 /**< primal CIP solution */
   );

/** output row to file stream */
extern
RETCODE SCIPprintRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/**@} */




/*
 * cutting plane methods
 */

/**@name Cutting Plane Methods */
/**@{ */

/** adds cut to separation storage;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
extern
RETCODE SCIPaddCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             cut,                /**< separated cut */
   Real             score               /**< separation score of cut (the larger, the better the cut) */
   );

/** if not already existing, adds row to global cut pool */
extern
RETCODE SCIPaddPoolCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< cutting plane to add */
   );

/** removes the row from the global cut pool */
extern
RETCODE SCIPdelPoolCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< cutting plane to add */
   );

/** gets current number of rows in the global cut pool */
extern
int SCIPgetNPoolCuts(
   SCIP*            scip                /**< SCIP data structure */
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
RETCODE SCIPgetLPBranchCands(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands,           /**< pointer to store the number of LP branching candidates, or NULL */
   int*             npriolpcands        /**< pointer to store the number of candidates with maximal priority, or NULL */
   );

/** gets number of branching candidates for LP solution branching (number of fractional variables) */
extern
int SCIPgetNLPBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of branching candidates with maximal priority for LP solution branching */
extern
int SCIPgetNPrioLPBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets branching candidates for pseudo solution branching (nonfixed variables) along with the number of candidates */
extern
RETCODE SCIPgetPseudoBranchCands(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands,       /**< pointer to store the number of pseudo branching candidates, or NULL */
   int*             npriopseudocands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   );

/** gets branching candidates for pseudo solution branching (nonfixed variables) */
extern
int SCIPgetNPseudoBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPgetNPrioPseudoBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of binary branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPgetNPrioPseudoBranchBins(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of integer branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPgetNPrioPseudoBranchInts(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of implicit integer branching candidates with maximal branch priority for pseudo solution branching */
extern
int SCIPgetNPrioPseudoBranchImpls(
   SCIP*            scip                /**< SCIP data structure */
   );

/** calculates the branching score out of the gain predictions for a binary branching */
extern
Real SCIPgetBranchScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   Real             downgain,           /**< prediction of objective gain for rounding downwards */
   Real             upgain              /**< prediction of objective gain for rounding upwards */
   );

/** calculates the branching score out of the gain predictions for a branching with arbitrary many children */
extern
Real SCIPgetBranchScoreMultiple(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   int              nchildren,          /**< number of children that the branching will create */
   Real*            gains               /**< prediction of objective gain for each child */
   );

/** creates a child node of the active node */
extern
RETCODE SCIPcreateChild(
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           node,               /**< pointer to node data structure */
   Real             nodeselprio         /**< node selection priority of new node */
   );

/** branches on a variable; if solution value x' is fractional, two child nodes are created
 *  (x <= floor(x'), x >= ceil(x')), if solution value is integral, three child nodes are created
 *  (x <= x'-1, x == x', x >= x'+1)
 */
extern
RETCODE SCIPbranchVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to branch on */
   );

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 */
extern
RETCODE SCIPbranchLP(
   SCIP*            scip,               /**< SCIP data structure */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   );

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN */
extern
RETCODE SCIPbranchPseudo(
   SCIP*            scip,               /**< SCIP data structure */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   );

/**@} */




/*
 * primal solutions
 */

/**@name Primal Solution Methods */
/**@{ */

/** creates a primal solution, initialized to zero */
extern
RETCODE SCIPcreateSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current LP solution */
extern
RETCODE SCIPcreateLPSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current pseudo solution */
extern
RETCODE SCIPcreatePseudoSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the current solution */
extern
RETCODE SCIPcreateCurrentSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** frees primal CIP solution */
extern
RETCODE SCIPfreeSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol                 /**< pointer to the solution */
   );

/** links a primal solution to the current LP solution */
extern
RETCODE SCIPlinkLPSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the current pseudo solution */
extern
RETCODE SCIPlinkPseudoSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the current LP or pseudo solution */
extern
RETCODE SCIPlinkCurrentSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** clears a primal solution */
extern
RETCODE SCIPclearSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** stores solution values of variables in solution's own array */
extern
RETCODE SCIPunlinkSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** sets value of variable in primal CIP solution */
extern
RETCODE SCIPsetSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   );

/** sets values of multiple variables in primal CIP solution */
extern
RETCODE SCIPsetSolVals(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   int              nvars,              /**< number of variables to set solution value for */
   VAR**            vars,               /**< array with variables to add to solution */
   Real*            vals                /**< array with solution values of variables */
   );

/** increases value of variable in primal CIP solution */
extern
RETCODE SCIPincSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   VAR*             var,                /**< variable to increase solution value for */
   Real             incval              /**< increment for solution value of variable */
   );

/** returns value of variable in primal CIP solution, or in current LP/pseudo solution */
extern
Real SCIPgetSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   VAR*             var                 /**< variable to get value for */
   );

/** gets values of multiple variables in primal CIP solution */
extern
RETCODE SCIPgetSolVals(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int              nvars,              /**< number of variables to get solution value for */
   VAR**            vars,               /**< array with variables to get value for */
   Real*            vals                /**< array to store solution values of variables */
   );

/** returns objective value of primal CIP solution w.r.t. original problem, or current LP/pseudo objective value */
extern
Real SCIPgetSolOrigObj(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   );

/** returns transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value */
extern
Real SCIPgetSolTransObj(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   );

/** maps original space objective value into transformed objective value */
extern
Real SCIPtransformObj(
   SCIP*            scip,               /**< SCIP data structure */
   Real             obj                 /**< original space objective value to transform */
   );

/** maps transformed objective value into original space */
extern
Real SCIPretransformObj(
   SCIP*            scip,               /**< SCIP data structure */
   Real             obj                 /**< transformed objective value to retransform in original space */
   );

/** gets clock time, when this solution was found */
extern
Real SCIPgetSolTime(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** gets branch and bound run number, where this solution was found */
extern
int SCIPgetSolRunnum(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** gets node number of the specific branch and bound run, where this solution was found */
extern
Longint SCIPgetSolNodenum(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
extern
HEUR* SCIPgetSolHeur(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** outputs non-zero original variables of solution to file stream */
extern
RETCODE SCIPprintSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** outputs non-zero transformed variables of solution to file stream */
extern
RETCODE SCIPprintTransSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** gets best feasible primal solution found so far, or NULL if no solution has been found */
extern
SOL* SCIPgetBestSol(
   SCIP*            scip                /**< SCIP data structure */
   );

/** outputs best feasible primal solution found so far to file stream */
extern
RETCODE SCIPprintBestSol(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** outputs best feasible primal solution found so far in transformed variables to file stream */
extern
RETCODE SCIPprintBestTransSol(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** try to round given solution */
extern
RETCODE SCIProundSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   Bool*            success             /**< pointer to store whether rounding was successful */
   );

/** adds feasible primal solution to solution storage by copying it */
extern
RETCODE SCIPaddSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal CIP solution */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds primal solution to solution storage, frees the solution afterwards */
extern
RETCODE SCIPaddSolFree(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** checks solution for feasibility; if possible, adds it to storage by copying */
extern
RETCODE SCIPtrySol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards */
extern
RETCODE SCIPtrySolFree(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   );

/**@} */




/*
 * event methods
 */

/**@name Event Methods */
/**@{ */

/** catches a global (not variable dependent) event */
extern
RETCODE SCIPcatchEvent(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   );

/** drops a global event (stops to track event) */
extern
RETCODE SCIPdropEvent(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   );

/** catches an objective value or domain change event on the given variable */
extern
RETCODE SCIPcatchVarEvent(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to catch event for */
   EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   );

/** drops an objective value or domain change event (stops to track event) on the given variable */
extern
RETCODE SCIPdropVarEvent(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to drop event for */
   EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   );

/**@} */




/*
 * tree methods
 */

/**@name Tree Methods */
/**@{ */

/** gets children of active node along with the number of children */
extern
RETCODE SCIPgetChildren(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          children,           /**< pointer to store children array, or NULL if not needed */
   int*             nchildren           /**< pointer to store number of children, or NULL if not needed */
   );

/** gets number of children of active node */
extern
int SCIPgetNChildren(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets siblings of active node along with the number of siblings */
extern
RETCODE SCIPgetSiblings(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          siblings,           /**< pointer to store siblings array, or NULL if not needed */
   int*             nsiblings           /**< pointer to store number of siblings, or NULL if not needed */
   );

/** gets number of siblings of active node */
extern
int SCIPgetNSiblings(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets leaves of the tree along with the number of leaves */
extern
RETCODE SCIPgetLeaves(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          leaves,             /**< pointer to store leaves array, or NULL if not needed */
   int*             nleaves             /**< pointer to store number of leaves, or NULL if not needed */
   );

/** gets number of leaves in the tree */
extern
int SCIPgetNLeaves(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the best child of the active node w.r.t. the node selection priority assigned by the branching rule */
extern
NODE* SCIPgetPrioChild(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the best sibling of the active node w.r.t. the node selection priority assigned by the branching rule */
extern
NODE* SCIPgetPrioSibling(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the best child of the active node w.r.t. the node selection strategy */
extern
NODE* SCIPgetBestChild(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the best sibling of the active node w.r.t. the node selection strategy */
extern
NODE* SCIPgetBestSibling(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the best leaf from the node queue w.r.t. the node selection strategy */
extern
NODE* SCIPgetBestLeaf(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy */
extern
NODE* SCIPgetBestNode(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the node with smallest lower bound from the tree (child, sibling, or leaf) */
extern
NODE* SCIPgetBestboundNode(
   SCIP*            scip                /**< SCIP data structure */
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
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of processed nodes, including the active node */
extern
Longint SCIPgetNNodes(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of processed nodes in all runs, including the active node */
extern
Longint SCIPgetNTotalNodes(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of nodes left in the tree (children + siblings + leaves) */
extern
int SCIPgetNNodesLeft(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far */
extern
int SCIPgetNLPs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in primal and dual simplex */
extern
Longint SCIPgetNLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far for node relaxations */
extern
int SCIPgetNNodeLPs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far for node relaxations */
extern
Longint SCIPgetNNodeLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of LPs solved so far during diving */
extern
int SCIPgetNDivingLPs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far during diving */
extern
Longint SCIPgetNDivingLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of times, strong branching was called (each call represents solving two LPs) */
extern
int SCIPgetNStrongbranchs(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of simplex iterations used so far in strong branching */
extern
Longint SCIPgetNStrongbranchLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of pricing rounds performed so far at the current node */
extern
int SCIPgetNPriceRounds(
   SCIP*            scip                /**< SCIP data structure */
   );

/** get current number of variables in the pricing store */
extern
int SCIPgetNPricevars(
   SCIP*            scip                /**< SCIP data structure */
   );

/** get total number of pricing variables found so far */
extern
int SCIPgetNPricevarsFound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** get total number of pricing variables applied to the LPs */
extern
int SCIPgetNPricevarsApplied(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of separation rounds performed so far at the current node */
extern
int SCIPgetNSepaRounds(
   SCIP*            scip                /**< SCIP data structure */
   );

/** get current number of cuts in the cut store */
extern
int SCIPgetNCuts(
   SCIP*            scip                /**< SCIP data structure */
   );

/** get total number of cuts found so far */
extern
int SCIPgetNCutsFound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** get total number of cuts applied to the LPs */
extern
int SCIPgetNCutsApplied(
   SCIP*            scip                /**< SCIP data structure */
   );

/** get total number of conflict constraints found */
extern
Longint SCIPgetNConflictsFound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets depth of active node, or -1 if no active node exists */
extern
int SCIPgetDepth(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets maximal depth of all processed nodes in current branch and bound run */
extern
int SCIPgetMaxDepth(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets maximal depth of all processed nodes over all branch and bound runs */
extern
int SCIPgetMaxTotalDepth(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of backtracks, i.e. number of times, the new node was selected from the leaves queue */
extern
Longint SCIPgetNBacktracks(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets current plunging depth (succ. times, a child was selected as next node) */
extern
int SCIPgetPlungeDepth(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of active constraints at the current node */
extern
int SCIPgetNActiveConss(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of enabled constraints at the current node */
extern
int SCIPgetNEnabledConss(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets total number of globally valid constraints currently in the problem */
extern
int SCIPgetNGlobalConss(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets average dual bound of all unprocessed nodes */
extern
Real SCIPgetAvgDualbound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets average lower (dual) bound of all unprocessed nodes in transformed problem */
extern
Real SCIPgetAvgLowerbound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets global dual bound */
extern
Real SCIPgetDualbound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets global lower (dual) bound in transformed problem */
extern
Real SCIPgetLowerbound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets global primal bound (objective value of best solution or user objective limit) */
extern
Real SCIPgetPrimalbound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets global upper (primal) bound in transformed problem (objective value of best solution or user objective limit) */
extern
Real SCIPgetUpperbound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets global cutoff bound in transformed problem: a sub problem with lower bound larger than the cutoff
 *  cannot contain a better feasible solution; usually, this bound is equal to the upper bound, but if the
 *  objective value is always integral, the cutoff bound is (nearly) one less than the upper bound
 */
extern
Real SCIPgetCutoffbound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets current gap |(primalbound - dualbound)/dualbound| */
extern
Real SCIPgetGap(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets current gap |(upperbound - lowerbound)/lowerbound| in transformed problem */
extern
Real SCIPgetTransGap(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of feasible primal solutions found so far */
extern
Longint SCIPgetNSolsFound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of primal solutions stored in the solution storage */
extern
Longint SCIPgetNSols(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 */
extern
Bool SCIPisPrimalboundSol(
   SCIP*            scip                /**< SCIP data structure */
   );

/** outputs original problem to file stream */
extern
RETCODE SCIPprintOrigProblem(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** outputs transformed problem to file stream */
extern
RETCODE SCIPprintTransProblem(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** outputs SCIP status */
extern
RETCODE SCIPprintStatus(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** outputs solving statistics */
extern
RETCODE SCIPprintStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** outputs history statistics about branchings on variables */
extern
RETCODE SCIPprintBranchingStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/**@} */




/*
 * timing methods
 */

/**@name Timing Methods */
/**@{ */

/** gets current time of day in seconds (standard time zone) */
extern
Real SCIPgetTimeOfDay(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates a clock using the default clock type */
extern
RETCODE SCIPcreateClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   );

/** creates a clock counting the CPU user seconds */
extern
RETCODE SCIPcreateCPUClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   );

/** creates a clock counting the wall clock seconds */
extern
RETCODE SCIPcreateWallClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   );

/** frees a clock */
extern
RETCODE SCIPfreeClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   );

/** resets the time measurement of a clock to zero and completely stops the clock */
extern
RETCODE SCIPresetClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   );

/** starts the time measurement of a clock */
extern
RETCODE SCIPstartClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   );

/** stops the time measurement of a clock */
extern
RETCODE SCIPstopClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   );

/** gets the measured time of a clock in seconds */
extern
Real SCIPgetClockTime(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   );

/** sets the measured time of a clock to the given value in seconds */
extern
RETCODE SCIPsetClockTime(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock,              /**< clock timer */
   Real             sec                 /**< time in seconds to set the clock's timer to */
   );

/** gets the current total SCIP time in seconds */
extern
Real SCIPgetTotalTime(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the current solving time in seconds */
extern
Real SCIPgetSolvingTime(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the current presolving time in seconds */
extern
Real SCIPgetPresolvingTime(
   SCIP*            scip                /**< SCIP data structure */
   );

/**@} */




/*
 * numeric values and comparisons
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity */
extern
Real SCIPinfinity(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns value treated as zero */
extern
Real SCIPepsilon(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns value treated as zero for sums of floating point values */
extern
Real SCIPsumepsilon(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns feasibility tolerance for constraints */
extern
Real SCIPfeastol(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns feasibility tolerance for reduced costs */
extern
Real SCIPdualfeastol(
   SCIP*            scip                /**< SCIP data structure */
   );

/** sets the feasibility tolerance for constraints */
extern
RETCODE SCIPsetFeastol(
   SCIP*            scip,               /**< SCIP data structure */
   Real             feastol             /**< new feasibility tolerance for constraints */
   );

/** sets the feasibility tolerance for reduced costs */
extern
RETCODE SCIPsetDualfeastol(
   SCIP*            scip,               /**< SCIP data structure */
   Real             dualfeastol         /**< new feasibility tolerance for reduced costs */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** checks, if values are in range of epsilon */
extern
Bool SCIPisEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) lower than val2 */
extern
Bool SCIPisLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) greater than val2 */
extern
Bool SCIPisLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) greater than val2 */
extern
Bool SCIPisGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) lower than val2 */
extern
Bool SCIPisGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range epsilon of 0.0 */
extern
Bool SCIPisZero(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than epsilon */
extern
Bool SCIPisPositive(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -epsilon */
extern
Bool SCIPisNegative(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if values are in range of sumepsilon */
extern
Bool SCIPisSumEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) lower than val2 */
extern
Bool SCIPisSumLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
extern
Bool SCIPisSumLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) greater than val2 */
extern
Bool SCIPisSumGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
extern
Bool SCIPisSumGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range sumepsilon of 0.0 */
extern
Bool SCIPisSumZero(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than sumepsilon */
extern
Bool SCIPisSumPositive(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -sumepsilon */
extern
Bool SCIPisSumNegative(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if values are in range of feasibility tolerance */
extern
Bool SCIPisFeasEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than feasibility tolerance) lower than val2 */
extern
Bool SCIPisFeasLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than feasibility tolerance) greater than val2 */
extern
Bool SCIPisFeasLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than feasibility tolerance) greater than val2 */
extern
Bool SCIPisFeasGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than feasibility tolerance) lower than val2 */
extern
Bool SCIPisFeasGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range feasibility tolerance of 0.0 */
extern
Bool SCIPisFeasZero(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than feasibility tolerance */
extern
Bool SCIPisFeasPositive(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -feasibility tolerance */
extern
Bool SCIPisFeasNegative(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if the first given lower bound is tighter (w.r.t. bound strengthening epsilon) than the second one */
extern
Bool SCIPisLbBetter(
   SCIP*            scip,               /**< SCIP data structure */
   Real             lb1,                /**< first lower bound to compare */
   Real             lb2                 /**< second lower bound to compare */
   );

/** checks, if the first given upper bound is tighter (w.r.t. bound strengthening epsilon) than the second one */
extern
Bool SCIPisUbBetter(
   SCIP*            scip,               /**< SCIP data structure */
   Real             ub1,                /**< first upper bound to compare */
   Real             ub2                 /**< second upper bound to compare */
   );

/** checks, if the cut's activity is more then cutvioleps larger than the given right hand side;
 *  both, the activity and the rhs, should be normed
 */
extern
Bool SCIPisCutViolated(
   SCIP*            scip,               /**< SCIP data structure */
   Real             cutactivity,        /**< activity of the cut */
   Real             cutrhs              /**< right hand side value of the cut */
   );

/** checks, if relative difference of values is in range of epsilon */
extern
Bool SCIPisRelEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than epsilon */
extern
Bool SCIPisRelLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
extern
Bool SCIPisRelLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than epsilon */
extern
Bool SCIPisRelGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
extern
Bool SCIPisRelGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if rel. difference of values is in range of sumepsilon */
extern
Bool SCIPisSumRelEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if rel. difference of val1 and val2 is lower than sumepsilon */
extern
Bool SCIPisSumRelLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if rel. difference of val1 and val2 is not greater than sumepsilon */
extern
Bool SCIPisSumRelLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if rel. difference of val1 and val2 is greater than sumepsilon */
extern
Bool SCIPisSumRelGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if rel. difference of val1 and val2 is not lower than -sumepsilon */
extern
Bool SCIPisSumRelGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is (positive) infinite */
extern
Bool SCIPisInfinity(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against infinity */
   );

/** checks, if value is non-negative within the LP feasibility bounds */
extern
Bool SCIPisFeasible(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is integral within the LP feasibility bounds */
extern
Bool SCIPisIntegral(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if given fractional part is smaller than feastol */
extern
Bool SCIPisFracIntegral(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** rounds value + feasibility tolerance down to the next integer */
extern
Real SCIPfloor(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** rounds value - feasibility tolerance up to the next integer */
extern
Real SCIPceil(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** returns fractional part of value, i.e. x - floor(x) */
extern
Real SCIPfrac(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to return fractional part for */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPisEQ(scip, val1, val2)       SCIPsetIsEQ((scip)->set, val1, val2)       
#define SCIPisLT(scip, val1, val2)       SCIPsetIsLT((scip)->set, val1, val2)       
#define SCIPisLE(scip, val1, val2)       SCIPsetIsLE((scip)->set, val1, val2)       
#define SCIPisGT(scip, val1, val2)       SCIPsetIsGT((scip)->set, val1, val2)       
#define SCIPisGE(scip, val1, val2)       SCIPsetIsGE((scip)->set, val1, val2)       
#define SCIPisZero(scip, val)            SCIPsetIsZero((scip)->set, val)            
#define SCIPisPositive(scip, val)        SCIPsetIsPositive((scip)->set, val)        
#define SCIPisNegative(scip, val)        SCIPsetIsNegative((scip)->set, val)        
                                                                           
#define SCIPisSumEQ(scip, val1, val2)    SCIPsetIsSumEQ((scip)->set, val1, val2)    
#define SCIPisSumLT(scip, val1, val2)    SCIPsetIsSumLT((scip)->set, val1, val2)    
#define SCIPisSumLE(scip, val1, val2)    SCIPsetIsSumLE((scip)->set, val1, val2)    
#define SCIPisSumGT(scip, val1, val2)    SCIPsetIsSumGT((scip)->set, val1, val2)    
#define SCIPisSumGE(scip, val1, val2)    SCIPsetIsSumGE((scip)->set, val1, val2)    
#define SCIPisSumZero(scip, val)         SCIPsetIsSumZero((scip)->set, val)         
#define SCIPisSumPositive(scip, val)     SCIPsetIsSumPositive((scip)->set, val)     
#define SCIPisSumNegative(scip, val)     SCIPsetIsSumNegative((scip)->set, val)     
                                                                           
#define SCIPisFeasEQ(scip, val1, val2)   SCIPsetIsFeasEQ((scip)->set, val1, val2)   
#define SCIPisFeasLT(scip, val1, val2)   SCIPsetIsFeasLT((scip)->set, val1, val2)   
#define SCIPisFeasLE(scip, val1, val2)   SCIPsetIsFeasLE((scip)->set, val1, val2)   
#define SCIPisFeasGT(scip, val1, val2)   SCIPsetIsFeasGT((scip)->set, val1, val2)   
#define SCIPisFeasGE(scip, val1, val2)   SCIPsetIsFeasGE((scip)->set, val1, val2)   
#define SCIPisFeasZero(scip, val)        SCIPsetIsFeasZero((scip)->set, val)        
#define SCIPisFeasPositive(scip, val)    SCIPsetIsFeasPositive((scip)->set, val)    
#define SCIPisFeasNegative(scip, val)    SCIPsetIsFeasNegative((scip)->set, val)    
                                                                           
#define SCIPisLbBetter(scip, lb1, lb2)   SCIPsetIsLbBetter(scip->set, lb1, lb2)
#define SCIPisUbBetter(scip, ub1, ub2)   SCIPsetIsUbBetter(scip->set, ub1, ub2)
#define SCIPisCutViolated(scip, act,rhs) SCIPsetIsCutViolated((scip)->set, \
                                         (SCIPnodeGetDepth((scip)->tree->actnode) == 0), act, rhs)
                                                                           
#define SCIPisRelEQ(scip, val1, val2)    SCIPsetIsRelEQ((scip)->set, val1, val2)    
#define SCIPisRelLT(scip, val1, val2)    SCIPsetIsRelLT((scip)->set, val1, val2)    
#define SCIPisRelLE(scip, val1, val2)    SCIPsetIsRelLE((scip)->set, val1, val2)    
#define SCIPisRelGT(scip, val1, val2)    SCIPsetIsRelGT((scip)->set, val1, val2)    
#define SCIPisRelGE(scip, val1, val2)    SCIPsetIsRelGE((scip)->set, val1, val2)    
                                                                           
#define SCIPisSumRelEQ(scip, val1, val2) SCIPsetIsSumRelEQ((scip)->set, val1, val2) 
#define SCIPisSumRelLT(scip, val1, val2) SCIPsetIsSumRelLT((scip)->set, val1, val2) 
#define SCIPisSumRelLE(scip, val1, val2) SCIPsetIsSumRelLE((scip)->set, val1, val2) 
#define SCIPisSumRelGT(scip, val1, val2) SCIPsetIsSumRelGT((scip)->set, val1, val2) 
#define SCIPisSumRelGE(scip, val1, val2) SCIPsetIsSumRelGE((scip)->set, val1, val2) 
                                                                           
#define SCIPisInfinity(scip, val)        SCIPsetIsInfinity((scip)->set, val)        
#define SCIPisFeasible(scip, val)        SCIPsetIsFeasible((scip)->set, val)        
#define SCIPisIntegral(scip, val)        SCIPsetIsIntegral((scip)->set, val)        
#define SCIPisFracIntegral(scip, val)    SCIPsetIsFracIntegral((scip)->set, val)    

#define SCIPfloor(scip, val)             SCIPsetFloor((scip)->set, val)             
#define SCIPceil(scip, val)              SCIPsetCeil((scip)->set, val)              
#define SCIPfrac(scip, val)              SCIPsetFrac((scip)->set, val)              

#endif

/** outputs a real number, or "+infinity", or "-infinity" to a file */
extern
void SCIPprintReal(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val,                /**< value to print */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/**@} */




/*
 * memory management
 */

/**@name Memory Management */
/**@{ */

#define SCIPallocMemory(scip,ptr)               ( (allocMemory((ptr)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocMemoryArray(scip,ptr,num)      ( (allocMemoryArray((ptr), (num)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocMemorySize(scip,ptr,size)      ( (allocMemorySize((ptr), (size)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocMemoryArray(scip,ptr,newnum) ( (reallocMemoryArray((ptr), (newnum)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocMemorySize(scip,ptr,newsize) ( (reallocMemorySize((ptr), (newsize)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateMemory(scip, ptr, source)  ( (duplicateMemory((ptr), (source)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateMemoryArray(scip, ptr, source, num) \
                                                ( (duplicateMemoryArray((ptr), (source), (num)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPfreeMemory(scip,ptr)                freeMemory(ptr)
#define SCIPfreeMemoryNull(scip,ptr)            freeMemoryNull(ptr)
#define SCIPfreeMemoryArray(scip,ptr)           freeMemoryArray(ptr)
#define SCIPfreeMemoryArrayNull(scip,ptr)       freeMemoryArrayNull(ptr)
#define SCIPfreeMemorySize(scip,ptr)            freeMemorySize(ptr)
#define SCIPfreeMemorySizeNull(scip,ptr)        freeMemorySizeNull(ptr)

#define SCIPallocBlockMemory(scip,ptr)          ( (allocBlockMemory(SCIPmemhdr(scip), (ptr)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocBlockMemoryArray(scip,ptr,num) ( (allocBlockMemoryArray(SCIPmemhdr(scip), (ptr), (num)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocBlockMemorySize(scip,ptr,size) ( (allocBlockMemorySize(SCIPmemhdr(scip), (ptr), (size)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocBlockMemoryArray(scip,ptr,oldnum,newnum) \
                                                ( (reallocBlockMemoryArray(SCIPmemhdr(scip), (ptr), (oldnum), (newnum)) \
                                                  == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocBlockMemorySize(scip,ptr,oldsize,newsize) \
                                                ( (reallocBlockMemorySize(SCIPmemhdr(scip), (ptr), (oldsize), (newsize)) \
                                                  == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateBlockMemory(scip, ptr, source) \
                                                ( (duplicateBlockMemory(SCIPmemhdr(scip), (ptr), (source)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateBlockMemoryArray(scip, ptr, source, num) \
                                                ( (duplicateBlockMemoryArray(SCIPmemhdr(scip), (ptr), (source), (num)) \
                                                  == NULL) ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPensureBlockMemoryArray(scip,ptr,arraysizeptr,minsize) \
                                                ( (SCIPensureBlockMemoryArray_call((scip), (void**)(ptr), sizeof(**(ptr)), \
                                                   (arraysizeptr), (minsize))) )
#define SCIPfreeBlockMemory(scip,ptr)           freeBlockMemory(SCIPmemhdr(scip), (ptr))
#define SCIPfreeBlockMemoryNull(scip,ptr)       freeBlockMemoryNull(SCIPmemhdr(scip), (ptr))
#define SCIPfreeBlockMemoryArray(scip,ptr,num)  freeBlockMemoryArray(SCIPmemhdr(scip), (ptr), (num))
#define SCIPfreeBlockMemoryArrayNull(scip,ptr,num) \
                                                freeBlockMemoryArrayNull(SCIPmemhdr(scip), (ptr), (num))
#define SCIPfreeBlockMemorySize(scip,ptr,size)  freeBlockMemorySize(SCIPmemhdr(scip), (ptr), (size))
#define SCIPfreeBlockMemorySizeNull(scip,ptr,size) \
                                                freeBlockMemorySizeNull(SCIPmemhdr(scip), (ptr), (size))

#define SCIPallocBufferArray(scip,ptr,num)      SCIPallocBuffer(scip, (void**)(ptr), (num)*(int)sizeof(**(ptr)))
#define SCIPduplicateBufferArray(scip,ptr,source,num) \
                                                SCIPduplicateBuffer(scip, (void**)(ptr), source, (num)*(int)sizeof(**(ptr)))
#define SCIPreallocBufferArray(scip,ptr,num)    SCIPreallocBuffer(scip, (void**)(ptr), (num)*(int)sizeof(**(ptr)))
#define SCIPfreeBufferArray(scip,ptr)           SCIPfreeBuffer(scip, (void**)(ptr), 0*(int)sizeof(**(ptr)))
#define SCIPallocBufferSize(scip,ptr,size)      SCIPallocBuffer(scip, (void**)(ptr), size)
#define SCIPduplicateBufferSize(scip,ptr,source,size) \
                                                SCIPduplicateBuffer(scip, (void**)(ptr), source, size)
#define SCIPreallocBufferSize(scip,ptr,size)    SCIPreallocBuffer(scip, (void**)(ptr), size)
#define SCIPfreeBufferSize(scip,ptr)            SCIPfreeBuffer(scip, (void**)(ptr), 0*(int)sizeof(**(ptr)))

/** returns block memory to use at the current time */
extern
MEMHDR* SCIPmemhdr(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns the total number of bytes used in block memory */
extern
Longint SCIPgetMemUsed(
   SCIP*            scip                /**< SCIP data structure */
   );

/** calculate memory size for dynamically allocated arrays */
extern
int SCIPcalcMemGrowSize(
   SCIP*            scip,               /**< SCIP data structure */
   int              num                 /**< minimum number of entries to store */
   );

/** extends a dynamically allocated block memory array to be able to store at least the given number of elements;
 *  use SCIPensureBlockMemoryArray() define to call this method!
 */
extern
RETCODE SCIPensureBlockMemoryArray_call(
   SCIP*            scip,               /**< SCIP data structure */
   void**           arrayptr,           /**< pointer to dynamically sized array */
   size_t           elemsize,           /**< size in bytes of each element in array */
   int*             arraysize,          /**< pointer to current array size */
   int              minsize             /**< required minimal array size */
   );

/** gets a memory buffer with at least the given size */
extern
RETCODE SCIPallocBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   int              size                /**< required size in bytes of buffer */
   );

/** allocates a memory buffer with at least the given size and copies the given memory into the buffer */
extern
RETCODE SCIPduplicateBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   void*            source,             /**< memory block to copy into the buffer */
   int              size                /**< required size in bytes of buffer */
   );

/** reallocates a memory buffer to at least the given size */
extern
RETCODE SCIPreallocBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to the buffer */
   int              size                /**< required size in bytes of buffer */
   );

/** frees a memory buffer */
extern
RETCODE SCIPfreeBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to the buffer */
   int              dummysize           /**< used to get a safer define for SCIPfreeBufferSize/Array */
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
RETCODE SCIPcreateRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY**      realarray           /**< pointer to store the real array */
   );

/** frees a dynamic array of real values */
extern
RETCODE SCIPfreeRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY**      realarray           /**< pointer to the real array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPextendRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic real array */
extern
RETCODE SCIPclearRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray           /**< dynamic real array */
   );

/** gets value of entry in dynamic array */
extern
Real SCIPgetRealarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPsetRealarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx,                /**< array index to set value for */
   Real             val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
RETCODE SCIPincRealarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx,                /**< array index to increase value for */
   Real             incval              /**< value to increase array index */
   );

/** creates a dynamic array of int values */
extern
RETCODE SCIPcreateIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY**       intarray            /**< pointer to store the int array */
   );

/** frees a dynamic array of int values */
extern
RETCODE SCIPfreeIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY**       intarray            /**< pointer to the int array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPextendIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic int array */
extern
RETCODE SCIPclearIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray            /**< dynamic int array */
   );

/** gets value of entry in dynamic array */
extern
int SCIPgetIntarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPsetIntarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx,                /**< array index to set value for */
   int              val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
RETCODE SCIPincIntarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx,                /**< array index to increase value for */
   int              incval              /**< value to increase array index */
   );

/** creates a dynamic array of bool values */
extern
RETCODE SCIPcreateBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY**      boolarray           /**< pointer to store the bool array */
   );

/** frees a dynamic array of bool values */
extern
RETCODE SCIPfreeBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY**      boolarray           /**< pointer to the bool array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPextendBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic bool array */
extern
RETCODE SCIPclearBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** gets value of entry in dynamic array */
extern
Bool SCIPgetBoolarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPsetBoolarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx,                /**< array index to set value for */
   Bool             val                 /**< value to set array index to */
   );

/** creates a dynamic array of pointers */
extern
RETCODE SCIPcreatePtrarray(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY**       ptrarray            /**< pointer to store the int array */
   );

/** frees a dynamic array of pointers */
extern
RETCODE SCIPfreePtrarray(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY**       ptrarray            /**< pointer to the int array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPextendPtrarray(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY*        ptrarray,           /**< dynamic int array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic pointer array */
extern
RETCODE SCIPclearPtrarray(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY*        ptrarray            /**< dynamic int array */
   );

/** gets value of entry in dynamic array */
extern
void* SCIPgetPtrarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY*        ptrarray,           /**< dynamic int array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPsetPtrarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY*        ptrarray,           /**< dynamic int array */
   int              idx,                /**< array index to set value for */
   void*            val                 /**< value to set array index to */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPcreateRealarray(scip, realarray) SCIPrealarrayCreate(realarray, SCIPmemhdr(scip))
#define SCIPfreeRealarray(scip, realarray)   SCIPrealarrayFree(realarray)
#define SCIPextendRealarray(scip, realarray, minidx, maxidx) \
                                             SCIPrealarrayExtend(realarray, (scip)->set, minidx, maxidx)
#define SCIPclearRealarray(scip, realarray)  SCIPrealarrayClear(realarray)
#define SCIPgetRealarrayVal(scip, realarray, idx) \
                                             SCIPrealarrayGetVal(realarray, idx)
#define SCIPsetRealarrayVal(scip, realarray, idx, val) \
                                             SCIPrealarraySetVal(realarray, (scip)->set, idx, val)
#define SCIPincRealarrayVal(scip, realarray, idx, incval) \
                                             SCIPrealarrayIncVal(realarray, scip->set, idx, incval)

#define SCIPcreateIntarray(scip, intarray)   SCIPintarrayCreate(intarray, SCIPmemhdr(scip))
#define SCIPfreeIntarray(scip, intarray)     SCIPintarrayFree(intarray)
#define SCIPextendIntarray(scip, intarray, minidx, maxidx) \
                                             SCIPintarrayExtend(intarray, (scip)->set, minidx, maxidx)
#define SCIPclearIntarray(scip, intarray)    SCIPintarrayClear(intarray)
#define SCIPgetIntarrayVal(scip, intarray, idx) \
                                             SCIPintarrayGetVal(intarray, idx)
#define SCIPsetIntarrayVal(scip, intarray, idx, val) \
                                             SCIPintarraySetVal(intarray, (scip)->set, idx, val)
#define SCIPincIntarrayVal(scip, intarray, idx, incval) \
                                             SCIPintarrayIncVal(intarray, scip->set, idx, incval)

#define SCIPcreateBoolarray(scip, boolarray) SCIPboolarrayCreate(boolarray, SCIPmemhdr(scip))
#define SCIPfreeBoolarray(scip, boolarray)   SCIPboolarrayFree(boolarray)
#define SCIPextendBoolarray(scip, boolarray, minidx, maxidx) \
                                             SCIPboolarrayExtend(boolarray, (scip)->set, minidx, maxidx)
#define SCIPclearBoolarray(scip, boolarray)  SCIPboolarrayClear(boolarray)
#define SCIPgetBoolarrayVal(scip, boolarray, idx) \
                                             SCIPboolarrayGetVal(boolarray, idx)
#define SCIPsetBoolarrayVal(scip, boolarray, idx, val) \
                                             SCIPboolarraySetVal(boolarray, (scip)->set, idx, val)

#define SCIPcreatePtrarray(scip, ptrarray)   SCIPptrarrayCreate(ptrarray, SCIPmemhdr(scip))
#define SCIPfreePtrarray(scip, ptrarray)     SCIPptrarrayFree(ptrarray)
#define SCIPextendPtrarray(scip, ptrarray, minidx, maxidx) \
                                             SCIPptrarrayExtend(ptrarray, (scip)->set, minidx, maxidx)
#define SCIPclearPtrarray(scip, ptrarray)    SCIPptrarrayClear(ptrarray)
#define SCIPgetPtrarrayVal(scip, ptrarray, idx) \
                                             SCIPptrarrayGetVal(ptrarray, idx)
#define SCIPsetPtrarrayVal(scip, ptrarray, idx, val) \
                                             SCIPptrarraySetVal(ptrarray, (scip)->set, idx, val)

#endif

/**@} */




#ifndef NDEBUG

/*
 * debugging methods
 */

/**@name Debugging Methods */
/**@{ */

/** prints output about used memory */
extern
void SCIPdebugMemory(
   SCIP*            scip                /**< SCIP data structure */
   );

/**@} */

#endif




#endif
