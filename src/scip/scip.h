/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip.h
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 */

/** Creating, capturing, releasing, and adding data objects.
 *
 *  Data objects (variables, constraints, rows) are subject to reference counting
 *  to avoid expensive copying operations. Creating such an object will set the
 *  reference count to one. Capturing an object increases the reference counter,
 *  releasing it decreases the counter. If the reference counter gets zero, the
 *  object is destroyed.
 *
 *  Remember that a created data object is automatically captured. If the user
 *  doesn't need the object anymore, he has to call the object's release() method.
 *
 *  When a data object is added to SCIP, it is captured again, such that a
 *  release() call does not destroy the object. If SCIP doesn't need the object
 *  anymore, it is automatically relased.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_H__
#define __SCIP_H__


/** SCIP operation stage */
enum Stage
{
   SCIP_STAGE_INIT       = 0,           /**< SCIP datastructures are initialized, no problem exists */
   SCIP_STAGE_PROBLEM    = 1,           /**< the problem is being created and modified */
   SCIP_STAGE_INITSOLVE  = 2,           /**< the solving process data is being initialized */
   SCIP_STAGE_PRESOLVING = 3,           /**< the problem is being presolved */
   SCIP_STAGE_SOLVING    = 4,           /**< the problem is being solved */
   SCIP_STAGE_SOLVED     = 5,           /**< the problem was solved */
   SCIP_STAGE_FREESOLVE  = 6            /**< the solving process data is being freed */
};
typedef enum Stage STAGE;


typedef struct Scip SCIP;               /**< SCIP main data structure */




#include <stdio.h>

#include "def.h"
#include "retcode.h"
#include "result.h"
#include "memory.h"
#include "message.h"
#include "reader.h"
#include "cons.h"
#include "var.h"
#include "lp.h"
#include "tree.h"
#include "nodesel.h"
#include "disp.h"
#include "branch.h"
#include "event.h"
#include "sepa.h"
#include "heur.h"
#include "misc.h"
#include "price.h"
#include "sepastore.h"
#include "cutpool.h"
#include "primal.h"




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
   const char*      msg                 /**< message to display */
   );

/** returns current stage of SCIP */
extern
STAGE SCIPstage(
   SCIP*            scip                /**< SCIP data structure */
   );

/**@} */




/*
 * parameter settings
 */

/** creates a Bool parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddBoolParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Bool             defaultvalue        /**< default value of the parameter */
   );

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddIntParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   int*             valueptr,           /**< pointer to store the current parameter value, or NULL */
   int              defaultvalue        /**< default value of the parameter */
   );

/** creates a Longint parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddLongintParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   Longint          defaultvalue        /**< default value of the parameter */
   );

/** creates a Real parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddRealParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Real             defaultvalue        /**< default value of the parameter */
   );

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddCharParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   char*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   char             defaultvalue        /**< default value of the parameter */
   );

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPaddStringParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   char**           valueptr,           /**< pointer to store the current parameter value, or NULL */
   const char*      defaultvalue        /**< default value of the parameter */
   );

/** gets the value of an existing Bool parameter */
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

/** creates a constraint handler and includes it in SCIP */
extern
RETCODE SCIPincludeConsHdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              chckpriority,       /**< priority of the constraint handler for checking infeasibility */
   int              sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int              propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   DECL_CONSINIT    ((*consinit)),      /**< initialise constraint handler */
   DECL_CONSEXIT    ((*consexit)),      /**< deinitialise constraint handler */
   DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA    ((*conssepa)),      /**< separate cutting planes */
   DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   );

/** returns the constraint handler of the given name, or NULL if not existing */
extern
CONSHDLR* SCIPfindConsHdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of constraint handler */
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
   DECL_SEPAINIT    ((*sepainit)),      /**< initialise separator */
   DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialise separator */
   DECL_SEPAEXEC    ((*sepaexec)),      /**< execution method of separator */
   SEPADATA*        sepadata            /**< separator data */
   );

/** returns the separator of the given name, or NULL if not existing */
extern
SEPA* SCIPfindSepa(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of separator */
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
   Bool             pseudonodes,        /**< call heuristic at nodes where only a pseudo solution exist? */
   DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   DECL_HEURINIT    ((*heurinit)),      /**< initialise primal heuristic */
   DECL_HEUREXIT    ((*heurexit)),      /**< deinitialise primal heuristic */
   DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   HEURDATA*        heurdata            /**< primal heuristic data */
   );

/** returns the primal heuristic of the given name, or NULL if not existing */
extern
HEUR* SCIPfindHeur(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of primal heuristic */
   );

/** creates an event handler and includes it in SCIP */
extern
RETCODE SCIPincludeEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of event handler */
   const char*      desc,               /**< description of event handler */
   DECL_EVENTFREE   ((*eventfree)),     /**< destructor of event handler */
   DECL_EVENTINIT   ((*eventinit)),     /**< initialise event handler */
   DECL_EVENTEXIT   ((*eventexit)),     /**< deinitialise event handler */
   DECL_EVENTDELETE ((*eventdelete)),   /**< free specific event data */
   DECL_EVENTEXEC   ((*eventexec)),     /**< execute event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   );

/** returns the event handler of the given name, or NULL if not existing */
extern
EVENTHDLR* SCIPfindEventHdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   );

/** creates a node selector and includes it in SCIP */
extern
RETCODE SCIPincludeNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of node selector */
   const char*      desc,               /**< description of node selector */
   DECL_NODESELFREE ((*nodeselfree)),   /**< destructor of node selector */
   DECL_NODESELINIT ((*nodeselinit)),   /**< initialise node selector */
   DECL_NODESELEXIT ((*nodeselexit)),   /**< deinitialise node selector */
   DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   NODESELDATA*     nodeseldata,        /**< node selector data */
   Bool             lowestboundfirst    /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
   );

/** returns the node selector of the given name, or NULL if not existing */
extern
NODESEL* SCIPfindNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   );

/** creates a branching rule and includes it in SCIP */
extern
RETCODE SCIPincludeBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   DECL_BRANCHINIT  ((*branchinit)),    /**< initialise branching rule */
   DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialise branching rule */
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

/** creates a display column and includes it in SCIP */
extern
RETCODE SCIPincludeDisp(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of display column */
   const char*      desc,               /**< description of display column */
   const char*      header,             /**< head line of display column */
   DECL_DISPFREE    ((*dispfree)),      /**< destructor of display column */
   DECL_DISPINIT    ((*dispinit)),      /**< initialise display column */
   DECL_DISPEXIT    ((*dispexit)),      /**< deinitialise display column */
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

/**@} */




/*
 * global problem methods
 */

/**@name Global Problem Methods */
/**@{ */

/** creates empty problem and initializes all solving data structures */
extern
RETCODE SCIPcreateProb(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< problem name */
   DECL_PROBDELETE  ((*probdelete)),    /**< frees user problem data */
   DECL_PROBTRANS   ((*probtrans)),     /**< transforms user problem data into data belonging to the transformed problem */
   PROBDATA*        probdata            /**< user problem data set by the reader */
   );

/** reads problem from file and initializes all solving data structures */
extern
RETCODE SCIPreadProb(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< problem file name */
   );

/** frees problem and branch-and-bound data structures */
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

/** adds variable to the problem */
extern
RETCODE SCIPaddVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to add */
   );

/** gets variables of the problem along with the numbers of different variable types */
extern
RETCODE SCIPgetVarsData(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*             nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*             nbin,               /**< pointer to store number of binary variables or NULL if not needed */
   int*             nint,               /**< pointer to store number of integer variables or NULL if not needed */
   int*             nimpl,              /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*             ncont               /**< pointer to store number of continous variables or NULL if not needed */
   );

/** gets array with active problem variables */
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

/** gets number of continous active problem variables */
extern
int SCIPgetNContVars(
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

/** adds global constraint to the problem */
extern
RETCODE SCIPaddCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   );

/** globally removes constraint from all subproblems; removes constraint from the subproblem of the node, where it
 *  was created, or from the global problem, if it was a globally valid problem constraint;
 *  the method must not be called for local check-constraint (i.e. constraints, that locally ensure feasibility);
 *  the constraint data is freed, and if the constraint is no longer used, it is freed completely
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

/** adds local constraint to the active node (and all of its subnodes) */
extern
RETCODE SCIPaddConsLocal(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   );

/** adds local constraint to the given node (and all of its subnodes) */
extern
RETCODE SCIPaddConsNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to add constraint to */
   CONS*            cons                /**< constraint to add */
   );

/** disables constraint's separation, enforcing, and propagation capabilities at the active node (and all subnodes);
 *  if the current node is the root node, or if the method is called during problem modification or presolving,
 *  the constraint is globally deleted from the problem
 */
extern
RETCODE SCIPdisableConsLocal(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to disable */
   );

/** disables constraint's separation, enforcing, and propagation capabilities at the given node (and all subnodes) */
extern
RETCODE SCIPdisableConsNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to disable constraint in */
   CONS*            cons                /**< constraint to disable */
   );

/**@} */




/*
 * solve methods
 */

/**@name Solve Methods */
/**@{ */

/** solves problem */
extern
RETCODE SCIPsolve(
   SCIP*            scip                /**< SCIP data structure */
   );

/** frees all solution process data, only original problem is kept */
extern
RETCODE SCIPfreeSolve(
   SCIP*            scip                /**< SCIP data structure */
   );

/**@} */




/*
 * variable methods
 */

/**@name Variable Methods */
/**@{ */

/** create and capture problem variable */
extern
RETCODE SCIPcreateVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var,                /**< pointer to variable object */
   const char*      name,               /**< name of column */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             removeable          /**< is var's column removeable from the LP (due to aging or cleanup)? */
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

/** gets solution value for variable in active node */
extern
Real SCIPgetVarSol(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get solution value for */
   );

/** gets strong branching information on COLUMN variable */
extern
RETCODE SCIPgetVarStrongbranch(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get solution value for */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching column down */
   Real*            up                  /**< stores dual bound after branching column up */
   );

/** changes lower bound of variable in the given node; if possible, adjust bound to integral value */
extern
RETCODE SCIPchgVarLbNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

/** changes upper bound of variable in the given node; if possible, adjust bound to integral value */
extern
RETCODE SCIPchgVarUbNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

/** depending on SCIP's stage, changes lower bound of variable in the problem, in preprocessing, or in active node;
 *  if possible, adjust bound to integral value
 */
extern
RETCODE SCIPchgVarLb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

/** depending on SCIP's stage, changes upper bound of variable in the problem, in preprocessing, or in active node;
 *  if possible, adjust bound to integral value
 */
extern
RETCODE SCIPchgVarUb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

/** changes type of variable in the problem */
extern
RETCODE SCIPchgVarType(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   VARTYPE          vartype             /**< new type of variable */
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
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate           /**< should the constraint be propagated during node processing? */
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

/** increases age of constraint; should be called in constraint separation, if no cut was found for this constraint,
 *  in constraint enforcing, if constraint was feasible, and in constraint propagation, if no domain reduction was
 *  deduced.
 */
extern
RETCODE SCIPincConsAge(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint */
   );

/** resets age of constraint to zero; should be called in constraint separation, if a cut was found for this constraint,
 *  in constraint enforcing, if the constraint was violated, and in constraint propagation, if a domain reduction was
 *  deduced.
 */
extern
RETCODE SCIPresetConsAge(
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
Bool SCIPhasActnodeLP(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets solution status of actual LP */
extern
LPSOLSTAT SCIPgetLPSolstat(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets objective value of actual LP, or SCIP_INVALID if LP is not solved yet */
extern
Real SCIPgetLPObjval(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets actual LP columns along with the actual number of LP columns */
extern
RETCODE SCIPgetLPColsData(
   SCIP*            scip,               /**< SCIP data structure */
   COL***           cols,               /**< pointer to store the array of LP columns, or NULL */
   int*             ncols               /**< pointer to store the number of LP columns, or NULL */
   );

/** gets actual LP columns */
extern
COL** SCIPgetLPCols(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets actual number of LP columns */
extern
int SCIPgetNLPCols(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets actual LP rows along with the actual number of LP rows */
extern
RETCODE SCIPgetLPRowsData(
   SCIP*            scip,               /**< SCIP data structure */
   ROW***           rows,               /**< pointer to store the array of LP rows, or NULL */
   int*             nrows               /**< pointer to store the number of LP rows, or NULL */
   );

/** gets actual LP rows */
extern
ROW** SCIPgetLPRows(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets actual number of LP rows */
extern
int SCIPgetNLPRows(
   SCIP*            scip                /**< SCIP data structure */
   );

/** returns TRUE iff all potential variables exist as columns in the LP, and FALSE, if there may be additional columns,
 *  that will be added in pricing and improve the objective value
 */
extern
Bool SCIPallVarsInLP(
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
   Real*            weights,            /**< row weights in row summation */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   Bool*            success             /**< pointer to store whether the returned coefficients are a valid MIR cut */
   );

/** writes actual LP to a file */
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

/** quits LP diving and resets bounds and objective values of columns to the actual node's values */
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
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
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

/** forbids roundings of variables in row that may violate row */
extern
RETCODE SCIPforbidRowRounding(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** allows roundings of variables in row that may violate row */
extern
RETCODE SCIPallowRowRounding(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
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

/** resolves variable to columns and adds them with the coefficient to the row */
extern
RETCODE SCIPaddVarToRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   VAR*             var,                /**< problem variable */
   Real             val                 /**< value of coefficient */
   );

/** tries to find a rational representation of the row and multiplies coefficients with common denominator */
extern
RETCODE SCIPmakeRowRational(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
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

/** returns the activity of a row in the last LP solution */
extern
Real SCIPgetRowLPActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the feasibility of a row in the last LP solution */
extern
Real SCIPgetRowLPFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the activity of a row for the actual pseudo solution */
extern
Real SCIPgetRowPseudoActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

/** returns the feasibility of a row for the actual pseudo solution */
extern
Real SCIPgetRowPseudoFeasibility(
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

/** adds cut to separation storage */
extern
RETCODE SCIPaddCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             cut,                /**< separated cut */
   Real             score               /**< separation score of cut (the larger, the better the cut) */
   );

/** if not already existing, adds row to global cut pool */
extern
RETCODE SCIPpoolCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< cutting plane to add */
   );

/** gets actual number of rows in the global cut pool */
extern
int SCIPgetPoolsize(
   SCIP*            scip                /**< SCIP data structure */
   );

/**@} */




/*
 * branching methods
 */

/**@name Branching Methods */
/**@{ */

/** gets branching candidates for LP solution branching (fractional variables) along with solution values,
 *  fractionalities, and number of branching candidates
 */
extern
RETCODE SCIPgetLPBranchCands(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands            /**< pointer to store the number of LP branching candidates, or NULL */
   );

/** gets number of branching candidates for LP solution branching (number of fractional variables) */
extern
int SCIPgetNLPBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets branching candidates for pseudo solution branching (nonfixed variables) along with the number of candidates */
extern
RETCODE SCIPgetPseudoBranchCands(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands        /**< pointer to store the number of pseudo branching candidates, or NULL */
   );

/** gets branching candidates for pseudo solution branching (nonfixed variables) */
extern
int SCIPgetNPseudoBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   );

/** calculates the branching score out of the downward and upward gain prediction */
extern
Real SCIPgetBranchScore(
   SCIP*            scip,               /**< SCIP data structure */
   Real             downgain,           /**< prediction of objective gain for branching downwards */
   Real             upgain              /**< prediction of objective gain for branching upwards */
   );

/** creates a child node of the active node */
extern
RETCODE SCIPcreateChild(
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           node                /**< pointer to node data structure */
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

/** calls branching rules to branch on an LP solution */
extern
RETCODE SCIPbranchLP(
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

/** creates a primal solution, initialized to the actual LP solution */
extern
RETCODE SCIPcreateLPSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the actual pseudo solution */
extern
RETCODE SCIPcreatePseudoSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a primal solution, initialized to the actual solution */
extern
RETCODE SCIPcreateActSol(
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

/** links a primal solution to the actual LP solution */
extern
RETCODE SCIPlinkLPSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the actual pseudo solution */
extern
RETCODE SCIPlinkPseudoSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

/** links a primal solution to the actual LP or pseudo solution */
extern
RETCODE SCIPlinkActSol(
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

/** increases value of variable in primal CIP solution */
extern
RETCODE SCIPincSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   VAR*             var,                /**< variable to increase solution value for */
   Real             incval              /**< increment for solution value of variable */
   );

/** returns value of variable in primal CIP solution, or in actual LP/pseudo solution */
extern
Real SCIPgetSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution, or NULL for actual LP/pseudo solution */
   VAR*             var                 /**< variable to get value for */
   );

/** returns objective value of primal CIP solution, or actual LP/pseudo objective value */
extern
Real SCIPgetSolObj(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution, or NULL for actual LP/pseudo objective value */
   );

/** returns transformed objective value of primal CIP solution, or transformed actual LP/pseudo objective value */
extern
Real SCIPgetSolTransObj(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution, or NULL for actual LP/pseudo objective value */
   );

/** maps transformed objective value into original space */
extern
Real SCIPretransformObj(
   SCIP*            scip,               /**< SCIP data structure */
   Real             obj                 /**< transformed objective value to retransform in original space */
   );

/** gets node number, where this solution was found */
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

/** adds feasible primal solution to solution storage by copying it */
extern
RETCODE SCIPaddSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal CIP solution */
   );

/** adds primal solution to solution storage, frees the solution afterwards */
extern
RETCODE SCIPaddSolFree(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol                 /**< pointer to primal CIP solution; is cleared in function call */
   );

/** checks solution for feasibility; if possible, adds it to storage by copying */
extern
RETCODE SCIPtrySol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal CIP solution */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards */
extern
RETCODE SCIPtrySolFree(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
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
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   );

/** catches a domain change event on the given variable */
extern
RETCODE SCIPcatchVarEvent(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to catch event for */
   EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   );

/** drops a domain change event (stops to track event) on the given variable */
extern
RETCODE SCIPdropVarEvent(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to drop event for */
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

/** gets the best child of the active node */
extern
NODE* SCIPgetBestChild(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the best sibling of the active node */
extern
NODE* SCIPgetBestSibling(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the best leaf from the node queue */
extern
NODE* SCIPgetBestLeaf(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets the best node from the tree (child, sibling, or leaf) */
extern
NODE* SCIPgetBestNode(
   SCIP*            scip                /**< SCIP data structure */
   );

/** if given value is larger than the node's lower bound, sets the node's lower bound to the new value */
extern
RETCODE SCIPupdateNodeLowerBound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to update lower bound for */
   Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   );

/**@} */




/*
 * statistic methods
 */

/**@name Statistic Methods */
/**@{ */

/** gets number of processed nodes, including the active node */
extern
Longint SCIPgetNodenum(
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

/** gets total number of simplex iterations used so far */
extern
int SCIPgetNLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets number of separation rounds performed so far at the current node */
extern
int SCIPgetNSepaRounds(
   SCIP*            scip                /**< SCIP data structure */
   );

/** get actual number of cuts in the cut store */
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

/** gets depth of active node, or -1 if no active node exists */
extern
int SCIPgetActDepth(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets maximal depth of all processed nodes */
extern
int SCIPgetMaxDepth(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets actual plunging depth (succ. times, a child was selected as next node) */
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

/** gets dual bound of active node */
extern
Real SCIPgetActDualBound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets lower (dual) bound of active node in transformed problem */
extern
Real SCIPgetActTransLowerBound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets average dual bound of all unprocessed nodes */
extern
Real SCIPgetAvgDualBound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets average lower (dual) bound of all unprocessed nodes in transformed problem */
extern
Real SCIPgetAvgTransLowerBound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets global dual bound */
extern
Real SCIPgetDualBound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets global lower (dual) bound in transformed problem */
extern
Real SCIPgetTransLowerBound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets global primal bound */
extern
Real SCIPgetPrimalBound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** gets global upper (primal) bound in transformed problem */
extern
Real SCIPgetTransUpperBound(
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
int SCIPgetNSolsFound(
   SCIP*            scip                /**< SCIP data structure */
   );

/** outputs solving statistics */
extern
RETCODE SCIPprintStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
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

/** returns feasibility tolerance */
extern
Real SCIPfeastol(
   SCIP*            scip                /**< SCIP data structure */
   );

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

/** checks, if value is integral within the LP feasibility bounds */
extern
Bool SCIPisIntegral(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is non-negative within the LP feasibility bounds */
extern
Bool SCIPisFeasible(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
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
#define SCIPfreeMemory(scip,ptr)                freeMemory((ptr))
#define SCIPfreeMemoryNull(scip,ptr)            freeMemoryNull((ptr))
#define SCIPfreeMemoryArray(scip,ptr)           freeMemoryArray((ptr))
#define SCIPfreeMemoryArrayNull(scip,ptr)       freeMemoryArrayNull((ptr))
#define SCIPfreeMemorySize(scip,ptr)            freeMemorySize((ptr))
#define SCIPfreeMemorySizeNull(scip,ptr)        freeMemorySizeNull((ptr))

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
#define SCIPfreeBlockMemory(scip,ptr)           freeBlockMemory(SCIPmemhdr(scip), (ptr))
#define SCIPfreeBlockMemoryNull(scip,ptr)       freeBlockMemoryNull(SCIPmemhdr(scip), (ptr))
#define SCIPfreeBlockMemoryArray(scip,ptr,num)  freeBlockMemoryArray(SCIPmemhdr(scip), (ptr), (num))
#define SCIPfreeBlockMemoryArrayNull(scip,ptr,num) \
                                                freeBlockMemoryArrayNull(SCIPmemhdr(scip), (ptr), (num))
#define SCIPfreeBlockMemorySize(scip,ptr,size)  freeBlockMemorySize(SCIPmemhdr(scip), (ptr), (size))
#define SCIPfreeBlockMemorySizeNull(scip,ptr,size) \
                                                freeBlockMemorySizeNull(SCIPmemhdr(scip), (ptr), (size))

#define SCIPcaptureBufferArray(scip,ptr,num)    SCIPcaptureBuffer(scip, (void**)(ptr), (num)*sizeof(**(ptr)))
#define SCIPreleaseBufferArray(scip,ptr)        SCIPreleaseBuffer(scip, (void**)(ptr), 0*sizeof(**(ptr)))
#define SCIPcaptureBufferSize(scip,ptr,size)    SCIPcaptureBuffer(scip, (void**)(ptr), size)
#define SCIPreleaseBufferSize(scip,ptr)         SCIPreleaseBuffer(scip, (void**)(ptr), 0*sizeof(**(ptr)))

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

/** gets a memory buffer with at least the given size */
extern
RETCODE SCIPcaptureBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   int              size                /**< required size in bytes of buffer */
   );

/** releases a memory buffer */
extern
RETCODE SCIPreleaseBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   int              dummysize           /**< used to get a safer define for SCIPreleaseBufferSize/Array */
   );

/**@} */




/*
 * dynamic arrays
 */

/**@name Dynamic Arrays */
/**@{ */

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
