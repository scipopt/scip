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



#define CHECK_SCIP(x) { RETCODE _retcode_; \
                        if( (_retcode_ = (x)) != SCIP_OKAY ) \
                          SCIPerror(stderr, _retcode_, __FILE__, __LINE__); \
                      }



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
#include "heur.h"
#include "sort.h"
#include "price.h"
#include "sepa.h"
#include "cutpool.h"
#include "primal.h"




/*
 * miscellaneous methods
 */

extern
Real SCIPversion(                       /**< returns scip version number */
   void
   );

extern
void SCIPprintVersion(                  /**< prints a version information line to a file stream */
   FILE*            file                /**< file stream to write version information, or NULL for stdout */
   );

extern
void SCIPerror(                         /**< prints error message and aborts program execution */
   FILE*            errout,             /**< file stream to write error message */
   RETCODE          retcode,            /**< SCIP return code causing the error */
   const char*      filename,           /**< source code file name */
   int              line                /**< source line */
   );




/*
 * general SCIP methods
 */

extern
RETCODE SCIPcreate(                     /**< creates and initializes SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   );

extern
RETCODE SCIPfree(                       /**< frees SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   );

extern
VERBLEVEL SCIPverbLevel(                /**< gets verbosity level for message output */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
RETCODE SCIPsetVerbLevel(               /**< sets verbosity level for message output */
   SCIP*            scip,               /**< SCIP data structure */
   VERBLEVEL        verblevel           /**< verbosity level for message output */
   );

extern
STAGE SCIPstage(                        /**< returns current stage of SCIP */
   SCIP*            scip                /**< SCIP data structure */
   );




/*
 * SCIP user functionality methods
 */

extern
RETCODE SCIPincludeReader(              /**< creates a reader and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
   DECL_READERFREE((*readerfree)),      /**< destructor of reader */
   DECL_READERINIT((*readerinit)),      /**< initialise reader */
   DECL_READEREXIT((*readerexit)),      /**< deinitialise reader */
   DECL_READERREAD((*readerread)),      /**< read method */
   READERDATA*      readerdata          /**< reader data */
   );

extern
RETCODE SCIPincludeConsHdlr(            /**< creates a constraint handler and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              chckpriority,       /**< priority of the constraint handler for checking infeasibility */
   int              propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE((*consfree)),          /**< destructor of constraint handler */
   DECL_CONSINIT((*consinit)),          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit)),          /**< deinitialise constraint handler */
   DECL_CONSDELE((*consdele)),          /**< free specific constraint data */
   DECL_CONSTRAN((*constran)),          /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA((*conssepa)),          /**< separate cutting planes */
   DECL_CONSENLP((*consenlp)),          /**< enforcing constraints for LP solutions */
   DECL_CONSENPS((*consenps)),          /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHCK((*conschck)),          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop)),          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   );

extern
RETCODE SCIPfindConsHdlr(               /**< finds the constraint handler of the given name */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint handler */
   CONSHDLR**       conshdlr            /**< pointer for storing the constraint handler (returns NULL, if not found) */
   );

extern
RETCODE SCIPincludeHeur(                /**< creates a primal heuristic and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of primal heuristic */
   const char*      desc,               /**< description of primal heuristic */
   int              priority,           /**< priority of the primal heuristic */
   int              freq,               /**< frequency for calling primal heuristic */
   DECL_HEURFREE((*heurfree)),          /**< destructor of primal heuristic */
   DECL_HEURINIT((*heurinit)),          /**< initialise primal heuristic */
   DECL_HEUREXIT((*heurexit)),          /**< deinitialise primal heuristic */
   DECL_HEUREXEC((*heurexec)),          /**< execution method of primal heuristic */
   HEURDATA*        heurdata            /**< primal heuristic data */
   );

extern
RETCODE SCIPfindHeur(                   /**< finds the primal heuristic of the given name */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of primal heuristic */
   HEUR**           heur                /**< pointer for storing the primal heuristic (returns NULL, if not found) */
   );

extern
RETCODE SCIPincludeEventhdlr(           /**< creates an event handler and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of event handler */
   const char*      desc,               /**< description of event handler */
   DECL_EVENTFREE((*eventfree)),        /**< destructor of event handler */
   DECL_EVENTINIT((*eventinit)),        /**< initialise event handler */
   DECL_EVENTEXIT((*eventexit)),        /**< deinitialise event handler */
   DECL_EVENTDELE((*eventdele)),        /**< free specific event data */
   DECL_EVENTEXEC((*eventexec)),        /**< execute event handler */
   EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   );

extern
RETCODE SCIPfindEventHdlr(              /**< finds the event handler of the given name */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of event handler */
   EVENTHDLR**      eventhdlr           /**< pointer for storing the event handler (returns NULL, if not found) */
   );

extern
RETCODE SCIPincludeNodesel(             /**< creates a node selector and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of node selector */
   const char*      desc,               /**< description of node selector */
   DECL_NODESELFREE((*nodeselfree)),    /**< destructor of node selector */
   DECL_NODESELINIT((*nodeselinit)),    /**< initialise node selector */
   DECL_NODESELEXIT((*nodeselexit)),    /**< deinitialise node selector */
   DECL_NODESELSLCT((*nodeselslct)),    /**< node selection method */
   DECL_NODESELCOMP((*nodeselcomp)),    /**< node comparison method */
   NODESELDATA*     nodeseldata,        /**< node selector data */
   Bool             lowestboundfirst    /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
   );

extern
RETCODE SCIPincludeBranchrule(          /**< creates a branching rule and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   DECL_BRANCHFREE((*branchfree)),      /**< destructor of branching rule */
   DECL_BRANCHINIT((*branchinit)),      /**< initialise branching rule */
   DECL_BRANCHEXIT((*branchexit)),      /**< deinitialise branching rule */
   DECL_BRANCHEXLP((*branchexlp)),      /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXPS((*branchexps)),      /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   );

extern
RETCODE SCIPincludeDisp(                /**< creates a display column and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of display column */
   const char*      desc,               /**< description of display column */
   const char*      header,             /**< head line of display column */
   DECL_DISPFREE((*dispfree)),          /**< destructor of display column */
   DECL_DISPINIT((*dispinit)),          /**< initialise display column */
   DECL_DISPEXIT((*dispexit)),          /**< deinitialise display column */
   DECL_DISPOUTP((*dispoutp)),          /**< output method */
   DISPDATA*        dispdata,           /**< display column data */
   int              width,              /**< width of display column (no. of chars used) */
   int              priority,           /**< priority of display column */
   int              position,           /**< relative position of display column */
   Bool             stripline           /**< should the column be separated with a line from its right neighbour? */
   );




/*
 * global problem methods
 */

extern
RETCODE SCIPcreateProb(                 /**< creates empty problem and initializes all solving data structures */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< problem name */
   );

extern
RETCODE SCIPreadProb(                   /**< reads problem from file and initializes all solving data structures */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< problem file name */
   );

extern
RETCODE SCIPfreeProb(                   /**< frees problem and branch-and-bound data structures */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
RETCODE SCIPsetObjsense(                /**< sets objective sense of problem */
   SCIP*            scip,               /**< SCIP data structure */
   OBJSENSE         objsense            /**< new objective sense */
   );

extern
RETCODE SCIPaddVar(                     /**< adds variable to the problem */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to add */
   );

extern
RETCODE SCIPgetVars(                    /**< gets variables of the problem */
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*             nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*             nbin,               /**< pointer to store number of binary variables or NULL if not needed */
   int*             nint,               /**< pointer to store number of integer variables or NULL if not needed */
   int*             nimpl,              /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*             ncont               /**< pointer to store number of continous variables or NULL if not needed */
   );

extern
RETCODE SCIPfindVar(                    /**< finds variable of given name in the problem */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of variable to find */
   VAR**            var                 /**< pointer to store the variable, returns NULL if not found */
   );

extern
RETCODE SCIPaddCons(                    /**< adds global constraint to the problem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   );

extern
RETCODE SCIPfindCons(                   /**< finds constraint of given name in the problem */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint to find */
   CONS**           cons                /**< pointer to store the constraint, returns NULL if not found */
   );




/*
 * local subproblem methods
 */

extern
RETCODE SCIPaddLocalCons(               /**< adds local constraint to the actual subproblem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   );




/*
 * solve methods
 */

extern
RETCODE SCIPsolve(                      /**< solves problem */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
RETCODE SCIPfreeSolve(                  /**< frees all solution process data, only original problem is kept */
   SCIP*            scip                /**< SCIP data structure */
   );




/*
 * variable methods
 */

extern
RETCODE SCIPcreateVar(                  /**< create and capture problem variable */
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var,                /**< pointer to variable object */
   const char*      name,               /**< name of column */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   );

extern
RETCODE SCIPcaptureVar(                 /**< increases usage counter of variable */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to capture */
   );

extern
RETCODE SCIPreleaseVar(                 /**< decreases usage counter of variable, and frees memory if necessary */
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var                 /**< pointer to variable */
   );

extern
RETCODE SCIPgetVarSol(                  /**< gets solution value for variable in active node */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get solution value for */
   Real*            solval              /**< pointer to store the solution value */
   );

extern
RETCODE SCIPchgVarLbNode(               /**< changes lower bound of variable in the given node; if possible, adjust bound
                                         *   to integral value */
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

extern
RETCODE SCIPchgVarUbNode(               /**< changes upper bound of variable in the given node; if possible, adjust bound
                                         *   to integral value */
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

extern
RETCODE SCIPchgVarLb(                   /**< depending on SCIP's stage, changes lower bound of variable in the problem,
                                         *   in preprocessing, or in active node; if possible, adjust bound to integral
                                         *   value */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

extern
RETCODE SCIPchgVarUb(                   /**< depending on SCIP's stage, changes upper bound of variable in the problem,
                                         *   in preprocessing, or in active node; if possible, adjust bound to integral
                                         *   value */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

extern
RETCODE SCIPchgVarType(                 /**< changes type of variable in the problem */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   VARTYPE          vartype             /**< new type of variable */
   );




/*
 * LP row methods
 */

extern
RETCODE SCIPcreateRow(                  /**< creates and captures an LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row,                /**< pointer to row */
   const char*      name,               /**< name of row */
   CONS*            cons,               /**< constraint, this row belongs to, or NULL if the row was separated from LP */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   );

extern
RETCODE SCIPcaptureRow(                 /**< increases usage counter of LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< row to capture */
   );

extern
RETCODE SCIPreleaseRow(                 /**< decreases usage counter of LP row, and frees memory if necessary */
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row                 /**< pointer to LP row */
   );

extern
RETCODE SCIPforbidRowRounding(          /**< forbids roundings of variables in row that may violate row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

extern
RETCODE SCIPallowRowRounding(           /**< allows roundings of variables in row that may violate row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

extern
RETCODE SCIPchgRowLhs(                  /**< changes left hand side of LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             lhs                 /**< new left hand side */
   );

extern
RETCODE SCIPchgRowRhs(                  /**< changes right hand side of LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             rhs                 /**< new right hand side */
   );

extern
RETCODE SCIPaddVarToRow(                /**< resolves variable to columns and adds them with the coefficient to the row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   VAR*             var,                /**< problem variable */
   Real             val                 /**< value of coefficient */
   );

extern
RETCODE SCIPgetRowActivityBounds(       /**< returns the minimal and maximal activity of a row w.r.t. the column's bounds */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            minactivity,        /**< pointer to store the minimal activity, or NULL */
   Real*            maxactivity         /**< pointer to store the maximal activity, or NULL */
   );

extern
RETCODE SCIPgetRowActivityResiduals(    /**< gets activity bounds for row after setting variable to zero */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   VAR*             var,                /**< variable to calculate activity residual for */
   Real             val,                /**< coefficient value of variable in linear constraint */
   Real*            minresactivity,     /**< pointer to store the minimal residual activity */
   Real*            maxresactivity      /**< pointer to store the maximal residual activity */
   );

extern
RETCODE SCIPinvalidRowActivityBounds(   /**< invalidates activity bounds, such that they are recalculated in next get */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   );

extern
RETCODE SCIPgetRowActivity(             /**< returns the activity of a row in the last LP solution */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            activity            /**< pointer to store the row's activity */
   );

extern
RETCODE SCIPgetRowFeasibility(          /**< returns the feasibility of a row in the last LP solution */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            feasibility         /**< pointer to store the row's feasibility */
   );

extern
RETCODE SCIPgetRowPseudoFeasibility(    /**< returns the feasibility of a row for the actual pseudo solution */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            pseudofeasibility   /**< pointer to store the row's pseudo feasibility */
   );

extern
RETCODE SCIPprintRow(                   /**< output row to file stream */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   FILE*            file                /**< output file (or NULL for standard output) */
   );




/*
 * LP methods
 */

extern
RETCODE SCIPhasActnodeLP(               /**< checks, whether the LP was solved in the active node */
   SCIP*            scip,               /**< SCIP data structure */
   Bool*            actnodehaslp        /**< pointer to store whether the active node has LP information */
   );

extern
RETCODE SCIPgetLPCols(                  /**< gets actual LP columns */
   SCIP*            scip,               /**< SCIP data structure */
   COL***           cols,               /**< pointer to store the array of LP columns, or NULL */
   int*             ncols               /**< pointer to store the number of LP columns, or NULL */
   );

extern
RETCODE SCIPgetLPRows(                  /**< gets actual LP rows */
   SCIP*            scip,               /**< SCIP data structure */
   ROW***           rows,               /**< pointer to store the array of LP rows, or NULL */
   int*             nrows               /**< pointer to store the number of LP rows, or NULL */
   );




/*
 * cutting plane methods
 */

extern
RETCODE SCIPaddCut(                     /**< adds cut to separation storage */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             cut,                /**< separated cut */
   Real             score               /**< separation score of cut (the larger, the better the cut) */
   );

extern
RETCODE SCIPpoolCut(                    /**< if not already existing, adds row to global cut pool */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< cutting plane to add */
   );

extern
RETCODE SCIPgetPoolsize(                /**< gets actual number of rows in the global cut pool */
   SCIP*            scip,               /**< SCIP data structure */
   int*             poolsize            /**< pointer to store the number of columns */
   );




/*
 * branching methods
 */

extern
RETCODE SCIPgetLPBranchCands(           /**< gets branching candidates for LP solution branching (fractional variables) */
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands            /**< pointer to store the number of LP branching candidates, or NULL */
   );

extern
RETCODE SCIPgetPseudoBranchCands(       /**< gets branching candidates for pseudo solution branching (nonfixed variables) */
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands        /**< pointer to store the number of pseudo branching candidates, or NULL */
   );

extern
RETCODE SCIPcreateChild(                /**< creates a child node of the active node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           node                /**< pointer to node data structure */
   );

extern
RETCODE SCIPbranchVar(                  /**< branches on a variable; if solution value x' is fractional, two child nodes
                                         *   are created (x <= floor(x'), x >= ceil(x')), if solution value is integral,
                                         *   three child nodes are created (x <= x'-1, x == x', x >= x'+1) */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to branch on */
   );

extern
RETCODE SCIPbranchLP(                   /**< calls branching rules to branch on an LP solution */
   SCIP*            scip,               /**< SCIP data structure */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   );




/*
 * constraint methods
 */

extern
RETCODE SCIPcreateCons(                 /**< creates and captures a constraint of the given constraint handler */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to constraint */
   const char*      name,               /**< name of constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             model               /**< is constraint necessary for feasibility? */
   );

extern
RETCODE SCIPcaptureCons(                /**< increases usage counter of constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to capture */
   );

extern
RETCODE SCIPreleaseCons(                /**< decreases usage counter of constraint, and frees memory if necessary */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons                /**< pointer to constraint */
   );




/*
 * primal solutions
 */

extern
RETCODE SCIPcreateSol(                  /**< creates a primal solution, initialized to zero */
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

extern
RETCODE SCIPcreateLPSol(                /**< creates a primal solution, initialized to the actual LP solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

extern
RETCODE SCIPcreatePseudoSol(            /**< creates a primal solution, initialized to the actual pseudo solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

extern
RETCODE SCIPcreateActSol(               /**< creates a primal solution, initialized to the actual solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

extern
RETCODE SCIPlinkLPSol(                  /**< links a primal solution to the actual LP solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

extern
RETCODE SCIPlinkPseudoSol(              /**< links a primal solution to the actual pseudo solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

extern
RETCODE SCIPlinkActSol(                 /**< links a primal solution to the actual LP or pseudo solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

extern
RETCODE SCIPclearSol(                   /**< clears a primal solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

extern
RETCODE SCIPunlinkSol(                  /**< stores solution values of variables in solution's own array */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   );

extern
RETCODE SCIPsetSolVal(                  /**< sets value of variable in primal CIP solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   );

extern
RETCODE SCIPincSolVal(                  /**< increases value of variable in primal CIP solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   VAR*             var,                /**< variable to increase solution value for */
   Real             incval              /**< increment for solution value of variable */
   );

extern
RETCODE SCIPgetSolVal(                  /**< returns value of variable in primal CIP solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   VAR*             var,                /**< variable to get value for */
   Real*            solval              /**< pointer to store the solution value */
   );

extern
RETCODE SCIPgetSolObj(                  /**< returns objective value of primal CIP solution */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   Real*            objval              /**< pointer to store the objective value */
   );

extern
RETCODE SCIPgetSolNodenum(              /**< gets node number, where this solution was found */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   Longint*         nodenum             /**< pointer to store the node number */
   );

extern
RETCODE SCIPgetSolHeur(                 /**< gets heuristic, that found this solution (or NULL if it's from the tree) */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   HEUR**           heur                /**< pointer to store the heuristic */
   );

extern
RETCODE SCIPprintSol(                   /**< outputs non-zero original variables of solution to file stream */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

extern
RETCODE SCIPprintSolTrans(              /**< outputs non-zero transformed variables of solution to file stream */
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

extern
RETCODE SCIPgetBestSol(                 /**< gets best feasible primal solution found so far */
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol                 /**< pointer to store the solution, returns NULL if no solution available */
   );

extern
RETCODE SCIPprintBestSol(               /**< outputs best feasible primal solution found so far to file stream */
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   );




/*
 * event methods
 */

extern
RETCODE SCIPcatchVarEvent(              /**< catches an event on the given variable */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to catch event for */
   EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   );

extern
RETCODE SCIPdropVarEvent(               /**< drops an event (stops to track event) on the given variable */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to drop event for */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   );




/*
 * tree methods
 */

extern
RETCODE SCIPgetChildren(                /**< gets children of active node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          children,           /**< pointer to store children array */
   int*             nchildren           /**< pointer to store number of children */
   );

extern
RETCODE SCIPgetSiblings(                /**< gets siblings of active node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          siblings,           /**< pointer to store siblings array */
   int*             nsiblings           /**< pointer to store number of siblings */
   );

extern
RETCODE SCIPgetBestChild(               /**< gets the best child of the active node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestchild           /**< pointer to store best child */
   );

extern
RETCODE SCIPgetBestSibling(             /**< gets the best sibling of the active node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestsibling         /**< pointer to store best sibling */
   );

extern
RETCODE SCIPgetBestLeaf(                /**< gets the best leaf from the node queue */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestleaf            /**< pointer to store best leaf */
   );

extern
RETCODE SCIPgetBestNode(                /**< gets the best node from the tree (child, sibling, or leaf) */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestnode            /**< pointer to store best leaf */
   );




/*
 * statistic methods
 */

extern
RETCODE SCIPgetNodenum(                 /**< gets number of processed nodes, including the active node */
   SCIP*            scip,               /**< SCIP data structure */
   Longint*         nodenum             /**< pointer to store the number of processed nodes */
   );

extern
RETCODE SCIPgetNNodesLeft(              /**< gets number of nodes left in the tree (children + siblings + leaves) */
   SCIP*            scip,               /**< SCIP data structure */
   int*             nnodes              /**< pointer to store the number of processed nodes */
   );

extern
RETCODE SCIPgetNLPIterations(           /**< gets total number of simplex iterations used so far */
   SCIP*            scip,               /**< SCIP data structure */
   int*             lpiterations        /**< pointer to store the iterations */
   );

extern
RETCODE SCIPgetActDepth(                /**< gets depth of active node */
   SCIP*            scip,               /**< SCIP data structure */
   int*             actdepth            /**< pointer to store the depth */
   );

extern
RETCODE SCIPgetMaxDepth(                /**< gets maximal depth of all processed nodes */
   SCIP*            scip,               /**< SCIP data structure */
   int*             maxdepth            /**< pointer to store the depth */
   );

extern
RETCODE SCIPgetPlungeDepth(             /**< gets actual plunging depth (succ. times, a child was selected as next node) */
   SCIP*            scip,               /**< SCIP data structure */
   int*             plungedepth         /**< pointer to store the depth */
   );

extern
RETCODE SCIPgetActDualBound(            /**< gets dual bound of active node */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            actdualbound        /**< pointer to store the dual bound */
   );

extern
RETCODE SCIPgetAvgDualBound(            /**< gets average dual bound of all unprocessed nodes */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            avgdualbound        /**< pointer to store the average dual bound */
   );

extern
RETCODE SCIPgetDualBound(               /**< gets global dual bound */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            dualbound           /**< pointer to store the dual bound */
   );

extern
RETCODE SCIPgetPrimalBound(             /**< gets global primal bound */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            primalbound         /**< pointer to store the primal bound */
   );




/*
 * numeric values and comparisons
 */

extern
Real SCIPinfinity(                      /**< returns value treated as infinity */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
Real SCIPepsilon(                       /**< returns value treated as zero */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
Real SCIPsumepsilon(                    /**< returns value treated as zero for sums of floating point values */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
Bool SCIPisEQ(                          /**< checks, if values are in range of epsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisL(                           /**< checks, if val1 is (more than epsilon) lower than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisLE(                          /**< checks, if val1 is not (more than epsilon) greater than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisG(                           /**< checks, if val1 is (more than epsilon) greater than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisGE(                          /**< checks, if val1 is not (more than epsilon) lower than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisZero(                        /**< checks, if value is in range epsilon of 0.0 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisPos(                         /**< checks, if value is greater than epsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisNeg(                         /**< checks, if value is lower than -epsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisSumEQ(                       /**< checks, if values are in range of sumepsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisSumL(                        /**< checks, if val1 is (more than sumepsilon) lower than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisSumLE(                       /**< checks, if val1 is not (more than sumepsilon) greater than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisSumG(                        /**< checks, if val1 is (more than sumepsilon) greater than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisSumGE(                       /**< checks, if val1 is not (more than sumepsilon) lower than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisSumZero(                     /**< checks, if value is in range sumepsilon of 0.0 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisSumPos(                      /**< checks, if value is greater than sumepsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisSumNeg(                      /**< checks, if value is lower than -sumepsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisRelEQ(                       /**< checks, if relative difference of values is in range of epsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisRelL(                        /**< checks, if relative difference of val1 and val2 is lower than epsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisRelLE(                       /**< checks, if relative difference of val1 and val2 is not greater than epsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisRelG(                        /**< checks, if relative difference of val1 and val2 is greater than epsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisRelGE(                       /**< checks, if relative difference of val1 and val2 is not lower than -epsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisSumRelEQ(                       /**< checks, if rel. difference of values is in range of sumepsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisSumRelL(                        /**< checks, if rel. difference of val1 and val2 is lower than sumepsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisSumRelLE(                       /**< checks, if rel. difference of val1 and val2 is not greater than sumepsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisSumRelG(                        /**< checks, if rel. difference of val1 and val2 is greater than sumepsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisSumRelGE(                       /**< checks, if rel. difference of val1 and val2 is not lower than -sumepsilon */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisInfinity(                    /**< checks, if value is (positive) infinite */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against infinity */
   );

extern
Real SCIPfloor(                         /**< rounds value down to the next integer */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Real SCIPceil(                          /**< rounds value up to the next integer */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisIntegral(                    /**< checks, if value is integral within the LP feasibility bounds */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisFeasible(                    /**< checks, if value is non-negative within the LP feasibility bounds */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );



/*
 * memory management
 */

#define SCIPallocBlockMemory(scip,ptr)          allocBlockMemory(SCIPmemhdr(scip), (ptr))
#define SCIPallocBlockMemoryArray(scip,ptr,num) allocBlockMemoryArray(SCIPmemhdr(scip), (ptr), (num))
#define SCIPallocBlockMemorySize(scip,ptr,size) allocBlockMemorySize(SCIPmemhdr(scip), (ptr), (size))
#define SCIPreallocBlockMemoryArray(scip,ptr,oldnum,newnum) \
                                                reallocBlockMemoryArray(SCIPmemhdr(scip), \
                                                  (ptr), (oldnum), (newnum))
#define SCIPreallocBlockMemorySize(scip,ptr,oldsize,newsize) \
                                                reallocBlockMemorySize(SCIPmemhdr(scip), \
                                                  (ptr), (oldsize), (newsize))
#define SCIPduplicateBlockMemory(scip, ptr, source) \
                                                duplicateBlockMemory(SCIPmemhdr(scip), (ptr), (source))
#define SCIPduplicateBlockMemoryArray(scip, ptr, source, num) \
                                                duplicateBlockMemoryArray(SCIPmemhdr(scip), (ptr), (source), (num))
#define SCIPfreeBlockMemory(scip,ptr)           freeBlockMemory(SCIPmemhdr(scip), (ptr))
#define SCIPfreeBlockMemoryNull(scip,ptr)       freeBlockMemoryNull(SCIPmemhdr(scip), (ptr))
#define SCIPfreeBlockMemoryArray(scip,ptr,num)  freeBlockMemoryArray(SCIPmemhdr(scip), (ptr), (num))
#define SCIPfreeBlockMemoryArrayNull(scip,ptr,num) \
                                                freeBlockMemoryArrayNull(SCIPmemhdr(scip), (ptr), (num))
#define SCIPfreeBlockMemorySize(scip,ptr,size)  freeBlockMemorySize(SCIPmemhdr(scip), (ptr), (size))
#define SCIPfreeBlockMemorySizeNull(scip,ptr,size) \
                                                freeBlockMemorySizeNull(SCIPmemhdr(scip), (ptr), (size))

#define SCIPcaptureBufferArray(scip,ptr,num)    SCIPcaptureBuffer(scip, (void**)(&(ptr)), (num)*sizeof(*(ptr)))
#define SCIPreleaseBufferArray(scip,ptr)        SCIPreleaseBuffer(scip, (void**)(&ptr))
#define SCIPcaptureBufferSize(scip,ptr,size)    SCIPcaptureBuffer(scip, (void**)(&(ptr)), size)
#define SCIPreleaseBufferSize(scip,ptr)         SCIPreleaseBuffer(scip, (void**)(&ptr))

extern
MEMHDR* SCIPmemhdr(                     /**< returns block memory to use at the current time */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
int SCIPcalcMemGrowSize(                /**< calculate memory size for dynamically allocated arrays */
   SCIP*            scip,               /**< SCIP data structure */
   int              num                 /**< minimum number of entries to store */
   );

extern
RETCODE SCIPcaptureBuffer(              /**< gets a memory buffer with at least the given size */
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   int              size                /**< required size in bytes of buffer */
   );

extern
RETCODE SCIPreleaseBuffer(              /**< releases a memory buffer */
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr                 /**< pointer to store the buffer */
   );




/*
 * dynamic arrays
 */

extern
RETCODE SCIPcreateRealarray(            /**< creates a dynamic array of real values */
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY**      realarray           /**< pointer to store the real array */
   );

extern
RETCODE SCIPfreeRealarray(              /**< frees a dynamic array of real values */
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY**      realarray           /**< pointer to the real array */
   );

extern
RETCODE SCIPextendRealarray(            /**< extends dynamic array to be able to store indices from minidx to maxidx */
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

extern
RETCODE SCIPclearRealarray(             /**< clears a dynamic real array */
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray           /**< dynamic real array */
   );

extern
Real SCIPgetRealarray(                  /**< gets value of entry in dynamic array */
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx                 /**< array index to get value for */
   );

extern
RETCODE SCIPsetRealarray(               /**< sets value of entry in dynamic array */
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx,                /**< array index to set value for */
   Real             val                 /**< value to set array index to */
   );

extern
RETCODE SCIPincRealarray(               /**< increases value of entry in dynamic array */
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx,                /**< array index to increase value for */
   Real             incval              /**< value to increase array index */
   );

extern
RETCODE SCIPcreateIntarray(             /**< creates a dynamic array of int values */
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY**       intarray            /**< pointer to store the int array */
   );

extern
RETCODE SCIPfreeIntarray(               /**< frees a dynamic array of int values */
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY**       intarray            /**< pointer to the int array */
   );

extern
RETCODE SCIPextendIntarray(             /**< extends dynamic array to be able to store indices from minidx to maxidx */
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

extern
RETCODE SCIPclearIntarray(              /**< clears a dynamic int array */
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray            /**< dynamic int array */
   );

extern
int SCIPgetIntarray(                    /**< gets value of entry in dynamic array */
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx                 /**< array index to get value for */
   );

extern
RETCODE SCIPsetIntarray(                /**< sets value of entry in dynamic array */
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx,                /**< array index to set value for */
   int              val                 /**< value to set array index to */
   );

extern
RETCODE SCIPincIntarray(                /**< increases value of entry in dynamic array */
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx,                /**< array index to increase value for */
   int              incval              /**< value to increase array index */
   );

extern
RETCODE SCIPcreateBoolarray(            /**< creates a dynamic array of bool values */
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY**      boolarray           /**< pointer to store the bool array */
   );

extern
RETCODE SCIPfreeBoolarray(              /**< frees a dynamic array of bool values */
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY**      boolarray           /**< pointer to the bool array */
   );

extern
RETCODE SCIPextendBoolarray(            /**< extends dynamic array to be able to store indices from minidx to maxidx */
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

extern
RETCODE SCIPclearBoolarray(             /**< clears a dynamic bool array */
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

extern
Bool SCIPgetBoolarray(                  /**< gets value of entry in dynamic array */
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx                 /**< array index to get value for */
   );

extern
RETCODE SCIPsetBoolarray(               /**< sets value of entry in dynamic array */
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx,                /**< array index to set value for */
   Bool             val                 /**< value to set array index to */
   );




/*
 * debugging methods
 */

#ifndef NDEBUG

extern
void SCIPdebugMemory(                   /**< prints output about used memory */
   SCIP*            scip                /**< SCIP data structure */
   );

#endif


#endif
