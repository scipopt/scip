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
 *  If the user wants to add a data object (variable, constraint, row) to SCIP
 *  without using it for himself, he has to call the object's create() method
 *  and the corresponding add() method. He doesn't have to take care of the
 *  destruction, because SCIP will free the object if no longer needed. The
 *  user should not use the object after it was added to SCIP, because he cannot
 *  assume, that the object is still living.
 *
 *  If the user wants to use the data object for himself, he has to call the
 *  object's create() method and the capture() method to tell SCIP, the object
 *  is in use by the user. The user may add the object to SCIP with a call to
 *  the add() method. If the user doesn't need the object any more, he must
 *  call the release() method to tell SCIP, that the object may be destroyed,
 *  if it is no longer in use.
 *  
 *  Don't call the release() method of an object, if you didn't call the capture()
 *  method before, because this can destroy an object used by SCIP.
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
#include "memory.h"
#include "message.h"
#include "reader.h"
#include "cons.h"
#include "var.h"
#include "lp.h"
#include "tree.h"
#include "nodesel.h"
#include "disp.h"


#define CHECK_SCIP(x) { RETCODE _retcode_; \
                        if( (_retcode_ = (x)) != SCIP_OKAY ) \
                          SCIPerror(stderr, _retcode_, __FILE__, __LINE__); \
                      }

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

extern
RETCODE SCIPcreate(                     /**< creates and initializes SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   );

extern
RETCODE SCIPfree(                       /**< frees SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   );

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
RETCODE SCIPsolve(                      /**< solves problem */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
RETCODE SCIPfreeSolve(                  /**< frees all solution process data, only original problem is kept */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
RETCODE SCIPcreateVar(                  /**< create problem variable */
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
RETCODE SCIPcreateRow(                  /**< creates an LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row,                /**< pointer to row */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs                 /**< right hand side of row */
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
RETCODE SCIPchgLhs(                     /**< changes left hand side of LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             lhs                 /**< new left hand side */
   );

extern
RETCODE SCIPchgRhs(                     /**< changes right hand side of LP row */
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
RETCODE SCIPgetRowFeasibility(          /**< returns the feasibility of a row in the last LP solution */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            feasibility         /**< pointer to store the row's feasibility */
   );

extern
RETCODE SCIPprintRow(                   /**< output row to file stream */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

extern
RETCODE SCIPaddCut(                     /**< adds cut to separation storage */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             cut,                /**< separated cut */
   Real             score,              /**< separation score of cut (the larger, the better the cut) */
   Bool             pool                /**< should the cut be used in the global cut pool? Cut must be global valid! */
   );

extern
RETCODE SCIPcreateChild(                /**< creates a child node of the active node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           node                /**< pointer to node data structure */
   );

extern
RETCODE SCIPincludeReader(              /**< creates a reader and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
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
   DECL_CONSINIT((*consinit)),          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit)),          /**< deinitialise constraint handler */
   DECL_CONSFREE((*consfree)),          /**< free specific constraint data */
   DECL_CONSTRAN((*constran)),          /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA((*conssepa)),          /**< separate cutting planes */
   DECL_CONSENFO((*consenfo)),          /**< enforcing constraints */
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
RETCODE SCIPincludeNodesel(             /**< creates a node selector and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of node selector */
   const char*      desc,               /**< description of node selector */
   DECL_NODESELINIT((*nodeselinit)),    /**< initialise node selector */
   DECL_NODESELEXIT((*nodeselexit)),    /**< deinitialise node selector */
   DECL_NODESELSLCT((*nodeselslct)),    /**< node selection method */
   DECL_NODESELCOMP((*nodeselcomp)),    /**< node comparison method */
   NODESELDATA*     nodeseldata,        /**< node selector data */
   Bool             lowestboundfirst    /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
   );

extern
RETCODE SCIPincludeDisp(                /**< creates a display column and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of display column */
   const char*      desc,               /**< description of display column */
   const char*      header,             /**< head line of display column */
   DECL_DISPINIT((*dispinit)),          /**< initialise display column */
   DECL_DISPEXIT((*dispexit)),          /**< deinitialise display column */
   DECL_DISPOUTP((*dispoutp)),          /**< output method */
   DISPDATA*        dispdata,           /**< display column data */
   int              width,              /**< width of display column (no. of chars used) */
   int              priority,           /**< priority of display column */
   int              position,           /**< relative position of display column */
   Bool             stripline           /**< should the column be separated with a line from its right neighbour? */
   );

extern
RETCODE SCIPcreateCons(                 /**< creates a constraint of the given constraint handler */
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

extern
RETCODE SCIPaddCons(                    /**< adds constraint to the problem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   );

extern
RETCODE SCIPaddLocalCons(               /**< adds local constraint to the actual subproblem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   );

extern
RETCODE SCIPfindCons(                   /**< finds constraint of given name in the problem */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint to find */
   CONS**           cons                /**< pointer to store the constraint, returns NULL if not found */
   );

extern
RETCODE SCIPchgNodeBd(                  /**< changes bound of variable at the given node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

extern
RETCODE SCIPchgNodeLb(                  /**< changes lower bound of variable in the given node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

extern
RETCODE SCIPchgNodeUb(                  /**< changes upper bound of variable in the given node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

extern
RETCODE SCIPchgLocalLb(                 /**< changes lower bound of variable in the active node */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

extern
RETCODE SCIPchgLocalUb(                 /**< changes upper bound of variable in the active node */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

extern
RETCODE SCIPchgLb(                      /**< changes lower bound of variable in the problem */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

extern
RETCODE SCIPchgUb(                      /**< changes upper bound of variable in the problem */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   );

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
RETCODE SCIPgetBestLeaf(                /**< gets the best leaf from the node queue */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestleaf            /**< pointer to store best leaf */
   );

extern
RETCODE SCIPgetBestNode(                /**< gets the best node from the tree (child, sibling, or leaf) */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestnode            /**< pointer to store best leaf */
   );

extern
RETCODE SCIPgetNodenum(                 /**< gets number of processed nodes, including the active node */
   SCIP*            scip,               /**< SCIP data structure */
   int*             nodenum             /**< pointer to store the number of processed nodes */
   );

extern
RETCODE SCIPgetNNodesLeft(              /**< gets number of nodes left in the tree (children + siblings + leaves) */
   SCIP*            scip,               /**< SCIP data structure */
   int*             nnodes              /**< pointer to store the number of processed nodes */
   );

extern
RETCODE SCIPgetNLPIterations(           /**< gets number of simplex iterations used so far */
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
RETCODE SCIPgetActNCols(                /**< gets number of columns in actual LP */
   SCIP*            scip,               /**< SCIP data structure */
   int*             actncols            /**< pointer to store the number of columns */
   );

extern
RETCODE SCIPgetActNRows(                /**< gets number of rows in actual LP */
   SCIP*            scip,               /**< SCIP data structure */
   int*             actnrows            /**< pointer to store the number of columns */
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
RETCODE SCIPgetDualBound(               /**< gets actual dual bound */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            dualbound           /**< pointer to store the dual bound */
   );

extern
RETCODE SCIPgetPrimalBound(             /**< gets actual primal bound */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            primalbound         /**< pointer to store the primal bound */
   );

extern
RETCODE SCIPgetBestSol(                 /**< gets best feasible primal solution found so far */
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol                 /**< pointer to store the solution, returns NULL if no solution available */
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
Real SCIPinfinity(                      /**< returns value treated as infinity */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
Real SCIPepsilon(                       /**< returns value treated as zero */
   SCIP*            scip                /**< SCIP data structure */
   );

extern
Bool SCIPisEQ(                          /**< checks, if values are in range of epsZero */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisL(                           /**< checks, if val1 is (more than epsZero) lower than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisLE(                          /**< checks, if val1 is not (more than epsZero) greater than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisG(                           /**< checks, if val1 is (more than epsZero) greater than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPisGE(                          /**< checks, if val1 is not (more than epsZero) lower than val2 */
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
Bool SCIPisZero(                        /**< checks, if value is in range epsZero of 0.0 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisPos(                         /**< checks, if value is greater than epsZero */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPisNeg(                         /**< checks, if value is lower than -epsZero */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
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
 * debug methods
 */

#ifndef NDEBUG

extern
void SCIPdebugMemory(                   /**< prints output about used memory */
   SCIP*            scip                /**< SCIP data structure */
   );

#endif


#endif
