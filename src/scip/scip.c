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

/**@file   scip.c
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "set.h"
#include "mem.h"
#include "prob.h"
#include "tree.h"
#include "stat.h"
#include "scip.h"
#include "solve.h"
#include "price.h"


/** SCIP main data structure */
struct Scip
{
   STAGE            stage;              /**< SCIP operation stage */
   SET*             set;                /**< global SCIP settings */
   MEM*             mem;                /**< block memory buffers */
   PROB*            prob;               /**< original problem data */
   STAT*            stat;               /**< dynamic problem statistics */
   PROB*            transprob;          /**< transformed problem after presolve */
   TREE*            tree;               /**< branch and bound tree */
   LP*              lp;                 /**< LP data */
   PRICE*           price;              /**< storage for priced variables */
};



Real SCIPversion(                       /**< returns scip version number */
   void
   )
{
   return (Real)(SCIP_VERSION)/100.0;
}

void SCIPerror(                         /**< prints error message and aborts program execution */
   FILE*            errout,             /**< file stream to write error message */
   RETCODE          retcode,            /**< SCIP return code causing the error */
   const char*      filename,           /**< source code file name */
   int              line                /**< source line */
   )
{
   fprintf(errout, "[%s:%d] SCIP Error (%d): ", filename, line, retcode);
   SCIPretcodePrint(errout, retcode);
   fprintf(errout, "\n");
   abort();
}

RETCODE SCIPcreate(                     /**< creates and initializes SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   ALLOC_OKAY( allocMemory(*scip) );

   (*scip)->stage = SCIP_STAGE_INIT;

   CHECK_OKAY( SCIPsetCreate(&(*scip)->set, *scip) );
   CHECK_OKAY( SCIPmemCreate(&(*scip)->mem) );
   (*scip)->prob = NULL;
   (*scip)->stat = NULL;
   (*scip)->transprob = NULL;
   (*scip)->tree = NULL;
   (*scip)->lp = NULL;
   (*scip)->price = NULL;

   return SCIP_OKAY;
}

RETCODE SCIPfree(                       /**< frees SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(*scip != NULL);

   CHECK_OKAY( SCIPfreeProb(*scip) );
   assert((*scip)->stage == SCIP_STAGE_INIT);

   CHECK_OKAY( SCIPmemFree(&(*scip)->mem) );
   CHECK_OKAY( SCIPsetFree(&(*scip)->set) );

   freeMemory(*scip);

   return SCIP_OKAY;
}

RETCODE SCIPcreateProb(                 /**< creates empty problem and initializes all solving data structures */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< problem name */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      scip->stage = SCIP_STAGE_PROBLEM;
      
      CHECK_OKAY( SCIPsetInitCallbacks(scip->set) );
      CHECK_OKAY( SCIPprobCreate(&scip->prob, name) );
      CHECK_OKAY( SCIPstatCreate(&scip->stat) );

      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      errorMessage("Problem already exists");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

RETCODE SCIPfreeProb(                   /**< frees problem and solution process data */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   CHECK_OKAY( SCIPfreeSolve(scip) );

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);
      
      CHECK_OKAY( SCIPprobFree(&scip->prob, scip->mem->probmem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPstatFree(&scip->stat) );
      CHECK_OKAY( SCIPsetExitCallbacks(scip->set) );

      scip->stage = SCIP_STAGE_INIT;

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      errorMessage("Cannot free problem while solving is in progress");
      return SCIP_ERROR;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot free problem while solution is at hand");
      return SCIP_ERROR;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPsolve(                      /**< solves problem */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      errorMessage("Cannot solve problem before problem was created");
      return SCIP_NOPROBLEM;

   case SCIP_STAGE_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      /* switch stage to SOLVING */
      scip->stage = SCIP_STAGE_SOLVING;

      /* init solve data structures */
      CHECK_OKAY( SCIPlpCreate(&scip->lp, SCIPprobGetName(scip->prob)) );
      CHECK_OKAY( SCIPtreeCreate(&scip->tree, scip->mem->solvemem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPpriceCreate(&scip->price) );

      /* copy problem in solve memory */
      CHECK_OKAY( SCIPprobTransform(scip->prob, scip->mem->solvemem, scip->set, &scip->transprob) );

      /* presolve problem */
      /* ??? */
      todoMessage("problem presolving");

      /* fallthrough */

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      /* continue solution process */
      CHECK_OKAY( SCIPsolveCIP(scip->set, scip->mem->solvemem, scip->stat, scip->transprob, scip->tree, 
                     scip->lp, scip->price) );

      /* detect, whether problem is solved */
      /* ??? if( ... ) */
      scip->stage = SCIP_STAGE_SOLVED;

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      return SCIP_OKAY;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPfreeSolve(                  /**< frees all solution process data, only original problem is kept */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      CHECK_OKAY( SCIPprobFree(&scip->transprob, scip->mem->solvemem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPpriceFree(&scip->price) );
      CHECK_OKAY( SCIPtreeFree(&scip->tree) );
      CHECK_OKAY( SCIPlpFree(&scip->lp) );

      clearBlockMemoryNull(scip->mem->solvemem);

      scip->stage = SCIP_STAGE_PROBLEM;

      return SCIP_OKAY;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPcreateVar(                  /**< create problem variable */
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var,                /**< pointer to variable object */
   const char*      name,               /**< name of column */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype             /**< type of variable */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(var != NULL);
   assert(name != NULL);
   assert(lb <= ub);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      errorMessage("Cannot create variables before problem was created");
      return SCIP_NOPROBLEM;

   case SCIP_STAGE_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      CHECK_OKAY( SCIPvarCreate(var, scip->mem->probmem, scip->set, name, lb, ub, obj, vartype) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      CHECK_OKAY( SCIPvarCreateTransformed(var, scip->mem->solvemem, scip->set, name, lb, ub, obj, vartype) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot create additional columns after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPaddVar(                     /**< adds variable to the problem */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to add */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(var != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      errorMessage("Cannot add variables before problem was created");
      return SCIP_NOPROBLEM;

   case SCIP_STAGE_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      if( var->varstatus != SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("Cannot add transformed variables to original problem");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->prob, scip->set, var) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("Cannot add original variables to transformed problem");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->transprob, scip->set, var) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot add additional variables after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPcreateRow(                  /**< creates an LP row and captures it */
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row,                /**< pointer to row */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             rhs,                /**< right hand side of row */
   Real             epsilon,            /**< maximal normed violation of row */
   Bool             equality            /**< is row an equality? otherwise, it is a lower or equal inequality */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(row != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      errorMessage("Cannot create rows before problem was created");
      return SCIP_NOPROBLEM;

   case SCIP_STAGE_PROBLEM:
      errorMessage("Cannot create rows before problem solution was started");
      return SCIP_INVALIDCALL;

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      CHECK_OKAY( SCIProwCreate(row, scip->mem->solvemem, scip->set, scip->lp, scip->stat, name, 
                     len, col, val, rhs, epsilon, equality) );
      CHECK_OKAY( SCIPcaptureRow(scip, *row) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot create additional rows after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPcaptureRow(                 /**< increases usage counter of LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< row to capture */
   )
{
   assert(scip != NULL);

   SCIProwCapture(row);

   return SCIP_OKAY;
}

RETCODE SCIPreleaseRow(                 /**< decreases usage counter of LP row, and frees memory if necessary */
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row                 /**< pointer to LP row */
   )
{
   assert(scip != NULL);
   assert(row != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      errorMessage("Cannot release rows outside problem solving process");
      return SCIP_INVALIDCALL;

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      SCIProwRelease(row, scip->mem->solvemem, scip->set, scip->lp);

      return SCIP_OKAY;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPincludeConsHdlr(            /**< creates a constraint handler and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   DECL_CONSINIT((*consinit)),          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit)),          /**< deinitialise constraint handler */
   DECL_CONSFREE((*consfree)),          /**< frees specific constraint data */
   DECL_CONSTRAN((*constran)),          /**< transforms constraint data into data belonging to the transformed problem */
   DECL_CONSCHCK((*conschck)),          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop)),          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   CONSHDLR* conshdlr;

   assert(scip != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      CHECK_OKAY( SCIPconshdlrCreate(&conshdlr, name, desc, 
                     consinit, consexit, consfree, constran, conschck, consprop, conshdlrdata) );
      CHECK_OKAY( SCIPsetIncludeConsHdlr(scip->set, conshdlr) );

      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      CHECK_OKAY( SCIPconshdlrCreate(&conshdlr, name, desc, 
                     consinit, consexit, consfree, constran, conschck, consprop, conshdlrdata) );
      CHECK_OKAY( SCIPsetIncludeConsHdlr(scip->set, conshdlr) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      errorMessage("Cannot create additional constraint handlers while problem is solved");
      return SCIP_INVALIDCALL;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot create additional constraint handlers after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPfindConsHdlr(               /**< finds the constraint handler of the given name */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint handler */
   CONSHDLR**       conshdlr            /**< pointer for storing the constraint handler (returns NULL, if not found) */
   )
{
   assert(scip != NULL);
   assert(name != NULL);
   assert(conshdlr != NULL);

   CHECK_OKAY( SCIPsetFindConsHdlr(scip->set, name, conshdlr) );

   return SCIP_OKAY;
}

const char* SCIPgetConsHdlrName(        /**< gets name of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handlert */
   )
{
   return SCIPconshdlrGetName(conshdlr);
}

RETCODE SCIPincludeNodesel(             /**< creates a node selector and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of node selector */
   const char*      desc,               /**< description of node selector */
   DECL_NODESELINIT((*nodeselinit)),    /**< initialise node selector */
   DECL_NODESELEXIT((*nodeselexit)),    /**< deinitialise node selector */
   DECL_NODESELSLCT((*nodeselslct)),    /**< node selection method */
   DECL_NODESELCOMP((*nodeselcomp)),    /**< node comparison method */
   NODESELDATA*     nodeseldata         /**< node selector data */
   )
{
   NODESEL* nodesel;

   assert(scip != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      CHECK_OKAY( SCIPnodeselCreate(&nodesel, name, desc,
                     nodeselinit, nodeselexit, nodeselslct, nodeselcomp, nodeseldata) );
      CHECK_OKAY( SCIPsetIncludeNodesel(scip->set, nodesel) );

      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      CHECK_OKAY( SCIPnodeselCreate(&nodesel, name, desc,
                     nodeselinit, nodeselexit, nodeselslct, nodeselcomp, nodeseldata) );
      CHECK_OKAY( SCIPsetIncludeNodesel(scip->set, nodesel) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      errorMessage("Cannot create additional node selectors while problem is solved");
      return SCIP_INVALIDCALL;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot create additional node selectors after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

const char* SCIPgetNodeselName(         /**< gets name of node selector */
   NODESEL*         nodesel             /**< node selector */
   )
{
   return SCIPnodeselGetName(nodesel);
}

RETCODE SCIPcreateCons(                 /**< creates a constraint of the given constraint handler */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             model               /**< is constraint necessary for feasibility? */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(cons != NULL);
   assert(conshdlr != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->probmem, model, conshdlr, consdata) );

      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->probmem, model, conshdlr, consdata) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->solvemem, model, conshdlr, consdata) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot create additional constraints after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPaddCons(                    /**< adds global constraint to the problem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(cons != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      errorMessage("Cannot add constraints before problem was created");
      return SCIP_NOPROBLEM;

   case SCIP_STAGE_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);

      CHECK_OKAY( SCIPprobAddCons(scip->prob, scip->mem->probmem, cons) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      CHECK_OKAY( SCIPtreeAddGlobalCons(scip->tree, scip->mem->solvemem, cons) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot add additional constraints after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPaddLocalCons(               /**< adds local constraint to the actual subproblem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(cons != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      errorMessage("Cannot add local constraints before problem was created");
      return SCIP_NOPROBLEM;

   case SCIP_STAGE_PROBLEM:
      errorMessage("Cannot add local constraints before problem solution process was started");
      return SCIP_INVALIDCALL;

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      CHECK_OKAY( SCIPtreeAddLocalCons(scip->tree, scip->mem->solvemem, cons) );

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot add additional constraints after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPgetChildren(                /**< gets children of actual node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          children,           /**< pointer to store children array */
   int*             nchildren           /**< pointer to store number of children */
   )
{
   assert(scip != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      errorMessage("Cannot get children before problem was created");
      return SCIP_NOPROBLEM;

   case SCIP_STAGE_PROBLEM:
      errorMessage("Cannot get children before problem solution process was started");
      return SCIP_INVALIDCALL;

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      *children = SCIPtreeGetChildren(scip->tree);
      *nchildren = SCIPtreeGetNChildren(scip->tree);

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot get children after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPgetSiblings(                /**< gets siblings of actual node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          siblings,           /**< pointer to store siblings array */
   int*             nsiblings           /**< pointer to store number of siblings */
   )
{
   assert(scip != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      errorMessage("Cannot get siblings before problem was created");
      return SCIP_NOPROBLEM;

   case SCIP_STAGE_PROBLEM:
      errorMessage("Cannot get siblings before problem solution process was started");
      return SCIP_INVALIDCALL;

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      *siblings = SCIPtreeGetSiblings(scip->tree);
      *nsiblings = SCIPtreeGetNSiblings(scip->tree);

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot get siblings after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPgetBestLeaf(                /**< gets the best leaf from the node queue */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestleaf            /**< pointer to store best leaf */
   )
{
   assert(scip != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      errorMessage("Cannot get best leaf before problem was created");
      return SCIP_NOPROBLEM;

   case SCIP_STAGE_PROBLEM:
      errorMessage("Cannot get best leaf before problem solution process was started");
      return SCIP_INVALIDCALL;

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      *bestleaf = SCIPtreeGetBestLeaf(scip->tree);
      
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot get best leaf after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPgetBestNode(                /**< gets the best node from the tree (child, sibling, or leaf) */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestnode            /**< pointer to store best node */
   )
{
   assert(scip != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      errorMessage("Cannot get best node before problem was created");
      return SCIP_NOPROBLEM;

   case SCIP_STAGE_PROBLEM:
      errorMessage("Cannot get best node before problem solution process was started");
      return SCIP_INVALIDCALL;

   case SCIP_STAGE_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);

      *bestnode = SCIPtreeGetBestNode(scip->tree, scip->set);

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      errorMessage("Cannot get best node after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}

VERBLEVEL SCIPverbLevel(                /**< gets verbosity level for message output */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   return scip->set->verblevel;
}

RETCODE SCIPsetVerbLevel(               /**< sets verbosity level for message output */
   SCIP*            scip,               /**< SCIP data structure */
   VERBLEVEL        verblevel           /**< verbosity level for message output */
   )
{
   return SCIPsetSetVerbLevel(scip->set, verblevel);
}

STAGE SCIPstage(                        /**< returns current stage of SCIP */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->stage;
}

MEMHDR* SCIPmemhdr(                     /**< returns block memory to use at the current time */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      return scip->mem->probmem;

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->mem->solvemem;

   default:
      errorMessage("Unknown SCIP stage");
      return NULL;
   }
}

Real SCIPinfinity(                      /**< returns value treated as infinity */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->set->infinity;
}

Real SCIPepsilon(                       /**< returns value treated as zero */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->set->epsilon;
}

Bool SCIPisEQ(                          /**< checks, if values are in range of epsZero */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);

   return SCIPsetIsEQ(scip->set, val1, val2);
}

Bool SCIPisL(                           /**< checks, if val1 is (more than epsZero) lower than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);

   return SCIPsetIsL(scip->set, val1, val2);
}

Bool SCIPisLE(                          /**< checks, if val1 is not (more than epsZero) greater than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);

   return SCIPsetIsLE(scip->set, val1, val2);
}

Bool SCIPisG(                           /**< checks, if val1 is (more than epsZero) greater than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);

   return SCIPsetIsG(scip->set, val1, val2);
}

Bool SCIPisGE(                          /**< checks, if val1 is not (more than epsZero) lower than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);

   return SCIPsetIsGE(scip->set, val1, val2);
}

Bool SCIPisInfinity(                    /**< checks, if value is infinite */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against infinity */
   )
{
   assert(scip != NULL);

   return SCIPsetIsInfinity(scip->set, val);
}

Bool SCIPisZero(                        /**< checks, if value is in range epsZero of 0.0 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);

   return SCIPsetIsZero(scip->set, val);
}

Bool SCIPisPos(                         /**< checks, if value is greater than epsZero */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);

   return SCIPsetIsPos(scip->set, val);
}

Bool SCIPisNeg(                         /**< checks, if value is lower than -epsZero */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);

   return SCIPsetIsNeg(scip->set, val);
}



/*
 * debug methods
 */

#ifndef NDEBUG

void SCIPdebugMemory(                   /**< prints output about used memory */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   memoryDiagnostic();

   printf("\nProblem Block Memory (%p):\n", scip->mem->probmem);
   blockMemoryDiagnostic(scip->mem->probmem);

   printf("\nSolution Block Memory (%p):\n", scip->mem->solvemem);
   blockMemoryDiagnostic(scip->mem->solvemem);

#if 0
   printf("\nTemporary Block Memory (%p):\n", scip->mem->tempmem);
   blockMemoryDiagnostic(scip->mem->tempmem);
#endif
}

#endif
