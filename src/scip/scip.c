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
#include "lp.h"
#include "stat.h"
#include "scip.h"


/** SCIP operation status */
enum Status
{
   SCIP_STATUS_INIT    = 0,             /**< SCIP datastructures are initialized, no problem exists */
   SCIP_STATUS_PROBLEM = 1,             /**< the problem is being created and modified */
   SCIP_STATUS_SOLVING = 2,             /**< the problem is being solved */
   SCIP_STATUS_SOLVED  = 3              /**< the problem is solved */
};
typedef enum Status STATUS;

/** SCIP main data structure */
struct Scip
{
   STATUS           status;             /**< SCIP operation status */
   SET*             set;                /**< global SCIP settings */
   MEM*             mem;                /**< block memory buffers */
   PROB*            prob;               /**< original problem data */
   STAT*            stat;               /**< dynamic problem statistics */
   PROB*            transprob;          /**< transformed problem after presolve */
   TREE*            tree;               /**< branch and bound tree */
   LP*              lp;                 /**< LP data */
};



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

   (*scip)->status = SCIP_STATUS_INIT;
   CHECK_OKAY( SCIPsetCreate(&(*scip)->set) );
   CHECK_OKAY( SCIPmemCreate(&(*scip)->mem) );
   (*scip)->prob = NULL;
   (*scip)->stat = NULL;
   (*scip)->transprob = NULL;
   (*scip)->tree = NULL;
   (*scip)->lp = NULL;

   return SCIP_OKAY;
}

RETCODE SCIPfree(                       /**< frees SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(*scip != NULL);

   CHECK_OKAY( SCIPfreeProb(*scip) );
   assert((*scip)->status == SCIP_STATUS_INIT);

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

   switch( scip->status )
   {
   case SCIP_STATUS_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      CHECK_OKAY( SCIPprobCreate(&scip->prob, name) );
      CHECK_OKAY( SCIPstatCreate(&scip->stat) );
      scip->status = SCIP_STATUS_PROBLEM;
      
      return SCIP_OKAY;

   case SCIP_STATUS_PROBLEM:
   case SCIP_STATUS_SOLVING:
   case SCIP_STATUS_SOLVED:
      errorMessage("Problem already exists");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP status");
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

   switch( scip->status )
   {
   case SCIP_STATUS_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      return SCIP_OKAY;

   case SCIP_STATUS_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      CHECK_OKAY( SCIPprobFree(&scip->prob) );
      CHECK_OKAY( SCIPstatFree(&scip->stat) );
      scip->status = SCIP_STATUS_INIT;

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVING:
      errorMessage("Cannot free problem while solving is in progress");
      return SCIP_ERROR;

   case SCIP_STATUS_SOLVED:
      errorMessage("Cannot free problem while solution is at hand");
      return SCIP_ERROR;

   default:
      errorMessage("Unknown SCIP status");
      return SCIP_ERROR;
   }
}

RETCODE SCIPsolve(                      /**< solves problem */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   switch( scip->status )
   {
   case SCIP_STATUS_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      return SCIP_NOPROBLEM;

   case SCIP_STATUS_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      /* init solve data structures */
      CHECK_OKAY( SCIPtreeCreate(&scip->tree, scip->set) );
      CHECK_OKAY( SCIPlpCreate(&scip->lp, SCIPprobGetName(scip->prob)) );

      /* copy problem in solve memory */
      CHECK_OKAY( SCIPprobDuplicate(&scip->transprob, scip->mem->solvemem, scip->prob) );

      /* preprocess problem */
      /* ??? */

      /* init tree with root relaxation */
      /* ??? */
      
      scip->status = SCIP_STATUS_SOLVING;
      /* fallthrough */

   case SCIP_STATUS_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);

      /* continue solution process */
      /* ??? */

      /* detect, whether problem is solved */
      /* ??? if( ... ) */
      scip->status = SCIP_STATUS_SOLVED;

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVED:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);

      return SCIP_OKAY;

   default:
      errorMessage("Unknown SCIP status");
      return SCIP_ERROR;
   }
}

RETCODE SCIPfreeSolve(                  /**< frees all solution process data, only original problem is kept */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   switch( scip->status )
   {
   case SCIP_STATUS_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      return SCIP_OKAY;

   case SCIP_STATUS_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVING:
   case SCIP_STATUS_SOLVED:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);

      CHECK_OKAY( SCIPprobFree(&scip->transprob) );
      CHECK_OKAY( SCIPtreeFree(&scip->tree) );
      CHECK_OKAY( SCIPlpFree(&scip->lp) );

      clearBlockMemoryNull(scip->mem->solvemem);

      scip->status = SCIP_STATUS_PROBLEM;

      return SCIP_OKAY;

   default:
      errorMessage("Unknown SCIP status");
      return SCIP_ERROR;
   }
}

RETCODE SCIPcreateCol(                  /**< create variable with empty column */
   SCIP*            scip,               /**< SCIP data structure */
   COL**            col,                /**< pointer to column object */
   const char*      name,               /**< name of column */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   COLTYPE          coltype             /**< type of variable */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(col != NULL);
   assert(name != NULL);
   assert(lb <= ub);

   switch( scip->status )
   {
   case SCIP_STATUS_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      return SCIP_NOPROBLEM;

   case SCIP_STATUS_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      CHECK_OKAY( SCIPcolCreate(col, scip->mem->probmem, scip->set, scip->lp, scip->stat, name, 
                     0, NULL, NULL, lb, ub, obj, coltype) );

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);

      CHECK_OKAY( SCIPcolCreate(col, scip->mem->solvemem, scip->set, scip->lp, scip->stat, name, 
                     0, NULL, NULL, lb, ub, obj, coltype) );

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVED:
      errorMessage("Cannot create additional columns after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP status");
      return SCIP_ERROR;
   }
}

RETCODE SCIPcreateRow(                  /**< creates an LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row,                /**< pointer to row */
   const char*      name,               /**< name of row */
   int              len,                /**< number of nonzeros in the row */
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             rhs,                /**< right hand side of row */
   Real             lhs,                /**< left hand side of row (for ranged rows) */
   Real             epsilon,            /**< maximal normed violation of row */
   ROWTYPE          rowtype             /**< type of row */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(row != NULL);

   switch( scip->status )
   {
   case SCIP_STATUS_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      return SCIP_NOPROBLEM;

   case SCIP_STATUS_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      CHECK_OKAY( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->lp, scip->stat, name, 
                     len, col, val, rhs, lhs, epsilon, rowtype) );

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);

      CHECK_OKAY( SCIProwCreate(row, scip->mem->solvemem, scip->set, scip->lp, scip->stat, name, 
                     len, col, val, rhs, lhs, epsilon, rowtype) );

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVED:
      errorMessage("Cannot create additional rows after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP status");
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

RETCODE SCIPcreateConstraint(           /**< creates a constraint of the given constraint handler */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             original,           /**< belongs constraint to the original problem formulation? */
   Bool             model               /**< is constraint necessary for feasibility? */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(cons != NULL);
   assert(conshdlr != NULL);

   switch( scip->status )
   {
   case SCIP_STATUS_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      return SCIP_NOPROBLEM;

   case SCIP_STATUS_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->probmem, original, model, conshdlr, consdata) );

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);

      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->solvemem, original, model, conshdlr, consdata) );

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVED:
      errorMessage("Cannot create additional constraints after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP status");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

RETCODE SCIPaddConstraint(              /**< adds global constraint to the problem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(cons != NULL);

   switch( scip->status )
   {
   case SCIP_STATUS_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      return SCIP_NOPROBLEM;

   case SCIP_STATUS_PROBLEM:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      CHECK_OKAY( SCIPprobAddConstraint(scip->prob, scip->mem->probmem, cons) );

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);

      CHECK_OKAY( SCIPtreeAddGlobalConstraint(scip->tree, scip->mem->solvemem, cons) );

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVED:
      errorMessage("Cannot add additional constraints after problem solution process was started");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP status");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

RETCODE SCIPaddLocalConstraint(         /**< adds local constraint to the actual subproblem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(cons != NULL);

   switch( scip->status )
   {
   case SCIP_STATUS_INIT:
      assert(scip->prob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);

      return SCIP_NOPROBLEM;

   case SCIP_STATUS_PROBLEM:
      errorMessage("Cannot add local constraints before problem solution process was started");
      return SCIP_INVALIDCALL;

   case SCIP_STATUS_SOLVING:
      assert(scip->prob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);

      CHECK_OKAY( SCIPtreeAddLocalConstraint(scip->tree, scip->mem->solvemem, cons) );

      return SCIP_OKAY;

   case SCIP_STATUS_SOLVED:
      errorMessage("Cannot add additional constraints after problem was solved");
      return SCIP_INVALIDCALL;

   default:
      errorMessage("Unknown SCIP status");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
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

   printf("\nTemporary Block Memory (%p):\n", scip->mem->tempmem);
   blockMemoryDiagnostic(scip->mem->tempmem);
}

#endif
