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
#include "sepa.h"
#include "primal.h"
#include "lpi.h"


/** SCIP main data structure */
struct Scip
{
   STAGE            stage;              /**< SCIP operation stage */
   SET*             set;                /**< global SCIP settings */
   MEM*             mem;                /**< block memory buffers */
   PROB*            origprob;           /**< original problem data */
   PROB*            transprob;          /**< transformed problem after presolve */
   STAT*            stat;               /**< dynamic problem statistics */
   TREE*            tree;               /**< branch and bound tree */
   LP*              lp;                 /**< LP data */
   PRICE*           price;              /**< storage for priced variables */
   SEPA*            sepa;               /**< storage for separated cuts */
   PRIMAL*          primal;             /**< primal data and solution storage */
};



static
RETCODE checkStage(                     /**< checks, if SCIP is in one of the feasible stages */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      method,             /**< method that was called */
   Bool             init,               /**< may method be called in the INIT stage? */
   Bool             problem,            /**< may method be called in the PROBLEM stage? */
   Bool             initsolve,          /**< may method be called in the INITSOLVE stage? */
   Bool             presolving,         /**< may method be called in the PRESOLVING stage? */
   Bool             solving,            /**< may method be called in the SOLVING stage? */
   Bool             solved,             /**< may method be called in the SOLVED stage? */
   Bool             freesolve           /**< may method be called in the FREESOLVE stage? */
   )
{
   char s[255];

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);
   assert(method != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      assert(scip->origprob == NULL);
      assert(scip->stat == NULL);
      assert(scip->transprob == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);
      assert(scip->sepa == NULL);
      assert(scip->primal == NULL);
      assert(scip->tree == NULL);

      if( !init )
      {
         sprintf(s, "Cannot call method <%s> in initialization stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob == NULL);
      assert(scip->lp == NULL);
      assert(scip->price == NULL);
      assert(scip->sepa == NULL);
      assert(scip->primal == NULL);
      assert(scip->tree == NULL);

      if( !problem )
      {
         sprintf(s, "Cannot call method <%s> in problem creation stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);
      assert(scip->tree == NULL);

      if( !initsolve )
      {
         sprintf(s, "Cannot call method <%s> in solve initialization stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);
      assert(scip->sepa != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree == NULL);

      if( !presolving )
      {
         sprintf(s, "Cannot call method <%s> in presolving stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);
      assert(scip->sepa != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree != NULL);

      if( !solving )
      {
         sprintf(s, "Cannot call method <%s> in solving stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);
      assert(scip->transprob != NULL);
      assert(scip->lp != NULL);
      assert(scip->price != NULL);
      assert(scip->sepa != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree != NULL);

      if( !solved )
      {
         sprintf(s, "Cannot call method <%s> in problem solved stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_FREESOLVE:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);

      if( !freesolve )
      {
         sprintf(s, "Cannot call method <%s> in solve deinitialization stage", method);
         errorMessage(s);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   default:
      errorMessage("Unknown SCIP stage");
      return SCIP_ERROR;
   }
}




Real SCIPversion(                       /**< returns scip version number */
   void
   )
{
   return (Real)(SCIP_VERSION)/100.0;
}

void SCIPprintVersion(                  /**< prints a version information line to a file stream */
   FILE*            file                /**< file stream to write version information, or NULL for stdout */
   )
{
   if( file == NULL )
      file = stdout;

   fprintf(file, "SCIP version %g [precision: %d byte]", SCIPversion(), (int)sizeof(Real));
#ifdef SCIP_BLOCKMEMORY
   fprintf(file, " [memory: block]");
#else
   fprintf(file, " [memory: standard]");
#endif
#ifdef NDEBUG
   fprintf(file, " [mode: optimized]");
#else
   fprintf(file, " [mode: debugging]");
#endif
   fprintf(file, " [LP solver: %s]\n", SCIPlpiGetSolverName());
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
   (*scip)->origprob = NULL;
   (*scip)->stat = NULL;
   (*scip)->transprob = NULL;
   (*scip)->tree = NULL;
   (*scip)->lp = NULL;
   (*scip)->price = NULL;
   (*scip)->sepa = NULL;
   (*scip)->primal = NULL;

   return SCIP_OKAY;
}

RETCODE SCIPfree(                       /**< frees SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   CHECK_OKAY( checkStage(*scip, "SCIPfree", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPcreateProb", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   scip->stage = SCIP_STAGE_PROBLEM;
   
   CHECK_OKAY( SCIPsetInitCallbacks(scip->set) );
   CHECK_OKAY( SCIPprobCreate(&scip->origprob, name) );
   CHECK_OKAY( SCIPstatCreate(&scip->stat) );
   
   return SCIP_OKAY;
}

RETCODE SCIPreadProb(                   /**< reads problem from file and initializes all solving data structures */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< problem file name */
   )
{
   RETCODE retcode;
   Bool read;
   int i;

   assert(filename != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreadProb", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* try all readers until one could read the file */
   read = FALSE;
   for( i = 0; i < scip->set->nreaders && !read; ++i )
   {
      CHECK_OKAY( retcode = SCIPreaderRead(scip->set->readers[i], scip, filename) );
      read = (retcode != SCIP_DIDNOTRUN);
   }

   if( !read )
   {
      printf("No reader for input file <%s> available\n", filename);
      return SCIP_READERR;
   }

   return SCIP_OKAY;
}

RETCODE SCIPfreeProb(                   /**< frees problem and solution process data */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeProb", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   CHECK_OKAY( SCIPfreeSolve(scip) );
   assert(scip->stage == SCIP_STAGE_INIT || scip->stage == SCIP_STAGE_PROBLEM);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPprobFree(&scip->origprob, scip->mem->probmem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPstatFree(&scip->stat) );
      CHECK_OKAY( SCIPsetExitCallbacks(scip->set) );

      scip->stage = SCIP_STAGE_INIT;

      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPsetObjsense(                /**< sets objective sense of problem */
   SCIP*            scip,               /**< SCIP data structure */
   OBJSENSE         objsense            /**< new objective sense */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetObjsense", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   if( objsense != SCIP_OBJSENSE_MAXIMIZE && objsense != SCIP_OBJSENSE_MINIMIZE )
   {
      errorMessage("Invalid objective sense");
      return SCIP_INVALIDDATA;
   }

   SCIPprobSetObjsense(scip->origprob, objsense);
   
   return SCIP_OKAY;
}
   
RETCODE SCIPsolve(                      /**< solves problem */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsolve", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      /* switch stage to INITSOLVE */
      scip->stage = SCIP_STAGE_INITSOLVE;

      /* mark statistics before solving */
      SCIPstatMark(scip->stat);

      /* init solve data structures */
      CHECK_OKAY( SCIPlpCreate(&scip->lp, scip->mem->solvemem, scip->set, SCIPprobGetName(scip->origprob)) );
      CHECK_OKAY( SCIPpriceCreate(&scip->price) );
      CHECK_OKAY( SCIPsepaCreate(&scip->sepa) );

      /* copy problem in solve memory */
      CHECK_OKAY( SCIPprobTransform(scip->origprob, scip->mem->solvemem, scip->set, scip->stat, &scip->transprob) );

      /* activate constraints in the problem */
      CHECK_OKAY( SCIPprobActivate(scip->transprob, scip->set) );

      /* create primal solution storage */
      CHECK_OKAY( SCIPprimalCreate(&scip->primal, scip->mem->solvemem, scip->set, scip->transprob, scip->lp) );

      /* switch stage to PRESOLVING */
      scip->stage = SCIP_STAGE_PRESOLVING;

      /* presolve problem */
      /* ??? */
      todoMessage("problem presolving");

      /* create branch-and-bound tree */
      CHECK_OKAY( SCIPtreeCreate(&scip->tree, scip->mem->solvemem, scip->set, scip->stat, scip->lp, scip->transprob) );
      
      /* switch stage to SOLVING */
      scip->stage = SCIP_STAGE_SOLVING;

      /* fallthrough */

   case SCIP_STAGE_SOLVING:
      /* continue solution process */
      CHECK_OKAY( SCIPsolveCIP(scip->set, scip->mem->solvemem, scip->stat, scip->transprob, scip->tree, 
                     scip->lp, scip->price, scip->sepa, scip->primal) );

      /* detect, whether problem is solved */
      todoMessage("detect, whether problem is solved");
      scip->stage = SCIP_STAGE_SOLVED;

      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPfreeSolve(                  /**< frees all solution process data, only original problem is kept */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeSolve", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      /* switch stage to FREESOLVE */
      scip->stage = SCIP_STAGE_FREESOLVE;

      /* free primal solution storage */
      CHECK_OKAY( SCIPprimalFree(&scip->primal, scip->mem->solvemem, scip->set, scip->lp) );

      /* deactivate constraints in the problem */
      CHECK_OKAY( SCIPprobDeactivate(scip->transprob) );

      /* deactivate the active node */
      CHECK_OKAY( SCIPnodeActivate(NULL, scip->mem->solvemem, scip->set, scip->lp, scip->tree) );

      CHECK_OKAY( SCIPlpClear(scip->lp, scip->mem->solvemem, scip->set) );
      CHECK_OKAY( SCIPtreeFree(&scip->tree, scip->mem->solvemem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPprobFree(&scip->transprob, scip->mem->solvemem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPsepaFree(&scip->sepa) );
      CHECK_OKAY( SCIPpriceFree(&scip->price) );
      CHECK_OKAY( SCIPlpFree(&scip->lp, scip->mem->solvemem, scip->set) );

      clearBlockMemoryNull(scip->mem->solvemem);

      /* reset statistics to the point before solving started */
      SCIPstatReset(scip->stat);

      /* switch stage to PROBLEM */
      scip->stage = SCIP_STAGE_PROBLEM;

      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
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
   assert(var != NULL);
   assert(name != NULL);
   assert(lb <= ub);

   CHECK_OKAY( checkStage(scip, "SCIPcreateVar", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarCreate(var, scip->mem->probmem, scip->set, scip->stat, name, lb, ub, obj, vartype) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPvarCreateTransformed(var, scip->mem->solvemem, scip->set, scip->stat, name, lb, ub, obj, vartype) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPcaptureVar(                 /**< increases usage counter of variable */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to capture */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIPvarCapture(var);

   return SCIP_OKAY;
}

RETCODE SCIPreleaseVar(                 /**< decreases usage counter of variable, and frees memory if necessary */
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var                 /**< pointer to variable */
   )
{
   assert(var != NULL);
   assert(*var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreleaseVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPvarRelease(var, scip->mem->probmem, scip->set, scip->lp);
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      if( (*var)->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("cannot release original variables while solving the problem");
         return SCIP_INVALIDCALL;
      }
      SCIPvarRelease(var, scip->mem->solvemem, scip->set, scip->lp);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPaddVar(                     /**< adds variable to the problem */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to add */
   )
{
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddVar", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( var->varstatus != SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("Cannot add transformed variables to original problem");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->origprob, scip->mem->probmem, scip->set, var) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
      if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("Cannot add original variables to transformed problem");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->transprob, scip->mem->solvemem, scip->set, var) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPgetVars(                    /**< gets variables of the problem */
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*             nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*             nbin,               /**< pointer to store number of binary variables or NULL if not needed */
   int*             nint,               /**< pointer to store number of integer variables or NULL if not needed */
   int*             nimpl,              /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*             ncont               /**< pointer to store number of continous variables or NULL if not needed */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( vars != NULL )
         *vars = scip->origprob->vars;
      if( nvars != NULL )
         *nvars = scip->origprob->nvars;
      if( nbin != NULL )
         *nbin = scip->origprob->nbin;
      if( nint != NULL )
         *nint = scip->origprob->nint;
      if( nimpl != NULL )
         *nimpl = scip->origprob->nimpl;
      if( ncont != NULL )
         *ncont = scip->origprob->ncont;
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      if( vars != NULL )
         *vars = scip->transprob->vars;
      if( nvars != NULL )
         *nvars = scip->transprob->nvars;
      if( nbin != NULL )
         *nbin = scip->transprob->nbin;
      if( nint != NULL )
         *nint = scip->transprob->nint;
      if( nimpl != NULL )
         *nimpl = scip->transprob->nimpl;
      if( ncont != NULL )
         *ncont = scip->transprob->ncont;
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPfindVar(                    /**< finds variable of given name in the problem */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of variable to find */
   VAR**            var                 /**< pointer to store the variable, returns NULL if not found */
   )
{
   assert(name != NULL);
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfindVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      *var = SCIPprobFindVar(scip->origprob, name);
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      *var = SCIPprobFindVar(scip->transprob, name);
      if( *var == NULL )
         *var = SCIPprobFindVar(scip->origprob, name);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
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
   Real             lhs,                /**< left hand side of row */
   Real             rhs                 /**< right hand side of row */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcreateRow", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwCreate(row, scip->mem->solvemem, scip->set, scip->lp, scip->stat, name, 
                  len, col, val, lhs, rhs) );

   return SCIP_OKAY;
}

RETCODE SCIPcaptureRow(                 /**< increases usage counter of LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< row to capture */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   SCIProwCapture(row);

   return SCIP_OKAY;
}

RETCODE SCIPreleaseRow(                 /**< decreases usage counter of LP row, and frees memory if necessary */
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row                 /**< pointer to LP row */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreleaseRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE) );

   SCIProwRelease(row, scip->mem->solvemem, scip->set, scip->lp);
   
   return SCIP_OKAY;
}

RETCODE SCIPchgLhs(                     /**< changes left hand side of LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             lhs                 /**< new left hand side */
   )
{
   assert(row != NULL);

   errorMessage("sides of row must not be changed");
   abort();
#if 0
   CHECK_OKAY( checkStage(scip, "SCIPchgLhs", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwChgLhs(row, scip->set, scip->lp, lhs) );
   return SCIP_OKAY;
#endif
}

RETCODE SCIPchgRhs(                     /**< changes right hand side of LP row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             rhs                 /**< new right hand side */
   )
{
   assert(row != NULL);

   errorMessage("sides of row must not be changed");
   abort();
#if 0
   CHECK_OKAY( checkStage(scip, "SCIPchgRhs", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwChgRhs(row, scip->set, scip->lp, rhs) );
   return SCIP_OKAY;
#endif
}

RETCODE SCIPaddVarToRow(                /**< resolves variable to columns and adds them with the coefficient to the row */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   VAR*             var,                /**< problem variable */
   Real             val                 /**< value of coefficient */
   )
{
   assert(row != NULL);
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddVarToRow", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarAddToRow(var, scip->mem->solvemem, scip->set, scip->lp, scip->stat, row, val) );
   
   return SCIP_OKAY;
}

RETCODE SCIPgetRowFeasibility(          /**< returns the feasibility of a row in the last LP solution */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real*            feasibility         /**< pointer to store the row's feasibility */
   )
{
   assert(row != NULL);
   assert(feasibility != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetRowFeasibility", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *feasibility = SCIProwGetFeasibility(row, scip->stat);
   
   return SCIP_OKAY;
}

RETCODE SCIPprintRow(                   /**< output row to file stream */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPprintRow", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIProwPrint(row, file);

   return SCIP_OKAY;
}

RETCODE SCIPaddCut(                     /**< adds cut to separation storage */
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             cut,                /**< separated cut */
   Real             score,              /**< separation score of cut (the larger, the better the cut) */
   Bool             pool                /**< should the cut be used in the global cut pool? Cut must be global valid! */
   )
{
   assert(cut != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddCut", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsepaAddCut(scip->sepa, scip->set, cut, score, pool) );
   
   return SCIP_OKAY;
}

RETCODE SCIPcreateChild(                /**< creates a child node of the active node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           node                /**< pointer to node data structure */
   )
{
   assert(node != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcreateChild", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPnodeCreate(node, scip->mem->solvemem, scip->set, scip->tree) );
   
   return SCIP_OKAY;
}

RETCODE SCIPincludeReader(              /**< creates a reader and includes it in SCIP */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
   DECL_READERINIT((*readerinit)),      /**< initialise reader */
   DECL_READEREXIT((*readerexit)),      /**< deinitialise reader */
   DECL_READERREAD((*readerread)),      /**< read method */
   READERDATA*      readerdata          /**< reader data */
   )
{
   READER* reader;

   CHECK_OKAY( checkStage(scip, "SCIPincludeReader", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPreaderCreate(&reader, name, desc, extension, readerinit, readerexit, readerread, readerdata) );
   CHECK_OKAY( SCIPsetIncludeReader(scip->set, reader) );
   
   return SCIP_OKAY;
}

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
   )
{
   CONSHDLR* conshdlr;

   CHECK_OKAY( checkStage(scip, "SCIPincludeConsHdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconshdlrCreate(&conshdlr, name, desc, sepapriority, enfopriority, chckpriority,
                  consinit, consexit, consfree, constran, conssepa, consenfo, conschck, consprop, conshdlrdata) );
   CHECK_OKAY( SCIPsetIncludeConsHdlr(scip->set, conshdlr) );
   
   return SCIP_OKAY;
}

RETCODE SCIPfindConsHdlr(               /**< finds the constraint handler of the given name */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint handler */
   CONSHDLR**       conshdlr            /**< pointer for storing the constraint handler (returns NULL, if not found) */
   )
{
   assert(name != NULL);
   assert(conshdlr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfindConsHdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetFindConsHdlr(scip->set, name, conshdlr) );

   return SCIP_OKAY;
}

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
   )
{
   NODESEL* nodesel;

   CHECK_OKAY( checkStage(scip, "SCIPincludeNodesel", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPnodeselCreate(&nodesel, name, desc,
                  nodeselinit, nodeselexit, nodeselslct, nodeselcomp, nodeseldata, lowestboundfirst) );
   CHECK_OKAY( SCIPsetIncludeNodesel(scip->set, nodesel) );
   
   return SCIP_OKAY;
}

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
   )
{
   DISP* disp;

   CHECK_OKAY( checkStage(scip, "SCIPincludeDisp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPdispCreate(&disp, name, desc, header, dispinit, dispexit, dispoutp, dispdata,
                  width, priority, position, stripline) );
   CHECK_OKAY( SCIPsetIncludeDisp(scip->set, disp) );
   
   return SCIP_OKAY;
}

RETCODE SCIPcreateCons(                 /**< creates a constraint of the given constraint handler */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to constraint */
   const char*      name,               /**< name of constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             model               /**< is constraint necessary for feasibility? */
   )
{
   assert(cons != NULL);
   assert(name != NULL);
   assert(conshdlr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcreateCons", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->probmem, name, conshdlr, consdata, model, TRUE) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->solvemem, name, conshdlr, consdata, model, FALSE) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPcaptureCons(                /**< increases usage counter of constraint */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to capture */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

RETCODE SCIPreleaseCons(                /**< decreases usage counter of constraint, and frees memory if necessary */
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons                /**< pointer to constraint */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreleaseCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPconsRelease(cons, scip->mem->probmem, scip->set);
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      if( SCIPconsIsOriginal(*cons) )
      {
         errorMessage("cannot release original constraint while solving the problem");
         return SCIP_INVALIDCALL;
      }
      SCIPconsRelease(cons, scip->mem->solvemem, scip->set);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPaddCons(                    /**< adds global constraint to the problem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddCons", FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPprobAddCons(scip->origprob, scip->mem->probmem, cons) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPprobAddCons(scip->transprob, scip->mem->solvemem, cons) );
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPtreeAddGlobalCons(scip->tree, scip->mem->solvemem, scip->set, cons) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPaddLocalCons(               /**< adds local constraint to the actual subproblem */
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddLocalCons", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPtreeAddLocalCons(scip->tree, scip->mem->solvemem, scip->set, cons) );
   
   return SCIP_OKAY;
}

RETCODE SCIPfindCons(                   /**< finds constraint of given name in the problem */
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of constraint to find */
   CONS**           cons                /**< pointer to store the constraint, returns NULL if not found */
   )
{
   assert(name != NULL);
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfindCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      *cons = SCIPprobFindCons(scip->origprob, name);
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      *cons = SCIPprobFindCons(scip->transprob, name);
      if( *cons == NULL )
         *cons = SCIPprobFindCons(scip->origprob, name);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPchgNodeBd(                  /**< changes bound of variable at the given node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPchgNodeBd", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( node == NULL )
      node = scip->tree->actnode;
   
   CHECK_OKAY( SCIPnodeAddBoundchg(node, scip->mem->solvemem, scip->set, scip->lp, scip->tree,
                  var, newbound, boundtype) );
   
   return SCIP_OKAY;
}

RETCODE SCIPchgNodeLb(                  /**< changes lower bound of variable in the given node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgNodeLb", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPchgNodeBd(scip, node, var, newbound, SCIP_BOUNDTYPE_LOWER);
}

RETCODE SCIPchgNodeUb(                  /**< changes upper bound of variable in the given node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgNodeUb", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPchgNodeBd(scip, node, var, newbound, SCIP_BOUNDTYPE_UPPER);
}

RETCODE SCIPchgLocalLb(                 /**< changes lower bound of variable in the active node */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgLocalLb", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPchgNodeBd(scip, NULL, var, newbound, SCIP_BOUNDTYPE_LOWER);
}

RETCODE SCIPchgLocalUb(                 /**< changes upper bound of variable in the active node */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgLocalUb", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPchgNodeBd(scip, NULL, var, newbound, SCIP_BOUNDTYPE_UPPER);
}

RETCODE SCIPchgLb(                      /**< changes lower bound of variable in the problem */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgLb", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarChgLb(var, scip->mem->probmem, scip->set, scip->lp, scip->tree, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
      if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("cannot change bounds of original variables while solving the problem");
         return SCIP_INVALIDCALL;
      }
      CHECK_OKAY( SCIPvarChgLb(var, scip->mem->solvemem, scip->set, scip->lp, scip->tree, newbound) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPchgUb(                      /**< changes upper bound of variable in the problem */
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPchgUb", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarChgUb(var, scip->mem->probmem, scip->set, scip->lp, scip->tree, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
      if( var->varstatus == SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("cannot change bounds of original variables while solving the problem");
         return SCIP_INVALIDCALL;
      }
      CHECK_OKAY( SCIPvarChgUb(var, scip->mem->solvemem, scip->set, scip->lp, scip->tree, newbound) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage");
      return SCIP_ERROR;
   }
}

RETCODE SCIPgetChildren(                /**< gets children of active node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          children,           /**< pointer to store children array */
   int*             nchildren           /**< pointer to store number of children */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetChildren", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *children = scip->tree->children;
   *nchildren = scip->tree->nchildren;
   
   return SCIP_OKAY;
}

RETCODE SCIPgetSiblings(                /**< gets siblings of active node */
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          siblings,           /**< pointer to store siblings array */
   int*             nsiblings           /**< pointer to store number of siblings */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetSiblings", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *siblings = scip->tree->siblings;
   *nsiblings = scip->tree->nsiblings;
   
   return SCIP_OKAY;
}

RETCODE SCIPgetBestLeaf(                /**< gets the best leaf from the node queue */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestleaf            /**< pointer to store best leaf */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetBestLeaf", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *bestleaf = SCIPtreeGetBestLeaf(scip->tree);
   
   return SCIP_OKAY;
}

RETCODE SCIPgetBestNode(                /**< gets the best node from the tree (child, sibling, or leaf) */
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           bestnode            /**< pointer to store best node */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetBestNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   *bestnode = SCIPtreeGetBestNode(scip->tree, scip->set);
   
   return SCIP_OKAY;
}

RETCODE SCIPgetNodenum(                 /**< gets number of processed nodes, including the active node */
   SCIP*            scip,               /**< SCIP data structure */
   int*             nodenum             /**< pointer to store the number of processed nodes */
   )
{
   assert(nodenum != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetNodenum", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   *nodenum = scip->stat->nnodes;

   return SCIP_OKAY;
}

RETCODE SCIPgetNNodesLeft(              /**< gets number of nodes left in the tree (children + siblings + leaves) */
   SCIP*            scip,               /**< SCIP data structure */
   int*             nnodes              /**< pointer to store the number of processed nodes */
   )
{
   assert(nnodes != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetNNodesLeft", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *nnodes = SCIPtreeGetNNodes(scip->tree);

   return SCIP_OKAY;
}

RETCODE SCIPgetNLPIterations(           /**< gets number of simplex iterations used so far */
   SCIP*            scip,               /**< SCIP data structure */
   int*             lpiterations        /**< pointer to store the iterations */
   )
{
   assert(lpiterations != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetLPIterations", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *lpiterations = scip->stat->nlpiterations;
   
   return SCIP_OKAY;
}

RETCODE SCIPgetActDepth(                /**< gets depth of active node */
   SCIP*            scip,               /**< SCIP data structure */
   int*             actdepth            /**< pointer to store the depth */
   )
{
   assert(actdepth != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetActDepth", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   if( scip->tree->actnode != NULL )
      *actdepth = scip->tree->actnode->depth;
   else
      *actdepth = -1;
   
   return SCIP_OKAY;
}

RETCODE SCIPgetMaxDepth(                /**< gets maximal depth of all processed nodes */
   SCIP*            scip,               /**< SCIP data structure */
   int*             maxdepth            /**< pointer to store the depth */
   )
{
   assert(maxdepth != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetMaxDepth", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *maxdepth = scip->stat->maxdepth;
   
   return SCIP_OKAY;
}

RETCODE SCIPgetActNCols(                /**< gets number of columns in actual LP */
   SCIP*            scip,               /**< SCIP data structure */
   int*             actncols            /**< pointer to store the number of columns */
   )
{
   assert(actncols != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetActNCols", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *actncols = scip->lp->ncols;
   
   return SCIP_OKAY;
}

RETCODE SCIPgetActNRows(                /**< gets number of rows in actual LP */
   SCIP*            scip,               /**< SCIP data structure */
   int*             actnrows            /**< pointer to store the number of columns */
   )
{
   assert(actnrows != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetActNRows", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *actnrows = scip->lp->nrows;
   
   return SCIP_OKAY;
}

RETCODE SCIPgetActDualBound(            /**< gets dual bound of active node */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            actdualbound        /**< pointer to store the dual bound */
   )
{
   assert(actdualbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetActDualBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *actdualbound = SCIPprobExternObjval(scip->origprob, 
      SCIPprobExternObjval(scip->transprob, SCIPtreeGetActLowerbound(scip->tree)));
   
   return SCIP_OKAY;
}

RETCODE SCIPgetAvgDualBound(            /**< gets average dual bound of all unprocessed nodes */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            avgdualbound        /**< pointer to store the average dual bound */
   )
{
   assert(avgdualbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetAvgDualBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *avgdualbound = SCIPprobExternObjval(scip->origprob, 
      SCIPprobExternObjval(scip->transprob, SCIPtreeGetAvgLowerbound(scip->tree)));
   
   return SCIP_OKAY;
}

RETCODE SCIPgetDualBound(               /**< gets actual dual bound */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            dualbound           /**< pointer to store the dual bound */
   )
{
   assert(dualbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetDualBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *dualbound = SCIPprobExternObjval(scip->origprob, 
      SCIPprobExternObjval(scip->transprob, SCIPtreeGetLowerbound(scip->tree, scip->set)));
   
   return SCIP_OKAY;
}

RETCODE SCIPgetPrimalBound(             /**< gets actual primal bound */
   SCIP*            scip,               /**< SCIP data structure */
   Real*            primalbound         /**< pointer to store the primal bound */
   )
{
   assert(primalbound != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetPrimalBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   *primalbound = SCIPprobExternObjval(scip->origprob, SCIPprobExternObjval(scip->transprob, scip->primal->upperbound));
   
   return SCIP_OKAY;
}

RETCODE SCIPgetBestSol(                 /**< gets best feasible primal solution found so far */
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol                 /**< pointer to store the solution, returns NULL if no solution available */
   )
{
   assert(sol != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetBestSol", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   if( scip->primal->nsols > 0 ) 
   {
      assert(scip->primal->sols != NULL);
      assert(scip->primal->sols[0] != NULL);
      *sol = scip->primal->sols[0];
   }
   else
      *sol = NULL;
   
   return SCIP_OKAY;
}

VERBLEVEL SCIPverbLevel(                /**< gets verbosity level for message output */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->verblevel;
}

RETCODE SCIPsetVerbLevel(               /**< sets verbosity level for message output */
   SCIP*            scip,               /**< SCIP data structure */
   VERBLEVEL        verblevel           /**< verbosity level for message output */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetVerbLevel", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      return scip->mem->probmem;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      return scip->mem->solvemem;

   default:
      errorMessage("Unknown SCIP stage");
      return NULL;
   }
}

int SCIPcalcMemGrowSize(                /**< calculate memory size for dynamically allocated arrays */
   SCIP*            scip,               /**< SCIP data structure */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);

   return SCIPsetCalcMemGrowSize(scip->set, num);
}

Real SCIPinfinity(                      /**< returns value treated as infinity */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->infinity;
}

Real SCIPepsilon(                       /**< returns value treated as zero */
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->epsilon;
}

Bool SCIPisEQ(                          /**< checks, if values are in range of epsZero */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsEQ(scip->set, val1, val2);
}

Bool SCIPisL(                           /**< checks, if val1 is (more than epsZero) lower than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsL(scip->set, val1, val2);
}

Bool SCIPisLE(                          /**< checks, if val1 is not (more than epsZero) greater than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsLE(scip->set, val1, val2);
}

Bool SCIPisG(                           /**< checks, if val1 is (more than epsZero) greater than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsG(scip->set, val1, val2);
}

Bool SCIPisGE(                          /**< checks, if val1 is not (more than epsZero) lower than val2 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsGE(scip->set, val1, val2);
}

Bool SCIPisInfinity(                    /**< checks, if value is (positive) infinite */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against infinity */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsInfinity(scip->set, val);
}

Bool SCIPisZero(                        /**< checks, if value is in range epsZero of 0.0 */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsZero(scip->set, val);
}

Bool SCIPisPos(                         /**< checks, if value is greater than epsZero */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsPos(scip->set, val);
}

Bool SCIPisNeg(                         /**< checks, if value is lower than -epsZero */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsNeg(scip->set, val);
}

Real SCIPfloor(                         /**< rounds value down to the next integer */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFloor(scip->set, val);
}

Real SCIPceil(                          /**< rounds value up to the next integer */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetCeil(scip->set, val);
}

Bool SCIPisIntegral(                    /**< checks, if value is integral within the LP feasibility bounds */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsIntegral(scip->set, val);
}

Bool SCIPisFeasible(                    /**< checks, if value is non-negative within the LP feasibility bounds */
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasible(scip->set, val);
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
   assert(scip->mem != NULL);

   memoryDiagnostic();

   printf("\nProblem Block Memory (%p):\n", scip->mem->probmem);
   blockMemoryDiagnostic(scip->mem->probmem);

   printf("\nSolution Block Memory (%p):\n", scip->mem->solvemem);
   blockMemoryDiagnostic(scip->mem->solvemem);
}

#endif
