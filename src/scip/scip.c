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
#pragma ident "@(#) $Id: scip.c,v 1.228 2004/11/19 14:45:12 bzfpfend Exp $"

/**@file   scip.c
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Kati Wolter
 */
/**@todo check all checkStage() calls, use bit flags instead of the Bool parameters */
/**@todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdarg.h>
#include <assert.h>

#include "def.h"
#include "retcode.h"
#include "set.h"
#include "stat.h"
#include "clock.h"
#include "vbc.h"
#include "interrupt.h"
#include "lpi.h"
#include "mem.h"
#include "misc.h"
#include "history.h"
#include "event.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
#include "sol.h"
#include "primal.h"
#include "tree.h"
#include "pricestore.h"
#include "sepastore.h"
#include "cutpool.h"
#include "solve.h"
#include "scip.h"

#include "branch.h"
#include "conflict.h"
#include "cons.h"
#include "dialog.h"
#include "disp.h"
#include "heur.h"
#include "nodesel.h"
#include "reader.h"
#include "presol.h"
#include "pricer.h"
#include "relax.h"
#include "sepa.h"
#include "prop.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "struct_scip.h"
#endif



#ifndef NDEBUG
/** checks, if SCIP is in one of the feasible stages */
static
RETCODE checkStage(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      method,             /**< method that was called */
   Bool             init,               /**< may method be called in the INIT stage? */
   Bool             problem,            /**< may method be called in the PROBLEM stage? */
   Bool             transforming,       /**< may method be called in the TRANSFORMING stage? */
   Bool             transformed,        /**< may method be called in the TRANSFORMED stage? */
   Bool             presolving,         /**< may method be called in the PRESOLVING stage? */
   Bool             presolved,          /**< may method be called in the PRESOLVED stage? */
   Bool             initsolve,          /**< may method be called in the INITSOLVE stage? */
   Bool             solving,            /**< may method be called in the SOLVING stage? */
   Bool             solved,             /**< may method be called in the SOLVED stage? */
   Bool             freesolve,          /**< may method be called in the FREESOLVE stage? */
   Bool             freetrans           /**< may method be called in the FREETRANS stage? */
   )
{
   assert(scip != NULL);
   assert(method != NULL);

   debugMessage("called method <%s> at stage %d ------------------------------------------------\n", method, scip->stage);

   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->interrupt != NULL);
   assert(scip->dialoghdlr != NULL);
   assert(scip->totaltime != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      assert(scip->stat == NULL);
      assert(scip->origprob == NULL);
      assert(scip->eventfilter == NULL);
      assert(scip->eventqueue == NULL);
      assert(scip->branchcand == NULL);
      assert(scip->lp == NULL);
      assert(scip->primal == NULL);
      assert(scip->transprob == NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->tree == NULL);

      if( !init )
      {
         errorMessage("cannot call method <%s> in initialization stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter == NULL);
      assert(scip->eventqueue == NULL);
      assert(scip->branchcand == NULL);
      assert(scip->lp == NULL);
      assert(scip->primal == NULL);
      assert(scip->transprob == NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->tree == NULL);

      if( !problem )
      {
         errorMessage("cannot call method <%s> in problem creation stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->tree == NULL);

      if( !transforming )
      {
         errorMessage("cannot call method <%s> in problem transformation stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->tree == NULL);

      if( !transformed )
      {
         errorMessage("cannot call method <%s> in problem transformed stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->tree == NULL);

      if( !presolving )
      {
         errorMessage("cannot call method <%s> in presolving stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->tree == NULL);

      if( !presolved )
      {
         errorMessage("cannot call method <%s> in problem presolved stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore != NULL);
      assert(scip->sepastore != NULL);
      assert(scip->cutpool != NULL);
      assert(scip->conflict != NULL);
      assert(scip->tree != NULL);

      if( !solving )
      {
         errorMessage("cannot call method <%s> in solving stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore != NULL);
      assert(scip->sepastore != NULL);
      assert(scip->cutpool != NULL);
      assert(scip->conflict != NULL);
      assert(scip->tree != NULL);

      if( !solved )
      {
         errorMessage("cannot call method <%s> in problem solved stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_FREESOLVE:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->transprob != NULL);

      if( !freesolve )
      {
         errorMessage("cannot call method <%s> in solve deinitialization stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_FREETRANS:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->tree == NULL);

      if( !freetrans )
      {
         errorMessage("cannot call method <%s> in free transformed problem stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }
}
#else
#define checkStage(scip,method,init,problem,transforming,transformed,presolving,presolved,initsolve,solving,solved, \
   freesolve,freetrans) SCIP_OKAY
#endif




/*
 * miscellaneous methods
 */

/** returns scip version number */
Real SCIPversion(
   void
   )
{
   return (Real)(SCIP_VERSION)/100.0;
}

/** prints a version information line to a file stream */
void SCIPprintVersion(
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   if( file == NULL )
      file = stdout;

   fprintf(file, "SCIP version %g [precision: %d byte]", SCIPversion(), (int)sizeof(Real));
#ifndef NOBLOCKMEM
   fprintf(file, " [memory: block]");
#else
   fprintf(file, " [memory: standard]");
#endif
#ifndef NDEBUG
   fprintf(file, " [mode: debug]");
#else
   fprintf(file, " [mode: optimized]");
#endif
   fprintf(file, " [LP solver: %s]\n", SCIPlpiGetSolverName());
}

/** prints error message for the given SCIP return code */
void SCIPprintError(
   RETCODE          retcode,            /**< SCIP return code causing the error */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   fprintf(file, "SCIP Error (%d): ", retcode);
   SCIPretcodePrint(file, retcode);
   fprintf(file, "\n");
}




/*
 * general SCIP methods
 */

/** creates and initializes SCIP data structures */
RETCODE SCIPcreate(
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   ALLOC_OKAY( allocMemory(scip) );

   (*scip)->stage = SCIP_STAGE_INIT;

   CHECK_OKAY( SCIPmemCreate(&(*scip)->mem) );
   CHECK_OKAY( SCIPsetCreate(&(*scip)->set, (*scip)->mem->setmem, *scip) );
   CHECK_OKAY( SCIPinterruptCreate(&(*scip)->interrupt) );
   CHECK_OKAY( SCIPdialoghdlrCreate(&(*scip)->dialoghdlr) );
   CHECK_OKAY( SCIPclockCreate(&(*scip)->totaltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIPclockStart((*scip)->totaltime, (*scip)->set);
   (*scip)->stat = NULL;
   (*scip)->origprob = NULL;
   (*scip)->eventfilter = NULL;
   (*scip)->eventqueue = NULL;
   (*scip)->branchcand = NULL;
   (*scip)->tree = NULL;
   (*scip)->lp = NULL;
   (*scip)->transprob = NULL;
   (*scip)->pricestore = NULL;
   (*scip)->sepastore = NULL;
   (*scip)->cutpool = NULL;
   (*scip)->conflict = NULL;
   (*scip)->primal = NULL;

   return SCIP_OKAY;
}

/** frees SCIP data structures */
RETCODE SCIPfree(
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   CHECK_OKAY( checkStage(*scip, "SCIPfree", TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPfreeProb(*scip) );
   assert((*scip)->stage == SCIP_STAGE_INIT);

   CHECK_OKAY( SCIPsetFree(&(*scip)->set, (*scip)->mem->setmem) );
   CHECK_OKAY( SCIPdialoghdlrFree(&(*scip)->dialoghdlr) );
   SCIPclockFree(&(*scip)->totaltime);
   SCIPinterruptFree(&(*scip)->interrupt);
   CHECK_OKAY( SCIPmemFree(&(*scip)->mem) );

   freeMemory(scip);

   return SCIP_OKAY;
}

/** prints a message depending on the verbosity level */
void SCIPmessage(
   SCIP*            scip,               /**< SCIP data structure */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   )
{
   va_list ap;

   CHECK_ABORT( checkStage(scip, "SCIPmessage", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   va_start(ap, formatstr);
   vinfoMessage(scip->set->disp_verblevel, msgverblevel, formatstr, ap);
   va_end(ap);
}

/** returns current stage of SCIP */
STAGE SCIPgetStage(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->stage;
}

/** outputs SCIP stage and solution status if applicable */
RETCODE SCIPprintStage(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintStage", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( file == NULL )
      file = stdout;

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      fprintf(file, "initialization");
      break;
   case SCIP_STAGE_PROBLEM:
      fprintf(file, "problem creation / modification");
      break;
   case SCIP_STAGE_TRANSFORMING:
      fprintf(file, "problem transformation");
      break;
   case SCIP_STAGE_TRANSFORMED:
      fprintf(file, "problem transformed");
      break;
   case SCIP_STAGE_PRESOLVING:
      fprintf(file, "presolving is running");
      break;
   case SCIP_STAGE_PRESOLVED:
      fprintf(file, "problem is presolved");
      break;
   case SCIP_STAGE_INITSOLVE:
      fprintf(file, "solving process initialization");
      break;
   case SCIP_STAGE_SOLVING:
      if( SCIPsolveIsStopped(scip->set, scip->stat) )
      {
         fprintf(file, "solving was interrupted [");
         CHECK_OKAY( SCIPprintStatus(scip, file) );
         fprintf(file, "]");
      }
      else
         fprintf(file, "solving process is running");
      break;
   case SCIP_STAGE_SOLVED:
      fprintf(file, "problem is solved [");
      CHECK_OKAY( SCIPprintStatus(scip, file) );
      fprintf(file, "]");
      break;
   case SCIP_STAGE_FREESOLVE:
      fprintf(file, "solving process deinitialization");
      break;
   case SCIP_STAGE_FREETRANS:
      fprintf(file, "freeing transformed problem");
      break;
   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** gets solution status */
STATUS SCIPgetStatus(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetStatus", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( scip->stage == SCIP_STAGE_INIT )
      return SCIP_STATUS_UNKNOWN;
   else
   {
      assert(scip->stat != NULL);

      return scip->stat->status;
   }
}

/** outputs solution status */
RETCODE SCIPprintStatus(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintStatus", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( file == NULL )
      file = stdout;

   switch( SCIPgetStatus(scip) )
   {
   case SCIP_STATUS_UNKNOWN:
      fprintf(file, "unknown");
      break;
   case SCIP_STATUS_USERINTERRUPT:
      fprintf(file, "user interrupt");
      break;
   case SCIP_STATUS_NODELIMIT:
      fprintf(file, "node limit reached");
      break;
   case SCIP_STATUS_TIMELIMIT:
      fprintf(file, "time limit reached");
      break;
   case SCIP_STATUS_MEMLIMIT:
      fprintf(file, "memory limit reached");
      break;
   case SCIP_STATUS_GAPLIMIT:
      fprintf(file, "gap limit reached");
      break;
   case SCIP_STATUS_SOLLIMIT:
      fprintf(file, "solution limit reached");
      break;
   case SCIP_STATUS_BESTSOLLIMIT:
      fprintf(file, "solution improvement limit reached");
      break;
   case SCIP_STATUS_OPTIMAL:
      fprintf(file, "optimal solution found");
      break;
   case SCIP_STATUS_INFEASIBLE:
      fprintf(file, "infeasible");
      break;
   case SCIP_STATUS_UNBOUNDED:
      fprintf(file, "unbounded");
      break;
   case SCIP_STATUS_INFORUNBD:
      fprintf(file, "infeasible or unbounded");
      break;
   default:
      errorMessage("invalid status code <%d>\n", SCIPgetStatus(scip));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** returns whether the current stage belongs to the transformed problem space */
Bool SCIPisTransformed(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return ((int)scip->stage >= (int)SCIP_STAGE_TRANSFORMING);
}

/** returns whether the solution process should be provably correct */
Bool SCIPisExactSolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return (scip->set->misc_exactsolve);
}

/** returns whether the user pressed CTRL-C to interrupt the solving process */
Bool SCIPpressedCtrlC(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPpressedCtrlC", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPinterrupted();
}




/*
 * parameter settings
 */

/** creates a Bool parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPaddBoolParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Bool             defaultvalue,       /**< default value of the parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddBoolParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetAddBoolParam(scip->set, scip->mem->setmem, name, desc, valueptr, defaultvalue,
         paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
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
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddIntParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetAddIntParam(scip->set, scip->mem->setmem, name, desc, valueptr, defaultvalue, minvalue, maxvalue,
         paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a Longint parameter, sets it to its default value, and adds it to the parameter set */
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
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddLongintParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetAddLongintParam(scip->set, scip->mem->setmem, name, desc,
         valueptr, defaultvalue, minvalue, maxvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a Real parameter, sets it to its default value, and adds it to the parameter set */
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
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddRealParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetAddRealParam(scip->set, scip->mem->setmem, name, desc, valueptr, defaultvalue, minvalue, maxvalue,
         paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPaddCharParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   char*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   char             defaultvalue,       /**< default value of the parameter */
   const char*      allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddCharParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetAddCharParam(scip->set, scip->mem->setmem, name, desc, valueptr, defaultvalue, allowedvalues,
         paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
RETCODE SCIPaddStringParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   char**           valueptr,           /**< pointer to store the current parameter value, or NULL */
   const char*      defaultvalue,       /**< default value of the parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddStringParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetAddStringParam(scip->set, scip->mem->setmem, name, desc, valueptr, defaultvalue,
         paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** gets the value of an existing Bool parameter */
RETCODE SCIPgetBoolParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Bool*            value               /**< pointer to store the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetBoolParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetGetBoolParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing Int parameter */
RETCODE SCIPgetIntParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   int*             value               /**< pointer to store the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetIntParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetGetIntParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing Longint parameter */
RETCODE SCIPgetLongintParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Longint*         value               /**< pointer to store the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLongintParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetGetLongintParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing Real parameter */
RETCODE SCIPgetRealParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Real*            value               /**< pointer to store the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetRealParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetGetRealParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing Char parameter */
RETCODE SCIPgetCharParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   char*            value               /**< pointer to store the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetCharParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetGetCharParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing String parameter */
RETCODE SCIPgetStringParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   char**           value               /**< pointer to store the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetStringParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetGetStringParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing Bool parameter */
RETCODE SCIPsetBoolParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Bool             value               /**< new value of the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetBoolParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetSetBoolParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing Int parameter */
RETCODE SCIPsetIntParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   int              value               /**< new value of the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetIntParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetSetIntParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing Longint parameter */
RETCODE SCIPsetLongintParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Longint          value               /**< new value of the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetLongintParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetSetLongintParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing Real parameter */
RETCODE SCIPsetRealParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   Real             value               /**< new value of the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetRealParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetSetRealParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing Char parameter */
RETCODE SCIPsetCharParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   char             value               /**< new value of the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetCharParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetSetCharParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing String parameter */
RETCODE SCIPsetStringParam(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of the parameter */
   const char*      value               /**< new value of the parameter */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetStringParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetSetStringParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** reads parameters from a file */
RETCODE SCIPreadParams(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< file name */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPreadParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetReadParams(scip->set, filename) );

   return SCIP_OKAY;
}

/** writes all parameters in the parameter set to a file */
RETCODE SCIPwriteParams(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename,           /**< file name, or NULL for stdout */
   Bool             comments,           /**< should parameter descriptions be written as comments? */
   Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPwriteParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetWriteParams(scip->set, filename, comments, onlychanged) );

   return SCIP_OKAY;
}

/** returns the array of all available SCIP parameters */
PARAM** SCIPgetParams(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPsetGetParams(scip->set);
}

/** returns the total number of all available SCIP parameters */
int SCIPgetNParams(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPsetGetNParams(scip->set);
}




/*
 * SCIP user functionality methods: managing plugins
 */

/** creates a reader and includes it in SCIP */
RETCODE SCIPincludeReader(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of reader */
   const char*      desc,               /**< description of reader */
   const char*      extension,          /**< file extension that reader processes */
   DECL_READERFREE  ((*readerfree)),    /**< destructor of reader */
   DECL_READERREAD  ((*readerread)),    /**< read method */
   READERDATA*      readerdata          /**< reader data */
   )
{
   READER* reader;

   CHECK_OKAY( checkStage(scip, "SCIPincludeReader", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPreaderCreate(&reader, name, desc, extension, readerfree, readerread, readerdata) );
   CHECK_OKAY( SCIPsetIncludeReader(scip->set, reader) );
   
   return SCIP_OKAY;
}

/** returns the reader of the given name, or NULL if not existing */
READER* SCIPfindReader(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of constraint handler */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindReader", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindReader(scip->set, name);
}

/** returns the array of currently available readers */
READER** SCIPgetReaders(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetReaders", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->readers;
}

/** returns the number of currently available readers */
int SCIPgetNReaders(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNReaders", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nreaders;
}

/** creates a variable pricer and includes it in SCIP
 *  To use the variable pricer for solving a problem, it first has to be activated with a call to SCIPactivatePricer().
 *  This should be done during the problem creation stage.
 */
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
   )
{
   PRICER* pricer;

   CHECK_OKAY( checkStage(scip, "SCIPincludePricer", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPpricerCreate(&pricer, scip->set, scip->mem->setmem,
         name, desc, priority,
         pricerfree, pricerinit, pricerexit, pricerredcost, pricerfarkas, pricerdata) );
   CHECK_OKAY( SCIPsetIncludePricer(scip->set, pricer) );
   
   return SCIP_OKAY;
}

/** returns the variable pricer of the given name, or NULL if not existing */
PRICER* SCIPfindPricer(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of variable pricer */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindPricer", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindPricer(scip->set, name);
}

/** returns the array of currently available variable pricers; active pricers are in the first slots of the array */
PRICER** SCIPgetPricers(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPricers", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortPricers(scip->set);

   return scip->set->pricers;
}

/** returns the number of currently available variable pricers */
int SCIPgetNPricers(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPricers", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->npricers;
}

/** returns the number of currently active variable pricers, that are used in the LP solving loop */
int SCIPgetNActivePricers(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNAcvitePricers", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nactivepricers;
}

/** sets the priority of a variable pricer */
RETCODE SCIPsetPricerPriority(
   SCIP*            scip,               /**< SCIP data structure */
   PRICER*          pricer,             /**< variable pricer */
   int              priority            /**< new priority of the variable pricer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetPricerPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPpricerSetPriority(pricer, scip->set, priority);

   return SCIP_OKAY;
}

/** activates pricer to be used for the current problem
 *  This method should be called during the problem creation stage for all pricers that are necessary to solve
 *  the problem model.
 *  The pricers are automatically deactivated when the problem is freed.
 */
RETCODE SCIPactivatePricer(
   SCIP*            scip,               /**< SCIP data structure */
   PRICER*          pricer              /**< variable pricer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPactivatePricer", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPpricerActivate(pricer, scip->set);

   return SCIP_OKAY;
}

/** creates a constraint handler and includes it in SCIP */
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
   int              maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
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
   DECL_CONSRESPROP ((*consresprop)),   /**< propagation conflict resolving method */
   DECL_CONSLOCK    ((*conslock)),      /**< variable rounding lock method */
   DECL_CONSUNLOCK  ((*consunlock)),    /**< variable rounding unlock method */
   DECL_CONSACTIVE  ((*consactive)),    /**< activation notification method */
   DECL_CONSDEACTIVE((*consdeactive)),  /**< deactivation notification method */
   DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   DECL_CONSPRINT   ((*consprint)),     /**< constraint display method */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   CONSHDLR* conshdlr;

   CHECK_OKAY( checkStage(scip, "SCIPincludeConshdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconshdlrCreate(&conshdlr, scip->set, scip->mem->setmem,
         name, desc, sepapriority, enfopriority, chckpriority, sepafreq, propfreq, eagerfreq, maxprerounds, needscons, 
         consfree, consinit, consexit, consinitpre, consexitpre, consinitsol, consexitsol, 
         consdelete, constrans, consinitlp, conssepa, consenfolp, consenfops, conscheck, consprop, conspresol,
         consresprop, conslock, consunlock, consactive, consdeactive, consenable, consdisable, consprint,
         conshdlrdata) );
   CHECK_OKAY( SCIPsetIncludeConshdlr(scip->set, conshdlr) );
   
   return SCIP_OKAY;
}

/** returns the constraint handler of the given name, or NULL if not existing */
CONSHDLR* SCIPfindConshdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of constraint handler */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindConshdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindConshdlr(scip->set, name);
}

/** returns the array of currently available constraint handlers */
CONSHDLR** SCIPgetConshdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetConshdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->conshdlrs;
}

/** returns the number of currently available constraint handlers */
int SCIPgetNConshdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNConshdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nconshdlrs;
}

/** creates a conflict handler and includes it in SCIP */
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
   )
{
   CONFLICTHDLR* conflicthdlr;

   CHECK_OKAY( checkStage(scip, "SCIPincludeConflicthdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconflicthdlrCreate(&conflicthdlr, scip->set, scip->mem->setmem, name, desc, priority, 
         conflictfree, conflictinit, conflictexit, conflictexec,
         conflicthdlrdata) );
   CHECK_OKAY( SCIPsetIncludeConflicthdlr(scip->set, conflicthdlr) );
   
   return SCIP_OKAY;
}

/** returns the conflict handler of the given name, or NULL if not existing */
CONFLICTHDLR* SCIPfindConflicthdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of conflict handler */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindConflicthdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindConflicthdlr(scip->set, name);
}

/** returns the array of currently available conflict handlers */
CONFLICTHDLR** SCIPgetConflicthdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetConflicthdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortConflicthdlrs(scip->set);

   return scip->set->conflicthdlrs;
}

/** returns the number of currently available conflict handlers */
int SCIPgetNConflicthdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNConflicthdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nconflicthdlrs;
}

/** sets the priority of a conflict handler */
RETCODE SCIPsetConflicthdlrPriority(
   SCIP*            scip,               /**< SCIP data structure */
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   int              priority            /**< new priority of the conflict handler */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetConflicthdlrPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPconflicthdlrSetPriority(conflicthdlr, scip->set, priority);

   return SCIP_OKAY;
}

/** creates a presolver and includes it in SCIP */
RETCODE SCIPincludePresol(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of presolver */
   const char*      desc,               /**< description of presolver */
   int              priority,           /**< priority of the presolver */
   int              maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
   DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver to free user data (called when SCIP is exiting) */
   DECL_PRESOLINIT  ((*presolinit)),    /**< initialization method of presolver (called after problem was transformed) */
   DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialization method of presolver (called before transformed problem is freed) */
   DECL_PRESOLINITPRE((*presolinitpre)),/**< presolving initialization method of presolver (called when presolving is about to begin) */
   DECL_PRESOLEXITPRE((*presolexitpre)),/**< presolving deinitialization method of presolver (called after presolving has been finished) */
   DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   PRESOL* presol;

   CHECK_OKAY( checkStage(scip, "SCIPincludePresol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPpresolCreate(&presol, scip->set, scip->mem->setmem, name, desc, priority, maxrounds,
         presolfree, presolinit, presolexit, presolinitpre, presolexitpre, presolexec, presoldata) );
   CHECK_OKAY( SCIPsetIncludePresol(scip->set, presol) );
   
   return SCIP_OKAY;
}

/** returns the presolver of the given name, or NULL if not existing */
PRESOL* SCIPfindPresol(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of presolver */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindPresol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindPresol(scip->set, name);
}

/** returns the array of currently available presolvers */
PRESOL** SCIPgetPresols(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPresols", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortPresols(scip->set);

   return scip->set->presols;
}

/** returns the number of currently available presolvers */
int SCIPgetNPresols(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPresols", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->npresols;
}

/** sets the priority of a presolver */
RETCODE SCIPsetPresolPriority(
   SCIP*            scip,               /**< SCIP data structure */
   PRESOL*          presol,             /**< presolver */
   int              priority            /**< new priority of the presolver */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetPresolPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPpresolSetPriority(presol, scip->set, priority);

   return SCIP_OKAY;
}

/** creates a relaxator and includes it in SCIP */
RETCODE SCIPincludeRelax(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of relaxator */
   const char*      desc,               /**< description of relaxator */
   int              priority,           /**< priority of the relaxator (negative: after LP, non-negative: before LP) */
   int              freq,               /**< frequency for calling relaxator */
   DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxator */
   DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxator */
   DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxator */
   DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxator */
   RELAXDATA*       relaxdata           /**< relaxator data */
   )
{
   RELAX* relax;

   CHECK_OKAY( checkStage(scip, "SCIPincludeRelax", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPrelaxCreate(&relax, scip->set, scip->mem->setmem,
         name, desc, priority, freq,
         relaxfree, relaxinit, relaxexit, relaxexec, relaxdata) );
   CHECK_OKAY( SCIPsetIncludeRelax(scip->set, relax) );
   
   return SCIP_OKAY;
}

/** returns the relaxator of the given name, or NULL if not existing */
RELAX* SCIPfindRelax(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of relaxator */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindRelax", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindRelax(scip->set, name);
}

/** returns the array of currently available relaxators */
RELAX** SCIPgetRelaxs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRelaxs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortRelaxs(scip->set);

   return scip->set->relaxs;
}

/** returns the number of currently available relaxators */
int SCIPgetNRelaxs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNRelaxs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nrelaxs;
}

/** sets the priority of a relaxator */
RETCODE SCIPsetRelaxPriority(
   SCIP*            scip,               /**< SCIP data structure */
   RELAX*           relax,              /**< primal relaxistic */
   int              priority            /**< new priority of the relaxator */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetRelaxPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPrelaxSetPriority(relax, scip->set, priority);

   return SCIP_OKAY;
}

/** creates a separator and includes it in SCIP */
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
   )
{
   SEPA* sepa;

   CHECK_OKAY( checkStage(scip, "SCIPincludeSepa", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsepaCreate(&sepa, scip->set, scip->mem->setmem,
         name, desc, priority, freq,
         sepafree, sepainit, sepaexit, sepaexec, sepadata) );
   CHECK_OKAY( SCIPsetIncludeSepa(scip->set, sepa) );
   
   return SCIP_OKAY;
}

/** returns the separator of the given name, or NULL if not existing */
SEPA* SCIPfindSepa(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of separator */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindSepa", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindSepa(scip->set, name);
}

/** returns the array of currently available separators */
SEPA** SCIPgetSepas(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSepas", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortSepas(scip->set);

   return scip->set->sepas;
}

/** returns the number of currently available separators */
int SCIPgetNSepas(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNSepas", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nsepas;
}

/** sets the priority of a separator */
RETCODE SCIPsetSepaPriority(
   SCIP*            scip,               /**< SCIP data structure */
   SEPA*            sepa,               /**< primal sepaistic */
   int              priority            /**< new priority of the separator */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetSepaPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsepaSetPriority(sepa, scip->set, priority);

   return SCIP_OKAY;
}

/** creates a propagator and includes it in SCIP */
RETCODE SCIPincludeProp(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of propagator */
   const char*      desc,               /**< description of propagator */
   int              priority,           /**< priority of the propagator */
   int              freq,               /**< frequency for calling propagator */
   DECL_PROPFREE    ((*propfree)),      /**< destructor of propagator */
   DECL_PROPINIT    ((*propinit)),      /**< initialize propagator */
   DECL_PROPEXIT    ((*propexit)),      /**< deinitialize propagator */
   DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   DECL_PROPRESPROP ((*propresprop)),   /**< propagation conflict resolving method */
   PROPDATA*        propdata            /**< propagator data */
   )
{
   PROP* prop;

   CHECK_OKAY( checkStage(scip, "SCIPincludeProp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPpropCreate(&prop, scip->set, scip->mem->setmem,
         name, desc, priority, freq,
         propfree, propinit, propexit, propexec, propresprop, propdata) );
   CHECK_OKAY( SCIPsetIncludeProp(scip->set, prop) );
   
   return SCIP_OKAY;
}

/** returns the propagator of the given name, or NULL if not existing */
PROP* SCIPfindProp(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of propagator */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindProp", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindProp(scip->set, name);
}

/** returns the array of currently available propagators */
PROP** SCIPgetProps(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetProps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortProps(scip->set);

   return scip->set->props;
}

/** returns the number of currently available propagators */
int SCIPgetNProps(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNProps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nprops;
}

/** sets the priority of a propagator */
RETCODE SCIPsetPropPriority(
   SCIP*            scip,               /**< SCIP data structure */
   PROP*            prop,               /**< primal propistic */
   int              priority            /**< new priority of the propagator */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetPropPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPpropSetPriority(prop, scip->set, priority);

   return SCIP_OKAY;
}

/** creates a primal heuristic and includes it in SCIP */
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
   )
{
   HEUR* heur;

   CHECK_OKAY( checkStage(scip, "SCIPincludeHeur", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPheurCreate(&heur, scip->set, scip->mem->setmem,
         name, desc, dispchar, priority, freq, freqofs, maxdepth, pseudonodes, duringplunging,
         heurfree, heurinit, heurexit, heurexec, heurdata) );
   CHECK_OKAY( SCIPsetIncludeHeur(scip->set, heur) );
   
   return SCIP_OKAY;
}

/** returns the primal heuristic of the given name, or NULL if not existing */
HEUR* SCIPfindHeur(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of primal heuristic */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindHeur", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindHeur(scip->set, name);
}

/** returns the array of currently available primal heuristics */
HEUR** SCIPgetHeurs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetHeurs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortHeurs(scip->set);

   return scip->set->heurs;
}

/** returns the number of currently available primal heuristics */
int SCIPgetNHeurs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNHeurs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nheurs;
}

/** sets the priority of a primal heuristic */
RETCODE SCIPsetHeurPriority(
   SCIP*            scip,               /**< SCIP data structure */
   HEUR*            heur,               /**< primal heuristic */
   int              priority            /**< new priority of the primal heuristic */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetHeurPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPheurSetPriority(heur, scip->set, priority);

   return SCIP_OKAY;
}

/** creates an event handler and includes it in SCIP */
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
   )
{
   EVENTHDLR* eventhdlr;

   CHECK_OKAY( checkStage(scip, "SCIPincludeEventhdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPeventhdlrCreate(&eventhdlr, name, desc, 
         eventfree, eventinit, eventexit, eventdelete, eventexec,
         eventhdlrdata) );
   CHECK_OKAY( SCIPsetIncludeEventhdlr(scip->set, eventhdlr) );
   
   return SCIP_OKAY;
}

/** returns the event handler of the given name, or NULL if not existing */
EVENTHDLR* SCIPfindEventhdlr(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindEventhdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindEventhdlr(scip->set, name);
}

/** returns the array of currently available event handlers */
EVENTHDLR** SCIPgetEventhdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetEventhdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->eventhdlrs;
}

/** returns the number of currently available event handlers */
int SCIPgetNEventhdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNEventhdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->neventhdlrs;
}

/** creates a node selector and includes it in SCIP */
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
   )
{
   NODESEL* nodesel;

   CHECK_OKAY( checkStage(scip, "SCIPincludeNodesel", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPnodeselCreate(&nodesel, scip->set, scip->mem->setmem, name, desc, stdpriority, memsavepriority, 
         lowestboundfirst, nodeselfree, nodeselinit, nodeselexit, nodeselselect, nodeselcomp, nodeseldata) );
   CHECK_OKAY( SCIPsetIncludeNodesel(scip->set, nodesel) );
   
   return SCIP_OKAY;
}

/** returns the node selector of the given name, or NULL if not existing */
NODESEL* SCIPfindNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindNodesel", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindNodesel(scip->set, name);
}

/** returns the array of currently available node selectors */
NODESEL** SCIPgetNodesels(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNodesels", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nodesels;
}

/** returns the number of currently available node selectors */
int SCIPgetNNodesels(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNNodesels", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nnodesels;
}

/** sets the priority of a node selector in standard mode */
RETCODE SCIPsetNodeselStdPriority(
   SCIP*            scip,               /**< SCIP data structure */
   NODESEL*         nodesel,            /**< node selector */
   int              priority            /**< new standard priority of the node selector */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetNodeselStdPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPnodeselSetStdPriority(nodesel, scip->set, priority);

   return SCIP_OKAY;
}

/** sets the priority of a node selector in memory saving mode */
RETCODE SCIPsetNodeselMemsavePriority(
   SCIP*            scip,               /**< SCIP data structure */
   NODESEL*         nodesel,            /**< node selector */
   int              priority            /**< new memory saving priority of the node selector */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetNodeselMemsavePriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPnodeselSetMemsavePriority(nodesel, scip->set, priority);

   return SCIP_OKAY;
}

/** returns the currently used node selector */
NODESEL* SCIPgetNodesel(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNodesel", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetGetNodesel(scip->set, scip->stat);
}

/** creates a branching rule and includes it in SCIP */
RETCODE SCIPincludeBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   int              maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
                                         *   compared to best node's dual bound for applying branching rule
                                         *   (0.0: only on current best node, 1.0: on all nodes) */
   DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   DECL_BRANCHINIT  ((*branchinit)),    /**< initialize branching rule */
   DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialize branching rule */
   DECL_BRANCHEXECLP((*branchexeclp)),  /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXECPS((*branchexecps)),  /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   BRANCHRULE* branchrule;

   CHECK_OKAY( checkStage(scip, "SCIPincludeBranchrule", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchruleCreate(&branchrule, scip->mem->setmem, scip->set, name, desc, priority, maxdepth, 
         maxbounddist, branchfree, branchinit, branchexit, branchexeclp, branchexecps, branchruledata) );
   CHECK_OKAY( SCIPsetIncludeBranchrule(scip->set, branchrule) );
   
   return SCIP_OKAY;
}

/** returns the branching rule of the given name, or NULL if not existing */
BRANCHRULE* SCIPfindBranchrule(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindBranchrule", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortBranchrules(scip->set);

   return SCIPsetFindBranchrule(scip->set, name);
}

/** returns the array of currently available branching rules */
BRANCHRULE** SCIPgetBranchrules(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBranchrules", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->branchrules;
}

/** returns the number of currently available branching rules */
int SCIPgetNBranchrules(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNBranchrules", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nbranchrules;
}

/** sets the priority of a branching rule */
RETCODE SCIPsetBranchrulePriority(
   SCIP*            scip,               /**< SCIP data structure */
   BRANCHRULE*      branchrule,         /**< branching rule */
   int              priority            /**< new priority of the branching rule */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetBranchrulePriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPbranchruleSetPriority(branchrule, scip->set, priority);

   return SCIP_OKAY;
}

/** sets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
RETCODE SCIPsetBranchruleMaxdepth(
   SCIP*            scip,               /**< SCIP data structure */
   BRANCHRULE*      branchrule,         /**< branching rule */
   int              maxdepth            /**< new maxdepth of the branching rule */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetBranchruleMaxdepth", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPbranchruleSetMaxdepth(branchrule, maxdepth);

   return SCIP_OKAY;
}

/** sets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
RETCODE SCIPsetBranchruleMaxbounddist(
   SCIP*            scip,               /**< SCIP data structure */
   BRANCHRULE*      branchrule,         /**< branching rule */
   Real             maxbounddist        /**< new maxbounddist of the branching rule */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetBranchruleMaxbounddist", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPbranchruleSetMaxbounddist(branchrule, maxbounddist);

   return SCIP_OKAY;
}

/** creates a display column and includes it in SCIP */
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
   )
{
   DISP* disp;

   CHECK_OKAY( checkStage(scip, "SCIPincludeDisp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPdispCreate(&disp, scip->set, scip->mem->setmem,
         name, desc, header, dispstatus, dispfree, dispinit, dispexit, dispoutput, dispdata,
         width, priority, position, stripline) );
   CHECK_OKAY( SCIPsetIncludeDisp(scip->set, disp) );
   
   return SCIP_OKAY;
}

/** returns the display column of the given name, or NULL if not existing */
DISP* SCIPfindDisp(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of event handler */
   )
{
   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindDisp", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindDisp(scip->set, name);
}

/** returns the array of currently available display columns */
DISP** SCIPgetDisps(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetDisps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->disps;
}

/** returns the number of currently available display columns */
int SCIPgetNDisps(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNDisps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->ndisps;
}

/** automatically selects display columns for being shown w.r.t. the display width parameter */
RETCODE SCIPautoselectDisps(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPselectDisps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPdispAutoActivate(scip->set) );

   return SCIP_OKAY;
}




/*
 * user interactive dialog methods
 */

/**@name User Interactive Dialog Methods */
/**@{ */

/** creates and captures a dialog */
RETCODE SCIPcreateDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG**         dialog,             /**< pointer to store the dialog */
   DECL_DIALOGEXEC  ((*dialogexec)),    /**< execution method of dialog */
   DECL_DIALOGDESC  ((*dialogdesc)),    /**< description output method of dialog, or NULL */
   const char*      name,               /**< name of dialog: command name appearing in parent's dialog menu */
   const char*      desc,               /**< description of dialog used if description output method is NULL */
   Bool             issubmenu,          /**< is the dialog a submenu? */
   DIALOGDATA*      dialogdata          /**< user defined dialog data */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPdialogCreate(dialog, dialogexec, dialogdesc, name, desc, issubmenu, dialogdata) );

   return SCIP_OKAY;
}

/** captures a dialog */
RETCODE SCIPcaptureDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          dialog              /**< dialog */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPdialogCapture(dialog);

   return SCIP_OKAY;
}

/** releases a dialog */
RETCODE SCIPreleaseDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG**         dialog              /**< pointer to the dialog */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPreleaseDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPdialogRelease(dialog) );

   return SCIP_OKAY;
}

/** makes given dialog the root dialog of SCIP's interactive user shell; captures dialog and releases former root dialog */
RETCODE SCIPsetRootDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          dialog              /**< dialog to be the root */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetRootDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPdialoghdlrSetRoot(scip->dialoghdlr, dialog) );

   return SCIP_OKAY;
}

/** returns the root dialog of SCIP's interactive user shell */
DIALOG* SCIPgetRootDialog(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRootDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPdialoghdlrGetRoot(scip->dialoghdlr);
}

/** adds a sub dialog to the given dialog as menu entry and captures it */
RETCODE SCIPaddDialogEntry(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          dialog,             /**< dialog to extend, or NULL for root dialog */
   DIALOG*          subdialog           /**< subdialog to add as menu entry in dialog */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddDialogEntry", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( dialog == NULL )
      dialog = SCIPdialoghdlrGetRoot(scip->dialoghdlr);

   CHECK_OKAY( SCIPdialogAddEntry(dialog, scip->set, subdialog) );

   return SCIP_OKAY;
}

/** starts interactive mode of SCIP by executing the root dialog */
RETCODE SCIPstartInteraction(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPstartInteraction", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPdialoghdlrExec(scip->dialoghdlr, scip) );

   return SCIP_OKAY;
}

/**@} */




/*
 * global problem methods
 */

/** creates empty problem and initializes all solving data structures (the objective sense is set to MINIMIZE)
 *  If the problem type requires the use of variable pricers, these pricers should be added to the problem with calls
 *  to SCIPactivatePricer(). These pricers are automatically deactivated, when the problem is freed.
 */
RETCODE SCIPcreateProb(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< problem name */
   DECL_PROBDELORIG ((*probdelorig)),   /**< frees user data of original problem */
   DECL_PROBTRANS   ((*probtrans)),     /**< creates user data of transformed problem by transforming original user data */
   DECL_PROBDELTRANS((*probdeltrans)),  /**< frees user data of transformed problem */
   DECL_PROBINITSOL ((*probinitsol)),   /**< solving process initialization method of transformed data */
   DECL_PROBEXITSOL ((*probexitsol)),   /**< solving process deinitialization method of transformed data */
   PROBDATA*        probdata            /**< user problem data set by the reader */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateProb", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* switch stage to PROBLEM */
   scip->stage = SCIP_STAGE_PROBLEM;
   
   CHECK_OKAY( SCIPstatCreate(&scip->stat, scip->mem->probmem, scip->set) );
   CHECK_OKAY( SCIPprobCreate(&scip->origprob, scip->mem->probmem, name, 
         probdelorig, probtrans, probdeltrans, probinitsol, probexitsol, probdata, FALSE) );

   return SCIP_OKAY;
}

/** reads problem from file and initializes all solving data structures */
RETCODE SCIPreadProb(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< problem file name */
   )
{
   RESULT result;
   int i;

   assert(filename != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreadProb", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* try all readers until one could read the file */
   result = SCIP_DIDNOTRUN;
   for( i = 0; i < scip->set->nreaders && result == SCIP_DIDNOTRUN; ++i )
   {
      CHECK_OKAY( SCIPreaderRead(scip->set->readers[i], scip->set, filename, &result) );
   }

   switch( result )
   {
   case SCIP_DIDNOTRUN:
      warningMessage("No reader for input file <%s> available\n", filename);
      return SCIP_READERROR;

   case SCIP_SUCCESS:
      if( scip->origprob != NULL )
      {      
         infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
            "original problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
            scip->origprob->nvars, scip->origprob->nbinvars, scip->origprob->nintvars,
            scip->origprob->nimplvars, scip->origprob->ncontvars,
            scip->origprob->nconss);
#if 0
         printf(" var names :  ");
         SCIPhashtablePrintStatistics(scip->origprob->varnames);
         printf(" cons names:  ");
         SCIPhashtablePrintStatistics(scip->origprob->consnames);
#endif
      }
      return SCIP_OKAY;

   default:
      assert(i < scip->set->nreaders);
      errorMessage("invalid result code <%d> from reader <%s> reading file <%s>\n", 
         result, SCIPreaderGetName(scip->set->readers[i]), filename);
      return SCIP_READERROR;
   }  /*lint !e788*/
}

/** frees problem and solution process data */
RETCODE SCIPfreeProb(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeProb", TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPfreeTransform(scip) );
   assert(scip->stage == SCIP_STAGE_INIT || scip->stage == SCIP_STAGE_PROBLEM);

   if( scip->stage == SCIP_STAGE_PROBLEM )
   {
      int p;

      /* deactivate all pricers */
      for( p = 0; p < scip->set->nactivepricers; ++p )
      {
         CHECK_OKAY( SCIPpricerDeactivate(scip->set->pricers[p], scip->set) );
      }
      assert(scip->set->nactivepricers == 0);

      /* free problem and problem statistics datastructures */
      CHECK_OKAY( SCIPprobFree(&scip->origprob, scip->mem->probmem, scip->set, scip->stat, scip->lp) );
      CHECK_OKAY( SCIPstatFree(&scip->stat, scip->mem->probmem) );

      /* switch stage to INIT */
      scip->stage = SCIP_STAGE_INIT;
   }
   assert(scip->stage == SCIP_STAGE_INIT);

   return SCIP_OKAY;
}

/** gets user problem data */
PROBDATA* SCIPgetProbData(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetProbData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobGetData(scip->origprob);

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      return SCIPprobGetData(scip->transprob);

   default:
      errorMessage("invalid SCIP stage\n");
      abort();
   }  /*lint !e788*/
}

/** sets user problem data */
RETCODE SCIPsetProbData(
   SCIP*            scip,               /**< SCIP data structure */
   PROBDATA*        probdata            /**< user problem data to use */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetProbData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetData(scip->origprob, probdata);
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      SCIPprobSetData(scip->transprob, probdata);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets objective sense of original problem */
OBJSENSE SCIPgetObjsense(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetObjsense", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->objsense;
}

/** sets objective sense of problem */
RETCODE SCIPsetObjsense(
   SCIP*            scip,               /**< SCIP data structure */
   OBJSENSE         objsense            /**< new objective sense */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetObjsense", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   if( objsense != SCIP_OBJSENSE_MAXIMIZE && objsense != SCIP_OBJSENSE_MINIMIZE )
   {
      errorMessage("invalid objective sense\n");
      return SCIP_INVALIDDATA;
   }

   SCIPprobSetObjsense(scip->origprob, objsense);
   
   return SCIP_OKAY;
}

/** sets limit on objective function, such that only solutions better than this limit are accepted */
RETCODE SCIPsetObjlimit(
   SCIP*            scip,               /**< SCIP data structure */
   Real             objlimit            /**< new primal objective limit */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetObjlimit", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetExternObjlim(scip->origprob, objlimit);
      break;
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      if( SCIPtransformObj(scip, objlimit) > SCIPprobGetInternObjlim(scip->transprob, scip->set) )
      {
         errorMessage("cannot relax objective limit from %g to %g after problem was transformed\n",
            SCIPprobGetExternObjlim(scip->origprob), objlimit);
         return SCIP_INVALIDDATA;
      }
      SCIPprobSetExternObjlim(scip->origprob, objlimit);
      SCIPprobSetExternObjlim(scip->transprob, SCIPprobGetInternObjlim(scip->origprob, scip->set));
      CHECK_OKAY( SCIPprimalUpdateUpperbound(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
            scip->tree, scip->lp) );
      break;
   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** gets current limit on objective function */
Real SCIPgetObjlimit(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetObjlimit", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPprobGetExternObjlim(scip->origprob);
}

/** informs SCIP, that the objective value is always integral in every feasible solution */
RETCODE SCIPsetObjIntegral(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetObjIntegral", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetObjIntegral(scip->origprob);
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIPprobSetObjIntegral(scip->transprob);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** returns whether the objective value is known to be integral in every feasible solution */
Bool SCIPisObjIntegral(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPisObjIntegral", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobIsObjIntegral(scip->origprob);
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIPprobIsObjIntegral(scip->transprob);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** adds variable to the problem */
RETCODE SCIPaddVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to add */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddVar", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* avoid inserting the same variable twice */
   if( SCIPvarGetProbindex(var) != -1 )
      return SCIP_OKAY;

   /* insert the negation variable x instead of the negated variable x' in x' = offset - x */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
   {
      assert(SCIPvarGetNegationVar(var) != NULL);
      CHECK_OKAY( SCIPaddVar(scip, SCIPvarGetNegationVar(var)) );
      return SCIP_OKAY;
   }

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL )
      {
         errorMessage("cannot add transformed variables to original problem\n");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->origprob, scip->mem->probmem, scip->set, scip->lp, scip->branchcand, 
            scip->eventfilter, scip->eventqueue, var) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
         {
            errorMessage("cannot add original variables to transformed problem\n");
         }
         else
         {
            errorMessage("cannot add fixed or aggregated variables to transformed problem\n");
         }
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->transprob, scip->mem->solvemem, scip->set, scip->lp, 
            scip->branchcand, scip->eventfilter, scip->eventqueue, var) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** adds variable to the problem and uses it as pricing candidate to enter the LP */
RETCODE SCIPaddPricedVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to add */
   Real             score               /**< pricing score of variable (the larger, the better the variable) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddPricedVar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* insert the negation variable x instead of the negated variable x' in x' = offset - x */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
   {
      assert(SCIPvarGetNegationVar(var) != NULL);
      CHECK_OKAY( SCIPaddPricedVar(scip, SCIPvarGetNegationVar(var), score) );
      return SCIP_OKAY;
   }

   /* add variable to problem if not yet inserted */
   if( SCIPvarGetProbindex(var) == -1 )
   {
      /* check variable's status */
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
         {
            errorMessage("cannot add original variables to transformed problem\n");
         }
         else
         {
            errorMessage("cannot add fixed or aggregated variables to transformed problem\n");
         }
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->transprob, scip->mem->solvemem, scip->set, scip->lp, 
            scip->branchcand, scip->eventfilter, scip->eventqueue, var) );
   }

   /* add variable to pricing storage */
   CHECK_OKAY( SCIPpricestoreAddVar(scip->pricestore, scip->mem->solvemem, scip->set, scip->lp, var, score,
         (SCIPtreeGetCurrentDepth(scip->tree) == 0)) );
   
   return SCIP_OKAY;
}

/** gets variables of the problem along with the numbers of different variable types; data may become invalid after
 *  calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 */
RETCODE SCIPgetVarsData(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*             nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*             nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*             nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*             nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*             ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetVarsData", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( vars != NULL )
         *vars = scip->origprob->vars;
      if( nvars != NULL )
         *nvars = scip->origprob->nvars;
      if( nbinvars != NULL )
         *nbinvars = scip->origprob->nbinvars;
      if( nintvars != NULL )
         *nintvars = scip->origprob->nintvars;
      if( nimplvars != NULL )
         *nimplvars = scip->origprob->nimplvars;
      if( ncontvars != NULL )
         *ncontvars = scip->origprob->ncontvars;
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      if( vars != NULL )
         *vars = scip->transprob->vars;
      if( nvars != NULL )
         *nvars = scip->transprob->nvars;
      if( nbinvars != NULL )
         *nbinvars = scip->transprob->nbinvars;
      if( nintvars != NULL )
         *nintvars = scip->transprob->nintvars;
      if( nimplvars != NULL )
         *nimplvars = scip->transprob->nimplvars;
      if( ncontvars != NULL )
         *ncontvars = scip->transprob->ncontvars;
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets array with active problem variables; data may become invalid after
 *  calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 */
VAR** SCIPgetVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->vars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->vars;

   default:
      errorMessage("invalid SCIP stage\n");
      abort();
   }  /*lint !e788*/
}

/** gets number of active problem variables */
int SCIPgetNVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nvars;

   default:
      errorMessage("invalid SCIP stage\n");
      abort();
   }  /*lint !e788*/
}

/** gets number of binary active problem variables */
int SCIPgetNBinVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNBinVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nbinvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nbinvars;

   default:
      errorMessage("invalid SCIP stage\n");
      abort();
   }  /*lint !e788*/
}

/** gets number of integer active problem variables */
int SCIPgetNIntVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNIntVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nintvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nintvars;

   default:
      errorMessage("invalid SCIP stage\n");
      abort();
   }  /*lint !e788*/
}

/** gets number of implicit integer active problem variables */
int SCIPgetNImplVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNImplVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nimplvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nimplvars;

   default:
      errorMessage("invalid SCIP stage\n");
      abort();
   }  /*lint !e788*/
}

/** gets number of continuous active problem variables */
int SCIPgetNContVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNContVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->ncontvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->ncontvars;

   default:
      errorMessage("invalid SCIP stage\n");
      abort();
   }  /*lint !e788*/
}

/** gets variables of the original problem along with the numbers of different variable types; data may become invalid
 *  after a call to SCIPchgVarType()
 */
RETCODE SCIPgetOrigVarsData(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*             nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*             nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*             nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*             nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*             ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetOrigVarsData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( vars != NULL )
      *vars = scip->origprob->vars;
   if( nvars != NULL )
      *nvars = scip->origprob->nvars;
   if( nbinvars != NULL )
      *nbinvars = scip->origprob->nbinvars;
   if( nintvars != NULL )
      *nintvars = scip->origprob->nintvars;
   if( nimplvars != NULL )
      *nimplvars = scip->origprob->nimplvars;
   if( ncontvars != NULL )
      *ncontvars = scip->origprob->ncontvars;

   return SCIP_OKAY;
}

/** gets array with original problem variables; data may become invalid after
 *  a call to SCIPchgVarType()
 */
VAR** SCIPgetOrigVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetOrigVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->vars;
}

/** gets number of original problem variables */
int SCIPgetNOrigVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNOrigVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nvars;
}

/** gets number of binary original problem variables */
int SCIPgetNOrigBinVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNOrigBinVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nbinvars;
}

/** gets number of integer original problem variables */
int SCIPgetNOrigIntVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNOrigIntVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nintvars;
}

/** gets number of implicit integer original problem variables */
int SCIPgetNOrigImplVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNOrigImplVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nimplvars;
}

/** gets number of continuous original problem variables */
int SCIPgetNOrigContVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNOrigContVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->ncontvars;
}

/** returns variable of given name in the problem, or NULL if not existing */
VAR* SCIPfindVar(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of variable to find */
   )
{
   VAR* var;

   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobFindVar(scip->origprob, name);

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      var = SCIPprobFindVar(scip->transprob, name);
      if( var == NULL )
         return SCIPprobFindVar(scip->origprob, name);
      else
         return var;

   default:
      errorMessage("invalid SCIP stage\n");
      abort();
   }  /*lint !e788*/
}

/** returns TRUE iff all potential variables exist in the problem, and FALSE, if there may be additional variables,
 *  that will be added in pricing
 */
Bool SCIPallVarsInProb(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPallVarsInProb", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return (scip->set->nactivepricers == 0);
}

/** adds constraint to the problem; if constraint is only valid locally, it is added to the local subproblem of the
 *  current node (and all of its subnodes); otherwise it is added to the global problem;
 *  if a local constraint is added at the root node, it is automatically upgraded into a global constraint
 */
RETCODE SCIPaddCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddCons", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPprobAddCons(scip->origprob, scip->set, scip->stat, cons) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPprobAddCons(scip->transprob, scip->set, scip->stat, cons) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      if( SCIPconsIsGlobal(cons) )
      {
         CHECK_OKAY( SCIPprobAddCons(scip->transprob, scip->set, scip->stat, cons) );
      }
      else
      {
         CHECK_OKAY( SCIPnodeAddCons(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
               scip->tree, cons) );
      }
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was added, or from the problem, if it was a problem constraint
 */
RETCODE SCIPdelCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to delete */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPdelCons", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(cons->addconssetchg == NULL);
      CHECK_OKAY( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->origprob) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      assert(cons->addconssetchg == NULL);
      /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPconsDelete(cons, scip->mem->solvemem, scip->set, scip->stat, scip->transprob) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** returns constraint of given name in the problem, or NULL if not existing */
CONS* SCIPfindCons(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of constraint to find */
   )
{
   CONS* cons;

   assert(name != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfindCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobFindCons(scip->origprob, name);

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      cons = SCIPprobFindCons(scip->transprob, name);
      if( cons == NULL )
         return SCIPprobFindCons(scip->origprob, name);
      else
         return cons;

   default:
      errorMessage("invalid SCIP stage\n");
      abort();
   }  /*lint !e788*/
}




/*
 * local subproblem methods
 */

/** adds constraint to the given node (and all of its subnodes), even if it is a global constraint;
 *  if a local constraint is added to the root node, it is automatically upgraded into a global constraint
 */
RETCODE SCIPaddConsNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to add constraint to */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);
   assert(node != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddConsNode", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPnodeGetDepth(node) == 0 )
   {
      assert(node == scip->tree->root);
      cons->local = FALSE;
      CHECK_OKAY( SCIPprobAddCons(scip->transprob, scip->set, scip->stat, cons) );
   }
   else
   {
      CHECK_OKAY( SCIPnodeAddCons(node, scip->mem->solvemem, scip->set, scip->stat, scip->tree, cons) );
   }

   return SCIP_OKAY;
}

/** adds constraint locally to the current node (and all of its subnodes), even if it is a global constraint;
 *  if a local constraint is added at the root node, it is automatically upgraded into a global constraint
 */
RETCODE SCIPaddConsLocal(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddConsLocal", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPaddConsNode(scip, SCIPtreeGetCurrentNode(scip->tree), cons) );
   
   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities at the given node (and all subnodes) */
RETCODE SCIPdisableConsNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to disable constraint in */
   CONS*            cons                /**< constraint to disable */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPdisableConsNode", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPnodeDisableCons(node, scip->mem->solvemem, scip->set, scip->stat, scip->tree, cons) );
   
   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities at the current node (and all subnodes);
 *  if the method is called during problem modification or presolving, the constraint is globally deleted from the problem
 */
RETCODE SCIPdisableConsLocal(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to disable */
   )
{
   NODE* node;

   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPdisableConsLocal", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(cons->addconssetchg == NULL);
      CHECK_OKAY( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->origprob) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      assert(cons->addconssetchg == NULL);
      CHECK_OKAY( SCIPconsDelete(cons, scip->mem->solvemem, scip->set, scip->stat, scip->transprob) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      node = SCIPtreeGetCurrentNode(scip->tree);
      CHECK_OKAY( SCIPnodeDisableCons(node, scip->mem->solvemem, scip->set, scip->stat, scip->tree, cons) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets dual bound of current node */
Real SCIPgetLocalDualbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   NODE* node;

   CHECK_ABORT( checkStage(scip, "SCIPgetLocalDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   node = SCIPtreeGetCurrentNode(scip->tree);
   return node != NULL
      ? SCIPprobExternObjval(scip->origprob, scip->set, 
         SCIPprobExternObjval(scip->transprob, scip->set, SCIPnodeGetLowerbound(node)))
      : SCIP_INVALID;
}

/** gets lower bound of current node in transformed problem */
Real SCIPgetLocalLowerbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   NODE* node;

   CHECK_ABORT( checkStage(scip, "SCIPgetLocalLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   node = SCIPtreeGetCurrentNode(scip->tree);

   return node != NULL ? SCIPnodeGetLowerbound(node) : SCIP_INVALID;
}

/** gets dual bound of given node */
Real SCIPgetNodeDualbound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node                /**< node to get dual bound for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNodeDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->origprob, scip->set, 
      SCIPprobExternObjval(scip->transprob, scip->set, SCIPnodeGetLowerbound(node)));
}

/** gets lower bound of given node in transformed problem */
Real SCIPgetNodeLowerbound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node                /**< node to get dual bound for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNodeLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPnodeGetLowerbound(node);
}

/** if given value is tighter (larger for minimization, smaller for maximization) than the current node's dual bound,
 *  sets the current node's dual bound to the new value
 */
RETCODE SCIPupdateLocalDualbound(
   SCIP*            scip,               /**< SCIP data structure */
   Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPupdateLocalDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(SCIPtreeGetCurrentNode(scip->tree), scip->stat,
      SCIPprobInternObjval(scip->transprob, scip->set, SCIPprobInternObjval(scip->origprob, scip->set, newbound)));

   return SCIP_OKAY;
}

/** if given value is larger than the current node's lower bound (in transformed problem), sets the current node's
 *  lower bound to the new value
 */
RETCODE SCIPupdateLocalLowerbound(
   SCIP*            scip,               /**< SCIP data structure */
   Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPupdateLocalLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(SCIPtreeGetCurrentNode(scip->tree), scip->stat, newbound);

   return SCIP_OKAY;
}

/** if given value is tighter (larger for minimization, smaller for maximization) than the node's dual bound,
 *  sets the node's dual bound to the new value
 */
RETCODE SCIPupdateNodeDualbound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to update dual bound for */
   Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPupdateNodeDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(node, scip->stat, SCIPprobInternObjval(scip->transprob, scip->set, 
         SCIPprobInternObjval(scip->origprob, scip->set, newbound)));

   return SCIP_OKAY;
}

/** if given value is larger than the node's lower bound (in transformed problem), sets the node's lower bound
 *  to the new value
 */
RETCODE SCIPupdateNodeLowerbound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to update lower bound for */
   Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPupdateNodeLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(node, scip->stat, newbound);

   return SCIP_OKAY;
}




/*
 * solve methods
 */

/** initializes solving data structures and transforms problem */
static
RETCODE transformProb(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   int h;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->stat != NULL);
   assert(scip->stat->status == SCIP_STATUS_UNKNOWN);
   assert(scip->stage == SCIP_STAGE_PROBLEM);

   /* check, if a node selector exists */
   if( SCIPsetGetNodesel(scip->set, scip->stat) == NULL )
   {
      errorMessage("no node selector available\n");
      return SCIP_PLUGINNOTFOUND;
   }
   
   /* remember number of constraints */
   SCIPprobMarkNConss(scip->origprob);
   
   /* switch stage to TRANSFORMING */
   scip->stage = SCIP_STAGE_TRANSFORMING;

   /* mark statistics before solving */
   SCIPstatMark(scip->stat);

   /* init solve data structures */
   CHECK_OKAY( SCIPeventfilterCreate(&scip->eventfilter, scip->mem->solvemem) );
   CHECK_OKAY( SCIPeventqueueCreate(&scip->eventqueue) );
   CHECK_OKAY( SCIPbranchcandCreate(&scip->branchcand) );
   CHECK_OKAY( SCIPlpCreate(&scip->lp, scip->set, SCIPprobGetName(scip->origprob)) );
   CHECK_OKAY( SCIPprimalCreate(&scip->primal) );

   /* copy problem in solve memory */
   CHECK_OKAY( SCIPprobTransform(scip->origprob, scip->mem->solvemem, scip->set, scip->stat, scip->lp, 
         scip->branchcand, scip->eventfilter, scip->eventqueue, &scip->transprob) );

   /* switch stage to TRANSFORMED */
   scip->stage = SCIP_STAGE_TRANSFORMED;

   /* update upper bound (e.g. objective limit) in primal data */
   CHECK_OKAY( SCIPprimalUpdateUpperbound(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
         scip->tree, scip->lp) );

   /* print transformed problem statistics */
   infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
      "transformed problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
      scip->transprob->nvars, scip->transprob->nbinvars, scip->transprob->nintvars, scip->transprob->nimplvars,
      scip->transprob->ncontvars, scip->transprob->nconss);
   
   for( h = 0; h < scip->set->nconshdlrs; ++h )
   {
      int nconss;
      
      nconss = SCIPconshdlrGetNConss(scip->set->conshdlrs[h]);
      if( nconss > 0 )
      {
         infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            " %5d constraints of type <%s>\n", nconss, SCIPconshdlrGetName(scip->set->conshdlrs[h]));
      }
   }
   infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL, "\n");

   /* init callback methods */
   CHECK_OKAY( SCIPsetInitCallbacks(scip->set) );

   return SCIP_OKAY;
}

/** loops through the included presolvers and constraint's presolve methods, until changes are too few */
static
RETCODE presolve(
   SCIP*            scip,               /**< SCIP data structure */
   Bool*            unbounded,          /**< pointer to store whether presolving detected unboundness */
   Bool*            infeasible          /**< pointer to store whether presolving detected infeasibility */
   )
{
   RESULT result;
   Bool aborted;
   Real abortfac;
   int maxnrounds;
   int nrounds;
   int nfixedvars;
   int naggrvars;
   int nchgvartypes;
   int nchgbds;
   int naddholes;
   int ndelconss;
   int nupgdconss;
   int nchgcoefs;
   int nchgsides;
   int lastnfixedvars;
   int lastnaggrvars;
   int lastnchgvartypes;
   int lastnchgbds;
   int lastnaddholes;
   int lastndelconss;
   int lastnupgdconss;
   int lastnchgcoefs;
   int lastnchgsides;
   int i;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->stat->status == SCIP_STATUS_UNKNOWN);
   assert(scip->transprob != NULL);
   assert(scip->stage == SCIP_STAGE_TRANSFORMED);
   assert(unbounded != NULL);
   assert(infeasible != NULL);

   *unbounded = FALSE;
   *infeasible = FALSE;

   /* update upper bound (e.g. objective limit) in primal data */
   CHECK_OKAY( SCIPprimalUpdateUpperbound(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
         scip->tree, scip->lp) );

   /* switch stage to PRESOLVING */
   scip->stage = SCIP_STAGE_PRESOLVING;
   
   /* start presolving timer */
   SCIPclockStart(scip->stat->presolvingtime, scip->set);

   nrounds = 0;
   nfixedvars = 0;
   naggrvars = 0;
   nchgvartypes = 0;
   nchgbds = 0;
   naddholes = 0;
   ndelconss = 0;
   nupgdconss = 0;
   nchgcoefs = 0;
   nchgsides = 0;

   maxnrounds = scip->set->presol_maxrounds;
   if( maxnrounds == -1 )
      maxnrounds = INT_MAX;

   abortfac = scip->set->presol_abortfac;

   aborted = SCIPsolveIsStopped(scip->set, scip->stat);
   result = SCIP_DIDNOTRUN;

   infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "presolving:\n");

   /* inform presolvers that the presolving is abound to begin */
   for( i = 0; i < scip->set->npresols; ++i )
   {
      CHECK_OKAY( SCIPpresolInitpre(scip->set->presols[i], scip, &result) );
      *unbounded = *unbounded || (result == SCIP_UNBOUNDED);
      *infeasible = *infeasible || (result == SCIP_CUTOFF);
   }

   /* inform constraint handlers that the presolving is abound to begin */
   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrInitpre(scip->set->conshdlrs[i], scip, &result) );
      *unbounded = *unbounded || (result == SCIP_UNBOUNDED);
      *infeasible =*infeasible ||  (result == SCIP_CUTOFF);
   }
   aborted = aborted || *unbounded || *infeasible;

   /* perform presolving rounds */
   while( nrounds < maxnrounds && !aborted )
   {
      lastnfixedvars = nfixedvars;
      lastnaggrvars = naggrvars;
      lastnchgvartypes = nchgvartypes;
      lastnchgbds = nchgbds;
      lastnaddholes = naddholes;
      lastndelconss = ndelconss;
      lastnupgdconss = nupgdconss;
      lastnchgcoefs = nchgcoefs;
      lastnchgsides = nchgsides;

      /* sort presolvers by priority */
      SCIPsetSortPresols(scip->set);

      /* call included presolvers */
      for( i = 0; i < scip->set->npresols && !(*unbounded) && !(*infeasible); ++i )
      {
         CHECK_OKAY( SCIPpresolExec(scip->set->presols[i], scip->set, nrounds, 
               &nfixedvars, &naggrvars, &nchgvartypes, &nchgbds, &naddholes,
               &ndelconss, &nupgdconss, &nchgcoefs, &nchgsides, &result) );
         *unbounded = *unbounded || (result == SCIP_UNBOUNDED);
         *infeasible = *infeasible || (result == SCIP_CUTOFF);
      }

      /* call presolve methods of constraint handlers */
      for( i = 0; i < scip->set->nconshdlrs && !(*unbounded) && !(*infeasible); ++i )
      {
         CHECK_OKAY( SCIPconshdlrPresolve(scip->set->conshdlrs[i], scip->mem->solvemem, scip->set, scip->stat, 
               scip->transprob, nrounds,
               &nfixedvars, &naggrvars, &nchgvartypes, &nchgbds, &naddholes,
               &ndelconss, &nupgdconss, &nchgcoefs, &nchgsides, &result) );
         *unbounded = *unbounded || (result == SCIP_UNBOUNDED);
         *infeasible = *infeasible || (result == SCIP_CUTOFF);
      }

      /* check, if we should abort presolving due to not enough changes in the last round */
      aborted = TRUE;

      /* don't abort, if enough changes were applied to the variables */
      aborted = aborted
         && (scip->transprob->nvars == 0
            || (nfixedvars - lastnfixedvars + naggrvars - lastnaggrvars + nchgvartypes - lastnchgvartypes
               + nchgbds - lastnchgbds + naddholes - lastnaddholes < abortfac * scip->transprob->nvars));

      /* don't abort, if enough changes were applied to the constraints */
      aborted = aborted
         && (scip->transprob->nconss == 0
            || (ndelconss - lastndelconss + nupgdconss - lastnupgdconss + nchgsides - lastnchgsides
               < abortfac * scip->transprob->nconss));

      /* don't abort, if enough changes were applied to the coefficients (assume a 20% density of non-zero elements) */
      aborted = aborted
         && (scip->transprob->nvars == 0 || scip->transprob->nconss == 0
            || (nchgcoefs - lastnchgcoefs < abortfac * 0.2 * scip->transprob->nvars * scip->transprob->nconss));

      /* abort if time limit was reached or user interrupted */
      aborted = aborted || SCIPsolveIsStopped(scip->set, scip->stat);

      /* abort if problem is infeasible or unbounded */
      aborted = aborted || *unbounded || *infeasible;

      /* increase round number */
      nrounds++;

      if( !aborted )
      {
         /* print presolving statistics */
         infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "(round %d)", nrounds);
         infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
            " %d deleted vars, %d deleted constraints, %d tightened bounds, %d added holes, %d changed sides, %d changed coefficients\n",
            nfixedvars + naggrvars, ndelconss, nchgbds, naddholes, nchgsides, nchgcoefs);
      }
   }
   
   /* inform constraint handlers that the presolving is finished, and perform final modifications */
   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      CHECK_OKAY( SCIPconshdlrExitpre(scip->set->conshdlrs[i], scip, &result) );
      *unbounded = *unbounded || (result == SCIP_UNBOUNDED);
      *infeasible = *infeasible || (result == SCIP_CUTOFF);
   }

   /* inform presolvers that the presolving is finished, and perform final modifications */
   for( i = 0; i < scip->set->npresols; ++i )
   {
      CHECK_OKAY( SCIPpresolExitpre(scip->set->presols[i], scip, &result) );
      *unbounded = *unbounded || (result == SCIP_UNBOUNDED);
      *infeasible = *infeasible || (result == SCIP_CUTOFF);
   }

   /* print presolving statistics */
   infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL, "presolving (%d rounds):\n", nrounds);
   infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      " %d deleted vars, %d deleted constraints, %d tightened bounds, %d added holes, %d changed sides, %d changed coefficients\n",
      nfixedvars + naggrvars, ndelconss, nchgbds, naddholes, nchgsides, nchgcoefs);

   /* replace variables in variable bounds with active problem variables, and 
    * check, whether the objective value is always integral
    */
   CHECK_OKAY( SCIPprobExitPresolve(scip->transprob, scip->set) );

   /* stop presolving time */
   SCIPclockStop(scip->stat->presolvingtime, scip->set);
   
   /* switch stage to PRESOLVED */
   scip->stage = SCIP_STAGE_PRESOLVED;

   return SCIP_OKAY;
}

/** initializes solution process data structures */
static
RETCODE initSolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   int h;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->stage == SCIP_STAGE_PRESOLVED);

   /* update upper bound (e.g. objective limit) in primal data */
   CHECK_OKAY( SCIPprimalUpdateUpperbound(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
         scip->tree, scip->lp) );

   /* switch stage to INITSOLVE */
   scip->stage = SCIP_STAGE_INITSOLVE;

   /* create VBC output file */
   CHECK_OKAY( SCIPvbcInit(scip->stat->vbc, scip->mem->solvemem, scip->set) );

   /* init solution process data structures */
   CHECK_OKAY( SCIPpricestoreCreate(&scip->pricestore) );
   CHECK_OKAY( SCIPsepastoreCreate(&scip->sepastore) );
   CHECK_OKAY( SCIPcutpoolCreate(&scip->cutpool, scip->mem->solvemem, scip->set->sepa_cutagelimit) );
   CHECK_OKAY( SCIPconflictCreate(&scip->conflict, scip->set) );
   CHECK_OKAY( SCIPtreeCreate(&scip->tree, scip->mem->solvemem, scip->set, scip->stat, scip->lp, 
         SCIPsetGetNodesel(scip->set, scip->stat)) );

   /* increase number of branch and bound runs */
   scip->stat->nruns++;

   /* switch stage to SOLVING */
   scip->stage = SCIP_STAGE_SOLVING;

   /* inform the transformed problem that the branch and bound process starts now */
   CHECK_OKAY( SCIPprobInitSolve(scip->transprob, scip->set) );

   /* inform constraint handlers that the branch and bound process starts now */
   for( h = 0; h < scip->set->nconshdlrs; ++h )
   {
      CHECK_OKAY( SCIPconshdlrInitsol(scip->set->conshdlrs[h], scip) );
   }
   
   /* inform primal heuristics that the branch and bound process starts now */
   for( h = 0; h < scip->set->nheurs; ++h )
   {
      CHECK_OKAY( SCIPheurInitsol(scip->set->heurs[h], scip->set) );
   }

   /* remember number of constraints */
   SCIPprobMarkNConss(scip->transprob);

   /* if all variables are known, calculate a trivial primal bound by setting all variables to their worst bound */
   if( scip->set->nactivepricers == 0 )
   {
      VAR* var;
      Real obj;
      Real objbound;
      Real bd;
      int v;
      
      objbound = 0.0;
      for( v = 0; v < scip->transprob->nvars && !SCIPsetIsInfinity(scip->set, objbound); ++v )
      {
         var = scip->transprob->vars[v];
         obj = SCIPvarGetObj(var);
         if( !SCIPsetIsZero(scip->set, obj) )
         {
            bd = SCIPvarGetWorstBound(var);
            if( SCIPsetIsInfinity(scip->set, REALABS(bd)) )
               objbound = SCIPsetInfinity(scip->set);
            else
               objbound += obj * bd;
         }
      }

      /* update primal bound (add 1.0 to primal bound, such that solution with worst bound may be found) */
      if( !SCIPsetIsInfinity(scip->set, objbound) && SCIPsetIsLT(scip->set, objbound + 1.0, scip->primal->upperbound) )
      {
         CHECK_OKAY( SCIPprimalSetUpperbound(scip->primal, scip->mem->solvemem, scip->set, scip->stat,
               scip->transprob, scip->tree, scip->lp, objbound + 1.0) );
      }
   }

   return SCIP_OKAY;
}

/** frees solution process data structures */
static
RETCODE freeSolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   int h;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->stage == SCIP_STAGE_SOLVING || scip->stage == SCIP_STAGE_SOLVED);

   /* remove focus from the current focus node */
   if( SCIPtreeGetFocusNode(scip->tree) != NULL )
   {
      NODE* node = NULL;
      Bool cutoff;

      CHECK_OKAY( SCIPnodeFocus(&node, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->primal,
            scip->tree, scip->lp, scip->branchcand, scip->eventfilter, scip->eventqueue, &cutoff) );
      assert(!cutoff);
   }

   /* switch stage to FREESOLVE */
   scip->stage = SCIP_STAGE_FREESOLVE;

   /* switch status to UNKNOWN */
   scip->stat->status = SCIP_STATUS_UNKNOWN;

   /* clear the LP, and flush the changes to clear the LP of the solver */
   CHECK_OKAY( SCIPlpClear(scip->lp, scip->mem->solvemem, scip->set) );
   CHECK_OKAY( SCIPlpFlush(scip->lp, scip->mem->solvemem, scip->set) );
   
   /* clear all row references in internal data structures */
   CHECK_OKAY( SCIPcutpoolClear(scip->cutpool, scip->mem->solvemem, scip->set, scip->lp) );

   /* deinitialize constraint handlers */
   for( h = 0; h < scip->set->nconshdlrs; ++h )
   {
      CHECK_OKAY( SCIPconshdlrExitsol(scip->set->conshdlrs[h], scip) );
   }

   /* we have to free the tree prior to the problem deinitialization, because the rows stored in the forks and
    * subroots have to be released
    */
   CHECK_OKAY( SCIPtreeFree(&scip->tree, scip->mem->solvemem, scip->set, scip->lp) );

   /* deinitialize transformed problem */
   CHECK_OKAY( SCIPprobExitSolve(scip->transprob, scip->mem->solvemem, scip->set, scip->lp) );

   /* free solution process data structures */
   CHECK_OKAY( SCIPconflictFree(&scip->conflict) );
   CHECK_OKAY( SCIPcutpoolFree(&scip->cutpool, scip->mem->solvemem, scip->set, scip->lp) );
   CHECK_OKAY( SCIPsepastoreFree(&scip->sepastore) );
   CHECK_OKAY( SCIPpricestoreFree(&scip->pricestore) );

   /* close VBC output file */
   SCIPvbcExit(scip->stat->vbc, scip->set);

   /* reset statistics for current branch and bound run */
   SCIPstatResetCurrentRun(scip->stat);

   /* switch stage to TRANSFORMED */
   scip->stage = SCIP_STAGE_TRANSFORMED;

   return SCIP_OKAY;
}

/** free transformed problem */
static
RETCODE freeTransform(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->stat != NULL);
   assert(scip->stat->status == SCIP_STATUS_UNKNOWN);
   assert(scip->stage == SCIP_STAGE_TRANSFORMED);

   /* exit callback methods */
   CHECK_OKAY( SCIPsetExitCallbacks(scip->set) );

   /* switch stage to FREETRANS */
   scip->stage = SCIP_STAGE_FREETRANS;

   /* free transformed problem data structures */
   CHECK_OKAY( SCIPprobFree(&scip->transprob, scip->mem->solvemem, scip->set, scip->stat, scip->lp) );
   CHECK_OKAY( SCIPprimalFree(&scip->primal, scip->mem->solvemem) );
   CHECK_OKAY( SCIPlpFree(&scip->lp, scip->mem->solvemem, scip->set) );
   CHECK_OKAY( SCIPbranchcandFree(&scip->branchcand) );
   CHECK_OKAY( SCIPeventfilterFree(&scip->eventfilter, scip->mem->solvemem, scip->set) );
   CHECK_OKAY( SCIPeventqueueFree(&scip->eventqueue) );

   /* free the transformed block memory */
#ifndef NDEBUG
   blockMemoryCheckEmpty(scip->mem->solvemem);
#endif
   freeAllBlockMemoryNull(scip->mem->solvemem);

   /* reset statistics to the point before solving started */
   SCIPstatReset(scip->stat);

   /* switch stage to PROBLEM */
   scip->stage = SCIP_STAGE_PROBLEM;

   return SCIP_OKAY;
}

/** presolves problem */
RETCODE SCIPpresolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   Bool unbounded;
   Bool infeasible;

   CHECK_OKAY( checkStage(scip, "SCIPpresolve", FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* capture the CTRL-C interrupt */
   if( scip->set->misc_catchctrlc )
      SCIPinterruptCapture(scip->interrupt);

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      /* initialize solving data structures and transform problem */
      CHECK_OKAY( transformProb(scip) );
      assert(scip->stage == SCIP_STAGE_TRANSFORMED);

      /*lint -fallthrough */

   case SCIP_STAGE_TRANSFORMED:
      /* presolve problem */
      CHECK_OKAY( presolve(scip, &unbounded, &infeasible) );
      assert(scip->stage == SCIP_STAGE_PRESOLVED);

      if( infeasible || unbounded )
      {
         /* initialize solving process data structures to be able to switch to SOLVED stage */
         CHECK_OKAY( initSolve(scip) );

         /* switch stage to SOLVED */
         scip->stage = SCIP_STAGE_SOLVED;

         /* print solution message */
         if( infeasible )
         {
            infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL, "presolving detected infeasibility.\n");

            /* switch status to INFEASIBLE */
            scip->stat->status = SCIP_STATUS_INFEASIBLE;
         }
         else if( scip->primal->nsols >= 1 )
         {
            infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL, "presolving detected unboundness.\n");

            /* switch status to UNBOUNDED */
            scip->stat->status = SCIP_STATUS_UNBOUNDED;
         }
         else
         {
            infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL, "presolving detected unboundness (or infeasibility).\n");

            /* switch status to INFORUNBD */
            scip->stat->status = SCIP_STATUS_INFORUNBD;
         }
      }
      else
      {
         int h;

         /* print presolved problem statistics */
         infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
            "presolved problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
            scip->transprob->nvars, scip->transprob->nbinvars, scip->transprob->nintvars, scip->transprob->nimplvars,
            scip->transprob->ncontvars, scip->transprob->nconss);
         
         for( h = 0; h < scip->set->nconshdlrs; ++h )
         {
            int nconss;
            
            nconss = SCIPconshdlrGetNConss(scip->set->conshdlrs[h]);
            if( nconss > 0 )
            {
               infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
                  " %5d constraints of type <%s>\n", nconss, SCIPconshdlrGetName(scip->set->conshdlrs[h]));
            }
         }

         if( SCIPprobIsObjIntegral(scip->transprob) )
         {
            infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "objective value is always integral\n");
         }
      }
      
      /* display timing statistics */
      infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
         "Presolving Time: %.2f\n", SCIPclockGetTime(scip->stat->presolvingtime));
      break;

   case SCIP_STAGE_PRESOLVED:
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }

   /* release the CTRL-C interrupt */
   if( scip->set->misc_catchctrlc )
      SCIPinterruptRelease(scip->interrupt);

   return SCIP_OKAY;
}

/** presolves and solves problem */
RETCODE SCIPsolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   Bool restart;

   CHECK_OKAY( checkStage(scip, "SCIPsolve", FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   /* check, if a node selector exists */
   if( SCIPsetGetNodesel(scip->set, scip->stat) == NULL )
   {
      errorMessage("no node selector available\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* automatic restarting loop */
   do
   {
      restart = FALSE;

      /* start solving timer */
      SCIPclockStart(scip->stat->solvingtime, scip->set);

      switch( scip->stage )
      {
      case SCIP_STAGE_PROBLEM:
      case SCIP_STAGE_TRANSFORMED:
         /* initialize solving data structures, transform and problem */
         CHECK_OKAY( SCIPpresolve(scip) );
         if( scip->stage == SCIP_STAGE_SOLVED )
            break;
         assert(scip->stage == SCIP_STAGE_PRESOLVED);

         /*lint -fallthrough*/

      case SCIP_STAGE_PRESOLVED:
         /* initialize solving process data structures */
         CHECK_OKAY( initSolve(scip) );
         assert(scip->stage == SCIP_STAGE_SOLVING);
         infoMessage(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL, "\n");
      
         /*lint -fallthrough*/

      case SCIP_STAGE_SOLVING:
         /* reset display */
         SCIPstatResetDisplay(scip->stat);

         /* capture the CTRL-C interrupt */
         if( scip->set->misc_catchctrlc )
            SCIPinterruptCapture(scip->interrupt);

         /* continue solution process */
         CHECK_OKAY( SCIPsolveCIP(scip->mem->solvemem, scip->set, scip->stat, scip->mem, scip->transprob,
               scip->primal, scip->tree, scip->lp, scip->pricestore, scip->sepastore, scip->cutpool,
               scip->branchcand, scip->conflict, scip->eventfilter, scip->eventqueue, &restart) );

         /* release the CTRL-C interrupt */
         if( scip->set->misc_catchctrlc )
            SCIPinterruptRelease(scip->interrupt);
         
         /* detect, whether problem is solved */
         if( SCIPtreeGetNNodes(scip->tree) == 0 && SCIPtreeGetCurrentNode(scip->tree) == NULL )
         {
            assert(scip->stat->status == SCIP_STATUS_OPTIMAL
               || scip->stat->status == SCIP_STATUS_INFEASIBLE
               || scip->stat->status == SCIP_STATUS_UNBOUNDED
               || scip->stat->status == SCIP_STATUS_INFORUNBD);

            /* tree is empty, and no current node exists -> problem is solved */
            scip->stage = SCIP_STAGE_SOLVED;
         }
         else if( restart )
         {
            /* free the solving process data in order to restart */
            SCIPmessage(scip, SCIP_VERBLEVEL_NORMAL,
               "(run %d) restarting after %d root node bound changes\n\n",
               scip->stat->nruns, scip->stat->nrootboundchgsrun);
            CHECK_OKAY( SCIPfreeSolve(scip) );
            assert(scip->stage == SCIP_STAGE_TRANSFORMED);
         }
         break;

      case SCIP_STAGE_SOLVED:
         assert(scip->stat->status == SCIP_STATUS_OPTIMAL
            || scip->stat->status == SCIP_STATUS_INFEASIBLE
            || scip->stat->status == SCIP_STATUS_UNBOUNDED
            || scip->stat->status == SCIP_STATUS_INFORUNBD);
         break;

      default:
         errorMessage("invalid SCIP stage\n");
         return SCIP_ERROR;
      }  /*lint !e788*/

      /* stop solving timer */
      SCIPclockStop(scip->stat->solvingtime, scip->set);
   }
   while( restart && !SCIPsolveIsStopped(scip->set, scip->stat) );

   /* display most relevant statistics */
   if( scip->set->disp_verblevel >= SCIP_VERBLEVEL_NORMAL )
   {
      printf("\n");
      printf("SCIP Status        : ");
      CHECK_OKAY( SCIPprintStage(scip, NULL) );
      printf("\n");
      printf("Solving Time (sec) : %.2f\n", SCIPclockGetTime(scip->stat->solvingtime));
      if( scip->stat->nruns > 1 )
         printf("Solving Nodes      : %lld (total of %lld nodes in %d runs)\n", 
            scip->stat->nnodes, scip->stat->ntotalnodes, scip->stat->nruns);
      else
         printf("Solving Nodes      : %lld\n", scip->stat->nnodes);
      if( scip->stage >= SCIP_STAGE_TRANSFORMED && scip->stage <= SCIP_STAGE_FREESOLVE )
         printf("Primal Bound       : %+.14e (%lld solutions)\n", SCIPgetPrimalbound(scip), scip->primal->nsolsfound);
      if( scip->stage >= SCIP_STAGE_SOLVING && scip->stage <= SCIP_STAGE_SOLVED )
      {
         printf("Dual Bound         : %+.14e\n", SCIPgetDualbound(scip));
         printf("Gap                : ");
         if( SCIPsetIsInfinity(scip->set, SCIPgetGap(scip)) )
            printf("infinite\n");
         else
            printf("%.2f %%\n", 100.0*SCIPgetGap(scip));
      }
   }

   return SCIP_OKAY;
}

/** frees branch and bound tree and all solution process data; statistics, presolving data and transformed problem is
 *  preserved
 */
RETCODE SCIPfreeSolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeSolve", TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMED:
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
      /* switch stage to TRANSFORMED */
      scip->stage = SCIP_STAGE_TRANSFORMED;
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      /* free solution process data structures */
      CHECK_OKAY( freeSolve(scip) );
      assert(scip->stage == SCIP_STAGE_TRANSFORMED);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** frees all solution process data including presolving and transformed problem, only original problem is kept */
RETCODE SCIPfreeTransform(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeTransform", TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      /* free solution process data */
      CHECK_OKAY( SCIPfreeSolve(scip) );
      assert(scip->stage == SCIP_STAGE_TRANSFORMED);

      /*lint -fallthrough */

   case SCIP_STAGE_TRANSFORMED:
      /* free transformed problem data structures */
      CHECK_OKAY( freeTransform(scip) );
      assert(scip->stage == SCIP_STAGE_PROBLEM);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}




/*
 * variable methods
 */

/** creates and captures problem variable; if variable is of integral type, fractional bounds are automatically rounded; 
 *  an integer variable with bounds zero and one is automatically converted into a binary variable
 */
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
   )
{
   assert(var != NULL);
   assert(lb <= ub);

   CHECK_OKAY( checkStage(scip, "SCIPcreateVar", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarCreateOriginal(var, scip->mem->probmem, scip->set, scip->stat, 
            name, lb, ub, obj, vartype, initial, removeable, vardelorig, vartrans, vardeltrans, vardata) );
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPvarCreateTransformed(var, scip->mem->solvemem, scip->set, scip->stat,
            name, lb, ub, obj, vartype, initial, removeable, NULL, NULL, vardeltrans, vardata) );
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** increases usage counter of variable */
RETCODE SCIPcaptureVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to capture */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIPvarCapture(var);

   return SCIP_OKAY;
}

/** decreases usage counter of variable, and frees memory if necessary */
RETCODE SCIPreleaseVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var                 /**< pointer to variable */
   )
{
   assert(var != NULL);
   assert(*var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreleaseVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarRelease(var, scip->mem->probmem, scip->set, scip->lp) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      if( !SCIPvarIsTransformed(*var) )
      {
         errorMessage("cannot release original variables while the transformed problem exists\n");
         return SCIP_INVALIDCALL;
      }
      CHECK_OKAY( SCIPvarRelease(var, scip->mem->solvemem, scip->set, scip->lp) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 */
RETCODE SCIPtransformVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get/create transformed variable for */
   VAR**            transvar            /**< pointer to store the transformed variable */
   )
{
   assert(transvar != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPtransformVar", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPvarIsTransformed(var) )
   {
      *transvar = var;
      SCIPvarCapture(*transvar);
   }
   else
   {
      CHECK_OKAY( SCIPvarTransform(var, scip->mem->solvemem, scip->set, scip->stat, scip->origprob->objsense, transvar) );
   }

   return SCIP_OKAY;
}

/** gets and captures transformed variables for an array of variables;
 *  if a variable of the array is not yet transformed, a new transformed variable for this variable is created;
 *  it is possible to call this method with vars == transvars
 */
RETCODE SCIPtransformVars(
   SCIP*            scip,               /**< SCIP data structure */
   int              nvars,              /**< number of variables to get/create transformed variables for */
   VAR**            vars,               /**< array with variables to get/create transformed variables for */
   VAR**            transvars           /**< array to store the transformed variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || transvars != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPtransformVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarIsTransformed(vars[v]) )
      {
         transvars[v] = vars[v];
         SCIPvarCapture(transvars[v]);
      }
      else
      {
         CHECK_OKAY( SCIPvarTransform(vars[v], scip->mem->solvemem, scip->set, scip->stat, scip->origprob->objsense,
               &transvars[v]) );
      }
   }

   return SCIP_OKAY;
}

/** gets corresponding transformed variable of a given variable;
 *  returns NULL as transvar, if transformed variable is not yet existing
 */
RETCODE SCIPgetTransformedVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get transformed variable for */
   VAR**            transvar            /**< pointer to store the transformed variable */
   )
{
   assert(transvar != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetTransformedVar", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPvarIsTransformed(var) )
      *transvar = var;
   else
   {
      CHECK_OKAY( SCIPvarGetTransformed(var, scip->mem->solvemem, scip->set, scip->stat, transvar) );
   }

   return SCIP_OKAY;
}

/** gets corresponding transformed variables for an array of variables;
 *  stores NULL in a transvars slot, if the transformed variable is not yet existing;
 *  it is possible to call this method with vars == transvars, but remember that variables that are not
 *  yet transformed will be replaced with NULL
 */
RETCODE SCIPgetTransformedVars(
   SCIP*            scip,               /**< SCIP data structure */
   int              nvars,              /**< number of variables to get transformed variables for */
   VAR**            vars,               /**< array with variables to get transformed variables for */
   VAR**            transvars           /**< array to store the transformed variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || transvars != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetTransformedVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarIsTransformed(vars[v]) )
         transvars[v] = vars[v];
      else
      {
         CHECK_OKAY( SCIPvarGetTransformed(vars[v], scip->mem->solvemem, scip->set, scip->stat, &transvars[v]) );
      }
   }

   return SCIP_OKAY;
}

/** gets negated variable x' = lb + ub - x of variable x; negated variable is created, if not yet existing */
RETCODE SCIPgetNegatedVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get negated variable for */
   VAR**            negvar              /**< pointer to store the negated variable */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetNegatedVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( !SCIPvarIsTransformed(var) )
   {
      CHECK_OKAY( SCIPvarNegate(var, scip->mem->probmem, scip->set, scip->stat, negvar) );
   }
   else
   {
      assert(scip->stage != SCIP_STAGE_PROBLEM);
      CHECK_OKAY( SCIPvarNegate(var, scip->mem->solvemem, scip->set, scip->stat, negvar) );
   }

   return SCIP_OKAY;
}

/** gets a binary variable that is equal to the given binary variable, and that is either active or the negated
 *  variable of an active binary variable; if the given variable is fixed, NULL is returned as representative,
 *  and *negated is TRUE iff the variable is fixed to TRUE
 */
RETCODE SCIPgetBinvarRepresentative(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< binary variable to get binary representative for */
   VAR**            repvar,             /**< pointer to store the binary representative */
   Bool*            negated             /**< pointer to store whether the negation of an active variable was returned */
   )
{
   assert(repvar != NULL);
   assert(negated != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetBinvarRepresentative", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   /* get the active representative of the given variable */
   *repvar = var;
   *negated = FALSE;
   CHECK_OKAY( SCIPvarGetProbvarBinary(repvar, negated) );

   /* negate the representative, if it corresponds to the negation of the given variable */
   if( *repvar != NULL && *negated )
   {
      CHECK_OKAY( SCIPvarNegate(*repvar, scip->mem->solvemem, scip->set, scip->stat, repvar) );
   }

   return SCIP_OKAY;
}

/** gets solution value for variable in current node */
Real SCIPgetVarSol(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get solution value for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPvarGetSol(var, SCIPtreeHasCurrentNodeLP(scip->tree));
}

/** gets solution values of multiple variables in current node */
RETCODE SCIPgetVarSols(
   SCIP*            scip,               /**< SCIP data structure */
   int              nvars,              /**< number of variables to get solution value for */
   VAR**            vars,               /**< array with variables to get value for */
   Real*            vals                /**< array to store solution values of variables */
   )
{
   int v;
   Bool uselp;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetVarSols", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      for( v = 0; v < nvars; ++v )
         vals[v] = SCIPvarGetLPSol(vars[v]);
   }
   else
   {
      for( v = 0; v < nvars; ++v )
         vals[v] = SCIPvarGetPseudoSol(vars[v]);
   }

   return SCIP_OKAY;
}

/** gets strong branching information on COLUMN variable */
RETCODE SCIPgetVarStrongbranch(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get strong branching values for */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching column down */
   Real*            up,                 /**< stores dual bound after branching column up */
   Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   Bool*            downconflict,       /**< pointer to store whether a conflict clause was created for an infeasible
                                         *   downwards branch, or NULL */
   Bool*            upconflict,         /**< pointer to store whether a conflict clause was created for an infeasible
                                         *   upwards branch, or NULL */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   COL* col;

   assert(lperror != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetVarStrongbranch", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( downinf != NULL )
      *downinf = FALSE;
   if( upinf != NULL )
      *upinf = FALSE;
   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      errorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   if( !SCIPcolIsInLP(col) )
   {
      errorMessage("cannot get strong branching information on variable <%s> not in current LP\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   /* call strong branching for column */
   CHECK_OKAY( SCIPcolGetStrongbranch(col, scip->set, scip->stat, scip->lp, itlim, down, up, lperror) );

   /* check, if the branchings are infeasible; in exact solving mode, we cannot trust the strong branching enough
    * declare the sub nodes infeasible
    */
   if( !(*lperror) && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && !scip->set->misc_exactsolve )
   {
      Bool downcutoff;
      Bool upcutoff;

      downcutoff = SCIPsetIsGE(scip->set, col->strongbranchdown, scip->lp->cutoffbound);
      upcutoff = SCIPsetIsGE(scip->set, col->strongbranchup, scip->lp->cutoffbound);
      if( downinf != NULL )
         *downinf = downcutoff;
      if( upinf != NULL )
         *upinf = upcutoff;
   
      /* analyze infeasible strong branching sub problems:
       * because the strong branching's bound change is necessary for infeasibility, it cannot be undone;
       * therefore, infeasible strong branchings on non-binary variables will not produce a valid conflict clause
       */
      if( scip->set->conf_usesb && scip->set->nconflicthdlrs > 0
         && SCIPvarGetType(var) == SCIP_VARTYPE_BINARY
         && SCIPtreeGetCurrentDepth(scip->tree) > 0 )
      {
         if( (downcutoff && SCIPsetCeil(scip->set, col->primsol-1.0) >= col->lb - 0.5)
            || (upcutoff && SCIPsetFloor(scip->set, col->primsol+1.0) <= col->ub + 0.5) )
         {
            CHECK_OKAY( SCIPconflictAnalyzeStrongbranch(scip->conflict, scip->mem->solvemem, scip->set, scip->stat,
                  scip->transprob, scip->tree, scip->lp, col, downconflict, upconflict) );
         }
      }
   }

   return SCIP_OKAY;
}

/** gets strong branching information on COLUMN variable of the last SCIPgetVarStrongbranch() call;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given variable;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 */
RETCODE SCIPgetVarStrongbranchLast(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get last strong branching values for */
   Real*            down,               /**< stores dual bound after branching column down */
   Real*            up,                 /**< stores dual bound after branching column up */
   Real*            solval,             /**< stores LP solution value of variable at the last strong branching call */
   Real*            lpobjval            /**< stores LP objective value at last strong branching call, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetVarStrongbranchLast", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      errorMessage("cannot get strong branching information on non-COLUMN variable\n");
      return SCIP_INVALIDDATA;
   }

   SCIPcolGetStrongbranchLast(SCIPvarGetCol(var), down, up, solval, lpobjval);

   return SCIP_OKAY;
}

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given column, or -1 if strong branching was never applied to the column in current run
 */
Longint SCIPgetVarStrongbranchNode(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get last strong branching node for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetVarStrongbranchNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      errorMessage("cannot get strong branching information on non-COLUMN variable\n");
      abort();
   }

   return SCIPcolGetStrongbranchNode(SCIPvarGetCol(var));
}

/** changes variable's objective value */
RETCODE SCIPchgVarObj(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the objective value for */
   Real             newobj              /**< new objective value */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarObj", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      CHECK_OKAY( SCIPvarChgObj(var, scip->mem->probmem, scip->set, scip->primal, scip->lp, scip->eventqueue, newobj) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarChgObj(var, scip->mem->solvemem, scip->set, scip->primal, scip->lp, scip->eventqueue, newobj) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/   
}

/** adds value to variable's objective value */
RETCODE SCIPaddVarObj(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the objective value for */
   Real             addobj              /**< additional objective value */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddVarObj", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      CHECK_OKAY( SCIPvarAddObj(var, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->primal, scip->lp,
            scip->eventqueue, addobj) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarAddObj(var, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->primal, scip->lp,
            scip->eventqueue, addobj) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/   
}

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) lower bound value;
 *  does not change the bounds of the variable
 */
Real SCIPadjustedVarLb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to adjust the bound for */
   Real             lb                  /**< lower bound value to adjust */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPadjustedVarLb", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPvarAdjustLb(var, scip->set, &lb);

   return lb;
}

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) lower bound value;
 *  does not change the bounds of the variable
 */
Real SCIPadjustedVarUb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to adjust the bound for */
   Real             ub                  /**< upper bound value to adjust */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPadjustedVarUb", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPvarAdjustUb(var, scip->set, &ub);

   return ub;
}

/** depending on SCIP's stage, changes lower bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 */
RETCODE SCIPchgVarLb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarLb", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarChgLbGlobal(var, scip->set, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarChgLbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarChgLbGlobal(var, scip->set, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** depending on SCIP's stage, changes upper bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 */
RETCODE SCIPchgVarUb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarUb", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarChgUbGlobal(var, scip->set, newbound) );
      CHECK_OKAY( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarChgUbGlobal(var, scip->set, newbound) );
      CHECK_OKAY( SCIPvarChgUbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** changes lower bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 */
RETCODE SCIPchgVarLbNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for current node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarLbNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarAdjustLb(var, scip->set, &newbound);

   CHECK_OKAY( SCIPnodeAddBoundchg(node, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, 
         scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
   
   return SCIP_OKAY;
}

/** changes upper bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 */
RETCODE SCIPchgVarUbNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for current node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarUbNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarAdjustUb(var, scip->set, &newbound);

   CHECK_OKAY( SCIPnodeAddBoundchg(node, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, 
         scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
   
   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  doesn't store any inference information in the bound change, such that in conflict analysis, this change
 *  is treated like a branching decision
 */
RETCODE SCIPtightenVarLb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new value for bound */
   Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   Real lb;
   Real ub;

   assert(infeasible != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPtightenVarLb", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));
   
   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);

   if( !SCIPsetIsLbBetter(scip->set, newbound, lb) )
   {
      if( tightened != NULL )
         *tightened = FALSE;
      return SCIP_OKAY;
   }

   switch( scip->stage )
   {
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarChgLbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarChgLbGlobal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  doesn't store any inference information in the bound change, such that in conflict analysis, this change
 *  is treated like a branching decision
 */
RETCODE SCIPtightenVarUb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new value for bound */
   Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   Real lb;
   Real ub;

   assert(infeasible != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPtightenVarUb", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));
   
   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   if( !SCIPsetIsUbBetter(scip->set, newbound, ub) )
   {
      if( tightened != NULL )
         *tightened = FALSE;
      return SCIP_OKAY;
   }

   switch( scip->stage )
   {
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarChgUbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarChgUbGlobal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 */
RETCODE SCIPinferVarLbCons(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new value for bound */
   CONS*            infercons,          /**< constraint that deduced the bound change */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   Real lb;
   Real ub;

   assert(infeasible != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPinferVarLbCons", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));
   
   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);
   
   if( !SCIPsetIsLbBetter(scip->set, newbound, lb) )
   {
      if( tightened != NULL )
         *tightened = FALSE;
      return SCIP_OKAY;
   }

   switch( scip->stage )
   {
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarChgLbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarChgLbGlobal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, 
            infercons, NULL, inferinfo, FALSE) );
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 */
RETCODE SCIPinferVarUbCons(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new value for bound */
   CONS*            infercons,          /**< constraint that deduced the bound change */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   Real lb;
   Real ub;

   assert(infeasible != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPinferVarUbCons", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));
   
   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   if( !SCIPsetIsUbBetter(scip->set, newbound, ub) )
   {
      if( tightened != NULL )
         *tightened = FALSE;
      return SCIP_OKAY;
   }

   switch( scip->stage )
   {
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarChgUbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarChgUbGlobal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, 
            infercons, NULL, inferinfo, FALSE) );
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** depending on SCIP's stage, fixes binary variable in the problem, in preprocessing, or in current node;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason for the
 *  deduction of the fixing
 */
RETCODE SCIPinferBinvarCons(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< binary variable to fix */
   Bool             fixedval,           /**< value to fix binary variable to */
   CONS*            infercons,          /**< constraint that deduced the fixing */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   Bool*            tightened           /**< pointer to store whether the fixing tightened the local bounds, or NULL */
   )
{
   Real lb;
   Real ub;

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(fixedval == TRUE || fixedval == FALSE);
   assert(infeasible != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPinferBinvarCons", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsEQ(scip->set, lb, 0.0) || SCIPsetIsEQ(scip->set, lb, 1.0));
   assert(SCIPsetIsEQ(scip->set, ub, 0.0) || SCIPsetIsEQ(scip->set, ub, 1.0));
   assert(SCIPsetIsLE(scip->set, lb, ub));

   /* check, if variable is already fixed */
   if( (lb > 0.5) || (ub < 0.5) )
   {
      *infeasible = (fixedval == (lb < 0.5));

      return SCIP_OKAY;
   }

   /* apply the fixing */
   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      if( fixedval == TRUE )
      {
         CHECK_OKAY( SCIPchgVarLb(scip, var, 1.0) );
      }
      else
      {
         CHECK_OKAY( SCIPchgVarUb(scip, var, 0.0) );
      }
      break;

   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarFix(var, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->primal, scip->lp,
                     scip->branchcand, scip->eventqueue, (Real)fixedval, infeasible) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      if( fixedval == TRUE )
      {
         CHECK_OKAY( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, 1.0, SCIP_BOUNDTYPE_LOWER, 
               infercons, NULL, inferinfo, FALSE) );
      }
      else
      {
         CHECK_OKAY( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, 0.0, SCIP_BOUNDTYPE_UPPER, 
               infercons, NULL, inferinfo, FALSE) );
      }
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 */
RETCODE SCIPinferVarLbProp(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new value for bound */
   PROP*            inferprop,          /**< propagator that deduced the bound change */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   Real lb;
   Real ub;

   assert(infeasible != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPinferVarLbProp", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));
   
   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);

   if( !SCIPsetIsLbBetter(scip->set, newbound, lb) )
   {
      if( tightened != NULL )
         *tightened = FALSE;
      return SCIP_OKAY;
   }

   switch( scip->stage )
   {
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarChgLbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarChgLbGlobal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, 
            NULL, inferprop, inferinfo, FALSE) );
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 */
RETCODE SCIPinferVarUbProp(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound,           /**< new value for bound */
   PROP*            inferprop,          /**< propagator that deduced the bound change */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   Real lb;
   Real ub;

   assert(infeasible != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPinferVarUbProp", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));
   
   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   if( !SCIPsetIsUbBetter(scip->set, newbound, ub) )
   {
      if( tightened != NULL )
         *tightened = FALSE;
      return SCIP_OKAY;
   }

   switch( scip->stage )
   {
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarChgUbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarChgUbGlobal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, 
            NULL, inferprop, inferinfo, FALSE) );
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** depending on SCIP's stage, fixes binary variable in the problem, in preprocessing, or in current node;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason for the
 *  deduction of the fixing
 */
RETCODE SCIPinferBinvarProp(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< binary variable to fix */
   Bool             fixedval,           /**< value to fix binary variable to */
   PROP*            inferprop,          /**< propagator that deduced the fixing */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   Bool*            tightened           /**< pointer to store whether the fixing tightened the local bounds, or NULL */
   )
{
   Real lb;
   Real ub;

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(fixedval == TRUE || fixedval == FALSE);
   assert(infeasible != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPinferBinvarProp", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsEQ(scip->set, lb, 0.0) || SCIPsetIsEQ(scip->set, lb, 1.0));
   assert(SCIPsetIsEQ(scip->set, ub, 0.0) || SCIPsetIsEQ(scip->set, ub, 1.0));
   assert(SCIPsetIsLE(scip->set, lb, ub));

   /* check, if variable is already fixed */
   if( (lb > 0.5) || (ub < 0.5) )
   {
      *infeasible = (fixedval == (lb < 0.5));

      return SCIP_OKAY;
   }

   /* apply the fixing */
   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      if( fixedval == TRUE )
      {
         CHECK_OKAY( SCIPchgVarLb(scip, var, 1.0) );
      }
      else
      {
         CHECK_OKAY( SCIPchgVarUb(scip, var, 0.0) );
      }
      break;

   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarFix(var, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->primal, scip->lp,
                     scip->branchcand, scip->eventqueue, (Real)fixedval, infeasible) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      if( fixedval == TRUE )
      {
         CHECK_OKAY( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, 1.0, SCIP_BOUNDTYPE_LOWER, 
               NULL, inferprop, inferinfo, FALSE) );
      }
      else
      {
         CHECK_OKAY( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, 0.0, SCIP_BOUNDTYPE_UPPER, 
               NULL, inferprop, inferinfo, FALSE) );
      }
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** informs variable x about a globally valid variable lower bound x >= b*z + d with integer variable z */
RETCODE SCIPaddVarVlb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   VAR*             vlbvar,             /**< variable z    in x >= b*z + d */
   Real             vlbcoef,            /**< coefficient b in x >= b*z + d */
   Real             vlbconstant         /**< constant d    in x >= b*z + d */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddVarVlb", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarAddVlb(var, scip->mem->solvemem, scip->set, vlbvar, vlbcoef, vlbconstant) );

   return SCIP_OKAY;
}

/** informs variable x about a globally valid variable upper bound x <= b*z + d with integer variable z */
RETCODE SCIPaddVarVub(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   VAR*             vubvar,             /**< variable z    in x <= b*z + d */
   Real             vubcoef,            /**< coefficient b in x <= b*z + d */
   Real             vubconstant         /**< constant d    in x <= b*z + d */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddVarVub", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarAddVub(var, scip->mem->solvemem, scip->set, vubvar, vubcoef, vubconstant) );

   return SCIP_OKAY;
}

/** informs variable x about a globally valid implication:  x >= b   =>   z <= c  or  z >= c */
RETCODE SCIPaddVarLbimpl(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             bound,              /**< bound b       bounding information    x >= b */
   VAR*             infervar,           /**< variable z    in inference            z <= c  or  z >= c */
   Bool             infertype,          /**< type          of inference    TRUE if z <= c, FALSE if z >= c */
   Real             inferbound          /**< bound c       in inference            z <= c  or  z >= c */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddVarLbimplic", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarAddLbimplic(var, scip->mem->solvemem, scip->set, bound, infervar, infertype, inferbound) );
   
   return SCIP_OKAY;
}

/** informs variable x about a globally valid implication:  x <= b   =>   z <= c  or  z >= c */
RETCODE SCIPaddVarUbimpl(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             bound,              /**< bound b       bounding information    x <= b */
   VAR*             infervar,           /**< variable z    in inference            z <= c  or  z >= c */
   Bool             infertype,          /**< type          of inference    TRUE if z <= c, FALSE if z >= c */
   Real             inferbound          /**< bound c       in inference            z <= c  or  z >= c */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddVarUbimplic", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarAddUbimplic(var, scip->mem->solvemem, scip->set, bound, infervar, infertype, inferbound) );

   return SCIP_OKAY;
}

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
RETCODE SCIPchgVarBranchFactor(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             branchfactor        /**< factor to weigh variable's branching score with */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchFactor(var, scip->set, branchfactor);

   return SCIP_OKAY;
}

/** scales the branch factor of the variable with the given value */
RETCODE SCIPscaleVarBranchFactor(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             scale               /**< factor to scale variable's branching factor with */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPscaleVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchFactor(var, scip->set, scale * SCIPvarGetBranchFactor(var));

   return SCIP_OKAY;
}

/** adds the given value to the branch factor of the variable */
RETCODE SCIPaddVarBranchFactor(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             addfactor           /**< value to add to the branch factor of the variable */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchFactor(var, scip->set, addfactor + SCIPvarGetBranchFactor(var));

   return SCIP_OKAY;
}

/** sets the branch priority of the variable; variables with higher branch priority are always prefered to variables
 *  with lower priority in selection of branching variable
 */
RETCODE SCIPchgVarBranchPriority(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   int              branchpriority      /**< branch priority of the variable */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchPriority(var, branchpriority);

   return SCIP_OKAY;
}

/** changes the branch priority of the variable to the given value, if it is larger than the current priority */
RETCODE SCIPupdateVarBranchPriority(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   int              branchpriority      /**< new branch priority of the variable, if it is larger than current priority */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPupdateVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( branchpriority > SCIPvarGetBranchPriority(var) )
      SCIPvarChgBranchPriority(var, branchpriority);

   return SCIP_OKAY;
}

/** adds the given value to the branch priority of the variable */
RETCODE SCIPaddVarBranchPriority(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   int              addpriority         /**< value to add to the branch priority of the variable */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchPriority(var, addpriority + SCIPvarGetBranchPriority(var));

   return SCIP_OKAY;
}

/** sets the branch direction of the variable (-1: prefer downwards branch, 0: automatic selection, +1: prefer upwards
 *  branch)
 */
RETCODE SCIPchgVarBranchDirection(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   BRANCHDIR        branchdirection     /**< preferred branch direction of the variable (downwards, upwards, auto) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarBranchDirection", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchDirection(var, branchdirection);

   return SCIP_OKAY;
}

/** changes type of variable in the problem; this changes the vars array returned from
 *  SCIPgetVars() and SCIPgetVarsData()
 */
RETCODE SCIPchgVarType(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   VARTYPE          vartype             /**< new type of variable */
   )
{
   assert(var != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPchgVarType", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      if( SCIPvarGetProbindex(var) >= 0 )
      {
         CHECK_OKAY( SCIPprobChgVarType(scip->origprob, scip->set, scip->branchcand, var, vartype) );
      }
      else
      {
         CHECK_OKAY( SCIPvarChgType(var, vartype) );
      }
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
      if( !SCIPvarIsTransformed(var) )
      {
         errorMessage("cannot change type of original variables while solving the problem\n");
         return SCIP_INVALIDCALL;
      }
      if( SCIPvarGetProbindex(var) >= 0 )
      {
         CHECK_OKAY( SCIPprobChgVarType(scip->transprob, scip->set, scip->branchcand, var, vartype) );
      }
      else
      {
         CHECK_OKAY( SCIPvarChgType(var, vartype) );
      }
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/   
}

/** in problem creation and solving stage, both bounds of the variable are set to the given value;
 *  in presolving stage, the variable is converted into a fixed variable, and bounds are changed respectively;
 *  conversion into a fixed variable changes the vars array returned from SCIPgetVars() and SCIPgetVarsData()
 */
RETCODE SCIPfixVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to fix */
   Real             fixedval,           /**< value to fix variable to */
   Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   Bool*            fixed               /**< pointer to store whether the fixing was performed (variable was unfixed) */
   )
{
   assert(var != NULL);
   assert(infeasible != NULL);
   assert(fixed != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfixVar", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *fixed = FALSE;

   if( (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !SCIPsetIsIntegral(scip->set, fixedval))
      || SCIPsetIsFeasLT(scip->set, fixedval, SCIPvarGetLbLocal(var))
      || SCIPsetIsFeasGT(scip->set, fixedval, SCIPvarGetUbLocal(var)) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
   {
      *infeasible = !SCIPsetIsFeasEQ(scip->set, fixedval, SCIPvarGetLbLocal(var));
      return SCIP_OKAY;
   }

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      *infeasible = FALSE;
      if( SCIPsetIsFeasGT(scip->set, fixedval, SCIPvarGetLbLocal(var)) )
      {
         CHECK_OKAY( SCIPchgVarLb(scip, var, fixedval) );
         *fixed = TRUE;
      }
      if( SCIPsetIsFeasLT(scip->set, fixedval, SCIPvarGetUbLocal(var)) )
      {
         CHECK_OKAY( SCIPchgVarUb(scip, var, fixedval) );
         *fixed = TRUE;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarFix(var, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->primal, scip->lp,
            scip->branchcand, scip->eventqueue, fixedval, infeasible) );
      *fixed = TRUE;
      return SCIP_OKAY;
      
   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** Tries to aggregate an equality a*x + b*y == c consisting of two integral active problem variables x and y.
 *  An integer aggregation (i.e. integral coefficients a' and b', such that a'*x + b'*y == c') is searched.
 *  This can lead to the detection of infeasibility (e.g. if c' is fractional), or to a rejection of the
 *  aggregation (denoted by aggregated == FALSE), if the resulting integer coefficients are too large and thus
 *  numerically instable.
 */
static
RETCODE aggregateActiveIntVars(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             varx,               /**< integral variable x in equality a*x + b*y == c */
   VAR*             vary,               /**< integral variable y in equality a*x + b*y == c */
   Real             scalarx,            /**< multiplier a in equality a*x + b*y == c */
   Real             scalary,            /**< multiplier b in equality a*x + b*y == c */
   Real             rhs,                /**< right hand side c in equality a*x + b*y == c */
   Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   )
{
   VAR* aggvar;
   char aggvarname[MAXSTRLEN];
   Longint scalarxn = 0;
   Longint scalarxd = 0;
   Longint scalaryn = 0;
   Longint scalaryd = 0;
   Longint a;
   Longint b;
   Longint c;
   Longint scm;
   Longint gcd;
   Longint currentclass;
   Longint classstep;
   Longint xsol;
   Longint ysol;
   Bool success;

#define MAXDNOM 1000000LL

   assert(scip->stage == SCIP_STAGE_PRESOLVING);
   assert(varx != NULL);
   assert(SCIPvarGetStatus(varx) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetType(varx) != SCIP_VARTYPE_CONTINUOUS);
   assert(vary != NULL);
   assert(SCIPvarGetStatus(vary) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetType(vary) != SCIP_VARTYPE_CONTINUOUS);
   assert(!SCIPsetIsZero(scip->set, scalarx));
   assert(!SCIPsetIsZero(scip->set, scalary));
   assert(infeasible != NULL);
   assert(aggregated != NULL);

   *infeasible = FALSE;
   *aggregated = FALSE;

   /* get rational representation of coefficients */
   success = SCIPrealToRational(scalarx, -SCIPsetEpsilon(scip->set), SCIPsetEpsilon(scip->set), MAXDNOM,
      &scalarxn, &scalarxd);
   if( success )
      success = SCIPrealToRational(scalary, -SCIPsetEpsilon(scip->set), SCIPsetEpsilon(scip->set), MAXDNOM,
         &scalaryn, &scalaryd);
   if( !success )
      return SCIP_OKAY;
   assert(scalarxd >= 1);
   assert(scalaryd >= 1);

   /* multiply equality with smallest common denominator */
   scm = SCIPcalcSmaComMul(scalarxd, scalaryd);
   a = (scm/scalarxd)*scalarxn;
   b = (scm/scalaryd)*scalaryn;
   rhs *= scm;
   
   /* divide equality by the greatest common divisor of a and b */
   gcd = SCIPcalcGreComDiv(ABS(a), ABS(b));
   a /= gcd;
   b /= gcd;
   rhs /= gcd;
   assert(a != 0);
   assert(b != 0);

   /* check, if right hand side is integral */
   if( !SCIPsetIsIntegral(scip->set, rhs) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   c = (Longint)(SCIPfloor(scip, rhs));

   /* check, if we are in an easy case with either |a| = 1 or |b| = 1 */
   if( a == 1 || a == -1 )
   {
      /* aggregate x = - b/a*y + c/a */
      /*lint --e{653}*/
      CHECK_OKAY( SCIPvarAggregate(varx, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
            scip->primal, scip->lp, scip->branchcand, scip->eventqueue, 
            vary, (Real)(-b/a), (Real)(c/a), infeasible) );
      *aggregated = TRUE;
      return SCIP_OKAY;
   }
   if( b == 1 || b == -1 )
   {
      /* aggregate y = - a/b*x + c/b */
      /*lint --e{653}*/
      CHECK_OKAY( SCIPvarAggregate(vary, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
            scip->primal, scip->lp, scip->branchcand, scip->eventqueue, 
            varx, (Real)(-a/b), (Real)(c/b), infeasible) );
      *aggregated = TRUE;
      return SCIP_OKAY;
   }

   /* Both variables are integers, their coefficients are not multiples of each other, and they don't have any
    * common divisor. Let (x',y') be a solution of the equality
    *   a*x + b*y == c    ->   a*x == c - b*y
    * Then x = -b*z + x', y = a*z + y' with z integral gives all solutions to the equality.
    */

   /* find initial solution (x',y'):
    *  - find y' such that c - b*y' is a multiple of a
    *    - start in equivalence class c%a
    *    - step through classes, where each step increases class number by (-b)%a, until class 0 is visited
    *    - if equivalence class 0 is visited, we are done: y' equals the number of steps taken
    *    - because a and b don't have a common divisor, each class is visited at most once, and at most a-1 steps are needed
    *  - calculate x' with x' = (c - b*y')/a (which must be integral)
    */

   /* search upwards from ysol = 0 */
   ysol = 0;
   currentclass = c%a;
   if( currentclass < 0 )
      currentclass += a;
   assert(0 <= currentclass && currentclass < a);
   classstep = (-b)%a;
   if( classstep < 0 )
      classstep += a;
   assert(0 < classstep && classstep < a);
   while( currentclass != 0 )
   {
      assert(0 <= currentclass && currentclass < a);
      currentclass += classstep;
      if( currentclass >= a )
         currentclass -= a;
      ysol++;
   }
   assert(ysol < a);
   assert(((c - b*ysol)%a) == 0);
   xsol = (c - b*ysol)/a;

   /* feasible solutions are (x,y) = (x',y') + z*(-b,a)
    * - create new integer variable z with infinite bounds
    * - aggregate variable x = -b*z + x'
    * - aggregate variable y =  a*z + y'
    * - the bounds of z are calculated automatically during aggregation
    */
   sprintf(aggvarname, "agg%d", scip->stat->nvaridx);
   CHECK_OKAY( SCIPvarCreateTransformed(&aggvar, scip->mem->solvemem, scip->set, scip->stat,
         aggvarname, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER,
         SCIPvarIsInitial(varx) || SCIPvarIsInitial(vary), SCIPvarIsRemoveable(varx) && SCIPvarIsRemoveable(vary),
         NULL, NULL, NULL, NULL) );
   CHECK_OKAY( SCIPprobAddVar(scip->transprob, scip->mem->solvemem, scip->set, scip->lp, 
         scip->branchcand, scip->eventfilter, scip->eventqueue, aggvar) );
   CHECK_OKAY( SCIPvarAggregate(varx, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
         scip->primal, scip->lp, scip->branchcand, scip->eventqueue,
         aggvar, (Real)(-b), (Real)xsol, infeasible) );
   if( !(*infeasible) )
   {
      CHECK_OKAY( SCIPvarAggregate(vary, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
            scip->primal, scip->lp, scip->branchcand, scip->eventqueue, 
            aggvar, (Real)a, (Real)ysol, infeasible) );
   }
   *aggregated = TRUE;

   /* release z */
   CHECK_OKAY( SCIPvarRelease(&aggvar, scip->mem->solvemem, scip->set, scip->lp) );

   return SCIP_OKAY;
}

/** performs second step of SCIPaggregateVars(): 
 *  the variable to be aggregated is chosen among active problem variables x' and y', prefering a less strict variable
 *  type as aggregation variable (i.e. continuous variables are prefered over implicit integers, implicit integers
 *  over integers, and integers over binaries). If none of the variables is continuous, it is tried to find an integer
 *  aggregation (i.e. integral coefficients a'' and b'', such that a''*x' + b''*y' == c''). This can lead to
 *  the detection of infeasibility (e.g. if c'' is fractional), or to a rejection of the aggregation (denoted by
 *  aggregated == FALSE), if the resulting integer coefficients are too large and thus numerically instable.
 */
static
RETCODE aggregateActiveVars(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             varx,               /**< variable x in equality a*x + b*y == c */
   VAR*             vary,               /**< variable y in equality a*x + b*y == c */
   Real             scalarx,            /**< multiplier a in equality a*x + b*y == c */
   Real             scalary,            /**< multiplier b in equality a*x + b*y == c */
   Real             rhs,                /**< right hand side c in equality a*x + b*y == c */
   Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   )
{
   int agg;

   assert(scip->stage == SCIP_STAGE_PRESOLVING);
   assert(varx != NULL);
   assert(SCIPvarGetStatus(varx) == SCIP_VARSTATUS_LOOSE);
   assert(vary != NULL);
   assert(SCIPvarGetStatus(vary) == SCIP_VARSTATUS_LOOSE);
   assert(!SCIPsetIsZero(scip->set, scalarx));
   assert(!SCIPsetIsZero(scip->set, scalary));
   assert(infeasible != NULL);
   assert(aggregated != NULL);

   *infeasible = FALSE;
   *aggregated = FALSE;

   /* figure out, which variable should be aggregated */
   agg = -1;
   
   /* a*x + b*y == c
    *  ->  x == -b/a * y + c/a  (agg=0)
    *  ->  y == -a/b * x + c/b  (agg=1)
    */
   if( SCIPvarGetType(varx) == SCIP_VARTYPE_CONTINUOUS )
      agg = 0;
   else if( SCIPvarGetType(vary) == SCIP_VARTYPE_CONTINUOUS )
      agg = 1;
   else if( SCIPvarGetType(varx) == SCIP_VARTYPE_IMPLINT )
      agg = 0;
   else if( SCIPvarGetType(vary) == SCIP_VARTYPE_IMPLINT )
      agg = 1;
   else if( SCIPsetIsIntegral(scip->set, scalary/scalarx) )
      agg = 0;
   else if( SCIPsetIsIntegral(scip->set, scalarx/scalary) )
      agg = 1;
   if( agg == 1 )
   {
      VAR* var;
      Real scalar;

      /* switch the variables, such that varx is the aggregated variable */
      var = vary;
      vary = varx;
      varx = var;
      scalar = scalary;
      scalary = scalarx;
      scalarx = scalar;
      agg = 0;
   }
   assert(agg == 0 || agg == -1);

   /* did we find an "easy" aggregation? */
   if( agg == 0 )
   {
      Real scalar;
      Real constant;

      /* calculate aggregation scalar and constant: a*x + b*y == c  =>  x == -b/a * y + c/a */
      scalar = -scalary/scalarx;
      constant = rhs/scalarx;

      /* check aggregation for integer feasibility */
      if( SCIPvarGetType(varx) != SCIP_VARTYPE_CONTINUOUS
         && SCIPvarGetType(vary) != SCIP_VARTYPE_CONTINUOUS
         && SCIPsetIsIntegral(scip->set, scalar) && !SCIPsetIsIntegral(scip->set, constant) )
      {
         *infeasible = TRUE;
         return SCIP_OKAY;
      }

      /* aggregate the variable */
      CHECK_OKAY( SCIPvarAggregate(varx, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, 
            scip->primal, scip->lp, scip->branchcand, scip->eventqueue, vary, scalar, constant, infeasible) );
      *aggregated = TRUE;
      return SCIP_OKAY;
   }

   /* the variables are both integral: we have to try to find an integer aggregation */
   CHECK_OKAY( aggregateActiveIntVars(scip, varx, vary, scalarx, scalary, rhs, infeasible, aggregated) );

   return SCIP_OKAY;
}

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
   )
{
   Real constantx;
   Real constanty;

   assert(infeasible != NULL);
   assert(redundant != NULL);
   assert(aggregated != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaggregateVars", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   *redundant = FALSE;
   *aggregated = FALSE;

   /* we cannot aggregate multiaggregated variables */
   if( SCIPvarGetStatus(varx) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(vary) == SCIP_VARSTATUS_MULTAGGR )
      return SCIP_OKAY;

   /* get the corresponding equality in active problem variable space:
    * transform both expressions "a*x + 0" and "b*y + 0" into problem variable space
    */
   constantx = 0.0;
   constanty = 0.0;
   CHECK_OKAY( SCIPvarGetProbvarSum(&varx, &scalarx, &constantx) );
   CHECK_OKAY( SCIPvarGetProbvarSum(&vary, &scalary, &constanty) );
   
   /* move the constant to the right hand side to acquire the form "a'*x' + b'*y' == c'" */
   rhs -= (constantx + constanty);

   /* if a scalar is zero, treat the variable as fixed-to-zero variable */
   if( SCIPsetIsZero(scip->set, scalarx) )
      varx = NULL;
   if( SCIPsetIsZero(scip->set, scalary) )
      vary = NULL;

   /* capture the special cases that less than two variables are left, due to resolutions to a fixed variable or
    * to the same active variable
    */
   if( varx == NULL && vary == NULL )
   {
      /* both variables were resolved to fixed variables */
      *infeasible = !SCIPsetIsZero(scip->set, rhs);
      *redundant = TRUE;
   }
   else if( varx == NULL )
   {
      assert(SCIPsetIsZero(scip->set, scalarx));
      assert(!SCIPsetIsZero(scip->set, scalary));
      
      /* variable x was resolved to fixed variable: variable y can be fixed to c'/b' */
      CHECK_OKAY( SCIPvarFix(vary, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->primal, scip->lp,
            scip->branchcand, scip->eventqueue, rhs/scalary, infeasible) );
      *aggregated = TRUE;
      *redundant = TRUE;
   }
   else if( vary == NULL )
   {
      assert(SCIPsetIsZero(scip->set, scalary));
      assert(!SCIPsetIsZero(scip->set, scalarx));
      
      /* variable y was resolved to fixed variable: variable x can be fixed to c'/a' */
      CHECK_OKAY( SCIPvarFix(varx, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->primal, scip->lp,
            scip->branchcand, scip->eventqueue, rhs/scalarx, infeasible) );
      *aggregated = TRUE;
      *redundant = TRUE;
   }
   else if( varx == vary )
   {
      /* both variables were resolved to the same active problem variable: this variable can be fixed */
      scalarx += scalary;
      if( SCIPsetIsZero(scip->set, scalarx) )
      {
         /* left hand side of equality is zero: equality is potentially infeasible */
         *infeasible = !SCIPsetIsZero(scip->set, rhs);
      }
      else
      {
         /* sum of scalars is not zero: fix variable x' == y' to c'/(a'+b') */
         CHECK_OKAY( SCIPvarFix(varx, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->primal, scip->lp,
               scip->branchcand, scip->eventqueue, rhs/scalarx, infeasible) );
         *aggregated = TRUE;
      }
      *redundant = TRUE;
   }
   else
   {
      /* both variables are different active problem variables, and both scalars are non-zero: try to aggregate them */
      CHECK_OKAY( aggregateActiveVars(scip, varx, vary, scalarx, scalary, rhs, infeasible, aggregated) );
      *redundant = *aggregated;
   }

   return SCIP_OKAY;
}

/** converts variable into multi-aggregated variable; this changes the vars array returned from
 *  SCIPgetVars() and SCIPgetVarsData(); Warning! The integrality condition is not checked anymore on
 *  the multiaggregated variable. You must not multiaggregate an integer variable without being sure,
 *  that integrality on the aggregation variables implies integrality on the aggregated variable.
 */
RETCODE SCIPmultiaggregateVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable x to aggregate */
   int              naggvars,           /**< number n of variables in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   VAR**            aggvars,            /**< variables y_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Real*            scalars,            /**< multipliers a_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Real             constant,           /**< constant shift c in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPmultiaggregateVar", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarMultiaggregate(var, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->primal, 
         scip->lp, scip->branchcand, scip->eventqueue, naggvars, aggvars, scalars, constant, infeasible) );

   return SCIP_OKAY;
}

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of "solvaldelta" in the
 *  variable's solution value and resulting change of "objdelta" in the in the LP's objective value
 */
RETCODE SCIPupdateVarPseudocost(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   Real             weight              /**< weight in (0,1] of this update in pseudo cost sum */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPupdateVarPseudocost", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarUpdatePseudocost(var, scip->set, scip->stat, solvaldelta, objdelta, weight) );

   return SCIP_OKAY;
}

/** gets the variable's pseudo cost value for the given direction */
Real SCIPgetVarPseudocost(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarPseudocost", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetPseudocost(var, scip->stat, solvaldelta);
}

/** gets the variable's pseudo cost value for the given direction,
 *  only using the pseudo cost information of the current run
 */
Real SCIPgetVarPseudocostCurrentRun(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarPseudocostCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetPseudocostCurrentRun(var, scip->stat, solvaldelta);
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
Real SCIPgetVarPseudocostCount(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarPseudocostCount", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetPseudocostCount(var, dir);
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 */
Real SCIPgetVarPseudocostCountCurrentRun(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarPseudocostCountCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetPseudocostCountCurrentRun(var, dir);
}

/** gets the variable's pseudo cost score value for the given LP solution value */
Real SCIPgetVarPseudocostScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             solval              /**< variable's LP solution value */
   )
{
   Real downsol;
   Real upsol;
   Real pscostdown;
   Real pscostup;

   CHECK_ABORT( checkStage(scip, "SCIPgetVarPseudocostScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   downsol = SCIPsetCeil(scip->set, solval-1.0);
   upsol = SCIPsetFloor(scip->set, solval+1.0);
   pscostdown = SCIPvarGetPseudocost(var, scip->stat, downsol-solval);
   pscostup = SCIPvarGetPseudocost(var, scip->stat, upsol-solval);

   return SCIPbranchGetScore(scip->set, var, pscostdown, pscostup);
}

/** gets the variable's pseudo cost score value for the given LP solution value,
 *  only using the pseudo cost information of the current run
 */
Real SCIPgetVarPseudocostScoreCurrentRun(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             solval              /**< variable's LP solution value */
   )
{
   Real downsol;
   Real upsol;
   Real pscostdown;
   Real pscostup;

   CHECK_ABORT( checkStage(scip, "SCIPgetVarPseudocostScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   downsol = SCIPsetCeil(scip->set, solval-1.0);
   upsol = SCIPsetFloor(scip->set, solval+1.0);
   pscostdown = SCIPvarGetPseudocostCurrentRun(var, scip->stat, downsol-solval);
   pscostup = SCIPvarGetPseudocostCurrentRun(var, scip->stat, upsol-solval);

   return SCIPbranchGetScore(scip->set, var, pscostdown, pscostup);
}

/** returns the average number of inferences found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 */
Real SCIPgetVarAvgInferences(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarAvgInferences", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetAvgInferences(var, scip->stat, dir);
}

/** returns the average number of inferences found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 */
Real SCIPgetVarAvgInferencesCurrentRun(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarAvgInferencesCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, dir);
}

/** returns the variable's average inference score value */
Real SCIPgetVarAvgInferenceScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< problem variable */
   )
{
   Real inferdown;
   Real inferup;

   CHECK_ABORT( checkStage(scip, "SCIPgetVarAvgInferenceScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   
   inferdown = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, inferdown, inferup);
}

/** returns the variable's average inference score value only using inferences of the current run */
Real SCIPgetVarAvgInferenceScoreCurrentRun(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< problem variable */
   )
{
   Real inferdown;
   Real inferup;

   CHECK_ABORT( checkStage(scip, "SCIPgetVarAvgInferenceScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   
   inferdown = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, inferdown, inferup);
}

/** returns the average number of cutoffs found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 */
Real SCIPgetVarAvgCutoffs(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarAvgCutoffs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetAvgCutoffs(var, scip->stat, dir);
}

/** returns the average number of cutoffs found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 */
Real SCIPgetVarAvgCutoffsCurrentRun(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarAvgCutoffsCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, dir);
}

/** returns the variable's average cutoff score value */
Real SCIPgetVarAvgCutoffScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< problem variable */
   )
{
   Real cutoffdown;
   Real cutoffup;

   CHECK_ABORT( checkStage(scip, "SCIPgetVarAvgCutoffScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   
   cutoffdown = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, cutoffdown, cutoffup);
}

/** returns the variable's average cutoff score value, only using cutoffs of the current run */
Real SCIPgetVarAvgCutoffScoreCurrentRun(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< problem variable */
   )
{
   Real cutoffdown;
   Real cutoffup;

   CHECK_ABORT( checkStage(scip, "SCIPgetVarAvgCutoffScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   
   cutoffdown = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, cutoffdown, cutoffup);
}

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor
 */
Real SCIPgetVarAvgInferenceCutoffScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   )
{
   Real avginferdown;
   Real avginferup;
   Real avginfer;
   Real inferdown;
   Real inferup;
   Real cutoffdown;
   Real cutoffup;

   CHECK_ABORT( checkStage(scip, "SCIPgetVarAvgInferenceCutoffScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   
   avginferdown = SCIPhistoryGetAvgInferences(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS);
   avginferup = SCIPhistoryGetAvgInferences(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS);
   avginfer = (avginferdown + avginferup)/2.0;
   inferdown = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);
   cutoffdown = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, 
      inferdown + cutoffweight * avginfer * cutoffdown, inferup + cutoffweight * avginfer * cutoffup);
}

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor, only using inferences and cutoffs of the current run
 */
Real SCIPgetVarAvgInferenceCutoffScoreCurrentRun(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   )
{
   Real avginferdown;
   Real avginferup;
   Real avginfer;
   Real inferdown;
   Real inferup;
   Real cutoffdown;
   Real cutoffup;

   CHECK_ABORT( checkStage(scip, "SCIPgetVarAvgInferenceCutoffScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   
   avginferdown = SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS);
   avginferup = SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS);
   avginfer = (avginferdown + avginferup)/2.0;
   inferdown = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);
   cutoffdown = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, 
      inferdown + cutoffweight * avginfer * cutoffdown, inferup + cutoffweight * avginfer * cutoffup);
}

/** gets user data for given variable */
VARDATA* SCIPgetVarData(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< problem variable */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPvarGetData(var);
}

/** outputs variable information to file stream */
RETCODE SCIPprintVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPvarPrint(var, scip->set, file);

   return SCIP_OKAY;
}




/*
 * conflict analysis methods
 */

/** initializes the conflict analysis by clearing the conflict candidate queue; this method must be called before
 *  you enter the conflict variables by calling SCIPaddConflictLb(), SCIPaddConflictUb(), or SCIPaddConflictBinvar();
 */
RETCODE SCIPinitConflictAnalysis(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPinitConflictAnalysis", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconflictInit(scip->conflict) );

   return SCIP_OKAY;
}

/** adds lower bound of variable at the time of the given bound change index to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictLb() should be called for each lower bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictLb() should be called
 *      for each lower bound, whose current assignment lead to the deduction of the given conflict bound.
 */
RETCODE SCIPaddConflictLb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable whose lower bound should be added to conflict candidate queue */
   BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                         *   conflicting bound was valid, NULL for current local bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddConflictLb", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconflictAddBound(scip->conflict, scip->mem->solvemem, scip->set, scip->stat,
         var, SCIP_BOUNDTYPE_LOWER, bdchgidx) );

   return SCIP_OKAY;
}

/** adds upper bound of variable at the time of the given bound change index to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictLb() should be called for each upper bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictUb() should be called
 *      for each upper bound, whose current assignment lead to the deduction of the given conflict bound.
 */
RETCODE SCIPaddConflictUb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                         *   conflicting bound was valid, NULL for current local bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddConflictUb", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconflictAddBound(scip->conflict, scip->mem->solvemem, scip->set, scip->stat,
         var, SCIP_BOUNDTYPE_UPPER, bdchgidx) );

   return SCIP_OKAY;
}

/** adds changed bound of fixed binary variable to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictBinvar() should be called for each fixed binary
 *      variable that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictBinvar() should be called
 *      for each binary variable, whose current fixing lead to the deduction of the given conflict bound.
 */
RETCODE SCIPaddConflictBinvar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< binary variable whose changed bound should be added to conflict queue */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddConflictBinvar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   if( SCIPvarGetLbLocal(var) > 0.5 )
   {
      CHECK_OKAY( SCIPconflictAddBound(scip->conflict, scip->mem->solvemem, scip->set, scip->stat,
            var, SCIP_BOUNDTYPE_LOWER, NULL) );
   }
   else if( SCIPvarGetUbLocal(var) < 0.5 )
   {
      CHECK_OKAY( SCIPconflictAddBound(scip->conflict, scip->mem->solvemem, scip->set, scip->stat,
            var, SCIP_BOUNDTYPE_UPPER, NULL) );
   }

   return SCIP_OKAY;
}

/** analyzes conflict bounds that were added after a call to SCIPinitConflictAnalysis() with calls to
 *  SCIPconflictAddLb(), SCIPconflictAddUb(), or SCIPaddConflictBinvar();
 *  on success, calls the conflict handlers to create a conflict constraint out of the resulting conflict set;
 *  the given valid depth must be a depth level, at which the conflict set defined by calls to SCIPaddConflictLb(),
 *  SCIPaddConflictUb() and SCIPaddConflictBinvar() is valid for the whole subtree; if the conflict was found by a
 *  violated constraint, use SCIPanalyzeConflictCons() instead of SCIPanalyzeConflict() to make sure, that the correct
 *  valid depth is used
 */
RETCODE SCIPanalyzeConflict(
   SCIP*            scip,               /**< SCIP data structure */
   int              validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPanalyzeConflict", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconflictAnalyze(scip->conflict, scip->mem->solvemem, scip->set, scip->stat, 
         scip->transprob, scip->tree, validdepth, success) );
   
   return SCIP_OKAY;
}

/** analyzes conflict bounds that were added with calls to SCIPconflictAddLb(), SCIPconflictAddUb(),
 *  or SCIPaddConflictBinvar(); on success, calls the conflict handlers to create a conflict constraint out of the
 *  resulting conflict set;
 *  the given constraint must be the constraint that detected the conflict, i.e. the constraint that is infeasible
 *  in the local bounds of the initial conflict set (defined by calls to SCIPaddConflictLb(), SCIPaddConflictUb(),
 *  and SCIPaddConflictBinvar())
 */
RETCODE SCIPanalyzeConflictCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint that detected the conflict */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPanalyzeConflictCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPconsIsGlobal(cons) )
   {
      CHECK_OKAY( SCIPconflictAnalyze(scip->conflict, scip->mem->solvemem, scip->set, scip->stat, 
            scip->transprob, scip->tree, 0, success) );
   }
   else if( SCIPconsIsActive(cons) )
   {
      CHECK_OKAY( SCIPconflictAnalyze(scip->conflict, scip->mem->solvemem, scip->set, scip->stat, 
            scip->transprob, scip->tree, SCIPconsGetActiveDepth(cons), success) );
   }
   
   return SCIP_OKAY;
}




/*
 * constraint methods
 */

/** creates and captures a constraint of the given constraint handler
 *  Warning! If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
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
   )
{
   assert(cons != NULL);
   assert(name != NULL);
   assert(conshdlr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcreateCons", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->probmem, name, conshdlr, consdata, 
            initial, separate, enforce, check, propagate, local, modifiable, removeable, TRUE) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->solvemem, name, conshdlr, consdata,
            initial, separate, enforce, check, propagate, local, modifiable, removeable, FALSE) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** increases usage counter of constraint */
RETCODE SCIPcaptureCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to capture */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

/** decreases usage counter of constraint, and frees memory if necessary */
RETCODE SCIPreleaseCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons                /**< pointer to constraint */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreleaseCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPconsRelease(cons, scip->mem->probmem, scip->set) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      if( SCIPconsIsOriginal(*cons) && (*cons)->nuses == 1 )
      {
         errorMessage("cannot release last use of original constraint while the transformed problem exists\n");
         return SCIP_INVALIDCALL;
      }
      CHECK_OKAY( SCIPconsRelease(cons, scip->mem->solvemem, scip->set) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets and captures transformed constraint of a given constraint; if the constraint is not yet transformed,
 *  a new transformed constraint for this constraint is created
 */
RETCODE SCIPtransformCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint to get/create transformed constraint for */
   CONS**           transcons           /**< pointer to store the transformed constraint */
   )
{
   assert(transcons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPtransformCons", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPconsIsTransformed(cons) )
   {
      *transcons = cons;
      SCIPconsCapture(*transcons);
   }
   else
   {
      CHECK_OKAY( SCIPconsTransform(cons, scip->mem->solvemem, scip->set, transcons) );
   }

   return SCIP_OKAY;
}

/** gets and captures transformed constraints for an array of constraints;
 *  if a constraint in the array is not yet transformed, a new transformed constraint for this constraint is created;
 *  it is possible to call this method with conss == transconss
 */
RETCODE SCIPtransformConss(
   SCIP*            scip,               /**< SCIP data structure */
   int              nconss,             /**< number of constraints to get/create transformed constraints for */
   CONS**           conss,              /**< array with constraints to get/create transformed constraints for */
   CONS**           transconss          /**< array to store the transformed constraints */
   )
{
   int c;

   assert(nconss == 0 || conss != NULL);
   assert(nconss == 0 || transconss != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPtransformConss", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsTransformed(conss[c]) )
      {
         transconss[c] = conss[c];
         SCIPconsCapture(transconss[c]);
      }
      else
      {
         CHECK_OKAY( SCIPconsTransform(conss[c], scip->mem->solvemem, scip->set, &transconss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** gets corresponding transformed constraint of a given constraint;
 *  returns NULL as transcons, if transformed constraint is not yet existing
 */
RETCODE SCIPgetTransformedCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint to get the transformed constraint for */
   CONS**           transcons           /**< pointer to store the transformed constraint */
   )
{
   assert(transcons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetTransformedCons", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPconsIsTransformed(cons) )
      *transcons = cons;
   else
      *transcons = SCIPconsGetTransformed(cons);

   return SCIP_OKAY;
}

/** gets corresponding transformed constraints for an array of constraints;
 *  stores NULL in a transconss slot, if the transformed constraint is not yet existing;
 *  it is possible to call this method with conss == transconss, but remember that constraints that are not
 *  yet transformed will be replaced with NULL
 */
RETCODE SCIPgetTransformedConss(
   SCIP*            scip,               /**< SCIP data structure */
   int              nconss,             /**< number of constraints to get the transformed constraints for */
   CONS**           conss,              /**< constraints to get the transformed constraints for */
   CONS**           transconss          /**< array to store the transformed constraints */
   )
{
   int c;

   assert(nconss == 0 || conss != NULL);
   assert(nconss == 0 || transconss != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetTransformedConss", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsTransformed(conss[c]) )
         transconss[c] = conss[c];
      else
         transconss[c] = SCIPconsGetTransformed(conss[c]);
   }

   return SCIP_OKAY;
}

/** adds given value to age of constraint, but age can never become negative;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 */
RETCODE SCIPaddConsAge(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint */
   Real             deltaage            /**< value to add to the constraint's age */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddConsAge", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsAddAge(cons, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, deltaage) );

   return SCIP_OKAY;
}

/** increases age of constraint by 1.0;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 */
RETCODE SCIPincConsAge(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPincConsAge", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsIncAge(cons, scip->mem->solvemem, scip->set, scip->stat, scip->transprob) );

   return SCIP_OKAY;
}

/** resets age of constraint to zero;
 *  should be called
 *   - in constraint separation, if a cut was found for this constraint,
 *   - in constraint enforcing, if the constraint was violated, and
 *   - in constraint propagation, if a domain reduction was deduced;
 */
RETCODE SCIPresetConsAge(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPresetConsAge", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsResetAge(cons, scip->set) );

   return SCIP_OKAY;
}

/** enables constraint's propagation capabilities */
RETCODE SCIPenableConsPropagation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPenableConsPropagation", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsEnablePropagation(cons, scip->set) );

   return SCIP_OKAY;
}

/** disables constraint's propagation capabilities s.t. the constraint is not propagated anymore until the propagation
 *  is enabled again with a call to SCIPenableConsPropagation()
 */
RETCODE SCIPdisableConsPropagation(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPdisableConsPropagation", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsDisablePropagation(cons, scip->set) );

   return SCIP_OKAY;
}

/** locks rounding of variables involved in the costraint */
RETCODE SCIPlockConsVars(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPlockConsVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsLockVars(cons, scip->set, nlockspos, nlocksneg) );

   return SCIP_OKAY;
}

/** unlocks rounding of variables involved in the costraint */
RETCODE SCIPunlockConsVars(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint */
   int              nunlockspos,        /**< decrease in number of rounding locks for constraint */
   int              nunlocksneg         /**< decrease in number of rounding locks for constraint's negation */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPunlockConsVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPconsUnlockVars(cons, scip->set, nunlockspos, nunlocksneg) );

   return SCIP_OKAY;
}

/** checks single constraint for feasibility of the given solution */
RETCODE SCIPcheckCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint to check */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcheckCons", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsCheck(cons, scip->set, sol, checkintegrality, checklprows, result) );

   return SCIP_OKAY;
}

/** marks the constraint to be essential for feasibility */
RETCODE SCIPsetConsChecked(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetConsChecked", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPconsSetChecked(cons, scip->set) );

   return SCIP_OKAY;
}

/** outputs constraint information to file stream */
RETCODE SCIPprintCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPconsPrint(cons, scip->set, file) );

   return SCIP_OKAY;
}




/*
 * LP methods
 */

/** returns, whether the LP was or is to be solved in the current node */
Bool SCIPhasCurrentNodeLP(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPhasCurrentNodeLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeHasCurrentNodeLP(scip->tree);
}

/** gets solution status of current LP */
LPSOLSTAT SCIPgetLPSolstat(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLPSolstat", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetSolstat(scip->lp);
}

/** gets objective value of current LP (which is the sum of column and loose objective value) */
Real SCIPgetLPObjval(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLPObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetObjval(scip->lp, scip->set);
}

/** gets part of objective value of current LP that results from COLUMN variables only */
Real SCIPgetLPColumnObjval(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLPColumnObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetColumnObjval(scip->lp);
}

/** gets part of objective value of current LP that results from LOOSE variables only */
Real SCIPgetLPLooseObjval(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLPLooseObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetLooseObjval(scip->lp, scip->set);
}

/** gets pseudo objective value of the current LP */
Real SCIPgetPseudoObjval(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPseudoObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   return SCIPlpGetPseudoObjval(scip->lp, scip->set);
}

/** gets current LP columns along with the current number of LP columns */
RETCODE SCIPgetLPColsData(
   SCIP*            scip,               /**< SCIP data structure */
   COL***           cols,               /**< pointer to store the array of LP columns, or NULL */
   int*             ncols               /**< pointer to store the number of LP columns, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPColsData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );   

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      if( cols != NULL )
         *cols = SCIPlpGetCols(scip->lp);
      if( ncols != NULL )
         *ncols = SCIPlpGetNCols(scip->lp);
   }
   else
   {
      if( cols != NULL )
         *cols = NULL;
      if( ncols != NULL )
         *ncols = scip->tree->pathnlpcols[scip->tree->pathlen-1];
   }

   return SCIP_OKAY;
}

/** gets current LP columns */
COL** SCIPgetLPCols(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLPCols", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );   

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIPlpGetCols(scip->lp);
   else
      return NULL;
}

/** gets current number of LP columns */
int SCIPgetNLPCols(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNLPCols", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );   

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIPlpGetNCols(scip->lp);
   else
      return scip->tree->pathnlpcols[scip->tree->pathlen-1];
}

/** gets current LP rows along with the current number of LP rows */
RETCODE SCIPgetLPRowsData(
   SCIP*            scip,               /**< SCIP data structure */
   ROW***           rows,               /**< pointer to store the array of LP rows, or NULL */
   int*             nrows               /**< pointer to store the number of LP rows, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPRowsData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );   
   
   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      if( rows != NULL )
         *rows = SCIPlpGetRows(scip->lp);
      if( nrows != NULL )
         *nrows = SCIPlpGetNRows(scip->lp);
   }
   else
   {
      if( rows != NULL )
         *rows = NULL;
      if( nrows != NULL )
         *nrows = scip->tree->pathnlprows[scip->tree->pathlen-1];
   }
   
   return SCIP_OKAY;
}

/** gets current LP rows */
ROW** SCIPgetLPRows(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );   
   
   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIPlpGetRows(scip->lp);
   else
      return NULL;
}

/** gets current number of LP rows */
int SCIPgetNLPRows(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );   

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIPlpGetNRows(scip->lp);
   else
      return scip->tree->pathnlprows[scip->tree->pathlen-1];
}

/** returns TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing
 */
Bool SCIPallColsInLP(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPallColsInLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp);
}

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1 */
RETCODE SCIPgetLPBasisInd(
   SCIP*            scip,               /**< SCIP data structure */
   int*             basisind            /**< pointer to store the basis indices */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPBasisInd", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPlpGetBasisInd(scip->lp, basisind) );

   return SCIP_OKAY;
}

/** gets a row from the inverse basis matrix B^-1 */
RETCODE SCIPgetLPBInvRow(
   SCIP*            scip,               /**< SCIP data structure */
   int              r,                  /**< row number */
   Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPBInvRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPlpGetBInvRow(scip->lp, r, coef) );

   return SCIP_OKAY;
}

/** gets a row from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A) */
RETCODE SCIPgetLPBInvARow(
   SCIP*            scip,               /**< SCIP data structure */
   int              r,                  /**< row number */
   Real*            binvrow,            /**< row in B^-1 from prior call to SCIPgetLPBInvRow(), or NULL */
   Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPBInvARow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPlpGetBInvARow(scip->lp, r, binvrow, coef) );

   return SCIP_OKAY;
}

/** calculates a weighted sum of all LP rows; for negative weights, the left and right hand side of the corresponding
 *  LP row are swapped in the summation
 */
RETCODE SCIPsumLPRows(
   SCIP*            scip,               /**< SCIP data structure */
   Real*            weights,            /**< row weights in row summation */
   REALARRAY*       sumcoef,            /**< array to store sum coefficients indexed by variables' probindex */
   Real*            sumlhs,             /**< pointer to store the left hand side of the row summation */
   Real*            sumrhs              /**< pointer to store the right hand side of the row summation */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsumLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPlpSumRows(scip->lp, scip->set, scip->transprob, weights, sumcoef, sumlhs, sumrhs) );

   return SCIP_OKAY;
}

/* calculates a MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
RETCODE SCIPcalcMIR(
   SCIP*            scip,               /**< SCIP data structure */
   Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   Bool*            cutislocal          /**< pointer to store whether the returned cut is only valid locally */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcalcMIR", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPlpCalcMIR(scip->lp, scip->set, scip->stat, scip->transprob,
         boundswitch, usevbds, allowlocal, maxweightrange, minfrac, weights, scale,
         mircoef, mirrhs, cutactivity, success, cutislocal) );

   return SCIP_OKAY;
}

/* calculates a strong CG cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
RETCODE SCIPcalcStrongCG(
   SCIP*            scip,               /**< SCIP data structure */
   Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real             scale,              /**< additional scaling factor multiplied to all rows */
   Real*            mircoef,            /**< array to store strong CG coefficients: must be of size SCIPgetNVars() */
   Real*            mirrhs,             /**< pointer to store the right hand side of the strong CG row */
   Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   Bool*            success,            /**< pointer to store whether the returned coefficients are a valid strong CG cut */
   Bool*            cutislocal          /**< pointer to store whether the returned cut is only valid locally */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcalcStrongCG", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPlpCalcStrongCG(scip->lp, scip->set, scip->stat, scip->transprob,
         boundswitch, usevbds, allowlocal, maxweightrange, minfrac, weights, scale,
         mircoef, mirrhs, cutactivity, success, cutislocal) );

   return SCIP_OKAY;
}

/** writes current LP to a file */
RETCODE SCIPwriteLP(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      fname               /**< file name */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPwriteLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   CHECK_OKAY( SCIPlpWrite(scip->lp, fname) );

   return SCIP_OKAY;
}




/*
 * LP row methods
 */

/** creates and captures an LP row */
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
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwCreate(row, scip->mem->solvemem, scip->set, scip->stat,
         name, len, cols, vals, lhs, rhs, local, modifiable, removeable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row without any coefficients */
RETCODE SCIPcreateEmptyRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row,                /**< pointer to row */
   const char*      name,               /**< name of row */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             local,              /**< is row only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateEmptyRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwCreate(row, scip->mem->solvemem, scip->set, scip->stat,
         name, 0, NULL, NULL, lhs, rhs, local, modifiable, removeable) );

   return SCIP_OKAY;
}

/** increases usage counter of LP row */
RETCODE SCIPcaptureRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< row to capture */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIProwCapture(row);

   return SCIP_OKAY;
}

/** decreases usage counter of LP row, and frees memory if necessary */
RETCODE SCIPreleaseRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row                 /**< pointer to LP row */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPreleaseRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE) );

   CHECK_OKAY( SCIProwRelease(row, scip->mem->solvemem, scip->set, scip->lp) );
   
   return SCIP_OKAY;
}

/** changes left hand side of LP row */
RETCODE SCIPchgRowLhs(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             lhs                 /**< new left hand side */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgRowLhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwChgLhs(row, scip->set, scip->lp, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of LP row */
RETCODE SCIPchgRowRhs(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             rhs                 /**< new right hand side */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgRowRhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwChgRhs(row, scip->set, scip->lp, rhs) );

   return SCIP_OKAY;
}

/** informs row, that all subsequent additions of variables to the row should be cached and not directly applied;
 *  after all additions were applied, SCIPflushRowExtensions() must be called;
 *  while the caching of row extensions is activated, information methods of the row give invalid results;
 *  caching should be used, if a row is build with SCIPaddVarToRow() calls variable by variable to increase
 *  the performance
 */
RETCODE SCIPcacheRowExtensions(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcacheRowExtension", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* delay the row sorting */
   SCIProwDelaySort(row);

   return SCIP_OKAY;
}

/** flushes all cached row extensions after a call of SCIPcacheRowExtensions() and merges coefficients with
 *  equal columns into a single coefficient
 */
RETCODE SCIPflushRowExtensions(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPflushRowExtension", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* force the row sorting, and merge equal column entries */
   SCIProwForceSort(row, scip->set);

   return SCIP_OKAY;
}

/** resolves variable to columns and adds them with the coefficient to the row */
RETCODE SCIPaddVarToRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   VAR*             var,                /**< problem variable */
   Real             val                 /**< value of coefficient */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddVarToRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarAddToRow(var, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp, row, val) );
   
   return SCIP_OKAY;
}

/** resolves variables to columns and adds them with the coefficients to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 */
RETCODE SCIPaddVarsToRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   int              nvars,              /**< number of variables to add to the row */
   VAR**            vars,               /**< problem variables to add */
   Real*            vals                /**< values of coefficients */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddVarsToRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* delay the row sorting */
   SCIProwDelaySort(row);

   /* add the variables to the row */
   for( v = 0; v < nvars; ++v )
   {
      CHECK_OKAY( SCIPvarAddToRow(vars[v], scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp,
            row, vals[v]) );
   }

   /* force the row sorting */
   SCIProwForceSort(row, scip->set);

   return SCIP_OKAY;
}

/** resolves variables to columns and adds them with the same single coefficient to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 */
RETCODE SCIPaddVarsToRowSameCoef(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   int              nvars,              /**< number of variables to add to the row */
   VAR**            vars,               /**< problem variables to add */
   Real             val                 /**< unique value of all coefficients */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddVarsToRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* resize the row to be able to store all variables (at least, if they are COLUMN variables) */
   CHECK_OKAY( SCIProwEnsureSize(row, scip->mem->solvemem, scip->set, nvars) );

   /* delay the row sorting */
   SCIProwDelaySort(row);

   /* add the variables to the row */
   for( v = 0; v < nvars; ++v )
   {
      CHECK_OKAY( SCIPvarAddToRow(vars[v], scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp,
            row, val) );
   }

   /* force the row sorting */
   SCIProwForceSort(row, scip->set);

   return SCIP_OKAY;
}

/** tries to find a value, such that all row coefficients, if scaled with this value become integral */
RETCODE SCIPcalcRowIntegralScalar(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   Real             maxscale,           /**< maximal allowed scalar */
   Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   Bool*            success             /**< stores whether returned value is valid */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcalcRowIntegralScalar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwCalcIntegralScalar(row, scip->set, mindelta, maxdelta, maxdnom, maxscale, 
         usecontvars, intscalar, success) );

   return SCIP_OKAY;
}

/** tries to scale row, s.t. all coefficients become integral */
RETCODE SCIPmakeRowIntegral(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   Real             maxscale,           /**< maximal value to scale row with */
   Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   Bool*            success             /**< stores whether row could be made rational */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPmakeRowIntegral", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwMakeIntegral(row, scip->set, scip->stat, scip->lp, mindelta, maxdelta, maxdnom, maxscale,
         usecontvars, success) );

   return SCIP_OKAY;
}

/** returns minimal absolute value of row vector's non-zero coefficients */
Real SCIPgetRowMinCoef(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowMinCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   return SCIProwGetMinval(row, scip->set);
}

/** returns maximal absolute value of row vector's non-zero coefficients */
Real SCIPgetRowMaxCoef(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowMaxCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   return SCIProwGetMaxval(row, scip->set);
}

/** returns the minimal activity of a row w.r.t. the column's bounds */
Real SCIPgetRowMinActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowMinActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   return SCIProwGetMinActivity(row, scip->set, scip->stat);
}

/** returns the maximal activity of a row w.r.t. the column's bounds */
Real SCIPgetRowMaxActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowMaxActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   return SCIProwGetMaxActivity(row, scip->set, scip->stat);
}

/** recalculates the activity of a row in the last LP solution */
RETCODE SCIPrecalcRowLPActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPrecalcRowLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   SCIProwRecalcLPActivity(row, scip->stat);

   return SCIP_OKAY;
}

/** returns the activity of a row in the last LP solution */
Real SCIPgetRowLPActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   return SCIProwGetLPActivity(row, scip->stat, scip->lp);
}

/** returns the feasibility of a row in the last LP solution: negative value means infeasibility */
Real SCIPgetRowLPFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowLPFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   return SCIProwGetLPFeasibility(row, scip->stat, scip->lp);
}

/** recalculates the activity of a row for the current pseudo solution */
RETCODE SCIPrecalcRowPseudoActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPrecalcRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   SCIProwRecalcPseudoActivity(row, scip->stat);

   return SCIP_OKAY;
}

/** returns the activity of a row for the current pseudo solution */
Real SCIPgetRowPseudoActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetPseudoActivity(row, scip->stat);
}

/** returns the feasibility of a row for the current pseudo solution: negative value means infeasibility */
Real SCIPgetRowPseudoFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowPseudoFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetPseudoFeasibility(row, scip->stat);
}

/** recalculates the activity of a row in the last LP or pseudo solution */
RETCODE SCIPrecalcRowActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPrecalcRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      SCIProwRecalcLPActivity(row, scip->stat);
   else
      SCIProwRecalcPseudoActivity(row, scip->stat);

   return SCIP_OKAY;
}

/** returns the activity of a row in the last LP or pseudo solution */
Real SCIPgetRowActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPActivity(row, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoActivity(row, scip->stat);
}

/** returns the feasibility of a row in the last LP or pseudo solution */
Real SCIPgetRowFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPFeasibility(row, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoFeasibility(row, scip->stat);
}

/** returns the activity of a row for the given primal solution */
Real SCIPgetRowSolActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   SOL*             sol                 /**< primal CIP solution */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowSolActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
   {
      Real activity;

      CHECK_ABORT( SCIProwGetSolActivity(row, scip->set, scip->stat, sol, &activity) );
      return activity;
   }
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPActivity(row, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoActivity(row, scip->stat);
}

/** returns the feasibility of a row for the given primal solution */
Real SCIPgetRowSolFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   SOL*             sol                 /**< primal CIP solution */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowSolFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
   {
      Real feasibility;

      CHECK_ABORT( SCIProwGetSolFeasibility(row, scip->set, scip->stat, sol, &feasibility) );
      return feasibility;
   }
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPFeasibility(row, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoFeasibility(row, scip->stat);
}

/** output row to file stream */
RETCODE SCIPprintRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPprintRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   SCIProwPrint(row, file);

   return SCIP_OKAY;
}




/*
 * cutting plane methods
 */

/** returns efficacy of the cut with respect to the current LP solution: e = -feasibility/norm */
Real SCIPgetCutEfficacy(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             cut                 /**< separated cut */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetCutEfficacy", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetEfficacy(cut, scip->set, scip->stat, scip->lp);
}

/** returns whether the cut's efficacy with respect to the current LP solution is greater than the minimal cut efficacy */
Bool SCIPisCutEfficacious(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             cut                 /**< separated cut */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPisCutEfficacious", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwIsEfficacious(cut, scip->set, scip->stat, scip->lp, (SCIPtreeGetCurrentDepth(scip->tree) == 0));
}

/** checks, if the given cut's efficacy is larger than the minimal cut efficacy */
Bool SCIPisEfficacious(
   SCIP*            scip,               /**< SCIP data structure */
   Real             efficacy            /**< efficacy of the cut */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPisCutEfficacious", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetIsEfficacious(scip->set, (SCIPtreeGetCurrentDepth(scip->tree) == 0), efficacy);
}

/** adds cut to separation storage;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
RETCODE SCIPaddCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             cut,                /**< separated cut */
   Bool             forcecut            /**< should the cut be forced to enter the LP? */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   assert(SCIPtreeGetCurrentNode(scip->tree) != NULL);

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      errorMessage("cannot add cuts, because node LP is not processed\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPsepastoreAddCut(scip->sepastore, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
         cut, forcecut, (SCIPtreeGetCurrentDepth(scip->tree) == 0)) );
   
   return SCIP_OKAY;
}

/** if not already existing, adds row to global cut pool */
RETCODE SCIPaddPoolCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< cutting plane to add */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddPoolCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPcutpoolAddRow(scip->cutpool, scip->mem->solvemem, scip->set, row) );

   return SCIP_OKAY;
}

/** removes the row from the global cut pool */
RETCODE SCIPdelPoolCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< cutting plane to add */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPdelPoolCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPcutpoolDelRow(scip->cutpool, scip->mem->solvemem, scip->set, scip->stat, scip->lp, row) );

   return SCIP_OKAY;
}

/** gets current number of rows in the global cut pool */
int SCIPgetNPoolCuts(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPoolCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPcutpoolGetNCuts(scip->cutpool);
}

/** creates a cut pool */
RETCODE SCIPcreateCutpool(
   SCIP*            scip,               /**< SCIP data structure */
   CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   int              agelimit            /**< maximum age a cut can reach before it is deleted from the pool */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateCutpool", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPcutpoolCreate(cutpool, scip->mem->solvemem, agelimit) );

   return SCIP_OKAY;
}

/** frees a cut pool */
RETCODE SCIPfreeCutpool(
   SCIP*            scip,               /**< SCIP data structure */
   CUTPOOL**        cutpool             /**< pointer to store cut pool */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeCutpool", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPcutpoolFree(cutpool, scip->mem->solvemem, scip->set, scip->lp) );

   return SCIP_OKAY;
}

/** if not already existing, adds row to a cut pool and captures it */
RETCODE SCIPaddRowCutpool(
   SCIP*            scip,               /**< SCIP data structure */
   CUTPOOL*         cutpool,            /**< cut pool */
   ROW*             row                 /**< cutting plane to add */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddRowCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPcutpoolAddRow(cutpool, scip->mem->solvemem, scip->set, row) );

   return SCIP_OKAY;
}

/** adds row to a cut pool and captures it; doesn't check for multiple cuts */
RETCODE SCIPaddNewRowCutpool(
   SCIP*            scip,               /**< SCIP data structure */
   CUTPOOL*         cutpool,            /**< cut pool */
   ROW*             row                 /**< cutting plane to add */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddNewRowCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPcutpoolAddNewRow(cutpool, scip->mem->solvemem, scip->set, row) );

   return SCIP_OKAY;
}

/** removes the LP row from a cut pool */
RETCODE SCIPdelRowCutpool(
   SCIP*            scip,               /**< SCIP data structure */
   CUTPOOL*         cutpool,            /**< cut pool */
   ROW*             row                 /**< row to remove */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPdelRowCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPcutpoolDelRow(cutpool, scip->mem->solvemem, scip->set, scip->stat, scip->lp, row) );
   return SCIP_OKAY;
}

/** separates cuts from a cut pool */
RETCODE SCIPseparateCutpool(
   SCIP*            scip,               /**< SCIP data structure */
   CUTPOOL*         cutpool,            /**< cut pool */
   RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   assert(SCIPtreeGetCurrentNode(scip->tree) != NULL);

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      errorMessage("cannot add cuts, because node LP is not processed\n");
      return SCIP_INVALIDCALL;
   }
   
   CHECK_OKAY( SCIPcutpoolSeparate(cutpool, scip->mem->solvemem, scip->set, scip->stat, scip->lp, scip->sepastore,
         (SCIPtreeGetCurrentDepth(scip->tree) == 0), result) );
   return SCIP_OKAY;
}




/*
 * LP diving methods
 */

/** initiates LP diving, making methods SCIPchgVarObjDive(), SCIPchgVarLbDive(), and SCIPchgVarUbDive() available */
RETCODE SCIPstartDive(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPstartDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPlpDiving(scip->lp) )
   {
      errorMessage("already in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      errorMessage("cannot start diving at a pseudo node\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPlpStartDive(scip->lp, scip->mem->solvemem, scip->set) );

   return SCIP_OKAY;
}

/** quits LP diving and resets bounds and objective values of columns to the current node's values */
RETCODE SCIPendDive(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPendDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      errorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* reset the probably changed LP's cutoff bound */
   CHECK_OKAY( SCIPlpSetCutoffbound(scip->lp, scip->set, scip->primal->cutoffbound) );

   /* unmark the diving flag in the LP and reset all variables' objective and bound values */
   CHECK_OKAY( SCIPlpEndDive(scip->lp, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
         scip->transprob->vars, scip->transprob->nvars) );

   /* if a new best solution was created, the cutoff of the tree was delayed due to diving;
    * the cutoff has to be done now.
    */
   if( scip->tree->cutoffdelayed )
   {
      CHECK_OKAY( SCIPtreeCutoff(scip->tree, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
            scip->primal->cutoffbound) );
   }

   return SCIP_OKAY;
}

/** changes variable's objective value in current dive */
RETCODE SCIPchgVarObjDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the objective value for */
   Real             newobj              /**< new objective value */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarObjDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      errorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPvarChgObjDive(var, scip->set, scip->lp, newobj) );

   /* invalidate the LP's cutoff bound, since this has nothing to do with the current objective value anymore;
    * the cutoff bound is reset in SCIPendDive()
    */
   CHECK_OKAY( SCIPlpSetCutoffbound(scip->lp, scip->set, SCIPsetInfinity(scip->set)) );

   return SCIP_OKAY;
}

/** changes variable's lower bound in current dive */
RETCODE SCIPchgVarLbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarLbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( !SCIPlpDiving(scip->lp) )
   {
      errorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPvarChgLbDive(var, scip->set, scip->lp, newbound) );

   return SCIP_OKAY;
}

/** changes variable's upper bound in current dive */
RETCODE SCIPchgVarUbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarUbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( !SCIPlpDiving(scip->lp) )
   {
      errorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPvarChgUbDive(var, scip->set, scip->lp, newbound) );

   return SCIP_OKAY;
}

/** gets variable's objective value in current dive */
Real SCIPgetVarObjDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get the bound for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarObjDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( !SCIPlpDiving(scip->lp) )
   {
      errorMessage("not in diving mode\n");
      abort();
   }

   return SCIPvarGetObjLP(var);
}

/** gets variable's lower bound in current dive */
Real SCIPgetVarLbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get the bound for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarLbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( !SCIPlpDiving(scip->lp) )
   {
      errorMessage("not in diving mode\n");
      abort();
   }

   return SCIPvarGetLbLP(var);
}

/** gets variable's upper bound in current dive */
Real SCIPgetVarUbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get the bound for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarUbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( !SCIPlpDiving(scip->lp) )
   {
      errorMessage("not in diving mode\n");
      abort();
   }

   return SCIPvarGetUbLP(var);
}

/** solves the LP of the current dive */
RETCODE SCIPsolveDiveLP(
   SCIP*            scip,               /**< SCIP data structure */
   int              itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsolveDiveLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      errorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* solve diving LP */
   CHECK_OKAY( SCIPlpSolveAndEval(scip->lp, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, 
         itlim, FALSE, FALSE, lperror) );

   /* analyze an infeasible LP (not necessary in the root node) */
   if( !scip->set->misc_exactsolve && SCIPtreeGetCurrentDepth(scip->tree) > 0
      && (SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_INFEASIBLE
         || SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_OBJLIMIT) )
   {
      CHECK_OKAY( SCIPconflictAnalyzeLP(scip->conflict, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
            scip->tree, scip->lp, NULL) );
   }

   return SCIP_OKAY;
}

/** returns the number of the node in the current branch and bound run, where the last LP diving was applied */
Longint SCIPgetLastDivenode(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLastDivenode", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return scip->stat->lastdivenode;
}




/*
 * probing methods
 */

/** initiates probing, making methods SCIPchgVarLbProbing(), SCIPchgVarUbProbing(), SCIPfixVarProbing() and
 *  SCIPpropagateProbing() available
 */
RETCODE SCIPstartProbing(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPstartProbing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeProbing(scip->tree) )
   {
      errorMessage("already in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPtreeStartProbing(scip->tree, scip->mem->solvemem, scip->set, scip->lp) );

   return SCIP_OKAY;
}

/** quits probing and resets bounds and constraints to the focus node's environment */
RETCODE SCIPendProbing(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPendProbing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      errorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   /* switch back from probing to normal operation mode and restore variables and constraints to focus node */
   CHECK_OKAY( SCIPtreeEndProbing(scip->tree, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
         scip->branchcand, scip->eventqueue) );

   return SCIP_OKAY;
}

/** injects a change of variable's lower bound into probing node; the same can also be achieved with a call to
 *  SCIPchgVarLb(), but in this case, the bound change would be treated like a deduction instead of a branching decision
 */
RETCODE SCIPchgVarLbProbing(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarLbProbing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( !SCIPtreeProbing(scip->tree) )
   {
      errorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPnodeAddBoundchg(SCIPtreeGetProbingNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
         scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, TRUE) );

   return SCIP_OKAY;
}

/** injects a change of variable's upper bound into probing node; the same can also be achieved with a call to
 *  SCIPchgVarUb(), but in this case, the bound change would be treated like a deduction instead of a branching decision
 */
RETCODE SCIPchgVarUbProbing(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarUbProbing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( !SCIPtreeProbing(scip->tree) )
   {
      errorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPnodeAddBoundchg(SCIPtreeGetProbingNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
         scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, TRUE) );

   return SCIP_OKAY;
}

/** injects a change of variable's bounds into probing node to fix the variable to the specified value; the same can also
 *  be achieved with a call to SCIPfixVar(), but in this case, the bound changes would be treated like deductions instead
 *  of branching decisions
 */
extern
RETCODE SCIPfixVarProbing(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             fixedval            /**< value to fix variable to */
   )
{
   Real lb;
   Real ub;

   CHECK_OKAY( checkStage(scip, "SCIPfixVarProbing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( !SCIPtreeProbing(scip->tree) )
   {
      errorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   if( SCIPsetIsGT(scip->set, fixedval, lb) )
   {
      CHECK_OKAY( SCIPnodeAddBoundchg(SCIPtreeGetProbingNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, fixedval, SCIP_BOUNDTYPE_LOWER, TRUE) );
   }
   if( SCIPsetIsLT(scip->set, fixedval, ub) )
   {
      CHECK_OKAY( SCIPnodeAddBoundchg(SCIPtreeGetProbingNode(scip->tree), scip->mem->solvemem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, fixedval, SCIP_BOUNDTYPE_UPPER, TRUE) );
   }

   return SCIP_OKAY;
}

/** applies domain propagation on the probing sub problem, that was changed after SCIPstartProbing() was called;
 *  the propagated domains of the variables can be accessed with the usual bound accessing calls SCIPvarGetLbLocal()
 *  and SCIPvarGetUbLocal(); the propagation is only valid locally, i.e. the local bounds as well as the changed
 *  bounds due to SCIPchgVarLbProbing(), SCIPchgVarUbProbing(), and SCIPfixVarProbing() are used for propagation
 */
RETCODE SCIPpropagateProbing(
   SCIP*            scip,               /**< SCIP data structure */
   Bool*            cutoff              /**< pointer to store whether the probing node can be cut off */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPpropagateProbing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( !SCIPtreeProbing(scip->tree) )
   {
      errorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPpropagateDomains(scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, cutoff) );

   return SCIP_OKAY;
}




/*
 * branching methods
 */

/** gets branching candidates for LP solution branching (fractional variables) along with solution values,
 *  fractionalities, and number of branching candidates;
 *  branching rules should always select the branching candidate among the first npriolpcands of the candidate
 *  list
 */
RETCODE SCIPgetLPBranchCands(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands,           /**< pointer to store the number of LP branching candidates, or NULL */
   int*             npriolpcands        /**< pointer to store the number of candidates with maximal priority, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_OPTIMAL && SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      errorMessage("LP not solved to optimality\n");
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPbranchcandGetLPCands(scip->branchcand, scip->set, scip->stat, scip->lp,
         lpcands, lpcandssol, lpcandsfrac, nlpcands, npriolpcands) );
   
   return SCIP_OKAY;
}

/** gets number of branching candidates for LP solution branching (number of fractional variables) */
int SCIPgetNLPBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   int nlpcands;

   CHECK_ABORT( checkStage(scip, "SCIPgetNLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      errorMessage("LP not solved to optimality\n");
      abort();
   }

   CHECK_ABORT( SCIPbranchcandGetLPCands(scip->branchcand, scip->set, scip->stat, scip->lp,
         NULL, NULL, NULL, &nlpcands, NULL) );
   
   return nlpcands;
}

/** gets number of branching candidates with maximal priority for LP solution branching */
int SCIPgetNPrioLPBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   int npriolpcands;

   CHECK_ABORT( checkStage(scip, "SCIPgetNPrioLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      errorMessage("LP not solved to optimality\n");
      abort();
   }

   CHECK_ABORT( SCIPbranchcandGetLPCands(scip->branchcand, scip->set, scip->stat, scip->lp,
         NULL, NULL, NULL, NULL, &npriolpcands) );
   
   return npriolpcands;
}

/** gets branching candidates for pseudo solution branching (nonfixed variables) along with the number of candidates */
RETCODE SCIPgetPseudoBranchCands(
   SCIP*            scip,               /**< SCIP data structure */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands,       /**< pointer to store the number of pseudo branching candidates, or NULL */
   int*             npriopseudocands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchcandGetPseudoCands(scip->branchcand, scip->set, scip->transprob, 
         pseudocands, npseudocands, npriopseudocands) );

   return SCIP_OKAY;
}

/** gets branching candidates for pseudo solution branching (nonfixed variables) */
int SCIPgetNPseudoBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPseudoCands(scip->branchcand);
}

/** gets number of branching candidates with maximal branch priority for pseudo solution branching */
int SCIPgetNPrioPseudoBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPrioPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoCands(scip->branchcand);
}

/** gets number of binary branching candidates with maximal branch priority for pseudo solution branching */
int SCIPgetNPrioPseudoBranchBins(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPrioPseudoBranchBins", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoBins(scip->branchcand);
}

/** gets number of integer branching candidates with maximal branch priority for pseudo solution branching */
int SCIPgetNPrioPseudoBranchInts(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPrioPseudoBranchInts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoInts(scip->branchcand);
}

/** gets number of implicit integer branching candidates with maximal branch priority for pseudo solution branching */
int SCIPgetNPrioPseudoBranchImpls(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPrioPseudoBranchImpls", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoImpls(scip->branchcand);
}

/** calculates the branching score out of the gain predictions for a binary branching */
Real SCIPgetBranchScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   Real             downgain,           /**< prediction of objective gain for rounding downwards */
   Real             upgain              /**< prediction of objective gain for rounding upwards */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBranchScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchGetScore(scip->set, var, downgain, upgain);
}

/** calculates the branching score out of the gain predictions for a branching with arbitrary many children */
Real SCIPgetBranchScoreMultiple(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   int              nchildren,          /**< number of children that the branching will create */
   Real*            gains               /**< prediction of objective gain for each child */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBranchScoreMultiple", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchGetScoreMultiple(scip->set, var, nchildren, gains);
}

/** creates a child node of the focus node */
RETCODE SCIPcreateChild(
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           node,               /**< pointer to node data structure */
   Real             nodeselprio         /**< node selection priority of new node */
   )
{
   assert(node != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcreateChild", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPnodeCreateChild(node, scip->mem->solvemem, scip->set, scip->stat, scip->tree, nodeselprio) );
   
   return SCIP_OKAY;
}

/** branches on a variable; if solution value x' is fractional, two child nodes are created
 *  (x <= floor(x'), x >= ceil(x')), if solution value is integral, three child nodes are created
 *  (x <= x'-1, x == x', x >= x'+1)
 */
RETCODE SCIPbranchVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to branch on */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPbranchVar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      errorMessage("cannot branch on continuous variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }
   if( SCIPsetIsEQ(scip->set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      errorMessage("cannot branch on variable <%s> with fixed domain [%g,%g]\n", 
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPtreeBranchVar(scip->tree, scip->mem->solvemem, scip->set, scip->stat, scip->lp, scip->branchcand,
         scip->eventqueue, var) );

   return SCIP_OKAY;
}

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 */
RETCODE SCIPbranchLP(
   SCIP*            scip,               /**< SCIP data structure */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPbranchLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchExecLP(scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, 
         scip->sepastore, scip->branchcand, scip->eventqueue, scip->primal->upperbound, TRUE, result) );

   return SCIP_OKAY;
}

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN */
RETCODE SCIPbranchPseudo(
   SCIP*            scip,               /**< SCIP data structure */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPbranchPseudo", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchExecPseudo(scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, 
         scip->branchcand, scip->eventqueue, scip->primal->upperbound, TRUE, result) );

   return SCIP_OKAY;
}




/*
 * primal solutions
 */

/** creates a primal solution, initialized to zero */
RETCODE SCIPcreateSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolCreate(sol, scip->mem->solvemem, scip->set, scip->stat, scip->primal, scip->tree, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current LP solution */
RETCODE SCIPcreateLPSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      errorMessage("LP solution does not exist\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPsolCreateLPSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->primal, scip->tree, scip->lp,
         heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current pseudo solution */
RETCODE SCIPcreatePseudoSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreatePseudoSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolCreatePseudoSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->primal, scip->tree, scip->lp,
         heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current solution */
RETCODE SCIPcreateCurrentSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolCreateCurrentSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->primal, scip->tree, scip->lp,
         heur) );

   return SCIP_OKAY;
}

/** frees primal CIP solution */
RETCODE SCIPfreeSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol                 /**< pointer to the solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolFree(sol, scip->mem->solvemem, scip->primal) );

   return SCIP_OKAY;
}

/** links a primal solution to the current LP solution */
RETCODE SCIPlinkLPSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPlinkLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      errorMessage("LP solution does not exist\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPsolLinkLPSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** links a primal solution to the current pseudo solution */
RETCODE SCIPlinkPseudoSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPlinkPseudoSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolLinkPseudoSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** links a primal solution to the current LP or pseudo solution */
RETCODE SCIPlinkCurrentSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPlinkCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolLinkCurrentSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** clears a primal solution */
RETCODE SCIPclearSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolClear(sol, scip->stat, scip->tree) );

   return SCIP_OKAY;
}

/** stores solution values of variables in solution's own array */
RETCODE SCIPunlinkSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPunlinkSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolUnlink(sol, scip->set, scip->transprob) );

   return SCIP_OKAY;
}

/** sets value of variable in primal CIP solution */
RETCODE SCIPsetSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetSolVal", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolSetVal(sol, scip->set, scip->stat, scip->tree, var, val) );

   return SCIP_OKAY;
}

/** sets values of multiple variables in primal CIP solution */
RETCODE SCIPsetSolVals(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   int              nvars,              /**< number of variables to set solution value for */
   VAR**            vars,               /**< array with variables to add to solution */
   Real*            vals                /**< array with solution values of variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPsetSolVals", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   for( v = 0; v < nvars; ++v )
   {
      CHECK_OKAY( SCIPsolSetVal(sol, scip->set, scip->stat, scip->tree, vars[v], vals[v]) );
   }

   return SCIP_OKAY;
}


/** increases value of variable in primal CIP solution */
RETCODE SCIPincSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   VAR*             var,                /**< variable to increase solution value for */
   Real             incval              /**< increment for solution value of variable */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPincSolVal", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolIncVal(sol, scip->set, scip->stat, scip->tree, var, incval) );

   return SCIP_OKAY;
}

/** returns value of variable in primal CIP solution, or in current LP/pseudo solution */
Real SCIPgetSolVal(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   VAR*             var                 /**< variable to get value for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolVal", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
      return SCIPsolGetVal(sol, scip->stat, var);
   else
   {
      CHECK_ABORT( checkStage(scip, "SCIPgetSolVal(sol==NULL)", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      return SCIPvarGetSol(var, SCIPtreeHasCurrentNodeLP(scip->tree));
   }
}

/** gets values of multiple variables in primal CIP solution */
RETCODE SCIPgetSolVals(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int              nvars,              /**< number of variables to get solution value for */
   VAR**            vars,               /**< array with variables to get value for */
   Real*            vals                /**< array to store solution values of variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetSolVals", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
   {
      for( v = 0; v < nvars; ++v )
         vals[v] = SCIPsolGetVal(sol, scip->stat, vars[v]);
   }
   else
   {
      CHECK_OKAY( SCIPgetVarSols(scip, nvars, vars, vals) );
   }

   return SCIP_OKAY;
}

/** returns objective value of primal CIP solution w.r.t. original problem, or current LP/pseudo objective value */
Real SCIPgetSolOrigObj(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolOrigObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
      return SCIPprobExternObjval(scip->origprob, scip->set, 
         SCIPprobExternObjval(scip->transprob, scip->set, SCIPsolGetObj(sol)));
   else
   {
      CHECK_ABORT( checkStage(scip, "SCIPgetSolOrigObj(sol==NULL)", 
            FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      if( SCIPtreeHasCurrentNodeLP(scip->tree) )
         return SCIPprobExternObjval(scip->origprob, scip->set,
            SCIPprobExternObjval(scip->transprob, scip->set, SCIPlpGetObjval(scip->lp, scip->set)));
      else
         return SCIPprobExternObjval(scip->origprob, scip->set, 
            SCIPprobExternObjval(scip->transprob, scip->set, SCIPlpGetPseudoObjval(scip->lp, scip->set)));
   }
}

/** returns transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value */
Real SCIPgetSolTransObj(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolTransObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
      return SCIPsolGetObj(sol);
   else
   {
      CHECK_ABORT( checkStage(scip, "SCIPgetSolTransObj(sol==NULL)", 
            FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      if( SCIPtreeHasCurrentNodeLP(scip->tree) )
         return SCIPlpGetObjval(scip->lp, scip->set);
      else
         return SCIPlpGetPseudoObjval(scip->lp, scip->set);
   }
}

/** maps original space objective value into transformed objective value */
Real SCIPtransformObj(
   SCIP*            scip,               /**< SCIP data structure */
   Real             obj                 /**< original space objective value to transform */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPtransformObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPprobInternObjval(scip->transprob, scip->set, SCIPprobInternObjval(scip->origprob, scip->set, obj));
}

/** maps transformed objective value into original space */
Real SCIPretransformObj(
   SCIP*            scip,               /**< SCIP data structure */
   Real             obj                 /**< transformed objective value to retransform in original space */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPretransformObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->origprob, scip->set, SCIPprobExternObjval(scip->transprob, scip->set, obj));
}

/** gets clock time, when this solution was found */
Real SCIPgetSolTime(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolTime", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetTime(sol);
}

/** gets branch and bound run number, where this solution was found */
int SCIPgetSolRunnum(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolRunnum", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetRunnum(sol);
}

/** gets node number of the specific branch and bound run, where this solution was found */
Longint SCIPgetSolNodenum(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolNodenum", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetNodenum(sol);
}

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
HEUR* SCIPgetSolHeur(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolHeur", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetHeur(sol);
}

/** outputs non-zero variables of solution in original problem space to file stream */
RETCODE SCIPprintSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( file == NULL )
      file = stdout;

   fprintf(file, "objective value:                 ");
   SCIPprintReal(scip, SCIPprobExternObjval(scip->origprob, scip->set, 
         SCIPprobExternObjval(scip->transprob, scip->set, SCIPsolGetObj(sol))), file);
   fprintf(file, "\n");

   CHECK_OKAY( SCIPsolPrint(sol, scip->set, scip->stat, scip->origprob, scip->transprob, file) );

   return SCIP_OKAY;
}

/** outputs non-zero variables of solution in transformed problem space to file stream */
RETCODE SCIPprintTransSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintSolTrans", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( file == NULL )
      file = stdout;

   fprintf(file, "objective value:                 ");
   SCIPprintReal(scip, SCIPsolGetObj(sol), file);
   fprintf(file, "\n");

   CHECK_OKAY( SCIPsolPrint(sol, scip->set, scip->stat, scip->transprob, NULL, file) );
   
   return SCIP_OKAY;
}

/** gets number of feasible primal solutions stored in the solution storage */
Longint SCIPgetNSols(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNSols", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->nsols;
}

/** gets array of feasible primal solutions stored in the solution storage */
SOL** SCIPgetSols(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSols", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->sols;
}

/** gets best feasible primal solution found so far, or NULL if no solution has been found */
SOL* SCIPgetBestSol(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( scip->primal->nsols > 0 ) 
   {
      assert(scip->primal->sols != NULL);
      assert(scip->primal->sols[0] != NULL);
      return scip->primal->sols[0];
   }
   return NULL;
}

/** outputs best feasible primal solution found so far to file stream */
RETCODE SCIPprintBestSol(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   SOL* sol;

   CHECK_OKAY( checkStage(scip, "SCIPprintBestSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   sol = SCIPgetBestSol(scip);

   if( file == NULL )
      file = stdout;

   if( sol == NULL )
      fprintf(file, "no solution available\n");
   else
   {
      CHECK_OKAY( SCIPprintSol(scip, sol, file) );
   }

   return SCIP_OKAY;
}

/** outputs best feasible primal solution found so far in transformed variables to file stream */
RETCODE SCIPprintBestTransSol(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   SOL* sol;

   CHECK_OKAY( checkStage(scip, "SCIPprintBestTransSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   sol = SCIPgetBestSol(scip);

   if( file == NULL )
      file = stdout;

   if( sol == NULL )
      fprintf(file, "no solution available\n");
   else
   {
      CHECK_OKAY( SCIPprintTransSol(scip, sol, file) );
   }

   return SCIP_OKAY;
}

/** try to round given solution */
RETCODE SCIProundSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   Bool*            success             /**< pointer to store whether rounding was successful */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIProundSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolRound(sol, scip->set, scip->stat, scip->transprob, scip->tree, success) );

   return SCIP_OKAY;
}

/** adds feasible primal solution to solution storage by copying it */
RETCODE SCIPaddSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal CIP solution */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalAddSol(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
         scip->lp, scip->eventfilter, sol, stored) );

   return SCIP_OKAY;
}

/** adds primal solution to solution storage, frees the solution afterwards */
RETCODE SCIPaddSolFree(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddSolFree", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalAddSolFree(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
         scip->lp, scip->eventfilter, sol, stored) );

   return SCIP_OKAY;
}

/** adds current LP/pseudo solution to solution storage */
RETCODE SCIPaddCurrentSol(
   SCIP*            scip,               /**< SCIP data structure */
   HEUR*            heur,               /**< heuristic that found the solution */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalAddCurrentSol(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, 
         scip->tree, scip->lp, scip->eventfilter, heur, stored) );

   return SCIP_OKAY;
}

/** checks solution for feasibility; if possible, adds it to storage by copying */
RETCODE SCIPtrySol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPtrySol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalTrySol(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
         scip->lp, scip->eventfilter, sol, checkintegrality, checklprows, stored) );

   return SCIP_OKAY;
}

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards */
RETCODE SCIPtrySolFree(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPtrySolFree", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalTrySolFree(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
         scip->lp, scip->eventfilter, sol, checkintegrality, checklprows, stored) );
   
   return SCIP_OKAY;
}

/** checks current LP/pseudo solution for feasibility; if possible, adds it to storage */
RETCODE SCIPtryCurrentSol(
   SCIP*            scip,               /**< SCIP data structure */
   HEUR*            heur,               /**< heuristic that found the solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPtryCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalTryCurrentSol(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, 
         scip->tree, scip->lp, scip->eventfilter, heur, checkintegrality, checklprows, stored) );

   return SCIP_OKAY;
}

/** checks solution for feasibility without adding it to the solution store */
RETCODE SCIPcheckSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            feasible            /**< stores whether given solution is feasible */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcheckSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   /* if we want to solve exactly, the constraint handlers cannot rely on the LP's feasibility */
   checklprows = checklprows || scip->set->misc_exactsolve;
   
   CHECK_OKAY( SCIPsolCheck(sol, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, 
         checkintegrality, checklprows, feasible) );

   return SCIP_OKAY;
}

/** checks solution for feasibility in original problem without adding it to the solution store;
 *  this method is used to double check a solution in order to validate the presolving process
 */
RETCODE SCIPcheckSolOrig(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal CIP solution */
   Bool*            feasible,           /**< stores whether given solution is feasible */
   CONSHDLR**       infeasconshdlr,     /**< pointer to store constraint handler of first infeasible constraint,
                                         *   or NULL if not needed */
   CONS**           infeascons          /**< pointer to store first infeasible constraint, or NULL if not needed;
                                         *   stores NULL, if a constraint handler that doesn't need constraints
                                         *   rejected the solution */
   )
{
   RESULT result;
   int c;
   int h;

   assert(feasible != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcheckSolOrig", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   if( infeasconshdlr != NULL )
      *infeasconshdlr = NULL;
   if( infeascons != NULL )
      *infeascons = NULL;

   /* check original constraints */
   *feasible = TRUE;
   for( c = 0; c < scip->origprob->nconss; ++c )
   {
      CHECK_OKAY( SCIPconsCheck(scip->origprob->conss[c], scip->set, sol, TRUE, TRUE, &result) );
      if( result != SCIP_FEASIBLE )
      {
         *feasible = FALSE;
         if( infeasconshdlr != NULL )
            *infeasconshdlr = SCIPconsGetHdlr(scip->origprob->conss[c]);
         if( infeascons != NULL )
            *infeascons = scip->origprob->conss[c];
         return SCIP_OKAY;
      }
   }

   /* call constraint handlers that don't need constraints */
   for( h = 0; h < scip->set->nconshdlrs; ++h )
   {
      if( !SCIPconshdlrNeedsCons(scip->set->conshdlrs[h]) )
      {
         CHECK_OKAY( SCIPconshdlrCheck(scip->set->conshdlrs[h], scip->mem->solvemem, scip->set, scip->stat, 
               scip->origprob, sol, TRUE, TRUE, &result) );
         if( result != SCIP_FEASIBLE )
         {
            *feasible = FALSE;
            if( infeasconshdlr != NULL )
               *infeasconshdlr = scip->set->conshdlrs[h];
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}




/*
 * event methods
 */

/** catches a global (not variable dependent) event */
RETCODE SCIPcatchEvent(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*             filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcatchEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPeventfilterAdd(scip->eventfilter, scip->mem->solvemem, scip->set, 
         eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** drops a global event (stops to track event) */
RETCODE SCIPdropEvent(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int              filterpos           /**< position of event filter entry returned by SCIPcatchEvent(), or -1 */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPdropEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPeventfilterDel(scip->eventfilter, scip->mem->solvemem, scip->set, 
         eventtype, eventhdlr, eventdata, filterpos) );
   
   return SCIP_OKAY;
}

/** catches an objective value or domain change event on the given variable */
RETCODE SCIPcatchVarEvent(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to catch event for */
   EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*             filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcatchVarEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( !SCIPvarIsTransformed(var) )
   {
      errorMessage("cannot catch events for original variables\n");
      return SCIP_INVALIDDATA;
   }

   if( (eventtype & SCIP_EVENTTYPE_VARCHANGED) == 0 )
   {
      errorMessage("event is neither an objective value nor a domain change event\n");
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPvarCatchEvent(var, scip->mem->solvemem, scip->set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** drops an objective value or domain change event (stops to track event) on the given variable */
RETCODE SCIPdropVarEvent(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to drop event for */
   EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int              filterpos           /**< position of event filter entry returned by SCIPcatchVarEvent(), or -1 */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPdropVarEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( !SCIPvarIsTransformed(var) )
   {
      errorMessage("cannot drop events for original variables\n");
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPvarDropEvent(var, scip->mem->solvemem, scip->set, eventtype, eventhdlr, eventdata, filterpos) );
   
   return SCIP_OKAY;
}




/*
 * tree methods
 */

/** gets children of focus node along with the number of children */
RETCODE SCIPgetChildren(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          children,           /**< pointer to store children array, or NULL if not needed */
   int*             nchildren           /**< pointer to store number of children, or NULL if not needed */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetChildren", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( children != NULL )
      *children = scip->tree->children;
   if( nchildren != NULL )
      *nchildren = scip->tree->nchildren;
   
   return SCIP_OKAY;
}

/** gets number of children of focus node */
int SCIPgetNChildren(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNChildren", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->tree->nchildren;
}

/** gets siblings of focus node along with the number of siblings */
RETCODE SCIPgetSiblings(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          siblings,           /**< pointer to store siblings array, or NULL if not needed */
   int*             nsiblings           /**< pointer to store number of siblings, or NULL if not needed */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetSiblings", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( siblings != NULL )
      *siblings = scip->tree->siblings;
   if( nsiblings != NULL )
      *nsiblings = scip->tree->nsiblings;
   
   return SCIP_OKAY;
}

/** gets number of siblings of focus node */
int SCIPgetNSiblings(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNSiblings", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->tree->nsiblings;
}

/** gets leaves of the tree along with the number of leaves */
RETCODE SCIPgetLeaves(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          leaves,             /**< pointer to store leaves array, or NULL if not needed */
   int*             nleaves             /**< pointer to store number of leaves, or NULL if not needed */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLeaves", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( leaves != NULL )
      *leaves = SCIPnodepqNodes(scip->tree->leaves);
   if( nleaves != NULL )
      *nleaves = SCIPnodepqLen(scip->tree->leaves);
   
   return SCIP_OKAY;
}

/** gets number of leaves in the tree */
int SCIPgetNLeaves(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNLeaves", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPnodepqLen(scip->tree->leaves);
}

/** gets the best child of the focus node w.r.t. the node selection priority assigned by the branching rule */
NODE* SCIPgetPrioChild(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPrioChild", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetPrioChild(scip->tree);
}

/** gets the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule */
NODE* SCIPgetPrioSibling(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPrioSibling", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetPrioSibling(scip->tree);
}

/** gets the best child of the focus node w.r.t. the node selection strategy */
NODE* SCIPgetBestChild(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestChild", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestChild(scip->tree, scip->set);
}

/** gets the best sibling of the focus node w.r.t. the node selection strategy */
NODE* SCIPgetBestSibling(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestSibling", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestSibling(scip->tree, scip->set);
}

/** gets the best leaf from the node queue w.r.t. the node selection strategy */
NODE* SCIPgetBestLeaf(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestLeaf", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestLeaf(scip->tree);
}

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy */
NODE* SCIPgetBestNode(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestNode(scip->tree, scip->set);
}

/** gets the node with smallest lower bound from the tree (child, sibling, or leaf) */
NODE* SCIPgetBestboundNode(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestboundNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetLowerboundNode(scip->tree, scip->set);
}

/** cuts off node and whole sub tree from branch and bound tree */
RETCODE SCIPcutoffNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node                /**< node that should be cut off */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcutoffNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodeCutoff(node, scip->set, scip->stat, scip->tree);

   return SCIP_OKAY;
}




/*
 * statistic methods
 */

/** gets number of branch and bound runs performed, including the current run */
int SCIPgetNRuns(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNRuns", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->stat->nruns;
}

/** gets number of processed nodes in current run, including the focus node */
Longint SCIPgetNNodes(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNNodes", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->stat->nnodes;
}

/** gets total number of processed nodes in all runs, including the focus node */
Longint SCIPgetNTotalNodes(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNTotalNodes", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->stat->ntotalnodes;
}

/** gets number of nodes left in the tree (children + siblings + leaves) */
int SCIPgetNNodesLeft(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNNodesLeft", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetNNodes(scip->tree);
}

/** gets total number of LPs solved so far */
int SCIPgetNLPs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nlps;
}

/** gets total number of simplex iterations used so far in primal and dual simplex */
Longint SCIPgetNLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nlpiterations;
}

/** gets total number of LPs solved so far that were resolved from an advanced start basis */
int SCIPgetNResolveLPs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNResolveLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nresolvelps;
}

/** gets total number of simplex iterations used so far in primal and dual simplex calls where an advanced start basis
 *  was available
 */
Longint SCIPgetNResolveLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNResolveLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nresolvelpiterations;
}

/** gets total number of LPs solved so far for node relaxations */
int SCIPgetNNodeLPs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNNodeLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nnodelps;
}

/** gets total number of simplex iterations used so far for node relaxations */
Longint SCIPgetNNodeLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNNodeLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nnodelpiterations;
}

/** gets total number of LPs solved so far for initial LP in node relaxations */
int SCIPgetNNodeInitLPs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNInitNodeLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->ninitlps;
}

/** gets total number of simplex iterations used so far for initial LP in node relaxations */
Longint SCIPgetNNodeInitLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNNodeInitLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->ninitlpiterations;
}

/** gets total number of LPs solved so far during diving */
int SCIPgetNDivingLPs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNDivingLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->ndivinglps;
}

/** gets total number of simplex iterations used so far during diving */
Longint SCIPgetNDivingLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNDivingLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->ndivinglpiterations;
}

/** gets total number of times, strong branching was called (each call represents solving two LPs) */
int SCIPgetNStrongbranchs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNStrongbranchs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nstrongbranchs;
}

/** gets total number of simplex iterations used so far in strong branching */
Longint SCIPgetNStrongbranchLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNStrongbranchLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nsblpiterations;
}

/** gets total number of times, strong branching was called at the root node (each call represents solving two LPs) */
int SCIPgetNRootStrongbranchs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNRootStrongbranchs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nrootstrongbranchs;
}

/** gets total number of simplex iterations used so far in strong branching at the root node */
Longint SCIPgetNRootStrongbranchLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNRootStrongbranchLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nrootsblpiterations;
}

/** gets number of pricing rounds performed so far at the current node */
int SCIPgetNPriceRounds(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPriceRounds", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->npricerounds;
}

/** get current number of variables in the pricing store */
int SCIPgetNPricevars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPricevars", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPpricestoreGetNVars(scip->pricestore);
}

/** get total number of pricing variables found so far */
int SCIPgetNPricevarsFound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPricevarsFound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPpricestoreGetNVarsFound(scip->pricestore);
}

/** get total number of pricing variables applied to the LPs */
int SCIPgetNPricevarsApplied(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPricevarsApplied", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPpricestoreGetNVarsApplied(scip->pricestore);
}

/** gets number of separation rounds performed so far at the current node */
int SCIPgetNSepaRounds(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNSepaRounds", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nseparounds;
}

/** get current number of cuts in the cut store */
int SCIPgetNCuts(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNCuts", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPsepastoreGetNCuts(scip->sepastore);
}

/** get total number of cuts found so far */
int SCIPgetNCutsFound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNCutsFound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPsepastoreGetNCutsFound(scip->sepastore);
}

/** get number of cuts found so far in current separation round */
int SCIPgetNCutsFoundRound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNCutsFoundRound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPsepastoreGetNCutsFoundRound(scip->sepastore);
}

/** get total number of cuts applied to the LPs */
int SCIPgetNCutsApplied(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNCutsApplied", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPsepastoreGetNCutsApplied(scip->sepastore);
}

/** get total number of clauses found in conflict analysis (conflict and reconvergence clauses) */
Longint SCIPgetNConflictClausesFound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNConflictClausesFound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPconflictGetNPropConflictClauses(scip->conflict)
      + SCIPconflictGetNPropReconvergenceClauses(scip->conflict)
      + SCIPconflictGetNLPConflictClauses(scip->conflict)
      + SCIPconflictGetNLPReconvergenceClauses(scip->conflict)
      + SCIPconflictGetNStrongbranchConflictClauses(scip->conflict)
      + SCIPconflictGetNStrongbranchReconvergenceClauses(scip->conflict)
      + SCIPconflictGetNPseudoConflictClauses(scip->conflict)
      + SCIPconflictGetNPseudoReconvergenceClauses(scip->conflict);
}

/** gets depth of current node, or -1 if no current node exists */
int SCIPgetDepth(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetDepth", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetCurrentDepth(scip->tree);
}

/** gets maximal depth of all processed nodes in current branch and bound run */
int SCIPgetMaxDepth(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetMaxDepth", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->maxdepth;
}

/** gets maximal depth of all processed nodes over all branch and bound runs */
int SCIPgetMaxTotalDepth(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetMaxTotalDepth", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->maxtotaldepth;
}

/** gets total number of backtracks, i.e. number of times, the new node was selected from the leaves queue */
Longint SCIPgetNBacktracks(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNBacktracks", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nbacktracks;
}

/** gets current plunging depth (succ. times, a child was selected as next node) */
int SCIPgetPlungeDepth(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPlungeDepth", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->plungedepth;
}

/** gets total number of active constraints at the current node */
int SCIPgetNActiveConss(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNActiveConss", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nactiveconss;
}

/** gets total number of enabled constraints at the current node */
int SCIPgetNEnabledConss(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNEnabledConss", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nenabledconss;
}

/** gets total number of globally valid constraints currently in the problem */
int SCIPgetNGlobalConss(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNGlobalConss", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nconss;

   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nconss;

   default:
      errorMessage("invalid SCIP stage\n");
      abort();
   }
}

/** gets average dual bound of all unprocessed nodes */
Real SCIPgetAvgDualbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetAvgDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->origprob, scip->set, 
      SCIPprobExternObjval(scip->transprob, scip->set, SCIPtreeGetAvgLowerbound(scip->tree, scip->primal->cutoffbound)));
}

/** gets average lower (dual) bound of all unprocessed nodes in transformed problem */
Real SCIPgetAvgLowerbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetAvgLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetAvgLowerbound(scip->tree, scip->primal->cutoffbound);
}

/** gets global dual bound */
Real SCIPgetDualbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   Real lowerbound;

   CHECK_ABORT( checkStage(scip, "SCIPgetDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   lowerbound = SCIPtreeGetLowerbound(scip->tree, scip->set);
   if( SCIPsetIsInfinity(scip->set, lowerbound) )
      return SCIPgetPrimalbound(scip);
   else
      return SCIPprobExternObjval(scip->origprob, scip->set, SCIPprobExternObjval(scip->transprob, scip->set, lowerbound));
}

/** gets global lower (dual) bound in transformed problem */
Real SCIPgetLowerbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetLowerbound(scip->tree, scip->set);
}

/** gets dual bound of the root node */
Real SCIPgetDualboundRoot(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetDualboundRoot", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( SCIPsetIsInfinity(scip->set, scip->stat->rootlowerbound) )
      return SCIPgetPrimalbound(scip);
   else
      return SCIPprobExternObjval(scip->origprob, scip->set, 
         SCIPprobExternObjval(scip->transprob, scip->set, scip->stat->rootlowerbound));
}

/** gets lower (dual) bound in transformed problem of the root node */
Real SCIPgetLowerboundRoot(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLowerboundRoot", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPnodeGetLowerbound(scip->tree->root);
}

/** gets global primal bound (objective value of best solution or user objective limit) */
Real SCIPgetPrimalbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPrimalbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPprobExternObjval(scip->origprob, scip->set, 
      SCIPprobExternObjval(scip->transprob, scip->set, scip->primal->upperbound));
}

/** gets global upper (primal) bound in transformed problem (objective value of best solution or user objective limit) */
Real SCIPgetUpperbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetUpperbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->upperbound;
}

/** gets global cutoff bound in transformed problem: a sub problem with lower bound larger than the cutoff
 *  cannot contain a better feasible solution; usually, this bound is equal to the upper bound, but if the
 *  objective value is always integral, the cutoff bound is (nearly) one less than the upper bound
 */
Real SCIPgetCutoffbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetCutoffbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->cutoffbound;
}

/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 */
Bool SCIPisPrimalboundSol(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPisPrimalboundSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPprimalUpperboundIsSol(scip->primal);
}

/** gets current gap |(primalbound - dualbound)/dualbound| if both bounds have same sign, or infinity, if they have
 *  opposite sign
 */
Real SCIPgetGap(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   Real primalbound;
   Real dualbound;

   CHECK_ABORT( checkStage(scip, "SCIPgetGap", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( SCIPsetIsInfinity(scip->set, SCIPgetLowerbound(scip)) )
      return 0.0;

   primalbound = SCIPgetPrimalbound(scip);
   dualbound = SCIPgetDualbound(scip);

   if( SCIPsetIsEQ(scip->set, primalbound, dualbound) )
      return 0.0;
   else if( SCIPsetIsZero(scip->set, dualbound)
      || SCIPsetIsInfinity(scip->set, REALABS(primalbound))
      || primalbound * dualbound < 0.0 )
      return SCIPsetInfinity(scip->set);
   else
      return REALABS((primalbound - dualbound)/dualbound);
}

/** gets current gap |(upperbound - lowerbound)/lowerbound| in transformed problem if both bounds have same sign, 
 *  or infinity, if they have opposite sign
 */
Real SCIPgetTransGap(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   Real upperbound;
   Real lowerbound;

   CHECK_ABORT( checkStage(scip, "SCIPgetTransGap", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   upperbound = SCIPgetUpperbound(scip);
   lowerbound = SCIPgetLowerbound(scip);

   if( SCIPsetIsInfinity(scip->set, lowerbound) )
      return 0.0;
   else if( SCIPsetIsEQ(scip->set, upperbound, lowerbound) )
      return 0.0;
   else if( SCIPsetIsZero(scip->set, lowerbound)
      || SCIPsetIsInfinity(scip->set, upperbound)
      || lowerbound * upperbound < 0.0 )
      return SCIPsetInfinity(scip->set);
   else
      return REALABS((upperbound - lowerbound)/lowerbound);
}

/** gets number of feasible primal solutions found so far */
Longint SCIPgetNSolsFound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNSolsFound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->nsolsfound;
}

/** gets number of feasible primal solutions found so far, that improved the primal bound at the time they were found */
Longint SCIPgetNBestSolsFound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNBestSolsFound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->nbestsolsfound;
}

/** outputs original problem to file stream */
RETCODE SCIPprintOrigProblem(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintOrigProblem", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( file == NULL )
      file = stdout;

   if( scip->origprob == NULL )
      fprintf(file, "no problem available\n");
   else
   {
      CHECK_OKAY( SCIPprobPrint(scip->origprob, scip->set, file) );
   }

   return SCIP_OKAY;
}

/** outputs transformed problem to file stream */
RETCODE SCIPprintTransProblem(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintTransProblem", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( file == NULL )
      file = stdout;

   if( scip->transprob == NULL )
      fprintf(file, "no transformed problem available\n");
   else
   {
      CHECK_OKAY( SCIPprobPrint(scip->transprob, scip->set, file) );
   }

   return SCIP_OKAY;
}

static
void printPresolverStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
#if 0
   int nfixedvars;
   int naggrvars;
   int nchgbds;
   int nholes;
   int ndelconss;
   int nchgsides;
   int nchgcoefs;
#endif
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(file != NULL);

   fprintf(file, "Presolvers         :       Time  FixedVars   AggrVars  ChgBounds   AddHoles    DelCons   ChgSides   ChgCoefs\n");

#if 0
   nfixedvars = 0;
   naggrvars = 0;
   nchgbds = 0;
   nholes = 0;
   ndelconss = 0;
   nchgsides = 0;
   nchgcoefs = 0;
#endif

   /* presolver statistics */
   for( i = 0; i < scip->set->npresols; ++i )
   {
      PRESOL* presol;
      
      presol = scip->set->presols[i];
      fprintf(file, "  %-17.17s:", SCIPpresolGetName(presol));
      fprintf(file, " %10.2f %10d %10d %10d %10d %10d %10d %10d\n",
         SCIPpresolGetTime(presol),
         SCIPpresolGetNFixedVars(presol),
         SCIPpresolGetNAggrVars(presol),
         SCIPpresolGetNChgBds(presol),
         SCIPpresolGetNAddHoles(presol),
         SCIPpresolGetNDelConss(presol),
         SCIPpresolGetNChgSides(presol),
         SCIPpresolGetNChgCoefs(presol));
#if 0
      nfixedvars += SCIPpresolGetNFixedVars(presol);
      naggrvars += SCIPpresolGetNAggrVars(presol);
      nchgbds += SCIPpresolGetNChgBds(presol);
      nholes += SCIPpresolGetNAddHoles(presol);
      ndelconss += SCIPpresolGetNDelConss(presol);
      nchgsides += SCIPpresolGetNChgSides(presol);
      nchgcoefs += SCIPpresolGetNChgCoefs(presol);
#endif
   }

   /* constraint handler presolving methods statistics */
   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      CONSHDLR* conshdlr;
      int maxnconss;
      
      conshdlr = scip->set->conshdlrs[i];
      maxnconss = SCIPconshdlrGetMaxNConss(conshdlr);
      if( SCIPconshdlrDoesPresolve(conshdlr)
         && (maxnconss > 0 || !SCIPconshdlrNeedsCons(conshdlr)
            || SCIPconshdlrGetNFixedVars(conshdlr) > 0
            || SCIPconshdlrGetNAggrVars(conshdlr) > 0
            || SCIPconshdlrGetNChgBds(conshdlr) > 0
            || SCIPconshdlrGetNAddHoles(conshdlr) > 0
            || SCIPconshdlrGetNDelConss(conshdlr) > 0
            || SCIPconshdlrGetNChgSides(conshdlr) > 0
            || SCIPconshdlrGetNChgCoefs(conshdlr) > 0) )
      {
         fprintf(file, "  %-17.17s:", SCIPconshdlrGetName(conshdlr));
         fprintf(file, " %10.2f %10d %10d %10d %10d %10d %10d %10d\n",
            SCIPconshdlrGetPresolTime(conshdlr),
            SCIPconshdlrGetNFixedVars(conshdlr),
            SCIPconshdlrGetNAggrVars(conshdlr),
            SCIPconshdlrGetNChgBds(conshdlr),
            SCIPconshdlrGetNAddHoles(conshdlr),
            SCIPconshdlrGetNDelConss(conshdlr),
            SCIPconshdlrGetNChgSides(conshdlr),
            SCIPconshdlrGetNChgCoefs(conshdlr));
#if 0         
         nfixedvars += SCIPconshdlrGetNFixedVars(conshdlr);
         naggrvars += SCIPconshdlrGetNAggrVars(conshdlr);
         nchgbds += SCIPconshdlrGetNChgBds(conshdlr);
         nholes += SCIPconshdlrGetNAddHoles(conshdlr);
         ndelconss += SCIPconshdlrGetNDelConss(conshdlr);
         nchgsides += SCIPconshdlrGetNChgSides(conshdlr);
         nchgcoefs += SCIPconshdlrGetNChgCoefs(conshdlr);
#endif
      }
   }

   /* root node bound changes */
   fprintf(file, "  root node        :          -          -          - %10lld          -          -          -          -\n",
      scip->stat->nrootboundchgs);

#if 0
   /* print total */
   fprintf(file, "  total            :");
   fprintf(file, " %10.2f %10d %10d %10d %10d %10d %10d %10d\n",
      SCIPclockGetTime(scip->stat->presolvingtime),
      nfixedvars,
      naggrvars,
      nchgbds,
      nholes,
      ndelconss,
      nchgsides,
      nchgcoefs);
#endif
}

static
void printConstraintStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(file != NULL);

   /**@todo add constraint statistics: how many constraints (instead of cuts) have been added? */
   fprintf(file, "Constraints        :     Number  #Separate #Propagate    #EnfoLP    #EnfoPS    Cutoffs    DomReds       Cuts      Conss   Children\n");

   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      CONSHDLR* conshdlr;
      int startnconss;
      int maxnconss;

      conshdlr = scip->set->conshdlrs[i];
      startnconss = SCIPconshdlrGetStartNConss(conshdlr);
      maxnconss = SCIPconshdlrGetMaxNConss(conshdlr);
      if( maxnconss > 0 || !SCIPconshdlrNeedsCons(conshdlr) )
      {
         fprintf(file, "  %-17.17s:", SCIPconshdlrGetName(conshdlr));
         fprintf(file, " %10d%c%10lld %10lld %10lld %10lld %10lld %10lld %10lld %10lld %10lld\n",
            startnconss,
            maxnconss > startnconss ? '+' : ' ',
            SCIPconshdlrGetNSepaCalls(conshdlr), 
            SCIPconshdlrGetNPropCalls(conshdlr), 
            SCIPconshdlrGetNEnfoLPCalls(conshdlr),
            SCIPconshdlrGetNEnfoPSCalls(conshdlr),
            SCIPconshdlrGetNCutoffs(conshdlr), 
            SCIPconshdlrGetNDomredsFound(conshdlr), 
            SCIPconshdlrGetNCutsFound(conshdlr), 
            SCIPconshdlrGetNConssFound(conshdlr), 
            SCIPconshdlrGetNChildren(conshdlr));
      }
   }
}

static
void printConstraintTimingStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(file != NULL);

   fprintf(file, "Constraint Timings :  TotalTime   Separate  Propagate     EnfoLP     EnfoPS\n");

   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      CONSHDLR* conshdlr;
      int maxnconss;
      
      conshdlr = scip->set->conshdlrs[i];
      maxnconss = SCIPconshdlrGetMaxNConss(conshdlr);
      if( maxnconss > 0 || !SCIPconshdlrNeedsCons(conshdlr) )
      {
         fprintf(file, "  %-17.17s:", SCIPconshdlrGetName(conshdlr));
         fprintf(file, " %10.2f %10.2f %10.2f %10.2f %10.2f\n",
            SCIPconshdlrGetSepaTime(conshdlr) + SCIPconshdlrGetPropTime(conshdlr) + SCIPconshdlrGetEnfoLPTime(conshdlr)
            + SCIPconshdlrGetEnfoPSTime(conshdlr),
            SCIPconshdlrGetSepaTime(conshdlr), 
            SCIPconshdlrGetPropTime(conshdlr), 
            SCIPconshdlrGetEnfoLPTime(conshdlr), 
            SCIPconshdlrGetEnfoPSTime(conshdlr));
      }
   }
}

static
void printPropagatorStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(file != NULL);

   fprintf(file, "Propagators        :       Time      Calls    Cutoffs    DomReds\n");

   fprintf(file, "  reduced cost str.: %10.2f %10lld          - %10lld\n",
      SCIPclockGetTime(scip->stat->redcoststrtime),
      scip->stat->nredcoststrcalls,
      scip->stat->nredcoststrfound);

   for( i = 0; i < scip->set->nprops; ++i )
      fprintf(file, "  %-17.17s: %10.2f %10lld %10lld %10lld\n",
         SCIPpropGetName(scip->set->props[i]),
         SCIPpropGetTime(scip->set->props[i]),
         SCIPpropGetNCalls(scip->set->props[i]),
         SCIPpropGetNCutoffs(scip->set->props[i]),
         SCIPpropGetNDomredsFound(scip->set->props[i]));
}

static
void printConflictStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   fprintf(file, "Conflict Analysis  :       Time      Calls  Conflicts   Literals    Reconvs ReconvLits   LP Iters\n");
   fprintf(file, "  propagation      : %10.2f %10lld %10lld %10.1f %10lld %10.1f          -\n",
      SCIPconflictGetPropTime(scip->conflict),
      SCIPconflictGetNPropCalls(scip->conflict),
      SCIPconflictGetNPropConflictClauses(scip->conflict),
      SCIPconflictGetNPropConflictClauses(scip->conflict) > 0
      ? (Real)SCIPconflictGetNPropConflictLiterals(scip->conflict)
      / (Real)SCIPconflictGetNPropConflictClauses(scip->conflict) : 0,
      SCIPconflictGetNPropReconvergenceClauses(scip->conflict),
      SCIPconflictGetNPropReconvergenceClauses(scip->conflict) > 0
      ? (Real)SCIPconflictGetNPropReconvergenceLiterals(scip->conflict)
      / (Real)SCIPconflictGetNPropReconvergenceClauses(scip->conflict) : 0);
   fprintf(file, "  infeasible LP    : %10.2f %10lld %10lld %10.1f %10lld %10.1f %10lld\n",
      SCIPconflictGetLPTime(scip->conflict),
      SCIPconflictGetNLPCalls(scip->conflict),
      SCIPconflictGetNLPConflictClauses(scip->conflict),
      SCIPconflictGetNLPConflictClauses(scip->conflict) > 0
      ? (Real)SCIPconflictGetNLPConflictLiterals(scip->conflict)
      / (Real)SCIPconflictGetNLPConflictClauses(scip->conflict) : 0,
      SCIPconflictGetNLPReconvergenceClauses(scip->conflict),
      SCIPconflictGetNLPReconvergenceClauses(scip->conflict) > 0
      ? (Real)SCIPconflictGetNLPReconvergenceLiterals(scip->conflict)
      / (Real)SCIPconflictGetNLPReconvergenceClauses(scip->conflict) : 0,
      SCIPconflictGetNLPIterations(scip->conflict));
   fprintf(file, "  strong branching : %10.2f %10lld %10lld %10.1f %10lld %10.1f %10lld\n",
      SCIPconflictGetStrongbranchTime(scip->conflict),
      SCIPconflictGetNStrongbranchCalls(scip->conflict),
      SCIPconflictGetNStrongbranchConflictClauses(scip->conflict),
      SCIPconflictGetNStrongbranchConflictClauses(scip->conflict) > 0
      ? (Real)SCIPconflictGetNStrongbranchConflictLiterals(scip->conflict)
      / (Real)SCIPconflictGetNStrongbranchConflictClauses(scip->conflict) : 0,
      SCIPconflictGetNStrongbranchReconvergenceClauses(scip->conflict),
      SCIPconflictGetNStrongbranchReconvergenceClauses(scip->conflict) > 0
      ? (Real)SCIPconflictGetNStrongbranchReconvergenceLiterals(scip->conflict)
      / (Real)SCIPconflictGetNStrongbranchReconvergenceClauses(scip->conflict) : 0,
      SCIPconflictGetNStrongbranchIterations(scip->conflict));
   fprintf(file, "  pseudo solution  : %10.2f %10lld %10lld %10.1f %10lld %10.1f          -\n",
      SCIPconflictGetPseudoTime(scip->conflict),
      SCIPconflictGetNPseudoCalls(scip->conflict),
      SCIPconflictGetNPseudoConflictClauses(scip->conflict),
      SCIPconflictGetNPseudoConflictClauses(scip->conflict) > 0
      ? (Real)SCIPconflictGetNPseudoConflictLiterals(scip->conflict)
      / (Real)SCIPconflictGetNPseudoConflictClauses(scip->conflict) : 0,
      SCIPconflictGetNPseudoReconvergenceClauses(scip->conflict),
      SCIPconflictGetNPseudoReconvergenceClauses(scip->conflict) > 0
      ? (Real)SCIPconflictGetNPseudoReconvergenceLiterals(scip->conflict)
      / (Real)SCIPconflictGetNPseudoReconvergenceClauses(scip->conflict) : 0);
}

static
void printSeparatorStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(file != NULL);

   fprintf(file, "Separators         :       Time      Calls       Cuts\n");
   fprintf(file, "  cut pool         : %10.2f %10lld %10lld   (maximal pool size: %d)\n",
      SCIPcutpoolGetTime(scip->cutpool), 
      SCIPcutpoolGetNCalls(scip->cutpool), 
      SCIPcutpoolGetNCutsFound(scip->cutpool),
      SCIPcutpoolGetMaxNCuts(scip->cutpool));

   for( i = 0; i < scip->set->nsepas; ++i )
      fprintf(file, "  %-17.17s: %10.2f %10lld %10lld\n",
         SCIPsepaGetName(scip->set->sepas[i]),
         SCIPsepaGetTime(scip->set->sepas[i]),
         SCIPsepaGetNCalls(scip->set->sepas[i]),
         SCIPsepaGetNCutsFound(scip->set->sepas[i]));
}

static
void printPricerStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(file != NULL);

   fprintf(file, "Pricers            :       Time      Calls       Vars\n");
   fprintf(file, "  problem variables: %10.2f %10d %10d\n",
      SCIPpricestoreGetProbPricingTime(scip->pricestore),
      SCIPpricestoreGetNProbPricings(scip->pricestore),
      SCIPpricestoreGetNProbvarsFound(scip->pricestore));

   for( i = 0; i < scip->set->nactivepricers; ++i )
      fprintf(file, "  %-17.17s: %10.2f %10d %10d\n",
         SCIPpricerGetName(scip->set->pricers[i]),
         SCIPpricerGetTime(scip->set->pricers[i]),
         SCIPpricerGetNCalls(scip->set->pricers[i]),
         SCIPpricerGetNVarsFound(scip->set->pricers[i]));
}

static
void printBranchruleStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(file != NULL);

   fprintf(file, "Branching Rules    :       Time      Calls    Cutoffs    DomReds       Cuts      Conss   Children\n");

   for( i = 0; i < scip->set->nbranchrules; ++i )
      fprintf(file, "  %-17.17s: %10.2f %10lld %10lld %10lld %10lld %10lld %10lld\n",
         SCIPbranchruleGetName(scip->set->branchrules[i]),
         SCIPbranchruleGetTime(scip->set->branchrules[i]),
         SCIPbranchruleGetNLPCalls(scip->set->branchrules[i]) + SCIPbranchruleGetNPseudoCalls(scip->set->branchrules[i]),
         SCIPbranchruleGetNCutoffs(scip->set->branchrules[i]),
         SCIPbranchruleGetNDomredsFound(scip->set->branchrules[i]),
         SCIPbranchruleGetNCutsFound(scip->set->branchrules[i]),
         SCIPbranchruleGetNConssFound(scip->set->branchrules[i]),
         SCIPbranchruleGetNChildren(scip->set->branchrules[i]));
}

static
void printHeuristicStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->tree != NULL);
   assert(file != NULL);

   fprintf(file, "Primal Heuristics  :       Time      Calls      Found\n");
   fprintf(file, "  LP solutions     : %10.2f          - %10lld\n",
      SCIPclockGetTime(scip->stat->lpsoltime),
      scip->stat->nlpsolsfound);
   fprintf(file, "  pseudo solutions : %10.2f          - %10lld\n",
      SCIPclockGetTime(scip->stat->pseudosoltime),
      scip->stat->npssolsfound);

   SCIPsetSortHeurs(scip->set);

   for( i = 0; i < scip->set->nheurs; ++i )
      fprintf(file, "  %-17.17s: %10.2f %10lld %10lld\n",
         SCIPheurGetName(scip->set->heurs[i]),
         SCIPheurGetTime(scip->set->heurs[i]),
         SCIPheurGetNCalls(scip->set->heurs[i]),
         SCIPheurGetNSolsFound(scip->set->heurs[i]));
}

static
void printLPStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->lp != NULL);
   assert(file != NULL);

   fprintf(file, "LP                 :       Time      Calls Iterations  Iter/call   Iter/sec\n");

   fprintf(file, "  primal LP        : %10.2f %10d %10lld %10.2f",
      SCIPclockGetTime(scip->stat->primallptime),
      scip->stat->nprimallps,
      scip->stat->nprimallpiterations,
      scip->stat->nprimallps > 0 ? (Real)scip->stat->nprimallpiterations/(Real)scip->stat->nprimallps : 0.0);
   if( SCIPclockGetTime(scip->stat->primallptime) >= 0.01 )
      fprintf(file, " %10.2f\n", (Real)scip->stat->nprimallpiterations/SCIPclockGetTime(scip->stat->primallptime));
   else
      fprintf(file, "          -\n");

   fprintf(file, "  dual LP          : %10.2f %10d %10lld %10.2f",
      SCIPclockGetTime(scip->stat->duallptime),
      scip->stat->nduallps, 
      scip->stat->nduallpiterations,
      scip->stat->nduallps > 0 ? (Real)scip->stat->nduallpiterations/(Real)scip->stat->nduallps : 0.0);
   if( SCIPclockGetTime(scip->stat->duallptime) >= 0.01 )
      fprintf(file, " %10.2f\n", (Real)scip->stat->nduallpiterations/SCIPclockGetTime(scip->stat->duallptime));
   else
      fprintf(file, "          -\n");

   fprintf(file, "  diving LP        : %10.2f %10d %10lld %10.2f",
      SCIPclockGetTime(scip->stat->divinglptime),
      scip->stat->ndivinglps, 
      scip->stat->ndivinglpiterations,
      scip->stat->ndivinglps > 0 ? (Real)scip->stat->ndivinglpiterations/(Real)scip->stat->ndivinglps : 0.0);
   if( SCIPclockGetTime(scip->stat->divinglptime) >= 0.01 )
      fprintf(file, " %10.2f\n", (Real)scip->stat->ndivinglpiterations/SCIPclockGetTime(scip->stat->divinglptime));
   else
      fprintf(file, "          -\n");

   fprintf(file, "  strong branching : %10.2f %10d %10lld %10.2f",
      SCIPclockGetTime(scip->stat->strongbranchtime),
      scip->stat->nstrongbranchs,
      scip->stat->nsblpiterations,
      scip->stat->nstrongbranchs > 0 ? (Real)scip->stat->nsblpiterations/(Real)scip->stat->nstrongbranchs : 0.0);
   if( SCIPclockGetTime(scip->stat->strongbranchtime) >= 0.01 )
      fprintf(file, " %10.2f\n", (Real)scip->stat->nsblpiterations/SCIPclockGetTime(scip->stat->strongbranchtime));
   else
      fprintf(file, "          -\n");

   fprintf(file, "    (at root node) :          - %10d %10lld %10.2f          -\n",
      scip->stat->nrootstrongbranchs,
      scip->stat->nrootsblpiterations,
      scip->stat->nrootstrongbranchs > 0
      ? (Real)scip->stat->nrootsblpiterations/(Real)scip->stat->nrootstrongbranchs : 0.0);

   fprintf(file, "  conflict analysis: %10.2f %10d %10lld %10.2f",
      SCIPclockGetTime(scip->stat->conflictlptime),
      scip->stat->nconflictlps, 
      scip->stat->nconflictlpiterations,
      scip->stat->nconflictlps > 0 ? (Real)scip->stat->nconflictlpiterations/(Real)scip->stat->nconflictlps : 0.0);
   if( SCIPclockGetTime(scip->stat->conflictlptime) >= 0.01 )
      fprintf(file, " %10.2f\n", (Real)scip->stat->nconflictlpiterations/SCIPclockGetTime(scip->stat->conflictlptime));
   else
      fprintf(file, "          -\n");

}

static
void printRelaxatorStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(file != NULL);

   if( scip->set->nrelaxs == 0 )
      return;

   fprintf(file, "Relaxators         :       Time      Calls\n");

   for( i = 0; i < scip->set->nrelaxs; ++i )
      fprintf(file, "  %-17.17s: %10.2f %10lld\n",
         SCIPrelaxGetName(scip->set->relaxs[i]),
         SCIPrelaxGetTime(scip->set->relaxs[i]),
         SCIPrelaxGetNCalls(scip->set->relaxs[i]));
}

static
void printTreeStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->tree != NULL);
   assert(file != NULL);

   fprintf(file, "B&B Tree           :\n");
   fprintf(file, "  number of runs   : %10d\n", scip->stat->nruns);
   fprintf(file, "  nodes            : %10lld\n", scip->stat->nnodes);
   fprintf(file, "  nodes (total)    : %10lld\n", scip->stat->ntotalnodes);
   fprintf(file, "  max depth        : %10d\n", scip->stat->maxdepth);
   fprintf(file, "  max depth (total): %10d\n", scip->stat->maxtotaldepth);
   fprintf(file, "  backtracks       : %10lld (%.1f%%)\n", scip->stat->nbacktracks, 
      scip->stat->nnodes > 0 ? 100.0 * (Real)scip->stat->nbacktracks / (Real)scip->stat->nnodes : 0.0);
   fprintf(file, "  delayed cutoffs  : %10lld\n", scip->stat->ndelayedcutoffs);
   fprintf(file, "  repropagations   : %10lld (%lld domain reductions)\n", 
      scip->stat->nreprops, scip->stat->nrepropboundchgs);
   fprintf(file, "  avg switch length: %10.2f\n", 
      scip->stat->nnodes > 0
      ? (Real)(scip->stat->nactivatednodes + scip->stat->ndeactivatednodes) / (Real)scip->stat->nnodes : 0.0);
   fprintf(file, "  switching time   : %10.2f\n", SCIPclockGetTime(scip->stat->nodeactivationtime));
}

static
void printSolutionStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   Real primalbound;
   Real dualbound;
   Real dualboundroot;
   Real bestsol;
   Real gap;

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->primal != NULL);
   assert(file != NULL);

   primalbound = SCIPgetPrimalbound(scip);
   dualbound = SCIPgetDualbound(scip);
   dualboundroot = SCIPgetDualboundRoot(scip);
   gap = SCIPgetGap(scip);

   fprintf(file, "Solution           :\n");
   fprintf(file, "  Solutions found  : %10lld\n", scip->primal->nsolsfound);
   if( SCIPsetIsInfinity(scip->set, REALABS(primalbound)) )
   {
      if( scip->stage == SCIP_STAGE_SOLVED )
      {
         if( scip->primal->nsols == 0 )
            fprintf(file, "  Primal Bound     : infeasible\n");
         else
            fprintf(file, "  Primal Bound     :  unbounded\n");
      }
      else
         fprintf(file, "  Primal Bound     :          -\n");
   }
   else
   {
      fprintf(file, "  Primal Bound     : %+21.14e", primalbound);
      if( scip->primal->nsols == 0 )
         fprintf(file, "   (user objective limit)\n");
      else
      {
         bestsol = SCIPsolGetObj(scip->primal->sols[0]);
         bestsol = SCIPretransformObj(scip, bestsol);
         if( SCIPsetIsGT(scip->set, bestsol, primalbound) )
         {
            fprintf(file, "   (user objective limit)\n");
            fprintf(file, "  Best Solution    : %+21.14e", bestsol);
         }
         fprintf(file, "   (after %lld nodes, %.2f seconds, depth %d, found by <%s>)\n", 
            SCIPsolGetNodenum(scip->primal->sols[0]), 
            SCIPsolGetTime(scip->primal->sols[0]),
            SCIPsolGetDepth(scip->primal->sols[0]),
            SCIPsolGetHeur(scip->primal->sols[0]) != NULL
            ? SCIPheurGetName(SCIPsolGetHeur(scip->primal->sols[0])) : "relaxation");
      }
   }
   if( SCIPsetIsInfinity(scip->set, REALABS(dualbound)) )
      fprintf(file, "  Dual Bound       :          -\n");
   else
      fprintf(file, "  Dual Bound       : %+21.14e\n", dualbound);
   if( SCIPsetIsInfinity(scip->set, gap) )
      fprintf(file, "  Gap              :   infinite\n");
   else
      fprintf(file, "  Gap              : %10.2f %%\n", 100.0 * gap);
   if( SCIPsetIsInfinity(scip->set, REALABS(dualboundroot)) )
      fprintf(file, "  Root Dual Bound  :          -\n");
   else
      fprintf(file, "  Root Dual Bound  : %+21.14e\n", dualboundroot);
}

/** outputs solving statistics */
RETCODE SCIPprintStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintStatistics", TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( file == NULL )
      file = stdout;

   fprintf(file, "SCIP Status        : ");
   CHECK_OKAY( SCIPprintStage(scip, file) );
   fprintf(file, "\n");

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      fprintf(file, "Original Problem   : no problem exists.\n");
      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      fprintf(file, "Original Problem   :\n");
      SCIPprobPrintStatistics(scip->origprob, file);
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVED:
      fprintf(file, "Solving Time       : %10.2f\n", SCIPclockGetTime(scip->stat->solvingtime));
      fprintf(file, "Original Problem   :\n");
      SCIPprobPrintStatistics(scip->origprob, file);
      fprintf(file, "Presolved Problem  :\n");
      SCIPprobPrintStatistics(scip->transprob, file);
      printPresolverStatistics(scip, file);
      printConstraintStatistics(scip, file);
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      fprintf(file, "Solving Time       : %10.2f\n", SCIPclockGetTime(scip->stat->solvingtime));
      fprintf(file, "Original Problem   :\n");
      SCIPprobPrintStatistics(scip->origprob, file);
      fprintf(file, "Presolved Problem  :\n");
      SCIPprobPrintStatistics(scip->transprob, file);
      printPresolverStatistics(scip, file);
      printConstraintStatistics(scip, file);
      printConstraintTimingStatistics(scip, file);
      printPropagatorStatistics(scip, file);
      printConflictStatistics(scip, file);
      printSeparatorStatistics(scip, file);
      printPricerStatistics(scip, file);
      printBranchruleStatistics(scip, file);
      printHeuristicStatistics(scip, file);
      printLPStatistics(scip, file);
      printRelaxatorStatistics(scip, file);
      printTreeStatistics(scip, file);
      printSolutionStatistics(scip, file);
      
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** outputs history statistics about branchings on variables */
RETCODE SCIPprintBranchingStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   VAR** vars;
   Real* depths;
   Real depth;
   int v;
   int i;

   CHECK_OKAY( checkStage(scip, "SCIPprintBranchingHistory", TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( file == NULL )
      file = stdout;

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      fprintf(file, "problem not yet solved. branching statistics not available.\n");
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      CHECK_OKAY( SCIPallocBufferArray(scip, &vars, scip->transprob->nvars) );
      CHECK_OKAY( SCIPallocBufferArray(scip, &depths, scip->transprob->nvars) );
      for( v = 0; v < scip->transprob->nvars; ++v )
      {
         depth = SCIPvarGetAvgBranchdepth(scip->transprob->vars[v], SCIP_BRANCHDIR_DOWNWARDS)
            + SCIPvarGetAvgBranchdepth(scip->transprob->vars[v], SCIP_BRANCHDIR_UPWARDS);
         for( i = v; i > 0 && depth < depths[i-1]; i-- )
         {
            vars[i] = vars[i-1];
            depths[i] = depths[i-1];
         }
         vars[i] = scip->transprob->vars[v];
         depths[i] = depth;
      }

      fprintf(file, "                                                  branchings        inferences         cutoffs               LP gain  \n");
      fprintf(file, " variable        priority   factor    depth     down       up     down       up     down       up       down         up\n");

      for( v = 0; v < scip->transprob->nvars; ++v )
      {
         if( SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_DOWNWARDS) > 0
            || SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_UPWARDS) > 0 )
         {
            fprintf(file, " %-16s %7d %8.1f %8.1f %8lld %8lld %8.1f %8.1f %7.1f%% %7.1f%% %10.1f %10.1f\n",
               SCIPvarGetName(vars[v]),
               SCIPvarGetBranchPriority(vars[v]),
               SCIPvarGetBranchFactor(vars[v]),
               (SCIPvarGetAvgBranchdepth(vars[v], SCIP_BRANCHDIR_DOWNWARDS)
                  + SCIPvarGetAvgBranchdepth(vars[v], SCIP_BRANCHDIR_UPWARDS))/2.0 - 1.0,
               SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_DOWNWARDS),
               SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_UPWARDS),
               SCIPvarGetAvgInferences(vars[v], scip->stat, SCIP_BRANCHDIR_DOWNWARDS),
               SCIPvarGetAvgInferences(vars[v], scip->stat, SCIP_BRANCHDIR_UPWARDS),
               100.0 * SCIPvarGetAvgCutoffs(vars[v], scip->stat, SCIP_BRANCHDIR_DOWNWARDS),
               100.0 * SCIPvarGetAvgCutoffs(vars[v], scip->stat, SCIP_BRANCHDIR_UPWARDS),
               SCIPvarGetPseudocost(vars[v], scip->stat, -1.0),
               SCIPvarGetPseudocost(vars[v], scip->stat, +1.0));
         }
      }
      SCIPfreeBufferArray(scip, &depths);
      SCIPfreeBufferArray(scip, &vars);

      return SCIP_OKAY;
      
   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** outputs node information display line */
RETCODE SCIPprintDisplayLine(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file,               /**< output file (or NULL for standard output) */
   VERBLEVEL        verblevel           /**< minimal verbosity level to actually display the information line */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintDisplayLine", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( (VERBLEVEL)scip->set->disp_verblevel >= verblevel )
   {
      CHECK_OKAY( SCIPdispPrintLine(scip->set, scip->stat, TRUE) );
   }

   return SCIP_OKAY;
}




/*
 * timing methods
 */

/** gets current time of day in seconds (standard time zone) */
Real SCIPgetTimeOfDay(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetTimeOfDay", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPclockGetTimeOfDay();
}

/** creates a clock using the default clock type */
RETCODE SCIPcreateClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPclockCreate(clock, SCIP_CLOCKTYPE_DEFAULT) );

   return SCIP_OKAY;
}

/** creates a clock counting the CPU user seconds */
RETCODE SCIPcreateCPUClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateCPUClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPclockCreate(clock, SCIP_CLOCKTYPE_CPU) );

   return SCIP_OKAY;
}

/** creates a clock counting the wall clock seconds */
RETCODE SCIPcreateWallClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateWallClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPclockCreate(clock, SCIP_CLOCKTYPE_WALL) );

   return SCIP_OKAY;
}

/** frees a clock */
RETCODE SCIPfreeClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockFree(clock);

   return SCIP_OKAY;
}

/** resets the time measurement of a clock to zero and completely stops the clock */
RETCODE SCIPresetClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPresetClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockReset(clock);

   return SCIP_OKAY;
}

/** starts the time measurement of a clock */
RETCODE SCIPstartClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPstartClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockStart(clock, scip->set);

   return SCIP_OKAY;
}

/** stops the time measurement of a clock */
RETCODE SCIPstopClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPstopClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockStop(clock, scip->set);

   return SCIP_OKAY;
}

/** gets the measured time of a clock in seconds */
Real SCIPgetClockTime(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetClockTime", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPclockGetTime(clock);
}

/** sets the measured time of a clock to the given value in seconds */
RETCODE SCIPsetClockTime(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock,              /**< clock timer */
   Real             sec                 /**< time in seconds to set the clock's timer to */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetClockTime", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockSetTime(clock, sec);

   return SCIP_OKAY;
}

/** gets the current total SCIP time in seconds */
Real SCIPgetTotalTime(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetTotalTime", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPclockGetTime(scip->totaltime);
}

/** gets the current solving time in seconds */
Real SCIPgetSolvingTime(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolvingTime", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPclockGetTime(scip->stat->solvingtime);
}

/** gets the current presolving time in seconds */
Real SCIPgetPresolvingTime(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPresolvingTime", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPclockGetTime(scip->stat->presolvingtime);
}




/*
 * numeric values and comparisons
 */

/** returns value treated as infinity */
Real SCIPinfinity(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetInfinity(scip->set);
}

/** returns value treated as zero */
Real SCIPepsilon(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetEpsilon(scip->set);
}

/** returns value treated as zero for sums of floating point values */
Real SCIPsumepsilon(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetSumepsilon(scip->set);
}

/** returns feasibility tolerance for constraints */
Real SCIPfeastol(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
 
   return SCIPsetFeastol(scip->set);
}

/** returns feasibility tolerance for reduced costs */
Real SCIPdualfeastol(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
 
   return SCIPsetDualfeastol(scip->set);
}

/** sets the feasibility tolerance for constraints */
RETCODE SCIPchgFeastol(
   SCIP*            scip,               /**< SCIP data structure */
   Real             feastol             /**< new feasibility tolerance for constraints */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgFeastol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* mark the LP unsolved, if the feasibility tolerance was tightened */
   if( scip->lp != NULL && feastol < SCIPsetFeastol(scip->set) )
      scip->lp->solved = FALSE;

   /* change the settings */
   CHECK_OKAY( SCIPsetSetFeastol(scip->set, feastol) );

   return SCIP_OKAY;
}

/** sets the feasibility tolerance for reduced costs */
RETCODE SCIPchgDualfeastol(
   SCIP*            scip,               /**< SCIP data structure */
   Real             dualfeastol         /**< new feasibility tolerance for reduced costs */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgDualfeastol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* mark the LP unsolved, if the dual feasibility tolerance was tightened */
   if( scip->lp != NULL && dualfeastol < SCIPsetDualfeastol(scip->set) )
      scip->lp->solved = FALSE;

   /* change the settings */
   CHECK_OKAY( SCIPsetSetDualfeastol(scip->set, dualfeastol) );

   return SCIP_OKAY;
}

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** checks, if values are in range of epsilon */
Bool SCIPisEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsEQ(scip->set, val1, val2);
}

/** checks, if val1 is (more than epsilon) lower than val2 */
Bool SCIPisLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsLT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than epsilon) greater than val2 */
Bool SCIPisLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsLE(scip->set, val1, val2);
}

/** checks, if val1 is (more than epsilon) greater than val2 */
Bool SCIPisGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsGT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than epsilon) lower than val2 */
Bool SCIPisGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsGE(scip->set, val1, val2);
}

/** checks, if value is in range epsilon of 0.0 */
Bool SCIPisZero(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsZero(scip->set, val);
}

/** checks, if value is greater than epsilon */
Bool SCIPisPositive(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsPositive(scip->set, val);
}

/** checks, if value is lower than -epsilon */
Bool SCIPisNegative(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsNegative(scip->set, val);
}

/** checks, if values are in range of sumepsilon */
Bool SCIPisSumEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumEQ(scip->set, val1, val2);
}

/** checks, if val1 is (more than sumepsilon) lower than val2 */
Bool SCIPisSumLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumLT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
Bool SCIPisSumLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumLE(scip->set, val1, val2);
}

/** checks, if val1 is (more than sumepsilon) greater than val2 */
Bool SCIPisSumGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumGT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
Bool SCIPisSumGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumGE(scip->set, val1, val2);
}

/** checks, if value is in range sumepsilon of 0.0 */
Bool SCIPisSumZero(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumZero(scip->set, val);
}

/** checks, if value is greater than sumepsilon */
Bool SCIPisSumPositive(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumPositive(scip->set, val);
}

/** checks, if value is lower than -sumepsilon */
Bool SCIPisSumNegative(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumNegative(scip->set, val);
}

/** checks, if values are in range of feasibility tolerance */
Bool SCIPisFeasEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasEQ(scip->set, val1, val2);
}

/** checks, if val1 is (more than feasibility tolerance) lower than val2 */
Bool SCIPisFeasLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasLT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than feasibility tolerance) greater than val2 */
Bool SCIPisFeasLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasLE(scip->set, val1, val2);
}

/** checks, if val1 is (more than feasibility tolerance) greater than val2 */
Bool SCIPisFeasGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasGT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than feasibility tolerance) lower than val2 */
Bool SCIPisFeasGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasGE(scip->set, val1, val2);
}

/** checks, if value is in range feasibility tolerance of 0.0 */
Bool SCIPisFeasZero(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasZero(scip->set, val);
}

/** checks, if value is greater than feasibility tolerance */
Bool SCIPisFeasPositive(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasPositive(scip->set, val);
}

/** checks, if value is lower than -feasibility tolerance */
Bool SCIPisFeasNegative(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasNegative(scip->set, val);
}

/** checks, if the first given lower bound is tighter (w.r.t. bound strengthening epsilon) than the second one */
Bool SCIPisLbBetter(
   SCIP*            scip,               /**< SCIP data structure */
   Real             lb1,                /**< first lower bound to compare */
   Real             lb2                 /**< second lower bound to compare */
   )
{
   assert(scip != NULL);

   return SCIPsetIsLbBetter(scip->set, lb1, lb2);
}

/** checks, if the first given upper bound is tighter (w.r.t. bound strengthening epsilon) than the second one */
Bool SCIPisUbBetter(
   SCIP*            scip,               /**< SCIP data structure */
   Real             ub1,                /**< first upper bound to compare */
   Real             ub2                 /**< second upper bound to compare */
   )
{
   assert(scip != NULL);

   return SCIPsetIsUbBetter(scip->set, ub1, ub2);
}

/** checks, if relative difference of values is in range of epsilon */
Bool SCIPisRelEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelEQ(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is lower than epsilon */
Bool SCIPisRelLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelLT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
Bool SCIPisRelLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelLE(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is greater than epsilon */
Bool SCIPisRelGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelGT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
Bool SCIPisRelGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelGE(scip->set, val1, val2);
}

/** checks, if rel. difference of values is in range of sumepsilon */
Bool SCIPisSumRelEQ(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelEQ(scip->set, val1, val2);
}

/** checks, if rel. difference of val1 and val2 is lower than sumepsilon */
Bool SCIPisSumRelLT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelLT(scip->set, val1, val2);
}

/** checks, if rel. difference of val1 and val2 is not greater than sumepsilon */
Bool SCIPisSumRelLE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelLE(scip->set, val1, val2);
}

/** checks, if rel. difference of val1 and val2 is greater than sumepsilon */
Bool SCIPisSumRelGT(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelGT(scip->set, val1, val2);
}

/** checks, if rel. difference of val1 and val2 is not lower than -sumepsilon */
Bool SCIPisSumRelGE(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelGE(scip->set, val1, val2);
}

/** checks, if value is (positive) infinite */
Bool SCIPisInfinity(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against infinity */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsInfinity(scip->set, val);
}

/** checks, if value is non-negative within the LP feasibility bounds */
Bool SCIPisFeasible(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasible(scip->set, val);
}

/** checks, if value is integral within the LP feasibility bounds */
Bool SCIPisIntegral(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   if( scip->set->misc_exactsolve )
      return (val == SCIPsetFloor(scip->set, val));
   else
      return SCIPsetIsIntegral(scip->set, val);
}

/** checks, if given fractional part is smaller than feastol */
Bool SCIPisFracIntegral(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   if( scip->set->misc_exactsolve )
      return (val == 0.0);
   else
      return SCIPsetIsFracIntegral(scip->set, val);
}

/** rounds value + feasibility tolerance down to the next integer */
Real SCIPfloor(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFloor(scip->set, val);
}

/** rounds value - feasibility tolerance up to the next integer */
Real SCIPceil(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to be compared against zero */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetCeil(scip->set, val);
}

/** returns fractional part of value, i.e. x - floor(x) */
Real SCIPfrac(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val                 /**< value to return fractional part for */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFrac(scip->set, val);
}

#endif

/** outputs a real number, or "+infinity", or "-infinity" to a file */
void SCIPprintReal(
   SCIP*            scip,               /**< SCIP data structure */
   Real             val,                /**< value to print */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   assert(scip != NULL);

   if( file == NULL )
      file = stdout;

   if( SCIPsetIsInfinity(scip->set, val) )
      fprintf(file, "+infinity");
   else if( SCIPsetIsInfinity(scip->set, -val) )
      fprintf(file, "-infinity");
   else
      fprintf(file, "%f", val);
}




/*
 * memory management
 */

/** returns block memory to use at the current time */
MEMHDR* SCIPmemhdr(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPmemhdr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      return scip->mem->probmem;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      return scip->mem->solvemem;

   default:
      errorMessage("invalid SCIP stage\n");
      return NULL;
   }  /*lint !e788*/
}

/** returns the total number of bytes used in block memory */
Longint SCIPgetMemUsed(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetMemUsed", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPmemGetUsed(scip->mem);
}

/** calculate memory size for dynamically allocated arrays */
int SCIPcalcMemGrowSize(
   SCIP*            scip,               /**< SCIP data structure */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);

   return SCIPsetCalcMemGrowSize(scip->set, num);
}

/** extends a dynamically allocated block memory array to be able to store at least the given number of elements;
 *  use SCIPensureBlockMemoryArray() define to call this method!
 */
RETCODE SCIPensureBlockMemoryArray_call(
   SCIP*            scip,               /**< SCIP data structure */
   void**           arrayptr,           /**< pointer to dynamically sized array */
   size_t           elemsize,           /**< size in bytes of each element in array */
   int*             arraysize,          /**< pointer to current array size */
   int              minsize             /**< required minimal array size */
   )
{
   assert(scip != NULL);
   assert(arrayptr != NULL);
   assert(elemsize > 0);
   assert(arraysize != NULL);

   if( minsize > *arraysize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(scip->set, minsize);
      ALLOC_OKAY( reallocBlockMemorySize(SCIPmemhdr(scip), arrayptr, *arraysize * elemsize, newsize * elemsize) );
      *arraysize = newsize;
   }

   return SCIP_OKAY;
}

/** gets a memory buffer with at least the given size */
RETCODE SCIPallocBufferSize(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   int              size                /**< required size in bytes of buffer */
   )
{
   assert(ptr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPallocBufferSize", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetAllocBufferSize(scip->set, ptr, size) );

   return SCIP_OKAY;
}

/** allocates a memory buffer with at least the given size and copies the given memory into the buffer */
RETCODE SCIPduplicateBufferSize(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   void*            source,             /**< memory block to copy into the buffer */
   int              size                /**< required size in bytes of buffer */
   )
{
   assert(ptr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPduplicateBufferSize", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetDuplicateBufferSize(scip->set, ptr, source, size) );

   return SCIP_OKAY;
}

/** reallocates a memory buffer to at least the given size */
RETCODE SCIPreallocBufferSize(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to the buffer */
   int              size                /**< required size in bytes of buffer */
   )
{
   assert(ptr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreallocBufferSize", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetReallocBufferSize(scip->set, ptr, size) );

   return SCIP_OKAY;
}

/** frees a memory buffer */
void SCIPfreeBufferSize(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to the buffer */
   int              dummysize           /**< used to get a safer define for SCIPfreeBuffer() and SCIPfreeBufferArray() */
   )
{  /*lint --e{715}*/
   assert(ptr != NULL);

   CHECK_ABORT( checkStage(scip, "SCIPfreeBufferSize", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetFreeBufferSize(scip->set, ptr);
}




/*
 * dynamic arrays
 */

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** creates a dynamic array of real values */
RETCODE SCIPcreateRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY**      realarray           /**< pointer to store the real array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayCreate(realarray, SCIPmemhdr(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
RETCODE SCIPfreeRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY**      realarray           /**< pointer to the real array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayFree(realarray) );
   
   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPextendRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPextendRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayExtend(realarray, scip->set, minidx, maxidx) );
   
   return SCIP_OKAY;
}

/** clears a dynamic real array */
RETCODE SCIPclearRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray           /**< dynamic real array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayClear(realarray) );
   
   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
Real SCIPgetRealarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx                 /**< array index to get value for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRealarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPrealarrayGetVal(realarray, idx);
}

/** sets value of entry in dynamic array */
RETCODE SCIPsetRealarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx,                /**< array index to set value for */
   Real             val                 /**< value to set array index to */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetRealarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarraySetVal(realarray, scip->set, idx, val) );
   
   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
RETCODE SCIPincRealarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx,                /**< array index to increase value for */
   Real             incval              /**< value to increase array index */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPincRealarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayIncVal(realarray, scip->set, idx, incval) );
   
   return SCIP_OKAY;
}

/** creates a dynamic array of int values */
RETCODE SCIPcreateIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY**       intarray            /**< pointer to store the int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayCreate(intarray, SCIPmemhdr(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of int values */
RETCODE SCIPfreeIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY**       intarray            /**< pointer to the int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayFree(intarray) );
   
   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPextendIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPextendIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayExtend(intarray, scip->set, minidx, maxidx) );
   
   return SCIP_OKAY;
}

/** clears a dynamic int array */
RETCODE SCIPclearIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray            /**< dynamic int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayClear(intarray) );
   
   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
int SCIPgetIntarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx                 /**< array index to get value for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetIntarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPintarrayGetVal(intarray, idx);
}

/** sets value of entry in dynamic array */
RETCODE SCIPsetIntarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx,                /**< array index to set value for */
   int              val                 /**< value to set array index to */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetIntarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarraySetVal(intarray, scip->set, idx, val) );
   
   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
RETCODE SCIPincIntarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx,                /**< array index to increase value for */
   int              incval              /**< value to increase array index */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPincIntarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayIncVal(intarray, scip->set, idx, incval) );
   
   return SCIP_OKAY;
}

/** creates a dynamic array of bool values */
RETCODE SCIPcreateBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY**      boolarray           /**< pointer to store the bool array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarrayCreate(boolarray, SCIPmemhdr(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of bool values */
RETCODE SCIPfreeBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY**      boolarray           /**< pointer to the bool array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarrayFree(boolarray) );
   
   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPextendBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPextendBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarrayExtend(boolarray, scip->set, minidx, maxidx) );
   
   return SCIP_OKAY;
}

/** clears a dynamic bool array */
RETCODE SCIPclearBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarrayClear(boolarray) );
   
   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
Bool SCIPgetBoolarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx                 /**< array index to get value for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBoolarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPboolarrayGetVal(boolarray, idx);
}

/** sets value of entry in dynamic array */
RETCODE SCIPsetBoolarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx,                /**< array index to set value for */
   Bool             val                 /**< value to set array index to */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetBoolarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarraySetVal(boolarray, scip->set, idx, val) );
   
   return SCIP_OKAY;
}

/** creates a dynamic array of pointers */
RETCODE SCIPcreatePtrarray(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY**       ptrarray            /**< pointer to store the int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreatePtrarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPptrarrayCreate(ptrarray, SCIPmemhdr(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of pointers */
RETCODE SCIPfreePtrarray(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY**       ptrarray            /**< pointer to the int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreePtrarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPptrarrayFree(ptrarray) );
   
   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
RETCODE SCIPextendPtrarray(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY*        ptrarray,           /**< dynamic int array */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPextendPtrarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPptrarrayExtend(ptrarray, scip->set, minidx, maxidx) );
   
   return SCIP_OKAY;
}

/** clears a dynamic pointer array */
RETCODE SCIPclearPtrarray(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY*        ptrarray            /**< dynamic int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearPtrarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPptrarrayClear(ptrarray) );
   
   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
void* SCIPgetPtrarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY*        ptrarray,           /**< dynamic int array */
   int              idx                 /**< array index to get value for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPtrarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPptrarrayGetVal(ptrarray, idx);
}

/** sets value of entry in dynamic array */
RETCODE SCIPsetPtrarrayVal(
   SCIP*            scip,               /**< SCIP data structure */
   PTRARRAY*        ptrarray,           /**< dynamic int array */
   int              idx,                /**< array index to set value for */
   void*            val                 /**< value to set array index to */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetPtrarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPptrarraySetVal(ptrarray, scip->set, idx, val) );
   
   return SCIP_OKAY;
}

#endif




/*
 * debugging methods
 */

#ifndef NDEBUG

/** prints output about used memory */
void SCIPdebugMemory(
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
