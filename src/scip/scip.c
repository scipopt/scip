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
#pragma ident "@(#) $Id: scip.c,v 1.147 2004/04/08 13:14:44 bzfwolte Exp $"

/**@file   scip.c
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 */

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
#include "sepa.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "struct_scip.h"
#endif



/** checks, if SCIP is in one of the feasible stages */
static
RETCODE checkStage(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      method,             /**< method that was called */
   Bool             init,               /**< may method be called in the INIT stage? */
   Bool             problem,            /**< may method be called in the PROBLEM stage? */
   Bool             initsolve,          /**< may method be called in the INITSOLVE stage? */
   Bool             presolving,         /**< may method be called in the PRESOLVING stage? */
   Bool             presolved,          /**< may method be called in the PRESOLVED stage? */
   Bool             solving,            /**< may method be called in the SOLVING stage? */
   Bool             solved,             /**< may method be called in the SOLVED stage? */
   Bool             freesolve           /**< may method be called in the FREESOLVE stage? */
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
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->transprob == NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->lpconflict == NULL);
      assert(scip->pseudoconflict == NULL);
      assert(scip->primal == NULL);

      if( !init )
      {
         errorMessage("Cannot call method <%s> in initialization stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter == NULL);
      assert(scip->eventqueue == NULL);
      assert(scip->branchcand == NULL);
      assert(scip->tree == NULL);
      assert(scip->lp == NULL);
      assert(scip->transprob == NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->lpconflict == NULL);
      assert(scip->pseudoconflict == NULL);
      assert(scip->primal == NULL);

      if( !problem )
      {
         errorMessage("Cannot call method <%s> in problem creation stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->lpconflict == NULL);
      assert(scip->pseudoconflict == NULL);
      assert(scip->primal == NULL);

      if( !initsolve )
      {
         errorMessage("Cannot call method <%s> in solve initialization stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      assert(scip->conflict == NULL);
      assert(scip->lpconflict == NULL);
      assert(scip->pseudoconflict == NULL);
      assert(scip->primal == NULL);

      if( !presolving )
      {
         errorMessage("Cannot call method <%s> in presolving stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore != NULL);
      assert(scip->sepastore != NULL);
      assert(scip->cutpool != NULL);
      assert(scip->conflict != NULL);
      assert(scip->lpconflict != NULL);
      assert(scip->pseudoconflict != NULL);
      assert(scip->primal != NULL);

      if( !presolved )
      {
         errorMessage("Cannot call method <%s> in problem presolved stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore != NULL);
      assert(scip->sepastore != NULL);
      assert(scip->cutpool != NULL);
      assert(scip->conflict != NULL);
      assert(scip->lpconflict != NULL);
      assert(scip->pseudoconflict != NULL);
      assert(scip->primal != NULL);

      if( !solving )
      {
         errorMessage("Cannot call method <%s> in solving stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->tree != NULL);
      assert(scip->lp != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore != NULL);
      assert(scip->sepastore != NULL);
      assert(scip->cutpool != NULL);
      assert(scip->conflict != NULL);
      assert(scip->lpconflict != NULL);
      assert(scip->pseudoconflict != NULL);
      assert(scip->primal != NULL);

      if( !solved )
      {
         errorMessage("Cannot call method <%s> in problem solved stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_FREESOLVE:
      assert(scip->origprob != NULL);
      assert(scip->stat != NULL);

      if( !freesolve )
      {
         errorMessage("Cannot call method <%s> in solve deinitialization stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   default:
      errorMessage("Unknown SCIP stage\n");
      return SCIP_ERROR;
   }
}




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
   (*scip)->lpconflict = NULL;
   (*scip)->pseudoconflict = NULL;
   (*scip)->primal = NULL;

   return SCIP_OKAY;
}

/** frees SCIP data structures */
RETCODE SCIPfree(
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   CHECK_OKAY( checkStage(*scip, "SCIPfree", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

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

   CHECK_ABORT( checkStage(scip, "SCIPmessage", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   va_start(ap, formatstr);
   vinfoMessage(scip->set->verblevel, msgverblevel, formatstr, ap);
   va_end(ap);
}

/** returns current stage of SCIP */
STAGE SCIPstage(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->stage;
}

/** returns whether the current stage belongs to the transformed problem space */
Bool SCIPisTransformed(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return (scip->stage != SCIP_STAGE_PROBLEM);
}

/** returns whether the solution process should be provably correct */
Bool SCIPisExactSolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return (scip->set->exactsolve);
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
   CHECK_OKAY( checkStage(scip, "SCIPaddBoolParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPaddIntParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPaddLongintParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPaddRealParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPaddCharParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPaddStringParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPgetBoolParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPgetIntParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPgetLongintParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPgetRealParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPgetCharParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPgetStringParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPsetBoolParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPsetIntParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPsetLongintParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPsetRealParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPsetCharParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPsetStringParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetSetStringParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** reads parameters from a file */
RETCODE SCIPreadParams(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< file name */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPreadParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPwriteParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetWriteParams(scip->set, filename, comments, onlychanged) );

   return SCIP_OKAY;
}

/** returns the array of all available SCIP parameters */
PARAM** SCIPgetParams(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPsetGetParams(scip->set);
}

/** returns the total number of all available SCIP parameters */
int SCIPgetNParams(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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

   CHECK_OKAY( checkStage(scip, "SCIPincludeReader", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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

   CHECK_ABORT( checkStage(scip, "SCIPfindReader", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindReader(scip->set, name);
}

/** returns the array of currently available readers */
READER** SCIPgetReaders(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetReaders", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->readers;
}

/** returns the number of currently available readers */
int SCIPgetNReaders(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNReaders", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPincludePricer", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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

   CHECK_ABORT( checkStage(scip, "SCIPfindPricer", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindPricer(scip->set, name);
}

/** returns the array of currently available variable pricers; active pricers are in the first slots of the array */
PRICER** SCIPgetPricers(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPricers", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortPricers(scip->set);

   return scip->set->pricers;
}

/** returns the number of currently available variable pricers */
int SCIPgetNPricers(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPricers", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->npricers;
}

/** returns the number of currently active variable pricers, that are used in the LP solving loop */
int SCIPgetNActivePricers(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNAcvitePricers", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nactivepricers;
}

/** sets the priority of a variable pricer */
RETCODE SCIPsetPricerPriority(
   SCIP*            scip,               /**< SCIP data structure */
   PRICER*          pricer,             /**< variable pricer */
   int              priority            /**< new priority of the variable pricer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetPricerPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPactivatePricer", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   DECL_CONSINIT    ((*consinit)),      /**< initialize constraint handler */
   DECL_CONSEXIT    ((*consexit)),      /**< deinitialize constraint handler */
   DECL_CONSSOLSTART((*conssolstart)),  /**< solving start notification method of constraint handler */
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
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   CONSHDLR* conshdlr;

   CHECK_OKAY( checkStage(scip, "SCIPincludeConshdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconshdlrCreate(&conshdlr, scip->set, scip->mem->setmem,
                  name, desc, sepapriority, enfopriority, chckpriority, sepafreq, propfreq, needscons, 
                  consfree, consinit, consexit, conssolstart, consdelete, constrans, consinitlp,
                  conssepa, consenfolp, consenfops, conscheck, consprop, conspresol, consrescvar,
                  conslock, consunlock, consactive, consdeactive, consenable, consdisable, conshdlrdata) );
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

   CHECK_ABORT( checkStage(scip, "SCIPfindConshdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindConshdlr(scip->set, name);
}

/** returns the array of currently available constraint handlers */
CONSHDLR** SCIPgetConshdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetConshdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->conshdlrs;
}

/** returns the number of currently available constraint handlers */
int SCIPgetNConshdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNConshdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPincludeConflicthdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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

   CHECK_ABORT( checkStage(scip, "SCIPfindConflicthdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindConflicthdlr(scip->set, name);
}

/** returns the array of currently available conflict handlers */
CONFLICTHDLR** SCIPgetConflicthdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetConflicthdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortConflicthdlrs(scip->set);

   return scip->set->conflicthdlrs;
}

/** returns the number of currently available conflict handlers */
int SCIPgetNConflicthdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNConflicthdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nconflicthdlrs;
}

/** sets the priority of a conflict handler */
RETCODE SCIPsetConflicthdlrPriority(
   SCIP*            scip,               /**< SCIP data structure */
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   int              priority            /**< new priority of the conflict handler */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetConflicthdlrPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPconflicthdlrSetPriority(conflicthdlr, scip->set, priority);

   return SCIP_OKAY;
}

/** creates a presolver and includes it in SCIP */
RETCODE SCIPincludePresol(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name,               /**< name of presolver */
   const char*      desc,               /**< description of presolver */
   int              priority,           /**< priority of the presolver */
   DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver */
   DECL_PRESOLINIT  ((*presolinit)),    /**< initialize presolver */
   DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialize presolver */
   DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   PRESOL* presol;

   CHECK_OKAY( checkStage(scip, "SCIPincludePresol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPpresolCreate(&presol, scip->set, scip->mem->setmem, name, desc, priority,
                  presolfree, presolinit, presolexit, presolexec, presoldata) );
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

   CHECK_ABORT( checkStage(scip, "SCIPfindPresol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindPresol(scip->set, name);
}

/** returns the array of currently available presolvers */
PRESOL** SCIPgetPresols(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPresols", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortPresols(scip->set);

   return scip->set->presols;
}

/** returns the number of currently available presolvers */
int SCIPgetNPresols(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPresols", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->npresols;
}

/** sets the priority of a presolver */
RETCODE SCIPsetPresolPriority(
   SCIP*            scip,               /**< SCIP data structure */
   PRESOL*          presol,             /**< presolver */
   int              priority            /**< new priority of the presolver */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetPresolPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPpresolSetPriority(presol, scip->set, priority);

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

   CHECK_OKAY( checkStage(scip, "SCIPincludeSepa", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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

   CHECK_ABORT( checkStage(scip, "SCIPfindSepa", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindSepa(scip->set, name);
}

/** returns the array of currently available separators */
SEPA** SCIPgetSepas(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSepas", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortSepas(scip->set);

   return scip->set->sepas;
}

/** returns the number of currently available separators */
int SCIPgetNSepas(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNSepas", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nsepas;
}

/** sets the priority of a separator */
RETCODE SCIPsetSepaPriority(
   SCIP*            scip,               /**< SCIP data structure */
   SEPA*            sepa,               /**< primal sepaistic */
   int              priority            /**< new priority of the separator */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetSepaPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsepaSetPriority(sepa, scip->set, priority);

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
   Bool             pseudonodes,        /**< call heuristic at nodes where only a pseudo solution exist? */
   DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   DECL_HEURINIT    ((*heurinit)),      /**< initialize primal heuristic */
   DECL_HEUREXIT    ((*heurexit)),      /**< deinitialize primal heuristic */
   DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   HEUR* heur;

   CHECK_OKAY( checkStage(scip, "SCIPincludeHeur", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPheurCreate(&heur, scip->set, scip->mem->setmem,
                  name, desc, dispchar, priority, freq, freqofs, pseudonodes,
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

   CHECK_ABORT( checkStage(scip, "SCIPfindHeur", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindHeur(scip->set, name);
}

/** returns the array of currently available primal heuristics */
HEUR** SCIPgetHeurs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetHeurs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortHeurs(scip->set);

   return scip->set->heurs;
}

/** returns the number of currently available primal heuristics */
int SCIPgetNHeurs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNHeurs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nheurs;
}

/** sets the priority of a primal heuristic */
RETCODE SCIPsetHeurPriority(
   SCIP*            scip,               /**< SCIP data structure */
   HEUR*            heur,               /**< primal heuristic */
   int              priority            /**< new priority of the primal heuristic */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetHeurPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPincludeEventhdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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

   CHECK_ABORT( checkStage(scip, "SCIPfindEventhdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindEventhdlr(scip->set, name);
}

/** returns the array of currently available event handlers */
EVENTHDLR** SCIPgetEventhdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetEventhdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->eventhdlrs;
}

/** returns the number of currently available event handlers */
int SCIPgetNEventhdlrs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNEventhdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPincludeNodesel", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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

   CHECK_ABORT( checkStage(scip, "SCIPfindNodesel", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindNodesel(scip->set, name);
}

/** returns the array of currently available node selectors */
NODESEL** SCIPgetNodesels(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNodesels", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nodesels;
}

/** returns the number of currently available node selectors */
int SCIPgetNNodesels(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNNodesels", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nnodesels;
}

/** sets the priority of a node selector in standard mode */
RETCODE SCIPsetNodeselStdPriority(
   SCIP*            scip,               /**< SCIP data structure */
   NODESEL*         nodesel,            /**< node selector */
   int              priority            /**< new standard priority of the node selector */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetNodeselStdPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPsetNodeselMemsavePriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPnodeselSetMemsavePriority(nodesel, scip->set, priority);

   return SCIP_OKAY;
}

/** returns the currently used node selector */
NODESEL* SCIPgetNodesel(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNodesel", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetGetNodesel(scip->set, scip->stat);
}

/** creates a branching rule and includes it in SCIP */
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
   )
{
   BRANCHRULE* branchrule;

   CHECK_OKAY( checkStage(scip, "SCIPincludeBranchrule", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchruleCreate(&branchrule, scip->mem->setmem, scip->set, name, desc, priority, maxdepth,
                  branchfree, branchinit, branchexit, branchexeclp, branchexecps, branchruledata) );
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

   CHECK_ABORT( checkStage(scip, "SCIPfindBranchrule", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortBranchrules(scip->set);

   return SCIPsetFindBranchrule(scip->set, name);
}

/** returns the array of currently available branching rules */
BRANCHRULE** SCIPgetBranchrules(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBranchrules", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->branchrules;
}

/** returns the number of currently available branching rules */
int SCIPgetNBranchrules(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNBranchrules", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nbranchrules;
}

/** sets the priority of a branching rule */
RETCODE SCIPsetBranchrulePriority(
   SCIP*            scip,               /**< SCIP data structure */
   BRANCHRULE*      branchrule,         /**< branching rule */
   int              priority            /**< new priority of the branching rule */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetBranchrulePriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPsetBranchruleMaxdepth", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPbranchruleSetMaxdepth(branchrule, scip->set, maxdepth);

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

   CHECK_OKAY( checkStage(scip, "SCIPincludeDisp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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

   CHECK_ABORT( checkStage(scip, "SCIPfindDisp", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindDisp(scip->set, name);
}

/** returns the array of currently available display columns */
DISP** SCIPgetDisps(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetDisps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->disps;
}

/** returns the number of currently available display columns */
int SCIPgetNDisps(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNDisps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->ndisps;
}

/** automatically selects display columns for being shown w.r.t. the display width parameter */
RETCODE SCIPautoselectDisps(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPselectDisps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPcreateDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPdialogCreate(dialog, dialogexec, dialogdesc, name, desc, issubmenu, dialogdata) );

   return SCIP_OKAY;
}

/** captures a dialog */
RETCODE SCIPcaptureDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          dialog              /**< dialog */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcaptureDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPdialogCapture(dialog);

   return SCIP_OKAY;
}

/** releases a dialog */
RETCODE SCIPreleaseDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG**         dialog              /**< pointer to the dialog */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPreleaseDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPdialogRelease(dialog) );

   return SCIP_OKAY;
}

/** makes given dialog the root dialog of SCIP's interactive user shell; captures dialog and releases former root dialog */
RETCODE SCIPsetRootDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          dialog              /**< dialog to be the root */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetRootDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPdialoghdlrSetRoot(scip->dialoghdlr, dialog) );

   return SCIP_OKAY;
}

/** returns the root dialog of SCIP's interactive user shell */
DIALOG* SCIPgetRootDialog(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRootDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return SCIPdialoghdlrGetRoot(scip->dialoghdlr);
}

/** adds a sub dialog to the given dialog as menu entry and captures it */
RETCODE SCIPaddDialogEntry(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          dialog,             /**< dialog to extend, or NULL for root dialog */
   DIALOG*          subdialog           /**< subdialog to add as menu entry in dialog */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddDialogEntry", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPstartInteraction", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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
   PROBDATA*        probdata            /**< user problem data set by the reader */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateProb", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   scip->stage = SCIP_STAGE_PROBLEM;
   
   CHECK_OKAY( SCIPstatCreate(&scip->stat, scip->mem->probmem, scip->set) );
   CHECK_OKAY( SCIPprobCreate(&scip->origprob, scip->mem->probmem, name, 
                  probdelorig, probtrans, probdeltrans, probdata, FALSE) );
   
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

   CHECK_OKAY( checkStage(scip, "SCIPreadProb", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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
         infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_NORMAL,
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
   CHECK_OKAY( checkStage(scip, "SCIPfreeProb", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   CHECK_OKAY( SCIPfreeSolve(scip) );
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
      CHECK_OKAY( SCIPprobFree(&scip->origprob, scip->mem->probmem, scip->set, scip->lp) );
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
   CHECK_ABORT( checkStage(scip, "SCIPgetProbData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobGetData(scip->origprob);

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
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
   CHECK_OKAY( checkStage(scip, "SCIPsetProbData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetData(scip->origprob, probdata);
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      SCIPprobSetData(scip->transprob, probdata);
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** sets objective sense of problem */
RETCODE SCIPsetObjsense(
   SCIP*            scip,               /**< SCIP data structure */
   OBJSENSE         objsense            /**< new objective sense */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetObjsense", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   if( objsense != SCIP_OBJSENSE_MAXIMIZE && objsense != SCIP_OBJSENSE_MINIMIZE )
   {
      errorMessage("Invalid objective sense\n");
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
   CHECK_OKAY( checkStage(scip, "SCIPsetObjlimit", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetExternObjlim(scip->origprob, objlimit);
      break;
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIPprobSetExternObjlim(scip->transprob, SCIPprobInternObjval(scip->origprob, scip->set, objlimit));
      objlimit = SCIPprobGetInternObjlim(scip->transprob, scip->set);
      if( SCIPsetIsLT(scip->set, objlimit, scip->primal->upperbound) )
      {
         CHECK_OKAY( SCIPprimalSetUpperbound(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
                        scip->tree, scip->lp, objlimit) );
      }
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
   CHECK_ABORT( checkStage(scip, "SCIPgetObjlimit", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPprobGetExternObjlim(scip->origprob);
}

/** informs SCIP, that the objective value is always integral in every feasible solution */
RETCODE SCIPsetObjIntegral(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetObjIntegral", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetObjIntegral(scip->origprob);
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
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
   CHECK_ABORT( checkStage(scip, "SCIPisObjIntegral", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobIsObjIntegral(scip->origprob);
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
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
   CHECK_OKAY( checkStage(scip, "SCIPaddVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

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
         errorMessage("Cannot add transformed variables to original problem\n");
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->origprob, scip->mem->probmem, scip->set, scip->lp, scip->branchcand, 
                     scip->eventfilter, scip->eventqueue, var) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
         {
            errorMessage("Cannot add original variables to transformed problem\n");
         }
         else
         {
            errorMessage("Cannot add fixed or aggregated variables to transformed problem\n");
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
   CHECK_OKAY( checkStage(scip, "SCIPaddPricedVar", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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
            errorMessage("Cannot add original variables to transformed problem\n");
         }
         else
         {
            errorMessage("Cannot add fixed or aggregated variables to transformed problem\n");
         }
         return SCIP_INVALIDDATA;
      }
      CHECK_OKAY( SCIPprobAddVar(scip->transprob, scip->mem->solvemem, scip->set, scip->lp, 
                     scip->branchcand, scip->eventfilter, scip->eventqueue, var) );
   }

   /* add variable to pricing storage */
   CHECK_OKAY( SCIPpricestoreAddVar(scip->pricestore, scip->mem->solvemem, scip->set, scip->lp, var, score,
                  (SCIPnodeGetDepth(scip->tree->actnode) == 0)) );
   
   return SCIP_OKAY;
}

/** gets variables of the problem along with the numbers of different variable types; data may become invalid after
 *  calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVar(), and SCIPmultiaggregateVar()
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
   CHECK_OKAY( checkStage(scip, "SCIPgetVarsData", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

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
 *  calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVar(), and SCIPmultiaggregateVar()
 */
VAR** SCIPgetVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->vars;

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
   CHECK_ABORT( checkStage(scip, "SCIPgetNVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nvars;

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
   CHECK_ABORT( checkStage(scip, "SCIPgetNBinVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nbinvars;

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
   CHECK_ABORT( checkStage(scip, "SCIPgetNIntVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nintvars;

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
   CHECK_ABORT( checkStage(scip, "SCIPgetNImplVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nimplvars;

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
   CHECK_ABORT( checkStage(scip, "SCIPgetNContVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->ncontvars;

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
   CHECK_OKAY( checkStage(scip, "SCIPgetOrigVarsData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_ABORT( checkStage(scip, "SCIPgetOrigVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->vars;
}

/** gets number of original problem variables */
int SCIPgetNOrigVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNOrigVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nvars;
}

/** gets number of binary original problem variables */
int SCIPgetNOrigBinVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNOrigBinVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nbinvars;
}

/** gets number of integer original problem variables */
int SCIPgetNOrigIntVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNOrigIntVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nintvars;
}

/** gets number of implicit integer original problem variables */
int SCIPgetNOrigImplVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNOrigImplVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nimplvars;
}

/** gets number of continuous original problem variables */
int SCIPgetNOrigContVars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNOrigContVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

   CHECK_ABORT( checkStage(scip, "SCIPfindVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobFindVar(scip->origprob, name);

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
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
   CHECK_ABORT( checkStage(scip, "SCIPallVarsInProb", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return (scip->set->nactivepricers == 0);
}

/** adds constraint to the problem; if constraint is only valid locally, it is added to the local subproblem of the
 *  active node (and all of its subnodes); otherwise it is added to the global problem;
 *  if a local constraint is added at the root node, it is automatically upgraded into a global constraint
 */
RETCODE SCIPaddCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddCons", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPprobAddCons(scip->origprob, scip->set, cons) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPprobAddCons(scip->transprob, scip->set, cons) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      if( SCIPconsIsGlobal(cons) )
      {
         CHECK_OKAY( SCIPprobAddCons(scip->transprob, scip->set, cons) );
      }
      else
      {
         CHECK_OKAY( SCIPnodeAddCons(scip->tree->actnode, scip->mem->solvemem, scip->set, scip->tree, cons) );
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

   CHECK_OKAY( checkStage(scip, "SCIPdelCons", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(cons->addconssetchg == NULL);
      CHECK_OKAY( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->origprob) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      assert(cons->addconssetchg == NULL);
      /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPconsDelete(cons, scip->mem->solvemem, scip->set, scip->transprob) );
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

   CHECK_ABORT( checkStage(scip, "SCIPfindCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobFindCons(scip->origprob, name);

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
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

   CHECK_OKAY( checkStage(scip, "SCIPaddConsNode", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( SCIPnodeGetDepth(node) == 0 )
   {
      assert(node == scip->tree->root);
      cons->local = FALSE;
      CHECK_OKAY( SCIPprobAddCons(scip->transprob, scip->set, cons) );
   }
   else
   {
      CHECK_OKAY( SCIPnodeAddCons(node, scip->mem->solvemem, scip->set, scip->tree, cons) );
   }

   return SCIP_OKAY;
}

/** adds constraint locally to the active node (and all of its subnodes), even if it is a global constraint;
 *  if a local constraint is added at the root node, it is automatically upgraded into a global constraint
 */
RETCODE SCIPaddConsLocal(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddConsLocal", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPaddConsNode(scip, scip->tree->actnode, cons) );
   
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

   CHECK_OKAY( checkStage(scip, "SCIPdisableConsNode", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPnodeDisableCons(node, scip->mem->solvemem, scip->set, scip->tree, cons) );
   
   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities at the active node (and all subnodes);
 *  if the method is called during problem modification or presolving, the constraint is globally deleted from the problem
 */
RETCODE SCIPdisableConsLocal(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint to disable */
   )
{
   assert(cons != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPdisableConsLocal", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(cons->addconssetchg == NULL);
      CHECK_OKAY( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->origprob) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      assert(cons->addconssetchg == NULL);
      CHECK_OKAY( SCIPconsDelete(cons, scip->mem->solvemem, scip->set, scip->transprob) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      assert(scip->tree->actnode != NULL);
      assert(SCIPnodeGetType(scip->tree->actnode) == SCIP_NODETYPE_ACTNODE);
      CHECK_OKAY( SCIPnodeDisableCons(scip->tree->actnode, scip->mem->solvemem, scip->set, scip->tree, cons) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets dual bound of active node */
Real SCIPgetLocalDualbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLocalDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return scip->tree->actnode != NULL
      ? SCIPprobExternObjval(scip->origprob, scip->set, 
         SCIPprobExternObjval(scip->transprob, scip->set, SCIPnodeGetLowerbound(scip->tree->actnode)))
      : SCIP_INVALID;
}

/** gets lower bound of active node in transformed problem */
Real SCIPgetLocalLowerbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLocalLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return scip->tree->actnode != NULL ? SCIPnodeGetLowerbound(scip->tree->actnode) : SCIP_INVALID;
}

/** gets dual bound of given node */
Real SCIPgetNodeDualbound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node                /**< node to get dual bound for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNodeDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->origprob, scip->set, 
      SCIPprobExternObjval(scip->transprob, scip->set, SCIPnodeGetLowerbound(node)));
}

/** gets lower bound of given node in transformed problem */
Real SCIPgetNodeLowerbound(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node                /**< node to get dual bound for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNodeLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPnodeGetLowerbound(node);
}

/** if given value is tighter (larger for minimization, smaller for maximization) than the active node's dual bound,
 *  sets the active node's dual bound to the new value
 */
RETCODE SCIPupdateLocalDualbound(
   SCIP*            scip,               /**< SCIP data structure */
   Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPupdateLocalDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(scip->tree->actnode, SCIPprobInternObjval(scip->transprob, scip->set, 
                               SCIPprobInternObjval(scip->origprob, scip->set, newbound)));

   return SCIP_OKAY;
}

/** if given value is larger than the active node's lower bound (in transformed problem), sets the active node's
 *  lower bound to the new value
 */
RETCODE SCIPupdateLocalLowerbound(
   SCIP*            scip,               /**< SCIP data structure */
   Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPupdateLocalLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(scip->tree->actnode, newbound);

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
   CHECK_OKAY( checkStage(scip, "SCIPupdateNodeDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(node, SCIPprobInternObjval(scip->transprob, scip->set, 
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
   CHECK_OKAY( checkStage(scip, "SCIPupdateNodeLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(node, newbound);

   return SCIP_OKAY;
}




/*
 * solve methods
 */

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

   assert(unbounded != NULL);
   assert(infeasible != NULL);

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

   maxnrounds = scip->set->maxpresolrounds;
   if( maxnrounds == -1 )
      maxnrounds = INT_MAX;

   abortfac = scip->set->presolabortfac;

   aborted = SCIPsolveIsStopped(scip->set, scip->stat);
   result = SCIP_DIDNOTRUN;

   infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_HIGH, "presolving:\n");

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
      for( i = 0; i < scip->set->npresols && result != SCIP_CUTOFF && result != SCIP_UNBOUNDED; ++i )
      {
         CHECK_OKAY( SCIPpresolExec(scip->set->presols[i], scip->set, nrounds, 
                        &nfixedvars, &naggrvars, &nchgvartypes, &nchgbds, &naddholes,
                        &ndelconss, &nupgdconss, &nchgcoefs, &nchgsides, &result) );
      }

      /* call presolve methods of constraint handlers */
      for( i = 0; i < scip->set->nconshdlrs && result != SCIP_CUTOFF && result != SCIP_UNBOUNDED; ++i )
      {
         CHECK_OKAY( SCIPconshdlrPresolve(scip->set->conshdlrs[i], scip->mem->solvemem, scip->set, scip->transprob, nrounds,
                        &nfixedvars, &naggrvars, &nchgvartypes, &nchgbds, &naddholes,
                        &ndelconss, &nupgdconss, &nchgcoefs, &nchgsides, &result) );
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
      aborted = aborted || (result == SCIP_CUTOFF) || (result == SCIP_UNBOUNDED);

      /* increase round number */
      nrounds++;

      if( !aborted )
      {
         /* print presolving statistics */
         infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_HIGH, "(round %d)", nrounds);
         infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_HIGH,
            " %d deleted vars, %d deleted constraints, %d tightened bounds, %d added holes, %d changed sides, %d changed coefficients\n",
            nfixedvars + naggrvars, ndelconss, nchgbds, naddholes, nchgsides, nchgcoefs);
      }
   }
   
   /* print presolving statistics */
   infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_NORMAL, "presolving (%d rounds):\n", nrounds);
   infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_NORMAL,
      " %d deleted vars, %d deleted constraints, %d tightened bounds, %d added holes, %d changed sides, %d changed coefficients\n",
      nfixedvars + naggrvars, ndelconss, nchgbds, naddholes, nchgsides, nchgcoefs);

   /* check for unboundness and infeasibility */
   *unbounded = (result == SCIP_UNBOUNDED);
   *infeasible = (result == SCIP_CUTOFF);

   return SCIP_OKAY;
}

/** presolves problem */
RETCODE SCIPpresolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   RESULT result;
   Bool unbounded;
   Bool infeasible;
   int h;

   CHECK_OKAY( checkStage(scip, "SCIPpresolve", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->stage == SCIP_STAGE_PRESOLVED )
      return SCIP_OKAY;

   /* check, if a node selector exists */
   if( SCIPsetGetNodesel(scip->set, scip->stat) == NULL )
   {
      errorMessage("no node selector available\n");
      return SCIP_PLUGINNOTFOUND;
   }
   
   /* inform original problem, that solving starts now */
   SCIPprobSolvingStarts(scip->origprob);

   /* switch stage to INITSOLVE */
   scip->stage = SCIP_STAGE_INITSOLVE;

   /* mark statistics before solving */
   SCIPstatMark(scip->stat);

   /* create VBC output file */
   CHECK_OKAY( SCIPvbcInit(scip->stat->vbc, scip->mem->solvemem, scip->set) );

   /* init solve data structures */
   CHECK_OKAY( SCIPeventfilterCreate(&scip->eventfilter, scip->mem->solvemem) );
   CHECK_OKAY( SCIPeventqueueCreate(&scip->eventqueue) );
   CHECK_OKAY( SCIPbranchcandCreate(&scip->branchcand) );
   CHECK_OKAY( SCIPlpCreate(&scip->lp, scip->set, SCIPprobGetName(scip->origprob)) );
   CHECK_OKAY( SCIPtreeCreate(&scip->tree, scip->mem->solvemem, scip->set, scip->stat, scip->lp, 
                  SCIPsetGetNodesel(scip->set, scip->stat)) );

   /* copy problem in solve memory */
   CHECK_OKAY( SCIPprobTransform(scip->origprob, scip->mem->solvemem, scip->set, scip->stat, scip->lp, 
                  scip->branchcand, scip->eventfilter, scip->eventqueue, &scip->transprob) );

   /* init callback methods */
   CHECK_OKAY( SCIPsetInitCallbacks(scip->set) );

   /* switch stage to PRESOLVING */
   scip->stage = SCIP_STAGE_PRESOLVING;

   /* start presolving timer */
   SCIPclockStart(scip->stat->presolvingtime, scip->set);

   /* presolve problem */
   CHECK_OKAY( presolve(scip, &unbounded, &infeasible) );

   /* inform constraint handlers that presolving is finished and the branch and bound process starts now */
   for( h = 0; h < scip->set->nconshdlrs; ++h )
   {
      CHECK_OKAY( SCIPconshdlrSolstart(scip->set->conshdlrs[h], scip, &result) );
      unbounded = unbounded || (result == SCIP_UNBOUNDED);
      infeasible = infeasible || (result == SCIP_INFEASIBLE);
   }

   /* check, whether the objective value is always integral */
   SCIPprobCheckObjIntegral(scip->transprob, scip->set);

   /* stop presolving time */
   SCIPclockStop(scip->stat->presolvingtime, scip->set);

   /* switch stage to PRESOLVED */
   scip->stage = SCIP_STAGE_PRESOLVED;

   /* init solution process data structures */
   CHECK_OKAY( SCIPpricestoreCreate(&scip->pricestore) );
   CHECK_OKAY( SCIPsepastoreCreate(&scip->sepastore) );
   CHECK_OKAY( SCIPcutpoolCreate(&scip->cutpool, scip->mem->solvemem, scip->set->cutagelimit) );
   CHECK_OKAY( SCIPconflictCreate(&scip->conflict, scip->set) );
   CHECK_OKAY( SCIPlpconflictCreate(&scip->lpconflict) );
   CHECK_OKAY( SCIPpseudoconflictCreate(&scip->pseudoconflict) );
   CHECK_OKAY( SCIPprimalCreate(&scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp) );

   /* inform problem, that presolving is finished and the branch and bound process starts now */
   SCIPprobSolvingStarts(scip->transprob);

   if( infeasible )
   {
      /* print solution message */
      infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_NORMAL, "presolving detected infeasibility.\n");

      /* switch stage to SOLVED */
      scip->stage = SCIP_STAGE_SOLVED;
   }
   else if( unbounded )
   {
      infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_NORMAL, "presolving detected unboundness.\n");

      /* switch stage to SOLVED */
      scip->stage = SCIP_STAGE_SOLVED;
   }
   else
   {
      /* print presolved problem statistics */
      infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_NORMAL,
         "presolved problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
         scip->transprob->nvars, scip->transprob->nbinvars, scip->transprob->nintvars, scip->transprob->nimplvars,
         scip->transprob->ncontvars, scip->transprob->nconss);

      for( h = 0; h < scip->set->nconshdlrs; ++h )
      {
         int nconss;

         nconss = SCIPconshdlrGetNConss(scip->set->conshdlrs[h]);
         if( nconss > 0 )
         {
            infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_HIGH,
               " %5d constraints of type <%s>\n", nconss, SCIPconshdlrGetName(scip->set->conshdlrs[h]));
         }
      }

      if( SCIPprobIsObjIntegral(scip->transprob) )
      {
         infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_HIGH, "objective value is always integral\n");
      }
   }

   /* display timing statistics */
   infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_HIGH,
      "Presolving Time: %.2f\n", SCIPclockGetTime(scip->stat->presolvingtime));

   /* if all variables are known, calculate a trivial primal bound by setting all variables to their
    * worst bound
    */
   if( scip->set->nactivepricers == 0 )
   {
      VAR* var;
      Real objbound;
      Real bd;
      int v;
      
      objbound = 0.0;
      for( v = 0; v < scip->transprob->nvars && !SCIPsetIsInfinity(scip->set, objbound); ++v )
      {
         var = scip->transprob->vars[v];
         bd = SCIPvarGetWorstBound(var);
         if( SCIPsetIsInfinity(scip->set, ABS(bd)) )
            objbound = scip->set->infinity;
         else
            objbound += SCIPvarGetObj(var) * bd;
      }
      if( !SCIPsetIsInfinity(scip->set, objbound) && SCIPsetIsLT(scip->set, objbound, scip->primal->upperbound) )
      {
         /* add 1.0 to primal bound, such that solution with worst bound may be found */
         CHECK_OKAY( SCIPprimalSetUpperbound(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
                        scip->tree, scip->lp, objbound + 1.0) );
      }
   }

   return SCIP_OKAY;
}

/** presolves and solves problem */
RETCODE SCIPsolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsolve", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   /* check, if a node selector exists */
   if( SCIPsetGetNodesel(scip->set, scip->stat) == NULL )
   {
      errorMessage("no node selector available\n");
      return SCIP_PLUGINNOTFOUND;
   }
   
   /* capture the CTRL-C interrupt */
   if( scip->set->catchctrlc )
      SCIPinterruptCapture(scip->interrupt);

   /* start solving timer */
   SCIPclockStart(scip->stat->solvingtime, scip->set);

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      /* presolve problem and set up data structures */
      CHECK_OKAY( SCIPpresolve(scip) );
      if( scip->stage == SCIP_STAGE_SOLVED )
         return SCIP_OKAY;
      assert(scip->stage == SCIP_STAGE_PRESOLVED);
      infoMessage(scip->set->verblevel, SCIP_VERBLEVEL_NORMAL, "\n");

      /*lint -fallthrough*/

   case SCIP_STAGE_PRESOLVED:
      /* switch stage to SOLVING */
      scip->stage = SCIP_STAGE_SOLVING;
      
      /*lint -fallthrough*/

   case SCIP_STAGE_SOLVING:
      /* reset display */
      SCIPstatResetDisplay(scip->stat);

      /* continue solution process */
      CHECK_OKAY( SCIPsolveCIP(scip->mem->solvemem, scip->set, scip->stat, scip->mem, scip->transprob, scip->tree, 
                     scip->lp, scip->pricestore, scip->sepastore, scip->branchcand, scip->cutpool, 
                     scip->conflict, scip->lpconflict, scip->pseudoconflict, scip->primal,
                     scip->eventfilter, scip->eventqueue) );

      /* detect, whether problem is solved */
      if( SCIPtreeGetNNodes(scip->tree) == 0 && scip->tree->actnode == NULL )
      {
         /* tree is empty, and no active node exists -> problem is solved */
         scip->stage = SCIP_STAGE_SOLVED;
      }

      /* display most relevant statistics */
      if( scip->set->verblevel >= SCIP_VERBLEVEL_NORMAL )
      {
         printf("\n");
         printf("SCIP Status        : ");
         CHECK_OKAY( SCIPprintStatus(scip, NULL) );
         printf("\n");
         printf("Solving Time (sec) : %.2f\n", SCIPclockGetTime(scip->stat->solvingtime));
         printf("Solving Nodes      : %lld\n", scip->stat->nnodes);
         printf("Primal Bound       : %+.14e (%lld solutions)\n", SCIPgetPrimalbound(scip), scip->primal->nsolsfound);
         printf("Dual Bound         : %+.14e\n", SCIPgetDualbound(scip));
         printf("Gap                : ");
         if( SCIPsetIsInfinity(scip->set, SCIPgetGap(scip)) )
            printf("infinite\n");
         else
            printf("%.2f %%\n", 100.0*SCIPgetGap(scip));
      }
      break;

   case SCIP_STAGE_SOLVED:
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   /* stop solving timer */
   SCIPclockStop(scip->stat->solvingtime, scip->set);

   /* release the CTRL-C interrupt */
   if( scip->set->catchctrlc )
      SCIPinterruptRelease(scip->interrupt);

   return SCIP_OKAY;
}

/** frees all solution process data, only original problem is kept */
RETCODE SCIPfreeSolve(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeSolve", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      /* switch stage to FREESOLVE */
      scip->stage = SCIP_STAGE_FREESOLVE;

      /* deactivate the active node */
      CHECK_OKAY( SCIPnodeActivate(NULL, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, 
                     scip->branchcand, scip->eventqueue, scip->primal->cutoffbound) );

      /* exit callback methods */
      CHECK_OKAY( SCIPsetExitCallbacks(scip->set) );

      /* clear the LP */
      CHECK_OKAY( SCIPlpClear(scip->lp, scip->mem->solvemem, scip->set) );

      /* free solution process data structures */
      CHECK_OKAY( SCIPprimalFree(&scip->primal, scip->mem->solvemem) );
      CHECK_OKAY( SCIPpseudoconflictFree(&scip->pseudoconflict) );
      CHECK_OKAY( SCIPlpconflictFree(&scip->lpconflict) );
      CHECK_OKAY( SCIPconflictFree(&scip->conflict) );
      CHECK_OKAY( SCIPcutpoolFree(&scip->cutpool, scip->mem->solvemem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPsepastoreFree(&scip->sepastore) );
      CHECK_OKAY( SCIPpricestoreFree(&scip->pricestore) );
      CHECK_OKAY( SCIPtreeFree(&scip->tree, scip->mem->solvemem, scip->set, scip->lp) );
      CHECK_OKAY( SCIPprobFree(&scip->transprob, scip->mem->solvemem, scip->set, scip->lp) );

      /* free solve data structures */
      CHECK_OKAY( SCIPlpFree(&scip->lp, scip->mem->solvemem, scip->set) );
      CHECK_OKAY( SCIPbranchcandFree(&scip->branchcand) );
      CHECK_OKAY( SCIPeventfilterFree(&scip->eventfilter, scip->mem->solvemem, scip->set) );
      CHECK_OKAY( SCIPeventqueueFree(&scip->eventqueue) );

      /* close VBC output file */
      SCIPvbcExit(scip->stat->vbc, scip->set);

      /* free the solve block memory */
#ifndef NDEBUG
      blockMemoryCheckEmpty(scip->mem->solvemem);
#endif
      freeAllBlockMemoryNull(scip->mem->solvemem);

      /* reset statistics to the point before solving started */
      SCIPstatReset(scip->stat);

      /* switch stage to PROBLEM */
      scip->stage = SCIP_STAGE_PROBLEM;

      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}




/*
 * variable methods
 */

/** create and capture problem variable; if variable is of integral type, fractional bounds are automatically rounded */
RETCODE SCIPcreateVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            var,                /**< pointer to variable object */
   const char*      name,               /**< name of variable, or NULL for automatic name creation */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             initial,            /**< should var's column be present in the initial root LP? */
   Bool             removeable          /**< is var's column removeable from the LP (due to aging or cleanup)? */
   )
{
   assert(var != NULL);
   assert(lb <= ub);

   CHECK_OKAY( checkStage(scip, "SCIPcreateVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarCreateOriginal(var, scip->mem->probmem, scip->set, scip->stat, 
                     name, lb, ub, obj, vartype, initial, removeable) );
      break;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPvarCreateTransformed(var, scip->mem->solvemem, scip->set, scip->stat,
                     name, lb, ub, obj, vartype, initial, removeable) );
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
   CHECK_OKAY( checkStage(scip, "SCIPcaptureVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPreleaseVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarRelease(var, scip->mem->probmem, scip->set, scip->lp) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      if( !SCIPvarIsTransformed(*var) )
      {
         errorMessage("cannot release original variables while solving the problem\n");
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

   CHECK_OKAY( checkStage(scip, "SCIPtransformVar", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPtransformVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPgetTransformedVar", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPgetTransformedVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPgetNegatedVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

/** gets solution value for variable in active node */
Real SCIPgetVarSol(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get solution value for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarSol", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetSol(var, scip->tree->actnodehaslp);
}

/** gets solution values of multiple variables in active node */
RETCODE SCIPgetVarSols(
   SCIP*            scip,               /**< SCIP data structure */
   int              nvars,              /**< number of variables to get solution value for */
   VAR**            vars,               /**< array with variables to get value for */
   Real*            vals                /**< array to store solution values of variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPgetVarSols", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   for( v = 0; v < nvars; ++v )
      vals[v] = SCIPvarGetSol(vars[v], scip->tree->actnodehaslp);

   return SCIP_OKAY;
}

/** gets strong branching information on COLUMN variable */
RETCODE SCIPgetVarStrongbranch(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to get strong branching values for */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching column down */
   Real*            up,                 /**< stores dual bound after branching column up */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetVarStrongbranch", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      errorMessage("cannot get strong branching information on non-COLUMN variable\n");
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPcolGetStrongbranch(SCIPvarGetCol(var), scip->set, scip->stat, scip->lp, itlim, down, up, lperror) );

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
   Real*            solval              /**< stores LP solution value of variable at the last strong branching call */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetVarStrongbranchLast", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      errorMessage("cannot get strong branching information on non-COLUMN variable\n");
      return SCIP_INVALIDDATA;
   }

   SCIPcolGetStrongbranchLast(SCIPvarGetCol(var), down, up, solval);

   return SCIP_OKAY;
}

/** gets number of the last node where strong branching was used on the given COLUMN variable,
 *  or -1 if strong branching was never applied to the variable
 */
Longint SCIPgetVarStrongbranchNode(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get last strong branching node for */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetVarStrongbranchNode", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPchgVarObj", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      CHECK_OKAY( SCIPvarChgObj(var, scip->mem->probmem, scip->set, scip->lp, scip->branchcand,
                     scip->eventqueue, newobj) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarChgObj(var, scip->mem->solvemem, scip->set, scip->lp, scip->branchcand,
                     scip->eventqueue, newobj) );
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
   CHECK_OKAY( checkStage(scip, "SCIPaddVarObj", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      CHECK_OKAY( SCIPvarAddObj(var, scip->mem->probmem, scip->set, scip->origprob, scip->lp,
                     scip->branchcand, scip->eventqueue, addobj) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarAddObj(var, scip->mem->solvemem, scip->set, scip->transprob, scip->lp,
                     scip->branchcand, scip->eventqueue, addobj) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/   
}

/** depending on SCIP's stage, changes lower bound of variable in the problem, in preprocessing, or in active node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 */
RETCODE SCIPchgVarLb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarLb", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarSetLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
                     scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarSetLbGlobal(var, scip->set, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarSetLbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
                     scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarSetLbGlobal(var, scip->set, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundchg(scip->tree->actnode, scip->mem->solvemem, scip->set, scip->stat, scip->tree,
                     scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER) );
      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** depending on SCIP's stage, changes upper bound of variable in the problem, in preprocessing, or in active node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 */
RETCODE SCIPchgVarUb(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarUb", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPvarSetUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
                     scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarSetUbGlobal(var, scip->set, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarSetUbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
                     scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarSetUbGlobal(var, scip->set, newbound) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundchg(scip->tree->actnode, scip->mem->solvemem, scip->set, scip->stat, scip->tree,
                     scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER) );
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
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarLbNode", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   SCIPvarAdjustLb(var, scip->set, &newbound);

   CHECK_OKAY( SCIPnodeAddBoundchg(node, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, 
                  scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER) );
   
   return SCIP_OKAY;
}

/** changes upper bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 */
RETCODE SCIPchgVarUbNode(
   SCIP*            scip,               /**< SCIP data structure */
   NODE*            node,               /**< node to change bound at, or NULL for active node */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarUbNode", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   SCIPvarAdjustUb(var, scip->set, &newbound);

   CHECK_OKAY( SCIPnodeAddBoundchg(node, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, 
                  scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER) );
   
   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the active node, if the new bound is tighter than the
 *  current bound; if possible, adjusts bound to integral value; doesn't store any inference information in the
 *  bound change, such that in conflict analysis, this change is treated like a branching decision
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

   CHECK_OKAY( checkStage(scip, "SCIPtightenVarLb", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   *infeasible = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   
   if( SCIPsetIsFeasLE(scip->set, newbound, lb) )
   {
      if( tightened != NULL )
         *tightened = FALSE;
      return SCIP_OKAY;
   }

   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   switch( scip->stage )
   {
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarSetLbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
                     scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarSetLbGlobal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundchg(scip->tree->actnode, scip->mem->solvemem, scip->set, scip->stat, scip->tree,
                     scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER) );
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the active node, if the new bound is tighter than the
 *  current bound; if possible, adjusts bound to integral value; doesn't store any inference information in the
 *  bound change, such that in conflict analysis, this change is treated like a branching decision
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

   CHECK_OKAY( checkStage(scip, "SCIPtightenVarUb", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   *infeasible = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   
   if( SCIPsetIsFeasGE(scip->set, newbound, ub) )
   {
      if( tightened != NULL )
         *tightened = FALSE;
      return SCIP_OKAY;
   }

   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   switch( scip->stage )
   {
   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarSetUbLocal(var, scip->mem->solvemem, scip->set, scip->stat, scip->lp,
                     scip->branchcand, scip->eventqueue, newbound) );
      CHECK_OKAY( SCIPvarSetUbGlobal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPnodeAddBoundchg(scip->tree->actnode, scip->mem->solvemem, scip->set, scip->stat, scip->tree,
                     scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER) );
      break;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** fixes binary variable to given value; in problem creation or preprocessing stage, the variable is converted
 *  into a fixed variable, and the given inference constraint is ignored; in solving stage, the variable is fixed
 *  locally at the given node, and the given inference constraint is stored, such that the conflict analysis is
 *  able to find out the reason for the deduction of the variable fixing
 */
RETCODE SCIPinferBinVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< binary variable, that is deduced to a fixed value */
   Bool             fixedval,           /**< value to fix binary variable to */
   CONS*            infercons,          /**< constraint that deduced the fixing */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   Real lb;
   Real ub;

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(fixedval == TRUE || fixedval == FALSE);
   assert(infeasible != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPinferBinVar", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

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
      CHECK_OKAY( SCIPvarFix(var, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp,
                     scip->branchcand, scip->eventqueue, (Real)fixedval, infeasible) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      if( fixedval == TRUE )
      {
         CHECK_OKAY( SCIPnodeAddBoundinfer(scip->tree->actnode, scip->mem->solvemem, scip->set, scip->stat, scip->tree,
                        scip->lp, scip->branchcand, scip->eventqueue, var, 1.0, SCIP_BOUNDTYPE_LOWER, 
                        infercons, inferinfo) );
      }
      else
      {
         CHECK_OKAY( SCIPnodeAddBoundinfer(scip->tree->actnode, scip->mem->solvemem, scip->set, scip->stat, scip->tree,
                        scip->lp, scip->branchcand, scip->eventqueue, var, 0.0, SCIP_BOUNDTYPE_UPPER, 
                        infercons, inferinfo) );
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

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
RETCODE SCIPchgVarBranchFactor(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             branchfactor        /**< factor to weigh variable's branching score with */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPscaleVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIPvarChgBranchFactor(var, scip->set, scale * SCIPvarGetBranchFactor(var));

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
   CHECK_OKAY( checkStage(scip, "SCIPchgVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIPvarChgBranchPriority(var, scip->set, branchpriority);

   return SCIP_OKAY;
}

/** changes the branch priority of the variable to the given value, if it is larger than the current priority */
RETCODE SCIPupdateVarBranchPriority(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   int              branchpriority      /**< new branch priority of the variable, if it is larger than current priority */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPupdateVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   if( branchpriority > SCIPvarGetBranchPriority(var) )
      SCIPvarChgBranchPriority(var, scip->set, branchpriority);

   return SCIP_OKAY;
}

/** adds the given value to the branch priority of the variable */
RETCODE SCIPaddVarBranchPriority(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   int              addpriority         /**< value to add to the branch priority of the variable */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIPvarChgBranchPriority(var, scip->set, addpriority + SCIPvarGetBranchPriority(var));

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

   CHECK_OKAY( checkStage(scip, "SCIPchgVarType", FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

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

   case SCIP_STAGE_INITSOLVE:
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
   Real             fixedval,           /**< value to fix variable at */
   Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   Bool*            fixed               /**< pointer to store whether the fixing was performed (variable was unfixed) */
   )
{
   assert(var != NULL);
   assert(infeasible != NULL);
   assert(fixed != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfixVar", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

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
      if( SCIPsetIsGT(scip->set, fixedval, SCIPvarGetLbLocal(var)) )
      {
         CHECK_OKAY( SCIPchgVarLb(scip, var, fixedval) );
         *fixed = TRUE;
      }
      if( SCIPsetIsLT(scip->set, fixedval, SCIPvarGetUbLocal(var)) )
      {
         CHECK_OKAY( SCIPchgVarUb(scip, var, fixedval) );
         *fixed = TRUE;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      CHECK_OKAY( SCIPvarFix(var, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp,
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
   success = SCIPrealToRational(scalarx, scip->set->epsilon, MAXDNOM, &scalarxn, &scalarxd);
   if( success )
      success = SCIPrealToRational(scalary, scip->set->epsilon, MAXDNOM, &scalaryn, &scalaryd);
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
                     scip->lp, scip->branchcand, scip->eventqueue, vary, (Real)(-b/a), (Real)(c/a), infeasible) );
      *aggregated = TRUE;
      return SCIP_OKAY;
   }
   if( b == 1 || b == -1 )
   {
      /* aggregate y = - a/b*x + c/b */
      /*lint --e{653}*/
      CHECK_OKAY( SCIPvarAggregate(vary, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
                     scip->lp, scip->branchcand, scip->eventqueue, varx, (Real)(-a/b), (Real)(c/b), infeasible) );
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
                  SCIPvarIsInitial(varx) || SCIPvarIsInitial(vary),
                  SCIPvarIsRemoveable(varx) && SCIPvarIsRemoveable(vary) ) );
   CHECK_OKAY( SCIPprobAddVar(scip->transprob, scip->mem->solvemem, scip->set, scip->lp, 
                  scip->branchcand, scip->eventfilter, scip->eventqueue, aggvar) );
   CHECK_OKAY( SCIPvarAggregate(varx, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
                  scip->lp, scip->branchcand, scip->eventqueue, aggvar, (Real)(-b), (Real)xsol, infeasible) );
   if( !(*infeasible) )
   {
      CHECK_OKAY( SCIPvarAggregate(vary, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
                     scip->lp, scip->branchcand, scip->eventqueue, aggvar, (Real)a, (Real)ysol, infeasible) );
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
                        scip->lp, scip->branchcand, scip->eventqueue, vary, scalar, constant, infeasible) );
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

   CHECK_OKAY( checkStage(scip, "SCIPaggregateVars", FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   *redundant = FALSE;
   *aggregated = FALSE;

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
      CHECK_OKAY( SCIPvarFix(vary, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp,
                     scip->branchcand, scip->eventqueue, rhs/scalary, infeasible) );
      *aggregated = TRUE;
      *redundant = TRUE;
   }
   else if( vary == NULL )
   {
      assert(SCIPsetIsZero(scip->set, scalary));
      assert(!SCIPsetIsZero(scip->set, scalarx));
      
      /* variable y was resolved to fixed variable: variable x can be fixed to c'/a' */
      CHECK_OKAY( SCIPvarFix(varx, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp,
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
         CHECK_OKAY( SCIPvarFix(varx, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp,
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
 *  SCIPgetVars() and SCIPgetVarsData()
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
   CHECK_OKAY( checkStage(scip, "SCIPmultiaggregateVar", FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarMultiaggregate(var, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp,
                  scip->branchcand, scip->eventqueue, naggvars, aggvars, scalars, constant, infeasible) );

   return SCIP_OKAY;
}




/*
 * conflict analysis methods
 */

/** initializes the conflict analysis by clearing the conflict variable candidate queue */
RETCODE SCIPinitConflictAnalysis(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPinitConflictAnalysis", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconflictInit(scip->conflict) );

   return SCIP_OKAY;
}

/** adds currently fixed binary variable to the conflict analysis' candidate storage; this method should be called in
 *  one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictVar() should be called for each variable,
 *      whose current assignment lead to the conflict (i.e. the infeasibility of a globally valid constraint).
 *   2. In the conflict variable resolution method of a constraint handler, SCIPaddConflictVar() should be called
 *      for each variable, whose current assignment lead to the deduction of the given conflict variable.
 */
RETCODE SCIPaddConflictVar(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< conflict variable to add to conflict candidate queue */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPaddConflictVar", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconflictAddVar(scip->conflict, scip->mem->solvemem, scip->set, scip->stat, var) );

   return SCIP_OKAY;
}

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set; the conflict analysis
 *  should only be called if a globally valid constraint was violated -- otherwise, the resulting conflict
 *  constraint wouldn't be globally valid
 */
RETCODE SCIPanalyzeConflict(
   SCIP*            scip,               /**< SCIP data structure */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPanalyzeConflict", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconflictAnalyze(scip->conflict, scip->set, scip->transprob, scip->tree, success) );
   
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

   CHECK_OKAY( checkStage(scip, "SCIPcreateCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPconsCreate(cons, scip->mem->probmem, name, conshdlr, consdata, 
                     initial, separate, enforce, check, propagate, local, modifiable, removeable, TRUE) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
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
   CHECK_OKAY( checkStage(scip, "SCIPcaptureCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPreleaseCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPconsRelease(cons, scip->mem->probmem, scip->set) );
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      if( SCIPconsIsOriginal(*cons) && (*cons)->nuses == 1 )
      {
         errorMessage("cannot release last use of original constraint while solving the problem\n");
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

   CHECK_OKAY( checkStage(scip, "SCIPtransformCons", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPtransformConss", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPgetTransformedCons", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPgetTransformedConss", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPaddConsAge", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsAddAge(cons, scip->mem->solvemem, scip->set, scip->transprob, deltaage) );

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
   CHECK_OKAY( checkStage(scip, "SCIPincConsAge", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsIncAge(cons, scip->mem->solvemem, scip->set, scip->transprob) );

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
   CHECK_OKAY( checkStage(scip, "SCIPresetConsAge", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsResetAge(cons, scip->set) );

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
   CHECK_OKAY( checkStage(scip, "SCIPlockConsVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPunlockConsVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPcheckCons", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPconsCheck(cons, scip->set, sol, checkintegrality, checklprows, result) );

   return SCIP_OKAY;
}

/** marks the constraint to be essential for feasibility */
RETCODE SCIPsetConsChecked(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetConsChecked", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPconsSetChecked(cons, scip->set) );

   return SCIP_OKAY;
}




/*
 * LP methods
 */

/** checks, whether the LP was solved in the active node */
Bool SCIPhasActNodeLP(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPhasActNodeLP", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return scip->tree->actnodehaslp;
}

/** gets solution status of current LP */
LPSOLSTAT SCIPgetLPSolstat(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLPSolstat", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPlpGetSolstat(scip->lp);
}

/** gets objective value of current LP */
Real SCIPgetLPObjval(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLPObjval", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPlpGetObjval(scip->lp, scip->set);
}

/** gets pseudo objective value of the current LP */
Real SCIPgetPseudoObjval(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPseudoObjval", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   return SCIPlpGetPseudoObjval(scip->lp, scip->set);
}

/** gets current LP columns along with the current number of LP columns */
RETCODE SCIPgetLPColsData(
   SCIP*            scip,               /**< SCIP data structure */
   COL***           cols,               /**< pointer to store the array of LP columns, or NULL */
   int*             ncols               /**< pointer to store the number of LP columns, or NULL */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPColsData", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );   

   if( scip->tree->actnodehaslp )
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
   CHECK_ABORT( checkStage(scip, "SCIPgetLPCols", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );   

   if( scip->tree->actnodehaslp )
      return SCIPlpGetCols(scip->lp);
   else
      return NULL;
}

/** gets current number of LP columns */
int SCIPgetNLPCols(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNLPCols", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );   

   if( scip->tree->actnodehaslp )
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
   CHECK_OKAY( checkStage(scip, "SCIPgetLPRowsData", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );   
   
   if( scip->tree->actnodehaslp )
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
   CHECK_ABORT( checkStage(scip, "SCIPgetLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );   
   
   if( scip->tree->actnodehaslp )
      return SCIPlpGetRows(scip->lp);
   else
      return NULL;
}

/** gets current number of LP rows */
int SCIPgetNLPRows(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );   

   if( scip->tree->actnodehaslp )
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
   CHECK_ABORT( checkStage(scip, "SCIPallColsInLP", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp);
}

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1 */
RETCODE SCIPgetLPBasisInd(
   SCIP*            scip,               /**< SCIP data structure */
   int*             basisind            /**< pointer to store the basis indices */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLPBasisInd", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPgetLPBInvRow", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPgetLPBInvARow", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPsumLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPlpSumRows(scip->lp, scip->set, scip->transprob->nvars, weights, sumcoef, sumlhs, sumrhs) );

   return SCIP_OKAY;
}

/* calculates a MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
RETCODE SCIPcalcMIR(
   SCIP*            scip,               /**< SCIP data structure */
   Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   Bool*            success             /**< pointer to store whether the returned coefficients are a valid MIR cut */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcalcMIR", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPlpCalcMIR(scip->lp, scip->set, scip->stat, scip->transprob->nvars, scip->transprob->vars,
                  minfrac, weights, mircoef, mirrhs, success) );

   return SCIP_OKAY;
}

/** writes current LP to a file */
RETCODE SCIPwriteLP(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      fname               /**< file name */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPwriteLP", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   CHECK_OKAY( SCIPlpWrite(scip->lp, fname) );

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
   CHECK_OKAY( checkStage(scip, "SCIPstartDive", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( scip->lp->diving )
   {
      errorMessage("already in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   if( !scip->tree->actnodehaslp )
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
   CHECK_OKAY( checkStage(scip, "SCIPendDive", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( !scip->lp->diving )
   {
      errorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* reset the probably changed LP's cutoff bound */
   CHECK_OKAY( SCIPlpSetCutoffbound(scip->lp, scip->primal->cutoffbound) );

   /* unmark the diving flag in the LP and reset all variables' objective and bound values */
   CHECK_OKAY( SCIPlpEndDive(scip->lp, scip->mem->solvemem, scip->set, scip->stat, scip->transprob,
                  scip->transprob->vars, scip->transprob->nvars) );

   /* if a new best solution was created, the cutoff of the tree was delayed due to diving;
    * the cutoff has to be done now.
    */
   if( scip->tree->cutoffdelayed )
   {
      CHECK_OKAY( SCIPtreeCutoff(scip->tree, scip->mem->solvemem, scip->set, scip->lp, scip->primal->cutoffbound) );
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
   CHECK_OKAY( checkStage(scip, "SCIPchgVarObjDive", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( !scip->lp->diving )
   {
      errorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPvarChgObjDive(var, scip->set, scip->lp, newobj) );

   /* invalidate the LP's cutoff bound, since this has nothing to do with the current objective value anymore;
    * the cutoff bound is reset in SCIPendDive()
    */
   CHECK_OKAY( SCIPlpSetCutoffbound(scip->lp, scip->set->infinity) );

   return SCIP_OKAY;
}

/** changes variable's lower bound in current dive */
RETCODE SCIPchgVarLbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to change the bound for */
   Real             newbound            /**< new value for bound */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPchgVarLbDive", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   if( !scip->lp->diving )
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
   CHECK_OKAY( checkStage(scip, "SCIPchgVarUbDive", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   if( !scip->lp->diving )
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
   CHECK_ABORT( checkStage(scip, "SCIPgetVarObjDive", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   if( !scip->lp->diving )
   {
      errorMessage("not in diving mode\n");
      abort();
   }

   return SCIPvarGetObjDive(var, scip->set);
}

/** gets variable's lower bound in current dive */
Real SCIPgetVarLbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get the bound for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarLbDive", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   if( !scip->lp->diving )
   {
      errorMessage("not in diving mode\n");
      abort();
   }

   return SCIPvarGetLbDive(var, scip->set);
}

/** gets variable's upper bound in current dive */
Real SCIPgetVarUbDive(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var                 /**< variable to get the bound for */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarUbDive", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   if( !scip->lp->diving )
   {
      errorMessage("not in diving mode\n");
      abort();
   }

   return SCIPvarGetUbDive(var, scip->set);
}

/** solves the LP of the current dive */
RETCODE SCIPsolveDiveLP(
   SCIP*            scip,               /**< SCIP data structure */
   int              itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   Bool*            lperror             /**< pointer to store whether an unresolved LP error occured */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsolveDiveLP", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( !scip->lp->diving )
   {
      errorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPlpSolveAndEval(scip->lp, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, 
                  itlim, FALSE, lperror) );

   return SCIP_OKAY;
}

/** returns the number of the node where the last LP diving was applied */
Longint SCIPgetLastDivenode(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetLastDivenode", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   return scip->stat->lastdivenode;
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
   COL**            col,                /**< array with columns of row entries */
   Real*            val,                /**< array with coefficients of row entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             local,              /**< is row only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateRow", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwCreate(row, scip->mem->solvemem, scip->set, scip->stat,
                  name, len, col, val, lhs, rhs, local, modifiable, removeable) );

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
   CHECK_OKAY( checkStage(scip, "SCIPcreateEmptyRow", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPcaptureRow", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   SCIProwCapture(row);

   return SCIP_OKAY;
}

/** decreases usage counter of LP row, and frees memory if necessary */
RETCODE SCIPreleaseRow(
   SCIP*            scip,               /**< SCIP data structure */
   ROW**            row                 /**< pointer to LP row */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPreleaseRow", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPchgRowLhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPchgRowRhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwChgRhs(row, scip->set, scip->lp, rhs) );

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
   CHECK_OKAY( checkStage(scip, "SCIPaddVarToRow", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPvarAddToRow(var, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->lp, row, val) );
   
   return SCIP_OKAY;
}

/** resolves variables to columns and adds them with the coefficients to the row;
 *  if you want to add more than one variable to a row, for performance reasons this method is highly
 *  preferable to many single calls to SCIPaddVarToRow()
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

   assert(vars != NULL);
   assert(vals != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddVarsToRow", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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
 *  if you want to add more than one variable to a row, for performance reasons this method is highly
 *  preferable to many single calls to SCIPaddVarToRow()
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

   assert(vars != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddVarsToRow", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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

/** tries to find a rational representation of the row and multiplies coefficients with common denominator */
RETCODE SCIPmakeRowRational(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row,                /**< LP row */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   Real             maxscale,           /**< maximal value to scale row with */
   Bool*            success             /**< stores whether row could be made rational */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPmakeRowRational", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIProwMakeRational(row, scip->set, scip->stat, scip->lp, maxdnom, maxscale, success) );

   return SCIP_OKAY;
}

/** returns the minimal activity of a row w.r.t. the column's bounds */
Real SCIPgetRowMinActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowMinActivity", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   return SCIProwGetMinActivity(row, scip->set, scip->stat);
}

/** returns the maximal activity of a row w.r.t. the column's bounds */
Real SCIPgetRowMaxActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowMaxActivity", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   return SCIProwGetMaxActivity(row, scip->set, scip->stat);
}

/** returns the activity of a row in the last LP solution */
Real SCIPgetRowLPActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   return SCIProwGetLPActivity(row, scip->stat, scip->lp);
}

/** returns the feasibility of a row in the last LP solution: negative value means infeasibility */
Real SCIPgetRowLPFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowLPFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
   
   return SCIProwGetLPFeasibility(row, scip->stat, scip->lp);
}

/** returns the activity of a row for the current pseudo solution */
Real SCIPgetRowPseudoActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIProwGetPseudoActivity(row, scip->stat);
}

/** returns the feasibility of a row for the current pseudo solution: negative value means infeasibility */
Real SCIPgetRowPseudoFeasibility(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowPseudoFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIProwGetPseudoFeasibility(row, scip->stat);
}

/** returns the activity of a row in the last LP or pseudo solution */
Real SCIPgetRowActivity(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< LP row */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( scip->tree->actnodehaslp )
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
   CHECK_ABORT( checkStage(scip, "SCIPgetRowFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( scip->tree->actnodehaslp )
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
   CHECK_ABORT( checkStage(scip, "SCIPgetRowSolActivity", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( sol != NULL )
   {
      Real activity;

      CHECK_ABORT( SCIProwGetSolActivity(row, scip->set, scip->stat, sol, &activity) );
      return activity;
   }
   else if( scip->tree->actnodehaslp )
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
   CHECK_ABORT( checkStage(scip, "SCIPgetRowSolFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( sol != NULL )
   {
      Real feasibility;

      CHECK_ABORT( SCIProwGetSolFeasibility(row, scip->set, scip->stat, sol, &feasibility) );
      return feasibility;
   }
   else if( scip->tree->actnodehaslp )
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

   CHECK_OKAY( checkStage(scip, "SCIPprintRow", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE) );

   SCIProwPrint(row, file);

   return SCIP_OKAY;
}




/*
 * cutting plane methods
 */

/** adds cut to separation storage;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
RETCODE SCIPaddCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             cut,                /**< separated cut */
   Real             score               /**< separation score of cut (the larger, the better the cut) */
   )
{
   assert(cut != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPaddCut", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   assert(scip->tree->actnode != NULL);

   if( !scip->tree->actnodehaslp )
   {
      errorMessage("cannot add cuts, because node LP is not processed\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPsepastoreAddCut(scip->sepastore, scip->mem->solvemem, scip->set, scip->lp,
                  cut, score, (SCIPnodeGetDepth(scip->tree->actnode) == 0)) );
   
   return SCIP_OKAY;
}

/** if not already existing, adds row to global cut pool */
RETCODE SCIPpoolCut(
   SCIP*            scip,               /**< SCIP data structure */
   ROW*             row                 /**< cutting plane to add */
   )
{
   assert(row != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPpoolCut", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPcutpoolAddRow(scip->cutpool, scip->mem->solvemem, scip->set, row) );

   return SCIP_OKAY;
}

/** gets current number of rows in the global cut pool */
int SCIPgetPoolsize(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPoolsize", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   return SCIPcutpoolGetNCuts(scip->cutpool);
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
   CHECK_OKAY( checkStage(scip, "SCIPgetLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_OPTIMAL
      && SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_UNBOUNDED )
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

   CHECK_ABORT( checkStage(scip, "SCIPgetNLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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

   CHECK_ABORT( checkStage(scip, "SCIPgetNPrioLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPgetPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchcandGetPseudoCands(scip->branchcand, scip->set, scip->transprob, 
                  pseudocands, npseudocands, npriopseudocands) );

   return SCIP_OKAY;
}

/** gets branching candidates for pseudo solution branching (nonfixed variables) */
int SCIPgetNPseudoBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPbranchcandGetNPseudoCands(scip->branchcand);
}

/** gets number of branching candidates with maximal branch priority for pseudo solution branching */
int SCIPgetNPrioPseudoBranchCands(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   int npriopseudocands;

   CHECK_ABORT( checkStage(scip, "SCIPgetNPrioPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoCands(scip->branchcand);
}

/** calculates the branching score out of the downward and upward gain prediction */
Real SCIPgetBranchScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   Real             downgain,           /**< prediction of objective gain for branching downwards */
   Real             upgain              /**< prediction of objective gain for branching upwards */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBranchScore", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPbranchGetScore(scip->set, var, downgain, upgain);
}

/** creates a child node of the active node */
RETCODE SCIPcreateChild(
   SCIP*            scip,               /**< SCIP data structure */
   NODE**           node                /**< pointer to node data structure */
   )
{
   assert(node != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPcreateChild", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPnodeCreate(node, scip->mem->solvemem, scip->set, scip->stat, scip->tree) );
   
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
   CHECK_OKAY( checkStage(scip, "SCIPbranchVar", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPbranchLP", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchExecLP(scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, 
                  scip->sepastore, scip->branchcand, scip->eventqueue, result) );

   return SCIP_OKAY;
}

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN */
RETCODE SCIPbranchPseudo(
   SCIP*            scip,               /**< SCIP data structure */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPbranchPseudo", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPbranchExecPseudo(scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, 
                  scip->branchcand, scip->eventqueue, result) );

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
   CHECK_OKAY( checkStage(scip, "SCIPcreateSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolCreate(sol, scip->mem->solvemem, scip->stat, scip->tree, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current LP solution */
RETCODE SCIPcreateLPSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( !scip->tree->actnodehaslp )
   {
      errorMessage("LP solution does not exist\n");
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( SCIPsolCreateLPSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current pseudo solution */
RETCODE SCIPcreatePseudoSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreatePseudoSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolCreatePseudoSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current solution */
RETCODE SCIPcreateCurrentSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol,                /**< pointer to store the solution */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolCreateCurrentSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp, heur) );

   return SCIP_OKAY;
}

/** frees primal CIP solution */
RETCODE SCIPfreeSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL**            sol                 /**< pointer to the solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolFree(sol, scip->mem->solvemem) );

   return SCIP_OKAY;
}

/** links a primal solution to the current LP solution */
RETCODE SCIPlinkLPSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPlinkLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( !scip->tree->actnodehaslp )
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
   CHECK_OKAY( checkStage(scip, "SCIPlinkPseudoSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolLinkPseudoSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** links a primal solution to the current LP or pseudo solution */
RETCODE SCIPlinkCurrentSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPlinkCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPsolLinkCurrentSol(sol, scip->mem->solvemem, scip->set, scip->stat, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** clears a primal solution */
RETCODE SCIPclearSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsolClear(sol, scip->stat, scip->tree) );

   return SCIP_OKAY;
}

/** stores solution values of variables in solution's own array */
RETCODE SCIPunlinkSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPunlinkSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPsetSolVal", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPsetSolVals", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPincSolVal", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   CHECK_ABORT( checkStage(scip, "SCIPgetSolVal", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
   {
      Real solval;

      CHECK_ABORT( SCIPsolGetVal(sol, scip->set, scip->stat, var, &solval) );
      return solval;
   }
   else
   {
      CHECK_ABORT( checkStage(scip, "SCIPgetSolVal(sol==NULL)", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
      return SCIPgetVarSol(scip, var);
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

   CHECK_OKAY( checkStage(scip, "SCIPgetSolVals", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
   {
      for( v = 0; v < nvars; ++v )
      {
         CHECK_OKAY( SCIPsolGetVal(sol, scip->set, scip->stat, vars[v], &vals[v]) );
      }
   }
   else
   {
      CHECK_OKAY( checkStage(scip, "SCIPgetSolVals(sol==NULL)", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

      for( v = 0; v < nvars; ++v )
         vals[v] = SCIPgetVarSol(scip, vars[v]);
   }

   return SCIP_OKAY;
}

/** returns objective value of primal CIP solution w.r.t. original problem, or current LP/pseudo objective value */
Real SCIPgetSolOrigObj(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolOrigObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
      return SCIPprobExternObjval(scip->origprob, scip->set, 
         SCIPprobExternObjval(scip->transprob, scip->set, SCIPsolGetObj(sol)));
   else
   {
      CHECK_ABORT( checkStage(scip, "SCIPgetSolOrigObj(sol==NULL)", 
                      FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
      if( scip->tree->actnodehaslp )
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
   CHECK_ABORT( checkStage(scip, "SCIPgetSolTransObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
      return SCIPsolGetObj(sol);
   else
   {
      CHECK_ABORT( checkStage(scip, "SCIPgetSolTransObj(sol==NULL)", 
                      FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );
      if( scip->tree->actnodehaslp )
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
   CHECK_ABORT( checkStage(scip, "SCIPtransformObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPprobInternObjval(scip->transprob, scip->set, SCIPprobInternObjval(scip->origprob, scip->set, obj));
}

/** maps transformed objective value into original space */
Real SCIPretransformObj(
   SCIP*            scip,               /**< SCIP data structure */
   Real             obj                 /**< transformed objective value to retransform in original space */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPretransformObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPprobExternObjval(scip->origprob, scip->set, SCIPprobExternObjval(scip->transprob, scip->set, obj));
}

/** gets clock time, when this solution was found */
Real SCIPgetSolTime(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolTime", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetTime(sol);
}

/** gets node number, where this solution was found */
Longint SCIPgetSolNodenum(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolNodenum", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetNodenum(sol);
}

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
HEUR* SCIPgetSolHeur(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol                 /**< primal solution */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolHeur", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetHeur(sol);
}

/** outputs non-zero original variables of solution to file stream */
RETCODE SCIPprintSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( file == NULL )
      file = stdout;

   fprintf(file, "objective value:                 ");
   SCIPprintReal(scip, SCIPprobExternObjval(scip->origprob, scip->set, 
                    SCIPprobExternObjval(scip->transprob, scip->set, SCIPsolGetObj(sol))), file);
   fprintf(file, "\n");

   /**@todo transformed variables without original counterpart (e.g. priced variables) must also be displayed */
   CHECK_OKAY( SCIPsolPrint(sol, scip->set, scip->stat, scip->origprob, file) );

   return SCIP_OKAY;
}

/** outputs non-zero transformed variables of solution to file stream */
RETCODE SCIPprintTransSol(
   SCIP*            scip,               /**< SCIP data structure */
   SOL*             sol,                /**< primal solution */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintSolTrans", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( file == NULL )
      file = stdout;

   fprintf(file, "objective value:                 ");
   SCIPprintReal(scip, SCIPsolGetObj(sol), file);
   fprintf(file, "\n");

   CHECK_OKAY( SCIPsolPrint(sol, scip->set, scip->stat, scip->transprob, file) );
   
   return SCIP_OKAY;
}

/** gets best feasible primal solution found so far, or NULL if no solution has been found */
SOL* SCIPgetBestSol(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestSol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_FREESOLVE:
      return NULL;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      if( scip->primal->nsols > 0 ) 
      {
         assert(scip->primal->sols != NULL);
         assert(scip->primal->sols[0] != NULL);
         return scip->primal->sols[0];
      }
      return NULL;

   default:
      errorMessage("Unknown SCIP stage\n");
      abort();
   }
}

/** outputs best feasible primal solution found so far to file stream */
RETCODE SCIPprintBestSol(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   SOL* sol;

   CHECK_OKAY( checkStage(scip, "SCIPprintBestSol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

   CHECK_OKAY( checkStage(scip, "SCIPprintBestTransSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIProundSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPaddSol", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPaddSolFree", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalAddSolFree(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
                  scip->lp, scip->eventfilter, sol, stored) );

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
   CHECK_OKAY( checkStage(scip, "SCIPtrySol", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

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
   CHECK_OKAY( checkStage(scip, "SCIPtrySolFree", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   CHECK_OKAY( SCIPprimalTrySolFree(scip->primal, scip->mem->solvemem, scip->set, scip->stat, scip->transprob, scip->tree, 
                  scip->lp, scip->eventfilter, sol, checkintegrality, checklprows, stored) );

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
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcatchEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPeventfilterAdd(scip->eventfilter, scip->mem->solvemem, scip->set, eventtype, eventhdlr, eventdata) );

   return SCIP_OKAY;
}

/** drops a global event (stops to track event) */
RETCODE SCIPdropEvent(
   SCIP*            scip,               /**< SCIP data structure */
   EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPdropEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPeventfilterDel(scip->eventfilter, scip->mem->solvemem, scip->set, eventtype, eventhdlr, eventdata) );
   
   return SCIP_OKAY;
}

/** catches an objective value or domain change event on the given variable */
RETCODE SCIPcatchVarEvent(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to catch event for */
   EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcatchVarEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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

   CHECK_OKAY( SCIPvarCatchEvent(var, scip->mem->solvemem, scip->set, eventtype, eventhdlr, eventdata) );

   return SCIP_OKAY;
}

/** drops an objective value or domain change event (stops to track event) on the given variable */
RETCODE SCIPdropVarEvent(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< variable to drop event for */
   EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler when processing this event */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPdropVarEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( !SCIPvarIsTransformed(var) )
   {
      errorMessage("cannot drop events for original variables\n");
      return SCIP_INVALIDDATA;
   }

   CHECK_OKAY( SCIPvarDropEvent(var, scip->mem->solvemem, scip->set, eventtype, eventhdlr, eventdata) );
   
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
   CHECK_OKAY( checkStage(scip, "SCIPupdateVarPseudocost", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   SCIPvarUpdatePseudocost(var, scip->set, scip->stat, solvaldelta, objdelta, weight);

   return SCIP_OKAY;
}

/** gets the variable's pseudo cost value for the given direction */
Real SCIPgetVarPseudocost(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarPseudocost", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   return SCIPvarGetPseudocost(var, scip->stat, solvaldelta);
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
Real SCIPgetVarPseudocostCount(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   int              dir                 /**< branching direction: 0 (down), or 1 (up) */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetVarPseudocostCount", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   return SCIPvarGetPseudocostCount(var, dir);
}

/** gets the variable's pseudo cost score value for the given LP solution value */
Real SCIPgetVarPseudocostScore(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             solval              /**< variable's LP solution value */
   )
{
   Real frac;
   Real pscostdown;
   Real pscostup;

   CHECK_ABORT( checkStage(scip, "SCIPgetVarPseudocostScore", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE) );

   frac = solval - SCIPsetFloor(scip->set, solval);
   pscostdown = SCIPvarGetPseudocost(var, scip->stat, -frac);
   pscostup = SCIPvarGetPseudocost(var, scip->stat, 1.0-frac);

   return SCIPbranchGetScore(scip->set, var, pscostdown, pscostup);
}




/*
 * tree methods
 */

/** gets children of active node along with the number of children */
RETCODE SCIPgetChildren(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          children,           /**< pointer to store children array, or NULL if not needed */
   int*             nchildren           /**< pointer to store number of children, or NULL if not needed */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetChildren", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( children != NULL )
      *children = scip->tree->children;
   if( nchildren != NULL )
      *nchildren = scip->tree->nchildren;
   
   return SCIP_OKAY;
}

/** gets number of children of active node */
int SCIPgetNChildren(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNChildren", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return scip->tree->nchildren;
}

/** gets siblings of active node along with the number of siblings */
RETCODE SCIPgetSiblings(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          siblings,           /**< pointer to store siblings array, or NULL if not needed */
   int*             nsiblings           /**< pointer to store number of siblings, or NULL if not needed */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetSiblings", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   if( siblings != NULL )
      *siblings = scip->tree->siblings;
   if( nsiblings != NULL )
      *nsiblings = scip->tree->nsiblings;
   
   return SCIP_OKAY;
}

/** gets number of siblings of active node */
int SCIPgetNSiblings(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNSiblings", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return scip->tree->nsiblings;
}

/** gets leaves of the tree along with the number of leaves */
RETCODE SCIPgetLeaves(
   SCIP*            scip,               /**< SCIP data structure */
   NODE***          leaves,             /**< pointer to store leaves array, or NULL if not needed */
   int*             nleaves             /**< pointer to store number of leaves, or NULL if not needed */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPgetLeaves", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

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
   CHECK_ABORT( checkStage(scip, "SCIPgetNLeaves", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPnodepqLen(scip->tree->leaves);
}

/** gets the best child of the active node */
NODE* SCIPgetBestChild(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestChild", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetBestChild(scip->tree, scip->set);
}

/** gets the best sibling of the active node */
NODE* SCIPgetBestSibling(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestSibling", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetBestSibling(scip->tree, scip->set);
}

/** gets the best leaf from the node queue */
NODE* SCIPgetBestLeaf(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestLeaf", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetBestLeaf(scip->tree);
}

/** gets the best node from the tree (child, sibling, or leaf) */
NODE* SCIPgetBestNode(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestNode", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetBestNode(scip->tree, scip->set);
}

/** gets the node with smallest lower bound from the tree (child, sibling, or leaf) */
NODE* SCIPgetBestboundNode(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetBestboundNode", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetLowerboundNode(scip->tree, scip->set);
}




/*
 * statistic methods
 */

/** gets number of processed nodes, including the active node */
Longint SCIPgetNodenum(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNodenum", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->stat->nnodes;
}

/** gets number of nodes left in the tree (children + siblings + leaves) */
int SCIPgetNNodesLeft(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNNodesLeft", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPtreeGetNNodes(scip->tree);
}

/** gets total number of LPs solved so far */
int SCIPgetNLPs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNLPs", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nlps;
}

/** gets total number of simplex iterations used so far in primal and dual simplex */
Longint SCIPgetNLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNLPIterations", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nlpiterations;
}

/** gets total number of LPs solved so far for node relaxations */
int SCIPgetNNodeLPs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNNodeLPs", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nnodelps;
}

/** gets total number of simplex iterations used so far for node relaxations */
Longint SCIPgetNNodeLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNNodeLPIterations", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nnodelpiterations;
}

/** gets total number of LPs solved so far during diving */
int SCIPgetNDivingLPs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNDivingLPs", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->ndivinglps;
}

/** gets total number of simplex iterations used so far during diving */
Longint SCIPgetNDivingLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNDivingLPIterations", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->ndivinglpiterations;
}

/** gets total number of times, strong branching was called (each call represents solving two LPs) */
int SCIPgetNStrongbranchs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNStrongbranchs", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nstrongbranchs;
}

/** gets total number of simplex iterations used so far in strong branching */
Longint SCIPgetNStrongbranchLPIterations(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNStrongbranchLPIterations", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nsblpiterations;
}

/** gets number of pricing rounds performed so far at the current node */
int SCIPgetNPriceRounds(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPriceRounds", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return scip->stat->npricerounds;
}

/** get current number of variables in the pricing store */
int SCIPgetNPricevars(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPricevars", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPpricestoreGetNVars(scip->pricestore);
}

/** get total number of pricing variables found so far */
int SCIPgetNPricevarsFound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPricevarsFound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPpricestoreGetNVarsFound(scip->pricestore);
}

/** get total number of pricing variables applied to the LPs */
int SCIPgetNPricevarsApplied(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNPricevarsApplied", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPpricestoreGetNVarsApplied(scip->pricestore);
}

/** gets number of separation rounds performed so far at the current node */
int SCIPgetNSepaRounds(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNSepaRounds", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE) );

   return scip->stat->nseparounds;
}

/** get current number of cuts in the cut store */
int SCIPgetNCuts(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNCuts", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPsepastoreGetNCuts(scip->sepastore);
}

/** get total number of cuts found so far */
int SCIPgetNCutsFound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNCutsFound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPsepastoreGetNCutsFound(scip->sepastore);
}

/** get total number of cuts applied to the LPs */
int SCIPgetNCutsApplied(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNCutsApplied", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPsepastoreGetNCutsApplied(scip->sepastore);
}

/** get total number of conflict constraints found */
Longint SCIPgetNConflictsFound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNConflictsFound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPconflictGetNConflicts(scip->conflict)
      + SCIPlpconflictGetNConflicts(scip->lpconflict)
      + SCIPpseudoconflictGetNConflicts(scip->pseudoconflict);
}

/** gets depth of active node, or -1 if no active node exists */
int SCIPgetDepth(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetDepth", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   if( scip->tree->actnode != NULL )
      return SCIPnodeGetDepth(scip->tree->actnode);
   else
      return -1;
}

/** gets maximal depth of all processed nodes */
int SCIPgetMaxDepth(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetMaxDepth", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->maxdepth;
}

/** gets current plunging depth (succ. times, a child was selected as next node) */
int SCIPgetPlungeDepth(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPlungeDepth", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->plungedepth;
}

/** gets total number of active constraints at the current node */
int SCIPgetNActiveConss(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   int nconss;
   int h;

   CHECK_ABORT( checkStage(scip, "SCIPgetNActiveConss", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   nconss = 0;
   for( h = 0; h < scip->set->nconshdlrs; ++h )
      nconss += SCIPconshdlrGetNConss(scip->set->conshdlrs[h]);

   return nconss;
}

/** gets total number of enabled constraints at the current node */
int SCIPgetNEnabledConss(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   int nenabledconss;
   int h;

   CHECK_ABORT( checkStage(scip, "SCIPgetNEnabledConss", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   nenabledconss = 0;
   for( h = 0; h < scip->set->nconshdlrs; ++h )
      nenabledconss += SCIPconshdlrGetNEnabledConss(scip->set->conshdlrs[h]);

   return nenabledconss;
}

/** gets total number of globally valid constraints currently in the problem */
int SCIPgetNGlobalConss(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNGlobalConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

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
      errorMessage("Unknown SCIP stage\n");
      abort();
   }
}

/** gets average dual bound of all unprocessed nodes */
Real SCIPgetAvgDualbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetAvgDualbound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPprobExternObjval(scip->origprob, scip->set, 
      SCIPprobExternObjval(scip->transprob, scip->set, SCIPtreeGetAvgLowerbound(scip->tree, scip->primal->cutoffbound)));
}

/** gets average lower (dual) bound of all unprocessed nodes in transformed problem */
Real SCIPgetAvgLowerbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetAvgLowerbound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPtreeGetAvgLowerbound(scip->tree, scip->primal->cutoffbound);
}

/** gets global dual bound */
Real SCIPgetDualbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   Real lowerbound;

   CHECK_ABORT( checkStage(scip, "SCIPgetDualbound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

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
   CHECK_ABORT( checkStage(scip, "SCIPgetLowerbound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPtreeGetLowerbound(scip->tree, scip->set);
}

/** gets global primal bound (objective value of best solution or user objective limit) */
Real SCIPgetPrimalbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPrimalbound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPprobExternObjval(scip->origprob, scip->set, 
      SCIPprobExternObjval(scip->transprob, scip->set, scip->primal->upperbound));
}

/** gets global upper (primal) bound in transformed problem (objective value of best solution or user objective limit) */
Real SCIPgetUpperbound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetUpperbound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

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
   CHECK_ABORT( checkStage(scip, "SCIPgetCutoffbound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->cutoffbound;
}

/** gets current gap |(primalbound - dualbound)/dualbound| */
Real SCIPgetGap(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   Real primalbound;
   Real dualbound;

   CHECK_ABORT( checkStage(scip, "SCIPgetGap", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   if( SCIPsetIsInfinity(scip->set, SCIPgetLowerbound(scip)) )
      return 0.0;

   primalbound = SCIPgetPrimalbound(scip);
   dualbound = SCIPgetDualbound(scip);
   if( SCIPsetIsZero(scip->set, dualbound) || SCIPsetIsInfinity(scip->set, ABS(primalbound)) )
      return scip->set->infinity;
   else
      return ABS((primalbound - dualbound)/dualbound);
}

/** gets current gap |(upperbound - lowerbound)/lowerbound| in transformed problem */
Real SCIPgetTransGap(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   Real upperbound;
   Real lowerbound;

   CHECK_ABORT( checkStage(scip, "SCIPgetTransGap", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   upperbound = SCIPgetUpperbound(scip);
   lowerbound = SCIPgetLowerbound(scip);

   if( SCIPsetIsInfinity(scip->set, lowerbound) )
      return 0.0;
   else if( SCIPsetIsZero(scip->set, lowerbound) || SCIPsetIsInfinity(scip->set, upperbound) )
      return scip->set->infinity;
   else
      return ABS((upperbound - lowerbound)/lowerbound);
}

/** gets number of feasible primal solutions found so far */
Longint SCIPgetNSolsFound(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetNSolsFound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->nsolsfound;
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

   fprintf(file, "Constraints        :     Number  #Separate #Propagate    #EnfoLP    #EnfoPS    Cutoffs       Cuts    DomReds   Children\n");

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
         fprintf(file, " %10d%c%10lld %10lld %10lld %10lld %10lld %10lld %10lld %10lld\n",
            startnconss,
            maxnconss > startnconss ? '+' : ' ',
            SCIPconshdlrGetNSepaCalls(conshdlr), 
            SCIPconshdlrGetNPropCalls(conshdlr), 
            SCIPconshdlrGetNEnfoLPCalls(conshdlr),
            SCIPconshdlrGetNEnfoPSCalls(conshdlr),
            SCIPconshdlrGetNCutoffs(conshdlr), 
            SCIPconshdlrGetNCutsFound(conshdlr), 
            SCIPconshdlrGetNDomredsFound(conshdlr), 
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
void printRedcostStrengtheningStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   fprintf(file, "Reduced Cost Str.  :       Time      Calls      Found\n");
   fprintf(file, "  total            : %10.2f %10lld %10lld\n",
      SCIPclockGetTime(scip->stat->redcoststrtime),
      scip->stat->nredcoststrcalls,
      scip->stat->nredcoststrfound);
}

static
void printConflictStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file */
   )
{
   fprintf(file, "Conflict Analysis  :       Time      Calls  Conflicts\n");
   fprintf(file, "  propagation      : %10.2f %10lld %10lld\n",
      SCIPconflictGetTime(scip->conflict),
      SCIPconflictGetNCalls(scip->conflict),
      SCIPconflictGetNConflicts(scip->conflict));
   fprintf(file, "  infeasible LP    : %10.2f %10lld %10lld   (%lld LP iterations)\n",
      SCIPlpconflictGetTime(scip->lpconflict),
      SCIPlpconflictGetNCalls(scip->lpconflict),
      SCIPlpconflictGetNConflicts(scip->lpconflict),
      SCIPlpconflictGetNLPIterations(scip->lpconflict));
   fprintf(file, "  pseudo solution  : %10.2f %10lld %10lld\n",
      SCIPpseudoconflictGetTime(scip->pseudoconflict),
      SCIPpseudoconflictGetNCalls(scip->pseudoconflict),
      SCIPpseudoconflictGetNConflicts(scip->pseudoconflict));
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

   fprintf(file, "Branching Rules    :       Time      Calls    Cutoffs       Cuts    DomReds   Children\n");

   for( i = 0; i < scip->set->nbranchrules; ++i )
      fprintf(file, "  %-17.17s: %10.2f %10lld %10lld %10lld %10lld %10lld\n",
         SCIPbranchruleGetName(scip->set->branchrules[i]),
         SCIPbranchruleGetTime(scip->set->branchrules[i]),
         SCIPbranchruleGetNLPCalls(scip->set->branchrules[i]) + SCIPbranchruleGetNPseudoCalls(scip->set->branchrules[i]),
         SCIPbranchruleGetNCutoffs(scip->set->branchrules[i]),
         SCIPbranchruleGetNCutsFound(scip->set->branchrules[i]),
         SCIPbranchruleGetNDomredsFound(scip->set->branchrules[i]),
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

   fprintf(file, "  strong branching : %10.2f %10d %10lld %10.2f",
      SCIPclockGetTime(scip->stat->strongbranchtime),
      scip->stat->nstrongbranchs,
      scip->stat->nsblpiterations,
      scip->stat->nstrongbranchs > 0 ? (Real)scip->stat->nsblpiterations/(Real)scip->stat->nstrongbranchs : 0.0);
   if( SCIPclockGetTime(scip->stat->strongbranchtime) >= 0.01 )
      fprintf(file, " %10.2f\n", (Real)scip->stat->nsblpiterations/SCIPclockGetTime(scip->stat->strongbranchtime));
   else
      fprintf(file, "          -\n");
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
   fprintf(file, "  nodes            : %10lld\n", scip->stat->nnodes);
   fprintf(file, "  max depth        : %10d\n", scip->stat->maxdepth);
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
   Real bestsol;
   Real gap;

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->primal != NULL);
   assert(file != NULL);

   primalbound = SCIPgetPrimalbound(scip);
   dualbound = SCIPgetDualbound(scip);
   gap = SCIPgetGap(scip);

   fprintf(file, "Solution           :\n");
   fprintf(file, "  Solutions found  : %10lld\n", scip->primal->nsolsfound);
   if( SCIPsetIsInfinity(scip->set, ABS(primalbound)) )
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
   if( SCIPsetIsInfinity(scip->set, ABS(dualbound)) )
      fprintf(file, "  Dual Bound       :          -\n");
   else
      fprintf(file, "  Dual Bound       : %+21.14e\n", dualbound);
   if( SCIPsetIsInfinity(scip->set, gap) )
      fprintf(file, "  Gap              :   infinite\n");
   else
      fprintf(file, "  Gap              : %10.2f %%\n", 100.0 * gap);
}

/** outputs SCIP status */
RETCODE SCIPprintStatus(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintStatus", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( file == NULL )
      file = stdout;

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      fprintf(file, "initialization");
      break;
   case SCIP_STAGE_PROBLEM:
      fprintf(file, "problem creation");
      break;
   case SCIP_STAGE_INITSOLVE:
      fprintf(file, "solving process initialization");
      break;
   case SCIP_STAGE_PRESOLVING:
      fprintf(file, "presolving is running");
      break;
   case SCIP_STAGE_PRESOLVED:
      fprintf(file, "problem is presolved");
      break;
   case SCIP_STAGE_SOLVING:
      if( SCIPsolveIsStopped(scip->set, scip->stat) )
      {
         fprintf(file, "solving was interrupted [");
         SCIPsolvePrintStopReason(scip->set, scip->stat, file);
         fprintf(file, "]");
      }
      else
         fprintf(file, "solving process is running");
      break;
   case SCIP_STAGE_SOLVED:
      fprintf(file, "problem is solved");
      if( scip->primal->nsols == 0 )
         fprintf(file, " [infeasible]");
      else if( SCIPsetIsInfinity(scip->set, -scip->primal->upperbound) )
         fprintf(file, " [unbounded]");
      else
         fprintf(file, " [optimal solution found]");
      break;
   case SCIP_STAGE_FREESOLVE:
      fprintf(file, "solving process deinitialization");
      break;
   }

   return SCIP_OKAY;
}

/** outputs solving statistics */
RETCODE SCIPprintStatistics(
   SCIP*            scip,               /**< SCIP data structure */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPprintStatistics", TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   if( file == NULL )
      file = stdout;

   switch( scip->stage )
   {
   case SCIP_STAGE_INIT:
      fprintf(file, "SCIP Status        : initialization\n");
      fprintf(file, "Original Problem   : no problem exists.\n");
      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
      fprintf(file, "SCIP Status        : problem creation / modification\n");
      fprintf(file, "Original Problem   :\n");
      SCIPprobPrintStatistics(scip->origprob, file);
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
      fprintf(file, "SCIP Status        : ");
      CHECK_OKAY( SCIPprintStatus(scip, file) );
      fprintf(file, "\n");
      fprintf(file, "Solving Time       : %10.2f\n", SCIPclockGetTime(scip->stat->solvingtime));
      fprintf(file, "Original Problem   :\n");
      SCIPprobPrintStatistics(scip->origprob, file);
      fprintf(file, "Transformed Problem:\n");
      SCIPprobPrintStatistics(scip->transprob, file);
      printPresolverStatistics(scip, file);
      printConstraintStatistics(scip, file);
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      fprintf(file, "SCIP Status        : ");
      CHECK_OKAY( SCIPprintStatus(scip, file) );
      fprintf(file, "\n");
      fprintf(file, "Solving Time       : %10.2f\n", SCIPclockGetTime(scip->stat->solvingtime));
      fprintf(file, "Original Problem   :\n");
      SCIPprobPrintStatistics(scip->origprob, file);
      fprintf(file, "Transformed Problem:\n");
      SCIPprobPrintStatistics(scip->transprob, file);
      printPresolverStatistics(scip, file);
      printConstraintStatistics(scip, file);
      printConstraintTimingStatistics(scip, file);
      printRedcostStrengtheningStatistics(scip, file);
      printConflictStatistics(scip, file);
      printSeparatorStatistics(scip, file);
      printPricerStatistics(scip, file);
      printBranchruleStatistics(scip, file);
      printHeuristicStatistics(scip, file);
      printLPStatistics(scip, file);
      printTreeStatistics(scip, file);
      printSolutionStatistics(scip, file);
      
#if 0
      { /*????????????????????????????????*/
         /**@todo remove this! */
         int v;
         for( v = 0; v < scip->transprob->nvars; ++v )
         {
            if( SCIPvarGetNBranchings(scip->transprob->vars[v]) > 0 )
            {
               char s[MAXSTRLEN];
               sprintf(s, "<%s>", SCIPvarGetName(scip->transprob->vars[v]));
               printf(" %16s: %6lld branchings, %8lld inferences  ->  %7.2f inferences/branching (avg. depth: %7.2f)\n",
                  s, SCIPvarGetNBranchings(scip->transprob->vars[v]),
                  SCIPvarGetNInferences(scip->transprob->vars[v]), SCIPvarGetAvgInferences(scip->transprob->vars[v]),
                  SCIPvarGetAvgBranchdepth(scip->transprob->vars[v]));
            }
         }
      }
#endif

      return SCIP_OKAY;

   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}




/*
 * timing methods
 */

/** creates a clock using the default clock type */
RETCODE SCIPcreateClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPclockCreate(clock, SCIP_CLOCKTYPE_DEFAULT) );

   return SCIP_OKAY;
}

/** creates a clock counting the CPU user seconds */
RETCODE SCIPcreateCPUClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateCPUClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPclockCreate(clock, SCIP_CLOCKTYPE_CPU) );

   return SCIP_OKAY;
}

/** creates a clock counting the wall clock seconds */
RETCODE SCIPcreateWallClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateWallClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPclockCreate(clock, SCIP_CLOCKTYPE_WALL) );

   return SCIP_OKAY;
}

/** frees a clock */
RETCODE SCIPfreeClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK**          clock               /**< pointer to clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockFree(clock);

   return SCIP_OKAY;
}

/** resets the time measurement of a clock to zero and completely stops the clock */
RETCODE SCIPresetClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPresetClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockReset(clock);

   return SCIP_OKAY;
}

/** starts the time measurement of a clock */
RETCODE SCIPstartClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPstartClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockStart(clock, scip->set);

   return SCIP_OKAY;
}

/** stops the time measurement of a clock */
RETCODE SCIPstopClock(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPstopClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockStop(clock, scip->set);

   return SCIP_OKAY;
}

/** gets the measured time of a clock in seconds */
Real SCIPgetClockTime(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock               /**< clock timer */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetClockTime", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPclockGetTime(clock);
}

/** sets the measured time of a clock to the given value in seconds */
RETCODE SCIPsetClockTime(
   SCIP*            scip,               /**< SCIP data structure */
   CLOCK*           clock,              /**< clock timer */
   Real             sec                 /**< time in seconds to set the clock's timer to */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetClockTime", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockSetTime(clock, sec);

   return SCIP_OKAY;
}

/** gets the current total SCIP time in seconds */
Real SCIPgetTotalTime(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetTotalTime", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPclockGetTime(scip->totaltime);
}

/** gets the current solving time in seconds */
Real SCIPgetSolvingTime(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetSolvingTime", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPclockGetTime(scip->stat->solvingtime);
}

/** gets the current presolving time in seconds */
Real SCIPgetPresolvingTime(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetPresolvingTime", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

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

   return scip->set->infinity;
}

/** returns value treated as zero */
Real SCIPepsilon(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->epsilon;
}

/** returns value treated as zero for sums of floating point values */
Real SCIPsumepsilon(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->sumepsilon;
}

/** returns feasibility tolerance for constraints */
Real SCIPfeastol(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
 
   return scip->set->feastol;
}

/** returns feasibility tolerance for reduced costs */
Real SCIPdualfeastol(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
 
   return scip->set->dualfeastol;
}

/** sets the feasibility tolerance for constraints */
RETCODE SCIPsetFeastol(
   SCIP*            scip,               /**< SCIP data structure */
   Real             feastol             /**< new feasibility tolerance for constraints */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetFeastol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* mark the LP unsolved, if the feasibility tolerance was tightened */
   if( scip->lp != NULL && feastol < scip->set->feastol )
      scip->lp->solved = FALSE;

   /* change the settings */
   CHECK_OKAY( SCIPsetSetFeastol(scip->set, feastol) );

   return SCIP_OKAY;
}

/** sets the feasibility tolerance for reduced costs */
RETCODE SCIPsetDualfeastol(
   SCIP*            scip,               /**< SCIP data structure */
   Real             dualfeastol         /**< new feasibility tolerance for reduced costs */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPsetDualfeastol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* mark the LP unsolved, if the dual feasibility tolerance was tightened */
   if( scip->lp != NULL && dualfeastol < scip->set->dualfeastol )
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

/** checks, if the cut's activity is more then cutvioleps larger than the given right hand side;
 *  both, the activity and the rhs, should be normed
 */
Bool SCIPisCutViolated(
   SCIP*            scip,               /**< SCIP data structure */
   Real             cutactivity,        /**< activity of the cut */
   Real             cutrhs              /**< right hand side value of the cut */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->tree != NULL);
   assert(scip->tree->actnode != NULL);

   return SCIPsetIsCutViolated(scip->set, (SCIPnodeGetDepth(scip->tree->actnode) == 0), cutactivity, cutrhs);
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

   if( scip->set->exactsolve )
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

   if( scip->set->exactsolve )
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

/** returns fractional part of value, i.e. ceil(x) - x */
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
   CHECK_ABORT( checkStage(scip, "SCIPmemhdr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      return scip->mem->solvemem;

   default:
      errorMessage("Unknown SCIP stage\n");
      return NULL;
   }  /*lint !e788*/
}

/** returns the total number of bytes used in block memory */
Longint SCIPgetMemUsed(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_ABORT( checkStage(scip, "SCIPgetMemUsed", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

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
RETCODE SCIPallocBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   int              size                /**< required size in bytes of buffer */
   )
{
   assert(ptr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPallocBuffer", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetAllocBufferSize(scip->set, ptr, size) );

   return SCIP_OKAY;
}

/** allocates a memory buffer with at least the given size and copies the given memory into the buffer */
RETCODE SCIPduplicateBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to store the buffer */
   void*            source,             /**< memory block to copy into the buffer */
   int              size                /**< required size in bytes of buffer */
   )
{
   assert(ptr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPduplicateBuffer", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetDuplicateBufferSize(scip->set, ptr, source, size) );

   return SCIP_OKAY;
}

/** reallocates a memory buffer to at least the given size */
RETCODE SCIPreallocBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to the buffer */
   int              size                /**< required size in bytes of buffer */
   )
{
   assert(ptr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPreallocBuffer", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   CHECK_OKAY( SCIPsetReallocBufferSize(scip->set, ptr, size) );

   return SCIP_OKAY;
}

/** frees a memory buffer */
RETCODE SCIPfreeBuffer(
   SCIP*            scip,               /**< SCIP data structure */
   void**           ptr,                /**< pointer to the buffer */
   int              dummysize           /**< used to get a safer define for SCIPfreeBufferSize/Array */
   )
{  /*lint --e{715}*/
   assert(ptr != NULL);

   CHECK_OKAY( checkStage(scip, "SCIPfreeBuffer", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetFreeBufferSize(scip->set, ptr);
   
   return SCIP_OKAY;
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
   CHECK_OKAY( checkStage(scip, "SCIPcreateRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayCreate(realarray, SCIPmemhdr(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
RETCODE SCIPfreeRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY**      realarray           /**< pointer to the real array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_OKAY( checkStage(scip, "SCIPextendRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayExtend(realarray, scip->set, minidx, maxidx) );
   
   return SCIP_OKAY;
}

/** clears a dynamic real array */
RETCODE SCIPclearRealarray(
   SCIP*            scip,               /**< SCIP data structure */
   REALARRAY*       realarray           /**< dynamic real array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_ABORT( checkStage(scip, "SCIPgetRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_OKAY( checkStage(scip, "SCIPsetRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_OKAY( checkStage(scip, "SCIPincRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPrealarrayIncVal(realarray, scip->set, idx, incval) );
   
   return SCIP_OKAY;
}

/** creates a dynamic array of int values */
RETCODE SCIPcreateIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY**       intarray            /**< pointer to store the int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayCreate(intarray, SCIPmemhdr(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of int values */
RETCODE SCIPfreeIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY**       intarray            /**< pointer to the int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_OKAY( checkStage(scip, "SCIPextendIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayExtend(intarray, scip->set, minidx, maxidx) );
   
   return SCIP_OKAY;
}

/** clears a dynamic int array */
RETCODE SCIPclearIntarray(
   SCIP*            scip,               /**< SCIP data structure */
   INTARRAY*        intarray            /**< dynamic int array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_ABORT( checkStage(scip, "SCIPgetIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_OKAY( checkStage(scip, "SCIPsetIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_OKAY( checkStage(scip, "SCIPincIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPintarrayIncVal(intarray, scip->set, idx, incval) );
   
   return SCIP_OKAY;
}

/** creates a dynamic array of bool values */
RETCODE SCIPcreateBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY**      boolarray           /**< pointer to store the bool array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPcreateBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarrayCreate(boolarray, SCIPmemhdr(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of bool values */
RETCODE SCIPfreeBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY**      boolarray           /**< pointer to the bool array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPfreeBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_OKAY( checkStage(scip, "SCIPextendBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarrayExtend(boolarray, scip->set, minidx, maxidx) );
   
   return SCIP_OKAY;
}

/** clears a dynamic bool array */
RETCODE SCIPclearBoolarray(
   SCIP*            scip,               /**< SCIP data structure */
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   CHECK_OKAY( checkStage(scip, "SCIPclearBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_ABORT( checkStage(scip, "SCIPgetBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
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
   CHECK_OKAY( checkStage(scip, "SCIPsetBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   CHECK_OKAY( SCIPboolarraySetVal(boolarray, scip->set, idx, val) );
   
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
