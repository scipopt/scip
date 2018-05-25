/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_message.c
 * @brief  public methods for message handling
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Robert Lion Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif


#include "lpi/lpi.h"
#include "nlpi/exprinterpret.h"
#include "nlpi/nlpi.h"
#include "scip/benders.h"
#include "scip/benderscut.h"
#include "scip/branch.h"
#include "scip/branch_nodereopt.h"
#include "scip/clock.h"
#include "scip/compr.h"
#include "scip/concsolver.h"
#include "scip/concurrent.h"
#include "scip/conflict.h"
#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/cons_linear.h"
#include "scip/cutpool.h"
#include "scip/cuts.h"
#include "scip/debug.h"
#include "scip/def.h"
#include "scip/dialog.h"
#include "scip/dialog_default.h"
#include "scip/disp.h"
#include "scip/event.h"
#include "scip/heur.h"
#include "scip/heur_ofins.h"
#include "scip/heur_reoptsols.h"
#include "scip/heur_trivialnegation.h"
#include "scip/heuristics.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/interrupt.h"
#include "scip/lp.h"
#include "scip/mem.h"
#include "scip/message_default.h"
#include "scip/misc.h"
#include "scip/nlp.h"
#include "scip/nodesel.h"
#include "scip/paramset.h"
#include "scip/presol.h"
#include "scip/presolve.h"
#include "scip/pricer.h"
#include "scip/pricestore.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/prop.h"
#include "scip/reader.h"
#include "scip/relax.h"
#include "scip/reopt.h"
#include "scip/retcode.h"
#include "scip/scipbuildflags.h"
#include "scip/scipcoreplugins.h"
#include "scip/scipgithash.h"
#include "scip/sepa.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/syncstore.h"
#include "scip/table.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/visual.h"
#include "xml/xml.h"

#include "scip/scip_message.h"

#include "scip/pub_message.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** installs the given message handler, such that all messages are passed to this handler. A messages handler can be
 *  created via SCIPmessagehdlrCreate().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note The currently installed messages handler gets freed if this SCIP instance is its last user (w.r.t. capture/release).
 */
SCIP_RETCODE SCIPsetMessagehdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   )
{
   int i;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetMessagehdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->set->nlpis != NULL || scip->set->nnlpis == 0);

   /* update message handler in NLP solver interfaces */
   for( i = 0; i < scip->set->nnlpis; ++i )
   {
      assert(scip->set->nlpis[i] != NULL);

      SCIP_CALL( SCIPnlpiSetMessageHdlr(scip->set->nlpis[i], messagehdlr) );
   }

   SCIPmessagehdlrCapture(messagehdlr);

   SCIP_CALL( SCIPmessagehdlrRelease(&scip->messagehdlr) );
   assert(scip->messagehdlr == NULL);

   scip->messagehdlr = messagehdlr;

   return SCIP_OKAY;
}

/** returns the currently installed message handler
 *
 *  @return the currently installed message handler, or NULL if messages are currently suppressed
 */
SCIP_MESSAGEHDLR* SCIPgetMessagehdlr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return scip->messagehdlr;
}

/** sets the log file name for the currently installed message handler */
void SCIPsetMessagehdlrLogfile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of log file, or NULL (no log) */
   )
{
   if( scip->messagehdlr != NULL )
   {
      SCIPmessagehdlrSetLogfile(scip->messagehdlr, filename);
   }
}

/** sets the currently installed message handler to be quiet (or not) */
void SCIPsetMessagehdlrQuiet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   )
{
   if( scip->messagehdlr != NULL )
   {
      SCIPmessagehdlrSetQuiet(scip->messagehdlr, quiet);
   }
}

/** prints a warning message via the message handler */
void SCIPwarningMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert(scip != NULL);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintWarning(scip->messagehdlr, formatstr, ap);
   va_end(ap);
}

/** prints a debug message */
void SCIPprintDebugMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline,         /**< line in the source file where the function was called */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   int subscipdepth = 0;
   va_list ap;

   assert( sourcefile != NULL );
   assert( scip != NULL );

   if ( scip->stat != NULL )
      subscipdepth = scip->stat->subscipdepth;
   if ( subscipdepth > 0 )
      SCIPmessageFPrintInfo(scip->messagehdlr, NULL, "%d: [%s:%d] debug: ", subscipdepth, sourcefile, sourceline);
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, NULL, "[%s:%d] debug: ", sourcefile, sourceline);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(scip->messagehdlr, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a debug message without precode */
void SCIPdebugMessagePrint(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert( scip != NULL );

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(scip->messagehdlr, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a dialog message that requests user interaction or is a direct response to a user interactive command */
void SCIPdialogMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert(scip != NULL);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintDialog(scip->messagehdlr, file, formatstr, ap);
   va_end(ap);
}

/** prints a message */
void SCIPinfoMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert(scip != NULL);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(scip->messagehdlr, file, formatstr, ap);
   va_end(ap);
}

/** prints a message depending on the verbosity level */
void SCIPverbMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert(scip != NULL);
   assert(scip->set != NULL);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, msgverblevel, file, formatstr, ap);
   va_end(ap);
}

/** returns the current message verbosity level
 *
 *  @return message verbosity level of SCIP
 *
 *  @see \ref SCIP_VerbLevel "SCIP_VERBLEVEL" for a list of all verbosity levels
 */
SCIP_VERBLEVEL SCIPgetVerbLevel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->disp_verblevel;
}
