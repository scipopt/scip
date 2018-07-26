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

/**@file   scip_message.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for message handling
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_MESSAGE_H__
#define __SCIP_SCIP_MESSAGE_H__


#include "scip/def.h"
#include "scip/type_message.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/var.h"
#include "scip/cons.h"
#include "scip/solve.h"
#include "scip/debug.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup MessageOutputMethods
 *
 * @{
 */

/* if we have a C99 compiler */
#ifdef SCIP_HAVE_VARIADIC_MACROS

/** prints a debugging message if SCIP_DEBUG flag is set */
#ifdef SCIP_DEBUG
#define SCIPdebugMsg(scip, ...)         SCIPprintDebugMessage(scip, __FILE__, __LINE__, __VA_ARGS__)
#define SCIPdebugMsgPrint(scip, ...)    SCIPdebugMessagePrint(scip, __VA_ARGS__)
#else
#define SCIPdebugMsg(scip, ...)         while ( FALSE ) SCIPprintDebugMessage(scip, __FILE__, __LINE__, __VA_ARGS__)
#define SCIPdebugMsgPrint(scip, ...)    while ( FALSE ) SCIPdebugMessagePrint(scip, __VA_ARGS__)
#endif

#else
/* if we do not have a C99 compiler, use a workaround that prints a message, but not the file and linenumber */

/** prints a debugging message if SCIP_DEBUG flag is set */
#ifdef SCIP_DEBUG
#define SCIPdebugMsg                    printf("debug: "), SCIPdebugMessagePrint
#define SCIPdebugMsgPrint               SCIPdebugMessagePrint
#else
#define SCIPdebugMsg                    while ( FALSE ) SCIPdebugMessagePrint
#define SCIPdebugMsgPrint               while ( FALSE ) SCIPdebugMessagePrint
#endif

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
EXTERN
SCIP_RETCODE SCIPsetMessagehdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   );

/** returns the currently installed message handler
 *
 *  @return the currently installed message handler, or NULL if messages are currently suppressed
 */
EXTERN
SCIP_MESSAGEHDLR* SCIPgetMessagehdlr(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the log file name for the currently installed message handler */
EXTERN
void SCIPsetMessagehdlrLogfile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of log file, or NULL (no log) */
   );

/** sets the currently installed message handler to be quiet (or not) */
EXTERN
void SCIPsetMessagehdlrQuiet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   );

/** prints a warning message via the message handler */
EXTERN
void SCIPwarningMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a debug message */
EXTERN
void SCIPprintDebugMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline,         /**< line in the source file where the function was called */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a debug message without precode */
EXTERN
void SCIPdebugMessagePrint(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a dialog message that requests user interaction or is a direct response to a user interactive command */
EXTERN
void SCIPdialogMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message */
EXTERN
void SCIPinfoMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message depending on the verbosity level */
EXTERN
void SCIPverbMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** returns the current message verbosity level
 *
 *  @return message verbosity level of SCIP
 *
 *  @see \ref SCIP_VerbLevel "SCIP_VERBLEVEL" for a list of all verbosity levels
 */
EXTERN
SCIP_VERBLEVEL SCIPgetVerbLevel(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**@} */

#ifdef __cplusplus
}
#endif

#endif
