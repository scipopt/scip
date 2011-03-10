/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   message.h
 * @brief  message output methods
 * @author Tobias Achterberg
 *
 * Because the message functions are implemented as defines with more than one
 * function call, they shouldn't be used as a single statement like in:
 *   if( error )
 *      SCIPerrorMessage("an error occured");
 * because this would produce the following macro extension:
 *   if( error )
 *      printf(("[%s:%d] ERROR: ", __FILE__, __LINE__);
 *   printf(("an error occured");
 * Instead, they should be protected with brackets:
 *   if( error )
 *   {
 *      SCIPerrorMessage("an error occured");
 *   }
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MESSAGE_H__
#define __SCIP_MESSAGE_H__


#include <stdio.h>
#include <stdarg.h>

#include "scip/def.h"
#include "scip/type_message.h"
#include "scip/struct_message.h"
#include "scip/pub_message.h"

#ifdef __cplusplus
extern "C" {
#endif

#define infoMessage                     SCIPmessagePrintInfo

#define printErrorHeader                SCIPmessagePrintErrorHeader
#define printError                      SCIPmessagePrintError
#define printInfo                       SCIPmessagePrintInfo


/** creates a message handler */
extern
SCIP_RETCODE SCIPmessagehdlrCreate(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   SCIP_Bool             bufferedoutput,     /**< should the output be buffered up to the next newline? */
   SCIP_DECL_MESSAGEERROR((*messageerror)),  /**< error message print method of message handler */
   SCIP_DECL_MESSAGEWARNING((*messagewarning)),/**< warning message print method of message handler */
   SCIP_DECL_MESSAGEDIALOG((*messagedialog)),/**< dialog message print method of message handler */
   SCIP_DECL_MESSAGEINFO ((*messageinfo)),   /**< info message print method of message handler */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< message handler data */
   );

/** frees message handler */
extern
void SCIPmessagehdlrFree(
   SCIP_MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   );

/** installs the given message handler */
extern
void SCIPmessageSetHandler(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   );

/** installs the default message handler that prints messages to stdout or stderr */
extern
void SCIPmessageSetDefaultHandler(
   void
   );

/** returns the currently installed message handler, or NULL if messages are currently suppressed */
extern
SCIP_MESSAGEHDLR* SCIPmessageGetHandler(
   void
   );

/** prints a message, acting like the printf() command */
extern
void SCIPmessagePrintInfo(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message, acting like the vprintf() command */
extern
void SCIPmessageVPrintInfo(
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a message into a file, acting like the fprintf() command */
extern
void SCIPmessageFPrintInfo(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message into a file, acting like the vfprintf() command */
extern
void SCIPmessageVFPrintInfo(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a dialog message that requests user interaction, acting like the printf() command */
extern
void SCIPmessagePrintDialog(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a dialog message that requests user interaction, acting like the vprintf() command */
extern
void SCIPmessageVPrintDialog(
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a dialog message that requests user interaction into a file, acting like the fprintf() command */
extern
void SCIPmessageFPrintDialog(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a dialog message that requests user interaction into a file, acting like the vfprintf() command */
extern
void SCIPmessageVFPrintDialog(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a message depending on the verbosity level, acting like the printf() command */
extern
void SCIPmessagePrintVerbInfo(
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message depending on the verbosity level, acting like the vprintf() command */
extern
void SCIPmessageVPrintVerbInfo(
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a message into a file depending on the verbosity level, acting like the fprintf() command */
extern
void SCIPmessageFPrintVerbInfo(
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message into a file depending on the verbosity level, acting like the vfprintf() command */
extern
void SCIPmessageVFPrintVerbInfo(
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

#ifdef __cplusplus
}
#endif

#endif

