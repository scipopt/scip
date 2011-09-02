/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_message.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for message output
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_MESSAGE_H__
#define __SCIP_PUB_MESSAGE_H__

#include "scip/def.h"
#include "scip/type_message.h"

#ifdef NDEBUG
#include "scip/struct_message.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** prints an error message */
#define SCIPerrorMessage                SCIPmessagePrintErrorHeader(__FILE__, __LINE__); \
                                        SCIPmessagePrintError

/** define used in blockmemshell/memory.c */
#define printErrorHeader                SCIPmessagePrintErrorHeader
#define printError                      SCIPmessagePrintError

#ifdef SCIP_DEBUG

/** executes command only if SCIP_DEBUG flag is set */
#define SCIPdebug(x)                        x

/** prints a debugging message if SCIP_DEBUG flag is set */
#define SCIPdebugMessage                printf("[%s:%d] debug: ", __FILE__, __LINE__); printf

/** executes printf command only if SCIP_DEBUG flag is set */
#define SCIPdebugPrintf                 printf

#else

/** executes command only if SCIP_DEBUG flag is set */
#define SCIPdebug(x)                        /**/

/** prints a debugging message if SCIP_DEBUG flag is set */
#define SCIPdebugMessage                while( FALSE ) printf

/** executes printf command only if SCIP_DEBUG flag is set */
#define SCIPdebugPrintf                 while( FALSE ) printf

#endif


/** Creates a message handler which deals with warning, information, and dialog (interactive shell) methods.
 *
 *  @note The message handler does not handle error messages. For that see SCIPmessageSetErrorPrinting()
 */
extern
SCIP_RETCODE SCIPmessagehdlrCreate(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   SCIP_Bool             bufferedoutput,     /**< should the output be buffered up to the next newline? */
   const char*           filename,           /**< name of log file, or NULL (stdout) */
   SCIP_Bool             quiet,              /**< should screen messages be suppressed? */
   SCIP_DECL_MESSAGEWARNING((*messagewarning)),/**< warning message print method of message handler */
   SCIP_DECL_MESSAGEDIALOG((*messagedialog)),/**< dialog message print method of message handler */
   SCIP_DECL_MESSAGEINFO ((*messageinfo)),   /**< info message print method of message handler */
   SCIP_DECL_MESSAGEHDLRFREE((*messagehdlrfree)), /**< destructor of message handler to free message handler data */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< message handler data */
   );

/** frees message handler */
extern
SCIP_RETCODE SCIPmessagehdlrFree(
   SCIP_MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   );

/** sets the user data of the message handler */
extern
SCIP_RETCODE SCIPmessagehdlrSetData(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler; must not be NULL */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< new message handler data to attach to the handler */
   );

/** sets the log file name for the message handler */
extern
void SCIPmessagehdlrSetLogfile(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           filename            /**< log file name where to copy messages into, or NULL */
   );

/** sets the messages handler to be quiet */
extern
void SCIPmessagehdlrSetQuiet(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   );

/** prints a message, acting like the printf() command */
extern
void SCIPmessagePrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message, acting like the vprintf() command */
extern
void SCIPmessageVPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a message into a file, acting like the fprintf() command */
extern
void SCIPmessageFPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message into a file, acting like the vfprintf() command */
extern
void SCIPmessageVFPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints the header with source file location for an error message */
extern
void SCIPmessagePrintWarningHeader(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline          /**< line in the source file where the function was called */
   );

/** prints a warning message, acting like the printf() command */
extern
void SCIPmessagePrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a warning message, acting like the vprintf() command */
extern
void SCIPmessageVPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a warning message into a file, acting like the fprintf() command */
extern
void SCIPmessageFPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a warning message into a file, acting like the vfprintf() command */
extern
void SCIPmessageVFPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a dialog message that requests user interaction, acting like the printf() command */
extern
void SCIPmessagePrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a dialog message that requests user interaction, acting like the vprintf() command */
extern
void SCIPmessageVPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a dialog message that requests user interaction into a file, acting like the fprintf() command */
extern
void SCIPmessageFPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a dialog message that requests user interaction into a file, acting like the vfprintf() command */
extern
void SCIPmessageVFPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a message depending on the verbosity level, acting like the printf() command */
extern
void SCIPmessagePrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message depending on the verbosity level, acting like the vprintf() command */
extern
void SCIPmessageVPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a message into a file depending on the verbosity level, acting like the fprintf() command */
extern
void SCIPmessageFPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message into a file depending on the verbosity level, acting like the vfprintf() command */
extern
void SCIPmessageVFPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints the header with source file location for an error message using the static message handler */
extern
void SCIPmessagePrintErrorHeader(
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline          /**< line in the source file where the function was called */
   );

/** prints an error message, acting like the printf() command using the static message handler */
extern
void SCIPmessagePrintError(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** Method to set the error printing method. Setting the error printing method to NULL will suspend all error methods.
 *
 *  @note The error printing method is static variable. That means all occurring errors are handled via that methods
 */
extern
void SCIPmessageSetErrorPrinting(
   SCIP_DECL_ERRORPRINTING((*errorPrinting)),/**< error message print method of message handler, or NULL */
   void*                 data                /**< data pointer which will be passed to the error printing method, or NULL */
   );

/** Method to set the error printing method to default version prints everything the stderr.
 *
 *  @note The error printing method is static variable. That means all occurring errors are handled via that methods
 */
extern
void SCIPmessageSetErrorPrintingDefault(
   void
   );


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the user data of the message handler */
extern
SCIP_MESSAGEHDLRDATA* SCIPmessagehdlrGetData(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** returns the log file or NULL for stdout */
extern
FILE* SCIPmessagehdlrGetLogfile(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** returns TRUE if the message handler is set to be quiet */
extern
SCIP_Bool SCIPmessagehdlrIsQuiet(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPmessagehdlrGetData(messagehdlr)     ((messagehdlr) != NULL) ? messagehdlr->messagehdlrdata : NULL
#define SCIPmessagehdlrGetLogfile(messagehdlr)  ((messagehdlr) == NULL ? NULL : (messagehdlr)->logfile)
#define SCIPmessagehdlrIsQuiet(messagehdlr)     ((messagehdlr) == NULL || (messagehdlr)->quiet)

#endif

#ifdef __cplusplus
}
#endif

#endif

