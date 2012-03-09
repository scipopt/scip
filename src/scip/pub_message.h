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

#ifdef __cplusplus
extern "C" {
#endif

/** prints an error message */
#define SCIPerrorMessage                SCIPmessagePrintErrorHeader(__FILE__, __LINE__); \
                                        SCIPmessagePrintError

/** prints a warning message */
#define SCIPwarningMessage              SCIPmessagePrintWarningHeader(__FILE__, __LINE__); \
                                        SCIPmessagePrintWarning

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


/** sets the user data of the message handler */
extern
SCIP_RETCODE SCIPmessagehdlrSetData(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler; must not be NULL */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< new message handler data to attach to the handler */
   );

/** returns the user data of the message handler */
extern
SCIP_MESSAGEHDLRDATA* SCIPmessagehdlrGetData(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** prints the header with source file location for an error message */
extern
void SCIPmessagePrintErrorHeader(
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline          /**< line in the source file where the function was called */
   );

/** prints an error message, acting like the printf() command */
extern
void SCIPmessagePrintError(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints the header with source file location for an error message */
extern
void SCIPmessagePrintWarningHeader(
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline          /**< line in the source file where the function was called */
   );

/** prints a warning message, acting like the printf() command */
extern
void SCIPmessagePrintWarning(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

#ifdef __cplusplus
}
#endif

#endif

