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

/**@file   message.h
 * @brief  message output methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MESSAGE_H__
#define __MESSAGE_H__


/** verbosity levels of output */
enum VerbLevel
{
   SCIP_VERBLEVEL_NONE   = 0,           /**< only error messages are displayed */
   SCIP_VERBLEVEL_NORMAL = 1,           /**< standard messages are displayed */
   SCIP_VERBLEVEL_HIGH   = 2,           /**< a lot of information is displayed */
   SCIP_VERBLEVEL_FULL   = 3            /**< all messages are displayed */
};
typedef enum VerbLevel VERBLEVEL;


#define errorMessage(msg)               errorMessage_call((msg), __FILE__, __LINE__)
#define todoMessage(msg)                todoMessage_call((msg), __FILE__, __LINE__)

#ifdef DEBUG
#define debugMessage                    printf("[%s:%d] debug: ", __FILE__, __LINE__); printf
#else
#define debugMessage                    if( FALSE ) printf
#endif


extern
void errorMessage_call(                 /**< prints an error message */
   const char*      msg,                /**< message to print */
   const char*      filename,           /**< name of the file, where the error occured */
   int              line                /**< line of the file, where the error occured */
   );

extern
void todoMessage_call(                  /**< prints a todo message */
   const char*      msg,                /**< message to print */
   const char*      filename,           /**< name of the file, where the error occured */
   int              line                /**< line of the file, where the error occured */
   );

extern
void infoMessage(                       /**< prints a message depending on the verbosity level */
   VERBLEVEL        verblevel,          /**< actual verbosity level */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*      msg                 /**< message to print */
   );


#endif
