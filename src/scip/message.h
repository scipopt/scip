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
#pragma ident "@(#) $Id: message.h,v 1.16 2004/07/20 14:34:41 bzfpfend Exp $"

/**@file   message.h
 * @brief  message output methods
 * @author Tobias Achterberg
 *
 * Because the message functions are implemented as defines with more than one
 * function call, they shouldn't be used as a single statement like in:
 *   if( error )
 *      errorMessage("an error occured");
 * because this would produce the following macro extension:
 *   if( error )
 *      printf(("[%s:%d] ERROR: ", __FILE__, __LINE__);
 *   printf(("an error occured");
 * Instead, they should be protected with brackets:
 *   if( error )
 *   {
 *      errorMessage("an error occured");
 *   }
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MESSAGE_H__
#define __MESSAGE_H__


/** verbosity levels of output */
enum VerbLevel
{
   SCIP_VERBLEVEL_NONE    = 0,           /**< only error and warning messages are displayed */
   SCIP_VERBLEVEL_MINIMAL = 1,           /**< a reduced number of messages are displayed */
   SCIP_VERBLEVEL_NORMAL  = 2,           /**< standard messages are displayed */
   SCIP_VERBLEVEL_HIGH    = 3,           /**< a lot of information is displayed */
   SCIP_VERBLEVEL_FULL    = 4            /**< all messages are displayed */
};
typedef enum VerbLevel VERBLEVEL;


#define errorMessage                    printf("[%s:%d] ERROR: ", __FILE__, __LINE__); printf
#define warningMessage                  printf("[%s:%d] Warning: ", __FILE__, __LINE__); printf

#ifdef DEBUG
#define debug(x)                        x
#define debugMessage                    printf("[%s:%d] debug: ", __FILE__, __LINE__); printf
#else
#define debug(x)                        /**/
#define debugMessage                    while( FALSE ) printf
#endif


#include <stdio.h>
#include <stdarg.h>



/** prints a message depending on the verbosity level */
extern
void infoMessage(
   VERBLEVEL        verblevel,          /**< current verbosity level */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   );

/** prints a message depending on the verbosity level, acting like the vprintf() command */
extern
void vinfoMessage(
   VERBLEVEL        verblevel,          /**< current verbosity level */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*      formatstr,          /**< format string like in printf() function */
   va_list          ap                  /**< variable argument list */
   );


#endif

