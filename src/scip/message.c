/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: message.c,v 1.14 2005/02/14 13:35:45 bzfpfend Exp $"

/**@file   message.c
 * @brief  message output methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>

#include "scip/message.h"


/** prints a message depending on the verbosity level, acting like the printf() command */
void infoMessage(
   VERBLEVEL        verblevel,          /**< current verbosity level */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert(msgverblevel > SCIP_VERBLEVEL_NONE);
   assert(msgverblevel <= SCIP_VERBLEVEL_FULL);
   assert(verblevel <= SCIP_VERBLEVEL_FULL);

   if( msgverblevel <= verblevel )
   {
      va_start(ap, formatstr);
      vprintf(formatstr, ap);
      va_end(ap);
   }
}

/** prints a message depending on the verbosity level, acting like the vprintf() command */
void vinfoMessage(
   VERBLEVEL        verblevel,          /**< current verbosity level */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*      formatstr,          /**< format string like in printf() function */
   va_list          ap                  /**< variable argument list */
   )
{
   assert(msgverblevel > SCIP_VERBLEVEL_NONE);
   assert(msgverblevel <= SCIP_VERBLEVEL_FULL);
   assert(verblevel <= SCIP_VERBLEVEL_FULL);

   if( msgverblevel <= verblevel )
      vprintf(formatstr, ap);
}

