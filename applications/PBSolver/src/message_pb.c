/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   message_pb.c
 * @brief  messagehdlr for the Pseudo-Boolean output format
 * @author Alexander Hoen
 * @author Gioni Mexi
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>

#include "scip/scip.h"
#include "message_pb.h"


static
void printMessage(
   SCIP_MESSAGEHDLR* messagehdlr,            /**< message handler */
   FILE* file,                               /**< file stream to print message into */
   const char* initial,                      /**< initial to prepend or NULL for no initial */
   const char* msg                           /**< message to print or NULL to flush only */
   )
{
   if( msg != NULL && msg[0] != '\0' )
   {
#ifdef PBSOLVER
      SCIP_MESSAGEHDLRDATA* messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);

      if( messagehdlrdata == NULL || messagehdlrdata->comment )
         fputs("c ", file);
#endif

      if( initial != NULL )
         fputs(initial, file);

      fputs(msg, file);
   }

   fflush(file);
}

/** default error printing method which is used to print all occurring errors */
static
SCIP_DECL_ERRORPRINTING(messageErrorPbSolver)
{ /*lint --e{715}*/
   printMessage(NULL, stderr, NULL, msg);
}

/** warning message print method of default message handler */
static
SCIP_DECL_MESSAGEWARNING(messageWarningPbSolver)
{ /*lint --e{715}*/
   printMessage(messagehdlr, file, "WARNING: ", msg);
}

/** dialog message print method of default message handler */
static
SCIP_DECL_MESSAGEDIALOG(messageDialogPbSolver)
{ /*lint --e{715}*/
   printMessage(messagehdlr, file, NULL, msg);
}

/** info message print method of default message handler */
static
SCIP_DECL_MESSAGEINFO(messageInfoPbSolver)
{ /*lint --e{715}*/
   printMessage(messagehdlr, file, NULL, msg);
}

static
SCIP_DECL_MESSAGEHDLRFREE(messagehdlrFreePbSolver)
{ /*lint --e{715}*/
   assert(messagehdlr != NULL);
   SCIP_MESSAGEHDLRDATA* messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);

   BMSfreeMemory(&messagehdlrdata);

   return SCIP_OKAY;
}

/** creates default pbsolver message handler */
SCIP_RETCODE SCIPcreateMessagehdlrPbSolver(
   SCIP_MESSAGEHDLR** messagehdlr,           /**< pointer to message handler */
   SCIP_Bool buffered,                       /**< should the output be buffered */
   const char* filename,                     /**< name of log file or NULL for stdout */
   SCIP_Bool quiet                           /**< should screen messages be suppressed? */
   )
{
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;

   SCIP_ALLOC( BMSallocMemory(&messagehdlrdata) );

#ifdef PBSOLVER
   messagehdlrdata->comment = TRUE;
#else
   messagehdlrdata->comment = FALSE;
#endif

   SCIP_CALL( SCIPmessagehdlrCreate(messagehdlr, buffered, filename, quiet, messageWarningPbSolver,
         messageDialogPbSolver, messageInfoPbSolver, messagehdlrFreePbSolver, messagehdlrdata) );

   SCIPmessageSetErrorPrinting(messageErrorPbSolver, NULL);

   return SCIP_OKAY;
}
