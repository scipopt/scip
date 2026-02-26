/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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
 * @author Dominik Kamp
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
      SCIP_MESSAGEHDLRDATA* messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);

      if( messagehdlrdata == NULL || messagehdlrdata->comment )
         fputs("c ", file);

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

   messagehdlrdata->comment = TRUE;

   SCIP_CALL( SCIPmessagehdlrCreate(messagehdlr, buffered, filename, quiet, messageWarningPbSolver,
         messageDialogPbSolver, messageInfoPbSolver, messagehdlrFreePbSolver, messagehdlrdata) );

   SCIPmessageSetErrorPrinting(messageErrorPbSolver, NULL);

   return SCIP_OKAY;
}

/** prints that the problem instance is unsupported */
SCIP_RETCODE SCIPprintUnsupportedPbSolver(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr = SCIPgetMessagehdlr(scip);
   assert(messagehdlr != NULL);
   SCIP_MESSAGEHDLRDATA* messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);
   SCIP_Bool quiet = SCIPmessagehdlrIsQuiet(messagehdlr);

   /* competition console log */
   if( messagehdlrdata->comment )
   {
      /* turn on message handler and remove comment declaration */
      SCIPmessagehdlrSetQuiet(messagehdlr, FALSE);
      messagehdlrdata->comment = FALSE;

      SCIPinfoMessage(scip, NULL, "s UNSUPPORTED\n");

      /* reset old values of message handler data */
      SCIPmessagehdlrSetQuiet(messagehdlr, quiet);
      messagehdlrdata->comment = TRUE;
   }

   return SCIP_OKAY;
}

/** prints the best primal solution */
SCIP_RETCODE SCIPprintSolutionPbSolver(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr = SCIPgetMessagehdlr(scip);
   assert(messagehdlr != NULL);
   SCIP_MESSAGEHDLRDATA* messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);
   SCIP_SOL* bestsol;
   SCIP_VAR** vars;
   SCIP_STATUS status;
   SCIP_Bool quiet = SCIPmessagehdlrIsQuiet(messagehdlr);
   SCIP_Bool hasobj;
   SCIP_Bool feasible;
   char* solutiontext;
   int strlength;
   int printlength;
   int nvars;
   int n;
   int v;

   /* enable message handler */
   SCIPmessagehdlrSetQuiet(messagehdlr, FALSE);

   /* competition console log */
   if( messagehdlrdata->comment )
   {
      /* disable comment prefix */
      messagehdlrdata->comment = FALSE;

      status = SCIPgetStatus(scip);
      if( status == SCIP_STATUS_INFEASIBLE )
         SCIPinfoMessage(scip, NULL, "s UNSATISFIABLE\n");
      else
      {
         feasible = FALSE;
         bestsol = SCIPgetBestSol(scip);

         if( bestsol != NULL )
         {
            /* check if solution is feasible in original problem */
            SCIP_CALL( SCIPcheckSolOrig(scip, bestsol, &feasible, FALSE, FALSE) );

            /* if solution is not feasible in original problem; post UNKNOWN */
            if( !feasible )
            {
               SCIPinfoMessage(scip, NULL, "c internal error\n");
               SCIPinfoMessage(scip, NULL, "s UNKNOWN\n");
            }
            else
            {
               vars = SCIPgetOrigVars(scip);
               nvars = SCIPgetNOrigVars(scip);
               assert(vars != NULL || nvars == 0);
               hasobj = FALSE;

               if( status == SCIP_STATUS_OPTIMAL )
               {
                  for( v = 0; v < nvars; ++v )
                  {
                     if( !SCIPisZero(scip, SCIPvarGetObj(vars[v])) )
                     {
                        hasobj = TRUE;
                        break;
                     }
                  }
               }

               if( !hasobj )
                  SCIPinfoMessage(scip, NULL, "s SATISFIABLE\n");
               else
                  SCIPinfoMessage(scip, NULL, "s OPTIMUM FOUND\n");

               strlength = 2 * SCIP_MAXSTRLEN;
               printlength = SCIP_MAXSTRLEN / 8;
               n = 0;
               v = 0;
               SCIP_CALL( SCIPallocBufferArray(scip, &solutiontext, strlength) );

               while( v < nvars )
               {
                  int printed = 0;

                  assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_ORIGINAL);

                  if( strstr(SCIPvarGetName(vars[v]), "andresultant_") == NULL
                     && strstr(SCIPvarGetName(vars[v]), "indslack_") == NULL
                     && strstr(SCIPvarGetName(vars[v]), "indicatorvar_") == NULL )
                  {
                     printed = SCIPsnprintf(solutiontext + n, strlength,
                           SCIPgetSolVal(scip, bestsol, vars[v]) > 0.5 ? " %s" : " -%s", SCIPvarGetName(vars[v]));
                     n += printed;
                     strlength -= printed;
                  }

                  if( n >= printlength || ( n >= 1 && v + 1 >= nvars ) )
                  {
                     if( strlength >= 1 )
                     {
                        strlength += n;
                        n = 0;
                        SCIPinfoMessage(scip, NULL, "v%s\n", solutiontext);
                     }
                     else
                     {
                        strlength += printed + SCIP_MAXSTRLEN;
                        n -= printed;
                        SCIP_CALL( SCIPreallocBufferArray(scip, &solutiontext, n + strlength) );

                        continue;
                     }
                  }
                  assert(strlength >= 1);

                  ++v;
               }

               SCIPfreeBufferArray(scip, &solutiontext);
            }
         }
         else
            SCIPinfoMessage(scip, NULL, "s UNKNOWN\n");
      }

      /* enable comment prefix */
      messagehdlrdata->comment = TRUE;
   }
   /* default console log */
   else
   {
      SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
   }

   /* reset message handler */
   SCIPmessagehdlrSetQuiet(messagehdlr, quiet);

   return SCIP_OKAY;
}
