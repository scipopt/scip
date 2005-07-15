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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: message.c,v 1.17 2005/07/15 17:20:11 bzfpfend Exp $"

/**@file   message.c
 * @brief  message output methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "scip/def.h"
#include "scip/memory.h"
#include "scip/message.h"


/** error message print method of default message handler */
static
DECL_MESSAGEERROR(messageErrorDefault)
{
   fputs(msg, file);
   fflush(file);
}

/** warning message print method of default message handler */
static
DECL_MESSAGEWARNING(messageWarningDefault)
{
   fputs(msg, file);
   fflush(file);
}

/** dialog message print method of default message handler */
static
DECL_MESSAGEDIALOG(messageDialogDefault)
{
   fputs(msg, file);
   fflush(file);
}

/** info message print method of default message handler */
static
DECL_MESSAGEINFO(messageInfoDefault)
{
   fputs(msg, file);
}

/** default message handler that prints messages to stdout or stderr */
static
MESSAGEHDLR messagehdlrDefault = {messageErrorDefault, messageWarningDefault, messageDialogDefault, messageInfoDefault,
                                  NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0};

/** static variable that contains the currently installed message handler;
 *  if the handler is set to NULL, messages are suppressed
 */
static MESSAGEHDLR* curmessagehdlr = &messagehdlrDefault;

/** stores the given message in the buffer and returns the part of the buffer and message that should be printed now */
static
void bufferMessage(
   char*            buffer,             /**< message buffer */
   int*             bufferlen,          /**< pointer to the currently used entries in the message buffer */
   const char*      msg,                /**< message to store in the buffer; NULL to flush the buffer */
   char*            outmsg              /**< array to store message that should be printed immediately */
   )
{
   assert(outmsg != NULL);
   assert(strlen(msg) < MAXSTRLEN);

   *outmsg = '\0';

   /* should the buffer be flushed? */
   if( msg == NULL )
   {
      strncpy(outmsg, buffer, MAXSTRLEN);
      (*bufferlen) = 0;
      assert(strlen(outmsg) < MAXSTRLEN);
      return;
   }

   /* do we have a buffer? */
   if( buffer == NULL )
   {
      /* no buffer exists -> just copy the message to the output */
      strncpy(outmsg, msg, MAXSTRLEN);
      assert(strlen(outmsg) < MAXSTRLEN);
   }
   else
   {
      assert(bufferlen != NULL);

      while( *msg != '\0' )
      {
         char c;

         assert(*bufferlen < MAXSTRLEN-1);
         c = *msg;
         msg++;
         buffer[*bufferlen] = c;
         (*bufferlen)++;

         /* if the buffer is full or we reached a newline, move the buffer to the output message and store the
          * remaining message in the buffer
          */
         if( *bufferlen >= MAXSTRLEN-1 || c == '\n' )
         {
            buffer[*bufferlen] = '\0';
            strncpy(outmsg, buffer, MAXSTRLEN);
            strncpy(buffer, msg, MAXSTRLEN);
            *bufferlen = strlen(msg);
            assert(*bufferlen < MAXSTRLEN-1);
            assert(strlen(outmsg) < MAXSTRLEN);
            break;
         }
      }
   }
}

/** prints error message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintError(
   const char*      msg                 /**< message to print; NULL to flush the output buffer */
   )
{
   if( curmessagehdlr != NULL && curmessagehdlr->messageerror != NULL )
   {
      char outmsg[MAXSTRLEN];

      bufferMessage(curmessagehdlr->errorbuffer, &curmessagehdlr->errorbufferlen, msg, outmsg);
      if( *outmsg != '\0' )
         curmessagehdlr->messageerror(curmessagehdlr, stderr, outmsg);
   }
}

/** prints warning message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintWarning(
   const char*      msg                 /**< message to print; NULL to flush the output buffer */
   )
{
   if( curmessagehdlr != NULL && curmessagehdlr->messagewarning != NULL )
   {
      char outmsg[MAXSTRLEN];

      bufferMessage(curmessagehdlr->warningbuffer, &curmessagehdlr->warningbufferlen, msg, outmsg);
      if( *outmsg != '\0' )
         curmessagehdlr->messagewarning(curmessagehdlr, stderr, outmsg);
   }
}

/** prints dialog message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintDialog(
   FILE*            file,               /**< file stream to print into, or NULL for stdout */
   const char*      msg                 /**< message to print; NULL to flush the output buffer */
   )
{
   if( curmessagehdlr != NULL && curmessagehdlr->messagedialog != NULL )
   {
      if( file == NULL || file == stdout )
      {
         char outmsg[MAXSTRLEN];
         
         bufferMessage(curmessagehdlr->dialogbuffer, &curmessagehdlr->dialogbufferlen, msg, outmsg);
         if( *outmsg != '\0' )
            curmessagehdlr->messagedialog(curmessagehdlr, stdout, outmsg);
      }
      else
      {
         /* file output cannot be buffered because the output file may change */
         if( *msg != '\0' )
            curmessagehdlr->messagedialog(curmessagehdlr, file, msg);
      }
   }
}

/** prints info message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintInfo(
   FILE*            file,               /**< file stream to print into, or NULL for stdout */
   const char*      msg                 /**< message to print; NULL to flush the output buffer */
   )
{
   if( curmessagehdlr != NULL && curmessagehdlr->messageinfo != NULL )
   {
      if( file == NULL || file == stdout )
      {
         char outmsg[MAXSTRLEN];
         
         bufferMessage(curmessagehdlr->infobuffer, &curmessagehdlr->infobufferlen, msg, outmsg);
         if( *outmsg != '\0' )
            curmessagehdlr->messageinfo(curmessagehdlr, stdout, outmsg);
      }
      else
      {
         /* file output cannot be buffered because the output file may change */
         if( *msg != '\0' )
            curmessagehdlr->messageinfo(curmessagehdlr, file, msg);
      }
   }
}

/** creates a message handler */
RETCODE SCIPmessagehdlrCreate(
   MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   Bool             bufferedoutput,     /**< should the output be buffered up to the next newline? */
   DECL_MESSAGEERROR((*messageerror)),  /**< error message print method of message handler */
   DECL_MESSAGEWARNING((*messagewarning)),/**< warning message print method of message handler */
   DECL_MESSAGEDIALOG((*messagedialog)),/**< dialog message print method of message handler */
   DECL_MESSAGEINFO ((*messageinfo)),   /**< info message print method of message handler */
   MESSAGEHDLRDATA* messagehdlrdata     /**< message handler data */
   )
{
   ALLOC_OKAY( allocMemory(messagehdlr) );
   (*messagehdlr)->messageerror = messageerror;
   (*messagehdlr)->messagewarning = messagewarning;
   (*messagehdlr)->messagedialog = messagedialog;
   (*messagehdlr)->messageinfo = messageinfo;
   (*messagehdlr)->messagehdlrdata = messagehdlrdata;
   (*messagehdlr)->errorbuffer = NULL;
   (*messagehdlr)->warningbuffer = NULL;
   (*messagehdlr)->dialogbuffer = NULL;
   (*messagehdlr)->infobuffer = NULL;
   (*messagehdlr)->errorbufferlen = 0;
   (*messagehdlr)->warningbufferlen = 0;
   (*messagehdlr)->dialogbufferlen = 0;
   (*messagehdlr)->infobufferlen = 0;

   /* allocate buffer for buffered output */
   if( bufferedoutput )
   {
      ALLOC_OKAY( allocMemoryArray(&(*messagehdlr)->errorbuffer, MAXSTRLEN) );
      ALLOC_OKAY( allocMemoryArray(&(*messagehdlr)->warningbuffer, MAXSTRLEN) );
      ALLOC_OKAY( allocMemoryArray(&(*messagehdlr)->dialogbuffer, MAXSTRLEN) );
      ALLOC_OKAY( allocMemoryArray(&(*messagehdlr)->infobuffer, MAXSTRLEN) );
      (*messagehdlr)->errorbuffer[0] = '\0';
      (*messagehdlr)->warningbuffer[0] = '\0';
      (*messagehdlr)->dialogbuffer[0] = '\0';
      (*messagehdlr)->infobuffer[0] = '\0';
   }

   return SCIP_OKAY;
}

/** frees message handler */
void SCIPmessagehdlrFree(
   MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   )
{
   assert(messagehdlr != NULL);

   if( *messagehdlr != NULL )
   {
      /* flush message buffers */
      messagePrintError(NULL);
      messagePrintWarning(NULL);
      messagePrintDialog(NULL, NULL);
      messagePrintInfo(NULL, NULL);

      freeMemoryArrayNull(&(*messagehdlr)->errorbuffer);
      freeMemoryArrayNull(&(*messagehdlr)->warningbuffer);
      freeMemoryArrayNull(&(*messagehdlr)->dialogbuffer);
      freeMemoryArrayNull(&(*messagehdlr)->infobuffer);
      freeMemory(messagehdlr);
   }
}

/** sets the user data of the message handler */
RETCODE SCIPmessagehdlrSetData(
   MESSAGEHDLR*     messagehdlr,        /**< message handler; must not be NULL */
   MESSAGEHDLRDATA* messagehdlrdata     /**< new message handler data to attach to the handler */
   )
{
   if( messagehdlr == NULL )
      return SCIP_INVALIDDATA;

   messagehdlr->messagehdlrdata = messagehdlrdata;

   return SCIP_OKAY;
}

/** returns the user data of the message handler */
MESSAGEHDLRDATA* SCIPmessagehdlrGetData(
   MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   if( messagehdlr != NULL )
      return messagehdlr->messagehdlrdata;
   else
      return NULL;
}

/** installs the given message handler */
void SCIPmessageSetHandler(
   MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   )
{
   curmessagehdlr = messagehdlr;
}

/** installs the default message handler that prints messages to stdout or stderr */
void SCIPmessageSetDefaultHandler(
   void
   )
{
   curmessagehdlr = &messagehdlrDefault;
}

/** returns the currently installed message handler, or NULL if messages are currently suppressed */
MESSAGEHDLR* SCIPmessageGetHandler(
   void
   )
{
   return curmessagehdlr;
}

/** prints the header with source file location for an error message */
void SCIPmessagePrintErrorHeader(
   const char*      sourcefile,         /**< name of the source file that called the function */
   int              sourceline          /**< line in the source file where the function was called */
   )
{
   char msg[MAXSTRLEN];

   snprintf(msg, MAXSTRLEN, "[%s:%d] ERROR: ", sourcefile, sourceline);
   messagePrintError(msg);
}

/** prints an error message, acting like the printf() command */
void SCIPmessagePrintError(
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   )
{
   char msg[MAXSTRLEN];
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   vsnprintf(msg, MAXSTRLEN, formatstr, ap);
   messagePrintError(msg);
   va_end(ap);
}

/** prints the header with source file location for an error message */
void SCIPmessagePrintWarningHeader(
   const char*      sourcefile,         /**< name of the source file that called the function */
   int              sourceline          /**< line in the source file where the function was called */
   )
{
   char msg[MAXSTRLEN];

   snprintf(msg, MAXSTRLEN, "[%s:%d] Warning: ", sourcefile, sourceline);
   messagePrintWarning(msg);
}

/** prints a warning message, acting like the printf() command */
void SCIPmessagePrintWarning(
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   )
{
   char msg[MAXSTRLEN];
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   vsnprintf(msg, MAXSTRLEN, formatstr, ap);
   messagePrintWarning(msg);
   va_end(ap);
}

/** prints a dialog message that requests user interaction, acting like the printf() command */
void SCIPmessagePrintDialog(
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintDialog(NULL, formatstr, ap);
   va_end(ap);
}

/** prints a dialog message that requests user interaction, acting like the vprintf() command */
void SCIPmessageVPrintDialog(
   const char*      formatstr,          /**< format string like in printf() function */
   va_list          ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintDialog(NULL, formatstr, ap);
}

/** prints a dialog message that requests user interaction into a file, acting like the fprintf() command */
void SCIPmessageFPrintDialog(
   FILE*            file,               /**< file stream to print into, or NULL for stdout */
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintDialog(file, formatstr, ap);
   va_end(ap);
}

/** prints a dialog message that requests user interaction into a file, acting like the vfprintf() command */
void SCIPmessageVFPrintDialog(
   FILE*            file,               /**< file stream to print into, or NULL for stdout */
   const char*      formatstr,          /**< format string like in printf() function */
   va_list          ap                  /**< variable argument list */
   )
{
   char msg[MAXSTRLEN];

   vsnprintf(msg, MAXSTRLEN, formatstr, ap);
   messagePrintDialog(file, msg);
}

/** prints a message, acting like the printf() command */
void SCIPmessagePrintInfo(
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintInfo(NULL, formatstr, ap);
   va_end(ap);
}

/** prints a message, acting like the vprintf() command */
void SCIPmessageVPrintInfo(
   const char*      formatstr,          /**< format string like in printf() function */
   va_list          ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintInfo(NULL, formatstr, ap);
}

/** prints a message into a file, acting like the fprintf() command */
void SCIPmessageFPrintInfo(
   FILE*            file,               /**< file stream to print into, or NULL for stdout */
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintInfo(file, formatstr, ap);
   va_end(ap);
}

/** prints a message into a file, acting like the vfprintf() command */
void SCIPmessageVFPrintInfo(
   FILE*            file,               /**< file stream to print into, or NULL for stdout */
   const char*      formatstr,          /**< format string like in printf() function */
   va_list          ap                  /**< variable argument list */
   )
{
   char msg[MAXSTRLEN];
   
   vsnprintf(msg, MAXSTRLEN, formatstr, ap);
   messagePrintInfo(file, msg);
}

/** prints a message depending on the verbosity level, acting like the printf() command */
void SCIPmessagePrintVerbInfo(
   VERBLEVEL        verblevel,          /**< current verbosity level */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintVerbInfo(verblevel, msgverblevel, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a message depending on the verbosity level, acting like the vprintf() command */
void SCIPmessageVPrintVerbInfo(
   VERBLEVEL        verblevel,          /**< current verbosity level */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*      formatstr,          /**< format string like in printf() function */
   va_list          ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintVerbInfo(verblevel, msgverblevel, NULL, formatstr, ap);
}

/** prints a message into a file depending on the verbosity level, acting like the fprintf() command */
void SCIPmessageFPrintVerbInfo(
   VERBLEVEL        verblevel,          /**< current verbosity level */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*            file,               /**< file stream to print into, or NULL for stdout */
   const char*      formatstr,          /**< format string like in printf() function */
   ...                                  /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintVerbInfo(verblevel, msgverblevel, file, formatstr, ap);
   va_end(ap);
}

/** prints a message into a file depending on the verbosity level, acting like the vfprintf() command */
void SCIPmessageVFPrintVerbInfo(
   VERBLEVEL        verblevel,          /**< current verbosity level */
   VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*            file,               /**< file stream to print into, or NULL for stdout */
   const char*      formatstr,          /**< format string like in printf() function */
   va_list          ap                  /**< variable argument list */
   )
{
   assert(msgverblevel > SCIP_VERBLEVEL_NONE);
   assert(msgverblevel <= SCIP_VERBLEVEL_FULL);
   assert(verblevel <= SCIP_VERBLEVEL_FULL);

   if( msgverblevel <= verblevel )
   {
      char msg[MAXSTRLEN];

      vsnprintf(msg, MAXSTRLEN, formatstr, ap);
      messagePrintInfo(file, msg);
   }
}

