/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This1 file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   message.c
 * @brief  message output methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "scip/type_message.h"
#include "scip/struct_message.h"
#include "scip/def.h"
#include "scip/pub_misc.h"
#include "blockmemshell/memory.h"


#ifndef va_copy
#define va_copy(dest, src) do { BMScopyMemory(&dest, &src); } while( 0 )
#endif

/* do defines for windows directly her to make the lpi more independent*/
#if defined(_WIN32) || defined(_WIN64)
#define snprintf _snprintf
#define vsnprintf _vsnprintf
#endif

/** stores the given message in the buffer and returns the part of the buffer and message that should be printed now */
static
SCIP_Bool bufferMessage(
   char*                 buffer,             /**< message buffer */
   int*                  bufferlen,          /**< pointer to the currently used entries in the message buffer */
   const char*           msg,                /**< message to store in the buffer; NULL to flush the buffer */
   char*                 outmsg,             /**< array to store message that should be printed immediately */
   int*                  outmsgsize          /**< pointer which stored allocated space in outmsg, size should be at
                                              *   least SCIP_MAXSTRLEN, if space would be not enough the needed sized is
                                              *   returned
                                              */
   )
{
   assert(outmsg != NULL);
   assert(outmsgsize != NULL);
   assert(*outmsgsize >= SCIP_MAXSTRLEN);

   *outmsg = '\0';

   /* should the buffer be flushed? */
   if( msg == NULL )
   {
      if( buffer != NULL )
         (void)strncpy(outmsg, buffer, SCIP_MAXSTRLEN);
      (*bufferlen) = 0;
      assert((int)strlen(outmsg) < SCIP_MAXSTRLEN);
      return TRUE;
   }

   /* do we have a buffer? */
   if( buffer == NULL )
   {
      /* check if the message fits into the output message */
      if( (int)strlen(msg) >= *outmsgsize )
      {
         *outmsgsize = (int)strlen(msg) + 1;
         return FALSE;
      }
      /* no buffer exists -> just copy the message to the output */
      (void)strncpy(outmsg, msg, strlen(msg));
      outmsg[strlen(msg)] = '\0';
      assert((int)strlen(outmsg) < *outmsgsize);
   }
   else
   {
      char* outmsgpos;

      assert(bufferlen != NULL);

      /* check if the message fits into the output message */
      if( (int)strlen(msg) + *bufferlen >= *outmsgsize + SCIP_MAXSTRLEN-1 )
      {
         *outmsgsize = (int)strlen(msg) + *bufferlen - SCIP_MAXSTRLEN + 2;
         return FALSE;
      }

      outmsgpos = outmsg;

      while( *msg != '\0' )
      {
         char c;

         assert(*bufferlen < SCIP_MAXSTRLEN-1);
         c = *msg;
         msg++;
         buffer[*bufferlen] = c;
         (*bufferlen)++;

         /* if the buffer is full or we reached a newline, move the buffer to the output message and store the
          * remaining message in the buffer
          */
         if( *bufferlen >= SCIP_MAXSTRLEN-1 || c == '\n' )
         {
            buffer[*bufferlen] = '\0';
            (void)strncpy(outmsgpos, buffer, SCIP_MAXSTRLEN);
            outmsgpos = &(outmsgpos[*bufferlen]);

            /* if rest fits into buffer copy all */
            if( (int)strlen(msg) < SCIP_MAXSTRLEN - 1 )
            {
               (void)strncpy(buffer, msg, SCIP_MAXSTRLEN);
               *bufferlen = (int)strlen(msg);
               assert(*bufferlen < SCIP_MAXSTRLEN-1);
               assert((int)strlen(outmsg) < *outmsgsize);
               break;
            }
            else
            {
               (void)strncpy(buffer, msg, SCIP_MAXSTRLEN - 2);
               *bufferlen = SCIP_MAXSTRLEN - 2;
               msg += (SCIP_MAXSTRLEN - 2);
            }
         }
      }
   }

   return TRUE;
}

/** default error printing method which is used to print all occurring errors */
static
SCIP_DECL_ERRORPRINTING(errorPrintingDefault)
{  /*lint --e{715}*/
   fputs(msg, stderr);
   fflush(stderr);
}

/** static variable which holds the error printing method */
static SCIP_DECL_ERRORPRINTING((*staticErrorPrinting)) = errorPrintingDefault;

/** static variable which holds a data pointer for the error prinint callback */
static void* staticErrorPrintingData = NULL;

/** prints error message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintError(
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   if( staticErrorPrinting != NULL )
      staticErrorPrinting(msg, msglength, staticErrorPrintingData);
}

/** prints warning message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   if( messagehdlr != NULL && messagehdlr->messagewarning != NULL && (!messagehdlr->quiet || messagehdlr->logfile != NULL) )
   {
      char* outmsg;
      int outmsgsize;

      outmsgsize = msglength + 1;

      if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
         return;

      if( !bufferMessage(messagehdlr->warningbuffer, &messagehdlr->warningbufferlen, msg, outmsg, &outmsgsize) )
      {
         assert(outmsgsize > msglength + 1);
         if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
         {
#ifndef NDEBUG
            SCIP_Bool ret;

            ret = bufferMessage(messagehdlr->warningbuffer, &messagehdlr->warningbufferlen, msg, outmsg, &outmsgsize);
            assert(ret);
#else
            bufferMessage(messagehdlr->warningbuffer, &messagehdlr->warningbufferlen, msg, outmsg, &outmsgsize);
#endif
         }
         else
         {
            BMSfreeMemory(&outmsg);
            return;
         }
      }

      if( *outmsg != '\0' )
      {
         if( !messagehdlr->quiet )
            messagehdlr->messagewarning(messagehdlr, stderr, outmsg);
         if( messagehdlr->logfile != NULL )
            messagehdlr->messagewarning(messagehdlr, messagehdlr->logfile, outmsg);
      }

      BMSfreeMemory(&outmsg);
   }
}

/** prints dialog message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{

   if( messagehdlr != NULL && messagehdlr->messagedialog != NULL )
   {
      if( (file == NULL || file == stdout) && !messagehdlr->quiet )
      {
         char* outmsg;
         int outmsgsize;

         outmsgsize = msglength + 1;

         if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
            return;

         if( !bufferMessage(messagehdlr->dialogbuffer, &messagehdlr->dialogbufferlen, msg, outmsg, &outmsgsize) )
         {
            assert(outmsgsize > msglength + 1);
            if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
            {
#ifndef NDEBUG
               SCIP_Bool ret;

               ret = bufferMessage(messagehdlr->dialogbuffer, &messagehdlr->dialogbufferlen, msg, outmsg, &outmsgsize);
               assert(ret);
#else
               bufferMessage(messagehdlr->dialogbuffer, &messagehdlr->dialogbufferlen, msg, outmsg, &outmsgsize);
#endif
            }
            else
            {
               BMSfreeMemory(&outmsg);
               return;
            }
         }

         if( *outmsg != '\0' )
         {
            messagehdlr->messagedialog(messagehdlr, stdout, outmsg);

            if( messagehdlr->logfile != NULL )
               messagehdlr->messagedialog(messagehdlr, messagehdlr->logfile, outmsg);
         }


         BMSfreeMemory(&outmsg);
      }
      else if( msg != NULL )
      {
         /* file output cannot be buffered because the output file may change */
         if( *msg != '\0' )
         {
            if( !messagehdlr->quiet || (file != NULL && file != stdout) )
            {
               assert(file != NULL);
               messagehdlr->messagedialog(messagehdlr, file, msg);
            }
            if( messagehdlr->logfile != NULL && (file == NULL || file == stdout) )
               messagehdlr->messagedialog(messagehdlr, messagehdlr->logfile, msg);
         }
      }
   }
}

/** prints info message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   if( messagehdlr != NULL && messagehdlr->messageinfo != NULL )
   {
      if( (file == NULL || file == stdout) && !messagehdlr->quiet && (msg == NULL || strlen(msg) < SCIP_MAXSTRLEN) )
      {
         char* outmsg;
         int outmsgsize;

         outmsgsize = msglength + 1;

         if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
            return;

         if( !bufferMessage(messagehdlr->infobuffer, &messagehdlr->infobufferlen, msg, outmsg, &outmsgsize) )
         {
            assert(outmsgsize > msglength + 1);
            if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
            {
#ifndef NDEBUG
               SCIP_Bool ret;

               ret = bufferMessage(messagehdlr->infobuffer, &messagehdlr->infobufferlen, msg, outmsg, &outmsgsize);
               assert(ret);
#else
               bufferMessage(messagehdlr->infobuffer, &messagehdlr->infobufferlen, msg, outmsg, &outmsgsize);
#endif
            }
            else
            {
               BMSfreeMemory(&outmsg);
               return;
            }
         }

         if( *outmsg != '\0' )
         {
            messagehdlr->messageinfo(messagehdlr, stdout, outmsg);

            if( messagehdlr->logfile != NULL )
               messagehdlr->messageinfo(messagehdlr, messagehdlr->logfile, outmsg);
         }

         BMSfreeMemory(&outmsg);
      }
      else if( msg != NULL )
      {
         /* file output cannot be buffered because the output file may change or the message is to long */
         if( *msg != '\0' )
         {
            if( !messagehdlr->quiet || (file != NULL && file != stdout) )
               messagehdlr->messageinfo(messagehdlr, file == NULL ? stdout : file, msg);
            if( messagehdlr->logfile != NULL && (file == NULL || file == stdout) )
               messagehdlr->messageinfo(messagehdlr, messagehdlr->logfile, msg);
         }
      }
   }
}

/** if the given file is not NULL a log file is opened */
static
void messagehdlrOpenLogfile(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           filename            /**< name of log file, or NULL (stdout) */
   )
{
   if( filename != NULL )
   {
      messagehdlr->logfile = fopen(filename, "a"); /* append to log file */

      if( messagehdlr->logfile == NULL )
      {
         SCIPerrorMessage("cannot open log file <%s> for writing\n", filename);
      }
   }
   else
      messagehdlr->logfile = NULL;
}

/** frees message handler */
static
SCIP_RETCODE messagehdlrFree(
   SCIP_MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   )
{
   assert(messagehdlr != NULL);

   if( *messagehdlr != NULL )
   {
      /* flush message buffers */
      messagePrintWarning(*messagehdlr, NULL, SCIP_MAXSTRLEN);
      messagePrintDialog(*messagehdlr, NULL, NULL, SCIP_MAXSTRLEN);
      messagePrintInfo(*messagehdlr, NULL, NULL, SCIP_MAXSTRLEN);

      if( (*messagehdlr)->messagehdlrfree != NULL )
      {
         /* call destructor method of message handler to free the message handler data */
         SCIP_CALL( (*messagehdlr)->messagehdlrfree(*messagehdlr) );
      }

      /* close the log file if one exists */
      if( (*messagehdlr)->logfile != NULL )
      {
         fclose((*messagehdlr)->logfile);
      }

      /* free buffer arrays */
      BMSfreeMemoryArrayNull(&(*messagehdlr)->warningbuffer);
      BMSfreeMemoryArrayNull(&(*messagehdlr)->dialogbuffer);
      BMSfreeMemoryArrayNull(&(*messagehdlr)->infobuffer);
      BMSfreeMemory(messagehdlr);
   }

   return SCIP_OKAY;
}

/** Creates a message handler which deals with warning, information, and dialog (interactive shell) methods.
 *
 *  @note The message handler does not handle error messages. For that see SCIPmessageSetErrorPrinting()
 *  @note Creating a message handler automatically captures it.
 */
SCIP_RETCODE SCIPmessagehdlrCreate(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   SCIP_Bool             bufferedoutput,     /**< should the output be buffered up to the next newline? */
   const char*           filename,           /**< name of log file, or NULL for no log */
   SCIP_Bool             quiet,              /**< should screen messages be suppressed? */
   SCIP_DECL_MESSAGEWARNING((*messagewarning)),/**< warning message print method of message handler */
   SCIP_DECL_MESSAGEDIALOG((*messagedialog)),/**< dialog message print method of message handler */
   SCIP_DECL_MESSAGEINFO ((*messageinfo)),   /**< info message print method of message handler */
   SCIP_DECL_MESSAGEHDLRFREE((*messagehdlrfree)), /**< destructor of message handler to free message handler data */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< message handler data */
   )
{
   SCIP_ALLOC( BMSallocMemory(messagehdlr) );
   (*messagehdlr)->messagewarning = messagewarning;
   (*messagehdlr)->messagedialog = messagedialog;
   (*messagehdlr)->messageinfo = messageinfo;
   (*messagehdlr)->messagehdlrfree = messagehdlrfree;
   (*messagehdlr)->messagehdlrdata = messagehdlrdata;
   (*messagehdlr)->warningbuffer = NULL;
   (*messagehdlr)->dialogbuffer = NULL;
   (*messagehdlr)->infobuffer = NULL;
   (*messagehdlr)->warningbufferlen = 0;
   (*messagehdlr)->dialogbufferlen = 0;
   (*messagehdlr)->infobufferlen = 0;
   (*messagehdlr)->nuses = 1;

   (*messagehdlr)->quiet = quiet;
   messagehdlrOpenLogfile(*messagehdlr, filename);

   /* allocate buffer for buffered output */
   if( bufferedoutput )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&(*messagehdlr)->warningbuffer, SCIP_MAXSTRLEN) ); /*lint !e506*/
      SCIP_ALLOC( BMSallocMemoryArray(&(*messagehdlr)->dialogbuffer, SCIP_MAXSTRLEN) ); /*lint !e506*/
      SCIP_ALLOC( BMSallocMemoryArray(&(*messagehdlr)->infobuffer, SCIP_MAXSTRLEN) ); /*lint !e506*/
      (*messagehdlr)->warningbuffer[0] = '\0';
      (*messagehdlr)->dialogbuffer[0] = '\0';
      (*messagehdlr)->infobuffer[0] = '\0';
   }

   return SCIP_OKAY;
}

/** captures message handler */
void SCIPmessagehdlrCapture(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler, or NULL */
   )
{
   if( messagehdlr != NULL )
      ++messagehdlr->nuses;
}

/** releases message handler */
SCIP_RETCODE SCIPmessagehdlrRelease(
   SCIP_MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   )
{
   assert(messagehdlr != NULL);

   if( *messagehdlr == NULL )
      return SCIP_OKAY;

   assert((*messagehdlr)->nuses >= 1);

   /* decrement usage counter */
   --(*messagehdlr)->nuses;

   /* the last one turns the light off */
   if( (*messagehdlr)->nuses == 0 )
   {
      SCIP_CALL( messagehdlrFree(messagehdlr) );
      assert(*messagehdlr == NULL);
   }
   else
   {
      *messagehdlr = NULL;
   }

   return SCIP_OKAY;
}

/** sets the user data of the message handler */
SCIP_RETCODE SCIPmessagehdlrSetData(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler; must not be NULL */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< new message handler data to attach to the handler */
   )
{
   if( messagehdlr == NULL )
      return SCIP_INVALIDDATA;

   messagehdlr->messagehdlrdata = messagehdlrdata;

   return SCIP_OKAY;
}

/** sets the log file name for the message handler */
void SCIPmessagehdlrSetLogfile(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           filename            /**< log file name where to copy messages into, or NULL */
   )
{
   /* close the old log file if one exists */
   if( messagehdlr->logfile != NULL )
   {
      fclose(messagehdlr->logfile);
   }

   /* opens the log file */
   messagehdlrOpenLogfile(messagehdlr, filename);
}

/** sets the messages handler to be quiet */
void SCIPmessagehdlrSetQuiet(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   )
{
   messagehdlr->quiet = quiet;
}

/** prints a warning message, acting like the printf() command */
void SCIPmessagePrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintWarning(messagehdlr, formatstr, ap);
   va_end(ap);
}

/** prints a warning message, acting like the vprintf() command */
void SCIPmessageVPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintWarning(messagehdlr, formatstr, ap);
}

/** prints a warning message, acting like the fprintf() command */
void SCIPmessageFPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintWarning(messagehdlr, formatstr, ap);
   va_end(ap);
}

/** prints a warning message, acting like the vfprintf() command */
void SCIPmessageVFPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   char msg[SCIP_MAXSTRLEN];
   int n;
   va_list aq;

   va_copy(aq, ap);

   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
      {
         va_end(aq);
         return;
      }

#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#endif
      assert(m == n);
      va_end(aq);
      messagePrintWarning(messagehdlr, bigmsg, n);
      BMSfreeMemory(&bigmsg);
      return;
   }

   messagePrintWarning(messagehdlr, msg, SCIP_MAXSTRLEN);
   va_end(aq);
}

/** prints a dialog message that requests user interaction, acting like the printf() command */
void SCIPmessagePrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintDialog(messagehdlr, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a dialog message that requests user interaction, acting like the vprintf() command */
void SCIPmessageVPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintDialog(messagehdlr, NULL, formatstr, ap);
}

/** prints a dialog message that requests user interaction into a file, acting like the fprintf() command */
void SCIPmessageFPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintDialog(messagehdlr, file, formatstr, ap);
   va_end(ap);
}

/** prints a dialog message that requests user interaction into a file, acting like the vfprintf() command */
void SCIPmessageVFPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   char msg[SCIP_MAXSTRLEN];
   int n;
   va_list aq;

   va_copy(aq, ap);

   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
      {
         va_end(aq);
         return;
      }

#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#endif
      assert(m == n);
      va_end(aq);
      messagePrintDialog(messagehdlr, file, bigmsg, n);
      BMSfreeMemory(&bigmsg);
      return;
   }
   messagePrintDialog(messagehdlr, file, msg, SCIP_MAXSTRLEN);
   va_end(aq);
}

/** prints a message, acting like the printf() command */
void SCIPmessagePrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintInfo(messagehdlr, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a message, acting like the vprintf() command */
void SCIPmessageVPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintInfo(messagehdlr, NULL, formatstr, ap);
}

/** prints a message into a file, acting like the fprintf() command */
void SCIPmessageFPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintInfo(messagehdlr, file, formatstr, ap);
   va_end(ap);
}

/** prints a message into a file, acting like the vfprintf() command */
void SCIPmessageVFPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   char msg[SCIP_MAXSTRLEN];
   int n;
   va_list aq;

   va_copy(aq, ap);

   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
      {
         va_end(aq);
         return;
      }

#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#endif
      assert(m == n);
      va_end(aq);
      messagePrintInfo(messagehdlr, file, bigmsg, n);
      BMSfreeMemory(&bigmsg);
      return;
   }
   messagePrintInfo(messagehdlr, file, msg, SCIP_MAXSTRLEN);
   va_end(aq);
}

/** prints a message depending on the verbosity level, acting like the printf() command */
void SCIPmessagePrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintVerbInfo(messagehdlr, verblevel, msgverblevel, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a message depending on the verbosity level, acting like the vprintf() command */
void SCIPmessageVPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintVerbInfo(messagehdlr, verblevel, msgverblevel, NULL, formatstr, ap);
}

/** prints a message into a file depending on the verbosity level, acting like the fprintf() command */
void SCIPmessageFPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintVerbInfo(messagehdlr, verblevel, msgverblevel, file, formatstr, ap);
   va_end(ap);
}

/** prints a message into a file depending on the verbosity level, acting like the vfprintf() command */
void SCIPmessageVFPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   assert(msgverblevel > SCIP_VERBLEVEL_NONE);
   assert(msgverblevel <= SCIP_VERBLEVEL_FULL);
   assert(verblevel <= SCIP_VERBLEVEL_FULL);

   if( msgverblevel <= verblevel )
   {
      char msg[SCIP_MAXSTRLEN];
      int n;
      va_list aq;

      va_copy(aq, ap);

      n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
      if( n < 0 )
         msg[SCIP_MAXSTRLEN-1] = '\0';
      else if( n >= SCIP_MAXSTRLEN )
      {
         char* bigmsg;
#ifndef NDEBUG
         int m;
#endif

         if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
         {
            va_end(aq);
            return;
         }

#ifndef NDEBUG
         m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#else
         vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#endif
         assert(m == n);
         va_end(aq);
         messagePrintInfo(messagehdlr, file, bigmsg, n);
         BMSfreeMemory(&bigmsg);
         return;
      }
      messagePrintInfo(messagehdlr, file, msg, SCIP_MAXSTRLEN);
      va_end(aq);
   }
}

/** prints the header with source file location for an error message using the static message handler */
void SCIPmessagePrintErrorHeader(
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline          /**< line in the source file where the function was called */
   )
{
   char msg[SCIP_MAXSTRLEN];

   /* safe string printing - do not use SCIPsnprintf() since message.c should be independent */
   (void) snprintf(msg, SCIP_MAXSTRLEN, "[%s:%d] ERROR: ", sourcefile, sourceline);
   msg[SCIP_MAXSTRLEN-1] = '\0';
   messagePrintError(msg, SCIP_MAXSTRLEN);
}

/** prints a error message, acting like the printf() command */
void SCIPmessagePrintError(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   char msg[SCIP_MAXSTRLEN];
   int n;
   va_list ap;

   va_start(ap, formatstr);  /*lint !e826*/

   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
      {
         va_end(ap);
         return;
      }

#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, ap);
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, ap);
#endif
      assert(m == n);
      va_end(ap);
      messagePrintError(bigmsg, n);
      BMSfreeMemory(&bigmsg);
      return;
   }

   messagePrintError(msg, SCIP_MAXSTRLEN);
   va_end(ap);
}

/** Method to set the error printing method. Setting the error printing method to NULL will suspend all error methods.
 *
 *  @note The error printing method is static variable. That means all occurring errors are handled via that methods
 */
void SCIPmessageSetErrorPrinting(
   SCIP_DECL_ERRORPRINTING((*errorPrinting)),/**< error message print method of message handler, or NULL */
   void*                 data                /**< data pointer which will be passed to the error printing method, or NULL */
   )
{
   staticErrorPrinting = errorPrinting;
   staticErrorPrintingData = data;
}

/** Method to set the error printing method to default version prints everything the stderr.
 *
 *  @note The error printing method is a static variable. This means that all occurring errors are handled via this method.
 */
void SCIPmessageSetErrorPrintingDefault(
   void
   )
{
   staticErrorPrinting = errorPrintingDefault;
   staticErrorPrintingData = NULL;
}

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPmessagehdlrGetData
#undef SCIPmessagehdlrGetLogfile
#undef SCIPmessagehdlrIsQuiet

/** returns the user data of the message handler */
SCIP_MESSAGEHDLRDATA* SCIPmessagehdlrGetData(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   if( messagehdlr != NULL )
      return messagehdlr->messagehdlrdata;
   else
      return NULL;
}


/** returns the log file or NULL for stdout */
FILE* SCIPmessagehdlrGetLogfile(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   if( messagehdlr == NULL )
      return NULL;

   return messagehdlr->logfile;
}

/** returns TRUE if the message handler is set to be quiet */
SCIP_Bool SCIPmessagehdlrIsQuiet(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   return (messagehdlr == NULL || messagehdlr->quiet);
}
