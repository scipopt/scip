/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dialog.c
 * @brief  user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "dialog.h"
#include "memory.h"

#ifdef WITH_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif



/** user interface dialog */
struct Dialog
{
   DECL_DIALOGEXEC  ((*dialogexec));    /**< execution method of dialog */
   DECL_DIALOGDESC  ((*dialogdesc));    /**< description output method of dialog, or NULL */
   char*            name;               /**< name of dialog: command name appearing in parent's dialog menu */
   char*            desc;               /**< description of dialog used if description output method is NULL */
   Bool             issubmenu;          /**< is the dialog a submenu? */
   DIALOG*          parent;             /**< parent dialog of dialog */
   DIALOG**         subdialogs;         /**< sub dialogs of dialog */
   int              nsubdialogs;        /**< number of sub dialogs */
   int              subdialogssize;     /**< size of subdialogs array */
   int              nuses;              /**< number of times, the dialog is used */
   DIALOGDATA*      dialogdata;         /**< user defined dialog data */
};

/** dialog handler */
struct Dialoghdlr
{
   DIALOG*          rootdialog;         /**< main (root) dialog */
   char*            buffer;             /**< command buffer */
   int              buffersize;         /**< size of command buffer */
   int              bufferpos;          /**< position of first unprocessed character in buffer */
};




/*
 * read line methods
 */

#ifdef WITH_READLINE

/** reads a line of input from stdin */
static
RETCODE readLine(
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   const char*      prompt              /**< prompt to display */
   )
{
   char* s;

   s = readline(prompt);
   (void)strncpy(&dialoghdlr->buffer[dialoghdlr->bufferpos], s, dialoghdlr->buffersize - dialoghdlr->bufferpos);
   free(s);

   return SCIP_OKAY;
}

/** puts the given string on the command history */
static
RETCODE addHistory(
   const char*      s                   /**< string to add to the command history */
   )
{
   add_history(s);
   
   return SCIP_OKAY;
}

#else

/** reads a line of input from stdin */
static
RETCODE readLine(
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   const char*      prompt              /**< prompt to display */
   )
{
   char* s;

   assert(dialoghdlr != NULL);
   assert(dialoghdlr->buffer != NULL);
   assert(dialoghdlr->bufferpos < dialoghdlr->buffersize);
   assert(dialoghdlr->buffer[dialoghdlr->bufferpos] == '\0');

   /* display prompt */
   printf(prompt);

   /* read line from stdin */
   (void)fgets(&dialoghdlr->buffer[dialoghdlr->bufferpos], dialoghdlr->buffersize - dialoghdlr->bufferpos, stdin);

   /* replace newline with \0 */
   s = strchr(&dialoghdlr->buffer[dialoghdlr->bufferpos], '\n');
   if( s != NULL )
      *s = '\0';
}

/** puts the given string on the command history */
static
RETCODE addHistory(
   const char*      s                   /**< string to add to the command history */
   )
{
   /* nothing to do here */
}

#endif




/*
 * dialog handler
 */

/** creates a dialog handler */
RETCODE SCIPdialoghdlrCreate(
   DIALOGHDLR**     dialoghdlr          /**< pointer to store dialog handler */
   )
{
   assert(dialoghdlr != NULL);
   
   ALLOC_OKAY( allocMemory(dialoghdlr) );
   (*dialoghdlr)->rootdialog = NULL;
   (*dialoghdlr)->buffersize = MAXSTRLEN;
   ALLOC_OKAY( allocMemoryArray(&(*dialoghdlr)->buffer, (*dialoghdlr)->buffersize) );

   SCIPdialoghdlrClearBuffer(*dialoghdlr);

   return SCIP_OKAY;
}

/** frees a dialog handler and it's dialog tree */
RETCODE SCIPdialoghdlrFree(
   DIALOGHDLR**     dialoghdlr          /**< pointer to dialog handler */
   )
{
   assert(dialoghdlr != NULL);
   
   CHECK_OKAY( SCIPdialoghdlrSetRoot(*dialoghdlr, NULL) );
   freeMemoryArray(&(*dialoghdlr)->buffer);
   freeMemory(dialoghdlr);

   return SCIP_OKAY;
}

/** executes the root dialog of the dialog handler */
RETCODE SCIPdialoghdlrExec(
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   DIALOG* dialog;

   assert(dialoghdlr != NULL);
   assert(dialoghdlr->buffer != NULL);

   /* clear the buffer, start with the root dialog */
   SCIPdialoghdlrClearBuffer(dialoghdlr);
   dialog = dialoghdlr->rootdialog;

   /* execute dialogs until a NULL is returned as nextdialog */
   while( dialog != NULL )
   {
      CHECK_OKAY( SCIPdialogExec(dialog, scip, dialoghdlr, &dialog) );
      
      /* reset buffer, it is was consumed completely */
      if( dialoghdlr->buffer[dialoghdlr->bufferpos] == '\0' )
         SCIPdialoghdlrClearBuffer(dialoghdlr);
   }

   return SCIP_OKAY;
}

/** makes given dialog the root dialog of dialog handler; captures dialog and releases former root dialog */
RETCODE SCIPdialoghdlrSetRoot(
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   DIALOG*          dialog              /**< dialog to be the root */
   )
{
   assert(dialoghdlr != NULL);

   if( dialoghdlr->rootdialog != NULL )
   {
      CHECK_OKAY( SCIPdialogRelease(&dialoghdlr->rootdialog) );
   }
   assert(dialoghdlr->rootdialog == NULL);

   dialoghdlr->rootdialog = dialog;

   if( dialog != NULL )
      SCIPdialogCapture(dialog);

   return SCIP_OKAY;
}

/** returns the root dialog of the dialog handler */
DIALOG* SCIPdialoghdlrGetRoot(
   DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   )
{
   assert(dialoghdlr != NULL);

   return dialoghdlr->rootdialog;
}

/** clears the input command buffer of the dialog handler */
void SCIPdialoghdlrClearBuffer(
   DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   )
{
   assert(dialoghdlr != NULL);

   dialoghdlr->buffer[0] = '\0';
   dialoghdlr->bufferpos = 0;
}

/** returns TRUE iff input command buffer is empty */
Bool SCIPdialoghdlrIsBufferEmpty(
   DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   )
{
   assert(dialoghdlr != NULL);
   assert(dialoghdlr->bufferpos < dialoghdlr->buffersize);

   return (dialoghdlr->buffer[dialoghdlr->bufferpos] == '\0');
}

/** returns the next word in the handler's command buffer; if the buffer is empty, displays the given prompt or the 
 *  current dialog's path and asks the user for further input; the user must not free or modify the returned string
 */
const char* SCIPdialoghdlrGetWord(
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   DIALOG*          dialog,             /**< current dialog */
   const char*      prompt              /**< prompt to display, or NULL to display the current dialog's path */
   )
{
   char* firstword;

   assert(dialoghdlr != NULL);
   assert(dialoghdlr->buffer != NULL);
   assert(dialoghdlr->bufferpos < dialoghdlr->buffersize);

   /* get input from the user, if the buffer is empty */
   if( SCIPdialoghdlrIsBufferEmpty(dialoghdlr) )
   {
      char path[MAXSTRLEN];
      char p[MAXSTRLEN];

      /* clear the buffer */
      SCIPdialoghdlrClearBuffer(dialoghdlr);

      if( prompt == NULL )
      {
         /* use current dialog's path as prompt */
         SCIPdialogGetPath(dialog, '/', path);
         snprintf(p, MAXSTRLEN, "%s> ", path);
         p[MAXSTRLEN-1] = '\0';
         prompt = p;
      }

      /* read command line from stdin */
      CHECK_ABORT( readLine(dialoghdlr, prompt) );

      /* insert command in command history */
      if( dialoghdlr->buffer[dialoghdlr->bufferpos] != '\0' )
      {
         CHECK_ABORT( SCIPdialoghdlrAddHistory(dialoghdlr, NULL, &dialoghdlr->buffer[dialoghdlr->bufferpos]) );
      }
   }

   /* the last character in the buffer must be a '\0' */
   dialoghdlr->buffer[dialoghdlr->buffersize-1] = '\0';

   /* skip leading spaces: find start of first word */
   while( isspace(dialoghdlr->buffer[dialoghdlr->bufferpos]) )
      dialoghdlr->bufferpos++;
   firstword = &dialoghdlr->buffer[dialoghdlr->bufferpos];

   /* find the next space in the buffer */
   while( dialoghdlr->buffer[dialoghdlr->bufferpos] != '\0' && !isspace(dialoghdlr->buffer[dialoghdlr->bufferpos]) )
      dialoghdlr->bufferpos++;

   /* truncate the command word in the buffer */
   if( dialoghdlr->buffer[dialoghdlr->bufferpos] != '\0' )
   {
      dialoghdlr->buffer[dialoghdlr->bufferpos] = '\0';
      dialoghdlr->bufferpos++;
   }

   /* remove additional spaces */
   while( isspace(dialoghdlr->buffer[dialoghdlr->bufferpos]) )
      dialoghdlr->bufferpos++;

   return firstword;
}

/** adds a command to the command history of the dialog handler; if a dialog is given, the command is preceeded
 *  by the dialog's command path; if no command is given, only the path to the dialog is added to the command history
 */
RETCODE SCIPdialoghdlrAddHistory(
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   DIALOG*          dialog,             /**< current dialog, or NULL */
   const char*      command             /**< command string to add to the command history, or NULL */
   )
{
   char s[MAXSTRLEN];
   char h[MAXSTRLEN];

   assert(dialoghdlr != NULL);

   s[MAXSTRLEN-1] = '\0';
   h[MAXSTRLEN-1] = '\0';

   if( command != NULL )
      strncpy(h, command, MAXSTRLEN-1);
   else
      h[0] = '\0';

   while( dialog != NULL && dialog != dialoghdlr->rootdialog )
   {
      snprintf(s, MAXSTRLEN-1, "%s %s", dialog->name, h);
      (void)strncpy(h, s, MAXSTRLEN-1);
      dialog = dialog->parent;
   }

   if( h[0] != '\0' )
   {
      CHECK_OKAY( addHistory(h) );
   }

   return SCIP_OKAY;
}




/*
 * dialog
 */

/** ensures, that subdialogs array can store at least the given number of sub dialogs */
static
RETCODE ensureSubdialogMem(
   DIALOG*          dialog,             /**< dialog */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal storage size for sub dialogs */
   )
{
   assert(dialog != NULL);

   if( num > dialog->subdialogssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&(dialog->subdialogs), newsize) );
      dialog->subdialogssize = newsize;
   }
   assert(num <= dialog->subdialogssize);

   return SCIP_OKAY;
}

/** creates and captures a user interface dialog */
RETCODE SCIPdialogCreate(
   DIALOG**         dialog,             /**< pointer to store the dialog */
   DECL_DIALOGEXEC  ((*dialogexec)),    /**< execution method of dialog */
   DECL_DIALOGDESC  ((*dialogdesc)),    /**< description output method of dialog, or NULL */
   const char*      name,               /**< name of dialog: command name appearing in parent's dialog menu */
   const char*      desc,               /**< description of dialog used if description output method is NULL */
   Bool             issubmenu,          /**< is the dialog a submenu? */
   DIALOGDATA*      dialogdata          /**< user defined dialog data */
   )
{
   assert(dialog != NULL);
   assert(name != NULL);

   ALLOC_OKAY( allocMemory(dialog) );
   (*dialog)->dialogexec = dialogexec;
   (*dialog)->dialogdesc = dialogdesc;

   ALLOC_OKAY( duplicateMemoryArray(&(*dialog)->name, name, strlen(name)+1) );
   if( desc != NULL )
   {
      ALLOC_OKAY( duplicateMemoryArray(&(*dialog)->desc, desc, strlen(desc)+1) );
   }
   else
      (*dialog)->desc = NULL;

   (*dialog)->issubmenu = issubmenu;
   (*dialog)->parent = NULL;
   (*dialog)->subdialogs = NULL;
   (*dialog)->nsubdialogs = 0;
   (*dialog)->subdialogssize = 0;
   (*dialog)->nuses = 0;
   (*dialog)->dialogdata = dialogdata;

   /* capture dialog */
   SCIPdialogCapture(*dialog);

   return SCIP_OKAY;
}

/** frees dialog and all of its sub dialogs */
static
RETCODE dialogFree(
   DIALOG**         dialog              /**< pointer to dialog */
   )
{
   int i;

   assert(dialog != NULL);
   assert(*dialog != NULL);
   assert((*dialog)->nuses == 0);

   /** release sub dialogs */
   for( i = 0; i < (*dialog)->nsubdialogs; ++i )
   {
      CHECK_OKAY( SCIPdialogRelease(&(*dialog)->subdialogs[i]) );
   }
   freeMemoryArrayNull(&(*dialog)->subdialogs);

   freeMemoryArrayNull(&(*dialog)->name);
   freeMemoryArrayNull(&(*dialog)->desc);
   freeMemory(dialog);

   return SCIP_OKAY;
}

/** captures a dialog */
void SCIPdialogCapture(
   DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   dialog->nuses++;
}

/** releases a dialog */
RETCODE SCIPdialogRelease(
   DIALOG**         dialog              /**< pointer to dialog */
   )
{
   assert(dialog != NULL);

   (*dialog)->nuses--;
   if( (*dialog)->nuses == 0 )
   {
      CHECK_OKAY( dialogFree(dialog) );
   }
   
   return SCIP_OKAY;
}

/** executes dialog */
RETCODE SCIPdialogExec(
   DIALOG*          dialog,             /**< dialog */
   SCIP*            scip,               /**< SCIP data structure */   
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   DIALOG**         nextdialog          /**< pointer to store the next dialog to process */
   )
{
   assert(dialog != NULL);
   assert(dialog->dialogexec != NULL);
   assert(nextdialog != NULL);

   CHECK_OKAY( dialog->dialogexec(scip, dialog, dialoghdlr, nextdialog) );

   return SCIP_OKAY;
}

/** adds a sub dialog to the given dialog as menu entry and captures the sub dialog */
RETCODE SCIPdialogAddEntry(
   DIALOG*          dialog,             /**< dialog */
   const SET*       set,                /**< global SCIP settings */
   DIALOG*          subdialog           /**< subdialog to add as menu entry in dialog */
   )
{
   assert(dialog != NULL);
   assert(subdialog != NULL);

   /* check, if subdialog already exists */
   if( SCIPdialogHasEntry(dialog, SCIPdialogGetName(subdialog)) )
   {
      char s[MAXSTRLEN];
      sprintf(s, "dialog entry with name <%s> already exists in dialog <%s>\n",
         SCIPdialogGetName(subdialog), SCIPdialogGetName(dialog));
      errorMessage(s);
      return SCIP_INVALIDDATA;
   }

   /* resize the subdialogs array */
   CHECK_OKAY( ensureSubdialogMem(dialog, set, dialog->nsubdialogs+1) );

   /* link the dialogs as parent-child pair */
   dialog->subdialogs[dialog->nsubdialogs] = subdialog;
   dialog->nsubdialogs++;
   subdialog->parent = dialog;

   /* capture sub dialog */
   SCIPdialogCapture(subdialog);

   return SCIP_OKAY;
}

/** returns TRUE iff a dialog entry matching exactly the given name is existing in the given dialog */
Bool SCIPdialogHasEntry(
   DIALOG*          dialog,             /**< dialog */
   const char*      entryname           /**< name of the dialog entry to find */
   )
{
   DIALOG** subdialogs;
   int nsubdialogs;
   int i;

   assert(dialog != NULL);
   assert(entryname != NULL);

   /* check entryname w.r.t. available dialog options */
   subdialogs = SCIPdialogGetSubdialogs(dialog);
   nsubdialogs = SCIPdialogGetNSubdialogs(dialog);
   for( i = 0; i < nsubdialogs; ++i )
   {
      /* check, if the sub dialog's name matches entryname */
      if( strcmp(entryname, SCIPdialogGetName(subdialogs[i])) == 0 )
         return TRUE;
   }

   return FALSE;
}

/** searches the dialog for entries corresponding to the given name;
 *  If a complete match is found, the entry is returned as "subdialog" and
 *  the return value is 1.
 *  If no dialog entry completely matches the given "entryname", the number
 *  of entries with names beginning with "entryname" is returned. If this
 *  number is 1, the single match is returned as "subdialog". Otherwise,
 *  "subdialog" is set to NULL.
 */
int SCIPdialogFindEntry(
   DIALOG*          dialog,             /**< dialog */
   const char*      entryname,          /**< name of the dialog entry to find */
   DIALOG**         subdialog           /**< pointer to store the found dialog entry */
   )
{
   DIALOG** subdialogs;
   int nsubdialogs;
   int namelen;
   int nfound;
   int i;

   assert(dialog != NULL);
   assert(entryname != NULL);
   assert(subdialog != NULL);

   *subdialog = NULL;

   /* check entryname w.r.t. available dialog options */
   subdialogs = SCIPdialogGetSubdialogs(dialog);
   nsubdialogs = SCIPdialogGetNSubdialogs(dialog);
   namelen = strlen(entryname);
   nfound = 0;
   for( i = 0; i < nsubdialogs; ++i )
   {
      /* check, if the beginning of the sub dialog's name matches entryname */
      if( strncmp(entryname, SCIPdialogGetName(subdialogs[i]), namelen) == 0 )
      {
         *subdialog = subdialogs[i];
         nfound++;

         /* if entryname exactly matches the subdialog's name, use this subdialog */
         if( namelen == (int)strlen(SCIPdialogGetName(subdialogs[i])) )
            return 1;
      }
   }

   if( nfound != 1 )
      *subdialog = NULL;

   return nfound;
}

/** displays the dialog's menu */
RETCODE SCIPdialogDisplayMenu(
   DIALOG*          dialog,             /**< dialog */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   int len;
   int i;

   assert(dialog != NULL);

   /* display the dialog's sub menus */
   for( i = 0; i < dialog->nsubdialogs; ++i )
   {
      if( SCIPdialogIsSubmenu(dialog->subdialogs[i]) )
      {
         CHECK_OKAY( SCIPdialogDisplayMenuEntry(dialog->subdialogs[i], scip) );
      }
   }

   /* display the dialog's menu options */
   for( i = 0; i < dialog->nsubdialogs; ++i )
   {
      if( !SCIPdialogIsSubmenu(dialog->subdialogs[i]) )
      {
         CHECK_OKAY( SCIPdialogDisplayMenuEntry(dialog->subdialogs[i], scip) );
      }
   }

   if( dialog->nsubdialogs == 0 )
      printf("<no options available>\n");

   return SCIP_OKAY;
}

/** displays the entry for the dialog in it's parent's menu */
RETCODE SCIPdialogDisplayMenuEntry(
   DIALOG*          dialog,             /**< dialog */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   char name[MAXSTRLEN];

   assert(dialog != NULL);

   /* display the dialog's name */
   if( dialog->issubmenu )
      sprintf(name, "<%s>", dialog->name);
   else
      sprintf(name, "%s", dialog->name);
   printf("  %-20s  ", name);
   if( strlen(name) > 20 )
   {
      /* break the line, and start the description in the next line */
      printf("\n                   -->  ");
   }

   /* display the dialog's description */
   if( dialog->dialogdesc != NULL )
   {
      CHECK_OKAY( dialog->dialogdesc(scip, dialog) );
   }
   else
      printf(dialog->desc);
   printf("\n");

   return SCIP_OKAY;
}

/** displays all dialog entries with names starting with the given "entryname" */
RETCODE SCIPdialogDisplayCompletions(
   DIALOG*          dialog,             /**< dialog */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      entryname           /**< name of the dialog entry to find */
   )
{
   DIALOG** subdialogs;
   int nsubdialogs;
   int namelen;
   int i;

   assert(dialog != NULL);
   assert(entryname != NULL);

   /* check entryname w.r.t. available dialog options */
   subdialogs = SCIPdialogGetSubdialogs(dialog);
   nsubdialogs = SCIPdialogGetNSubdialogs(dialog);
   namelen = strlen(entryname);
   for( i = 0; i < nsubdialogs; ++i )
   {
      /* check, if the beginning of the sub dialog's name matches entryname */
      if( strncmp(entryname, SCIPdialogGetName(subdialogs[i]), namelen) == 0 )
      {
         CHECK_OKAY( SCIPdialogDisplayMenuEntry(subdialogs[i], scip) );
      }
   }

   return SCIP_OKAY;
}

/** gets the name of the current path in the dialog tree, separated by the given character */
void SCIPdialogGetPath(
   DIALOG*          dialog,             /**< dialog */
   const char       sepchar,            /**< separation character to insert in path */
   char*            path                /**< string buffer to store the path */
   )
{
   char s[MAXSTRLEN];

   assert(dialog != NULL);

   (void)strcpy(path, dialog->name);
   dialog = dialog->parent;
   while( dialog != NULL )
   {
      snprintf(s, MAXSTRLEN, "%s%c%s", dialog->name, sepchar, path);
      s[MAXSTRLEN-1] = '\0';
      (void)strcpy(path, s);
      dialog = dialog->parent;
   }
}

/** gets the command name of the dialog */
const char* SCIPdialogGetName(
   DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->name;
}

/** gets the description of the dialog */
const char* SCIPdialogGetDesc(
   DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->desc;
}

/** returns whether the dialog is a sub menu */
Bool SCIPdialogIsSubmenu(
   DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->issubmenu;
}

/** gets the parent dialog of the given dialog */
DIALOG* SCIPdialogGetParent(
   DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->parent;
}

/** gets the array of subdialogs associated with the given dialog */
DIALOG** SCIPdialogGetSubdialogs(
   DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->subdialogs;
}

/** gets the number of subdialogs associated with the given dialog */
int SCIPdialogGetNSubdialogs(
   DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->nsubdialogs;
}

/** gets the user defined data associated with the given dialog */
DIALOGDATA* SCIPdialogGetData(
   DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->dialogdata;
}
