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

/**@file   disp.c
 * @brief  datastructures and methods for displaying runtime statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "disp.h"



/** display column */
struct Disp
{
   char*            name;               /**< name of display column */
   char*            desc;               /**< description of display column */
   char*            header;             /**< head line of display column */
   DECL_DISPFREE    ((*dispfree));      /**< destructor of display column */
   DECL_DISPINIT    ((*dispinit));      /**< initialise display column */
   DECL_DISPEXIT    ((*dispexit));      /**< deinitialise display column */
   DECL_DISPOUTPUT  ((*dispoutput));    /**< output method */
   DISPDATA*        dispdata;           /**< display column data */
   int              width;              /**< width of display column (no. of chars used) */
   int              priority;           /**< priority of display column */
   int              position;           /**< relative position of display column */
   unsigned int     stripline:1;        /**< should the column be separated with a line from its right neighbour? */
   unsigned int     initialized:1;      /**< is display column initialized? */
   unsigned int     active:1;           /**< should column be displayed to the screen? */
};



/* display column methods */

/** creates a display column */
RETCODE SCIPdispCreate(
   DISP**           disp,               /**< pointer to store display column */
   const char*      name,               /**< name of display column */
   const char*      desc,               /**< description of display column */
   const char*      header,             /**< head line of display column */
   DECL_DISPFREE    ((*dispfree)),      /**< destructor of display column */
   DECL_DISPINIT    ((*dispinit)),      /**< initialise display column */
   DECL_DISPEXIT    ((*dispexit)),      /**< deinitialise display column */
   DECL_DISPOUTPUT  ((*dispoutput)),    /**< output method */
   DISPDATA*        dispdata,           /**< display column data */
   int              width,              /**< width of display column (no. of chars used) */
   int              priority,           /**< priority of display column */
   int              position,           /**< relative position of display column */
   Bool             stripline           /**< should the column be separated with a line from its right neighbour? */
   )
{
   assert(disp != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(header != NULL);
   assert(dispoutput != NULL);
   assert(width >= 0);

   ALLOC_OKAY( allocMemory(disp) );
   ALLOC_OKAY( duplicateMemoryArray(&(*disp)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*disp)->desc, desc, strlen(desc)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*disp)->header, header, strlen(header)+1) );
   (*disp)->dispfree = dispfree;
   (*disp)->dispinit = dispinit;
   (*disp)->dispexit = dispexit;
   (*disp)->dispoutput = dispoutput;
   (*disp)->dispdata = dispdata;
   (*disp)->width = width;
   (*disp)->priority = priority;
   (*disp)->position = position;
   (*disp)->stripline = stripline;
   (*disp)->initialized = FALSE;
   (*disp)->active = FALSE;

   return SCIP_OKAY;
}
   
/** frees memory of display column */
RETCODE SCIPdispFree(
   DISP**           disp,               /**< pointer to display column data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(disp != NULL);
   assert(*disp != NULL);
   assert(!(*disp)->initialized);

   /* call destructor of display column */
   if( (*disp)->dispfree != NULL )
   {
      CHECK_OKAY( (*disp)->dispfree(scip, *disp) );
   }

   freeMemoryArray(&(*disp)->name);
   freeMemoryArray(&(*disp)->desc);
   freeMemoryArray(&(*disp)->header);
   freeMemory(disp);

   return SCIP_OKAY;
}

/** initializes display column */
RETCODE SCIPdispInit(
   DISP*            disp,               /**< display column */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(disp != NULL);
   assert(scip != NULL);

   if( disp->initialized )
   {
      char s[MAXSTRLEN];
      sprintf(s, "display column <%s> already initialized", disp->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( disp->dispinit != NULL )
   {
      CHECK_OKAY( disp->dispinit(scip, disp) );
   }
   disp->initialized = TRUE;

   return SCIP_OKAY;
}

/** deinitializes display column */
RETCODE SCIPdispExit(
   DISP*            disp,               /**< display column */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(disp != NULL);
   assert(scip != NULL);

   if( !disp->initialized )
   {
      char s[MAXSTRLEN];
      sprintf(s, "display column <%s> not initialized", disp->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( disp->dispexit != NULL )
   {
      CHECK_OKAY( disp->dispexit(scip, disp) );
   }
   disp->initialized = FALSE;

   return SCIP_OKAY;
}

/** output display column to screen */
RETCODE SCIPdispOutput(
   DISP*            disp,               /**< display column */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(disp != NULL);
   assert(disp->dispoutput != NULL);
   assert(scip != NULL);

   CHECK_OKAY( disp->dispoutput(scip, disp, stdout) );

   return SCIP_OKAY;
}

/** gets name of display column */
const char* SCIPdispGetName(
   DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->name;
}

/** gets user data of display column */
DISPDATA* SCIPdispGetData(
   DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->dispdata;
}

/** sets user data of display column; user has to free old data in advance! */
void SCIPdispSetData(
   DISP*            disp,               /**< display column */
   DISPDATA*        dispdata            /**< new display column user data */
   )
{
   assert(disp != NULL);

   disp->dispdata = dispdata;
}

/** gets position of display column */
int SCIPdispGetPosition(
   DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->position;
}

/** is display column initialized? */
Bool SCIPdispIsInitialized(
   DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->initialized;
}

/** prints one line of output with the active display columns */
RETCODE SCIPdispPrintLine(
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   Bool             forcedisplay        /**< should the line be printed without regarding frequency? */
   )
{
   assert(set != NULL);
   assert(set->dispfreq >= -1);
   assert(set->dispheaderfreq >= -1);
   assert(stat != NULL);

   if( set->verblevel < SCIP_VERBLEVEL_NORMAL || set->dispfreq == -1 )
      return SCIP_OKAY;

   if( forcedisplay
      || (stat->nnodes != stat->lastdispnode
         && set->dispfreq > 0
         && (stat->nnodes % set->dispfreq == 0 || stat->nnodes == 1)) )
   {
      int i;
      int j;
      Bool stripline;

      /* display header line */
      if( (set->dispheaderfreq == 0 && stat->ndisplines == 0)
         || (set->dispheaderfreq > 0 && stat->ndisplines % set->dispheaderfreq == 0) )
      {
         int fillspace;

         stripline = FALSE;
         for( i = 0; i < set->ndisps; ++i )
         {
            assert(set->disps[i] != NULL);
            if( set->disps[i]->active )
            {
               if( stripline )
                  printf("|");
               fillspace = set->disps[i]->width - strlen(set->disps[i]->header);
               for( j = 0; j < (fillspace)/2; ++j )
                  printf(" ");
               printf(set->disps[i]->header);
               for( j = 0; j < (fillspace+1)/2; ++j )
                  printf(" ");
               stripline = set->disps[i]->stripline;
            }
         }
         printf("\n");
      }

      /* display node information line */
      stripline = FALSE;
      for( i = 0; i < set->ndisps; ++i )
      {
         assert(set->disps[i] != NULL);
         if( set->disps[i]->active )
         {
            if( stripline )
               printf("|");
            CHECK_OKAY( SCIPdispOutput(set->disps[i], set->scip) );
            stripline = set->disps[i]->stripline;
         }
      }
      printf("\n");

      stat->lastdispnode = stat->nnodes;
      stat->ndisplines++;
   }

   return SCIP_OKAY;
}

static
DECL_SORTPTRCOMP(dispCmp)
{
   return ((DISP*)elem2)->priority - ((DISP*)elem1)->priority;
}

/** activates all display lines fitting in the display w.r. to priority */
RETCODE SCIPdispAutoActivate(
   const SET*       set                 /**< global SCIP settings */
   )
{
   DISP** disps;
   int width;
   int actwidth;
   int i;

   assert(set != NULL);

   /* sort display columns w.r. to their priority */
   ALLOC_OKAY( duplicateMemoryArray(&disps, set->disps, set->ndisps) );
   SCIPbsortPtr((void**)disps, set->ndisps, dispCmp);

   /* beginning with highest priority display column, activate columns as long as it fits into display width */
   width = 0;
   for( i = 0; i < set->ndisps; ++i )
   {
      actwidth = disps[i]->width;
      if( disps[i]->stripline )
         actwidth++;
      if( width + actwidth <= set->dispwidth )
      {
         disps[i]->active = TRUE;
         width += actwidth;
      }
      else
         disps[i]->active = FALSE;
   }

   /* free temporary memory */
   freeMemoryArray(&disps);

   return SCIP_OKAY;
}

static
const char powerchar[] = {' ', 'k', 'M', 'G', 'T', 'P', 'E'};
#define MAXPOWER 6

/** displays an integer in decimal form fitting in a given width */
void SCIPdispDecimal(
   FILE*            file,               /**< output stream */
   Longint          val,                /**< value to display */
   int              width               /**< width to fit into */
   )
{
   assert(width >= 1);

   if( width == 1 )
   {
      if( val < 0 )
         fprintf(file, "-");
      else if( val < 10 )
         fprintf(file, "%lld", val);
      else
         fprintf(file, "+");
   }
   else
   {
      char format[MAXSTRLEN];
      Longint maxval;
      int power;
      int i;

      maxval = 1;
      for( i = 0; i < width-1; ++i )
         maxval *= 10;
      if( val < 0 )
         maxval /= 10;
      power = 0;
      while( ABS(val) >= maxval && power < MAXPOWER )
      {
         power++;
         val /= 1000;
      }
      sprintf(format, "%%%dlld%c", width-1, powerchar[power]);

      if( width == 2 && val < 0 )
         fprintf(file, "-%c", powerchar[power]);
      else
         fprintf(file, format, val);
   }
}
