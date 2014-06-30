#ident "@(#) $Id: grphsave.c,v 1.5 1999/09/30 16:02:10 bzfkocht Exp $"
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: grphsave.c                                                    */
/*   Name....: Graph File Saving Routines                                    */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "portab.h"
#include "grph.h"


void SCIPwriteStp(
   SCIP* scip,
   const GRAPH* g,
   FILE*        fp,
   SCIP_Real    offset
   )
{
   int i;

   assert(g  != NULL);
   assert(fp != NULL);

   fprintf(fp, "%8x STP File, STP Format Version %2d.%02d\n",
      STP_MAGIC, VERSION_MAJOR, VERSION_MINOR);

   fprintf(fp, "Section Comment\n");
   fprintf(fp, "Name \"%s\"\n", "noname");
   fprintf(fp, "End\n\n");

   if( !SCIPisEQ(scip, offset, 0.0) )
   {
      fprintf(fp, "Section Presolve\n");
      fprintf(fp, "Fixed %f \n", offset);
      fprintf(fp, "End\n\n");
   }

   fprintf(fp, "Section Graph\n");
   fprintf(fp, "Nodes %d\n", g->knots);
   fprintf(fp, "Edges %d\n", g->edges / 2);

   for(i = 0; i < g->edges; i += 2)
   {
      if (g->ieat[i] != EAT_FREE)
      {
         assert(g->oeat[i] != EAT_FREE);

         fprintf(fp, "E %d %d %g\n",
            g->tail[i] + 1,
            g->head[i] + 1,
            g->cost[i]);
      }
   }
   fprintf(fp, "End\n\n");
   fprintf(fp, "Section Terminals\n");
   fprintf(fp, "Terminals %d\n", g->terms);

   for(i = 0; i < g->knots; i++)
   {
      if (Is_term(g->term[i]))
         fprintf(fp, "T %d\n", i + 1);
   }
   fprintf(fp, "End\n\n");

   if (g->flags & GRAPH_HAS_COORDINATES)
   {
      fprintf(fp, "Section Coordinates\n");

      for(i = 0; i < g->knots; i++)
         fprintf(fp, "DD %d %d %d\n", i + 1, g->xpos[i], g->ypos[i]);

      fprintf(fp, "End\n\n");
   }

   fprintf(fp, "EOF\n");
}




/*---------------------------------------------------------------------------*/
/*--- Name     : STP Save Graph                                           ---*/
/*--- Function : Writes a graph to a file in STP format.                  ---*/
/*--- Parameter: Graph and Filepointer.                                   ---*/
/*--- Returns  : Nothing.                                                 ---*/
/*---------------------------------------------------------------------------*/
static void stp_save(
   const GRAPH* g,
   FILE*        fp
   )
{
   int   i;

   assert(g  != NULL);
   assert(fp != NULL);

   fprintf(fp, "%8x STP File, STP Format Version %2d.%02d\n",
      STP_MAGIC, VERSION_MAJOR, VERSION_MINOR);

   fprintf(fp, "Section Comment\n");
   fprintf(fp, "Name \"%s\"\n", "noname");
   fprintf(fp, "End\n\n");

   fprintf(fp, "Section Graph\n");
   fprintf(fp, "Nodes %d\n", g->knots);
   fprintf(fp, "Edges %d\n", g->edges / 2);

   for(i = 0; i < g->edges; i += 2)
   {
      if (g->ieat[i] != EAT_FREE)
      {
         assert(g->oeat[i] != EAT_FREE);

         fprintf(fp, "E %d %d %g\n",
            g->tail[i] + 1,
            g->head[i] + 1,
            g->cost[i]);
      }
   }
   fprintf(fp, "End\n\n");
   fprintf(fp, "Section Terminals\n");
   fprintf(fp, "Terminals %d\n", g->terms);

   for(i = 0; i < g->knots; i++)
   {
      if (Is_term(g->term[i]))
         fprintf(fp, "T %d\n", i + 1);
   }
   fprintf(fp, "End\n\n");

   if (g->flags & GRAPH_HAS_COORDINATES)
   {
      fprintf(fp, "Section Coordinates\n");

      for(i = 0; i < g->knots; i++)
         fprintf(fp, "DD %d %d %d\n", i + 1, g->xpos[i], g->ypos[i]);

      fprintf(fp, "End\n\n");
   }

   fprintf(fp, "EOF\n");
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Beasley Save Graph                                       ---*/
/*--- Function : Writes a graph to a file in Beasleys format.             ---*/
/*--- Parameter: Graph and Filepointer.                                   ---*/
/*--- Returns  : Nothing.                                                 ---*/
/*---------------------------------------------------------------------------*/
static void bea_save(
   const GRAPH* g,
   FILE*        fp)
{
   int i;

   assert(g  != NULL);
   assert(fp != NULL);

   fprintf(fp, "%d %d\n",
      g->knots,
      g->edges / 2);

   for(i = 0; i < g->edges; i += 2)
   {
      if (g->ieat[i] != EAT_FREE)
      {
         assert(g->oeat[i] != EAT_FREE);

         fprintf(fp, "%d %d %g\n",
            g->tail[i] + 1,
            g->head[i] + 1,
            g->cost[i]);
      }
   }
   fprintf(fp, "%d\n", g->terms);

   for(i = 0; i < g->knots; i++)
   {
      if (Is_term(g->term[i]))
         fprintf(fp, "%d\n", i + 1);
   }
}

/*---------------------------------------------------------------------------*/
/*--- Name     : Graph Save                                               ---*/
/*--- Function : Write a graph to a file.                                 ---*/
/*--- Parameter: Graph, Filename, Filetype (Fileformat).                  ---*/
/*--- Returns  : Nothing.                                                 ---*/
/*---------------------------------------------------------------------------*/
void graph_save(
   const GRAPH* g,
   const char*  filename,
   FILETYPE     type)
{
   const char* msg_writing_s = "Writing Graph to File %s:\n";

   FILE* fp;

   assert(g                != NULL);
   assert(graph_valid(g));
   assert(filename         != NULL);
   assert(strlen(filename) > 0);
   assert((type == FF_BEA) || (type == FF_STP));

   if ((fp = fopen(filename, "w")) == NULL)
      perror(filename);
   else
   {
      printf(msg_writing_s, filename);

      switch(type)
      {
      case FF_BEA :
         bea_save(g, fp);
         break;
      case FF_STP :
         stp_save(g, fp);
         break;
      case FF_PRB :
      case FF_GRD :
      default     :
         /* CONSTCOND */
         assert(FALSE);
      }
      fclose(fp);
   }
}
