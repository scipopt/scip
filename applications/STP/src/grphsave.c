/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   grphsave.c
 * @brief  Saving routines for Steiner problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 * This file includes several saving routines for Steiner problems
 *
 * A list of all interface methods can be found in grph.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

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
   int e;
   int root = g->source;
   int hopfactor;

   assert(g  != NULL);
   assert(fp != NULL);

   fprintf(fp, "%8x STP File, STP Format Version %2d.%02d\n",
      STP_MAGIC, VERSION_MAJOR, VERSION_MINOR);

   fprintf(fp, "Section Comment\n");
   fprintf(fp, "Name ");
   switch( g->stp_type )
   {
      case STP_SPG:
         fprintf(fp, "\"%s\"\n", "SPG");
         break;

      case STP_PCSPG:
         fprintf(fp, "\"%s\"\n", "UNKNOWN");
         break;

      case STP_RPCSPG:
         fprintf(fp, "\"%s\"\n", "RPCST");
         break;

      case STP_NWSPG:
         fprintf(fp, "\"%s\"\n", "NWSPG");
         break;

      case STP_DCSTP:
         fprintf(fp, "\"%s\"\n", "UNKNOWN");
         break;

      case STP_RSMT:
         fprintf(fp, "\"%s\"\n", "RSMT");
         break;

      case STP_OARSMT:
         fprintf(fp, "\"%s\"\n", "OARSMT");
         break;

      case STP_MWCSP:
         fprintf(fp, "\"%s\"\n", "UNKNOWN");
         break;

      case STP_DHCSTP:
         fprintf(fp, "\"%s\"\n", "UNKNOWN");
         break;

      default:
         fprintf(fp, "\"%s\"\n", "UNKNOWN");
   }
   fprintf(fp, "Remark \"Transformed\"\n");
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

   for( i = 0; i < g->edges; i += 2 )
   {
      if (g->ieat[i] != EAT_FREE)
      {
         assert(g->oeat[i] != EAT_FREE);

         if( g->stp_type == STP_SPG || g->stp_type == STP_DCSTP || g->stp_type == STP_RSMT
           || g->stp_type == STP_OARSMT )
            fprintf(fp, "E ");
         else
            fprintf(fp, "AA ");

         fprintf(fp, "%d %d ", g->tail[i] + 1, g->head[i] + 1);

         if( g->stp_type == STP_SPG || g->stp_type == STP_DCSTP || g->stp_type == STP_RSMT
           || g->stp_type == STP_OARSMT )
            fprintf(fp, "%f\n", g->cost[i]);
         else
            fprintf(fp, "%f %f\n", g->cost[i], g->cost[Edge_anti(i)]);
      }
   }
   fprintf(fp, "End\n\n");
   fprintf(fp, "Section Terminals\n");
   fprintf(fp, "Terminals %d\n", g->terms);

   if( g->stp_type == STP_RPCSPG )
      fprintf(fp, "Root %d\n", g->source + 1);

   for(i = 0; i < g->knots; i++)
   {
      if (Is_term(g->term[i]) && (g->stp_type != STP_RPCSPG || i != g->source))
         fprintf(fp, "T %d\n", i + 1);
   }
   fprintf(fp, "End\n\n");
/*
   if (g->flags & GRAPH_HAS_COORDINATES)
   {
      fprintf(fp, "Section Coordinates\n");

      for(i = 0; i < g->knots; i++)
         fprintf(fp, "DD %d %d %d\n", i + 1, g->xpos[i], g->ypos[i]);

      fprintf(fp, "End\n\n");
   }
*/
   /* Hop-Constrained STP */
   if( g->stp_type == STP_DHCSTP )
   {
      fprintf(fp, "Section Hop Constraint\n");
      fprintf(fp, "limit %d\n", g->hoplimit);
      for( e = 0; e < g->edges; e++ )
      {
         /* TODO: When presolving is used: MODIFY */
         hopfactor = 1;
         fprintf(fp, "HC %d %d\n", e + 1, hopfactor);
      }
      fprintf(fp, "End\n\n");
   }

   /* Degree-Constrained STP */
   /* It is just enough to know that the degree constraint is required */
   if( g->stp_type == STP_DCSTP )
   {
      fprintf(fp, "Section Degree Constraint\n");
      fprintf(fp, "End\n\n");
   }

	 /* PRIZECOLLECTING STP */
   if( g->stp_type == STP_PCSPG || g->stp_type == STP_MWCSP )
   {
      fprintf(fp, "Section Prize Collecting Constraint\n");

      for( e = g->outbeg[root]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( !Is_term(g->term[g->head[e]]) )
            fprintf(fp, "PC %d\n", e + 1);
      }

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
/*
   if (g->flags & GRAPH_HAS_COORDINATES)
   {
      fprintf(fp, "Section Coordinates\n");

      for(i = 0; i < g->knots; i++)
         fprintf(fp, "DD %d %d %d\n", i + 1, g->xpos[i], g->ypos[i]);

      fprintf(fp, "End\n\n");
   }
*/
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
