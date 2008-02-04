/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_ppm.c,v 1.1 2008/02/04 12:32:17 bzfpfets Exp $"

/**@file   reader_ppm.c
 * @brief  PPM file reader
 * @author Tobias Achterberg
 * @author Michael Winkler
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h>
#endif
#include <ctype.h>

#include "scip/reader_ppm.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"

#define READER_NAME             "ppmreader"
#define READER_DESC             "file writer for ppm file format"
#define READER_EXTENSION        "ppm"


/*
 * Data structures
 */
#define PPM_MAX_LINELEN       70


/*
 * Local methods (for writing)
 */


/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< vars array to get active variables for */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n inrc/scip/reader_ppm.c linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( scalars != NULL );
   assert( nvars != NULL );
   assert( constant != NULL );

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize) );

      if( requiredsize > *nvars )
      {
         *nvars = requiredsize;
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, *nvars ) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &scalars, *nvars ) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize) );
         assert( requiredsize <= *nvars );
      }
   }
   else
   {
      for( v = 0; v < *nvars; ++v )
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );
   }
   return SCIP_OKAY;
}

/** clears the given line buffer */
static
void clearLine(
   char*                 linebuffer,         /**< line */
   int*                  linecnt             /**< number of charaters in line */
   )
{
   assert( linebuffer != NULL );
   assert( linecnt != NULL );

   (*linecnt) = 0;
   linebuffer[0] = '\0';
}

/** ends the given line with '\\0' and prints it to the given file stream */
static
void endLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt             /**< number of charaters in line */
   )
{
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );

   if( (*linecnt) > 0 )
   {
      if (linebuffer[(*linecnt)-1] == '\n')
         linebuffer[(*linecnt)-1] = '\0';
      else
         linebuffer[(*linecnt)] = '\0';
      SCIPinfoMessage(scip, file, "%s\n", linebuffer);
      clearLine(linebuffer, linecnt);
   }
}

/** appends extension to line and prints it to the give file stream if the line exceeded
    PPM_PRINTLEN */
static
void appendLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt,            /**< number of charaters in line */
   const char*           extension           /**< string to extent the line */
   )
{
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );
   assert( extension != NULL );

   sprintf(linebuffer, "%s%s", linebuffer, extension);
   (*linecnt) += strlen(extension) + 1;

   if( (*linecnt) > PPM_MAX_LINELEN )
      endLine(scip, file, linebuffer, linecnt);
}


/* calculates the color value for a given coefficient */
static
unsigned short calcColorValue(
   SCIP_Real             coef,               /**< Coefficient to scale */
   SCIP_Real             scale,              /**< scaling coefficient */
   unsigned short        maxvalue            /**< maximum value for output */
   )
{
   return (maxvalue - (unsigned short) (coef/scale*maxvalue));
}


/* print row in PPM format to file stream */
static
void printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VAR**            consvars,           /**< array of constraint variables */
   SCIP_Real*            consvals,           /**< array of constraint values */
   int                   nconsvars,          /**< number of constraint variables */
   int                   nvars,              /**< number of variables */
   SCIP_Real             maxcoef             /**< maximal coefficient */
   )
{
   assert (nconsvars > 0);
   int v, i = 0, j;
   char linebuffer[PPM_MAX_LINELEN + 1];
   int linecnt;
   int varindex = -1, actvarindex, maxvarindex = 0;
   int indexconsvar = 0;

   char buffer[PPM_MAX_LINELEN];

   assert( scip != NULL );

   clearLine(linebuffer, &linecnt);

   for( v = 0; v < nconsvars; ++v )
   {
      if (maxvarindex < SCIPvarGetProbindex(consvars[v]))
         maxvarindex = SCIPvarGetProbindex(consvars[v]);
   }

   assert(maxvarindex < nvars);

   /* print coefficients */
   for (v = 0; v < nconsvars; ++v)
   {
      actvarindex = maxvarindex;
      for (j = 0; j < nconsvars; ++j)
      {
         if ( SCIPvarGetProbindex(consvars[j]) <= actvarindex && SCIPvarGetProbindex(consvars[j]) > varindex )
         {
            actvarindex = SCIPvarGetProbindex(consvars[j]);
            indexconsvar = j;
         }
      }
      varindex = actvarindex;

      for( ; i < varindex; ++i )
      {
         appendLine(scip, file, linebuffer, &linecnt, "\t255\t255\t255\t");
         if ( (i+1)%5 == 0)
            endLine(scip, file, linebuffer, &linecnt);
      }

      sprintf(buffer, "\t255\t%d\t%d\t", calcColorValue(consvals[indexconsvar], maxcoef, 255), calcColorValue(consvals[indexconsvar], maxcoef, 255));

      appendLine(scip, file, linebuffer, &linecnt, buffer);
      i = varindex+1;
      if ( (i)%5 == 0)
         endLine(scip, file, linebuffer, &linecnt);
   }

   for( ; i < nvars; ++i )
   {
      appendLine(scip, file, linebuffer, &linecnt, "\t255\t255\t255\t");
      if ( (i+1)%5 == 0)
         endLine(scip, file, linebuffer, &linecnt);
   }

   endLine(scip, file, linebuffer, &linecnt);
}



/** prints given linear constraint information in PPM format to file stream */
static
SCIP_RETCODE printLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   int                   ncompletevars,      /**< number of variables in whole problem */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Real*            maxcoef,            /**< maximal coefficient */
   SCIP_Bool             printbool           /**< print row or calculate maximum coefficient */
   )
{
   int v;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   //assert( vars != NULL );
   assert( nvars > 0 );

   /* duplicate variable and value array */
   nactivevars = nvars;
   SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars );
   if( vals != NULL )
      SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars );
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

      for( v = 0; v < nactivevars; ++v )
         activevals[v] = 1.0;
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );


   if (!printbool)
      for (v = 0; v < nactivevars; ++v)
      {
         if ( REALABS(activevals[v]) > *maxcoef)
            *maxcoef = REALABS(activevals[v]);
      }
   else
   {
      assert (*maxcoef > 0);
      /* print constraint */
      printRow(scip, file, activevars, activevals, nactivevars, ncompletevars, *maxcoef);
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreePpm NULL

/** problem reading method of reader */
#define readerReadPpm NULL

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWritePpm)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPwritePpm(scip, file, name, transformed, vars, nvars, conss, nconss, result) );

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the ppm file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderPpm(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create ppm reader data */
   readerdata = NULL;

   /* include ppm reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerFreePpm, readerReadPpm, readerWritePpm, readerdata) );

   return SCIP_OKAY;
}



/* writes problem to file */
SCIP_RETCODE SCIPwritePpm(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**< output file, or NULL if standard output should be used */
   const char*        name,               /**< problem name */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_VAR**         vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                nvars,              /**< number of mutable variables in the problem */
   SCIP_CONS**        conss,              /**< array with constraints of the problem */
   int                nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*       result              /**< pointer to store the result of the file writing call */
   )
{
   int c,v,i;

   int linecnt;
   char linebuffer[PPM_MAX_LINELEN + 1];

   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;

   SCIP_Real maxcoef = 0;
   SCIP_Bool printbool = FALSE;

   assert( scip != NULL );

   /* print statistics as comment to file */
   SCIPinfoMessage(scip, file, "P3\n");
   SCIPinfoMessage(scip, file, "# %s\n", name);
   SCIPinfoMessage(scip, file, "%d %d\n", nvars, nconss);
   SCIPinfoMessage(scip, file, "255\n");

   clearLine(linebuffer, &linecnt);

   for (i = 0; i < 2; ++i)
   {
      if (i)
      {
         printbool = TRUE;
         SCIPdebugPrintf("Maximal coefficient = %g\n", maxcoef);
      }

      for (c = 0; c < nconss; ++c)
      {
         cons = conss[c];
         assert( cons != NULL);

         /* in case the transformed is written only constraint are posted which are enabled in the current node */
         if( transformed && !SCIPconsIsEnabled(cons) )
            continue;

         conshdlr = SCIPconsGetHdlr(cons);
         assert( conshdlr != NULL );

         conshdlrname = SCIPconshdlrGetName(conshdlr);
         assert( transformed == SCIPconsIsTransformed(cons) );

         if( strcmp(conshdlrname, "linear") == 0 )
         {
            SCIP_CALL( printLinearCons(scip, file, SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
                  SCIPgetNVarsLinear(scip, cons), nvars, transformed, &maxcoef, printbool) );
         }
         else if( strcmp(conshdlrname, "setppc") == 0 )
         {
            consvars = SCIPgetVarsSetppc(scip, cons);
            nconsvars = SCIPgetNVarsSetppc(scip, cons);

            SCIP_CALL( printLinearCons(scip, file, consvars, NULL, nconsvars, nvars, transformed, &maxcoef, printbool) );
         }
         else if ( strcmp(conshdlrname, "logicor") == 0 )
         {
            SCIP_CALL( printLinearCons(scip, file, SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons),
                  nvars, transformed, &maxcoef, printbool) );
         }
         else if ( strcmp(conshdlrname, "knapsack") == 0 )
         {
            SCIP_Longint* weights;

            consvars = SCIPgetVarsKnapsack(scip, cons);
            nconsvars = SCIPgetNVarsKnapsack(scip, cons);

            /* copy Longint array to SCIP_Real array */
            weights = SCIPgetWeightsKnapsack(scip, cons);
            SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
            for( v = 0; v < nconsvars; ++v )
               consvals[v] = weights[v];

            SCIP_CALL( printLinearCons(scip, file, consvars, consvals, nconsvars, nvars, transformed, &maxcoef, printbool) );

            SCIPfreeBufferArray(scip, &consvals);
         }
         else if ( strcmp(conshdlrname, "varbound") == 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
            SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

            consvars[0] = SCIPgetVarVarbound(scip, cons);
            consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

            consvals[0] = 1.0;
            consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

            SCIP_CALL( printLinearCons(scip, file, consvars, consvals, 2, nvars, transformed, &maxcoef, printbool) );

            SCIPfreeBufferArray(scip, &consvars);
            SCIPfreeBufferArray(scip, &consvals);
         }
         else
         {
            SCIPwarningMessage("constraint handler <%s> can not print requested format\n", conshdlrname );
            SCIPinfoMessage(scip, file, "\\ ");
            SCIP_CALL( SCIPprintCons(scip, cons, file) );
         }
      }
   }
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
