/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_gms.c,v 1.3 2009/07/21 06:04:15 bzfgleix Exp $"

/**@file   reader_gms.c
 * @ingroup FILEReaders 
 * @brief  GAMS file reader
 * @author Ambros Gleixner
 *
 * @todo Test for uniqueness of variable names (after cutting down).
 * @todo Check for words reserved for GAMS.
 * @todo Routines for reading.
 * @todo Can general SOS constraints be modelled in GAMS?
 * @todo Can indicator constraints be modelled within one file?
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

#include "scip/reader_gms.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/pub_misc.h"

#define READER_NAME             "gmsreader"
#define READER_DESC             "file reader for MI(NL)Ps in GAMS file format"
#define READER_EXTENSION        "gms"


/*
 * Data structures
 */
#define GMS_MAX_LINELEN      256
#define GMS_MAX_PUSHEDTOKENS 2
#define GMS_INIT_COEFSSIZE   8192
#define GMS_MAX_PRINTLEN     256       /**< the maximum length of any line is 255 + '\\0' = 256*/
#define GMS_MAX_NAMELEN      64        /**< the maximum length for any name is 63 + '\\0' = 64 */
#define GMS_PRINTLEN         100



/*
 * Local methods (for writing)
 */

#if 0
/* prints variable name LP format conform; always use this method to stay consistent
 *
 * 1) variable names should not start with a digit
 * 2) avoid variable name starting with an 'e' or 'E' since this notation is reserved for exponential entries
 */
static
void printVarName(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file,               /**< output file (or NULL for standard output) */
   SCIP_VAR*            var,                /**< variable */
   SCIP_Bool            genericnames        /**< use generic variable names? */
   )
{
   const char* name;

   assert( scip != NULL );
   assert( var != NULL );

   name = SCIPvarGetName(var);
   assert( name != NULL );

   if( genericnames || name[0] == '\0' )
      SCIPinfoMessage(scip, file, "x%d", SCIPvarGetProbindex(var) + 1);
   else
   {
      if( isdigit(name[0]) || name[0] == 'e' || name[0] == 'E' )
         SCIPinfoMessage(scip, file, "_%s", name);
      else
         SCIPinfoMessage(scip, file, "%s", name);
   }
}
#endif

/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< vars array to get active variables for */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
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
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &scalars, requiredsize) );
         
         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
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
      linebuffer[(*linecnt)] = '\0';
      SCIPinfoMessage(scip, file, "%s\n", linebuffer);
      clearLine(linebuffer, linecnt);
   }
}

/** appends extension to line and prints it to the give file stream if the
 *  line exceeded the length given in the define GMS_PRINTLEN */
static
void appendLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt,            /**< number of characters in line */
   const char*           extension           /**< string to extend the line */
   )
{
   int len;
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );
   assert( extension != NULL );
   assert( strlen(linebuffer) + strlen(extension) < GMS_MAX_PRINTLEN );

   /* NOTE: avoid
    *   sprintf(linebuffer, "%s%s", linebuffer, extension); 
    * because of overlapping memory areas in memcpy used in sprintf.
    */
   len = strlen(linebuffer);
   strncat(linebuffer, extension, GMS_MAX_PRINTLEN - len);

   (*linecnt) += (int) strlen(extension);

   SCIPdebugMessage("linebuffer <%s>, length = %zd\n", linebuffer, strlen(linebuffer));
   
   if( (*linecnt) > GMS_PRINTLEN )
      endLine(scip, file, linebuffer, linecnt);
}



/* print row in GAMS format to file stream */
static
void printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=e=", "=l=", or "=g=") */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of values */
   int                   nvars,              /**< number of variables */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   int v;
   char linebuffer[GMS_MAX_PRINTLEN] = { '\0' };
   int linecnt;

   SCIP_VAR* var;
   char varname[GMS_MAX_NAMELEN];
   char consname[GMS_MAX_NAMELEN + 4]; /* four extra characters for ' .. ' */
   char buffer[GMS_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strcmp(type, "=e=") == 0 || strcmp(type, "=l=") || strcmp(type, "=g=") );
   assert( nvars == 0 || (vars != NULL && vals != NULL) );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   if ( strlen(rowname) > 0 || strlen(rownameextension) > 0 )
   {
      (void) SCIPsnprintf(consname, GMS_MAX_NAMELEN + 4, "%s%s .. ", rowname, rownameextension);
      appendLine(scip, file, linebuffer, &linecnt, consname);
   }

   /* print coefficients */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );

      /* we start a new line; therefore we tab this line */
      if (linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, "     ");

      (void) SCIPsnprintf(varname, GMS_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %+.15g*%s", vals[v], varname);
      
      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   /* print right hand side */
   if( SCIPisZero(scip, rhs) )
      rhs = 0.0;

   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s %+.15g;", type, rhs);

   /* we start a new line; therefore we tab this line */
   if (linecnt == 0 )
      appendLine(scip, file, linebuffer, &linecnt, "     ");
   appendLine(scip, file, linebuffer, &linecnt, buffer);

   endLine(scip, file, linebuffer, &linecnt);
}



/** prints given linear constraint information in GAMS format to file stream */
static
SCIP_RETCODE printLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< name of the row */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int v;
   SCIP_VAR** activevars = NULL;
   SCIP_Real* activevals = NULL;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( rowname != NULL );

   /* The GAMS format does not forbid that the variable array is empty */
   assert( nvars == 0 || vars != NULL );
   assert( nvars > 0 || vars == NULL );

   assert( lhs <= rhs );
   
   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   nactivevars = nvars;
   if( nvars > 0 ) 
   {
      /* duplicate variable and value array */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars) );
      if( vals != NULL )
         SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars) );
      else
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );
         
         for( v = 0; v < nactivevars; ++v )
            activevals[v] = 1.0;
      }

      /* retransform given variables to active variables */
      SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );
   }
   
   /* print row(s) in GAMS format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* print equality constraint */
      printRow(scip, file, rowname, "", "=e=", activevars, activevals, nactivevars, rhs - activeconstant);
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         printRow(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", "=g=",
            activevars, activevals, nactivevars, lhs - activeconstant);
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         printRow(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "=l=",
            activevars, activevals, nactivevars, rhs - activeconstant);
      }
   }

   if( nvars > 0 )
   {
      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &activevars);
      SCIPfreeBufferArray(scip, &activevals);
   }
   
   return SCIP_OKAY;
}



/** method check if the variable names are not longer than GMS_MAX_NAMELEN */
static
void checkVarnames(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_VAR**         vars,               /**< array of variables */
   int                nvars               /**< number of variables */
   )
{
   int v;
   SCIP_VAR* var;

   assert( scip != NULL );
   assert( vars != NULL );

   /* check if the variable names are not to long */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );

      if( strlen(SCIPvarGetName(var)) > GMS_MAX_NAMELEN )
      {
         SCIPwarningMessage("there is a variable name which has to be cut down to %d characters; GAMS model might be corrupted\n", 
            GMS_MAX_NAMELEN - 1);
         return;
      }
   }
}

/** method check if the constraint names are not longer than GMS_MAX_NAMELEN */
static
void checkConsnames(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONS**        conss,              /**< array of constraints */
   int                nconss,             /**< number of constraints */
   SCIP_Bool          transformed         /**< TRUE iff problem is the transformed problem */
   )
{
   int c;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   assert( scip != NULL );
   assert( conss != NULL );

   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL );

      /* in case the transformed is written, only constraints are posted which are enabled in the current node */
      if( transformed && !SCIPconsIsEnabled(cons) )
         continue;

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 ) /* ambros: or quadratic */
      {
         SCIP_Real lhs = SCIPgetLhsLinear(scip, cons);
         SCIP_Real rhs = SCIPgetLhsLinear(scip, cons);

         if( (SCIPisEQ(scip, lhs, rhs) && strlen(SCIPconsGetName(conss[c])) > GMS_MAX_NAMELEN)
            || ( !SCIPisEQ(scip, lhs, rhs) && strlen(SCIPconsGetName(conss[c])) > GMS_MAX_NAMELEN - 4) )
         {
            SCIPwarningMessage("there is a constraint name which has to be cut down to %d characters;\n",
               GMS_MAX_NAMELEN - 1);
            return;
         }
      }
      else if( strlen(SCIPconsGetName(conss[c])) > GMS_MAX_NAMELEN )
      {
         SCIPwarningMessage("there is a constraint name which has to be cut down to %d characters;\n",
            GMS_MAX_NAMELEN - 1);
         return;
      }
   }
}



/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeGms NULL

/** problem reading method of reader */
#define readerReadGms NULL

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteGms)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPwriteGms(scip, file, name, transformed, objsense, objscale, objoffset, vars,
			  nvars, nbinvars, nintvars, nimplvars, ncontvars, conss, nconss, result) );

   return SCIP_OKAY;
}



/*
 * reader specific interface methods
 */

/** includes the gms file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderGms(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create gms reader data */
   readerdata = NULL;

   /* include gms reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerFreeGms, readerReadGms, readerWriteGms, readerdata) );

   /* add gms reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/gmsreader/freeints", "are integer variables free by default (depending on GAMS version)?",
         NULL, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}



/* writes problem to gms file */
SCIP_RETCODE SCIPwriteGms(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**< output file, or NULL if standard output should be used */
   const char*        name,               /**< problem name */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_OBJSENSE      objsense,           /**< objective sense */
   SCIP_Real          objscale,           /**< scalar applied to objective function; external objective value is
   					       extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real          objoffset,          /**< objective offset from bound shifting and fixing */
   SCIP_VAR**         vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                nvars,              /**< number of mutable variables in the problem */
   int                nbinvars,           /**< number of binary variables */
   int                nintvars,           /**< number of general integer variables */
   int                nimplvars,          /**< number of implicit integer variables */
   int                ncontvars,          /**< number of continuous variables */
   SCIP_CONS**        conss,              /**< array with constraints of the problem */
   int                nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*       result              /**< pointer to store the result of the file writing call */
   )
{
   int c,v;

   int linecnt;
   char linebuffer[GMS_MAX_PRINTLEN];

   char varname[GMS_MAX_NAMELEN];
   char buffer[GMS_MAX_PRINTLEN];

   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;

   char consname[GMS_MAX_NAMELEN];

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;

   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool freeints;
   SCIP_Bool nondefbounds;

   assert( scip != NULL );
   assert( nvars > 0 );

#if 0
   /* find indicator constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "indicator");
   assert( conshdlr != NULL );

   /* create hashtable storing linear constraints that should not be output */
   SCIP_CALL( SCIPhashmapCreate(&consHidden, SCIPblkmem(scip), 1000) );

   /* loop through indicator constraints (works only in transformed problem) */
   if ( transformed )
   {
      SCIP_CONS** consInd;
      int nConsInd = 0;
      
      consInd = SCIPconshdlrGetConss(conshdlr);
      nConsInd = SCIPconshdlrGetNConss(conshdlr);
      SCIPdebugMessage("Number of indicator constraints: %d\n", nConsInd);

      for( c = 0; c < nConsInd; ++c )
      {
	 assert( consInd[c] != NULL );
	 cons = SCIPgetLinearConsIndicator(consInd[c]);
	 
	 assert( !SCIPhashmapExists(consHidden, (void*) cons) );
	 SCIP_CALL( SCIPhashmapSetImage(consHidden, (void*) cons, (void*) TRUE) );
	 SCIPdebugMessage("Marked linear constraint <%s> as hidden.\n", SCIPconsGetName(cons));
      }
   }
   else
   {
      /* otherwise we have to pass through all constraints */
      for (c = 0; c < nconss; ++c)
      {
	 cons = conss[c];
	 assert( cons != NULL);

	 conshdlr = SCIPconsGetHdlr(cons);
	 assert( conshdlr != NULL );
	 conshdlrname = SCIPconshdlrGetName(conshdlr);

	 if( strcmp(conshdlrname, "indicator") == 0 )
	 {
	    SCIP_CONS* lincons;

	    lincons = SCIPgetLinearConsIndicator(cons);
	    assert( lincons != NULL );

	    assert( !SCIPhashmapExists(consHidden, (void*) lincons) );
	    SCIP_CALL( SCIPhashmapSetImage(consHidden, (void*) lincons, (void*) TRUE) );
	    SCIPdebugMessage("Marked linear constraint <%s> as hidden.\n", SCIPconsGetName(lincons));
	 }
      }
   }
#endif


   /* check if the variable names are not too long */
   checkVarnames(scip, vars, nvars);
   /* check if the constraint names are too long */
   checkConsnames(scip, conss, nconss, transformed);

   /* print statistics as comment to file */
   SCIPinfoMessage(scip, file, "* SCIP STATISTICS\n");
   SCIPinfoMessage(scip, file, "*   Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "*   Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   SCIPinfoMessage(scip, file, "*   Constraints      : %d\n", nconss);
   SCIPinfoMessage(scip, file, "*   Obj. scale       : %.15g\n", objscale);
   SCIPinfoMessage(scip, file, "*   Obj. offset      : %.15g\n\n", objoffset);

   /* print flags */
   SCIPinfoMessage(scip, file, "$MAXCOL %d\n", GMS_MAX_LINELEN - 1);
   SCIPinfoMessage(scip, file, "$OFFDIGIT\n\n");

   /* print variable section */
   SCIPinfoMessage(scip, file, "Variables\n");
   clearLine(linebuffer, &linecnt);

   /* auxiliary objective variable */
   SCIPinfoMessage(scip, file, " objvar,\n");

   /* "model" variables */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );

      (void) SCIPsnprintf(varname, GMS_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s%s", varname, (v < nvars - 1) ? "," : ";");
      appendLine(scip, file, linebuffer, &linecnt, buffer);

      if( (linecnt > 0 && (v == nbinvars - 1 || v == nbinvars + nintvars - 1 ||
	       v == nbinvars + nintvars + nimplvars - 1)) || v == nvars - 1 )
      {
	 endLine(scip, file, linebuffer, &linecnt);
	 clearLine(linebuffer, &linecnt);
      }
   }

   SCIPinfoMessage(scip, file, "\n");

   /* declare binary variables if present */
   if( nbinvars > 0 )
   {
      SCIPinfoMessage(scip, file, "Binary variables\n");
      clearLine(linebuffer, &linecnt);
   
      for( v = 0; v < nbinvars; ++v )
      {
	 var = vars[v];

	 (void) SCIPsnprintf(varname, GMS_MAX_NAMELEN, "%s", SCIPvarGetName(var));
	 (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s%s", varname, (v < nbinvars - 1) ? "," : ";");

	 appendLine(scip, file, linebuffer, &linecnt, buffer);
      }

      endLine(scip, file, linebuffer, &linecnt);
      SCIPinfoMessage(scip, file, "\n");
   }

   /* declare integer variables if present */
   if( nintvars > 0 )
   {
      SCIPinfoMessage(scip, file, "Integer variables\n");
      clearLine(linebuffer, &linecnt);
   
      for( v = 0; v < nintvars; ++v )
      {
	 var = vars[nbinvars + v];

	 (void) SCIPsnprintf(varname, GMS_MAX_NAMELEN, "%s", SCIPvarGetName(var));
	 (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s%s", varname, (v < nintvars - 1) ? "," : ";");

	 appendLine(scip, file, linebuffer, &linecnt, buffer);
      }

      endLine(scip, file, linebuffer, &linecnt);
      SCIPinfoMessage(scip, file, "\n");
   }

   /* print variable bounds */
   SCIPinfoMessage(scip, file, "* Variable bounds\n");
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/gmsreader/freeints", &freeints) );
   nondefbounds = FALSE;

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );
      (void) SCIPsnprintf(varname, GMS_MAX_NAMELEN, "%s", SCIPvarGetName(var));

      if( transformed )
      {
         /* in case the transformed is written only local bounds are posted which are valid in the current node */
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
      }
      else
      {
         lb = SCIPvarGetLbOriginal(var);
         ub = SCIPvarGetUbOriginal(var);
      }
      assert( lb <= ub );

      /* fixed */
      if( SCIPisEQ(scip, lb, ub) )
      {
	 if( v < nintvars )
	    SCIPinfoMessage(scip, file, " %s.fx = %g;\n", varname, SCIPfloor(scip, lb + 0.5));
	 else
	    SCIPinfoMessage(scip, file, " %s.fx = %.15g;\n", varname, lb);
	 nondefbounds = TRUE;
      }

      /* lower bound */
      if( v < nbinvars || (v < nintvars && !freeints) )
      {
	 /* default lower bound is 0 */
	 if( !SCIPisZero(scip, lb) )
	 {
	    if( !SCIPisInfinity(scip, -lb) )
	       SCIPinfoMessage(scip, file, " %s.lo = %g;\n", varname, SCIPceil(scip, lb));
	    else
	       SCIPinfoMessage(scip, file, " %s.lo = %g;\n", varname, -SCIPinfinity(scip)); /* sorry, -inf not allowed in gams file here */
	    nondefbounds = TRUE;
	 }
      }
      else if( v < nintvars && !SCIPisInfinity(scip, -lb) )
      {
	 /* freeints == TRUE: integer variables are free by default */
	 SCIPinfoMessage(scip, file, " %s.lo = %g;\n", varname, SCIPceil(scip, lb));
	 nondefbounds = TRUE;
      }
      else if( v >= nintvars && !SCIPisInfinity(scip, -lb) )
      {
	 /* continuous variables are free by default */
	 SCIPinfoMessage(scip, file, " %s.lo = %.15g;\n", varname, lb);
	 nondefbounds = TRUE;
      }

      /* upper bound */
      if( v < nbinvars )
      {
	 if( !SCIPisEQ(scip, ub, 1.0) )
	 {
	    SCIPinfoMessage(scip, file, " %s.up = %g;\n", varname, SCIPfloor(scip, ub));
	    nondefbounds = TRUE;
	 }
      }
      else if( v < nintvars && !freeints )
      {
	 /* freeints == FALSE: integer variables have upper bound 100 by default */
	 if( !SCIPisEQ(scip, ub, 100.0) )
	 {
	    if( !SCIPisInfinity(scip, ub) )
	       SCIPinfoMessage(scip, file, " %s.up = %g;\n", varname, SCIPfloor(scip, ub));
	    else
	       SCIPinfoMessage(scip, file, " %s.up = %g;\n", varname, SCIPinfinity(scip)); /* sorry, +inf not allowed in gams file here */
	    nondefbounds = TRUE;
	 }
      }
      else if( v < nintvars && !SCIPisInfinity(scip, ub) )
      {
	 /* freeints == TRUE: integer variables are free by default */
	 SCIPinfoMessage(scip, file, " %s.up = %g;\n", varname, SCIPfloor(scip, ub));
	 nondefbounds = TRUE;
      }
      else if( v >= nintvars && !SCIPisInfinity(scip, ub) )
      {
	 /* continuous variables are free by default */
	 SCIPinfoMessage(scip, file, " %s.up = %.15g;\n", varname, ub);
	 nondefbounds = TRUE;
      }
   }

   if( !nondefbounds )
      SCIPinfoMessage(scip, file, "* (All other bounds at default value: binary [0,1], integer [%s], continuous [-inf,+inf].)\n", freeints ? "-inf,+inf" : "0,100");
   SCIPinfoMessage(scip, file, "\n");

   /* print equations section */
   SCIPinfoMessage(scip, file, "Equations\n");
   clearLine(linebuffer, &linecnt);

   SCIPinfoMessage(scip, file, " objequ%s\n", (nconss > 0) ? "," : ";");

   /* declare equations */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL );

      /* we declare all constraints, although we might not define each of them later */
      (void) SCIPsnprintf(consname, GMS_MAX_NAMELEN, "%s", SCIPconsGetName(cons));
      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s%s", consname, (c < nconss - 1) ? "," : ";");
      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }
      
   endLine(scip, file, linebuffer, &linecnt);
   SCIPinfoMessage(scip, file, "\n");
   
   /* print objective function equation */
   clearLine(linebuffer, &linecnt);
   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " objequ .. objvar =e= ");
   appendLine(scip, file, linebuffer, &linecnt, buffer);

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );

#ifndef NDEBUG
      /* in case the original problem has to be posted the variables have to be either "original" or "negated" */
      if ( !transformed )
         assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL ||
            SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED );
#endif
      
      if( SCIPisZero(scip, SCIPvarGetObj(var)) )
         continue;

      if( linecnt == 0 )
         /* we start a new line; therefore we tab this line */
         appendLine(scip, file, linebuffer, &linecnt, "     ");

      (void) SCIPsnprintf(varname, GMS_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %+.15g*%s%s", SCIPvarGetObj(var), varname, v == nvars - 1 ? ";" : "");

      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   endLine(scip, file, linebuffer, &linecnt);
   SCIPinfoMessage(scip, file, "\n");

   /* print constraints */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL );

      /* in case the transformed is written, only constraints are posted which are enabled in the current node */
      if( transformed && !SCIPconsIsEnabled(cons) )
         continue;

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      (void) SCIPsnprintf(consname, GMS_MAX_NAMELEN, "%s", SCIPconsGetName(cons));
      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_CALL( printLinearCons(scip, file, consname,
               SCIPgetNVarsLinear(scip, cons), SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
               SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), transformed) );
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch ( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( printLinearCons(scip, file, consname,
                  nconsvar, consvars, NULL, 1.0, 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( printLinearCons(scip, file, consname,
                  nconsvars, consvars, NULL, -SCIPinfinity(scip), 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( printLinearCons(scip, file, consname,
                  nconsvars, consvars, NULL, 1.0, SCIPinfinity(scip), transformed) );
            break;
         }
      }
      else if ( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( printLinearCons(scip, file, consname,
               SCIPgetNVarsLogicor(scip, cons), SCIPgetVarsLogicor(scip, cons), NULL,
               1.0, SCIPinfinity(scip), transformed) );
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

         SCIP_CALL( printLinearCons(scip, file, consname, nconsvars, consvars, consvals,
	       -SCIPinfinity(scip), (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), transformed) );

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

         SCIP_CALL( printLinearCons(scip, file, consname,
               2, consvars, consvals,
               SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons), transformed) );

         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else
      {
         SCIPwarningMessage("constraint handler <%s> can not print requested format\n", conshdlrname );
         SCIPinfoMessage(scip, file, "* ");
         SCIP_CALL( SCIPprintCons(scip, cons, file) );
         SCIPinfoMessage(scip, file, "\n");
      }

      SCIPinfoMessage(scip, file, "\n");
   }

   /* print model creation */
   SCIPinfoMessage(scip, file, "Model m / all /;\n\n");

   /* print solve command */
   SCIPinfoMessage(scip, file, "$if not set MIP $set MIP MIP\n");
   SCIPinfoMessage(scip, file, "Solve m using %%MIP%% %simizing objvar;\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "min" : "max");

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
