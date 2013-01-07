/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
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

/**@file   probdata_lop.h
 * @brief  handling of data needed for solving linear ordering problems
 * @author Marc Pfetsch
 */

#include "probdata_lop.h"

#include "cons_linearordering.h"
#include <scip/misc.h>


struct SCIP_ProbData
{
   int                   n;                  /**< number of elements */
   SCIP_Real**           W;                  /**< weight matrix */
   SCIP_VAR***           vars;               /**< variables */
};



/* ----------------- SCIP interface functions ------------------------ */

/** delete problem data */
static
SCIP_DECL_PROBDELORIG(probdelorigLOP)
{  /*lint --e{831} */
   int i, j;

   assert( probdata != NULL );
   assert( *probdata != NULL );

   /* free matrix and release and free variables */
   assert( (*probdata)->W != NULL );
   assert( (*probdata)->vars != NULL );
   for (i = 0; i < (*probdata)->n; ++i)
   {
      for (j = 0; j < (*probdata)->n; ++j)
      {
	 if (j != i)
	    SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i][j]) );
      }
      SCIPfreeMemoryArray(scip, &(*probdata)->vars[i]);
      SCIPfreeMemoryArray(scip, &((*probdata)->W[i]));
   }
   SCIPfreeMemoryArray(scip, &(*probdata)->vars);
   SCIPfreeMemoryArray(scip, &((*probdata)->W));

   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


/** copies user data of source SCIP for the target SCIP */
static
SCIP_DECL_PROBCOPY(probcopyLOP)
{
   int n;
   int i;
   int j;

   assert( scip != NULL );
   assert( sourcescip != NULL );
   assert( sourcedata != NULL );
   assert( targetdata != NULL );

   /* set up data */
   SCIP_CALL( SCIPallocMemory(scip, targetdata) );

   n = sourcedata->n;
   (*targetdata)->n = n;

   /* set matrices */
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->W), n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->vars), n) );

   for( i = 0; i < n; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->W[i]), n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->vars[i]), n) );

      for( j = 0; j < n; ++j )
      {
         if( i != j )
         {
            SCIP_VAR* var;
            SCIP_Bool success;
            
	    SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->vars[i][j], &var) );
            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->vars[i][j]), varmap, consmap, global, &success) );
            assert(success);
            assert((*targetdata)->vars[i][j] != NULL);

            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->vars[i][j]) );
         }
         else
            (*targetdata)->vars[i][j] = NULL;
      }
   }

   *result = SCIP_SUCCESS;
   
   return SCIP_OKAY;
}




/* ----------------- auxiliary functions ------------------------ */

/** read weight matrix from file (in LOLIB format)
 *
 *  Format:
 *  comment line
 *  # of elements
 *  weight matrix (doubles)
 */
static
SCIP_RETCODE LOPreadFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of file to read */
   SCIP_PROBDATA*        probdata            /**< problem data to be filled */
   )
{
   int i, j;
   FILE *file;
   int status;
   int n;
   SCIP_Real** W;
   char s[SCIP_MAXSTRLEN];

   /* open file */
   file = fopen(filename, "r");
   if ( file == NULL )
   {
      SCIPerrorMessage("Could not open file %s.\n", filename);
      return SCIP_NOFILE;
   }

   /* skip one line */
   fgets(s, SCIP_MAXSTRLEN, file);

   /* read number of elements */
   status = fscanf(file, "%d", &n);
   if ( ! status )
   {
      SCIPerrorMessage("Reading failed.\n");
      return SCIP_READERROR;
   }
   assert( 0 < n );
   SCIPinfoMessage(scip, NULL, "Number of elements:\t%d\n\n", n);
   probdata->n = n;

   /* set up matrix */
   SCIP_CALL( SCIPallocMemoryArray(scip, &W, n) );
   for (i = 0; i < n; ++i)
      SCIP_CALL( SCIPallocMemoryArray(scip, &(W[i]), n) );
   probdata->W = W;

   /* read matrix */
   for (i = 0; i < n; ++i)
   {
      for (j = 0; j < n; ++j)
      {
	 SCIP_Real val;
	 status = fscanf(file, "%lf", &val);
	 if ( ! status )
	 {
	    SCIPerrorMessage("Reading failed.\n");
	    return SCIP_READERROR;
	 }
	 W[i][j] = val;
      }
   }
   fclose( file );

   return SCIP_OKAY;
}


/** get problem name
 *
 *  Returns NULL on error
 */
static
SCIP_RETCODE getProblemName(
   const char*           filename,           /**< input filename */
   char*                 probname,           /**< output problemname */
   int                   maxSize             /**< maximum size of probname */
   )
{
   int i = 0;
   int j = 0;
   int l;

   /* first find end of string */
   while ( filename[i] != 0)
      ++i;
   l = i;

   /* go back until '.' or '/' or '\' appears */
   while ((i > 0) && (filename[i] != '.') && (filename[i] != '/') && (filename[i] != '\\'))
      --i;

   /* if we found '.', search for '/' or '\\' */
   if (filename[i] == '.')
   {
      l = i;
      while ((i > 0) && (filename[i] != '/') && (filename[i] != '\\'))
	 --i;
   }

   /* correct counter */
   if ((filename[i] == '/') || (filename[i] == '\\'))
      ++i;

   /* copy name */
   while ( (i < l) && (filename[i] != 0) )
   {
      probname[j++] = filename[i++];
      if (j > maxSize-1)
	 return SCIP_ERROR;
   }
   probname[j] = 0;

   return SCIP_OKAY;
}



/* ----------------- public interface functions ------------------------ */


/** create linear ordering problem instance */
SCIP_RETCODE LOPcreateProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of file to read */
   )
{
   SCIP_PROBDATA* probdata = NULL;
   char probname[SCIP_MAXSTRLEN];

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, &probdata) );

   /* take filename as problem name */
   SCIP_CALL( getProblemName(filename, probname, SCIP_MAXSTRLEN) );

   SCIPinfoMessage(scip, NULL, "File name:\t\t%s\n", filename);
   SCIPinfoMessage(scip, NULL, "Problem name:\t\t%s\n", probname);

   /* read file */
   SCIP_CALL( LOPreadFile(scip, filename, probdata) );
   probdata->vars = NULL;

   SCIP_CALL( SCIPcreateProb(scip, probname, probdelorigLOP, NULL, NULL,
	 NULL, NULL, probcopyLOP, probdata) );

   return SCIP_OKAY;
}


/** create linear ordering problem model */
SCIP_RETCODE LOPgenerateModel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS* cons;
   int i, j;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   /* generate variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->vars, probdata->n) );
   for (i = 0; i < probdata->n; ++i)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->vars[i]), probdata->n) );
      for (j = 0; j < probdata->n; ++j)
      {
	 if (j != i)
	 {
	    char s[SCIP_MAXSTRLEN];
	    (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "x#%d#%d", i, j);
	    SCIP_CALL( SCIPcreateVar(scip, &(probdata->vars[i][j]), s, 0.0, 1.0, probdata->W[i][j], SCIP_VARTYPE_BINARY,
		  TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
	    SCIP_CALL( SCIPaddVar(scip, probdata->vars[i][j]) );
	 }
	 else
	    probdata->vars[i][j] = NULL;
      }
   }

   /* generate linear ordering constraint */
   SCIP_CALL( SCIPcreateConsLinearOrdering(scip, &cons, "LOP", probdata->n, probdata->vars, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
	 FALSE, FALSE, FALSE, FALSE));
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* set maximization */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   return SCIP_OKAY;
}


/** evalutate solution */
SCIP_RETCODE LOPevalSolution(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_VAR*** vars;
   SCIP_SOL* sol;
   int* outDegree;
   int* indices;
   int i;
   int j;
   int n;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   assert( probdata->vars != NULL );

   n = probdata->n;
   vars = probdata->vars;
   sol = SCIPgetBestSol(scip);

   if ( sol == NULL )
      printf("\nNo solution found.\n");
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &outDegree, n) );
      SCIP_CALL( SCIPallocBufferArray(scip, &indices, n) );

      /* compute out-degree */
      for (i = 0; i < n; ++i)
      {
	 int deg = 0;
	 for (j = 0; j < n; ++j)
	 {
	    SCIP_Real val = 0.0;
	    if (j == i)
	       continue;

	    val = SCIPgetSolVal(scip, sol, vars[i][j]);
	    assert( SCIPisIntegral(scip, val) );
	    if ( val < 0.5 )
	       ++deg;
	 }
	 outDegree[i] = deg;
	 indices[i] = i;
      }

      /* sort such that degrees are non-decreasing */
      SCIPsortIntInt(outDegree, indices, n);

      /* output */
      printf("\nFinal order:\n");
      for (i = 0; i < n; ++i)
	 printf("%d ", indices[i]);
      printf("\n");

      SCIPfreeBufferArray(scip, &indices);
      SCIPfreeBufferArray(scip, &outDegree);
   }

   return SCIP_OKAY;
}



/** return the number of elements */
int LOPgetNElements(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   return probdata->n;
}
