/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: probdata_LOP.c,v 1.3 2007/10/01 20:03:08 bzfpfets Exp $"

#include "probdata_LOP.h"

#include "cons_LO.h"
#include <string.h>
#include <libgen.h>


struct SCIP_ProbData
{
   int n;             /**< number of elements */
   SCIP_Real** W;     /**< weight matrix */
   SCIP_VAR*** Vars;  /**< variables */
};



/* ----------------- SCIP interface functions ------------------------ */

/** delete problem data */
static
SCIP_DECL_PROBDELORIG(probdelorigLOP)
{
   int i, j;

   assert(probdata != NULL);
   assert(*probdata != NULL);

   /* free matrix and release and free variables */
   assert( (*probdata)->W != NULL );
   assert( (*probdata)->Vars != NULL );
   for (i = 0; i < (*probdata)->n; ++i)
   {
      for (j = 0; j < (*probdata)->n; ++j)
      {
	 if (j != i)
	    SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->Vars[i][j]) );
      }
      SCIPfreeMemoryArray(scip, &(*probdata)->Vars[i]);
      SCIPfreeMemoryArray(scip, &((*probdata)->W[i]));
   }
   SCIPfreeMemoryArray(scip, &(*probdata)->Vars);
   SCIPfreeMemoryArray(scip, &((*probdata)->W));

   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


#define probtransLOP NULL
#define probdeltransLOP NULL
#define probinitsolLOP NULL
#define probexitsolLOP NULL




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
   SCIP*        scip,          /**< SCIP data structure */
   const char*  filename,      /**< name of file to read */
   SCIP_PROBDATA* probdata     /**< problem data to be filled */
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
   SCIPmessagePrintInfo("Number of elements: %d\n\n", n);
   probdata->n = n;

   /* set up matrix */
   SCIP_CALL( SCIPallocMemoryArray(scip, &W, n) );
   for (i = 0; i < n; ++i)
      SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &(W[i]), n) );
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





/* ----------------- outside interface functions ------------------------ */

/** create linear ordering problem instance */
SCIP_RETCODE LOPcreateProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of file to read */
   )
{
   SCIP_PROBDATA* probdata = NULL;
   char* filenameCopy;
   char* probname;

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, &probdata) );

   /* take filename as problem name */
   filenameCopy = strdup(filename);     /* need copy since some implementations of basename modify filename */
   probname = basename(filenameCopy);

   SCIPmessagePrintInfo("Problem name: %s\n\n", probname);

   /* read file */
   SCIP_CALL( LOPreadFile(scip, filename, probdata) );
   probdata->Vars = NULL;

   SCIP_CALL( SCIPcreateProb(scip, probname, probdelorigLOP, probtransLOP, probdeltransLOP,
			     probinitsolLOP, probexitsolLOP, probdata) );

   free((char*) filenameCopy);

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
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->Vars, probdata->n) );
   for (i = 0; i < probdata->n; ++i)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->Vars[i]), probdata->n) );
      for (j = 0; j < probdata->n; ++j)
      {
	 if (j != i)
	 {
	    char s[SCIP_MAXSTRLEN];
	    sprintf(s, "x#%d#%d", i, j);
	    SCIP_CALL( SCIPcreateVar(scip, &(probdata->Vars[i][j]), s, 0.0, 1.0, probdata->W[i][j], SCIP_VARTYPE_BINARY,
				     TRUE, FALSE, NULL, NULL, NULL, NULL));
	    SCIP_CALL( SCIPaddVar(scip, probdata->Vars[i][j]) );
	 }
	 else
	    probdata->Vars[i][j] = NULL;
      }
   }

   /* generate linear ordering constraint */
   SCIP_CALL( LOcreateCons(scip, &cons, "LOP", probdata->n, probdata->Vars, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE));
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* set maximization */
   SCIP_CALL_ABORT( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   return SCIP_OKAY;
}


/** evalutate solution */
SCIP_RETCODE LOPevalSolution(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   int i, j, n;
   SCIP_VAR*** Vars;
   SCIP_SOL* sol;
   SCIP_Real* inDegree;
   int* indices;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   assert( probdata->Vars != NULL );

   n = probdata->n;
   Vars = probdata->Vars;
   sol = SCIPgetBestSol(scip);

   if ( sol == NULL )
      printf("No solution found.\n");
   else
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &inDegree, n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &indices, n) );
      
      /* compute in-degree */
      for (i = 0; i < n; ++i)
      {
	 int deg = 0;
	 for (j = 0; j < n; ++j)
	 {
	    if (j == i)
	       continue;
	    
	    SCIP_Real val = 0.0;
	    val = SCIPgetSolVal(scip, sol, Vars[i][j]);
	    assert( SCIPisIntegral(scip, val) );
	    if ( val > 0.5 )
	       ++deg;
	 }
	 inDegree[i] = (SCIP_Real) deg;
	 indices[i] = i;
      }

      /* sort such that degrees are non-decreasing */
      SCIPbsortRealPtr(inDegree, (void**) indices, n);
      
      /* output */
      printf("\nFinal order:\n");
      for (i = 0; i < n; ++i)
	 printf("%d ", indices[i]);
      printf("\n");
      
      SCIPfreeMemoryArray(scip, &inDegree);
      SCIPfreeMemoryArray(scip, &indices);
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
