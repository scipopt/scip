/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  main file for linear ordering example
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include "probdata_lop.h"
#include "cons_linearordering.h"

/** define macro to print error message and exit */
#define SCIP_CALL_ERROR(x)   do                                                                               \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             SCIPprintError(_restat_);                                                        \
                             return -1;                                                                       \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )


/** read comand line arguments */
static
SCIP_Bool readArguments(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char**          filename,           /**< file name from arguments */
   const char**          settingsname,       /**< name of settings file */
   SCIP_Real*            timelimit,          /**< time limit read from arguments */
   SCIP_Real*            memlimit,           /**< memory limit read from arguments */
   SCIP_Longint*         nodelimit,          /**< node limit read from arguments */
   int*                  dispfreq            /**< display frequency */
   )
{
   int i;
   char usage[256];
#ifndef NDEBUG
   int status;

   assert( argc > 0 );
   assert( argv != NULL );
   assert( filename != NULL );
   assert( settingsname != NULL );
   assert( timelimit != NULL );
   assert( memlimit != NULL );
   assert( nodelimit != NULL );
   assert( dispfreq != NULL );

   /* init usage text */
   status = snprintf(usage, 255, "usage: %s <LOP file> [-s <setting file>] [-t <time limit>] [-m <mem limit>] [-n <node limit>] [-d <display frequency>]", argv[0]);
   assert( 0 <= status && status < 256 );
#else
   (void) snprintf(usage, 255, "usage: %s <LOP file> [-s <setting file>] [-t <time limit>] [-m <mem limit>] [-n <node limit>] [-d <display frequency>]", argv[0]);
#endif

   /* init arguments */
   *timelimit = 1e20;
   *memlimit = 1e20;
   *nodelimit = SCIP_LONGINT_MAX;
   *filename = NULL;
   *settingsname = NULL;
   *dispfreq = -1;

   /* check all arguments */
   for (i = 1; i < argc; ++i)
   {
      if ( ! strcmp(argv[i], "-s") )
      {
         if ( *settingsname != NULL )
         {
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         if ( i == argc-1 )
         {
            SCIPinfoMessage(scip, NULL, "Error: No setting file name supplied.\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         ++i;
         *settingsname = argv[i];
         assert( i < argc );
      }
      /* check for time limit */
      else if ( ! strcmp(argv[i], "-t") )
      {
         if ( i == argc-1 )
         {
            SCIPinfoMessage(scip, NULL, "No time limit supplied.\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         ++i;
         *timelimit = atof(argv[i]);
         assert( i < argc );
      }
      /* check for memory limit */
      else if ( ! strcmp(argv[i], "-m") )
      {
         if ( i == argc-1 )
         {
            SCIPinfoMessage(scip, NULL, "No memory limit supplied.\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         ++i;
         *memlimit = atof(argv[i]);
         assert( i < argc );
      }
      /* check for node limit */
      else if ( ! strcmp(argv[i], "-n") )
      {
         if ( i == argc-1 )
         {
            SCIPinfoMessage(scip, NULL, "No node limit supplied.\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         ++i;
         *nodelimit = atol(argv[i]);
         assert( i < argc );
      }
      /* check for display frequency */
      else if ( ! strcmp(argv[i], "-d") )
      {
         if ( i == argc-1 )
         {
            SCIPinfoMessage(scip, NULL, "No display frequency supplied.\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         ++i;
         *dispfreq = atoi(argv[i]);
         assert( i < argc );
      }
      else
      {
         /* if filename is already specified */
         if ( *filename != NULL )
         {
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         *filename = argv[i];
      }
   }

   if ( *filename == NULL ) /*lint !e850*/
   {
      SCIPinfoMessage(scip, NULL, "No filename supplied.\n");
      SCIPinfoMessage(scip, NULL, "%s\n", usage);
      return FALSE;
   }

   return TRUE;
}


/** main function, which starts the solution of the linear ordering problem */
int main(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )
{
   SCIP* scip = NULL;
   const char* filename;
   const char* settingsname;
   SCIP_Real timelimit;
   SCIP_Real memlimit;
   SCIP_Longint nodelimit;
   int dispfreq;
   int errorcode;

   /* initialize SCIP */
   SCIP_CALL_ERROR( SCIPcreate(&scip) );

   /* parse command line arguments */
   if( readArguments(scip, argc, argv, &filename, &settingsname, &timelimit, &memlimit, &nodelimit, &dispfreq) )
   {
      assert( filename != NULL );

      errorcode = 0;

      /* output version information */
      SCIPinfoMessage(scip, NULL, "Solving the linear ordering problem using SCIP.\n");
      SCIPinfoMessage(scip, NULL, "\n");

      SCIPprintVersion(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");

      /* include default SCIP plugins */
      SCIP_CALL_ERROR( SCIPincludeDefaultPlugins(scip) );

      /* include linear ordering constraint handler */
      SCIP_CALL_ERROR( SCIPincludeConshdlrLinearOrdering(scip) );

      /* set time, node, and memory limit */
      if ( ! SCIPisInfinity(scip, timelimit) )
      {
         SCIPinfoMessage(scip, NULL, "parameter <limits/time> set to %g\n", timelimit);
         SCIP_CALL_ERROR( SCIPsetRealParam(scip, "limits/time", timelimit) );
      }
      if ( ! SCIPisInfinity(scip, memlimit) )
      {
         SCIP_CALL_ERROR( SCIPsetRealParam(scip, "limits/memory", memlimit) );
      }
      if ( nodelimit < SCIP_LONGINT_MAX )
      {
         SCIP_CALL_ERROR( SCIPsetLongintParam(scip, "limits/nodes", nodelimit) );
      }
      if ( dispfreq >= 0 )
      {
         SCIP_CALL_ERROR( SCIPsetIntParam(scip, "display/freq", dispfreq) );
      }

      /* check for parameters */
      if ( settingsname != NULL )
      {
         if ( SCIPfileExists(settingsname) )
         {
            SCIPinfoMessage(scip, NULL, "reading parameter file <%s> ...\n\n", settingsname);
            SCIP_CALL_ERROR( SCIPreadParams(scip, settingsname) );
         }
         else
         {
            SCIPwarningMessage(scip, "parameter file <%s> not found - using default parameters.\n", settingsname);
         }
      }

      /* read problem data */
      SCIP_CALL_ERROR( LOPcreateProb(scip, argv[1]) );

      /* generate linear ordering model */
      SCIP_CALL_ERROR( LOPgenerateModel(scip) );

      /* print model */
      if ( LOPgetNElements(scip) <= 10 )
      {
         SCIP_CALL_ERROR( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );
         SCIPinfoMessage(scip, NULL, "\n");
      }

      /* solve the model */
      SCIP_CALL_ERROR( SCIPsolve(scip) );

      /* print statistics */
      SCIP_CALL_ERROR( SCIPprintStatistics(scip, NULL) );

      /* SCIP_CALL_ERROR( SCIPprintBestSol(scip, NULL, FALSE) );*/
      SCIP_CALL_ERROR( LOPevalSolution(scip) );
   }
   else
      errorcode = 1;

   SCIP_CALL_ERROR( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return errorcode;
}
