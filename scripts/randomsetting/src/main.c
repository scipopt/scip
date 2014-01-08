/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic Licence.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   main.c
 * @brief  random settings generator
 * @author Ambros Gleixner
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "scip/scip.h"
#include "scip/scipgithash.h"
#include "scip/scipdefplugins.h"

/** checks whether parameter is in the list of excluded parameters
 *
 *  @note this should be parameters for which SCIP most of the time aborts with a proper error in optimized mode;
 *        ideally there should be no random parameter settings which cause an assert in debug mode and undefined
 *        behaviour in optimized mode
 */
#define MAX_PARAMNAME_LEN 255
static
SCIP_Bool paramIsExcluded(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param               /**< pointer to parameter */
   )
{
   const char* paramname;
   SCIP_Bool excluded;

   assert(scip != NULL);
   assert(param != NULL);

   paramname = SCIPparamGetName(param);
   excluded = FALSE;

   /* otherwise the readers do not work */
   excluded = excluded
      || strncmp(paramname, "misc/usevartable", MAX_PARAMNAME_LEN) == 0
      || strncmp(paramname, "misc/useconstable", MAX_PARAMNAME_LEN) == 0;

   /* too large initial values easily yield a memory out */
   excluded = excluded
      || strncmp(paramname, "memory/arraygrowinit", MAX_PARAMNAME_LEN) == 0
      || strncmp(paramname, "memory/treegrowinit", MAX_PARAMNAME_LEN) == 0
      || strncmp(paramname, "memory/pathgrowinit", MAX_PARAMNAME_LEN) == 0;

   /* already larger values inside the feasible interval [1,10] yield many memory outs */
   excluded = excluded
      || strncmp(paramname, "memory/arraygrowfac", MAX_PARAMNAME_LEN) == 0
      || strncmp(paramname, "memory/treegrowfac", MAX_PARAMNAME_LEN) == 0
      || strncmp(paramname, "memory/pathgrowfac", MAX_PARAMNAME_LEN) == 0;

   if( excluded )
   {
      SCIPdebugMessage("excluding parameter <%s>\n", SCIPparamGetName(param));
      return TRUE;
   }

   return FALSE;
}

/** sets each parameter to a random value within its domain */
static
SCIP_RETCODE generateRandomSettings(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int*         seedp,              /**< pointer to seed value */
   SCIP_Bool             advanced            /**< should advanced parameters also be changed? */
   )
{
   SCIP_PARAM** params = SCIPgetParams(scip);
   int nparams = SCIPgetNParams(scip);
   int i;

   assert(scip != NULL);
   assert(seedp != NULL);
   assert(params != NULL);
   assert(nparams >= 0);

   for( i = 0; i < nparams; ++i )
   {
      SCIP_PARAM* param = params[i];

      assert(param != NULL);

      if( advanced == FALSE && SCIPparamIsAdvanced(param) )
         continue;

      if( paramIsExcluded(scip, param) )
         continue;

      switch( SCIPparamGetType(param) )
      {
      case SCIP_PARAMTYPE_BOOL:
         {
            SCIP_RETCODE retcode;
            int value = SCIPgetRandomInt(0, 1, seedp);
            assert(value >= 0);
            assert(value <= 1);

            SCIPdebugMessage("changing bool parameter <%s> to %d\n", SCIPparamGetName(param), value);

            retcode = SCIPchgBoolParam(scip, param, (SCIP_Bool)value);
            if( retcode != SCIP_OKAY && retcode != SCIP_PARAMETERWRONGVAL )
            {
               SCIP_CALL( retcode );
            }

            break;
         }

      case SCIP_PARAMTYPE_INT:
         {
            SCIP_RETCODE retcode;
            int value = SCIPgetRandomInt(SCIPparamGetIntMin(param), SCIPparamGetIntMax(param), seedp);
            assert(value >= SCIPparamGetIntMin(param));
            assert(value <= SCIPparamGetIntMax(param));

            SCIPdebugMessage("changing int parameter <%s> to %d\n", SCIPparamGetName(param), value);

            retcode = SCIPchgIntParam(scip, param, value);
            if( retcode != SCIP_OKAY && retcode != SCIP_PARAMETERWRONGVAL )
            {
               SCIP_CALL( retcode );
            }

            break;
         }

      case SCIP_PARAMTYPE_LONGINT:
         {
            SCIP_RETCODE retcode;
            int minvalue = (SCIPparamGetLongintMin(param) >= (SCIP_Longint)INT_MIN)
               ? (int)SCIPparamGetLongintMin(param)
               : INT_MIN;
            int maxvalue = (SCIPparamGetLongintMax(param) <= (SCIP_Longint)INT_MAX)
               ? (int)SCIPparamGetLongintMax(param)
               : INT_MAX;
            SCIP_Longint value = (SCIP_Longint)SCIPgetRandomInt(minvalue, maxvalue, seedp);
            assert(value >= SCIPparamGetLongintMin(param));
            assert(value <= SCIPparamGetLongintMax(param));

            SCIPdebugMessage("changing longint parameter <%s> to %"SCIP_LONGINT_FORMAT"\n", SCIPparamGetName(param), value);

            retcode = SCIPchgLongintParam(scip, param, value);
            if( retcode != SCIP_OKAY && retcode != SCIP_PARAMETERWRONGVAL )
            {
               SCIP_CALL( retcode );
            }

            break;
         }

      case SCIP_PARAMTYPE_REAL:
         {
            SCIP_RETCODE retcode;
            SCIP_Real value = SCIPgetRandomReal(SCIPparamGetRealMin(param), SCIPparamGetRealMax(param), seedp);
            assert(value >= SCIPparamGetRealMin(param));
            assert(value <= SCIPparamGetRealMax(param));

            SCIPdebugMessage("changing real parameter <%s> to %g\n", SCIPparamGetName(param), value);

            retcode = SCIPchgRealParam(scip, param, value);
            if( retcode != SCIP_OKAY && retcode != SCIP_PARAMETERWRONGVAL )
            {
               SCIP_CALL( retcode );
            }
            break;
         }

      case SCIP_PARAMTYPE_CHAR:
         {
            SCIP_RETCODE retcode;
            char* allowedvalues;
            char value;

            allowedvalues = SCIPparamGetCharAllowedValues(param);
            if( allowedvalues == NULL )
            {
               value = (char)SCIPgetRandomInt((int)CHAR_MIN, (int)CHAR_MAX, seedp);
            }
            else
            {
               int length;
               int pos;

               length = strlen(allowedvalues);
               assert(length >= 1);

               pos = SCIPgetRandomInt(0, length - 1, seedp);
               value = allowedvalues[pos];
            }

            SCIPdebugMessage("changing char parameter <%s> to %c\n", SCIPparamGetName(param), value);

            retcode = SCIPchgCharParam(scip, param, value);
            if( retcode != SCIP_OKAY && retcode != SCIP_PARAMETERWRONGVAL )
            {
               SCIP_CALL( retcode );
            }
            break;
         }
      case SCIP_PARAMTYPE_STRING:
         /* we currently do net test string parameters; this would need to be done manually */
         break;

      default:
         SCIPABORT();
         break;
      }
   }

   return SCIP_OKAY;
}

/** creates a settings file with random parameters */
static
SCIP_RETCODE createSettingsFile(
   unsigned int          seed                /**< random seed */
   )
{
   SCIP* scip = NULL;
   char filename[256];
   int status;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIPdebug( SCIPprintVersion(scip, NULL) );

   SCIP_CALL( generateRandomSettings(scip, &seed, TRUE) );

   /* write settings file */
   status = snprintf(filename, 255, "%s-%u.set", SCIPgetGitHash(), seed);
   if( status < 0 || status >= 256 )
   {
      SCIPerrorMessage("error creating filename: snprintf() returned %d\n", status);
      return SCIP_ERROR;
   }
   else
   {
      SCIP_CALL( SCIPwriteParams(scip, filename, TRUE, TRUE) );
      SCIPinfoMessage(scip, NULL, "-> saved random settings with seed %ld to file <%s>\n",
         seed, filename);
   }

   /* close SCIP */
   SCIP_CALL( SCIPfree(&scip) );
   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method  */
int main(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )
{
   SCIP_RETCODE retcode = SCIP_ERROR;

   if( argc >= 2 )
   {
      /* because of the filename we require a non-negative seed */
      if( strtol(argv[1], NULL, 10) < 0 )
      {
         SCIPerrorMessage("error creating filename: seed value is negative\n");
         return -1;
      }

      retcode = createSettingsFile((unsigned int)strtol(argv[1], NULL, 10));
   }
   else
      retcode = createSettingsFile((unsigned int)time(NULL));

   if( retcode != SCIP_OKAY )
   {
      printf("usage: %s [<seed>]\n\n", argv[0]);
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
