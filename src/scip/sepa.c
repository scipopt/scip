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

/**@file   sepa.c
 * @brief  methods and datastructures for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "sepa.h"
#include "clock.h"


/** separators data */
struct Sepa
{
   char*            name;               /**< name of separator */
   char*            desc;               /**< description of separator */
   int              priority;           /**< priority of the separator */
   int              freq;               /**< frequency for calling separator */
   DECL_SEPAFREE    ((*sepafree));      /**< destructor of separator */
   DECL_SEPAINIT    ((*sepainit));      /**< initialise separator */
   DECL_SEPAEXIT    ((*sepaexit));      /**< deinitialise separator */
   DECL_SEPAEXEC    ((*sepaexec));      /**< execution method of separator */
   SEPADATA*        sepadata;           /**< separators local data */
   CLOCK*           clock;              /**< separation time */
   Longint          ncalls;             /**< number of times, this separator was called */
   Longint          ncutsfound;         /**< number of cutting planes found so far by this separator */
   unsigned int     initialized:1;      /**< is separator initialized? */
};



/** creates a separator */
RETCODE SCIPsepaCreate(
   SEPA**           sepa,               /**< pointer to separator data structure */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of separator */
   const char*      desc,               /**< description of separator */
   int              priority,           /**< priority of the separator */
   int              freq,               /**< frequency for calling separator */
   DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   DECL_SEPAINIT    ((*sepainit)),      /**< initialise separator */
   DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialise separator */
   DECL_SEPAEXEC    ((*sepaexec)),      /**< execution method of separator */
   SEPADATA*        sepadata            /**< separator data */
   )
{
   char paramname[MAXSTRLEN];

   assert(sepa != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(sepaexec != NULL);

   ALLOC_OKAY( allocMemory(sepa) );
   ALLOC_OKAY( duplicateMemoryArray(&(*sepa)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*sepa)->desc, desc, strlen(desc)+1) );
   (*sepa)->priority = priority;
   (*sepa)->freq = freq;
   (*sepa)->sepafree = sepafree;
   (*sepa)->sepainit = sepainit;
   (*sepa)->sepaexit = sepaexit;
   (*sepa)->sepaexec = sepaexec;
   (*sepa)->sepadata = sepadata;
   CHECK_OKAY( SCIPclockCreate(&(*sepa)->clock, SCIP_CLOCKTYPE_DEFAULT) );
   (*sepa)->ncalls = 0;
   (*sepa)->ncutsfound = 0;
   (*sepa)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "separator/%s/freq", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, 
                  paramname, "frequency for calling separator (-1: never, 0: only in root node)",
                  &(*sepa)->freq, freq, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of separator */
RETCODE SCIPsepaFree(
   SEPA**           sepa,               /**< pointer to separator data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(sepa != NULL);
   assert(*sepa != NULL);
   assert(!(*sepa)->initialized);
   assert(scip != NULL);

   /* call destructor of separator */
   if( (*sepa)->sepafree != NULL )
   {
      CHECK_OKAY( (*sepa)->sepafree(scip, *sepa) );
   }

   SCIPclockFree(&(*sepa)->clock);
   freeMemoryArray(&(*sepa)->name);
   freeMemoryArray(&(*sepa)->desc);
   freeMemory(sepa);

   return SCIP_OKAY;
}

/** initializes separator */
RETCODE SCIPsepaInit(
   SEPA*            sepa,               /**< separator */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(sepa != NULL);
   assert(scip != NULL);

   if( sepa->initialized )
   {
      char s[MAXSTRLEN];
      sprintf(s, "separator <%s> already initialized", sepa->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( sepa->sepainit != NULL )
   {
      CHECK_OKAY( sepa->sepainit(scip, sepa) );
   }

   SCIPclockReset(sepa->clock);

   sepa->ncalls = 0;
   sepa->ncutsfound = 0;
   sepa->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of separator */
RETCODE SCIPsepaExit(
   SEPA*            sepa,               /**< separator */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(sepa != NULL);
   assert(scip != NULL);

   if( !sepa->initialized )
   {
      char s[MAXSTRLEN];
      sprintf(s, "separator <%s> not initialized", sepa->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( sepa->sepaexit != NULL )
   {
      CHECK_OKAY( sepa->sepaexit(scip, sepa) );
   }
   sepa->initialized = FALSE;

   return SCIP_OKAY;
}

/** calls execution method of separator */
RETCODE SCIPsepaExec(
   SEPA*            sepa,               /**< separator */
   const SET*       set,                /**< global SCIP settings */
   SEPASTORE*       sepastore,          /**< separation storage */
   int              actdepth,           /**< depth of active node */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(sepa != NULL);
   assert(sepa->sepaexec != NULL);
   assert(sepa->freq >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(actdepth >= 0);
   assert(result != NULL);

   if( (actdepth == 0 && sepa->freq == 0) || (sepa->freq > 0 && actdepth % sepa->freq == 0) )
   {
      int oldncutsfound;

      debugMessage("executing separator <%s>\n", sepa->name);

      oldncutsfound = SCIPsepastoreGetNCutsFound(sepastore);

      /* start timing */
      SCIPclockStart(sepa->clock, set);

      /* call external separation method */
      CHECK_OKAY( sepa->sepaexec(set->scip, sepa, result) );

      /* stop timing */
      SCIPclockStop(sepa->clock, set);

      /* evaluate result */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_SEPARATED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_CONSADDED
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         char s[MAXSTRLEN];
         sprintf(s, "execution method of separator <%s> returned invalid result <%d>", 
            sepa->name, *result);
         errorMessage(s);
         return SCIP_INVALIDRESULT;
      }
      if( *result != SCIP_DIDNOTRUN )
         sepa->ncalls++;

      sepa->ncutsfound += SCIPsepastoreGetNCutsFound(sepastore) - oldncutsfound;
   }
   else
      *result = SCIP_DIDNOTRUN;
   
   return SCIP_OKAY;
}


/** gets name of separator */
const char* SCIPsepaGetName(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->name;
}

/** gets user data of separator */
SEPADATA* SCIPsepaGetData(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->sepadata;
}

/** sets user data of separator; user has to free old data in advance! */
void SCIPsepaSetData(
   SEPA*            sepa,               /**< separator */
   SEPADATA*        sepadata            /**< new separator user data */
   )
{
   assert(sepa != NULL);

   sepa->sepadata = sepadata;
}

/** gets frequency of separator */
int SCIPsepaGetFreq(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->freq;
}

/** gets time in seconds used in this separator */
Real SCIPsepaGetTime(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return SCIPclockGetTime(sepa->clock);
}

/** gets the number of times, the separator was called and tried to find a solution */
Longint SCIPsepaGetNCalls(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncalls;
}

/** gets the number of cutting planes found by this separator */
Longint SCIPsepaGetNCutsFound(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncutsfound;
}

/** is separator initialized? */
Bool SCIPsepaIsInitialized(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->initialized;
}
