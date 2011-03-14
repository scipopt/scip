/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "probdata_healthcare.h"
#include "scip/cons_setppc.h"
#include "scip/misc.h"


struct SCIP_ProbData
{
   SCIP_CONS**           cons_servejobs;
   SCIP_CONS**           cons_workers;
   int                   njobs;
   int                   nworkers;
};


/** generates problem data */
static
SCIP_RETCODE createProbdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, probdata) );

   /* read problem data */
   (*probdata)->njobs = 17;
   (*probdata)->nworkers = 4;

   /* initialize variables and constraints */
   (*probdata)->cons_servejobs = NULL;
   (*probdata)->cons_workers = NULL;

   return SCIP_OKAY;
}


static
SCIP_DECL_PROBDELORIG(probdelorigHealthcare)
{
   int i;

   assert(probdata != NULL);
   assert(*probdata != NULL);

   /* release and free constraints */
   for( i = 0; i < (*probdata)->njobs; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->cons_servejobs[i]) );
   }
   SCIPfreeMemoryArray(scip, &(*probdata)->cons_servejobs);

   for( i = 0; i < (*probdata)->nworkers; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->cons_workers[i]) );
   }
   SCIPfreeMemoryArray(scip, &(*probdata)->cons_workers);

   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


#define probtransHealthcare NULL
#define probdeltransHealthcare NULL
#define probinitsolHealthcare NULL
#define probexitsolHealthcare NULL
#define probcopyHealthcare NULL


SCIP_RETCODE HCPcreateProbHealthcare(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< problem name */
   )
{
   SCIP_PROBDATA* probdata = NULL;

   /* generate problem data */
   SCIP_CALL( createProbdata(scip, &probdata) );

   SCIP_CALL( SCIPcreateProb(scip, name, probdelorigHealthcare, probtransHealthcare, probdeltransHealthcare, 
         probinitsolHealthcare, probexitsolHealthcare, probcopyHealthcare, probdata) );

   return SCIP_OKAY;
}

SCIP_RETCODE HCPgenerateModel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   int i;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* generate constraints */
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->cons_servejobs, probdata->njobs) );
   for( i = 0; i < probdata->njobs; ++i )
   {
      char consname[SCIP_MAXSTRLEN];

      SCIPsnprintf(consname, SCIP_MAXSTRLEN, "servejobs_%d", i);
      SCIP_CALL( SCIPcreateConsSetpart(scip, &probdata->cons_servejobs[i], consname, 
            0, NULL, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, probdata->cons_servejobs[i]) );
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->cons_workers, probdata->nworkers) );
   for( i = 0; i < probdata->nworkers; ++i )
   {
      char consname[SCIP_MAXSTRLEN];

      SCIPsnprintf(consname, SCIP_MAXSTRLEN, "workers_%d", i);
      SCIP_CALL( SCIPcreateConsSetpack(scip, &probdata->cons_workers[i], consname, 
            0, NULL, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, probdata->cons_workers[i]) );
   }

   return SCIP_OKAY;
}

SCIP_CONS** HCPgetConsServejobs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->cons_servejobs;
}

SCIP_CONS** HCPgetConsWorkers(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->cons_workers;
}

int HCPgetNJobs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->njobs;
}

int HCPgetNWorkers(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->nworkers;
}
