/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_minor.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  minor separator
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_minor.h"


#define SEPA_NAME              "minor"
#define SEPA_DESC              "separator template"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXMINORS           100 /**< default maximum number for minors (0: no limit) */

/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_VAR**            minors;             /**< variables of 2x2 minors; each minor is stored like (auxvar_x^2,auxvar_y^2,auxvar_xy) */
   int                   nminors;            /**< total number of minors */
   int                   minorssize;         /**< size of minors array */
   int                   maxminors;          /**< maximum number for minors (0: no limit) */
   SCIP_Bool             detectedminors;     /**< has the minor detection beeing called? */
};

/*
 * Local methods
 */

/** helper method to clear separation data */
static
SCIP_RETCODE sepadataClear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   int i;

   assert(sepadata != NULL);

   /* release captured variables */
   for( i = 0; i < 3 * sepadata->nminors; ++i )
   {
      assert(sepadata->minors[i] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &sepadata->minors[i]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &sepadata->minors, sepadata->minorssize);

   /* reset counters */
   sepadata->nminors = 0;
   sepadata->minorssize = 0;

   return SCIP_OKAY;
}

/** method to detect and store minors */
static
SCIP_RETCODE detectMinors(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   SCIP_SEPADATA* sepadata;

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* check whether minor detection has been called already */
   if( sepadata->detectedminors )
      return SCIP_OKAY;

   assert(sepadata->minors == NULL);
   assert(sepadata->nminors == 0);

   /* TODO check whether auxiliary variables have been already generated */

   sepadata->detectedminors = TRUE;

   return SCIP_OKAY;
}

/** separates cuts for stored minors */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_Bool             allowlocal,         /**< should local cuts be allowed */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_SEPADATA* sepadata;

   assert(sepa != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* check whether there are some minors available */
   if( sepadata->nminors == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyMinor)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaMinor(scip) );

   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->minors == NULL);
   assert(sepadata->nminors == 0);
   assert(sepadata->minorssize == 0);

   /* free separator data */
   SCIPfreeBlockMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitMinor NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitMinor NULL
#endif


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
static
SCIP_DECL_SEPAINITSOL(sepaInitsolMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* TODO */

   return SCIP_OKAY;
}


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* clear separation data */
   SCIP_CALL( sepadataClear(scip, sepadata) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMinor)
{  /*lint --e{715}*/

   /* try to detect minors */
   SCIP_CALL( detectMinors(scip, sepa) );

   /* call separation method */
   SCIP_CALL( separatePoint(scip, sepa, NULL, allowlocal, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolMinor)
{  /*lint --e{715}*/

   /* try to detect minors */
   SCIP_CALL( detectMinors(scip, sepa) );

   /* call separation method */
   SCIP_CALL( separatePoint(scip, sepa, sol, allowlocal, result) );

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the minor separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMinor(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata = NULL;
   SCIP_SEPA* sepa = NULL;

   /* create minor separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   BMSclearMemory(sepadata);

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpMinor, sepaExecsolMinor,
         sepadata) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyMinor) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeMinor) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitMinor) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitMinor) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolMinor) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolMinor) );

   /* add minor separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxminors",
         "maximum number for minors (0: no limit)",
         &sepadata->maxminors, TRUE, DEFAULT_MAXMINORS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
