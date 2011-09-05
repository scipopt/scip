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

/**@file   sepa_closecuts.c
 * @ingroup SEPARATORS
 * @brief  closecuts meta separator
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_closecuts.h"


#define SEPA_NAME              "closecuts"
#define SEPA_DESC              "closecuts meta separator"
#define SEPA_PRIORITY           1000000
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */


/* default values for parameters */
#define SCIP_DEFAULT_SEPARELINT              TRUE /**< generate close cuts w.r.t. relative interior point (best solution otherwise)? */
#define SCIP_DEFAULT_SEPACOMBVALUE           0.30 /**< convex combination value for close cuts */
#define SCIP_DEFAULT_SEPAROOTONLY            TRUE /**< generate close cuts in the root only? */
#define SCIP_DEFAULT_SEPATHRESHOLD             50 /**< threshold on number of generated cuts below which the ordinary separation is started */
#define SCIP_DEFAULT_INCLOBJCUTOFF          FALSE /**< include the objective cutoff when computing the relative interior? */
#define SCIP_DEFAULT_RECOMPUTERELINT        FALSE /**< recompute relative interior in each separation call? */
#define SCIP_DEFAULT_RELINTNORMTYPE           'o' /**< type of norm to use when computing relative interior */
#define SCIP_DEFAULT_MAXUNSUCCESSFULL           0 /**< turn off separation in current node after unsuccessful calls (-1 never turn off) */


/** separator data */
struct SCIP_SepaData
{
   SCIP_Bool             separelint;         /**< generate close cuts w.r.t. relative interior point (best solution otherwise)? */
   SCIP_Bool             separootonly;       /**< generate close cuts in the root only? */
   SCIP_Real             sepacombvalue;      /**< convex combination value for close cuts */
   int                   sepathreshold;      /**< threshold on number of generated cuts below which the ordinary separation is started */
   SCIP_Bool             inclobjcutoff;      /**< include the objective cutoff when computing the relative interior? */
   SCIP_Bool             recomputerelint;    /**< recompute relative interior in each separation call? */
   char                  relintnormtype;     /**< type of norm to use when computing relative interior */
   int                   maxunsuccessful;    /**< turn off separation in current node after unsuccessful calls (-1 never turn off) */
   SCIP_SOL*             sepasol;            /**< solution that can be used for generating close cuts */
   SCIP_Longint          discardnode;        /**< number of node for which separation is discarded */
   int                   nunsuccessful;      /**< number of consecutive unsuccessful calls */
};




/** generate point for close cut separation
 *
 *  The constructed point is the convex combination of the point stored in set->closesol and the
 *  current LP solution. The convexity parameter is set->sepa_closecombvalue. If this parameter is
 *  0, the point coincides with the LP solution.
 */
static
SCIP_RETCODE generateCloseCutPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL**            point               /**< point to be generated (or NULL if unsuccessful) */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real alpha;
   SCIP_Real onealpha;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( point != NULL );

   *point = NULL;
   if ( sepadata->sepasol == NULL )
      return SCIP_OKAY;

   alpha = sepadata->sepacombvalue;
   if ( alpha < 0.001 )
      return SCIP_OKAY;
   onealpha = 1.0 - alpha;

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, point, NULL) );

   /* generate convex combination */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for (i = 0; i < nvars; ++i)
   {
      var = vars[i];
      val = alpha * SCIPgetSolVal(scip, sepadata->sepasol, var) + onealpha * SCIPvarGetLPSol(var);

      if ( ! SCIPisZero(scip, val) )
      {
         SCIP_CALL( SCIPsetSolVal(scip, *point, var, val) );
      }
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of separator
 */


/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyClosecuts)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaClosecuts(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeClosecuts)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** initialization method of separator (called after problem was transformed) */
#define sepaInitClosecuts NULL

/** deinitialization method of separator (called before transformed problem is freed) */
#define sepaExitClosecuts NULL

/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolClosecuts NULL

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolClosecuts)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   if ( sepadata->separelint && sepadata->sepasol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &sepadata->sepasol) );
   }

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpClosecuts)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_Longint currentnodenumber;
   SCIP_Bool isroot;

   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   currentnodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
   if ( sepadata->discardnode == currentnodenumber )
      return SCIP_OKAY;

   isroot = FALSE;
   if (SCIPgetNNodes(scip) == 0)
      isroot = TRUE;

   /* only separate close cuts in the root if required */
   if ( sepadata->separootonly || isroot )
   {
      SCIP_SOL* point = NULL;

      SCIPdebugMessage("Separation method of closecuts separator.\n");
      *result = SCIP_DIDNOTFIND;

      /* check whether we have to compute a relative interior point */
      if ( sepadata->separelint )
      {
         /* check if previous relative interior point should be forgotten,
          * otherwise it is computed only once and the same point is used for all nodes */
         if ( sepadata->recomputerelint && sepadata->sepasol != NULL )
         {
            SCIP_CALL( SCIPfreeSol(scip, &sepadata->sepasol) );
         }
         if ( sepadata->sepasol == NULL )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Computing relative interior point (norm type: %c) ...\n", sepadata->relintnormtype);
            assert(sepadata->relintnormtype == 'o' || sepadata->relintnormtype == 'i');
            if( sepadata->relintnormtype == 'o' )
            {
               SCIP_CALL( SCIPcomputeLPRelIntPointOneNorm(scip, TRUE, sepadata->inclobjcutoff, &sepadata->sepasol) );
            }
            else
            {
               SCIP_CALL( SCIPcomputeLPRelIntPointSupNorm(scip, sepadata->inclobjcutoff, &sepadata->sepasol) );
            }
         }
      }
      else
      {
         /* get best solution (NULL if not present) */
         sepadata->sepasol = SCIPgetBestSol(scip);
      }

      /* separate close cuts */
      if ( sepadata->sepasol != NULL )
      {
         SCIPdebugMessage("Generating close cuts ... (combination value: %f)\n", sepadata->sepacombvalue);

         /* generate point to be separated */
         SCIP_CALL( generateCloseCutPoint(scip, sepadata, &point) );

         /* apply a separation round to generated point */
         if ( point != NULL )
         {
            int noldcuts;
            SCIP_Bool delayed;
            SCIP_Bool cutoff;

            noldcuts = SCIPgetNCuts(scip);

            SCIP_CALL( SCIPseparateSol(scip, point, isroot, FALSE, &delayed, &cutoff) );

            SCIP_CALL( SCIPfreeSol(scip, &point) );
            assert( point == NULL );

            /* the cuts can be not violated by the current LP if the computed point is strange */
            SCIP_CALL( SCIPremoveInefficaciousCuts(scip) );

            if ( cutoff )
               *result = SCIP_CUTOFF;
            else
            {
               if ( SCIPgetNCuts(scip) - noldcuts > sepadata->sepathreshold )
               {
                  sepadata->nunsuccessful = 0;
                  *result = SCIP_NEWROUND;
               }
               else
               {
                  if ( SCIPgetNCuts(scip) > noldcuts )
                  {
                     sepadata->nunsuccessful = 0;
                     *result = SCIP_SEPARATED;
                  }
                  else
                     ++sepadata->nunsuccessful;
               }
            }

            SCIPdebugMessage("Separated close cuts: %d (enoughcuts: %d, unsuccessfull: %d).\n", SCIPgetNCuts(scip) - noldcuts,
               SCIPgetNCuts(scip) - noldcuts > sepadata->sepathreshold, sepadata->nunsuccessful);

            if ( sepadata->maxunsuccessful >= 0 && sepadata->nunsuccessful > sepadata->maxunsuccessful )
            {
               SCIPdebugMessage("Turn off close cut separation, because of %d unsuccessfull calls.\n", sepadata->nunsuccessful);
               sepadata->discardnode = currentnodenumber;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** arbitrary primal solution separation method of separator */
#define sepaExecsolClosecuts NULL




/*
 * separator specific interface methods
 */

/** creates the closecuts separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaClosecuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create closecuts separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );
   sepadata->sepasol = NULL;
   sepadata->discardnode = -1;
   sepadata->nunsuccessful = 0;

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaCopyClosecuts, sepaFreeClosecuts, sepaInitClosecuts, sepaExitClosecuts,
         sepaInitsolClosecuts, sepaExitsolClosecuts, sepaExeclpClosecuts, sepaExecsolClosecuts,
         sepadata) );

   /* add closecuts separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/closecuts/separelint",
         "generate close cuts w.r.t. relative interior point (best solution otherwise)?",
         &sepadata->separelint, TRUE, SCIP_DEFAULT_SEPARELINT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/closecuts/sepacombvalue",
         "convex combination value for close cuts",
         &sepadata->sepacombvalue, TRUE, SCIP_DEFAULT_SEPACOMBVALUE, 0.0, 1.0,
         NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/closecuts/separootonly",
         "generate close cuts in the root only?",
         &sepadata->separootonly, TRUE, SCIP_DEFAULT_SEPAROOTONLY, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/closecuts/closethres",
         "threshold on number of generated cuts below which the ordinary separation is started",
         &sepadata->sepathreshold, TRUE, SCIP_DEFAULT_SEPATHRESHOLD, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/closecuts/inclobjcutoff",
         "include an objective cutoff when computing the relative interior?",
         &sepadata->inclobjcutoff, TRUE, SCIP_DEFAULT_INCLOBJCUTOFF, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/closecuts/recomputerelint",
         "recompute relative interior point in each separation call?",
         &sepadata->recomputerelint, TRUE, SCIP_DEFAULT_RECOMPUTERELINT, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip,
         "separating/closecuts/relintnormtype",
         "type of norm to use when computing relative interior: 'o'ne norm, 'i'nfinity norm",
         &sepadata->relintnormtype, TRUE, SCIP_DEFAULT_RELINTNORMTYPE, "oi", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/closecuts/maxunsuccessful",
         "turn off separation in current node after unsuccessful calls (-1 never turn off)",
         &sepadata->maxunsuccessful, TRUE, SCIP_DEFAULT_MAXUNSUCCESSFULL, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
