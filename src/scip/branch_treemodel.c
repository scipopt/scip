/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_treemodel.c
 * @brief  Branching rules based on the Single-Variable-Branching (SVB) model
 * @author Daniel Anderson
 * @author Pierre Le Bodic
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/branch_treemodel.h"

#include <limits.h>

#define DEFAULT_ENABLE         TRUE     /**< should candidate branching variables be scored using the Treemodel rule? */
#define DEFAULT_HIGHRULE       'r'      /**< scoring function to use at nodes predicted to be high in the tree. ('d'efault, 's'vts, 'r'atio, 't'ree sample) */
#define DEFAULT_LOWRULE        'd'      /**< scoring function to use at nodes predicted to be low in the tree ('d'efault, 's'vts, 'r'atio, 't'ree sample) */
#define DEFAULT_HEIGHT         10       /**< estimated tree height at which we switch from using the low rule to the high rule */
#define DEFAULT_FILTERHIGH     'a'      /**< should dominated candidates be filtered before using the high scoring function? ('a'uto, 't'rue, 'f'alse) */
#define DEFAULT_FILTERLOW      'a'      /**< should dominated candidates be filtered before using the low scoring function? ('a'uto, 't'rue, 'f'alse) */
#define DEFAULT_MAXFPITER      24       /**< maximum number of fixed-point iterations when computing the ratio */
#define DEFAULT_MAXSVTSHEIGHT  100      /**< maximum height to compute the SVTS score exactly before approximating */
#define DEFAULT_FALLBACKINF    'r'      /**< which method should be used as a fallback if the tree size estimates are infinite? ('d'efault, 'r'atio) */
#define DEFAULT_FALLBACKNOPRIM 'r'      /**< which method should be used as a fallback if there is no primal bound available? ('d'efault, 'r'atio) */

/** Parameters required by the Treemodel branching rules */
struct SCIP_BranchTreemodel
{
   SCIP_Bool             enabled;            /**< should candidate branching variables be scored using the Treemodel rule? */
   char                  highrule;           /**< scoring function to use at nodes predicted to be high in the tree. ('d'efault, 's'vts, 'r'atio, 't'ree sample) */
   char                  lowrule;            /**< scoring function to use at nodes predicted to be low in the tree ('d'efault, 's'vts, 'r'atio, 't'ree sample) */
   int                   height;             /**< estimated tree height at which we switch from using the low rule to the high rule */
   char                  filterhigh;         /**< should dominated candidates be filtered before using the high scoring function? ('a'uto, 't'rue, 'f'alse) */
   char                  filterlow;          /**< should dominated candidates be filtered before using the low scoring function? ('a'uto, 't'rue, 'f'alse) */
   int                   maxfpiter;          /**< maximum number of fixed-point iterations when computing the ratio */
   int                   maxsvtsheight;      /**< maximum height to compute the SVTS score exactly before approximating */
   char                  fallbackinf;        /**< which method should be used as a fallback if the tree size estimates are infinite? ('d'efault, 'r'atio) */
   char                  fallbacknoprim;     /**< which method should be used as a fallback if there is no primal bound available? ('d'efault, 'r'atio) */
};

/** Initialises the Treemodel parameter data structure */
SCIP_RETCODE SCIPtreemodelInit(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_BRANCHTREEMODEL**  treemodel   /**< Treemodel parameter data structure */
   )
{

   SCIP_CALL( SCIPallocBlockMemory(scip, treemodel) );

   assert(treemodel != NULL && *treemodel != NULL);

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/treemodel/enable",
         "should candidate branching variables be scored using the Treemodel branching rules?", &(*treemodel)->enabled, FALSE, DEFAULT_ENABLE,
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/highrule",
         "scoring function to use at nodes predicted to be high in the tree ('d'efault, 's'vts, 'r'atio, 't'ree sample)", &(*treemodel)->highrule, FALSE, DEFAULT_HIGHRULE, "dsrt",
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/lowrule",
         "scoring function to use at nodes predicted to be low in the tree ('d'efault, 's'vts, 'r'atio, 't'ree sample)", &(*treemodel)->lowrule, FALSE, DEFAULT_LOWRULE, "dsrt",
         NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/treemodel/height",
         "estimated tree height at which we switch from using the low rule to the high rule", &(*treemodel)->height, FALSE, DEFAULT_HEIGHT, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/filterhigh",
         "should dominated candidates be filtered before using the high scoring function? ('a'uto, 't'rue, 'f'alse)", &(*treemodel)->filterhigh, TRUE, DEFAULT_FILTERHIGH, "atf",
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/filterlow",
         "should dominated candidates be filtered before using the low scoring function? ('a'uto, 't'rue, 'f'alse)", &(*treemodel)->filterlow, TRUE, DEFAULT_FILTERLOW, "atf",
         NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/treemodel/maxfpiter",
         "maximum number of fixed-point iterations when computing the ratio", &(*treemodel)->maxfpiter, TRUE, DEFAULT_MAXFPITER, 1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/treemodel/maxsvtsheight",
         "maximum height to compute the SVTS score exactly before approximating", &(*treemodel)->maxsvtsheight, TRUE, DEFAULT_MAXSVTSHEIGHT, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/fallbackinf",
         "which method should be used as a fallback if the tree size estimates are infinite? ('d'efault, 'r'atio)", &(*treemodel)->fallbackinf, TRUE, DEFAULT_FALLBACKINF, "dr",
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/fallbacknoprim",
         "which method should be used as a fallback if there is no primal bound available? ('d'efault, 'r'atio)", &(*treemodel)->fallbacknoprim, TRUE, DEFAULT_FALLBACKNOPRIM, "dr",
         NULL, NULL) );
   
   return SCIP_OKAY;
}

/** Frees the Treemodel parameter data structure */
SCIP_RETCODE SCIPtreemodelFree(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_BRANCHTREEMODEL**  treemodel   /**< Treemodel parameter data structure */
   )
{
   assert(treemodel != NULL && *treemodel != NULL);

   SCIPfreeBlockMemory(scip, treemodel);
   return SCIP_OKAY;
}