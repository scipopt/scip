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

/**@file   JniScipHeurRens.c
 * @ingroup PUBLICMETHODS
 * @brief  JNI SCIP heur rens callable library
 * @author Stefan Heinz
 * @author Alexandra Kraft
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "JniScipHeurRens.h"
#include "def.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_rens.h"

#include <string.h>

/** creates RENS primal heuristic and includes it in SCIP */
JNIEXPORT
void JNISCIPHEURRENS(includeHeurRens)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip               /**< SCIP data structure */
   )
{
   SCIP* scip;

   /* convert JNI pointer into C pointer */
   scip = (SCIP*) (size_t) jscip;
   assert(scip != NULL);

   JNISCIP_CALL( SCIPincludeHeurRens(scip) );
}

/** main procedure of the RENS heuristic, creates and solves a sub-SCIP */
JNIEXPORT
void JNISCIPHEURRENS(applyRens)(
   JNIEnv*               env,                /**< JNI environment variable */
   jobject               jobj,               /**< JNI class pointer */
   jlong                 jscip,              /**< SCIP data structure */
   jlong                 jheur,              /**< heuristic data structure                                            */
   jintArray             jresult,            /**< result data structure                                               */
   jdouble               jminfixingrate,     /**< minimum percentage of integer variables that have to be fixed       */
   jdouble               jminimprove,        /**< factor by which RENS should at least improve the incumbent          */
   jlong                 jmaxnodes,          /**< maximum number of  nodes for the subproblem                         */
   jlong                 jnstallnodes,       /**< number of stalling nodes for the subproblem                         */
   jchar                 jstartsol,          /**< solution used for fixing values ('l'p relaxation, 'n'lp relaxation) */
   jboolean              jbinarybounds,      /**< should general integers get binary bounds [floor(.),ceil(.)]?       */
   jboolean              juselprows          /**< should subproblem be created out of the rows in the LP rows?        */
   )
{
   SCIPerrorMessage("method applyRens is not implemented yet (deprecated)\n");
   JNISCIP_CALL( SCIP_ERROR );
}
